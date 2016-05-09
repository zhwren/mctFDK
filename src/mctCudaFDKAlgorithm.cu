
#include "mctCudaFDKAlgorithm.h"



texture<float, cudaTextureType3D,cudaReadModeElementType> Prj_tex;

__constant__ float cSrcAxisDist,cSrcDetDist;

__constant__ int cDetColumn,cDetRow;
__constant__ float cBinWidth,cBinHeight;

__constant__ int cReconWidth,cReconHeight,cReconSlice;
__constant__ float cPixSizeX,cPixSizeY, cPixSizeZ, cFOV,cCoef;

__constant__ int cPrjNum;

__constant__ float cInvBinWidth, cInvBinHeight;
__constant__ float cReconWidthMid, cReconHeightMid, cReconSliceMid;
__constant__ float cDetColumnRayCenter, cDetRowRayCenter; 
__constant__ float cSDDInvBinHeight,cReconZminSDDInvBinHeight, cSDDInvBinHeightSliceSpace; 

static const unsigned g_MaxAngles = 1440;
__constant__ float gC_angle_sin[g_MaxAngles];
__constant__ float gC_angle_cos[g_MaxAngles];

__global__ void Correction(float *dSinoData,int iPaddedDetCount,float *dCorrectionMatrix,int DetColumn,int DetRow,int PrjSize)
{
	int i = threadIdx.x + mul24(blockDim.x,blockIdx.x);
	int j = threadIdx.y + mul24(blockDim.y,blockIdx.y);

	if (i < DetColumn && j < DetRow*PrjSize)
	{
		dSinoData[j*iPaddedDetCount + i] *= dCorrectionMatrix[(j%DetRow)*DetColumn + i];
	}
}

__global__ void Filter_FFT(cufftComplex *dFilteredSinoData,cufftComplex *dRampfilter,int DetRow,int iHalfFFTSize)
{
	int i = threadIdx.x + mul24(blockDim.x,blockIdx.x);
	int j = threadIdx.y + mul24(blockDim.y,blockIdx.y);

	if (i < iHalfFFTSize && j < DetRow)
	{
		dFilteredSinoData[j*iHalfFFTSize + i].x *= dRampfilter[i].x*cCoef;
		dFilteredSinoData[j*iHalfFFTSize + i].y *= dRampfilter[i].x*cCoef;
	}
}

__global__ void BP (float *dReconData, size_t Pitch, int iBatch)
{
	int i = threadIdx.x + mul24(blockDim.x,blockIdx.x);
	int j = threadIdx.y + mul24(blockDim.y,blockIdx.y);
	int k = threadIdx.z + mul24(blockDim.z,blockIdx.z);

	float x = (i-cReconWidthMid)*cPixSizeX;
	float y = (j-cReconHeightMid)*cPixSizeY;

	int p = iBatch*cPrjNum + k;

	if ( x*x + y*y < cFOV*cFOV)
	{
		float t = x*gC_angle_cos[p]+y*gC_angle_sin[p];
		float s = cSrcAxisDist + x*gC_angle_sin[p]-y*gC_angle_cos[p];

		float L2 = t*t+s*s;

		float m = atanf(t/s)*cInvBinWidth + cDetColumnRayCenter;
		float n = rsqrtf(L2)*cReconZminSDDInvBinHeight + cDetRowRayCenter;

		float dn = rsqrtf(L2)*cSDDInvBinHeightSliceSpace;

	#pragma unroll
		for (int iz = 0; iz < cReconSlice; iz++)
		{ 
			atomicAdd(dReconData+ iz*Pitch*cReconHeight +j*Pitch +i, tex3D(Prj_tex,m,n,k+0.5f)/L2);
			n += dn;
		}
	}
}

int calcZeropadFactorPowerOfTwo(int n, int iZeropad)
{
	if (iZeropad > 0)
	{
		double logBase2 = log((double)n) / log((double)2);
		int nextPowerOf2 = static_cast<int>(floor(logBase2));

		if (logBase2 != floor(logBase2))
			nextPowerOf2++;

		nextPowerOf2 += (iZeropad - 1);
		n = 1 << nextPowerOf2;
	}

	return n;
}

namespace mct
{
const int CFDK::m_nBatch = 5;//4;

CFDK::CFDK()
{
	this->InitParams();
}

CFDK::~CFDK(void)
{
	this->FreeObjects();
}

void CFDK::InitParams()
{
	this->m_PrjArray = NULL;
	this->m_PrjArrayLen = 0;

	this->m_DCorrectionMatrix = NULL;

	this->m_DReconData = NULL;
	this->m_DReconDataLenSize = 0;
	this->m_DReconDataLenCount = 0;

	this->m_DsinoData = NULL;
	this->m_DsinoDataLen = 0;

	this->m_DFilteredsinoData = NULL;
	this->m_DFilteredsinoDataLen = 0;

	this->m_RampFilter = NULL;
	this->m_RampFilterLen = 0;

	this->hReconData = NULL;
	this->hReconDataLen = 0;

	this->m_ProjectionAngleCountMatix = 0;
	this->m_ProjectionAngleCountRampFilter = 0;
	this->m_ProjectionAngleCountAngle = 0;
	this->m_DetectorSpacingMatix = 0;
	this->m_DetectorSpacingRampFilter = 0;

	this->m_FFTLen = 0;
	this->m_iPaddedDetCountOld =0;
}

void CFDK::FreeObjects()
{
	cudaFreeArray(m_PrjArray);
	m_PrjArray = NULL;

	cudaFree(m_DCorrectionMatrix);
	m_DCorrectionMatrix = NULL;

	cudaFree(m_DReconData); 
	m_DReconData = NULL;
	
	cudaFree(m_DsinoData);
	m_DsinoData = NULL;

	cudaFree(m_DFilteredsinoData);
	m_DFilteredsinoData = NULL;

	cudaFree(m_RampFilter);
	m_RampFilter = NULL;

	cufftDestroy(m_FwdFFT);
	cufftDestroy(m_BwdFFT);
	this->m_FFTLen = 0;
	this->m_iPaddedDetCountOld =0;

	if(this->hReconData != NULL)
	{
		delete[] this->hReconData;
		this->hReconData = NULL;
		this->hReconDataLen = 0;
	}
}

bool CFDK::SetParams(ScannerGeometry scanGeometry, ProjectionParams prjParams, ReconstructionParams reconParams, float* hSinoData, int iGPUIndex)
{
	SetPrjGeometry(scanGeometry, prjParams, reconParams);

	if(!setGPUIndex(iGPUIndex))
	{
		return false;
	}

	if(!CpyToSymbol())
	{
		return false;
	}
	if(!allocateBuffers())
	{
		return false;
	}
	
	if(!genRampFilter())
	{
		return false;
	}
	if(!caculateCorrectMatix())
	{
		return false;
	}

	m_hSinoData = hSinoData;

	return true;
}

void CFDK::SetPrjGeometry(ScannerGeometry scanGeometry, ProjectionParams prjParams, ReconstructionParams reconParams)
{
	m_SourceToIsoCenterDistance = scanGeometry.m_SourceToIsoCenterDistance;
	m_SourceToDetectorDistance = scanGeometry.m_SourceToDetectorDistance;

	m_DetectorSpacingX = scanGeometry.m_DetectorSpacingX;
	m_DetectorSpacingY = scanGeometry.m_DetectorSpacingY;

	m_DetectorColumnCount = scanGeometry.m_DetectorColumnCount + scanGeometry.m_DetectorCount - 1;
	m_DetectorRowCount = scanGeometry.m_DetectorRowCount;

	m_DetectorColumnRayCenter = scanGeometry.m_DetectorColumnRayCenter;
	m_DetectorRowRayCenter = scanGeometry.m_DetectorRowCount/2-0.5f; //scanGeometry.m_DetectorRowRayCenter;
	
	m_ProjectionAngleCount = prjParams.m_ProjectionAngleCount;

	m_ProjectionAngleStart = prjParams.m_ProjectionAngleStart;
	m_ProjectionAngleStep = -2*PI/prjParams.m_ProjectionAnglesPerRotation;

	m_DetectorLengthX = scanGeometry.m_DetectorSpacingX*m_DetectorColumnCount;
	m_DetectorLengthY = scanGeometry.m_DetectorSpacingY*m_DetectorRowCount;

	m_fFOV = m_SourceToIsoCenterDistance*sin(0.5f*(m_DetectorLengthX-m_DetectorSpacingX));

	m_ReconColumnCount = reconParams.m_ReconColumnCount;
	m_ReconRowCount = reconParams.m_ReconRowCount;
	m_ReconSliceCount = reconParams.m_ReconSliceCount*reconParams.m_MergedNum; //调整为重建所有层，输出根据MergedNum进行输出  2015.12.16

	m_nMergeNum = reconParams.m_MergedNum;

	m_ReconWindowMidColumn = reconParams.m_ReconWindowMidColumn; 
	m_ReconWindowMidRow = reconParams.m_ReconWindowMidRow;
	m_ReconWindowMidSlice = reconParams.m_ReconWindowMidSlice;

	m_PixelSpacingX = reconParams.m_PixelSpacingX;
	m_PixelSpacingY = reconParams.m_PixelSpacingY;
	m_PixelSpacingZ = reconParams.m_PixelSpacingZ;
	
	m_nPrjBatchSize = m_ProjectionAngleCount/m_nBatch;

	m_iPaddedDetCount = calcZeropadFactorPowerOfTwo(2*m_DetectorColumnCount-1, 1); 
	m_iHalfFFTSize = (m_iPaddedDetCount/2 + 1);
}


bool CFDK::CpyToSymbol()
{
	//可优化  2015.11.26 
	cudaMemcpyToSymbol(cSrcAxisDist,&m_SourceToIsoCenterDistance,sizeof(float));
	cudaMemcpyToSymbol(cSrcDetDist,&m_SourceToDetectorDistance,sizeof(float));

	cudaMemcpyToSymbol(cBinWidth,&m_DetectorSpacingX,sizeof(float));
	cudaMemcpyToSymbol(cBinHeight,&m_DetectorSpacingY,sizeof(float));

	cudaMemcpyToSymbol(cPixSizeX,&m_PixelSpacingX,sizeof(float));
	cudaMemcpyToSymbol(cPixSizeY,&m_PixelSpacingY,sizeof(float));
	cudaMemcpyToSymbol(cPixSizeZ,&m_PixelSpacingZ,sizeof(float));

	cudaMemcpyToSymbol(cFOV,&m_fFOV,sizeof(float));
	cudaMemcpyToSymbol(cDetColumn,&m_DetectorColumnCount,sizeof(int));
	cudaMemcpyToSymbol(cDetRow,&m_DetectorRowCount,sizeof(int));

	cudaMemcpyToSymbol(cReconWidth,&m_ReconColumnCount,sizeof(int));
	cudaMemcpyToSymbol(cReconHeight,&m_ReconRowCount,sizeof(int));
	cudaMemcpyToSymbol(cReconSlice,&m_ReconSliceCount,sizeof(int));

	cudaMemcpyToSymbol(cPrjNum,&m_nPrjBatchSize,sizeof(int));

	cudaMemcpyToSymbol(cReconWidthMid,&m_ReconWindowMidColumn,sizeof(float));
	cudaMemcpyToSymbol(cReconHeightMid,&m_ReconWindowMidRow,sizeof(float));
	cudaMemcpyToSymbol(cReconSliceMid,&m_ReconWindowMidSlice,sizeof(float));

	cudaMemcpyToSymbol(cDetColumnRayCenter,&m_DetectorColumnRayCenter,sizeof(float));
	cudaMemcpyToSymbol(cDetRowRayCenter,&m_DetectorRowRayCenter,sizeof(float));

	float InvBinWidth = 1.f/m_DetectorSpacingX;
	float InvBinHeight = 1.f/m_DetectorSpacingY;

	cudaMemcpyToSymbol(cInvBinWidth,&InvBinWidth,sizeof(float));
	cudaMemcpyToSymbol(cInvBinHeight,&InvBinHeight,sizeof(float)); 

	float coef = 1.f/m_iPaddedDetCount*abs(m_DetectorSpacingX*m_ProjectionAngleStep*m_SourceToIsoCenterDistance);
	cudaMemcpyToSymbol(cCoef,&coef,sizeof(float)); 

	float SDDInvBinHeight = m_SourceToDetectorDistance/m_DetectorSpacingY;
	float ReconZminSDDInvBinHeight = (-m_ReconWindowMidSlice)*m_PixelSpacingZ*SDDInvBinHeight;
	float SDDInvBinHeightSliceSpace = m_SourceToDetectorDistance/m_DetectorSpacingY*m_PixelSpacingZ; 

	cudaMemcpyToSymbol(cSDDInvBinHeight,&SDDInvBinHeight,sizeof(float)); 
	cudaMemcpyToSymbol(cReconZminSDDInvBinHeight,&ReconZminSDDInvBinHeight,sizeof(float)); 
	cudaMemcpyToSymbol(cSDDInvBinHeightSliceSpace,&SDDInvBinHeightSliceSpace,sizeof(float)); 

	if(this->m_ProjectionAngleCountAngle != m_ProjectionAngleCount)
	{
		float* angle_sin = new float[m_ProjectionAngleCount];
		if(angle_sin == NULL)
			return false;

		float* angle_cos = new float[m_ProjectionAngleCount];
		if(angle_cos == NULL)
		{
			delete []angle_sin;
			return false;
		}
	
		float angles = m_ProjectionAngleStart;		//TODO 目前起始角度都一致

		for (unsigned int i = 0; i < m_ProjectionAngleCount; ++i) 
		{
			angle_sin[i] = sinf(angles);
			angle_cos[i] = cosf(angles);

			angles += m_ProjectionAngleStep;
		}

		cudaMemcpyToSymbol(gC_angle_sin, angle_sin, m_ProjectionAngleCount*sizeof(float), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(gC_angle_cos, angle_cos, m_ProjectionAngleCount*sizeof(float), 0, cudaMemcpyHostToDevice);

		delete []angle_sin;
		delete []angle_cos;

		this->m_ProjectionAngleCountAngle = this->m_ProjectionAngleCount;
	}

	return true;
}

bool CFDK::setGPUIndex(int iGPUIndex)
{
	cudaSetDevice(iGPUIndex);
	cudaError_t err = cudaGetLastError();

	// Ignore errors caused by calling cudaSetDevice multiple times
	if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
		return false;

	return true;
}

bool CFDK::allocateBuffers()
{
	cudaError_t err;
	int dsinoDataLen;

	if((this->m_DReconData == NULL) || (this->m_DReconDataLenSize != m_ReconRowCount) ||(this->m_DReconDataLenCount!=m_ReconSliceCount))
	{
		if(this->m_DReconData != NULL)
		{
			cudaFree(this->m_DReconData);
			this->m_DReconData == NULL;
			this->m_DReconDataLenSize = 0;
			this->m_DReconDataLenCount = 0;
		}
		/// <summary>
		/// Allocate GPU Memory for Reconstruction Output
		/// </summary>
		err = cudaMallocPitch((void**)&this->m_DReconData, &m_nReconDataPitch, sizeof(float)*m_ReconColumnCount, m_ReconRowCount*m_ReconSliceCount);		
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);	
		this->m_DReconDataLenSize = m_ReconRowCount;
		this->m_DReconDataLenCount = m_ReconSliceCount;
	}
	else
	{
		m_nReconDataPitch = sizeof(float)*m_ReconColumnCount;
	}
	err = cudaMemset2D(this->m_DReconData, m_nReconDataPitch, 0, sizeof(float)*m_ReconColumnCount, m_ReconRowCount*m_ReconSliceCount);
	if(cudaSuccess != err)
	{
		return false;
	}
	CUDA_ASSERT(err);
	m_nReconDataPitch = m_nReconDataPitch/sizeof(float);

	dsinoDataLen = sizeof(cufftReal)*m_DetectorRowCount*m_iPaddedDetCount*m_nPrjBatchSize;
	if((this->m_DCorrectionMatrix == NULL) || (this->m_DsinoDataLen != dsinoDataLen))
	{
		if(this->m_DCorrectionMatrix != NULL)
		{
			cudaFree(this->m_DCorrectionMatrix);
			this->m_DCorrectionMatrix = NULL;
		}
		/// <summary>
		/// Allocate GPU Memory for FDK correction mattrix
		/// </summary>
		err = cudaMalloc((void**)&this->m_DCorrectionMatrix, sizeof(float)*m_DetectorColumnCount*m_DetectorRowCount);
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);
	}
	
	if((this->m_DsinoData == NULL) || (this->m_DsinoDataLen != dsinoDataLen))
	{
		if(this->m_DsinoData != NULL)
		{
			cudaFree(this->m_DsinoData);
			this->m_DsinoData = NULL;
		}
		/// <summary>
		/// Allocate GPU Memory for Sinogram Data
		/// 由于在滤波是需要补领，因此在分配内存时就特意分配大的内存
		/// </summary>
		err = cudaMalloc((void**)&this->m_DsinoData, dsinoDataLen);
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);
		this->m_DsinoDataLen = dsinoDataLen;
	}

	if((this->m_DFilteredsinoData ==NULL) || (this->m_DFilteredsinoDataLen != sizeof(cufftComplex)*m_DetectorRowCount*m_iHalfFFTSize*m_nPrjBatchSize))
	{
		if(this->m_DFilteredsinoData !=NULL)
		{			
			cudaFree(this->m_DFilteredsinoData);
			this->m_DFilteredsinoData = NULL;
			this->m_DFilteredsinoDataLen = 0;
		}
		/// <summary>
		/// Allocate Memory for Filtered Sinogram Data
		/// </summary>
		err = cudaMalloc((void**)&this->m_DFilteredsinoData, sizeof(cufftComplex)*m_DetectorRowCount*m_iHalfFFTSize*m_nPrjBatchSize);
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);
		this->m_DFilteredsinoDataLen = sizeof(cufftComplex)*m_DetectorRowCount*m_iHalfFFTSize*m_nPrjBatchSize;
	}
	if((this->m_RampFilter ==NULL) || (this->m_RampFilterLen != sizeof(cufftComplex)*m_DetectorRowCount*m_iHalfFFTSize*m_nPrjBatchSize))
	{
		if(this->m_RampFilter !=NULL)
		{			
			cudaFree(this->m_RampFilter);
			this->m_RampFilter = NULL;
			this->m_RampFilterLen = 0;
		}
		/// <summary>
		/// Allocate GPU Memory for Ramp Filter
		/// </summary>
		err = cudaMalloc((void**)&this->m_RampFilter, sizeof(cufftComplex)*m_iHalfFFTSize);
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);
		this->m_RampFilterLen = sizeof(cufftComplex)*m_DetectorRowCount*m_iHalfFFTSize*m_nPrjBatchSize;
	}

	if((this->m_PrjArray == NULL) || (this->m_PrjArrayLen != dsinoDataLen))
	{
		if(this->m_PrjArray != NULL)
		{
			cudaFreeArray(this->m_PrjArray);
			this->m_PrjArray = NULL;
			this->m_PrjArrayLen = 0;
		}
		/// <summary>
		/// Allocate GPU Memory for BackProjection Sinogram Array
		/// </summary>
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
		cudaExtent extent = make_cudaExtent(m_DetectorColumnCount,m_DetectorRowCount,m_nPrjBatchSize);
		err = cudaMalloc3DArray(&this->m_PrjArray,&channelDesc,extent);
		if(cudaSuccess != err)
		{
			return false;
		}
		CUDA_ASSERT(err);
		this->m_PrjArrayLen = dsinoDataLen;
	}
	
	if((this->m_iPaddedDetCountOld != m_iPaddedDetCount) || (this->m_FFTLen != m_DetectorRowCount*m_nPrjBatchSize))
	{
		if((this->m_iPaddedDetCountOld != 0) ||(this->m_FFTLen != 0))
		{
			//Frees all GPU resources associated with a cuFFT plan and destroys the internal plan data structure. 
			//This function should be called once a plan is no longer needed, to avoid wasting GPU memory
			cufftDestroy(m_FwdFFT);
			cufftDestroy(m_BwdFFT);
		}

		cufftPlan1d(&m_FwdFFT, m_iPaddedDetCount, CUFFT_R2C, m_DetectorRowCount*m_nPrjBatchSize);
		cufftPlan1d(&m_BwdFFT, m_iPaddedDetCount, CUFFT_C2R, m_DetectorRowCount*m_nPrjBatchSize);

		this->m_iPaddedDetCountOld = m_iPaddedDetCount;
		this->m_FFTLen = m_DetectorRowCount*m_nPrjBatchSize;
	}

	return true;
}


bool CFDK::caculateCorrectMatix()
{
	if((this->m_DetectorSpacingMatix != m_DetectorSpacingX) || (this->m_ProjectionAngleCountMatix!=this->m_ProjectionAngleCount))  //多次重建只有m_DetectorSpacingX和m_DetectorSpacingY，角度数变化，其它参数一致
	{
		float *hCorrectionMatrix = new float[m_DetectorColumnCount * m_DetectorRowCount ];
		if(hCorrectionMatrix == NULL)
			return false;

		for (size_t j = 0; j < m_DetectorRowCount; j++)
		{
			float y = (j-m_DetectorRowRayCenter)*m_DetectorSpacingY;
			float cosa = m_SourceToDetectorDistance/sqrt(m_SourceToDetectorDistance*m_SourceToDetectorDistance + y*y);

			for (size_t i = 0; i < m_DetectorColumnCount; i++)
			{
				float x = (i-m_DetectorColumnRayCenter)*m_DetectorSpacingX;
				hCorrectionMatrix[j*m_DetectorColumnCount+i] = cosa*cos(x);
			}
		}

		cudaError_t err;
		err = cudaMemcpy(m_DCorrectionMatrix, hCorrectionMatrix, sizeof(float)*m_DetectorColumnCount*m_DetectorRowCount, cudaMemcpyHostToDevice);
		CUDA_ASSERT(err);

		delete []hCorrectionMatrix;

		this->m_DetectorSpacingMatix = m_DetectorSpacingX;
		this->m_ProjectionAngleCountMatix =this->m_ProjectionAngleCount;
	}

	return true;
}

bool CFDK::genRampFilter()
{
	if((this->m_DetectorSpacingRampFilter != m_DetectorSpacingX) ||(this->m_ProjectionAngleCountRampFilter!=this->m_ProjectionAngleCount))
	{
		float *rampFilter = new float[m_iPaddedDetCount];
		if(rampFilter == NULL)
			return false;

		/// <summary>
		/// Step 1: Caculate RampFilter Spatial Domain Response
		/// </summary>
		memset(rampFilter,0,sizeof(float)*m_iPaddedDetCount);
		for (size_t i = 1;i < m_DetectorColumnCount;i += 2)
		{
			rampFilter[i] = rampFilter[m_iPaddedDetCount-i] = -1.f/(2*PI*PI*sin(i*m_DetectorSpacingX)*sin(i*m_DetectorSpacingX));
		}
		rampFilter[0] = 0.125f/(m_DetectorSpacingX*m_DetectorSpacingX);

		/// <summary>
		/// Step 2: Copy to GPU Memory
		/// </summary>
		float *DrampFilter;
		cudaError_t err;
		err = cudaMalloc((void**)&DrampFilter, sizeof(cufftReal)*m_iPaddedDetCount);
		CUDA_ASSERT(err);

		err = cudaMemcpy(DrampFilter, rampFilter, sizeof(cufftReal)*m_iPaddedDetCount, cudaMemcpyHostToDevice);
		CUDA_ASSERT(err);

		/// <summary>
		/// Step 3: FFT and get RampFilter's frequency domain spectrumn
		/// </summary>
		cufftHandle RampFilterFFT;
		cufftPlan1d(&RampFilterFFT, m_iPaddedDetCount, CUFFT_R2C, 1);
		cufftExecR2C(RampFilterFFT,DrampFilter, m_RampFilter);

		delete []rampFilter;
		cudaFree(DrampFilter);
		cufftDestroy(RampFilterFFT);

		this->m_DetectorSpacingRampFilter = m_DetectorSpacingX;
		this->m_ProjectionAngleCountRampFilter = this->m_ProjectionAngleCount;
	//genFilter(E_FBPFILTER::FILTER_RAMLAK, m_fBinWidth, m_nDetColumn, 
	//	0, E_FILTER_GENERATION::FILTER_GENERATION_INVERSE_FOURIER, m_RampFilter,
	//	1, E_SCANNER::GEOMETRY_EQUIANGULAR, m_fSrcAxisDist, m_fSrcDetDist);
	}

	return true;
}

bool CFDK::bindProjDataTexture(const cudaArray* array)
{
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	Prj_tex.addressMode[0] = cudaAddressModeClamp;
	Prj_tex.addressMode[1] = cudaAddressModeClamp;
	Prj_tex.addressMode[2] = cudaAddressModeClamp;
	Prj_tex.filterMode = cudaFilterModeLinear;
	Prj_tex.normalized = false;

	cudaBindTextureToArray(Prj_tex, array, channelDesc);

	// TODO: error value?

	return true;
}

bool CFDK::getReconstruction(float* hMergeReconData)
{
	int reconDataLen = m_ReconColumnCount*m_ReconRowCount*m_ReconSliceCount;
	if((this->hReconData == NULL) || (this->hReconDataLen != reconDataLen))
	{
		if(this->hReconData != NULL)
		{
			delete[] this->hReconData;
			this->hReconData = NULL;
			this->hReconDataLen = 0;
		}

		this->hReconData = new float[reconDataLen]; 
		memset(hReconData, 0, sizeof(float)*reconDataLen);

		this->hReconDataLen = reconDataLen;
	}

	cudaMemcpy2D(hReconData, 
		sizeof(float)*m_ReconColumnCount, 
		m_DReconData, 
		sizeof(float)*m_nReconDataPitch, 
		sizeof(float)*m_ReconColumnCount, 
		m_ReconRowCount*m_ReconSliceCount, 
		cudaMemcpyDeviceToHost );

	float* pReconData; 
	float* pMergeData; 

	int dataCountPerSlice =m_ReconColumnCount*m_ReconRowCount;
	for( int iSlice=0; iSlice<m_ReconSliceCount; iSlice=iSlice+m_nMergeNum)
	{
		int iSRIndex = iSlice*dataCountPerSlice;
		int iSMIndex = dataCountPerSlice*iSlice/m_nMergeNum;
		for(int iMatrix=0; iMatrix<dataCountPerSlice; iMatrix++ )
		{			 
			pReconData = hReconData + iSRIndex + iMatrix ;
			pMergeData = hMergeReconData + iSMIndex + iMatrix;

			float all=0;
			for(int iM =0; iM < m_nMergeNum; ++iM)
			{
				all +=  *(pReconData + iM*dataCountPerSlice);
			}
			*pMergeData = all/m_nMergeNum;
		}
	}
			
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
	{
		return false;
	}

	return true;
}

bool CFDK::CallRecon()
{
	bool isOK = true;

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
	{
		isOK = false;
	}

	if(isOK)
	{
		dim3 Block(16,16,1);
		dim3 DetGrid((m_DetectorColumnCount+15)/16,(m_DetectorRowCount*m_nPrjBatchSize+15)/16,1);
		dim3 FFTGrid((m_iHalfFFTSize+15)/16,(m_DetectorRowCount*m_nPrjBatchSize+15)/16,1);
		dim3 BPGrid((m_ReconColumnCount+15)/16,(m_ReconRowCount+15)/16,m_nPrjBatchSize);
	
		cudaMemcpy3DParms cpy3d = {0};
		cpy3d.srcPtr = make_cudaPitchedPtr(m_DsinoData, sizeof(float)*m_iPaddedDetCount, sizeof(float)*m_DetectorColumnCount, m_DetectorRowCount);//
		cpy3d.dstArray = m_PrjArray;
		cpy3d.extent = make_cudaExtent(m_DetectorColumnCount, m_DetectorRowCount, m_nPrjBatchSize);
		cpy3d.kind = cudaMemcpyDeviceToDevice;

		for (int j = 0 ; j < m_nBatch; j++)
		{
			/// <summary>
			/// Copy host Sinogram to GPU Memory
			/// </summary>
			cudaMemset(m_DsinoData, 0, sizeof(cufftReal)*m_DetectorRowCount*m_iPaddedDetCount* m_nPrjBatchSize);
			cudaMemcpy2D(m_DsinoData, 
						sizeof(cufftReal)*m_iPaddedDetCount, 
						m_hSinoData + j*m_DetectorColumnCount*m_DetectorRowCount*m_nPrjBatchSize,  
						sizeof(float)*m_DetectorColumnCount, 
						sizeof(float)*m_DetectorColumnCount,
						m_DetectorRowCount*m_nPrjBatchSize, 
						cudaMemcpyHostToDevice);
	
			err = cudaGetLastError();
			if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
			{
				isOK = false;
				break;
			}
			/// <summary>
			/// Step1. Do FDK Preweighting
			/// </summary>
			Correction<<<DetGrid,Block>>>(m_DsinoData, m_iPaddedDetCount, m_DCorrectionMatrix, m_DetectorColumnCount, m_DetectorRowCount, m_nPrjBatchSize);
	
			err = cudaGetLastError();
			if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
			{
				isOK = false;
				break;
			}
			/// <summary>
			/// Step2. do Filter
			/// </summary>
			cufftExecR2C(m_FwdFFT, m_DsinoData, m_DFilteredsinoData);
			Filter_FFT<<<FFTGrid,Block>>>(m_DFilteredsinoData, m_RampFilter, m_DetectorRowCount*m_nPrjBatchSize, 
				m_iHalfFFTSize);
			cufftExecC2R(m_BwdFFT, m_DFilteredsinoData, m_DsinoData);

			err = cudaGetLastError();
			if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
			{
				isOK = false;
				break;
			}
			/// <summary>
			/// Step3. Bind Filtered Sinogram to Array and do Backprojection
			/// </summary>
			cudaMemcpy3D(&cpy3d);
			bindProjDataTexture(m_PrjArray);
		
			BP<<<BPGrid,Block>>>(m_DReconData, m_nReconDataPitch, j);
		
			cudaUnbindTexture(Prj_tex);

			err = cudaGetLastError();
			if (err != cudaSuccess && err != cudaErrorSetOnActiveProcess)
			{
				isOK = false;
				break;
			}
		}	
	}

	////Frees all GPU resources associated with a cuFFT plan and destroys the internal plan data structure. 
	////This function should be called once a plan is no longer needed, to avoid wasting GPU memory
	//cufftDestroy(m_FwdFFT);
	//cufftDestroy(m_BwdFFT);

	return isOK;
}

}