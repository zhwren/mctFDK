#ifndef _INC_MCTCUDAFDKALGORITHM_H
#define _INC_MCTCUDAFDKALGORITHM_H


#include <cuda_runtime.h>
#include <cufft.h>

#include "mctCudaUtilities.h"
#include "mctGlobals.h"
#include "mctSupport.h"

namespace mct
{

class CFDK
{
public:
	CFDK(void);
	~CFDK(void);

	void FreeObjects();
	
	void InitParams();
	bool SetParams(ScannerGeometry scanGeometry, ProjectionParams prjParams, ReconstructionParams reconParams, float* hSinoData, int iGPUIndex);
	bool CallRecon();
	bool getReconstruction(float* hReconData);

private:
	bool setGPUIndex(int iGPUIndex);

	void SetPrjGeometry(ScannerGeometry scanGeometry, ProjectionParams prjParams, ReconstructionParams reconParams);
	bool CpyToSymbol();
	bool allocateBuffers();
	bool bindProjDataTexture(const cudaArray* array);

	/// <summary>
	/// 计算校正矩阵
	/// </summary>
	bool caculateCorrectMatix();

	/// <summary>
	/// 生产RampFilter
	/// </summary>
	bool genRampFilter();

	/// <summary>
	/// 球管焦点到旋转中心距离，球管焦点到探测器距离
	/// </summary>
	float m_SourceToIsoCenterDistance;
	float m_SourceToDetectorDistance;

	/// <summary>
	/// 探测器行列数目，524*16,由于每两个探测器模块之间插值了一列，因此为504+20 = 524
	/// </summary>
	int m_DetectorColumnCount;
	int m_DetectorRowCount;
	
	float m_DetectorSpacingRampFilter;
	float m_DetectorSpacingMatix;
	/// <summary>
	/// 探测器每个像素尺寸大小
	/// </summary>
	float m_DetectorSpacingX;
	float m_DetectorSpacingY;

	/// <summary>
	/// 探测器总的尺寸大小
	/// </summary>
	float m_DetectorLengthX;
	float m_DetectorLengthY;

	/// <summary>
	/// 探测器中心位置
	/// </summary>
	float m_DetectorColumnRayCenter;
	float m_DetectorRowRayCenter;

	int m_ProjectionAngleCountRampFilter;
	int m_ProjectionAngleCountMatix;
	int m_ProjectionAngleCountAngle;
	/// <summary>
	/// 投影角度数目，间隔，起始位置
	/// </summary>
	int m_ProjectionAngleCount;
	float m_ProjectionAngleStep;
	float m_ProjectionAngleStart;

	/// <summary>
	/// 重建长宽高
	/// </summary>
	int m_ReconColumnCount;
	int m_ReconRowCount;
	int m_ReconSliceCount;

	/// 探测器排数合并
	int m_nMergeNum;

	/// <summary>
	/// 重建中心的长宽高位置
	/// </summary>
	float m_ReconWindowMidColumn;
	float m_ReconWindowMidRow;
	float m_ReconWindowMidSlice;

	/// <summary>
	/// 重建像素大小
	/// </summary>
	float m_PixelSpacingX;
	float m_PixelSpacingY;
	float m_PixelSpacingZ;

	float m_fFOV;

	static const int m_nBatch;
	int m_nPrjBatchSize;

	/// <summary>
	/// 重建滤波核
	/// </summary>
	ReconFilter m_ReconFilter;

	/// <summary>
	/// 进行滤波补零后的长度
	/// </summary>
	int m_iPaddedDetCount; 
	int m_iHalfFFTSize;

	/// <summary>
	/// Host Sinogram Memory Address
	/// </summary>
	float *m_hSinoData;
	
	int m_DReconDataLenSize;
	int m_DReconDataLenCount;
	/// <summary>
	/// GPU Memory Address for Reconstruction Output, 分配成2D线性内存的形式
	/// </summary>
	float *m_DReconData; 
	size_t m_nReconDataPitch;

	/// <summary>
	/// GPU Memory Address for FDK Preweighting Matrix, 分配成1D线性内存的形式
	/// </summary>
	float *m_DCorrectionMatrix;

	int m_PrjArrayLen;
	/// <summary>
	/// GPU Array Address for FDK BackProjection
	/// </summary>
	cudaArray *m_PrjArray;

	int m_DsinoDataLen;
	/// <summary>
	/// GPU Memory Address for ZeroPadded Sinogram
	/// </summary>
	cufftReal *m_DsinoData;
	
	int m_RampFilterLen;
	/// <summary>
	/// GPU Memory Address for FDK RampFilter
	/// </summary>
	cufftComplex *m_RampFilter;

	int m_DFilteredsinoDataLen;
	/// <summary>
	/// GPU Memory Address for Filtered Sinogram
	/// </summary>
	cufftComplex *m_DFilteredsinoData;

	int m_iPaddedDetCountOld;
	int m_FFTLen;
	/// <summary>
	/// cufftHandle for Sinogram FFT, IFFT
	/// </summary>
	cufftHandle m_FwdFFT,m_BwdFFT;

	int hReconDataLen;
	float* hReconData;
};

}
#endif
