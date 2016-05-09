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
	/// ����У������
	/// </summary>
	bool caculateCorrectMatix();

	/// <summary>
	/// ����RampFilter
	/// </summary>
	bool genRampFilter();

	/// <summary>
	/// ��ܽ��㵽��ת���ľ��룬��ܽ��㵽̽��������
	/// </summary>
	float m_SourceToIsoCenterDistance;
	float m_SourceToDetectorDistance;

	/// <summary>
	/// ̽����������Ŀ��524*16,����ÿ����̽����ģ��֮���ֵ��һ�У����Ϊ504+20 = 524
	/// </summary>
	int m_DetectorColumnCount;
	int m_DetectorRowCount;
	
	float m_DetectorSpacingRampFilter;
	float m_DetectorSpacingMatix;
	/// <summary>
	/// ̽����ÿ�����سߴ��С
	/// </summary>
	float m_DetectorSpacingX;
	float m_DetectorSpacingY;

	/// <summary>
	/// ̽�����ܵĳߴ��С
	/// </summary>
	float m_DetectorLengthX;
	float m_DetectorLengthY;

	/// <summary>
	/// ̽��������λ��
	/// </summary>
	float m_DetectorColumnRayCenter;
	float m_DetectorRowRayCenter;

	int m_ProjectionAngleCountRampFilter;
	int m_ProjectionAngleCountMatix;
	int m_ProjectionAngleCountAngle;
	/// <summary>
	/// ͶӰ�Ƕ���Ŀ���������ʼλ��
	/// </summary>
	int m_ProjectionAngleCount;
	float m_ProjectionAngleStep;
	float m_ProjectionAngleStart;

	/// <summary>
	/// �ؽ������
	/// </summary>
	int m_ReconColumnCount;
	int m_ReconRowCount;
	int m_ReconSliceCount;

	/// ̽���������ϲ�
	int m_nMergeNum;

	/// <summary>
	/// �ؽ����ĵĳ����λ��
	/// </summary>
	float m_ReconWindowMidColumn;
	float m_ReconWindowMidRow;
	float m_ReconWindowMidSlice;

	/// <summary>
	/// �ؽ����ش�С
	/// </summary>
	float m_PixelSpacingX;
	float m_PixelSpacingY;
	float m_PixelSpacingZ;

	float m_fFOV;

	static const int m_nBatch;
	int m_nPrjBatchSize;

	/// <summary>
	/// �ؽ��˲���
	/// </summary>
	ReconFilter m_ReconFilter;

	/// <summary>
	/// �����˲������ĳ���
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
	/// GPU Memory Address for Reconstruction Output, �����2D�����ڴ����ʽ
	/// </summary>
	float *m_DReconData; 
	size_t m_nReconDataPitch;

	/// <summary>
	/// GPU Memory Address for FDK Preweighting Matrix, �����1D�����ڴ����ʽ
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
