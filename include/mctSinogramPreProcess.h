#ifndef _INC_MCT_RAWCORRECTION
#define _INC_MCT_RAWCORRECTION

#include "mctGlobals.h"
#include <map>

namespace mct
{
	class SinogramPreProcess
	{
	public:
		SinogramPreProcess();
		~SinogramPreProcess();

		void FreeRams();

		int SetPreParams(ScannerGeometry scanGeometry, ReconstructionParams reconParams, int prjViews, float *pDarkImg, float *pAirscanImg, float *pSinogram);
		/// <summary>
		/// ԭʼ����Ԥ����������������log��ÿ����̽����ģ����Ҫ��ֵһ�У�
		/// ���԰���Щ������GPU������������ЩУ��������������ؽ�ǰ������ʱ��CPU����Ԥ����
		/// </summary>
		int CallPreProcess(std::map<int,float>);

		float* GetPreProcessedSinogram();

	private:
		/// <summary>
		/// ̽����������Ŀ504*16
		/// </summary>
		int m_DetectorColumns;
		int m_DetectorRows;

		/// <summary>
		/// ̽����ģ����Ŀ
		/// </summary>
		int m_DetectorCounts;

		/// <summary>
		/// ÿ��̽����ģ������Ŀ
		/// </summary>
		int m_ColumnsPerDetector;

		/// <summary>
		/// ͶӰ�Ƕ���Ŀ
		/// </summary>
		int m_ProjectionCounts;
		//ÿȦ�ؽ���Ӱ�����
		int m_ReconSliceCount;

		/// <summary>
		/// ���������ġ�ͶӰ����
		/// </summary>
		float *m_DarkImg;
		float *m_AirScanImg;
		float *m_Sinogram;
				
		/// <summary>
		/// ������ͶӰ����
		/// </summary>
		float *m_ProcessedSinogram;
		int m_ProcessedSinogramLen;
		float *m_ProcessedSinogram_1;

		std::map<int,float> m_Correction;
	};
}

#endif
