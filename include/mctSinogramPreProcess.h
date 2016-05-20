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
		/// 原始数据预处理，包括减暗场，log，每两个探测器模块需要插值一列；
		/// 可以把这些工作用GPU做，但后期有些校正处理部分需放在重建前做，暂时用CPU进行预处理
		/// </summary>
		int CallPreProcess(std::map<int,float>);

		float* GetPreProcessedSinogram();

	private:
		/// <summary>
		/// 探测器行列数目504*16
		/// </summary>
		int m_DetectorColumns;
		int m_DetectorRows;

		/// <summary>
		/// 探测器模块数目
		/// </summary>
		int m_DetectorCounts;

		/// <summary>
		/// 每个探测器模块行数目
		/// </summary>
		int m_ColumnsPerDetector;

		/// <summary>
		/// 投影角度数目
		/// </summary>
		int m_ProjectionCounts;
		//每圈重建的影像个数
		int m_ReconSliceCount;

		/// <summary>
		/// 暗场、空拍、投影数据
		/// </summary>
		float *m_DarkImg;
		float *m_AirScanImg;
		float *m_Sinogram;
				
		/// <summary>
		/// 处理后的投影数据
		/// </summary>
		float *m_ProcessedSinogram;
		int m_ProcessedSinogramLen;
		float *m_ProcessedSinogram_1;

		std::map<int,float> m_Correction;
	};
}

#endif
