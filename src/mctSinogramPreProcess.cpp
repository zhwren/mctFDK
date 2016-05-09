#include <cmath> 
#include <stdio.h>
#include <memory.h>
using namespace std;

//#include "mctGlobals.h"
#include "mctSinogramPreProcess.h"


namespace mct
{
	SinogramPreProcess::SinogramPreProcess()
	{
		this->m_ProcessedSinogram = NULL;
		this->m_ProcessedSinogramLen = 0;
	}

	SinogramPreProcess::~SinogramPreProcess()
	{
		this->FreeRams();
	}
	
	void SinogramPreProcess::FreeRams()
	{
		if(m_ProcessedSinogram!=NULL)
		{
			delete[] m_ProcessedSinogram;
			m_ProcessedSinogram = NULL;
			this->m_ProcessedSinogramLen = 0;
		}
		if(m_ProcessedSinogram_1!=NULL)
		{
			delete[] m_ProcessedSinogram_1;
			m_ProcessedSinogram_1 = NULL;
		}
	}

	int SinogramPreProcess::SetPreParams(ScannerGeometry scanGeometry, ReconstructionParams reconParams, 
		int prjViews, float *pDarkImg, float *pAirscanImg, float *pSinogram)
	{
		this-> m_DetectorColumns = scanGeometry.m_DetectorColumnCount;
		this-> m_DetectorRows = scanGeometry.m_DetectorRowCount;
		this-> m_DetectorCounts = scanGeometry.m_DetectorCount;
		this-> m_ColumnsPerDetector = scanGeometry.m_ColumnsPerDetector;

		this-> m_ProjectionCounts = prjViews;

		this->m_ReconSliceCount = reconParams.m_ReconSliceCount;


		this->m_DarkImg = pDarkImg;
		this->m_AirScanImg = pAirscanImg;
		this->m_Sinogram = pSinogram;


		return 0;
	}
	
	int SinogramPreProcess::CallPreProcess()
	{
		MCT_ASSERT(m_DarkImg);
		MCT_ASSERT(m_AirScanImg);
		MCT_ASSERT(m_Sinogram);
		
		int processedSinogramLen = (m_DetectorColumns+m_DetectorCounts-1)*m_DetectorRows*m_ProjectionCounts;
		if((this->m_ProcessedSinogram == NULL) || (this->m_ProcessedSinogramLen < processedSinogramLen))
		{
			if(this->m_ProcessedSinogram != NULL)
			{
				delete[] this->m_ProcessedSinogram;
				this->m_ProcessedSinogram = NULL;
			}
			this->m_ProcessedSinogram = new float[processedSinogramLen];
			this->m_ProcessedSinogramLen = processedSinogramLen;

			this->m_ProcessedSinogram_1 = new float[this->m_ProcessedSinogramLen];
		}	
		memset(this->m_ProcessedSinogram, 0, sizeof(float)*this->m_ProcessedSinogramLen);		
		memset(this->m_ProcessedSinogram_1, 0, sizeof(float)*this->m_ProcessedSinogramLen);

		float *pDark = m_DarkImg;
		float *pAir = m_AirScanImg;
		float *pSinoProcessed, *pSinoProcessed_1;	//计算校正系数所用指针 
		float *pSino = m_Sinogram;
		float Sino1,Sino2,Sino3,Sino4,iView1,iView2,iRow1,iRow2,Air_aver,Sino_aver;
		float Views = 3;
		float Rows = 1;

		//转换成电流值
		/*for(int iNum = 0; iNum < m_DetectorRows*m_DetectorColumns; iNum ++)
		{
			*pDark = 0.15492*127*2.048*(*pDark-1024)*0.00002;
			pDark++;
		}
		for(int iPrj = 0; iPrj < m_ProjectionCounts; iPrj ++)
		{
			pDark = m_DarkImg;
			for(int iNum = 0; iNum < m_DetectorRows*m_DetectorColumns; iNum ++)
			{
				*pAir = 0.15492*127*2.048*(*pAir-1024)*0.00002;
				*pSino = 0.15492*127*2.048*(*pSino-1024)*0.00002;
				*pAir = *pAir - *pDark;
				*pSino = *pSino - *pDark;
				pDark++; pAir++; pSino++;
			}
		}*/

		
		//束流一致性校正		
		for(int iPrj = 0; iPrj < m_ProjectionCounts; iPrj ++)
		{
			pDark = m_DarkImg;
			pAir = m_AirScanImg + iPrj*m_DetectorRows*m_DetectorColumns;
			pSino = m_Sinogram + iPrj*m_DetectorRows*m_DetectorColumns;
			Air_aver = 0;
			Sino_aver = 0;
			for(int iRows = 0; iRows < m_DetectorRows; iRows++)
			{
				Air_aver += pAir[iRows*m_DetectorColumns]-pDark[iRows*m_DetectorColumns];
				Air_aver += pAir[iRows*m_DetectorColumns + m_DetectorColumns-1]-pDark[iRows*m_DetectorColumns + m_DetectorColumns-1];
				Sino_aver += pSino[iRows*m_DetectorColumns]-pDark[iRows*m_DetectorColumns];
				Sino_aver += pSino[iRows*m_DetectorColumns + m_DetectorColumns-1]-pDark[iRows*m_DetectorColumns + m_DetectorColumns-1];
			}
			float fator =  Air_aver/Sino_aver;
			if((fator < 0.9) || (fator > 1.1))   //如果超出范围，则认为边缘被遮挡或其它异常，直接使用原始值
			{
				fator = 1;
			}
			for(int Count = 0; Count < m_DetectorRows*m_DetectorColumns; Count++)
			{
				*(pSino+Count) = (*(pSino+Count)-*pDark) *fator;
				pDark++;
			}
		}
		
		pAir = m_AirScanImg;
		pSino = m_Sinogram;
		for(int iPrj = 0; iPrj < m_ProjectionCounts; iPrj ++)
		{		
			pDark = m_DarkImg;
			for(int iRows = 0; iRows < m_DetectorRows; iRows++)
			{
				for(int iCols = 0; iCols < m_DetectorColumns; iCols++)
				{
					pSinoProcessed = m_ProcessedSinogram + iPrj*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) + iRows*(m_DetectorColumns+m_DetectorCounts-1) + iCols + iCols/m_ColumnsPerDetector;
					pSinoProcessed_1 = m_ProcessedSinogram_1 + iPrj*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) + iRows*(m_DetectorColumns+m_DetectorCounts-1) + iCols + iCols/m_ColumnsPerDetector;					
					if(*pSino > 0.01)
					{
						float logP = (*pAir-*pDark)/(*pSino);
						if(logP>1)
						{
							*pSinoProcessed = log(logP);
						}
						else
						{
							*pSinoProcessed = 0;
						}
					}
					else
					{
						*pSinoProcessed = 0;
					}
					*pSinoProcessed_1 = *pSinoProcessed;

					pAir++; pDark++; pSino++;
				}
			}
		}

		// normal_tech修正
		//for(int iPrj = 0; iPrj < m_ProjectionCounts; iPrj ++)
		//{			
		//	for(int iRows = 0; iRows < m_DetectorRows; iRows++)
		//	{
		//		for(int iDets = 0; iDets < m_DetectorCounts-1; iDets++)
		//		{
		//			Sino1 = 0;
		//			Sino2 = 0;
		//			Sino3 = 0;
		//			Sino4 = 0;
		//			iView1 = (iPrj-Views)>0?(iPrj-Views):0;
		//			iView2 = (iPrj+Views)<(m_ProjectionCounts-1)?(iPrj+Views):(m_ProjectionCounts-1);
		//			iRow1 = (iRows-Rows)>0?(iRows-Rows):0;
		//			iRow2 = (iRows+Rows)<(m_DetectorRows-1)?(iRows+Rows):(m_DetectorRows-1);
		//			for(int iView = iView1; iView <= iView2; iView++)
		//			{
		//				for(int iRow = iRow1; iRow <= iRow2; iRow++)
		//				{
		//					pSinoProcessed = m_ProcessedSinogram_1 + iView*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) 
		//						+ iRow*(m_DetectorColumns+m_DetectorCounts-1) + (iDets+1)*(m_ColumnsPerDetector+1)-1;
		//					Sino1 = Sino1 + *(pSinoProcessed-2) + 0.25*(*(pSinoProcessed+2)-*(pSinoProcessed-2));
		//					Sino2 = Sino2 + *(pSinoProcessed-1);
		//					Sino3 = Sino3 + *(pSinoProcessed-2) + 0.75*(*(pSinoProcessed+2)-*(pSinoProcessed-2));
		//					Sino4 = Sino4 + *(pSinoProcessed+1);
		//				}
		//			}
		//			pSinoProcessed = m_ProcessedSinogram + iPrj*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) 
		//				+ iRows*(m_DetectorColumns+m_DetectorCounts-1) + (iDets+1)*(m_ColumnsPerDetector+1)-1;	
		//			if(((Sino1/Sino2)<1.005)&&((Sino1/Sino2)>0.99))
		//			{
		//				*(pSinoProcessed-1) = *(pSinoProcessed-1) * Sino1 / Sino2;
		//			}
		//			if(((Sino3/Sino4)<1.005)&&((Sino3/Sino4)>0.99))
		//			{
		//				*(pSinoProcessed+1) = *(pSinoProcessed+1) * Sino3 / Sino4;
		//			}
		//		}
		//	}
		//}

		for(int iRows = 0; iRows < m_ProjectionCounts*m_DetectorRows; iRows ++)
		{
			//插值
			for(int iDets = 0; iDets < m_DetectorCounts-1; iDets++)
			{
				pSinoProcessed = m_ProcessedSinogram + (m_DetectorColumns+m_DetectorCounts-1)*iRows + (iDets+1)*(m_ColumnsPerDetector+1)-1;
				
				*pSinoProcessed = (*(pSinoProcessed-1) + *(pSinoProcessed+1))/2;				
			}			
		}

		return 0;
	}

	float* SinogramPreProcess::GetPreProcessedSinogram()
	{
		return this->m_ProcessedSinogram;
	}

}// end namespace mct
