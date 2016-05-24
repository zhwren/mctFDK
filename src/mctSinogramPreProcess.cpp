#include <cmath> 
#include <stdio.h>
#include <memory.h>
using namespace std;

//#include "mctGlobals.h"
#include "mctSinogramPreProcess.h"
#include <iostream>

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
	
	int SinogramPreProcess::CallPreProcess(float* correction)
	{
		MCT_ASSERT(m_DarkImg);
		MCT_ASSERT(m_AirScanImg);
		MCT_ASSERT(m_Sinogram);
		m_Correction = correction;
		
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
		float Air_aver,Sino_aver;
		float iMax;
		int   nColumn;
		float a_0,a_1,a_2,a_3,a_4,Len;
		a_0 = 0.0225; a_1 = -0.0503; a_2 = 0.4400; a_3 = -0.6497; a_4 = 44.6151;

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
			pDark = m_DarkImg; nColumn = 0; iMax = 0;
			for(int iRows = 0; iRows < m_DetectorRows; iRows++)
			{
			        for(int iCols = 0; iCols < m_DetectorColumns; iCols++)
				{
				        if( *pSino < 0.01 ) *pSino = 0.01;
				        *pSino = (*pAir - *pDark) / (*pSino);
					if( *pSino > iMax )
					{
					      nColumn = iCols;
					      iMax    = *pSino;
					}
					*pSino = 1./(*pSino);
					pDark++; pAir++; pSino++;
				}
			}
                        int lCentralCol = nColumn / 24 * 24 + 11;
                        int rCentralCol = lCentralCol + 1;
			int key;
			pSino = pSino - m_DetectorColumns*m_DetectorRows;

			for(int iRows = 0; iRows < m_DetectorRows; iRows++)
			{
				for(int iCols = 0; iCols < m_DetectorColumns; iCols++)
				{
				        std::cout << *pSino << "     ";
				        if( iCols<rCentralCol ) key = lCentralCol - iCols;
					if( iCols>lCentralCol ) key = iCols - rCentralCol;
					if( key==11  ) *pSino *= m_Correction[0];
					if( key==12  ) *pSino *= m_Correction[1];
					if( key<=12  ) *pSino *= m_Correction[2];
					if( key==35  ) *pSino *= m_Correction[3];
					if( key==36  ) *pSino *= m_Correction[4];
					if( key<=36  ) *pSino *= m_Correction[5];
					if( key==59  ) *pSino *= m_Correction[6];
					if( key==60  ) *pSino *= m_Correction[7];
					if( key<=60  ) *pSino *= m_Correction[8];
					if( key==83  ) *pSino *= m_Correction[9];
					if( key==84  ) *pSino *= m_Correction[10];
					if( key<=84  ) *pSino *= m_Correction[11];
					if( key==107 ) *pSino *= m_Correction[12];
					if( key==108 ) *pSino *= m_Correction[13];
					if( key<=108 ) *pSino *= m_Correction[14];
					if( key==131 ) *pSino *= m_Correction[15];
					if( key==132 ) *pSino *= m_Correction[16];
					if( key<=132 ) *pSino *= m_Correction[17];
					if( key==155 ) *pSino *= m_Correction[18];
					if( key==156 ) *pSino *= m_Correction[19];
					if( key>156  ) *pSino *= m_Correction[20];
				        std::cout << *pSino << std::endl;;

					pSinoProcessed = m_ProcessedSinogram + iPrj*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) + iRows*(m_DetectorColumns+m_DetectorCounts-1) + iCols + iCols/m_ColumnsPerDetector;
					pSinoProcessed_1 = m_ProcessedSinogram_1 + iPrj*m_DetectorRows*(m_DetectorColumns+m_DetectorCounts-1) + iRows*(m_DetectorColumns+m_DetectorCounts-1) + iCols + iCols/m_ColumnsPerDetector;					
					float logP = 1./(*pSino);
					if(logP>1)
					{
						*pSinoProcessed = log(logP);
					}
					else
					{
						*pSinoProcessed = 0;
					}
					Len = a_1*pow(*pSinoProcessed,4) + a_2*pow(*pSinoProcessed,3) + a_3*pow(*pSinoProcessed,2) + a_4*pow(*pSinoProcessed,1);
					*pSinoProcessed = a_0*Len;
					*pSinoProcessed_1 = *pSinoProcessed;

					pSino++;
				}
			}
		}

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
