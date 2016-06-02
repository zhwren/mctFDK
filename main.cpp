#include "mct.h"
#include <stdio.h>

using namespace mct;

/// <summary>
/// Init paras to fit compute condition
/// </summary>
/// <returns></returns>
ScannerGeometry InitScannerGeometry()
{
	ScannerGeometry scannerGeometry;

	scannerGeometry.m_SourceToIsoCenterDistance = 240.0f;
	scannerGeometry.m_SourceToDetectorDistance = 440.0f;


	scannerGeometry.m_DetectorColumnCount = 504;
	scannerGeometry.m_DetectorRowCount = 16;

	scannerGeometry.m_DetectorSpacingX = 0.002443;//0.14167*PI/180;;
	scannerGeometry.m_DetectorSpacingY = 2.005f;//1.915f;

	scannerGeometry.m_ColumnsPerDetector = 24;
	scannerGeometry.m_DetectorCount = 21;

	scannerGeometry.m_DetectorColumnRayCenter = 262;//264;
	scannerGeometry.m_DetectorRowRayCenter = 7.5f;

	return scannerGeometry;
}

ProjectionParams InitProjectionParams()
{
	ProjectionParams prjParams;

	prjParams.m_ProjectionAngleCount = 1250;//1000;

	prjParams.m_ProjectionAngleStart = 0.0f;

	prjParams.m_ProjectionAnglesPerRotation = 1250;//1000;

	prjParams.RotNum = 1;

	prjParams.m_Pitch = 10;//0.0f;

	return prjParams;
}

ReconstructionParams InitReconParams()
{
	ReconstructionParams reconParams;

	reconParams.m_ReconColumnCount = 512;
	reconParams.m_ReconRowCount = 512;
	reconParams.m_ReconSliceCount = 16;

	reconParams.m_PixelSpacingX = 0.5f;
	reconParams.m_PixelSpacingY = 0.5f;
	reconParams.m_PixelSpacingZ = 0.5f;


	reconParams.m_ReconWindowMidColumn = (reconParams.m_ReconColumnCount-1)/2.0;
	reconParams.m_ReconWindowMidRow = (reconParams.m_ReconRowCount-1)/2.0;
	reconParams.m_ReconWindowMidSlice = (reconParams.m_ReconSliceCount-1)/2.0;

	reconParams.m_MergedNum=1;
	reconParams.m_ReconSliceSpaceNum = 0;
	//reconParams.m_ReconFilter = 0 ;

	return reconParams;
}


int main(int argc, char* argv[])
{
        clock_t start = clock();
	ScannerGeometry scannerGeometry = InitScannerGeometry();
	ProjectionParams prjParams = InitProjectionParams();
	ReconstructionParams reconParams = InitReconParams();


	float *pDarkImg = new float[scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount];
	float *pAirscanImg = new float[scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount*prjParams.m_ProjectionAngleCount];
	float *pSinogram = new float[scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount*prjParams.m_ProjectionAngleCount];
	float *pRecon = new float[reconParams.m_ReconColumnCount*reconParams.m_ReconRowCount*reconParams.m_ReconSliceCount];

	FILE *fp = fopen("/home/zhwren/Workfs/CTRecons/Datas/20160601/Dark/Dark.bin","rb");
	fread(pDarkImg,sizeof(float),scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount,fp);
	fclose(fp);

	fp = fopen("/home/zhwren/Workfs/CTRecons/Datas/20160601/Air/AirScan.bin","rb");
	fread(pAirscanImg,sizeof(float),scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount*prjParams.m_ProjectionAngleCount,fp);
	fclose(fp);

	fp = fopen("/home/zhwren/Workfs/CTRecons/Datas/20160601/Density/ScanData.bin","rb");
	fread(pSinogram,sizeof(float),scannerGeometry.m_DetectorColumnCount*scannerGeometry.m_DetectorRowCount*prjParams.m_ProjectionAngleCount,fp);
	fclose(fp);


	float *correction = new float[21];
	for(int i=0; i<21; i++) correction[i] = 1.;

	std::string Correction = argc>1?argv[1]:"";
	std::ifstream file( Correction.c_str() );
	if( file.is_open() )
	{
	  int key=0;
	  while( file >> correction[key] )
	    key++;
	}
	file.close();

	SinogramPreProcess process;
	process.SetPreParams(scannerGeometry, reconParams, prjParams.m_ProjectionAngleCount, pDarkImg, pAirscanImg, pSinogram);;
	process.CallPreProcess(correction);

	CFDK recon;
	recon.SetParams(scannerGeometry, prjParams, reconParams, process.GetPreProcessedSinogram(), 0);
	//CFDK recon(scannerGeometry, prjParams, reconParams, pSinogram, 0);
	recon.CallRecon();
	recon.getReconstruction(pRecon);

	std::string ReconData = argc>2?argv[2]:"ReconData.rcn";
	fp = fopen( ReconData.c_str(), "wb");
	fwrite(pRecon,sizeof(float),reconParams.m_ReconColumnCount*reconParams.m_ReconRowCount*reconParams.m_ReconSliceCount,fp);
	fclose(fp);

	delete[] pDarkImg;
	delete[] pAirscanImg;
	delete[] pSinogram;
	delete[] pRecon;
	clock_t finish = clock();
	std::cerr << double(finish-start)/CLOCKS_PER_SEC << std::endl;
}
