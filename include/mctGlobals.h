#ifndef _INC_GLOBALS
#define _INC_GLOBALS

#include <assert.h>

//----------------------------------------------------------------------------------------
// DEFINE
#ifndef NULL
#define NULL 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FLASE
#define FLASE 0
#endif

//----------------------------------------------------------------------------------------
#define MCT_ASSERT(a) assert(a)

#define MCT_DELETE(a) if (a) { delete a; a = NULL; }
#define MCT_DELETE_ARRAY(a) if (a) { delete[] a; a = NULL; }


//----------------------------------------------------------------------------------------
// errors
namespace mct 
{
typedef enum {MCT_SUCCESS, 
			  MCT_ERROR_NOT_INITIALIZED,
			  MCT_ERROR_INVALID_FILE,
			  MCT_ERROR_OUT_OF_RANGE,
			  MCT_ERROR_DIMENSION_MISMATCH,
			  MCT_ERROR_EXTERNAL_LIBRARY,
			  MCT_ERROR_ALLOCATION,
			  MCT_ERROR_NOT_IMPLEMENTED} MctError;
}

//----------------------------------------------------------------------------------------
// typedefs
namespace mct 
{
	typedef float float32;
	typedef double float64;
	typedef unsigned short int uint16;
	typedef signed short int sint16;
	typedef unsigned char uchar8;
	typedef signed char schar8;

	typedef int int32;
	typedef short int int16;
}

//----------------------------------------------------------------------------------------
// variables
namespace mct 
{
	const float32 PI = 3.14159265358979323846264338328f;
	const float32 PI32 = 3.14159265358979323846264338328f;
	const float32 PIdiv2 = PI / 2;
	const float32 PIdiv4 = PI / 4;
	const float32 eps = 1e-7f;

	#define TWOPI   6.28318530717958647692
	#define SQRT2   1.414213562373095049

	#define F_EPSILON       1.0E-6
	#define D_EPSILON       1.0E-10

	#define ASSUMEDZERO		1E-10
	#define eu				-1.f
}

//----------------------------------------------------------------------------------------
// structures
namespace mct 
{
	/// <summary>
	///  struct ScannerGeometry
	///  CT扫描仪几何参数
	/// </summary>
	struct ScannerGeometry
	{
		/// <summary>
		/// 球管焦点到旋转中心距离，球管焦点到探测器距离
		/// </summary>
		float m_SourceToIsoCenterDistance;
		float m_SourceToDetectorDistance;

		/// <summary>
		/// 探测器行列数目，504*16
		/// </summary>
		int m_DetectorColumnCount;
		int m_DetectorRowCount;

		/// <summary>
		/// 探测器每个像素尺寸大小
		/// </summary>
		float m_DetectorSpacingX;
		float m_DetectorSpacingY;

		/// <summary>
		/// 每个探测器模块的行数24，探测器模块数目21
		/// </summary>
		int m_ColumnsPerDetector;
		int m_DetectorCount;

		/// <summary>
		/// 探测器中心位置
		/// </summary>
		float m_DetectorColumnRayCenter;//524/2=262,由于球管位置偏差，出于矫正作用，将其改为261
		float m_DetectorRowRayCenter;
	};

	/// <summary>
	///  struct ProjectionParams
	///  投影参数
	/// </summary>
	struct ProjectionParams
	{
		/// <summary>
		/// 投影角度数目
		/// </summary>
		int m_ProjectionAngleCount;

		/// <summary>
		/// 投影角度起始位置
		/// 定义投影起始的焦点的Z坐标为0
		/// </summary>
		float m_ProjectionAngleStart;

		/// <summary>
		/// 每一圈的投影角度数目
		/// </summary>
		int m_ProjectionAnglesPerRotation;

		/// <summary>
		/// 旋转一圈Z方向移动距离,人体平面是x-y平面
		/// </summary>
		float m_Pitch;		
		/// <summary>
		/// 圈个数
		/// </summary>
		int RotNum;
	};

	/// <summary>
	///  ReconFilter 
	///  重建滤波核。
	///  
	/// </summary>
	typedef enum 
	{

	}ReconFilter;


	/// <summary>
	///  struct ReconstructionParams
	///  重建参数设置
	/// </summary>
	struct ReconstructionParams
	{
		/// <summary>
		/// 重建长宽高
		/// </summary>
		int m_ReconColumnCount;
		int m_ReconRowCount;
		int m_ReconSliceCount;				

		/// 探测器排数合并
		int m_MergedNum;
		int m_ReconSliceSpaceNum;//指定重建的层间距,单位层
		//float m_SliceSpace;

		/// <summary>
		/// 重建像素大小
		/// </summary>
		float m_PixelSpacingX;
		float m_PixelSpacingY;
		float m_PixelSpacingZ;

		/// <summary>
		/// 重建中心的位置
		/// </summary>
		float m_ReconWindowMidColumn;
		float m_ReconWindowMidRow;
		float m_ReconWindowMidSlice;

		/// <summary>
		/// 重建滤波核
		/// </summary>
		ReconFilter m_ReconFilter;
		float ZMin;//光源起始轴向位置（使用时不用赋值，动态库中自行计算了）
	};
}

#endif