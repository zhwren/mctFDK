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
	///  CTɨ���Ǽ��β���
	/// </summary>
	struct ScannerGeometry
	{
		/// <summary>
		/// ��ܽ��㵽��ת���ľ��룬��ܽ��㵽̽��������
		/// </summary>
		float m_SourceToIsoCenterDistance;
		float m_SourceToDetectorDistance;

		/// <summary>
		/// ̽����������Ŀ��504*16
		/// </summary>
		int m_DetectorColumnCount;
		int m_DetectorRowCount;

		/// <summary>
		/// ̽����ÿ�����سߴ��С
		/// </summary>
		float m_DetectorSpacingX;
		float m_DetectorSpacingY;

		/// <summary>
		/// ÿ��̽����ģ�������24��̽����ģ����Ŀ21
		/// </summary>
		int m_ColumnsPerDetector;
		int m_DetectorCount;

		/// <summary>
		/// ̽��������λ��
		/// </summary>
		float m_DetectorColumnRayCenter;//524/2=262,�������λ��ƫ����ڽ������ã������Ϊ261
		float m_DetectorRowRayCenter;
	};

	/// <summary>
	///  struct ProjectionParams
	///  ͶӰ����
	/// </summary>
	struct ProjectionParams
	{
		/// <summary>
		/// ͶӰ�Ƕ���Ŀ
		/// </summary>
		int m_ProjectionAngleCount;

		/// <summary>
		/// ͶӰ�Ƕ���ʼλ��
		/// ����ͶӰ��ʼ�Ľ����Z����Ϊ0
		/// </summary>
		float m_ProjectionAngleStart;

		/// <summary>
		/// ÿһȦ��ͶӰ�Ƕ���Ŀ
		/// </summary>
		int m_ProjectionAnglesPerRotation;

		/// <summary>
		/// ��תһȦZ�����ƶ�����,����ƽ����x-yƽ��
		/// </summary>
		float m_Pitch;		
		/// <summary>
		/// Ȧ����
		/// </summary>
		int RotNum;
	};

	/// <summary>
	///  ReconFilter 
	///  �ؽ��˲��ˡ�
	///  
	/// </summary>
	typedef enum 
	{

	}ReconFilter;


	/// <summary>
	///  struct ReconstructionParams
	///  �ؽ���������
	/// </summary>
	struct ReconstructionParams
	{
		/// <summary>
		/// �ؽ������
		/// </summary>
		int m_ReconColumnCount;
		int m_ReconRowCount;
		int m_ReconSliceCount;				

		/// ̽���������ϲ�
		int m_MergedNum;
		int m_ReconSliceSpaceNum;//ָ���ؽ��Ĳ���,��λ��
		//float m_SliceSpace;

		/// <summary>
		/// �ؽ����ش�С
		/// </summary>
		float m_PixelSpacingX;
		float m_PixelSpacingY;
		float m_PixelSpacingZ;

		/// <summary>
		/// �ؽ����ĵ�λ��
		/// </summary>
		float m_ReconWindowMidColumn;
		float m_ReconWindowMidRow;
		float m_ReconWindowMidSlice;

		/// <summary>
		/// �ؽ��˲���
		/// </summary>
		ReconFilter m_ReconFilter;
		float ZMin;//��Դ��ʼ����λ�ã�ʹ��ʱ���ø�ֵ����̬�������м����ˣ�
	};
}

#endif