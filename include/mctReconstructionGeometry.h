#ifndef _INC_MCT_RECONSTRUCTIONGEOMETRY
#define _INC_MCT_RECONSTRUCTIONGEOMETRY

namespace mct
{
	class ReconstructionGeometry
	{
	public:
		ReconstructionGeometry(int reconRows, int reconCols, int reconSlices, float pixelSpacingX, float pixelSpacingY, float pixelSpacingZ, float ReconMidX, float ReconMidY, float ReconMidZ);
		ReconstructionGeometry(const ReconstructionGeometry& reconGeometry);
		void operator=(const ReconstructionGeometry& reconGeometry);
		~ReconstructionGeometry();

		bool isInitialized() const;

		float GetPixelSpacingX() const;
		float GetPixelSpacingY() const;
		float GetPixelSpacingZ() const;

		float GetReconWindowMinX() const;
		float GetReconWindowMaxX() const;
		float GetReconWindowMinY() const;
		float GetReconWindowMaxY() const;
		float GetReconWindowMinZ() const;
		float GetReconWindowMaxZ() const;

		float GetReconWindowLengthX() const;
		float GetReconWindowLengthY() const;
		float GetReconWindowLengthZ() const;

		int GetReconRowCount() const;
		int GetReconColumnCount() const;
		int GetReconSliceCount() const;

	private:	
		bool m_bInitialized;

		float m_PixelSpacingX;
		float m_PixelSpacingY;
		float m_PixelSpacingZ;

		float m_ReconWindowMinX;
		float m_ReconWindowMaxX;
		float m_ReconWindowMinY;
		float m_ReconWindowMaxY;
		float m_ReconWindowMinZ;
		float m_ReconWindowMaxZ;

		float m_ReconWindowLengthX;
		float m_ReconWindowLengthY;
		float m_ReconWindowLengthZ;

		int m_ReconRowCount;
		int m_ReconColumnCount;
		int m_ReconSliceCount;
	};
}


#endif