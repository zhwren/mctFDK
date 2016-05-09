#ifndef __mctCudaUtilities_hcu
#define __mctCudaUtilities_hcu

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cuda_runtime.h>


namespace mct
{

#define CUDA_CHECK_ERROR \
    do{ \
    cudaError_t err = cudaGetLastError(); \
    if (cudaSuccess != err) \
	{                                               \
	fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
	err, __FILE__, __LINE__, cudaGetErrorString( err) );\
	exit(EXIT_FAILURE);} \
    }while(0)

#define CUFFT_CHECK_ERROR(result) \
    { \
    if (result) \
	fprinf(stderr, "CUFFT error: %s in file '%s' in line %i.\n",    \
	result, __FILE__, __LINE__ );    \
    }


#define CUDA_ASSERT(errorMessage) do {                                     \
	cudaError_t err = cudaThreadSynchronize();                               \
	if( cudaSuccess != err) {                                                \
	fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
	errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
	exit(EXIT_FAILURE);                                                  \
	} } while (0)

std::vector<int> GetListOfCudaDevices();

std::pair<int, int> GetCudaComputeCapability(int device);

/** Get the devices that are available */
int CudaGetAvailableDevices(std::vector<cudaDeviceProp> &devices);

/** Get the device that has the maximum FLOPS in the current context */
int CudaGetMaxFlopsDev();

/** Print device name and info */
void CudaPrintDeviceInfo(int device, bool verbose = false);

/** Find the Cuda platform that matches the "name" */
int CudaSelectPlatform(const char* name);

int CudaGetAvailableDevices(std::vector<cudaDeviceProp> &devices);

/** Check Cuda error */
void CudaCheckError(cudaError_t error, const char* filename = "", int lineno = 0, const char* location = "");

/** Check if Cuda-enabled Cuda is present. */
bool IsCudaAvailable();

/** Get Typename */
std::string GetTypename(const std::type_info& intype);

/** Get Typename in String if a valid type */
bool GetValidTypename(const std::type_info& intype, const std::vector<std::string>& validtypes, std::string& retTypeName);

/** Get 64-bit pragma */
std::string Get64BitPragma();

/** Get Typename in String */
void GetTypenameInString(const std::type_info& intype, std::ostringstream& ret);

/** Get pixel dimension (number of channels).
 * For high-dimensional pixel format, only itk::Vector< type, 2/3 > is acceptable. */
int GetPixelDimension(const std::type_info& intype);

}

#endif
