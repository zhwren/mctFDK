#include "mctCudaUtilities.h"
#include <assert.h>
#include <iostream>

using namespace std;

namespace mct
{

//
// Get the devices that are available.
//
int CudaGetAvailableDevices(std::vector<cudaDeviceProp> &devices)
{
  int numAvailableDevices = 0;
  cudaGetDeviceCount(&numAvailableDevices);

  if (numAvailableDevices == 0)
    {
    return 0;
    }

  devices.resize(numAvailableDevices);

  for (int i = 0; i < numAvailableDevices; ++i)
    {
    cudaGetDeviceProperties(&devices[i], i);
    }

  return numAvailableDevices;
}

//
// Get the device that has the maximum FLOPS
//
int CudaGetMaxFlopsDev()
{
	std::vector<cudaDeviceProp> devices;

	int numAvailableDevices = CudaGetAvailableDevices(devices);
	if (numAvailableDevices == 0)
	{

    return -1;
    }

  int max_flops = 0;
  int max_flops_device = 0;
  for (int i = 0; i < numAvailableDevices; ++i)
    {
    int flops = devices[i].multiProcessorCount * devices[i].clockRate;
    if (flops > max_flops)
      {
      max_flops = flops;
      max_flops_device = i;
      }
    }

  return max_flops_device;
}


std::vector<int> GetListOfCudaDevices()
{
  std::vector<int>      deviceList;
  int                   deviceCount;
  struct cudaDeviceProp properties;

  if (cudaGetDeviceCount(&deviceCount) == cudaSuccess)
    {
    for (int device = 0; device < deviceCount; ++device) {
      cudaGetDeviceProperties(&properties, device);
      if (properties.major != 9999)   /* 9999 means emulation only */
        deviceList.push_back(device);
      }
    }

  if(deviceList.size()<1)
      std::cout << "No CUDA device available" << std::endl;

  return deviceList;
}


std::pair<int, int> GetCudaComputeCapability(int device)
{
  struct cudaDeviceProp properties;

  if (cudaGetDeviceProperties(&properties, device) != cudaSuccess)
    {
        std::cout << "Unvalid CUDA device" << std::endl;
    }
  return std::make_pair(properties.major, properties.minor);
}

//
// Print device name & info
//
void CudaPrintDeviceInfo(int device)
{
  cudaDeviceProp prop;

  if (cudaGetDeviceProperties(&prop, device) != cudaSuccess)
  {
    std::cout << "Cuda Error : no device found!" << std::endl;
    return;
  }

  std::cout << prop.name << std::endl;
  std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
  std::cout << "Clockrate: " << prop.clockRate << std::endl;
  std::cout << "Global memory: " << prop.totalGlobalMem << std::endl;
  std::cout << "Constant memory: " << prop.totalConstMem << std::endl;
  std::cout << "Number of Multi Processors: " << prop.multiProcessorCount << std::endl;
  std::cout << "Maximum Thread Dim: { " << prop.maxThreadsDim[0] << ", " << prop.maxThreadsDim[1] << ", " << prop.maxThreadsDim[2] << " }" << std::endl;
  std::cout << "Maximum Threads per Block: " << prop.maxThreadsPerBlock << std::endl;
  std::cout << "Maximum Grid Size: { " << prop.maxGridSize[0] << ", " << prop.maxGridSize[1] << ", " << prop.maxGridSize[2] << " }" << std::endl;
}


//
// Find the Cuda platform that matches the "name"
//
int CudaSelectPlatform(const char* name)
{
  int numAvailableDevices = 0;

  std::vector<cudaDeviceProp> devices;
  numAvailableDevices = CudaGetAvailableDevices(devices);

  if (numAvailableDevices == 0)
  {
    std::cout << "Cuda Error : no device found!" << std::endl;
    return -1;
  }

  for (int i = 0; i < numAvailableDevices; ++i)
  {
    if (!strcmp(devices[i].name, name))
    {
      return i;
    }
  }

  return -1;
}

void CudaCheckError(cudaError_t error, const char* filename, int lineno, const char* location)
{
  if (error != cudaSuccess)
  {
    // print error message
    std::ostringstream errorMsg;
    errorMsg << "Cuda Error : " << cudaGetErrorString(error) << std::endl;
    std::cerr << filename << ":" << lineno << " @ " << location << " : " << errorMsg.str() << std::endl;
    //throw e_;
  }
}


/** Check if OpenCL-enabled Cuda is present. */
bool IsCudaAvailable()
{
  int count = 0;
  cudaError_t err = cudaGetDeviceCount(&count);
  CUDA_ASSERT(err);
  return count >= 1;
}


} // end namespace mct