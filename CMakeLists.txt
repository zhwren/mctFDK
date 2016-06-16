#************************************************************
#* Author       : ZhuHaiWen                                 *
#* Email        : zhuhw@ihep.ac.cn/zhwren0211@whu.edu.cn    *
#* Last modified: 2015-10-27 16:41:1445935314               *
#* Filename     : CMakeLists.txt                            *
#* Phone Number : 18625272373                               *
#* Discription  :                                           *
#************************************************************
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
PROJECT(Mct.exe)
FIND_PACKAGE(CUDA REQUIRED)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include
		    ${CUDA_INCLUDE_DIR})
FILE(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp 
		  ${PROJECT_SOURCE_DIR}/src/*.cu
		  ${PROJECT_SOURCE_DIR}/src/*.cc)
FILE(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h
                  ${PROJECT_SOURCE_DIR}/include/*.hh)
CUDA_ADD_EXECUTABLE(Mct.exe main.cpp ${sources} ${headers})
SET(CUDA_LIB /usr/local/cuda-7.5/lib64/libcufft.so)
TARGET_LINK_LIBRARIES(Mct.exe ${CUDA_LIB})
SET(EXECUTABLE_OUTPUT_PATH /home/zhwren/Workfs/geant4_workdir/bin/Linux-g++/)
