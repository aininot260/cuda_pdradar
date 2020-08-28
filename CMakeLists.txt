# author: aininot260

cmake_minimum_required (VERSION 2.8)
project (PD_Radar)

# -----------------------------------------------------------------------------
# ----- BUILD LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

option(GPUACCELERATED "Whether to program will run on GPU. Set to 1 to enable" 1)

if (GPUACCELERATED)
    find_package(CUDA REQUIRED)
    # Stop nvcc sending c compile flags through using -Xcompiler and breaking
    # on compilation of a cpp file receiving -std=c99. In long term should figure 
    # out why CMAKE_C_FLAGS and not CMAKE_CXX_FLAGS are being sent through to a cpp file
    set(CUDA_PROPAGATE_HOST_FLAGS FALSE)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=compute_50 -code=sm_50 -O2")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

add_subdirectory(src)

if (WIN32 AND GPUACCELERATED)
    message("WIN32 AND GPUACCELERATED")
    include_directories(${PROJECT_SOURCE_DIR}/fftw/)
    link_directories(${PROJECT_SOURCE_DIR}/fftw/)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler \" /MD /W0 /MP /Gm- \" ")
    cuda_add_library(PD_Radar SHARED ${PD_RADAR_SRC} ${PD_RADAR_GPU_SRC})
    target_link_libraries(PD_Radar libfftw3-3 libfftw3f-3 libfftw3l-3 cufft cublas)
elseif(WIN32)
    message("WIN32")
    include_directories(${PROJECT_SOURCE_DIR}/fftw/)
    link_directories(${PROJECT_SOURCE_DIR}/fftw/)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /MP /Gm-")
    add_library(PD_Radar SHARED ${PD_RADAR_SRC})
    target_link_libraries(PD_Radar libfftw3-3 libfftw3f-3 libfftw3l-3)
elseif(GPUACCELERATED)
    message("Linux AND GPUACCELERATED")
    set(CUDA_CUFFT_LIBRARIES "/usr/local/cuda/lib64/libcufft.so")
    set(CUDA_CUBLAS_LIBRARIES "/usr/local/cuda/lib64/libcublas.so")
    set(CUDA_cublas_device_LIBRARY "")

    cuda_add_library(PD_Radar SHARED ${PD_RADAR_SRC} ${PD_RADAR_GPU_SRC})
    target_link_libraries(PD_Radar fftw3 fftw3f)
    CUDA_ADD_CUFFT_TO_TARGET(PD_Radar)
    CUDA_ADD_CUBLAS_TO_TARGET(PD_Radar)
else()
    message("Linux")
    add_library(PD_Radar SHARED ${PD_RADAR_SRC})
    target_link_libraries(PD_Radar fftw3 fftw3f)
endif()

if (GPUACCELERATED)
    target_compile_definitions(PD_Radar PUBLIC USE_GPU)
endif()

target_include_directories(PD_Radar 
    PUBLIC include
)

# -----------------------------------------------------------------------------
# ----- USER EXECUTABLE -------------------------------------------------------
# -----------------------------------------------------------------------------

add_executable(main main.cpp)
target_link_libraries(main PD_Radar)