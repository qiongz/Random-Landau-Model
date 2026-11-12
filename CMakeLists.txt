cmake_minimum_required(VERSION 3.16)
project(rlm_fp64 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MKLROOT "/opt/intel/oneapi/mkl/2025.2" CACHE PATH "oneMKL root")
set(RLM_GPU "OFF" CACHE STRING "OFF, CUDA, or HIP")
set(RLM_HIP_ARCH "gfx1201" CACHE STRING "HIP offload arch")
set_property(CACHE RLM_GPU PROPERTY STRINGS OFF CUDA HIP)

add_executable(rlm
    src/dispatch.cpp
    src/solver_mkl.cpp
    src/chern.cpp
    src/pipeline.cpp
    src/continuum.cpp
    src/lattice.cpp
    src/cli.cpp
    src/check.cpp
)

target_include_directories(rlm PRIVATE include ${MKLROOT}/include)
target_compile_definitions(rlm PRIVATE RLM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_compile_options(rlm PRIVATE -O3 -fopenmp -march=native)
target_link_directories(rlm PRIVATE ${MKLROOT}/lib)
target_link_options(rlm PRIVATE -fopenmp -Wl,-rpath,${MKLROOT}/lib -Wl,--no-as-needed)
target_link_libraries(rlm PRIVATE mkl_intel_lp64 mkl_gnu_thread mkl_core gomp pthread m dl)

if(RLM_GPU STREQUAL "CUDA")
    enable_language(CUDA)
    target_sources(rlm PRIVATE src/solver_gpu.cu)
    set_source_files_properties(src/solver_gpu.cu PROPERTIES LANGUAGE CUDA)
    target_compile_definitions(rlm PRIVATE RLM_HAS_GPU=1 RLM_CUDA=1)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(rlm PRIVATE CUDA::cusolver CUDA::cublas CUDA::cudart)
elseif(RLM_GPU STREQUAL "HIP")
    set(RLM_GPU_OBJ ${CMAKE_BINARY_DIR}/solver_gpu.o)
    add_custom_command(
        OUTPUT ${RLM_GPU_OBJ}
        COMMAND hipcc -O3 -std=c++17 -fPIE -DRLM_HAS_GPU=1 -DRLM_HIP=1
                -I${CMAKE_SOURCE_DIR}/include
                --offload-arch=${RLM_HIP_ARCH}
                -c ${CMAKE_SOURCE_DIR}/src/solver_gpu.cu
                -o ${RLM_GPU_OBJ}
        DEPENDS ${CMAKE_SOURCE_DIR}/src/solver_gpu.cu ${CMAKE_SOURCE_DIR}/include/rlm/rlm.hpp
        VERBATIM
    )
    add_custom_target(rlm_gpu_obj DEPENDS ${RLM_GPU_OBJ})
    add_dependencies(rlm rlm_gpu_obj)
    target_compile_definitions(rlm PRIVATE RLM_HAS_GPU=1 RLM_HIP=1)
    target_link_directories(rlm PRIVATE /opt/rocm/lib)
    target_link_options(rlm PRIVATE -Wl,-rpath,/opt/rocm/lib)
    target_link_libraries(rlm PRIVATE ${RLM_GPU_OBJ} rocsolver rocblas amdhip64)
elseif(NOT RLM_GPU STREQUAL "OFF")
    message(FATAL_ERROR "RLM_GPU must be OFF, CUDA, or HIP (got ${RLM_GPU})")
endif()
