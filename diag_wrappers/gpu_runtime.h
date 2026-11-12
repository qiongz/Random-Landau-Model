#pragma once

#ifdef USE_ROCM
#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

typedef hipFloatComplex GpuComplex;

#define GPU_CX_MUL(a, b) hipCmulf((a), (b))
#define GPU_CX_ADD(a, b) hipCaddf((a), (b))
#define GPU_CX_CONJ(a) hipConjf(a)

#define gpuMalloc hipMalloc
#define gpuFree hipFree
#define gpuMemcpy hipMemcpy
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuHostMalloc hipHostMalloc
#define gpuFreeHost hipHostFree
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuHostMallocDefault hipHostMallocDefault

#define gpuEvent_t hipEvent_t
#define gpuEventCreate hipEventCreate
#define gpuEventRecord hipEventRecord
#define gpuEventSynchronize hipEventSynchronize
#define gpuEventElapsedTime hipEventElapsedTime
#define gpuEventDestroy hipEventDestroy

#else
#include <cuda_runtime.h>

typedef cuFloatComplex GpuComplex;

#define GPU_CX_MUL(a, b) cuCmulf((a), (b))
#define GPU_CX_ADD(a, b) cuCaddf((a), (b))
#define GPU_CX_CONJ(a) cuConjf(a)

#define gpuMalloc cudaMalloc
#define gpuFree cudaFree
#define gpuMemcpy cudaMemcpy
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuHostMalloc cudaHostAlloc
#define gpuFreeHost cudaFreeHost
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuHostMallocDefault cudaHostAllocDefault

#define gpuEvent_t cudaEvent_t
#define gpuEventCreate cudaEventCreate
#define gpuEventRecord cudaEventRecord
#define gpuEventSynchronize cudaEventSynchronize
#define gpuEventElapsedTime cudaEventElapsedTime
#define gpuEventDestroy cudaEventDestroy

#endif
