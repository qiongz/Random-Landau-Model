#include "rlm/rlm.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

namespace rlm {
namespace {

[[noreturn]] void gpu_fail(const char* step, int status) {
    throw std::runtime_error(std::string(step) + " failed: status " + std::to_string(status));
}

[[noreturn]] void info_fail(const char* step, int info) {
    throw std::runtime_error(std::string(step) + " info != 0: " + std::to_string(info));
}

}  // namespace
}  // namespace rlm

#if defined(RLM_HIP)

#include <hip/hip_runtime.h>
#include <rocsolver/rocsolver.h>

namespace rlm {
namespace {

template <typename T>
class DeviceMem {
public:
    DeviceMem() = default;
    ~DeviceMem() { release(); }

    DeviceMem(const DeviceMem&) = delete;
    DeviceMem& operator=(const DeviceMem&) = delete;

    T* get() const { return ptr_; }

    void ensure(std::size_t count, const char* step) {
        if (count <= cap_) return;
        void* raw = nullptr;
        hipError_t st = hipMalloc(&raw, count * sizeof(T));
        if (st != hipSuccess) gpu_fail(step, static_cast<int>(st));
        if (ptr_ != nullptr) {
            hipError_t freed = hipFree(ptr_);
            if (freed != hipSuccess) {
                hipError_t drop = hipFree(raw);
                (void)drop;
                gpu_fail("hipFree", static_cast<int>(freed));
            }
        }
        ptr_ = static_cast<T*>(raw);
        cap_ = count;
    }

private:
    void release() {
        if (ptr_ == nullptr) return;
        hipError_t st = hipFree(ptr_);
        (void)st;
        ptr_ = nullptr;
        cap_ = 0;
    }

    T* ptr_ = nullptr;
    std::size_t cap_ = 0;
};

class RocHandle {
public:
    RocHandle() {
        rocblas_status st = rocblas_create_handle(&handle_);
        if (st != rocblas_status_success) {
            handle_ = nullptr;
            gpu_fail("rocblas_create_handle", static_cast<int>(st));
        }
    }

    ~RocHandle() {
        if (handle_ == nullptr) return;
        rocblas_status st = rocblas_destroy_handle(handle_);
        (void)st;
        handle_ = nullptr;
    }

    RocHandle(const RocHandle&) = delete;
    RocHandle& operator=(const RocHandle&) = delete;

    rocblas_handle get() const { return handle_; }

private:
    rocblas_handle handle_ = nullptr;
};

static_assert(sizeof(rocblas_double_complex) == 16,
              "rocblas_double_complex must be 16 bytes");
static_assert(sizeof(cx) == 16, "std::complex<double> must be 16 bytes");

class GpuSolver final : public Solver {
public:
    bool on_device() const override { return true; }

    void diagonalize(cx* H, double* eval, int n) override {
        if (n <= 0) throw std::runtime_error("diagonalize: n <= 0");

        const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
        a_.ensure(nn, "hipMalloc A");
        d_.ensure(static_cast<std::size_t>(n), "hipMalloc D");
        e_.ensure(static_cast<std::size_t>(n), "hipMalloc E");
        info_.ensure(1, "hipMalloc info");

        // Both types are two doubles. Byte-copy; do not convert element by element.
        hipError_t copy_in = hipMemcpy(a_.get(), H, nn * sizeof(cx), hipMemcpyHostToDevice);
        if (copy_in != hipSuccess) gpu_fail("hipMemcpy H2D", static_cast<int>(copy_in));

        rocblas_status solved = rocsolver_zheevd(handle_.get(),
                                                  rocblas_evect_original,
                                                  rocblas_fill_upper,
                                                  static_cast<rocblas_int>(n),
                                                  a_.get(),
                                                  static_cast<rocblas_int>(n),
                                                  d_.get(),
                                                  e_.get(),
                                                  info_.get());
        if (solved != rocblas_status_success) {
            gpu_fail("rocsolver_zheevd", static_cast<int>(solved));
        }

        hipError_t synced = hipDeviceSynchronize();
        if (synced != hipSuccess) gpu_fail("hipDeviceSynchronize", static_cast<int>(synced));

        rocblas_int info = 0;
        hipError_t copy_info = hipMemcpy(&info, info_.get(), sizeof(info), hipMemcpyDeviceToHost);
        if (copy_info != hipSuccess) gpu_fail("hipMemcpy info", static_cast<int>(copy_info));
        if (info != 0) info_fail("rocsolver_zheevd", static_cast<int>(info));

        hipError_t copy_vec = hipMemcpy(H, a_.get(), nn * sizeof(cx), hipMemcpyDeviceToHost);
        if (copy_vec != hipSuccess) gpu_fail("hipMemcpy eigenvectors", static_cast<int>(copy_vec));

        hipError_t copy_eval =
            hipMemcpy(eval, d_.get(), static_cast<std::size_t>(n) * sizeof(double), hipMemcpyDeviceToHost);
        if (copy_eval != hipSuccess) gpu_fail("hipMemcpy eigenvalues", static_cast<int>(copy_eval));
    }

private:
    RocHandle handle_;
    DeviceMem<rocblas_double_complex> a_;
    DeviceMem<double> d_;
    DeviceMem<double> e_;
    DeviceMem<rocblas_int> info_;
};

}  // namespace

std::unique_ptr<Solver> make_gpu_solver() {
    return std::make_unique<GpuSolver>();
}

}  // namespace rlm

#elif defined(RLM_CUDA)

#include <cuda_runtime.h>
#include <cusolverDn.h>

namespace rlm {
namespace {

template <typename T>
class DeviceMem {
public:
    DeviceMem() = default;
    ~DeviceMem() { release(); }

    DeviceMem(const DeviceMem&) = delete;
    DeviceMem& operator=(const DeviceMem&) = delete;

    T* get() const { return ptr_; }

    void ensure(std::size_t count, const char* step) {
        if (count <= cap_) return;
        void* raw = nullptr;
        cudaError_t st = cudaMalloc(&raw, count * sizeof(T));
        if (st != cudaSuccess) gpu_fail(step, static_cast<int>(st));
        if (ptr_ != nullptr) {
            cudaError_t freed = cudaFree(ptr_);
            if (freed != cudaSuccess) {
                cudaError_t drop = cudaFree(raw);
                (void)drop;
                gpu_fail("cudaFree", static_cast<int>(freed));
            }
        }
        ptr_ = static_cast<T*>(raw);
        cap_ = count;
    }

private:
    void release() {
        if (ptr_ == nullptr) return;
        cudaError_t st = cudaFree(ptr_);
        (void)st;
        ptr_ = nullptr;
        cap_ = 0;
    }

    T* ptr_ = nullptr;
    std::size_t cap_ = 0;
};

class CusolverHandle {
public:
    CusolverHandle() {
        cusolverStatus_t st = cusolverDnCreate(&handle_);
        if (st != CUSOLVER_STATUS_SUCCESS) {
            handle_ = nullptr;
            gpu_fail("cusolverDnCreate", static_cast<int>(st));
        }
    }

    ~CusolverHandle() {
        if (handle_ == nullptr) return;
        cusolverStatus_t st = cusolverDnDestroy(handle_);
        (void)st;
        handle_ = nullptr;
    }

    CusolverHandle(const CusolverHandle&) = delete;
    CusolverHandle& operator=(const CusolverHandle&) = delete;

    cusolverDnHandle_t get() const { return handle_; }

private:
    cusolverDnHandle_t handle_ = nullptr;
};

static_assert(sizeof(cuDoubleComplex) == 16, "cuDoubleComplex must be 16 bytes");
static_assert(sizeof(cx) == 16, "std::complex<double> must be 16 bytes");

class GpuSolver final : public Solver {
public:
    bool on_device() const override { return true; }

    void diagonalize(cx* H, double* eval, int n) override {
        if (n <= 0) throw std::runtime_error("diagonalize: n <= 0");

        const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
        a_.ensure(nn, "cudaMalloc A");
        w_.ensure(static_cast<std::size_t>(n), "cudaMalloc W");
        info_.ensure(1, "cudaMalloc info");

        cudaError_t copy_in = cudaMemcpy(a_.get(), H, nn * sizeof(cx), cudaMemcpyHostToDevice);
        if (copy_in != cudaSuccess) gpu_fail("cudaMemcpy H2D", static_cast<int>(copy_in));

        int lwork = 0;
        cusolverStatus_t sized = cusolverDnZheevd_bufferSize(handle_.get(),
                                                             CUSOLVER_EIG_MODE_VECTOR,
                                                             CUBLAS_FILL_MODE_UPPER,
                                                             n,
                                                             a_.get(),
                                                             n,
                                                             w_.get(),
                                                             &lwork);
        if (sized != CUSOLVER_STATUS_SUCCESS) {
            gpu_fail("cusolverDnZheevd_bufferSize", static_cast<int>(sized));
        }
        if (lwork < 0) throw std::runtime_error("cusolverDnZheevd_bufferSize failed: status -1");
        work_.ensure(static_cast<std::size_t>(lwork), "cudaMalloc work");

        cusolverStatus_t solved = cusolverDnZheevd(handle_.get(),
                                                    CUSOLVER_EIG_MODE_VECTOR,
                                                    CUBLAS_FILL_MODE_UPPER,
                                                    n,
                                                    a_.get(),
                                                    n,
                                                    w_.get(),
                                                    work_.get(),
                                                    lwork,
                                                    info_.get());
        if (solved != CUSOLVER_STATUS_SUCCESS) {
            gpu_fail("cusolverDnZheevd", static_cast<int>(solved));
        }

        cudaError_t synced = cudaDeviceSynchronize();
        if (synced != cudaSuccess) gpu_fail("cudaDeviceSynchronize", static_cast<int>(synced));

        int info = 0;
        cudaError_t copy_info = cudaMemcpy(&info, info_.get(), sizeof(info), cudaMemcpyDeviceToHost);
        if (copy_info != cudaSuccess) gpu_fail("cudaMemcpy info", static_cast<int>(copy_info));
        if (info != 0) info_fail("cusolverDnZheevd", info);

        cudaError_t copy_vec = cudaMemcpy(H, a_.get(), nn * sizeof(cx), cudaMemcpyDeviceToHost);
        if (copy_vec != cudaSuccess) gpu_fail("cudaMemcpy eigenvectors", static_cast<int>(copy_vec));

        cudaError_t copy_eval =
            cudaMemcpy(eval, w_.get(), static_cast<std::size_t>(n) * sizeof(double), cudaMemcpyDeviceToHost);
        if (copy_eval != cudaSuccess) gpu_fail("cudaMemcpy eigenvalues", static_cast<int>(copy_eval));
    }

private:
    CusolverHandle handle_;
    DeviceMem<cuDoubleComplex> a_;
    DeviceMem<double> w_;
    DeviceMem<cuDoubleComplex> work_;
    DeviceMem<int> info_;
};

}  // namespace

std::unique_ptr<Solver> make_gpu_solver() {
    return std::make_unique<GpuSolver>();
}

}  // namespace rlm

#else
#error "Define RLM_HIP or RLM_CUDA"
#endif
