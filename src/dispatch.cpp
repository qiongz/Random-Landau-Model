#include "rlm/rlm.hpp"

#include <memory>
#include <stdexcept>

namespace rlm {

std::unique_ptr<Solver> make_solver(Backend backend) {
    switch (backend) {
        case Backend::Mkl:
            return make_mkl_solver();
        case Backend::Gpu:
#ifdef RLM_HAS_GPU
            return make_gpu_solver();
#else
            throw std::runtime_error(
                "GPU solver is not built; rebuild with -DRLM_GPU=HIP or -DRLM_GPU=CUDA");
#endif
    }
    throw std::runtime_error("unknown solver backend");
}

}  // namespace rlm
