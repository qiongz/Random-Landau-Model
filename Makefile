CC            = g++
ICC           = icpc
LD            = g++
NVCC          = nvcc
CFLAGS        = -Wall -std=c++11
ICC_LDFLAGS   = -Wall -qopenmp 
FPIC          = -fPIC

MAGMADIR     = /usr/local/magma
CUDADIR      = /usr/local/cuda
OPENBLASDIR  = /usr/local/openblas

MKL_CFLAGS       := -m64 -I$(MKLROOT)/include 

CUDA_INC  := -I$(CUDADIR)/include 

MAGMA_INC        := -DADD_ \
	            -I$(MAGMADIR)/include \
		    -I$(CUDADIR)/include

SEQ_MKL_LIBS     := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 \
       	  	     -lmkl_sequential -lmkl_core -lpthread -lm -ldl 

PARAL_MKL_LIBS   := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 \
	             -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

CUSOLVER_LIBS    := -L$(CUDADIR)/lib64 -lcublas -lcudart -lcusolver

MAGMA_LIBS       := -L$(MAGMADIR)/lib -lmagma -lmagma_sparse \
                    -L$(CUDADIR)/lib64 -lcublas -lcudart -lcusolver

NVCC_CFLAGS      := -gencode arch=compute_30,code=compute_30 \
                    -gencode arch=compute_32,code=compute_32 \
                    -gencode arch=compute_35,code=compute_35 \
                    -gencode arch=compute_50,code=compute_50 \
                    -gencode arch=compute_60,code=compute_60 \
                    -gencode arch=compute_61,code=compute_61 \
                    --device-c  

NVCC_LDFLAGS     := -gencode arch=compute_30,code=compute_30 \
                    -gencode arch=compute_32,code=compute_32 \
                    -gencode arch=compute_35,code=compute_35 \
                    -gencode arch=compute_50,code=compute_50 \
                    -gencode arch=compute_60,code=compute_60 \
                    -gencode arch=compute_61,code=compute_61 --device-link

# ----------------------------------------
default:
	@echo "Available make targets are:"
	@echo "  make cpu gpu       # compiles main_gpu.cu main_cpu.cpp"
	@echo "  make gpu           # compiles main_gpu.cu"
	@echo "  make cpu           # compiles main_cpu.cpp"
	@echo "  make clean         # deletes executables and object files"

all:cpu gpu clean

cpu:rlm_cpu rlm_cpu_dp clean

gpu:rlm_gpu_cusolver rlm_gpu_magma rlm_gpu_mkl rlm_gpu_magma_wfs clean

.PHONY:clean remove
clean:
	rm -f *.o 

remove:
	rm -f  rlm_gpu_* rlm_cpu

# ------ compile gpu version ------
# ---- object files -----
init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c -o $@ $<

wfs_file.o:wfs_file.cpp wfs_file.h
	$(CC) $(CFLAGS) -c -o $@ $<

matrix_coefficients.o:matrix_coefficients.cpp matrix_coefficients.h
	$(CC) $(CFLAGS) -c -o $@ $<

disorder_potential.o:disorder_potential.cpp disorder_potential.h
	$(CC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

potential_coeff.o:potential_coeff.cu
	$(NVCC) $(NVCC_CFLAGS) $(CUDA_INC) -c -o $@ $<

chern.o:chern.cpp chern.h wfs_file.h
	$(CC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

hamiltonian_gpu.o:hamiltonian_gpu.cu hamiltonian_gpu.h
	$(NVCC) $(NVCC_CFLAGS) $(CUDA_INC) -c -o $@ $<

cusolver_diag.o:cusolver_diag.cu cusolver_diag.h
	$(NVCC) $(NVCC_CFLAG) $(CUDA_INC) -c -o $@ $<

magma_diag.o:magma_diag.cu magma_diag.h 
	$(NVCC) $(MAGMA_CFLAGS) $(MAGMA_INC) $(NVCC_CFLAGS)  -c -o $@ $<

mkl_diag.o:mkl_diag.cpp mkl_diag.h
	$(CC) $(CFLAGS)  $(MKL_CFLAGS) -c -o $@ $<

main_gpu_cusolver.o:main_gpu.cu 
	$(NVCC) -Dcusolver $(NVCC_CFLAGS) $(CUDA_INC) -c -o $@ $<	

main_gpu_magma.o:main_gpu.cu
	$(NVCC) -Dmagma $(MAGMA_CFLAGS) $(MAGMA_INC) $(NVCC_CFLAGS) -c -o $@ $<

main_gpu_magma_wfs.o:main_gpu.cu
	$(NVCC) -Dmagma -DwfsIO $(MAGMA_CFLAGS) $(MAGMA_INC) $(NVCC_CFLAGS) -c -o $@ $<

main_gpu_mkl.o:main_gpu.cu 
	$(NVCC) -Dmkl $(NVCC_CFLAGS) $(CUDA_INC) -c -o $@ $<	

gpu_linker_cusolver.o:main_gpu_cusolver.o hamiltonian_gpu.o cusolver_diag.o potential_coeff.o
	$(NVCC) $(NVCC_LDFLAGS)  $^ -o $@

gpu_linker_magma.o:main_gpu_magma.o hamiltonian_gpu.o  magma_diag.o potential_coeff.o
	$(NVCC) $(NVCC_LDFLAGS)  $^ -o $@

gpu_linker_magma_wfs.o:main_gpu_magma_wfs.o hamiltonian_gpu.o  magma_diag.o potential_coeff.o
	$(NVCC) $(NVCC_LDFLAGS)  $^ -o $@

gpu_linker_mkl.o:main_gpu_mkl.o hamiltonian_gpu.o potential_coeff.o
	$(NVCC) $(NVCC_LDFLAGS)  $^ -o $@


# ----- executables -----
rlm_gpu_cusolver:main_gpu_cusolver.o gpu_linker_cusolver.o init.o matrix_coefficients.o disorder_potential.o hamiltonian_gpu.o chern.o cusolver_diag.o  wfs_file.o potential_coeff.o
	$(NVCC) -Dcusolver $^ $(PARAL_MKL_LIBS) ${CUSOLVER_LIBS} ${MAGMA_LIBS} -O3 -o $@

rlm_gpu_magma:main_gpu_magma.o gpu_linker_magma.o init.o matrix_coefficients.o disorder_potential.o hamiltonian_gpu.o chern.o magma_diag.o wfs_file.o potential_coeff.o
	$(NVCC) -Dmagma $^ $(PARAL_MKL_LIBS) ${MAGMA_LIBS} -O3 -o $@

rlm_gpu_mkl:main_gpu_mkl.o gpu_linker_mkl.o init.o matrix_coefficients.o disorder_potential.o hamiltonian_gpu.o chern.o mkl_diag.o wfs_file.o potential_coeff.o
	$(NVCC) -Dmkl $^ $(PARAL_MKL_LIBS) ${CUSOLVER_LIBS} ${MAGMA_LIBS} -O3 -o $@

rlm_gpu_magma_wfs:main_gpu_magma_wfs.o gpu_linker_magma_wfs.o init.o matrix_coefficients.o disorder_potential.o hamiltonian_gpu.o chern.o magma_diag.o wfs_file.o potential_coeff.o
	$(NVCC) -Dmagma -DwfsIO $^ $(PARAL_MKL_LIBS) ${MAGMA_LIBS} -O3 -o $@



# ------ compile cpu version ------

init_icc.o:init.cpp init.h
	$(ICC) $(CFLAGS) -c -o $@ $<
	
wfs_file_icc.o:wfs_file.cpp wfs_file.h
	$(ICC) $(CFLAGS) -c -o $@ $<

matrix_coefficients_icc.o:matrix_coefficients.cpp matrix_coefficients.h
	$(ICC) $(CFLAGS) -c -o $@ $<

disorder_potential_icc.o:disorder_potential.cpp disorder_potential.h
	$(ICC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

hamiltonian_cpu.o:hamiltonian_cpu.cpp hamiltonian_cpu.h
	$(ICC) $(MKL_CFLAGS) -c -o $@ $<

chern_icc.o:chern.cpp chern.h
	$(ICC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

mkl_diag_icc.o:mkl_diag.cpp mkl_diag.h
	$(ICC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

main_cpu.o:main_cpu.cpp 
	$(ICC) $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $< 

rlm_cpu:main_cpu.o hamiltonian_cpu.o init_icc.o disorder_potential_icc.o chern_icc.o matrix_coefficients_icc.o mkl_diag_icc.o wfs_file_icc.o
	$(ICC) $(ICC_LDFLAGS) $^ $(SEQ_MKL_LIBS) -O3 -o $@ 



# ------ compile cpu double precision version ------

init_icc_dp.o:init.cpp init.h
	$(ICC) -DDP $(CFLAGS) -c -o $@ $<
	
wfs_file_icc_dp.o:wfs_file.cpp wfs_file.h
	$(ICC) -DDP $(CFLAGS) -c -o $@ $<

matrix_coefficients_icc_dp.o:matrix_coefficients.cpp matrix_coefficients.h
	$(ICC) -DDP $(CFLAGS) -c -o $@ $<

disorder_potential_icc_dp.o:disorder_potential.cpp disorder_potential.h
	$(ICC) -DDP $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

hamiltonian_cpu_dp.o:hamiltonian_cpu.cpp hamiltonian_cpu.h
	$(ICC) -DDP $(MKL_CFLAGS) -c -o $@ $<

chern_icc_dp.o:chern.cpp chern.h
	$(ICC) -DDP $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

mkl_diag_icc_dp.o:mkl_diag.cpp mkl_diag.h
	$(ICC) -DDP $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $<

main_cpu_dp.o:main_cpu.cpp 
	$(ICC) -DDP $(CFLAGS) $(MKL_CFLAGS) -c -o $@ $< 

rlm_cpu_dp:main_cpu_dp.o hamiltonian_cpu_dp.o init_icc_dp.o disorder_potential_icc_dp.o chern_icc_dp.o matrix_coefficients_icc_dp.o mkl_diag_icc_dp.o wfs_file_icc_dp.o
	$(ICC) -DDP $(ICC_LDFLAGS) $^ $(SEQ_MKL_LIBS) -O3 -o $@ 
