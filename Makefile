CUDA_PATH	=	/usr/local/apps/cuda/cuda-9.2
CUDA_BIN_PATH	=	$(CUDA_PATH)/bin
CUDA_NVCC	=	$(CUDA_BIN_PATH)/nvcc

project6:	project6.cu
		$(CUDA_NVCC) ${CXXFLAGS} -o project6  project6.cu randomNumberFunctions.cpp -lm &> compile_errors.txt 
