# libcudapairhmm.so: the CUDA PairHMM variant-rs calls (C ABI).
#   make                      all GPU architectures below (fat binary + PTX for newer ones)
#   make ARCH="-gencode arch=compute_120,code=sm_120"   one architecture only (faster build)
NVCC     ?= nvcc
# the conda toolchain's g++ when there is one (matching its libstdc++), else the system's
CXX_HOST ?= $(shell command -v x86_64-conda-linux-gnu-g++ || echo g++)
ARCH     ?= $(foreach a,70 75 80 86 89 90 100 120,-gencode arch=compute_$(a),code=sm_$(a)) -gencode arch=compute_120,code=compute_120
# -ftz=true -fmad=false: GKL's numerics (flush-to-zero; separate multiply and add, as its AVX code)
NVFLAGS   = -O3 -std=c++17 $(ARCH) -Xcompiler -fPIC -ccbin $(CXX_HOST) -ftz=true -fmad=false -lineinfo $(EXTRA)

build/libcudapairhmm.so: native/cuda_pairhmm.cu
	@mkdir -p build
	$(NVCC) $(NVFLAGS) -shared -o $@ $< -cudart static

clean:
	rm -rf build
