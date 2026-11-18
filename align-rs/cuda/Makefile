# libalngpu.so: align-rs CUDA kernels, for Turing (sm_75) through Blackwell (sm_120), plus sm_120 PTX that the
# driver compiles for newer GPUs. A GPU outside these runs everything on the CPU (align-rs says so at start).
# One architecture only, for a faster build:  make -C cuda ARCH="-gencode arch=compute_120,code=sm_120"
NVCC     ?= nvcc
# the conda toolchain's g++ when there is one (matching its libstdc++), else the system's
CXX_HOST ?= $(shell command -v x86_64-conda-linux-gnu-g++ || echo g++)
ARCH     ?= $(foreach a,75 80 86 89 90 100 120,-gencode arch=compute_$(a),code=sm_$(a)) -gencode arch=compute_120,code=compute_120
NVFLAGS   = -O3 -std=c++17 $(ARCH) -Xcompiler -fPIC -ccbin $(CXX_HOST) -lineinfo --default-stream per-thread

build/libalngpu.so: aln_ext.cu
	@mkdir -p build
	$(NVCC) $(NVFLAGS) -shared -o $@ $< -cudart static

clean:
	rm -rf build
