# RELION SYCL/DPC++ version

This is SYCL/DPC++ version for [RELION](https://github.com/3dem/relion)

## Build & Running

+ Additional requirements for RELION SYCL/DPC++ version
	+ Intel(R) oneAPI Base Toolkit and HPC Toolkit (Installing all components are recommended. - https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
	+ Intel(R) software for general purpose GPU capabilities (https://dgpu-docs.intel.com)
	+ Intel(R) CPU Runtime for OpenCL(TM) Applications (optional - https://www.intel.com/content/www/us/en/developer/articles/technical/intel-cpu-runtime-for-opencl-applications-with-sycl-support.html)
	+ Codeplay(R) oneAPI for NVIDIA(R) GPU (optional - https://developer.codeplay.com/products/oneapi/nvidia/)
	+ Codeplay(R) oneAPI for AMD GPU (optional - https://developer.codeplay.com/products/oneapi/amd/)

+ SYCL specific command line arguments
	+ `--gpu`, `--sycl` or `--sycl-levelzero`: This is for Intel(R) GPU specific Level Zero backend (Recommended for Intel GPU)
	+ `--sycl-opencl`: This is for OpenCL(TM) GPU backend
	+ `--sycl-cpu`: This is for OpenCL(TM) CPU backend
	+ `--sycl-cuda`: This is for CUDA(R) backend over SYCL (Not tested)
	+ `--sycl-hip`: This is for HIP backend over SYCL (Not tested)
	+ If you need finer control of SYCL devices, you can set `ONEAPI_DEVICE_SELECTOR` environment variable. For detailed description, please look at https://intel.github.io/llvm-docs/EnvironmentVariables.html#oneapi-device-selector


```bash
$ git clone https://github.com/3dem/relion-devel.git relion_sycl -b sycl-merge
$ cd relion_sycl; mkdir build_sycl; cd build_sycl
$ {Load Intel oneAPI toolkit and SYCL/Level Zero/OpenCL runtime environment}
$ sycl-ls					# This will display available SYCL devices
$ cmake \
-DCMAKE_C_COMPILER=mpiicx \
-DCMAKE_CXX_COMPILER=mpiicpx \
-DMPI_C_COMPILER=mpiicx \
-DMPI_CXX_COMPILER=mpiicpx \
-DMKLFFT=ON \
-DSYCL=ON \
-DCUDA=OFF \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_FLAGS="-O3 -fsigned-zeros" \
-DCMAKE_CXX_FLAGS="-O3 -fsigned-zeros -march=native" \
-DCMAKE_EXE_LINKER_FLAGS="-O3 -march=native" \
..

$ #### This is Intel GPU Level Zero backend specific #####
$ export ZE_AFFINITY_MASK=0 # Use only the first available Level Zero device
$ export ZEX_NUMBER_OF_CCS=0:4,1:4
$ export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=2 # Don't use this with Intel Arc GPUs. Only for Max GPUs
$ export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
$ export SYCL_PI_LEVEL_ZERO_USM_ALLOCATOR="1;4G;host:16M,512,64K;device:16M,1024,64K;shared:0,0,64K"
$ #### End of Intel GPU Level Zero backend specific  #####
$ # For finer control of SYCL devcices, please see the above descrpition on ONEAPI_DEVICE_SELECTOR
$ 
$ ulimit -n 512000 # this is necessary for multi-GPU jobs
$ {Run 2D/3D/refinement application by replacing --gpu/--cpu with --gpu/--sycl/--sycl-opencl/--sycl-cpu/--sycl-cuda/--sycl-hip}
```


## Added macros

+ For CMake configuration
	+ `SYCL`(=ON/OFF): Enable SYCL based acceleration build
	+ `SyclForceOneDPL`(=ON/OFF): Use oneDPL(https://github.com/oneapi-src/oneDPL) if it can be used. This has the same effect as setting "-DSYCL_OFFLOAD_SORT" for CMAKE_CXX_FLAGS below. (experimental)
	+ `SYCL_AOT_COMPILE`(=ON/OFF): Enable AOT(Ahead-Of-Time) compilation for SPIR64 target. Default target is pvc. (experimental)
	+ `SYCL_AOT_TARGET`(=ON/OFF): Specify AOT(Ahead-Of-Time) SPIR64 target. Possible list can be checked using "ocloc compile --help" command. (experimental)
	+ `SYCL_CUDA_COMPILE`(=ON/OFF): Enable SYCL compilation for CUDA target (Not tested)
	+ `SYCL_CUDA_ARCH`: SYCL CUDA arch target (Not tested)
	+ `SYCL_HIP_COMPILE`(=ON/OFF): Enable SYCL compilation for HIP target (Not tested)
	+ `SYCL_HIP_ARCH`: SYCL HIP arch target (Not tested)
	+ `SYCL_ALL_COMPILE`(=ON/OFF): Enable SYCL compilation for SPIR64, CUDA, and HIP target (Not tested)
	+ `SYCL_HOST_FLAGS`(=list of flags with space as separator): Additional flags for host compiler (for future use)
	+ `SYCL_COMPILER_NAME`: SYCL compiler command name (for future use)
	+ `SYCL_COMPILE_FLAGS`: Additional SYCL compile flags (for future use)
	+ `SYCL_LINK_FLAGS`: SYCL link flags except "-lsycl -lOpenCL" if needed (for future use)

+ Others for testing purpose (just defining with -D* is needed in cmake -DCMAKE_CXX_FLAGS)
	+ `SYCL_OFFLOAD_SORT`: Use SYCL kernel for weight sorting routines. If this is set, oneDPL(https://github.com/oneapi-src/oneDPL) is used when it is beneficial. (experimental)
		+ `USE_LESS_ONEDPL`: If this is set, oneDPL is used only when there is no other implementation.
		+ `PREFER_ONEDPL`: If this is set, oneDPL is used everywhere when it is applicable. If this is set, you SHOULD NOT SET DisableIndirectAccess=1.
	+ `SYCL_OFFLOAD_FFT`: Use SYCL kernel for the current FFTW routines. (Not implemented)
	+ `INTEL_EXPLICIT_SIMD`: Use Explicit SIMD extension for SYCL kernels. (Not implemented)
	+ `INTEL_SG_SIZE`: Used for Intel sub-group size in SYCL kernel. 32 is recommended for PVC and 16 is for ATS. (Not tested well)
	+ `USE_IPP`: Use Intel IPP library's RadixSort for sortOnDevice instead of std::sort. Enabled by default if IPP library exists.
	+ `USE_MPI_COLLECTIVE`: Use MPI collective whenever possible. Enabled by default for ALTCPU and SYCL.
	+ `USE_INORDER_QUEUE`: Use in-order SYCL queue. Without this, out-of-order SYCL queue is used by default. (experimental)
	+ `USE_ASYNC_SYCL_SUBMIT`: Remove wait() for each SYCL kernel submission. (experimental)
	+ `USE_SYCL_STREAM`: Create new in-order SYCL queue for each cudaStream. (experimental)
	+ `USE_SUBSUB_DEVICE`: Create separate SYCL queue for each CCS. (experimental)
	+ `USE_EXISTING_SYCL_DEVICE`: This will copy and use created SYCL device pointer instead of creating new SYCL device for each thread. Not recommended.
	+ `USE_SINCOS_TABLE`: Pre-calculate sine/cosine table before main loop in some kernels (Not implemented)

+ Macros defined in source code (No need to specify in CMake configuration and it is defined by -DSYCL=ON above)
	+ `_SYCL_ENABLED`: source codes with SYCL which may be compiled with any SYCL compiler
	+ `_DPCPP_ENABLED`: source codes with DPC++ which should be compiled with Intel oneAPI compiler only

+ Prefined macros by Intel(R) oneAPI compiler
	+ `__INTEL_LLVM_COMPILER`: For icx and icpx

