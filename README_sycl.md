# RELION SYCL/DPC++ version

This is SYCL/DPC++ version of [RELION](https://github.com/3dem/relion)

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
$ mkdir build_sycl; cd build_sycl
$ {Load Intel oneAPI toolkit and SYCL/Level Zero/OpenCL runtime environment}
$ sycl-ls                        # This will display available SYCL devices
$ cmake \
-DCMAKE_C_COMPILER=mpiicx \
-DCMAKE_CXX_COMPILER=mpiicpx \
-DMPI_C_COMPILER=mpiicx \
-DMPI_CXX_COMPILER=mpiicpx \
-DMKLFFT=ON \
-DSYCL=ON \
-DCUDA=OFF \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_FLAGS="-O3" \
-DCMAKE_CXX_FLAGS="-O3 -march=native" \
..

$ #### This is Intel GPU Level Zero backend specific #####
$ export ZE_AFFINITY_MASK=0         # Use only the first available Level Zero device. This can be replaced by --gpu 0 syntax.
$ export ZEX_NUMBER_OF_CCS=0:4,1:4  # Set this only if you are putting more than one MPI ranks per GPU. 0:4 means 4 MPI ranks running on card 0
$ #### End of Intel GPU Level Zero backend specific  #####
$ # For finer control of SYCL devcices, please see the above descrpition on ONEAPI_DEVICE_SELECTOR
$ {Run 2D/3D/refinement application by replacing --gpu/--cpu with --gpu/--sycl/--sycl-opencl/--sycl-cpu/--sycl-cuda/--sycl-hip}
```


## Optional runtime environment variables

+ The below shell environment variables can be tested for more potential SYCL specific tuning. Setting it to "1" or "on" will enable these features.
	+ `relionSyclUseCuda`: This with --gpu have the same meaning as --sycl-cuda. --sycl-cuda will be used even if --gpu/--sycl is specified in command lines.
	+ `relionSyclUseHip`: This with --gpu have the same meaning as --sycl-hip. --sycl-hip will be used even if --gpu/--sycl is specified in command lines
	+ `relionSyclUseInOrderQueue`: Use in-order SYCL queue. Without this, out-of-order SYCL queue is used by default. (experimental)
	+ `relionSyclUseAsyncSubmission`: Remove wait() for each SYCL kernel submission. (experimental)
	+ `relionSyclUseStream`: Create new in-order SYCL queue for each cudaStream. (experimental)
	+ `relionSyclUseSubSubDevice`: Create separate SYCL queue for each CCS. (experimental)
	+ `relionSyclBlockSize`: SYCL memory pool block size. This takes precedence over relionSyclHostBlockSize and relionSyclDeviceBlockSize.
	+ `relionSyclHostBlockSize`: SYCL memory pool block size for sycl::malloc_host. Default is 256MB
	+ `relionSyclDeviceBlockSize`: SYCL memory pool block size for sycl::malloc_device. Default is 256MB
	+ `MAX_MPI_BLOCK`: Maximum MPI message size per single MPI API call for point-to-point and collective communication. This takes precedence over MAX_MPI_P2P_BLOCK and MAX_MPI_COLL_BLOCK.
	+ `MAX_MPI_P2P_BLOCK`: Maximum MPI message size per single MPI API call for point-to-point communication. Default is 4GB.
	+ `MAX_MPI_COLL_BLOCK`: Maximum MPI message size per single MPI API call for collective communication. Default is 64MB.


## Added macros

+ For CMake configuration
	+ `SYCL`(=ON/OFF): Enable SYCL based acceleration build
	+ `SYCL_CUDA_COMPILE`(=ON/OFF): Enable SYCL compilation for CUDA target
	+ `SYCL_CUDA_TARGET`: SYCL CUDA arch target (i.e. 80)
	+ `SYCL_HIP_COMPILE`(=ON/OFF): Enable SYCL compilation for HIP target
	+ `SYCL_HIP_TARGET`: SYCL HIP arch target (i.e gfx90a)
	+ `SyclForceOneDPL`(=ON/OFF): Use oneDPL(https://github.com/oneapi-src/oneDPL) if it can be used. This has the same effect as setting "-DUSE_ONEDPL" for CMAKE_CXX_FLAGS below. (experimental)
	+ `SYCL_AOT_COMPILE`(=ON/OFF): Enable AOT(Ahead-Of-Time) compilation for SPIR64 target. Default target is pvc. (for future use)
	+ `SYCL_AOT_TARGET`(=ON/OFF): Specify AOT(Ahead-Of-Time) SPIR64 target. Possible list can be checked using "ocloc compile --help" command. (for future use)
	+ `SYCL_HOST_FLAGS`(=list of flags with space as separator): Additional flags for host compiler (for future use)
	+ `SYCL_COMPILER_NAME`: SYCL compiler command name (for future use)
	+ `SYCL_COMPILE_FLAGS`: Additional SYCL compile flags (for future use)
	+ `SYCL_LINK_FLAGS`: SYCL link flags except "-lsycl -lOpenCL" if needed (for future use)

+ Others for testing purpose (just defining with -D* is needed in cmake -DCMAKE_CXX_FLAGS)
	+ `USE_ONEDPL`: Use SYCL kernel for weight sorting routines. If this is set, oneDPL(https://github.com/oneapi-src/oneDPL) is used when it is beneficial. (experimental)
		+ `USE_LESS_ONEDPL`: If this is set, oneDPL is used only when there is no other implementation.
		+ `USE_MORE_ONEDPL`: If this is set, oneDPL is used everywhere when it is applicable. If this is set, you SHOULD NOT SET DisableIndirectAccess=1.
	+ `SYCL_OFFLOAD_FFT`: Use SYCL kernel for the current FFTW routines. (Not implemented)
	+ `INTEL_EXPLICIT_SIMD`: Use Explicit SIMD extension for SYCL kernels. (Not implemented)
	+ `INTEL_SG_SIZE`: Used for Intel sub-group size in SYCL kernel. 32 is recommended for PVC and 16 is for ATS. (Not tested well)
	+ `USE_IPP`: Use Intel IPP library's RadixSort for sortOnDevice instead of std::sort. Enabled by default if IPP library exists.
	+ `USE_MPI_COLLECTIVE`: Use MPI collective whenever possible. Enabled by default for ALTCPU and SYCL.
	+ `USE_EXISTING_SYCL_DEVICE`: This will copy and use created SYCL device pointer instead of creating new SYCL device for each thread. Not recommended.
	+ `USE_SINCOS_TABLE`: Pre-calculate sine/cosine table before main loop in some kernels (Not implemented)

+ Macros defined in source code (No need to specify in CMake configuration and it is defined by -DSYCL=ON above)
	+ `_SYCL_ENABLED`: source codes with SYCL which may be compiled with any SYCL compiler
	+ `_DPCPP_ENABLED`: source codes with DPC++ which should be compiled with Intel oneAPI compiler only

+ Prefined macros by Intel(R) oneAPI compiler
	+ `__INTEL_LLVM_COMPILER`: For icx and icpx

