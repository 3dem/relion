RELION
======


RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

The underlying theory of MAP refinement is given in a [scientific publication](https://www.ncbi.nlm.nih.gov/pubmed/22100448)
. If RELION is useful in your work, please cite this paper.


The more comprehensive documentation of RELION is stored on the [Wiki](http://www2.mrc-lmb.cam.ac.uk/relion)

## Installation


More extensive options and configurations are available [here](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install), but the outlines to clone and install relion for typical use are made easy through [cmake](https://en.wikipedia.org/wiki/CMake).

On ubuntu machines, installing cmake, the compiler, and additional dependencies (mpi, fftw) is as easy as:

```
sudo apt install cmake build-essential mpi-default-bin mpi-default-dev libfftw3-dev
```

On other systems it is typically just as easy, you simply have to modify "apt" to
the appropriate package manager. You will also need [git](https://en.wikipedia.org/wiki/Git), which is just as easy;

```
sudo apt install git
```


Once git and cmake are installed, relion can be easily installed through
```
git clone https://github.com/3dem/relion.git
cd relion
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/where/to/install/ ..
make -j4
make install
```
(NOTES: "/where/to/install/.." is typically "/usr/local/relion".
 Make sure you create this directory beforehand.
 Installing to that location requires sudo, so in this case be sure to use
 "sudo make install" instead of "make install" in the last step.)

These steps will download the source-code, create a build-directory,
then configure and build relion, and lastly install it to be generally
available on the system.
## Updating


RELION is intermittently updated, with both minor and major features.
To update an existing installation, simply use the following commands

```
cd relion
git pull
cd build
make -j4
make install    # (or "sudo make install")
```



## Options for accelerated versions

Parts of the cryo-EM processing pipeline can be very computationally demanding, and in some cases special hardware can be used to make these faster. There are two such cases at the moment; 

* Since RELION-2: Use one or more PGUs, or graphics cards. RELION only supports CUDA-capable GPUs of compute capabilty 3.5 or higher.  
* Since RELION-3: Use the vectorized version. RELION only supports GCC and ICC 2018.3 or later.

There are more benefits than speed; the accelearated versions also have a decreased memory footprint. Details about how to enable either of these options is listed below. 

### Note that... 
If you are using intels compiler (icc), you cannot have both these features in the same binary at the moment. Not all tools are GPU-accelerated or vectorized. 

## GPU-acceleration
Tools that are GPU-accelerated: 
* relion\_refine
* relion\_autopick

To build with support for CUDA-accelerated kernels in addition to the original CPU version, build by setting `CUDA=ON` and `CudaTexture=ON/OFF`:

```
cd build
rm -r *
cmake -DCUDA=ON -DCudaTexture=ON ...
make -j4
make install
```
### Use
If you run relion\_refine with a the "`--gpu`" flag, you will run the accelerated CUDA version of the kernels.   If you leave out the "`--gpu`" flag, it will run the original CPU version.

## CPU-acceleration
Tools that are CPU-accelerated (vectorized): 
* relion\_refine

To build with support for CPU-accelerated kernels in addition to the original
CPU version, build by setting `ALTCPU=ON`
```
cd build
rm -r *
cmake -DALTCPU=ON ..
make -j4
make install
```
This will require intels TBB (thread building blocks) library. RELION will ook for TBB and fetch+install it if it cannot find it on your system. You can force this behaviour (and make sure you are using the latest version) by adding 
```
-DFORCE_OWN_TBB=ON
```
In addition, you can make use if intels MKL (Math Kernel Library). This is optional ( but will scale better with increased threads), and added through 
```
-DMKLFFT=ON
```





### Use
If you run relion\_refine with a the "`--cpu`" flag, you will run the accelerated version. If you leave it the original CPU version will be run. YOu should use this flag if you can, unless you want to verify old runs or behaviour.

With the Plasmodium ribosome benchmark noted on the [RELION website](https://www2.mrc-lmb.cam.ac.uk/relion/index.php?title=Benchmarks_%26_computer_hardware), GCC 7.3 hardware-optimized builds appear to run about 3x slower than those built with Intel(R) Parallel Studio XE 2018 Cluster Edition (for reasons still under investigation).

### Further considerations and hardware-specific compilation 

Before running cmake, you may need to set up the compile and build environment by sourcing the appropriate files:
```
source /opt/intel/compilers_and_libraries_<version>/linux/bin/compilervars.sh intel64
source /opt/intel/compilers_and_libraries_<version>/linux/mkl/bin/mklvars.sh intel64
source /opt/intel/impi/<version>/intel64/bin/mpivars.sh intel64
```

## AVX2
Best performance on AVX2 platforms has been measured using the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX2 -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX2 -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
Note that this binary will only run on AXV2 platforms.

## AVX512
Best performance on AVX512 platforms has been measured using the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
Note that the resulting binary will only work on AVX512 platforms.
## Xeon Phi 
Best performance on Intel Xeon Phi Processors (x700 Family) has been measuredusing the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xMIC-AVX512 -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xMIC-AVX512 -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
Note that this binary will only run on Xeon Phi Processor family platforms.

## Getting the intel tools

The intel compiler (icc) is not free. However, if you want to *run* a version of RELION built with the Intel compiler, you can download the runtime libraries for free:
  *  the [Intel C++ Compiler libraries](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2018-compilers-for-linux) (including TBB)
    * the [Intel Math Kernel Libraries and Intel MPI](https://software.intel.com/en-us/mkl)
 
The above nneds to be installed in a location that is accessible to all machines running the icc-compiled binary.

## Optimizing run-settings for CPU-acceleration
Best performance seems to be seen when the pool size (`--pool`) is roughly the same as the number of threads (`--j`).  Lowering the pool size may also decrease the memory used by a process or rank.  When running multi-node, 4 MPI ranks per node seems to also work well (so `--j` and `--pool` should be set  to the total number of threads in the machine divided by 4), although this may also depend on the data set. Systems with 256GB or 512GB are recommended for the CPU-accelerated kernels.

