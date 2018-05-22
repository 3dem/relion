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


More extensive options and configurations are available
[here](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install),
but the outlines to clone and install relion for typical use are made easy
through [cmake](https://en.wikipedia.org/wiki/CMake).

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

## Options for Accelerated versions

### CUDA-accelerated kernels
To build with support for CUDA-accelerated kernels in addition to the original CPU version,
build by setting `CUDA=ON` and `CudaTexture=ON/OFF`:

```
cd build
rm -r *
cmake -DCUDA=ON -DCudaTexture=ON ...
make -j4
make install
```
NOTE:  If you run relion\_refine\* with a the "`--gpu`" flag, you will run the
accelerated CUDA version of the kernels.   If you leave out the "`--gpu`" flag,
it will run the original CPU version.

### CPU-accelerated kernels
To build with support for CPU-accelerated kernels in addition to the original
CPU version, build by setting `CUDA=OFF`, `ALTCPU=ON`, `CudaTexture=OFF`, and
optionally `MKLFFT=ON` and/or `FORCE_OWN_TBB=ON`:
```
cd build
rm -r *
cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -DFORCE_OWN_TBB=ON ...
make -j4
make install
```
NOTES:  
* If you run relion\_refine\* with a the "`--cpu`" flag, you will run the
accelerated CPU version of the kernels.   If you leave out the "`--cpu`" flag,
it will run the original CPU version.
* If building with GCC, make sure you source some MPI installation first, such as
[Open MPI](https://www.open-mpi.org/) or
[Intel MPI](https://software.intel.com/en-us/intel-mpi-library).  You may also
minimize difficulties by specifying `FORCE_OWN_TBB=ON`
* Best performance of the CPU-accelerated kernels enabled with "`--cpu`" will be
obtained using [Intel Parallel Studio XE 2018 Cluster Edition](https://software.intel.com/en-us/parallel-studio-xe).  
  * Before running cmake, set up the compile and build environment by sourcing the appropriate files:
```
source /opt/intel/compilers_and_libraries_2018.2.199/linux/bin/compilervars.sh intel64
source /opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/bin/mklvars.sh intel64
source /opt/intel/impi/2018.2.199/intel64/bin/mpivars.sh intel64
```
  * Now, set up the project explicitly specifying the Intel compiler, any Intel-specific
compiler flags, and all the CMAKE flags for the accelerated version
(`FORCE_OWN_TBB=ON` is not needed since the build will pick up the version of
TBB that ships with the Intel(R) C++ Compiler).   For example:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info" ...
```
    * Best performance on AVX2 platforms has been measured using the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX2 -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX2 -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
    * Best performance on AVX512 platforms has been measured using the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
  * If you want to run a version of RELION built with the Intel(R) C++ Compiler
  on another machine or cluster, you can download the runtime libraries for free
  for use on those machines from the following sources:
    *  [Intel (R) C++ Compiler libraries](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2018-compilers-for-linux) including Threading Building Blocks
      * Decompress the .tgz file
      * Create a directory to hold all the Intel components somewhere that all
      machines using the binary can access
      * Install by running `install.sh` and selecting the above holding directory
    *  [Intel(R) Math Kernel Libraries and Intel(R) MPI](https://software.intel.com/en-us/mkl)
      *  Register
      *  Now download at least "Intel(R) Math Kernel Library (Intel MKL)" and "Intel
      MPI Library (Linux Package)" (make sure the product "Intel(R) Performance
      Libraries for Linux*" is selected at the top)
      *  Install each by running `install.sh` and select a user-level installation
      that goes to your holding directory
    * Once done, `source <holding_dir>/compilers_and_libraries_2018.2.199/linux/bin/compilervars.sh intel64`
    should set up the runtime environment on your machine without a local Intel software installation.


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
