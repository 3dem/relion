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

## Options for Accelerated versions of relion\_refine\*

Depending on your build options, you can build two different RELION binaries.
One will contain the original CPU version of RELION plus the CUDA-accelerated
version, invoked with the "`--gpu`" flag.   The other will contain the original
CPU version of RELION plus the CPU-accelerated version, invoked with the "`--cpu`"
flag.  Unfortunately, all three versions cannot be built into the same binaries
at this time.  So, if you wish to have all three available at
your site, you will need to copy the binaries at `<RELION_root>/build/bin` and the
shared library `librelion_lib.so` at `<RELION_root>/build/lib` into a holding directory,
clean out the build area, set up and build with a different set of build options, and then
copy the new binaries and `librelion_lib.so` into a different holding directory.

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
* The `FORCE_OWN_TBB=ON` option ensures you are building with the latest version of
Threading Building Blocks (TBB), an open-source library originally development by
Intel Corporation to implement a task-based parallelism.  It has been found that
the version of TBB installed on some systems is missing functionality found in
later versions.   Using this flag should remove any build or link issues related
to TBB.
* The `MKLFFT=ON` option builds relion\_refine\* with the Intel(R) Math Kernel Library
(Intel(R) MKL), which has an FFT implementation that was found to scale better with
increased threads than the one found in FFTW.
* Best performance of the CPU-accelerated kernels enabled with "`--cpu`" will be
obtained using [Intel(R) Parallel Studio XE 2018 Cluster Edition](https://software.intel.com/en-us/parallel-studio-xe).  With the Plasmodium ribosome benchmark noted on the [RELION website](https://www2.mrc-lmb.cam.ac.uk/relion/index.php?title=Benchmarks_%26_computer_hardware),
GCC 7.3 hardware-optimized builds appear to run about 3x slower than those built
with Intel(R) Parallel Studio XE 2018 Cluster Edition (for reasons still under
investigation).
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
Note that this binary will only run on AXV2 platforms.
    * Best performance on AVX512 platforms has been measured using the following CMAKE command:
```
CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xCORE-AVX512 -qopt-zmm-usage=high -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
```
Note that the resulting binary will only work on AVX512 platforms.
    * Best performance on Intel(R) Xeon Phi(TM) Processors (x700 Family) has been measured
     using the following CMAKE command:
     ```
     CC=mpiicc CXX=mpiicpc cmake -DCUDA=OFF -DALTCPU=ON -DCudaTexture=OFF -DMKLFFT=ON -D CMAKE_C_FLAGS="-O3 -ip -g -debug inline-debug-info -xMIC-AVX512 -restrict " -D CMAKE_CXX_FLAGS="-O3 -ip -g -debug inline-debug-info -xMIC-AVX512 -restrict " -DGUI=OFF -D CMAKE_BUILD_TYPE=Release ..
     ```
     Note that this binary will only run on Intel(R) Xeon Phi(TM) Processor family platforms.
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
* Best performance seems to be seen when the pool size (`--pool`) is about twice
the number of threads (`--j`), to about half the number of threads.  Lowering the
pool size may also decrease the memory used by a process or rank.  When running
multi-node, 4 MPI ranks per node seems to also work well (so `--j` should be set
to the total number of threads in the machine divided by 4, and `--pool` adjusted
accordingly), although how many ranks you use per machine may also depend on the data set.
* Systems with 256GB or 512GB are recommended for the CPU-accelerated kernels, with
the most memory utilization being seen at the last iteration of autorefinement.
For example, during a multi-rank MPI run, a single process/rank working on an
800x800 image stack was seen to have a peak memory use of about 170GB during this
final iteration.  By contrast, the previous iterations of the same data set have a
peak memory utilization of about 90 GB per rank (the same was seen during 3D
classification).   This data set would be run one rank per node on a 256GB system
due to the memory requirements of this final iteration.
By contrast, a data set with 300x300 images (Plasmodium ribosome) can run with
4 ranks per 256GB system without swapping.
* To run the Plasmodium ribosome data set workload found at the [RELION benchmark page](https://www2.mrc-lmb.cam.ac.uk/relion/index.php?title=Benchmarks_%26_computer_hardware)
on a machine with two sockets populated by AVX512 CPUs containing 20 cores each supporting
Intel(R) Hyper-Threading Technology (so two threads per core for a total of 80 threads -
Intel(R) Xeon(R) Gold 6148 processors):
  * 1 machine - completes in about 4.7 hours when built with Intel(R) Parallel
Studio XE 2018 Cluster Edition and the options listed above:
```
export OMP_SCHEDULE="dynamic"
export KMP_BLOCKTIME=0
relion_refine --o Results/skx --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --ctf --ctf_corrected_ref --tau2_fudge 4 --particle_diameter 360 --K 6 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --random_seed 0 --dont_combine_weights_via_disc --preread_images --iter 25 --pool 160 --j 80 --cpu
```
  * 2 machines - completes in about 2.6 hours when built with Intel(R) Parallel
Studio XE 2018 Cluster Edition and the options listed above:
  ```
export OMP_SCHEDULE="dynamic"
export KMP_BLOCKTIME=0
mpirun -perhost 4 -np 8 relion_refine_mpi --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --dont_combine_weights_via_disc --ctf --ctf_corrected_ref --tau2_fudge 4 --particle_diameter 360 --K 6 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --random_seed 0 --o Results/cpu --pool 40 --j 20 --iter 25 --cpu
  ```
  * 4 machines - completes in about 1.4 hours when built with Intel(R) Parallel
Studio XE 2018 Cluster Edition and the options listed above:
```
export OMP_SCHEDULE="dynamic"
export KMP_BLOCKTIME=0
mpirun -perhost 4 -np 16 relion_refine_mpi --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --dont_combine_weights_via_disc --ctf --ctf_corrected_ref --tau2_fudge 4 --particle_diameter 360 --K 6 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --random_seed 0 --o Results/cpu --pool 40 --j 20 --iter 25 --cpu
```
  * 8 machines - completes in about 1 hour when built with Intel(R) Parallel
Studio XE 2018 Cluster Edition and the options listed above:
```
export OMP_SCHEDULE="dynamic"
export KMP_BLOCKTIME=0
mpirun -perhost 4 -np 32 relion_refine_mpi --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --dont_combine_weights_via_disc --ctf --ctf_corrected_ref --tau2_fudge 4 --particle_diameter 360 --K 6 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --random_seed 0 --o Results/cpu --pool 40 --j 20 --iter 25 --cpu
```


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
