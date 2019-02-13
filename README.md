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

If FLTK related errors are reported, please add `-DFORCE_OWN_FLTK=ON` to
`cmake`. For FFTW related errors, try `-DFORCE_OWN_FFTW=ON`.

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
This will require the TBB (Threading Building Blocks) library. RELION will look for TBB and fetch+install it if it cannot find it on your system. You can force this behaviour (and make sure you are using the latest version) by adding
```
-DFORCE_OWN_TBB=ON
```
In addition, you can make use the Intel(R) Math Kernel Library (Intel(R) MKL). This is optional ( but will scale better with increased threads).   Add this by:
```
-DMKLFFT=ON
```

### Use
If you run relion\_refine with a the "`--cpu`" flag, you will run the accelerated version. If you leave it the original CPU version will be run. You should use this flag if you can, unless you want to verify old runs or behaviour.


### Building with the Intel(R) Compiler
With the Plasmodium ribosome benchmark noted on the [RELION website](https://www2.mrc-lmb.cam.ac.uk/relion/index.php?title=Benchmarks_%26_computer_hardware), GCC 7.3 hardware-optimized builds appear to run about 3x slower than those built with Intel(R) Parallel Studio XE 2018 Cluster Edition (for reasons still under investigation).

To build with the Intel compiler:
```
source /opt/intel/impi/<version>/intel64/bin/mpivars.sh intel64cd build
source /opt/intel/compilers_and_libraries_<version>/linux/bin/compilervars.sh intel64
source /opt/intel/compilers_and_libraries_<version>/linux/mkl/bin/mklvars.sh intel64
rm -r *
CC=mpiicc CXX=mpiicpc cmake -DALTCPU=ON -DMKLFFT=ON ..
make -j4
make install
```

### Getting the Intel runtime libraries

The Intel(R) Compiler (icc) is not free. However, if you want to *run* a version of RELION built with the Intel compiler, you can download the runtime libraries for free:
  *  The [Intel C++ Compiler libraries](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2018-compilers-for-linux) (including TBB)
  * The [Intel Math Kernel Libraries and Intel MPI](https://software.intel.com/en-us/mkl)
      *  Download at least "Intel(R) Math Kernel Library (Intel MKL)" and "Intel
      MPI Library (Linux Package)" (make sure the product "Intel(R) Performance
      Libraries for Linux*" is selected at the top)

The above needs to be installed in a location that is accessible to all machines running the icc-compiled binary.

## Optimized run-settings for CPU-acceleration
Best performance seems to be seen when the pool size (`--pool`) is roughly the same as the number of threads (`--j`).  Lowering the pool size may also decrease the memory used by a process or rank.  When running multi-node, 4 MPI ranks per node seems to also work well (so `--j` and `--pool` should be set  to the total number of threads in the machine divided by 4), although this may also depend on the data set. Systems with 256GB or 512GB are recommended for the CPU-accelerated kernels.
