RELION 3.1 (+ EER support)
===========================

RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

The underlying theory of MAP refinement is given in a [scientific publication](https://www.ncbi.nlm.nih.gov/pubmed/22100448).
If RELION is useful in your work, please cite this paper.

The more comprehensive documentation of RELION is stored on the [Wiki](http://www2.mrc-lmb.cam.ac.uk/relion)

## Notes specific to this EER branch: PLEASE READ

This is a special branch for processing movies in the EER format.
Once the development is complete, this will be merged into the main branch.

Please read instructions in [our wiki](https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Image_compression#Falcon4_EER) carefully before starting processing.

## Installation

More extensive options and configurations are available [here](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install),
but the outlines to clone and install relion for typical use are made easy through [cmake](https://en.wikipedia.org/wiki/CMake).

On Debian or Ubuntu machines, installing cmake, the compiler, and additional dependencies (mpi, fftw) is as easy as:

```
sudo apt install cmake git build-essential mpi-default-bin mpi-default-dev libfftw3-dev libtiff-dev
```

On other systems it is typically just as easy, you simply have to modify "apt" to
the appropriate package manager (e.g. yum).

Once git and cmake are installed, relion can be easily installed through:

```
git clone https://github.com/3dem/relion.git
cd relion
git checkout ver3.1
mkdir build
cd build
cmake ..
make
```

The binaries will be produced in the `build/bin` directory. If you want to copy binaries
into somewhere else, run `cmake` with `-DCMAKE_INSTALL_PREFIX=/where/to/install/` and
perform `make install` as the final step. Do not specify the build directory itself
as `CMAKE_INSTALL_PREFIX`! This will not work.

Also note that the MPI library used for compilation must be the one you intend to use RELION with.
Compiling RELION with one version of MPI and running the resulting binary with mpirun from another
version can cause crash. See our wiki below for details.

In any case, you have to make sure your PATH environmental variable points to the directory
containing relion binaries. Launching RELION as `/path/to/relion` is NOT a right way; this
starts the right GUI, but the GUI might invoke other versions of RELION in the PATH.

If FLTK related errors are reported, please add `-DFORCE_OWN_FLTK=ON` to
`cmake`. For FFTW related errors, try `-DFORCE_OWN_FFTW=ON`.

RELION also requires libtiff. Most Linux distributions have packages like `libtiff-dev` or `libtiff-devel`.
Note that you need a developer package. You need version 4.0.x to read BigTIFF files. If you installed
libtiff in a non-standard location, specify the location by
`-DTIFF_INCLUDE_DIR=/path/to/include -DTIFF_LIBRARY=/path/to/libtiff.so.5`.

See [our wiki](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install) for more
options, troubleshooting and useful environmental variables (especially in HPC clusters).

## Updating

RELION is intermittently updated, with both minor and major features.
To update an existing installation, simply use the following commands

```
cd relion
git pull
cd build
make
make install # Only when you have specified CMAKE_INSTALL_PREFIX in the cmake step
```

If something went wrong, remove the `build` directory and try again from `cmake`.

## Options for accelerated versions

Parts of the cryo-EM processing pipeline can be very computationally demanding, and in some cases special
hardware can be used to make these faster. There are two such cases at the moment;

* GPU acceleration: RELION only supports CUDA-capable GPUs of compute capabilty 3.5 or higher.
* Vectorized CPU code path: RELION only supports GCC and ICC 2018.3 or later.

Note that you cannot have both acceleration in the same binary at the moment.

There are more benefits than speed; the accelearated versions also have a decreased memory footprint.
Details about how to enable either of these options is listed below.

## GPU-acceleration

Tools that are GPU-accelerated:
* relion\_refine (i.e. Class2D, Class3D, Refine3D, Multibody refinement)
* relion\_autopick

Classification without alignment is not accelerated.

When CUDA SDK is available, GPU support is automatically compiled.

### Use

If you run relion\_refine with a the "`--gpu`" flag, you will run the accelerated CUDA version of the kernels.
If you leave out the "`--gpu`" flag, it will run the original CPU version.

## CPU-acceleration

Tools that are CPU-accelerated (vectorized):
* relion\_refine (i.e. Class2D, Class3D, Refine3D, Multibody refinement)

Classification without alignment is not accelerated.

To build with support for CPU-accelerated kernels in addition to the original CPU version, build by setting `ALTCPU=ON`

```
cd build
rm -r *
cmake -DALTCPU=ON ..
make
make install
```

This will require the Intel TBB (Threading Building Blocks) library. RELION will look for TBB,
and fetch and install it when it is missing on your system. You can force this behaviour (and make sure
you are using the latest version) by adding:

```
-DFORCE_OWN_TBB=ON
```

In addition, you can make use the Intel Math Kernel Library (Intel MKL).
This is optional (but will scale better with increased threads). Add this by:
```
-DMKLFFT=ON
```

### Use

If you run relion\_refine with a the "`--cpu`" flag, you will run the accelerated version.
If you leave it the original CPU version will be run. You should use this flag if you can, unless you want to verify old runs or behaviour.

For details on how to compile with Intel compilers and optimal runtime configulations,
please look at our [wiki](https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Benchmarks_%26_computer_hardware#Accelerated_RELION.2C_using_GPUs_or_CPU-vectorization).
