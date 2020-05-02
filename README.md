RELION 3.1 beta (EER branch)
=============================

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

To process EER dataset, proceed as follows.

1. Decide how many (internal) frames to group into a fraction.

   For example, if you have 1000 internal frames and group them by 30,
   you will get 33 fractions. The remaining 10 (= 1000 - 30 * 33) frames will be ignored.
2. Calculate how many e/A2 each fraction has.
3. Import EER files as usual.

   The pixel size should be **half of the physical size**, because RELION renders
   electrons in a 8K x 8K super resolution grid by default.
4. Run motion correction.
   - Specify the dose rate calculated in step 2.
   - Specify the gain reference.
   - **Grouping in the GUI must be 1** regardless of what you choose in step 1.
   - **Binning should be 2** to bring a 8K super resolution grid into a 4K physical grid by Fourier cropping.
   - Add `--eer_grouping 30` as additional arguments to specify the value decided in step 1.

If you specify `--eer_upsampling 1`, RELION renders electrons in a 4K x 4K grid from the beginning,
thus saving memory and increasing the processing speed. In this case, you should specify the physical
pixel size in step 3 and the binning should be 1 in step 4. This might reduce the data quality due with
noise beyond physical Nyquist aliasing back, but in our tests, the difference was tiny (only one
or two shells) if any.

The gain reference for EER is different from multiplicative gain references for K2/K3. When the
movie is in EER format, RELION will *divide* raw pixel values with the provided gain. When the gain
is zero, the pixel is considered as defective. With `--eer_upsamling 2` (default), the gain reference
can be 8K x 8K or 4K x 4K. In the latter case, the gain is upsampled. With `--eer_upsamling 1`, the
gain reference must be 4K x 4K.

If memory usage is a concern, consider building RELION in CPU single precision (`cmake -DDoublePrec_CPU=OFF`).

In future, we will make it possible to change `eer_upsampling` and `eer_grouping` during Polish.
This way, you can start processing at 4K and coarse slicing and then switch to 8K and finer slicing
in Polish to save processing time. Currently you have to manually modify trajectory STAR files to
do this (not officially supported).

Another useful tool is `relion_eer_to_tiff`, which renders an EER movie into a compressed integer TIFF.
Due to the different meanings of the gain reference for EER and TIFF, you have to take the inverse
of the EER gain reference yourself before processing the resulting TIFF files.

## Installation

More extensive options and configurations are available [here](http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install),
but the outlines to clone and install relion for typical use are made easy through [cmake](https://en.wikipedia.org/wiki/CMake).

On Debian or Ubuntu machines, installing cmake, the compiler, and additional dependencies (mpi, fftw) is as easy as:

```
sudo apt install cmake git build-essential mpi-default-bin mpi-default-dev libfftw3-dev
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

RELION requires libtiff to read movies in TIFF and EER. Most Linux distributions have packages like
`libtiff-dev` or `libtiff-devel`. Note that you need a developer package. You need version 4.0.x
to read BigTIFF files. If you installed libtiff in a non-standard location, specify the location by
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
