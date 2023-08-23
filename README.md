RELION 5.0-alpha
============

RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

The underlying theory of MAP refinement is given in a [scientific publication](https://www.ncbi.nlm.nih.gov/pubmed/22100448).
If RELION is useful in your work, please cite this paper.

The more comprehensive documentation of RELION is stored [here](https://relion.readthedocs.io/en/release-4.0/).

## Installation

More extensive options and configurations are available [here](https://relion.readthedocs.io/en/release-4.0/Installation.html),
but the outlines to clone and install RELION for typical use are made easy through [cmake](https://en.wikipedia.org/wiki/CMake).

On Debian or Ubuntu machines, installing cmake, the compiler, and additional dependencies (mpi, fftw) is as easy as:

```
sudo apt install cmake git build-essential mpi-default-bin mpi-default-dev libfftw3-dev libtiff-dev libpng-dev ghostscript libxft-dev
```

RedHat-like systems (CentOS, RHEL, Scientific Linux etc) use `yum` package manager:

```
sudo yum install cmake git gcc gcc-c++ openmpi-devel fftw-devel libtiff-devel libpng-devel ghostscript libXft-devel libX11-devel
```

To add support for Python modules (e.g. Blush, ModelAngelo and DynaMight) you will have to install a Python environment.
We recommend installing miniconda3. You can find that here: [miniconda](https://docs.conda.io/en/latest/miniconda.html)

Once you have Conda setup, you can install all the RELION Python dependencies into a new environment by running:
```conda env create -f environment.yml```

<details>
<summary><b>Using a Custom Conda/Python Environment with RELION (Advanced)</b></summary>
To enforce RELION to utilize a particular Python interpreter, incorporate the following flag during the cmake call: 

```-DPYTHON_EXE_PATH=path/to/python```

Additionally, if you intend to make use of automatically downloaded pretrained model weights (used in e.g. Blush, ModelAngelo and Classranker),
it's recommended to set the TORCH_HOME directory. To do this, include the following flag:

```-DTORCH_HOME_PATH=path/to/torch/home```
</details>

Before installing RELION make sure Conda is not active. If it is, you can deactivate it by running:
```conda deactivate```

RELION can now be installed through:

```
git clone https://github.com/3dem/relion.git
cd relion
git checkout master # or ver4.0; see below
mkdir build
cd build
cmake ..
make
```

By performing `git checkout ver5.0` instead of `git checkout master`, you can access the latest
(developmental) updates for RELION 5.0.x. The code there is not tested as throughfully as that in
the master branch and not generally recommended.

The binaries will be produced in the `build/bin` directory. If you want to copy binaries
into somewhere else, run `cmake` with `-DCMAKE_INSTALL_PREFIX=/where/to/install/` and
perform `make install` as the final step. Do not specify the build directory itself
as `CMAKE_INSTALL_PREFIX`! This will not work.

Also note that the MPI library used for compilation must be the one you intend to use RELION with.
Compiling RELION with one version of MPI and running the resulting binary with mpirun from another
version can cause crash. See our wiki below for details.

In any case, you have to make sure your PATH environmental variable points to the directory
containing RELION binaries. Launching RELION as `/path/to/relion` is NOT a right way; this
starts the right GUI, but the GUI might invoke other versions of RELION in the PATH.

If FLTK related errors are reported, please add `-DFORCE_OWN_FLTK=ON` to
`cmake`. For FFTW related errors, try `-DFORCE_OWN_FFTW=ON`.

RELION also requires libtiff. Most Linux distributions have packages like `libtiff-dev` or `libtiff-devel`.
Note that you need a developer package. You need version 4.0.x to read BigTIFF files. If you installed
libtiff in a non-standard location, specify the location by
`-DTIFF_INCLUDE_DIR=/path/to/include -DTIFF_LIBRARY=/path/to/libtiff.so.5`.

## Building with HIP

Currently the AMD GPU ROCm port and HIP optimisation of RELION are housed in an internal repo.
Once the pre-requistes are loaded (follow the steps in previous sections), the configuration with HIP can be easily done through:
```
git clone git@github.com:AMD-HPC/RELION.git relion
cd relion
git checkout suyash/ver4.0-hip
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/where/to/install/ \
      -DCMAKE_BUILD_TYPE=release   \
      -DHIP=on -DHIP_ARCH="gfx90a" \
      -DGUI=off                    \
      -DTIFF_INCLUDE_DIR=/path/to/include  \
      -DTIFF_LIBRARY=/path/to/libtiff.so.5 \
      -DAMDFFTW=on .. # only on AMD systems to build an optimized version of FFTW lib
make -j
```

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


## Class Ranker
The default model for the class ranker has been trained and tested in Python 3.9.12 with Pytorch 1.10.0 and Numpy 1.20.0.
If you wish to retrain the class ranker model with your own data, please refer to [this repo](https://github.com/3dem/relion-classranker).
