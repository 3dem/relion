RELION 5.0-beta
===============

RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

If RELION is useful in your work, please cite our papers.

Comprehensive documentation of RELION and tutorials are stored [here](https://relion.readthedocs.io/).

## Installation

See our relion.readthedocs.io page for [Installation instructions](https://relion.readthedocs.io/en/release-5.0/Installation.html),

Compared to previous release of RELION, please note that to add support for Python modules (e.g. Blush, ModelAngelo and DynaMight) you will have to install a Python environment. We recommend installing miniconda3. You can find that here: [miniconda](https://docs.conda.io/en/latest/miniconda.html) 

Also note that the MPI library used for compilation must be the one you intend to use RELION with.
Compiling RELION with one version of MPI and running the resulting binary with mpirun from another
version can cause crash. 

In any case, you have to make sure your PATH environmental variable points to the directory
containing RELION binaries. Launching RELION as `/path/to/relion` is NOT a right way; this
starts the right GUI, but the GUI might invoke other versions of RELION in the PATH.

If FLTK related errors are reported, please add `-DFORCE_OWN_FLTK=ON` to `cmake`. For FFTW related errors, try `-DFORCE_OWN_FFTW=ON`.

RELION also requires libtiff. Most Linux distributions have packages like `libtiff-dev` or `libtiff-devel`.
Note that you need a developer package. You need version 4.0.x to read BigTIFF files. If you installed
libtiff in a non-standard location, specify the location by
`-DTIFF_INCLUDE_DIR=/path/to/include -DTIFF_LIBRARY=/path/to/libtiff.so.5`.


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
