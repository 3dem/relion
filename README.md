RELION 5.0.0
============

RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

If RELION is useful in your work, please cite our papers.

Comprehensive documentation of RELION and tutorials are stored [here](https://relion.readthedocs.io/).

## Installation

See our [installation instructions](https://relion.readthedocs.io/en/release-5.0/Installation.html).

You will have to set up a Python environment to use Python modules (e.g. Blush, ModelAngelo and DynaMight).
Thus, please read the above instructions carefully even if you are familiar with earlier versions.

## Class Ranker

The default model for the class ranker has been trained and tested in Python 3.9.12 with Pytorch 1.10.0 and Numpy 1.20.0.
If you wish to retrain the class ranker model with your own data, please refer to [this repo](https://github.com/3dem/relion-classranker).
