import os

import mrcfile
import numpy as np


def read_mrc(filename: os.PathLike) -> np.ndarray:
    with mrcfile.open(filename, permissive=True) as mrc:
        data = mrc.data
    return data


def get_image_dimensions(filename: os.PathLike) -> np.ndarray:
    """Get array of image dimensions (xyz) from an MRC file header."""
    with mrcfile.open(filename, header_only=True) as mrc:
        return np.array([mrc.header.nx, mrc.header.ny, mrc.header.nz])
