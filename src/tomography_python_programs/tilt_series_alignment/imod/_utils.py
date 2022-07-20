import os
from functools import lru_cache
from pathlib import Path
from typing import Tuple

import numpy as np


def read_xf(file: os.PathLike) -> np.ndarray:
    """Read an IMOD xf file into an (n, 6) numpy array.

    The xf file with alignment transforms contains one
    line per view, each with a linear transformation specified by six numbers:
        A11 A12 A21 A22 DX DY
    where the coordinate (X, Y) is transformed to (X', Y') by:
        X' = A11 * X + A12 * Y + DX
        Y' = A21 * X + A22 * Y + DY
    """
    return np.loadtxt(fname=file, dtype=float).reshape((-1, 6))


def read_tlt(file: os.PathLike) -> np.ndarray:
    """Read an IMOD tlt file into an (n, ) numpy array."""
    return np.loadtxt(fname=file, dtype=float).reshape(-1)


def get_xf_shifts(xf: np.ndarray) -> np.ndarray:
    """Extract XY shifts from an IMOD xf file.

    Output is an (n, 2) numpy array of XY shifts. Shifts in an xf file are
    applied after rotations. IMOD xf files contain linear transformations. In
    the context of tilt-series alignment they contain transformations which are
    applied to 'align' a tilt-series such that images represent a fixed body rotating
    around the Y-axis.
    """
    return xf[:, -2:]


def get_xf_transformation_matrices(xf: np.ndarray) -> np.ndarray:
    """Extract the 2D transformation matrix from IMOD xf data.

    Output is an (n, 2, 2) numpy array of matrices.
    """
    return xf[:, :4].reshape((-1, 2, 2))


def get_in_plane_rotations(xf: np.ndarray) -> np.ndarray:
    """Extract the in plane rotation angle from IMOD xf data.

    Output is an (n, ) numpy array of angles in degrees. This function assumes
    that the transformation in the xf file is a simple 2D rotation.
    """
    transformation_matrices = get_xf_transformation_matrices(xf)
    cos_theta = transformation_matrices[:, 0, 0]
    return np.rad2deg(np.arccos(cos_theta))


def calculate_specimen_shifts(xf: np.ndarray) -> np.ndarray:
    """Extract specimen shifts from IMOD xf data.

    Output is an (n, 2) numpy array of XY shifts. Specimen shifts are shifts in
    the camera plane applied to the projected image of a specimen to align it
    with a tilt-image.

    This function relies on the fact that:
    RSP = S'RP where S' = R_invS

    R is a rotation
    R_inv is the inverse rotation
    S is a shift vector
    P is a point
    S' is a rotated shift vector

    In words: Shifting then rotating is equivalent to rotating then shifting by a
    rotated shift-vector.

    In IMOD, the rotation center is at (N - 1) / 2. In RELION, the rotation center
    is at N / 2. An extra half-pixel shift is applied to account for these
    differences.
    """
    rotation_matrices = get_xf_transformation_matrices(xf)
    inverse_rotation_matrices = rotation_matrices.transpose((0, 2, 1))

    image_shifts = xf[:, -2:].reshape((-1, 2, 1))
    specimen_shifts = inverse_rotation_matrices @ -image_shifts
    return specimen_shifts.reshape((-1, 2)) - 0.5


@lru_cache(maxsize=100)
def get_tilt_series_id(alignment_directory: Path) -> str:
    """Get the tilt-series stack basename from an Etomo directory."""
    edf_files = list(alignment_directory.glob('*.edf'))
    if len(edf_files) != 1:
        raise RuntimeError('Multiple Etomo directive files found')
    return edf_files[0].stem


def get_edf_file(imod_directory: Path, tilt_series_id: str) -> Path:
    """Get the Etomo directive file from an IMOD directory."""
    return imod_directory / f'{tilt_series_id}.edf'


def get_xf_file(alignment_directory: Path, tilt_series_id: str) -> Path:
    """Get the xf file containing image transforms from an IMOD directory."""
    return alignment_directory / f'{tilt_series_id}.xf'


def get_tlt_file(alignment_directory: Path, tilt_series_id: str) -> Path:
    """Get the tlt file containing tilt angles from an IMOD directory."""
    return alignment_directory / f'{tilt_series_id}.tlt'


def get_tilt_series_alignment_data(
        alignment_directory: Path, tilt_series_id: str
) -> Tuple[np.ndarray, np.ndarray]:
    """Get IMOD tilt-series alignment data."""
    xf = read_xf(get_xf_file(alignment_directory, tilt_series_id))
    tlt = read_tlt(get_tlt_file(alignment_directory, tilt_series_id))
    return xf, tlt


def get_xyz_extrinsic_euler_angles(
        imod_xf_data: np.ndarray, imod_tlt_data: np.ndarray
) -> np.ndarray:
    """Get XYZ extrinsic Euler angles which rotate the specimen."""
    euler_angles = np.zeros(shape=(len(imod_xf_data), 3))
    euler_angles[:, 1] = imod_tlt_data
    euler_angles[:, 2] = get_in_plane_rotations(imod_xf_data)
    return euler_angles


