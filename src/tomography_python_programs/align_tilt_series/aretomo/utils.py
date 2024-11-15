from pathlib import Path

import numpy as np

from lil_aretomo.utils import read_aln


def get_specimen_shifts(aln_file: Path) -> np.ndarray:
    """Get specimen shifts from AreTomo alignments file."""
    df = read_aln(aln_file)
    return np.array(df[['TX', 'TY']])


def get_xyz_extrinsic_euler_angles(aln_file: Path) -> np.ndarray:
    """Get XYZ-extrinsic Euler angles froom AreTomo alignments file."""
    df = read_aln(aln_file)
    n_images = len(df)
    euler_angles = np.empty(shape=(n_images, 3))
    euler_angles[:, 0] = 0
    euler_angles[:, 1] = df['TILT']
    euler_angles[:, 2] = df['ROT']
    return euler_angles
