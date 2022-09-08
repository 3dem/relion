import einops
import numpy as np


def Rx(angles_degrees: np.ndarray) -> np.ndarray:
    """Affine matrix for a rotation around the X-axis."""
    angles_degrees = np.asarray(angles_degrees).reshape(-1)
    c = np.cos(np.deg2rad(angles_degrees))
    s = np.sin(np.deg2rad(angles_degrees))
    matrices = einops.repeat(
        np.eye(4), 'i j -> n i j', n=len(angles_degrees)
    )
    matrices[:, 1, 1] = c
    matrices[:, 1, 2] = -s
    matrices[:, 2, 1] = s
    matrices[:, 2, 2] = c
    return np.squeeze(matrices)


def Ry(angles_degrees: np.ndarray) -> np.ndarray:
    """Affine matrix for a rotation around the Y-axis."""
    angles_degrees = np.asarray(angles_degrees).reshape(-1)
    c = np.cos(np.deg2rad(angles_degrees))
    s = np.sin(np.deg2rad(angles_degrees))
    matrices = einops.repeat(
        np.eye(4), 'i j -> n i j', n=len(angles_degrees)
    )
    matrices[:, 0, 0] = c
    matrices[:, 0, 2] = s
    matrices[:, 2, 0] = -s
    matrices[:, 2, 2] = c
    return np.squeeze(matrices)


def Rz(angles_degrees: float) -> np.ndarray:
    """Affine matrix for a rotation around the Z-axis."""
    angle_degrees = np.asarray(angles_degrees).reshape(-1)
    c = np.cos(np.deg2rad(angle_degrees))
    s = np.sin(np.deg2rad(angle_degrees))
    matrices = einops.repeat(
        np.eye(4), 'i j -> n i j', n=len(angle_degrees)
    )
    matrices[:, 0, 0] = c
    matrices[:, 0, 1] = -s
    matrices[:, 1, 0] = s
    matrices[:, 1, 1] = c
    return np.squeeze(matrices)


def S(shifts: np.ndarray) -> np.ndarray:
    """Affine matrices for shifts.

    Shifts supplied can be 2D or 3D.
    """
    shifts = np.asarray(shifts, dtype=float)
    if shifts.shape[-1] == 2:
        shifts = promote_2d_to_3d(shifts)
    shifts = np.array(shifts).reshape((-1, 3))
    matrices = einops.repeat(np.eye(4), 'i j -> n i j', n=shifts.shape[0])
    matrices[:, 0:3, 3] = shifts
    return np.squeeze(matrices)


def promote_2d_to_3d(shifts: np.ndarray) -> np.ndarray:
    """Promote 2D vectors to 3D with zeros in the last dimension."""
    shifts = np.asarray(shifts).reshape(-1, 2)
    shifts = np.c_[shifts, np.zeros(shifts.shape[0])]
    return np.squeeze(shifts)