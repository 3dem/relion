import numpy as np

from tomography_preprocessing.utils.transformations import Rx, Ry, Rz, S


def test_single_x_rotation():
    """Rotation around X should be right handed.
    """
    xyzw = np.array([1, 1, 1, 1]).reshape((4, 1))
    result = Rx(90) @ xyzw
    expected = np.array([1, -1, 1, 1]).reshape((4, 1))
    assert result.shape == (4, 1)
    assert np.allclose(result, expected)


def test_single_y_rotation():
    """Rotation around Y should be right handed.
    """
    xyzw = np.array([1, 1, 1, 1]).reshape((4, 1))
    result = Ry(90) @ xyzw
    expected = np.array([1, 1, -1, 1]).reshape((4, 1))
    assert result.shape == (4, 1)
    assert np.allclose(result, expected)


def test_single_z_rotation():
    """Rotation around Z should be right handed.
    """
    xyzw = np.array([1, 1, 1, 1]).reshape((4, 1))
    result = Rz(90) @ xyzw
    expected = np.array([-1, 1, 1, 1]).reshape((4, 1))
    assert result.shape == (4, 1)
    assert np.allclose(result, expected)


def test_multiple_x_rotations():
    """Should be able to generate multiple matrices from an array of angles."""
    angles = np.linspace(0, 90, 10)
    matrices = Rx(angles)
    assert matrices.shape == (10, 4, 4)
    assert np.allclose(matrices[0], np.eye(4))
    assert np.allclose(matrices[-1], Rx(90))


def test_multiple_y_rotations():
    """Should be able to generate multiple matrices from an array of angles."""
    angles = np.linspace(0, 90, 10)
    matrices = Ry(angles)
    assert matrices.shape == (10, 4, 4)
    assert np.allclose(matrices[0], np.eye(4))
    assert np.allclose(matrices[-1], Ry(90))


def test_multiple_z_rotations():
    """Should be able to generate multiple matrices from an array of angles."""
    angles = np.linspace(0, 90, 10)
    matrices = Rz(angles)
    assert matrices.shape == (10, 4, 4)
    assert np.allclose(matrices[0], np.eye(4))
    assert np.allclose(matrices[-1], Rz(90))


def test_single_shift_matrix_3d():
    """Generate a single (4, 4) affine matrix from a 3D shift vector."""
    shifts = [1, 2, 3]
    matrix = S(shifts)
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix[:3, 3], shifts)
    assert np.allclose(matrix[:3, :3], np.eye(3))


def test_single_shift_matrix_2d():
    """Generate a single (4, 4) affine matrix from a 2D shift vector."""
    shifts = [1, 2]
    matrix = S(shifts)
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix[:2, 3], shifts)
    assert matrix[2, 3] == 0
    assert np.allclose(matrix[:3, :3], np.eye(3))


def test_multiple_shift_matrices_3d():
    """Generate multiple (4, 4) affine matrices from an array of 3D shift vectors."""
    shifts = np.arange(18).reshape((6, 3))
    matrices = S(shifts)
    assert matrices.shape == (6, 4, 4)
    assert np.allclose(matrices[:, :3, 3], shifts)
    assert np.allclose(matrices[:, :3, :3], np.eye(3))


def test_multiple_shift_matrices_3d():
    """Generate multiple (4, 4) affine matrices from an array of 2D shift vectors."""
    shifts = np.arange(18).reshape((9, 2))
    matrices = S(shifts)
    assert matrices.shape == (9, 4, 4)
    assert np.allclose(matrices[:, :2, 3], shifts)
    assert np.allclose(matrices[:, 2, 3], 0)
    assert np.allclose(matrices[:, :3, :3], np.eye(3))
