import numpy as np

from yet_another_imod_wrapper.utils.etomo import \
    EtomoOutput, get_tilt_angle_offset, get_input_tilt_axis_rotation_angle
from yet_another_imod_wrapper.utils.xf import XF
from yet_another_imod_wrapper.utils.io import read_tlt


def get_specimen_shifts(etomo_output: EtomoOutput) -> np.ndarray:
    """Get specimen shifts (2D, in pixels) from etomo output.

    These 2D shifts are applied after rotating a specimen to align the projection
    of the specimen with the experimental image.
    """
    xf = XF.from_file(etomo_output.xf_file)
    return xf.specimen_shifts


def get_xyz_extrinsic_euler_angles(etomo_output: EtomoOutput) -> np.ndarray:
    """Get XYZ extrinsic Euler angles which rotate the specimen from etomo."""
    in_plane_rotation = get_input_tilt_axis_rotation_angle(
        etomo_output.align_log_file
    )
    xf = XF.from_file(etomo_output.xf_file, initial_tilt_axis_rotation_angle=in_plane_rotation)
    tilt_angles = read_tlt(etomo_output.tlt_file)
    n_images = len(tilt_angles)
    tilt_angle_offset = get_tilt_angle_offset(etomo_output.align_log_file)
    if tilt_angle_offset is not None:
        tilt_angles -= tilt_angle_offset
    euler_angles = np.empty(shape=(n_images, 3))
    euler_angles[:, 0] = 0
    euler_angles[:, 1] = tilt_angles
    euler_angles[:, 2] = xf.in_plane_rotations
    return euler_angles
