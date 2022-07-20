from pathlib import Path
from typing import Tuple, Optional

import numpy as np
import pandas as pd
import starfile
import typer

from .imod._utils import get_tilt_series_alignment_data, calculate_specimen_shifts, \
    get_xyz_extrinsic_euler_angles
from ..utils.transformations import S, Rx, Ry, Rz
from ._cli import cli


def create_alignment_job_directory_structure(output_directory: Path) -> Tuple[Path, Path, Path]:
    """Create directory structure for a tilt-series alignment job."""
    stacks_directory = output_directory / 'stacks'
    stacks_directory.mkdir(parents=True, exist_ok=True)

    external_directory = output_directory / 'external'
    external_directory.mkdir(parents=True, exist_ok=True)

    metadata_directory = output_directory / 'tilt_series'
    metadata_directory.mkdir(parents=True, exist_ok=True)
    return stacks_directory, external_directory, metadata_directory


def tilt_series_alignment_parameters_to_relion_projection_matrices(
        specimen_shifts: pd.DataFrame,
        euler_angles: pd.DataFrame,
        tilt_image_dimensions: np.ndarray,
        tomogram_dimensions: np.ndarray,
):
    """Generate affine matrices transforming points in 3D to 2D in tilt-images.

    Projection model:
    3D specimen is rotated about its center then translated such that the projection
    of points onto the XY-plane gives their position in a tilt-image.

    More specifically
    - 3D specimen is rotated about its center by
        - shifting the origin to the specimen center
        - rotated extrinsically about the Y-axis by the tilt angle
        - rotated extrinsically about the Z-axis by the in plane rotation angle
    - 3D specimen is translated to align coordinate system with tilt-image
        - move center-of-rotation of specimen to center of tilt-image
        - move center-of-rotation of specimen to rotation center in tilt-image

    Parameters
    ----------
    specimen_shifts: XY-shifts which align the projected specimen with tilt-images
    euler_angles: YZX intrinsic Euler angles which transform the specimen
    tilt_image_dimensions: XY-dimensions of tilt-series.
    tomogram_dimensions: size of tomogram in XYZ
    """
    tilt_image_center = tilt_image_dimensions / 2
    specimen_center = tomogram_dimensions / 2

    # Transformations, defined in order of application
    s0 = S(-specimen_center)  # put specimen center-of-rotation at the origin
    r0 = Rx(euler_angles['rlnTomoXTilt'])  # rotate specimen around X-axis
    r1 = Ry(euler_angles['rlnTomoYTilt'])  # rotate specimen around Y-axis
    r2 = Rz(euler_angles['rlnTomoZRot'])  # rotate specimen around Z-axis
    s1 = S(specimen_shifts)  # shift projected specimen in xy (camera) plane
    s2 = S(tilt_image_center)  # move specimen back into tilt-image coordinate system

    # compose matrices
    transformations = s2 @ s1 @ r2 @ r1 @ r0 @ s0
    return np.squeeze(transformations)


def write_single_tilt_series_alignment_output(
        tilt_image_df: pd.DataFrame,
        tilt_series_id: str,
        pixel_size: float,
        alignment_directory: Path,
        output_star_file: Path,
):
    """Write metadata from a tilt-series alignment experiment."""
    xf, tlt = get_tilt_series_alignment_data(alignment_directory, tilt_series_id)

    shifts_px = calculate_specimen_shifts(xf)
    tilt_image_df[['rlnTomoXShiftAngst', 'rlnTomoYShiftAngst']] = shifts_px * pixel_size

    euler_angles = get_xyz_extrinsic_euler_angles(xf, tlt)
    tilt_image_df[['rlnTomoXTilt', 'rlnTomoYTilt', 'rlnTomoZRot']] = euler_angles

    starfile.write({tilt_series_id: tilt_image_df}, output_star_file)


@cli.command(name='write-global-output')
def write_aligned_tilt_series_star_file(
        original_tilt_series_star_file: Path = typer.Option(...),
        job_directory: Path = typer.Option(...)
):
    """Write output from a batch of tilt-series alignment experiments."""
    df = starfile.read(original_tilt_series_star_file, always_dict=True)['global']

    # update individual tilt series star files
    df['rlnTomoTiltSeriesStarFile'] = [
        job_directory / 'tilt_series' / f'{tilt_series_id}.star'
        for tilt_series_id in df['rlnTomoName']
    ]
    df['rlnEtomoDirectiveFile'] = [
        job_directory / 'external' / tilt_series_id / f'{tilt_series_id}.edf'
        for tilt_series_id in df['rlnTomoName']
    ]
    etomo_directives_exist = any(df['rlnEtomoDirectiveFile'].apply(lambda x: Path(x).exists()))
    if etomo_directives_exist is False:
        df = df.drop(columns=['rlnEtomoDirectiveFile'])

    # check which output files were succesfully generated, take only those
    df = df[df['rlnTomoTiltSeriesStarFile'].apply(lambda x: x.exists())]
    starfile.write({'global': df}, job_directory / 'aligned_tilt_series.star')
