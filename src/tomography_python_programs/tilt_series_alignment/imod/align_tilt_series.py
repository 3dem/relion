from pathlib import Path
from typing import Callable, Dict, Any

import mrcfile
import numpy as np
import pandas as pd
from rich.console import Console

from .._job_utils import (
    create_alignment_job_directories,
    write_single_tilt_series_alignment_output
)
from ._utils import get_xyz_extrinsic_euler_angles, get_specimen_shifts


def align_single_tilt_series(
        tilt_series_id: str,
        tilt_series_df: pd.DataFrame,
        tilt_image_df: pd.DataFrame,
        alignment_function: Callable,
        alignment_function_kwargs: Dict[str, Any],
        output_directory: Path,
):
    """Align a single tilt-series in IMOD using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series_id: 'rlnTomoName' in RELION tilt-series metadata.
    tilt_series_df: master file for tilt-series metadata.
    tilt_image_df: file containing information for images in a single tilt-series.
    alignment_function: alignment function from yet_another_imod_wrapper.
    alignment_function_kwargs: keyword arguments specific to the alignment function.
    output_directory: directory in which results will be stored.
    """
    console = Console(record=True)

    # Create output directory structure
    imod_directory, metadata_directory = \
        create_alignment_job_directories(output_directory, tilt_series_id)

    # Order is important in IMOD, sort by tilt angle
    tilt_image_df = tilt_image_df.sort_values(by='rlnTomoNominalStageTiltAngle', ascending=True)

    # Align tilt-series using IMOD
    # implicit assumption - one tilt-axis angle per tilt-series
    console.log('Running IMOD alignment')
    etomo_output = alignment_function(
        tilt_series=np.stack([mrcfile.read(f) for f in tilt_image_df['rlnMicrographName']]),
        tilt_angles=tilt_image_df['rlnTomoNominalStageTiltAngle'],
        pixel_size=tilt_series_df['rlnTomoTiltSeriesPixelSize'],
        nominal_rotation_angle=tilt_image_df['rlnTomoNominalTiltAxisAngle'][0],
        basename=tilt_series_id,
        output_directory=imod_directory,
        **alignment_function_kwargs,
    )

    # Write out alignment metadata
    console.log('Converting metadata and writing STAR file for aligned tilt-series')
    euler_angles = get_xyz_extrinsic_euler_angles(etomo_output)
    specimen_shifts_px = get_specimen_shifts(etomo_output)
    write_single_tilt_series_alignment_output(
        tilt_series_df=tilt_image_df,
        tilt_series_id=tilt_series_id,
        euler_angles=euler_angles,
        specimen_shifts=specimen_shifts_px * tilt_series_df['rlnTomoTiltSeriesPixelSize'],
        output_star_file=metadata_directory / f'{tilt_series_id}.star',
    )
