from pathlib import Path

import pandas as pd
import numpy as np
import mrcfile 

from lil_aretomo import align_tilt_series_with_aretomo
from rich.console import Console
from typing import Optional, Tuple

from ._utils import coerce_gpu_ids

from .._job_utils import (
    create_alignment_job_directory_structure,
    write_single_tilt_series_alignment_output
)
from ... import utils

def align_single_tilt_series(
        tilt_series_id: str,
        tilt_series_df: pd.DataFrame,
        tilt_image_df: pd.DataFrame,
        do_local_alignments: bool,
        output_resolution: float,
        n_patches_xy: Tuple[int, int],
        alignment_thickness_px: float,
        tilt_angle_offset_correction: bool,
        gpu_ids: Optional[str],
        job_directory: Path,
):
    """Align a single tilt-series in AreTomo using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series_id: 'rlnTomoName' in RELION tilt-series metadata.
    tilt_series_df: master file for tilt-series metadata.
    tilt_image_df: file containing information for images in a single tilt-series.
    do_local_alignments: flag to enable local alignments.
    output_resolution: resolution for output tilt series in angstroms. Has no bearing on final resolution of tomogram.
    n_patches_xy: number of patches in x and y for local alignments
    alignment_thickness_px: thickness of intermediate reconstruction during alignments.
    tilt_angle_offset_correction: flag to enable/disable stage tilt offset correction (-TiltCor) in AreTomo
    gpu_ids: string to specify GPUs. GPU identifiers should be separated by colons e.g. 0:1:2:3
    job_directory: directory in which results will be stored.
    """
    console = Console(record=True)

    # Create output directory structure
    external_directory, metadata_directory = \
        create_alignment_job_directory_structure(job_directory)
    aretomo_directory = external_directory / tilt_series_id
    aretomo_directory.mkdir(parents=True, exist_ok=True)

    # Establish filenames
    tilt_series_filename = f'{tilt_series_id}.mrc'
    tilt_image_metadata_filename = f'{tilt_series_id}.star'

    # Order is important in IMOD, sort by tilt angle
    tilt_image_df = tilt_image_df.sort_values(by='rlnTomoNominalStageTiltAngle', ascending=True)
    
    if gpu_ids is not None:
        gpu_ids = coerce_gpu_ids(
            gpu_ids=gpu_ids
        )
    console.log('Running AreTomo')
    align_tilt_series_with_aretomo(
        tilt_series=np.stack([mrcfile.read(f) for f in tilt_image_df['rlnMicrographName']]),
        tilt_angles=tilt_image_df['rlnTomoNominalStageTiltAngle'],
        pixel_size=tilt_series_df['rlnTomoTiltSeriesPixelSize'],
        basename=tilt_series_id,
        output_directory=aretomo_directory,
        expected_sample_thickness=alignment_thickness_px,
        do_local_alignments=do_local_alignments,
        output_pixel_size=output_resolution / 2,
        nominal_rotation_angle=tilt_image_df['rlnTomoNominalTiltAxisAngle'][0],
        n_patches_xy=n_patches_xy,
        correct_tilt_angle_offset=tilt_angle_offset_correction,
        gpu_ids=gpu_ids
    )
    console.log('Writing STAR file for aligned tilt-series')
    write_single_tilt_series_alignment_output(
        tilt_image_df=tilt_image_df,
        tilt_series_id=tilt_series_id,
        pixel_size=tilt_series_df['rlnTomoTiltSeriesPixelSize'],
        alignment_directory=aretomo_directory,
        output_star_file=metadata_directory / tilt_image_metadata_filename,
    )
