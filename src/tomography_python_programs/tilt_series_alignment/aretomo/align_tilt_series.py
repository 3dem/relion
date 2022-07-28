from pathlib import Path

import pandas as pd
import numpy as np
import mrcfile

from lil_aretomo import align_tilt_series as align_tilt_series_with_aretomo
from lil_aretomo.utils import read_aln
from rich.console import Console
from typing import Optional, Tuple

from .._job_utils import (
    create_alignment_job_directories,
    write_single_tilt_series_alignment_output
)
from .utils import get_specimen_shifts, get_xyz_extrinsic_euler_angles


def align_single_tilt_series(
        tilt_series_id: str,
        global_df: pd.DataFrame,
        tilt_series_df: pd.DataFrame,
        sample_thickness_nanometers: float,
        tilt_angle_offset_correction: bool,
        gpu_ids: Optional[str],
        job_directory: Path,
):
    """Align a single tilt-series in AreTomo using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series_id: 'rlnTomoName' in RELION tilt-series metadata.
    global_df: data from global tilt-series metadata.
    tilt_series_df: file containing information for images in a single tilt-series.
    sample_thickness_nanometers: thickness of intermediate reconstruction during alignments.
    tilt_angle_offset_correction: flag to enable/disable stage tilt offset correction (-TiltCor) in AreTomo
    gpu_ids: string to specify GPUs. GPU identifiers should be separated by colons e.g. 0:1:2:3
    job_directory: directory in which results will be stored.
    """
    console = Console(record=True)

    # Create output directory structure
    aretomo_directory, metadata_directory = \
        create_alignment_job_directories(job_directory, tilt_series_id)

    # Order is important in IMOD, sort by tilt angle
    tilt_series_df = tilt_series_df.sort_values(by='rlnTomoNominalStageTiltAngle', ascending=True)

    console.log(f'Running AreTomo on {tilt_series_id}')
    aretomo_output = align_tilt_series_with_aretomo(
        tilt_series=np.stack([mrcfile.read(f) for f in tilt_series_df['rlnMicrographName']]),
        tilt_angles=tilt_series_df['rlnTomoNominalStageTiltAngle'],
        pixel_size=global_df['rlnTomoTiltSeriesPixelSize'],
        nominal_rotation_angle=tilt_series_df['rlnTomoNominalTiltAxisAngle'][0],
        basename=tilt_series_id,
        output_directory=aretomo_directory,
        sample_thickness_nanometers=sample_thickness_nanometers,
        do_local_alignments=False,
        correct_tilt_angle_offset=tilt_angle_offset_correction,
        output_pixel_size=20,
        gpu_ids=gpu_ids
    )
    console.log('Writing STAR file for aligned tilt-series')
    write_single_tilt_series_alignment_output(
        tilt_series_df=tilt_series_df,
        tilt_series_id=tilt_series_id,
        euler_angles=get_xyz_extrinsic_euler_angles(aretomo_output.aln_file),
        specimen_shifts=get_specimen_shifts(aretomo_output.aln_file) * global_df['rlnTomoTiltSeriesPixelSize'],
        output_star_file=metadata_directory / f'{tilt_series_id}.star',
    )
