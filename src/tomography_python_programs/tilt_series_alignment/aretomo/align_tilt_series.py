from pathlib import Path

import numpy as np
import mrcfile

from lil_aretomo import align_tilt_series as align_tilt_series_with_aretomo
from rich.console import Console
from typing import Optional, List

from .._job_utils import create_alignment_job_directories
from .utils import get_specimen_shifts, get_xyz_extrinsic_euler_angles
from ...metadata_model import RelionTiltSeries


def align_single_tilt_series(
        tilt_series: RelionTiltSeries,
        pixel_spacing_angstroms: float,
        sample_thickness_nanometers: float,
        tilt_angle_offset_correction: bool,
        gpu_ids: Optional[List[int]],
        job_directory: Path,
):
    """Align a single tilt-series in AreTomo using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series: RELION tilt-series metadata.
    pixel_spacing_angstroms: pixel spacing in Angstroms.
    sample_thickness_nanometers: thickness of intermediate reconstruction during alignments.
    tilt_angle_offset_correction: flag to enable/disable stage tilt offset correction (-TiltCor) in AreTomo
    gpu_ids: List of integers to specify zero-indexed GPU IDs. 
    job_directory: directory in which results will be stored.
    """
    console = Console(record=True)

    # Create output directory structure
    aretomo_directory, metadata_directory = \
        create_alignment_job_directories(job_directory, tilt_series.name)

    # Order is important in IMOD, sort by tilt angle
    df = tilt_series.data.sort_values(by='rlnTomoNominalStageTiltAngle', ascending=True)

    console.log(f'Running AreTomo on {tilt_series.name}')
    aretomo_output = align_tilt_series_with_aretomo(
        tilt_series=np.stack([mrcfile.read(f) for f in df['rlnMicrographName']]),
        tilt_angles=df['rlnTomoNominalStageTiltAngle'],
        pixel_size=pixel_spacing_angstroms,
        nominal_rotation_angle=df['rlnTomoNominalTiltAxisAngle'][0],
        basename=tilt_series.name,
        output_directory=aretomo_directory,
        sample_thickness_nanometers=sample_thickness_nanometers,
        do_local_alignments=False,
        correct_tilt_angle_offset=tilt_angle_offset_correction,
        output_pixel_size=20,
        gpu_ids=gpu_ids
    )

    console.log('Writing STAR file for aligned tilt-series')
    df[['rlnTomoXTilt', 'rlnTomoYTilt', 'rlnTomoZRot']] = \
        get_xyz_extrinsic_euler_angles(aretomo_output.aln_file)
    df[['rlnTomoXShiftAngst', 'rlnTomoYShiftAngst']] = \
        get_specimen_shifts(aretomo_output.aln_file) * pixel_spacing_angstroms
    RelionTiltSeries(name=tilt_series.name, data=df).write_star_file(
        metadata_directory / f'{tilt_series.name}.star')
