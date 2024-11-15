from pathlib import Path
from typing import Callable, Dict, Any

import mrcfile
import numpy as np
from rich.console import Console

from .._job_utils import (
    create_alignment_job_directories,
)
from ._utils import get_xyz_extrinsic_euler_angles, get_specimen_shifts
from ..._metadata_models.relion.tilt_series import RlnTiltSeries


def align_single_tilt_series(
        tilt_series: RlnTiltSeries,
        pixel_spacing_angstroms: float,
        alignment_function: Callable,
        alignment_function_kwargs: Dict[str, Any],
        output_directory: Path,
):
    """Align a single tilt-series in IMOD using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series: RELION metadata for a tilt-series.
    pixel_spacing_angstroms: spacing between pixels in Angstroms.
    alignment_function: alignment function from yet_another_imod_wrapper.
    alignment_function_kwargs: keyword arguments specific to the alignment function.
    output_directory: directory in which results will be stored.
    """
    console = Console(record=True)
    error_console = Console(stderr=True, style="bold red")

    # Create output directory structure
    imod_directory, metadata_directory = \
        create_alignment_job_directories(output_directory, tilt_series.name)

    # Order is important in IMOD, sort by tilt angle
    df = tilt_series.data.sort_values(by='rlnTomoNominalStageTiltAngle', ascending=True)

    # Align tilt-series using IMOD
    # implicit assumption - one tilt-axis angle per tilt-series
    console.log('Running IMOD alignment')
    try:
        etomo_output = alignment_function(
            tilt_series=np.stack([mrcfile.read(f) for f in df['rlnMicrographName']]),
            tilt_angles=df['rlnTomoNominalStageTiltAngle'],
            pixel_size=pixel_spacing_angstroms,
            nominal_rotation_angle=df['rlnTomoNominalTiltAxisAngle'][0],
            basename=tilt_series.name,
            output_directory=imod_directory,
            **alignment_function_kwargs,
        )

        # Write out alignment metadata
        console.log('Converting metadata and writing STAR file for aligned tilt-series')
        df[['rlnTomoXTilt', 'rlnTomoYTilt', 'rlnTomoZRot']] = get_xyz_extrinsic_euler_angles(
            etomo_output)
        df[['rlnTomoXShiftAngst', 'rlnTomoYShiftAngst']] = get_specimen_shifts(
            etomo_output) * pixel_spacing_angstroms
        RlnTiltSeries(name=tilt_series.name, data=df).write_star_file(
            metadata_directory / f'{tilt_series.name}.star')
    except RuntimeError as e:
        error_console.log(f'ERROR: alignment for {tilt_series.name} failed with error:{str(e)}')
