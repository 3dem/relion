from pathlib import Path
from typing import Optional, List

import typer
from rich.console import Console

from .align_tilt_series import align_single_tilt_series
from .._cli import cli
from .._job_utils import write_global_output
from ..._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from ..._utils.relion import relion_pipeline_job

console = Console(record=True)


@cli.command(name='AreTomo')
@relion_pipeline_job
def aretomo_cli(
        tilt_series_star_file: Path = typer.Option(..., help='RELION tilt-series STAR file'),
        output_directory: Path = typer.Option(..., help='directory in which results will be stored'),
        sample_thickness_nanometers: float = typer.Option(150, help='estimated sample thickness in nanometers'),
        do_tilt_angle_offset_correction: bool = typer.Option(False, help='enable/disable stage tilt offset correction (-TiltCor in AreTomo)'),
        tomogram_name: Optional[str] = typer.Option(None, help="'rlnTomoName' for a specific tilt-series"),
        gpu: Optional[List[int]] = typer.Option(None, help="zero-indexed GPU identifiers as integers"),
):
    """Align one or multiple tilt-series in AreTomo."""
    if not tilt_series_star_file.exists():
        raise RuntimeError('Could not find tilt series star file')
    console.log('Extracting metadata from STAR file')

    tilt_series_set = RlnTiltSeriesSet.from_star_file(
        filename=tilt_series_star_file, tilt_series_id=tomogram_name
    )
    for global_data, tilt_series in tilt_series_set:
        console.log(f'Aligning {tilt_series.name}...')
        align_single_tilt_series(
            tilt_series=tilt_series,
            pixel_spacing_angstroms=global_data['rlnTomoTiltSeriesPixelSize'],
            sample_thickness_nanometers=sample_thickness_nanometers,
            tilt_angle_offset_correction=do_tilt_angle_offset_correction,
            gpu_ids=gpu,
            job_directory=output_directory,
        )
    if tomogram_name is None:  # write out STAR file for set of tilt-series
        console.log('Writing aligned_tilt_series.star')
        write_global_output(
            input_tilt_series_set=tilt_series_set,
            job_directory=output_directory
        )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
