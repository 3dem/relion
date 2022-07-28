from pathlib import Path
from typing import Optional, Tuple, List

import typer
from rich.console import Console
from rich.progress import track

from .align_tilt_series import align_single_tilt_series
from .._cli import cli
from .._job_utils import write_global_output
from ... import utils
from ...utils.relion import relion_pipeline_job

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
    console.log('Extracting metadata for tilt series.')
    tilt_series_metadata = list(utils.star.iterate_tilt_series_metadata(
        tilt_series_star_file=tilt_series_star_file,
        tilt_series_id=tomogram_name
    ))
    for tilt_series_id, tilt_series_df, tilt_image_df in tilt_series_metadata:
        console.log(f'Aligning {tilt_series_id}...')
        align_single_tilt_series(
            tilt_series_id=tilt_series_id,
            global_df=tilt_series_df,
            tilt_series_df=tilt_image_df,
            sample_thickness_nanometers=sample_thickness_nanometers,
            tilt_angle_offset_correction=do_tilt_angle_offset_correction,
            gpu_ids=gpu,
            job_directory=output_directory,
        )
    if tomogram_name is None:  # write out STAR file for set of tilt-series
        console.log('Writing aligned_tilt_series.star')
        write_global_output(
            original_tilt_series_star_file=tilt_series_star_file,
            job_directory=output_directory
        )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
