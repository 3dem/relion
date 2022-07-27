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
        tilt_series_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        sample_thickness_nanometers: float = typer.Option(150),
        n_patches_xy: Optional[Tuple[int, int]] = typer.Option((None, None)),
        tilt_angle_offset_correction: bool = typer.Option(False),
        tomogram_name: Optional[str] = typer.Option(None),
        gpu_ids: Optional[List[int]] = typer.Option(None),
):
    """Align one or multiple tilt-series in AreTomo using RELION tilt-series metadata.

    Parameters
    ----------
    tilt_series_star_file: RELION tilt-series STAR file.
    output_directory: directory in which results will be stored.
    n_patches_xy: number of patches in x and y used in local alignments.
    sample_thickness_nanometers: thickness of intermediate reconstructions in angstroms.
    tomogram_name: 'rlnTomoName' for a specific tilt-series (optional).
    tilt_angle_offset_correction: enable/disable stage tilt offset correction (-TiltCor) in AreTomo.
    gpu_ids: zero-indexed GPU identifiers as integers.
    """
    if not tilt_series_star_file.exists():
        e = 'Could not find tilt series star file'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)
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
            n_patches_xy=n_patches_xy,
            sample_thickness_nanometers=sample_thickness_nanometers,
            tilt_angle_offset_correction=tilt_angle_offset_correction,
            gpu_ids=gpu_ids,
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
