from pathlib import Path
from typing import Optional

import typer
from yet_another_imod_wrapper import align_tilt_series_using_patch_tracking
from rich.console import Console
from rich.progress import track

from .align_tilt_series import align_single_tilt_series
from .._cli import cli
from .._job_utils import write_global_output
from ... import utils
from ...utils.relion import relion_pipeline_job

console = Console(record=True)


@cli.command(name='IMOD:patch-tracking')
@relion_pipeline_job
def patch_tracking_cli(
        tilt_series_star_file: Path = typer.Option(..., help='RELION tilt-series STAR file'),
        output_directory: Path = typer.Option(..., help='directory in which to store results'),
        patch_size_nanometers: float = typer.Option(...),
        patch_overlap_percentage: float = typer.Option(...),
        tomogram_name: Optional[str] = typer.Option(None),
):
    """Align one or multiple tilt-series with patch-tracking in IMOD.

    Parameters
    ----------
    tilt_series_star_file: RELION tilt-series STAR file.
    output_directory: directory in which to store results.
    tomogram_name: 'rlnTomoName' in tilt-series STAR file.
    patch_size_nanometers: size of 2D patches used for alignment.
    patch_overlap_percentage: percentage of overlap between tracked patches.
    """
    if not tilt_series_star_file.exists():
        e = 'Could not find tilt series star file'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)
    console.log('Extracting metadata for tilt series.')
    tilt_series_metadata = utils.star.iterate_tilt_series_metadata(
        tilt_series_star_file=tilt_series_star_file,
        tilt_series_id=tomogram_name
    )
    for tilt_series_id, tilt_series_df, tilt_image_df in tilt_series_metadata:
        console.log(f'Aligning {tilt_series_id}...')
        align_single_tilt_series(
            tilt_series_id=tilt_series_id,
            tilt_series_df=tilt_series_df,
            tilt_image_df=tilt_image_df,
            alignment_function=align_tilt_series_using_patch_tracking,
            alignment_function_kwargs={
                'patch_size': int(patch_size_nanometers * 10),
                'patch_overlap_percentage': patch_overlap_percentage,
            },
            output_directory=output_directory,
        )
    if tomogram_name is None:  # write global output for set of tilt-series
        console.log('Writing aligned_tilt_series.star')
        write_global_output(
            original_tilt_series_star_file=tilt_series_star_file,
            job_directory=output_directory
        )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
