from pathlib import Path
from typing import Optional

import typer
from yet_another_imod_wrapper import align_tilt_series_using_fiducials
from rich.console import Console

from .._job_utils import write_aligned_tilt_series_star_file
from .align_tilt_series import align_single_tilt_series
from .._cli import cli
from ... import utils
from ...utils.relion import relion_pipeline_job

console = Console(record=True)


@cli.command(name='IMOD:fiducials')
@relion_pipeline_job
def batch_fiducials(
        tilt_series_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        tomogram_name: Optional[str] = typer.Option(None),
        nominal_fiducial_diameter_nanometres: float = typer.Option(...),
):
    """Align one or multiple tilt-series with fiducials in IMOD.

    Parameters
    ----------
    tilt_series_star_file: RELION tilt-series STAR file.
    output_directory: directory in which results will be stored.
    tomogram_name: 'rlnTomoName' in RELION tilt-series metadata.
    nominal_fiducial_diameter_nanometres: nominal fiducial diameter in nanometers.
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
            alignment_function=align_tilt_series_using_fiducials,
            alignment_function_kwargs={
                'fiducial_size': nominal_fiducial_diameter_nanometres
            },
            output_directory=output_directory,
        )
    if tomogram_name is None:  # write out STAR file for set of tilt-series
        console.log('Writing aligned_tilt_series.star')
        write_aligned_tilt_series_star_file(
            original_tilt_series_star_file=tilt_series_star_file,
            job_directory=output_directory
        )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
