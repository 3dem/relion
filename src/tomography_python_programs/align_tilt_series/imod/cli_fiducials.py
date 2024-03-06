from pathlib import Path
from typing import Optional

import typer
from yet_another_imod_wrapper import align_tilt_series_using_fiducials
from rich.console import Console

from .._job_utils import write_global_output
from .align_tilt_series import align_single_tilt_series
from .._cli import cli
from ..._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from ..._utils.relion import relion_pipeline_job

console = Console(record=True)


@cli.command(name='IMOD:fiducials')
@relion_pipeline_job
def fiducials_cli(
        tilt_series_star_file: Path = typer.Option(..., help='RELION tilt-series STAR file'),
        output_directory: Path = typer.Option(..., help='directory in which results will be stored'),
        nominal_fiducial_diameter_nanometers: float = typer.Option(..., help='nominal fiducial diameter in nanometers'),
        tomogram_name: Optional[str] = typer.Option(None, help="'rlnTomoName' in RELION tilt-series metadata")
):
    """Align one or multiple tilt-series with fiducials in IMOD."""
    if not tilt_series_star_file.exists():
        raise RuntimeError('Could not find tilt series star file')

    console.log('Extracting metadata for tilt series.')
    tilt_series_set = RlnTiltSeriesSet.from_star_file(
        filename=tilt_series_star_file, tilt_series_id=tomogram_name
    )
    for global_data, tilt_series in tilt_series_set:
        aligned_tilt_series_star_file = output_directory / 'tilt_series' / f'{tilt_series.name}.star'
        if aligned_tilt_series_star_file.exists():
            console.log(f'tilt_series/{tilt_series.name}.star already exists, skipping...')
            continue

        console.log(f'Aligning {tilt_series.name}...')
        align_single_tilt_series(
            tilt_series=tilt_series,
            pixel_spacing_angstroms=global_data['rlnTomoTiltSeriesPixelSize'],
            alignment_function=align_tilt_series_using_fiducials,
            alignment_function_kwargs={
                'fiducial_size': nominal_fiducial_diameter_nanometers,
                'skip_if_completed': True
            },
            output_directory=output_directory,
        )
    if tomogram_name is None:  # write out STAR file for set of tilt-series
        console.log('Writing aligned_tilt_series.star')
        write_global_output(
            input_tilt_series_set=tilt_series_set,
            job_directory=output_directory
        )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
