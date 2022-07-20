from pathlib import Path

import starfile

from .._cli import cli
from .._job_utils import create_alignment_job_directory_structure
from ...utils.relion import relion_pipeline_job


@cli.command(name='IMOD:etomo-launcher')
@relion_pipeline_job
def etomo_launcher(
        tilt_series_star_file: Path,
        output_directory: Path,
):
    tilt_series_directory, imod_alignments_directory = \
        create_alignment_job_directory_structure(output_directory)
    tilt_series_df = starfile.read(tilt_series_star_file, always_dict=True)['global']


class EtomoLauncher:
    def __init__(
            self, tilt_series_star_file: Path, output_directory: Path
    ):
        self.tilt_series_directory, self.imod_alignments_directory = \
            create_alignment_job_directory_structure(output_directory)
        self.tilt_series_df = starfile.read(tilt_series_star_file, always_dict=True)['global']
