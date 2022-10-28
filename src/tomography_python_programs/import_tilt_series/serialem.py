from pathlib import Path
from typing import Optional, List

import mdocfile
import numpy as np
import pandas as pd
import rich
import starfile
import typer
from rich.progress import track

from ._cli import cli
from .. import utils
from ..utils.mdoc import calculate_pre_exposure_dose, basename_from_sub_frame_path
from ..utils.relion import relion_pipeline_job
from ..utils.file import match_filenames


console = rich.console.Console(record=True)


@cli.command(name='SerialEM')
@relion_pipeline_job
def import_tilt_series_from_serial_em(
        tilt_image_movie_pattern: str = typer.Option(...),
        mdoc_file_pattern: str = typer.Option(...),
        output_directory: Path = typer.Option(...),
        nominal_tilt_axis_angle: float = typer.Option(...),
        nominal_pixel_size: float = typer.Option(...),
        voltage: float = typer.Option(...),
        spherical_aberration: float = typer.Option(...),
        amplitude_contrast: float = typer.Option(...),
        invert_defocus_handedness: Optional[bool] = typer.Option(None),
        images_are_motion_corrected: bool = False,
        dose_per_tilt_image: Optional[float] = None,
        dose_per_movie_frame: Optional[float] = None,
        prefix: str = '',
        mtf_file: Optional[Path] = None,
) -> Path:
    """Import tilt-series data using SerialEM metadata.

    Parameters
    ----------
    tilt_image_movie_pattern : Filename pattern with wildcard characters
        for tilt image movie files.
    mdoc_file_pattern : Filename pattern wildcard characters for mdoc files
        containing tilt-series metadata.
    output_directory : Directory in which output files are written.
    nominal_tilt_axis_angle : Nominal tilt axis angle in the images.
    nominal_pixel_size : Nominal pixel spacing of the tilt-images in angstroms.
    voltage : Acceleration voltage (keV)
    spherical_aberration : Spherical aberration (mm)
    amplitude_contrast : Amplitude contrast fraction (e.g. 0.1)
    invert_defocus_handedness: Set this to flip the handedness of the defocus geometry (default=1).
        The value of this parameter is either +1 or -1, and it describes whether the focus 
        increases or decreases as a function of Z distance. It has to be determined experimentally. 
    images_are_motion_corrected : Set this if your raw images have already been motion corrected and/or are not movies.
    dose_per_tilt_image : dose applied in electrons per square angstrom to each tilt image.
        If set, this will override the values from the mdoc file.
        Only one of `dose_per_tilt_image` and `dose_per_movie_frame` can be set.
    dose_per_movie_frame : dose in electrons per square angstrom in each movie to each movie frame.
        If set, this will override the values from the mdoc file.
        Only one of `dose_per_tilt_image` and `dose_per_movie_frame` can be set.
    prefix : a prefix which will be added to the tilt-series name.
        This avoids name collisions when processing data from multiple collection
        sessions.
    mtf_file: a text file containing the MTF of the detector.

    Returns
    -------
    tilt_series_star_file : a text file pointing to STAR files containing per tilt-series metadata
    """
    console.log('Started import_tilt_series job.')

    # Create output directory structure
    tilt_series_directory = Path(output_directory) / 'tilt_series'
    console.log(f'Creating output directory at {tilt_series_directory}')
    tilt_series_directory.mkdir(parents=True, exist_ok=True)

    # Get lists of necessary files
    tilt_image_files = list(Path().glob(tilt_image_movie_pattern))
    mdoc_files = sorted(list(Path().glob(mdoc_file_pattern)))
    console.log(f'Found {len(mdoc_files)} mdoc files and {len(tilt_image_files)} image files.')

    # Raise error if no tilt images or mdocs found
    if len(tilt_image_files) == 0:
        raise RuntimeError('Could not find any image files')
    if len(mdoc_files) == 0:
        raise RuntimeError('Could not find any mdoc files')

    # Get tomogram ids and construct paths for per-tilt-series STAR files
    tomogram_ids = [
        utils.mdoc.construct_tomogram_id(mdoc_file, prefix)
        for mdoc_file in mdoc_files
    ]
    tilt_series_star_files = [
        tilt_series_directory / f"{tomogram_id}.star"
        for tomogram_id in tomogram_ids
    ]

    # Construct a global dataframe pointing at per-tilt-series STAR files
    global_star_file = Path(output_directory) / 'tilt_series.star'
    global_df = pd.DataFrame({
        'rlnTomoName': tomogram_ids,
        'rlnTomoTiltSeriesStarFile': tilt_series_star_files,
    })
    global_df['rlnVoltage'] = voltage
    global_df['rlnSphericalAberration'] = spherical_aberration
    global_df['rlnAmplitudeContrast'] = amplitude_contrast
    global_df['rlnMicrographOriginalPixelSize'] = nominal_pixel_size
    if images_are_motion_corrected:
        global_df['rlnTomoTiltSeriesPixelSize'] = nominal_pixel_size
    global_df['rlnTomoHand'] = -1 if invert_defocus_handedness else 1
    if mtf_file is not None:
        global_df['rlnMtfFileName'] = mtf_file

    starfile.write(
        data={'global': global_df},
        filename=global_star_file,
        overwrite=True
    )
    console.log(f'Wrote tilt-series data STAR file {global_star_file}')

    # Write out per tilt-series STAR files
    console.log('writing per tilt-series STAR files...')
    for tomogram_id, mdoc_file, output_filename in \
            track(list(zip(tomogram_ids, mdoc_files, tilt_series_star_files))):
        tilt_image_df = _generate_tilt_image_dataframe(
            mdoc_file=mdoc_file,
            tilt_image_files=tilt_image_files,
            dose_per_tilt_image=dose_per_tilt_image,
            dose_per_movie_frame=dose_per_movie_frame,
            nominal_tilt_axis_angle=nominal_tilt_axis_angle,
            images_are_motion_corrected=images_are_motion_corrected
        )
        starfile.write(
            data={f'{tomogram_id}': tilt_image_df},
            filename=output_filename,
            overwrite=True
        )
    console.log(f'Wrote STAR files for {len(tilt_series_star_files)} tilt-series.')
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
    return global_star_file


def _generate_tilt_image_dataframe(
        mdoc_file: Path,
        tilt_image_files: List[Path],
        nominal_tilt_axis_angle: float,
        dose_per_tilt_image: Optional[float],
        dose_per_movie_frame: Optional[float],
        images_are_motion_corrected: Optional[bool],
) -> pd.DataFrame:
    """Generate a dataframe containing data about images in a tilt-series."""
    df = mdocfile.read(mdoc_file, camel_to_snake=True)
    df['date_time'] = pd.to_datetime(df['date_time'], infer_datetime_format=True)
    df = df.sort_values(by="date_time", ascending=True)
    df['image_index'] = np.arange(len(df)) + 1  # 0 -> 1 indexing
    df['pre_exposure_dose'] = calculate_pre_exposure_dose(
        df, dose_per_tilt_image=dose_per_tilt_image, dose_per_movie_frame=dose_per_movie_frame,
    )
    df['sub_frame_basename'] = df['sub_frame_path'].apply(basename_from_sub_frame_path)
    df['tilt_image_file'] = match_filenames(df['sub_frame_basename'], to_match=tilt_image_files)
    df = df[df['tilt_image_file'] != None]
    df['nominal_tilt_axis_angle'] = nominal_tilt_axis_angle

    output_df = pd.DataFrame({
        'rlnTomoNominalStageTiltAngle': df['tilt_angle'],
        'rlnTomoNominalTiltAxisAngle': df['nominal_tilt_axis_angle'],
        'rlnMicrographPreExposure': df['pre_exposure_dose'],
    })
    
    if not images_are_motion_corrected:
        movie_information_df=pd.DataFrame({
            'rlnMicrographMovieName': df['tilt_image_file'],
            'rlnTomoTiltMovieFrameCount': df['num_sub_frames'],
            'rlnTomoTiltMovieIndex': df['image_index']
        })
        output_df=pd.concat([movie_information_df, output_df], axis=1)
    if 'target_defocus' in df.columns:
        output_df['rlnTomoNominalDefocus'] = df['target_defocus']
    if images_are_motion_corrected:
        output_df['rlnMicrographName'] = df['tilt_image_file']
    return output_df
