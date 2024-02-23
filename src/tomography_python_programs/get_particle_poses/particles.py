from pathlib import Path

import pandas as pd
import numpy as np
import starfile
import typer
from rich.console import Console
from scipy.spatial.transform import Rotation

from ._cli import cli
from .._utils.relion import relion_pipeline_job

console = Console(record=True)

@cli.command(name='particles', no_args_is_help=True)
@relion_pipeline_job
def combine_particle_annotations(
    tilt_series_star_file: Path = typer.Option(
        ..., help='tilt-series STAR file containing tomogram'
    ),
    annotations_directory: Path = typer.Option(
        ..., help='directory containing annotations in each tomogram'
    ),
    output_directory: Path = typer.Option(
        ..., help="directory into which 'particles.star' will be written."
    )
):
    console.log("Running get_particle_poses combine-particles.") 

    global_df = starfile.read(tilt_series_star_file)
    global_df = global_df.set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_particles.star')
    dfs = []
    for file in annotation_files:
        df = starfile.read(file)
        tilt_series_id = '_'.join(file.name.split('_')[:-1])
        scale_factor = float(global_df.loc[tilt_series_id, 'rlnTomoTomogramBinning'])
        xyz = df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']]
        xyz = xyz.to_numpy() * scale_factor
        df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']] = xyz
        dfs.append(df)
    df = pd.concat(dfs)
    output_file = output_directory / 'particles.star'
    starfile.write({'particles': df}, output_file, overwrite=True)
    console.log(f'  Wrote {output_file}')

    df2 = pd.DataFrame({'rlnTomoParticlesFile' : [output_file],
                        'rlnTomoTomogramsFile' : [tilt_series_star_file]})
    opt_file = output_directory / 'optimisation_set.star'
    starfile.write({'optimisation_set': df2}, opt_file, overwrite=True)
    console.log(f'  Wrote {opt_file}')


@cli.command(name='particles-from-star', no_args_is_help=True)
@relion_pipeline_job
def create_annotations_from_previous_star_file(
    tomograms_file: Path = typer.Option(
        ..., help='STAR file with tomograms data'
    ),
    annotations_directory: Path = typer.Option(
        ..., help='directory containing annotations in each tomogram' 
    ),
    in_star_file: Path = typer.Option(
        None, help='STAR file with particles to annotate on the tomogram'
    )
):
    console.log("Running get_particle_poses split-particles.") 

    annotations_directory.mkdir(parents=True, exist_ok=True)

    star_data = starfile.read(in_star_file)
    tomo_data = starfile.read(tomograms_file)

    if isinstance(star_data, pd.DataFrame):
        particles_df = star_data
    else:
        particles_df = star_data['particles']

    tomo_names = particles_df.rlnTomoName.unique()
    console.log(f'  Tomograms found in the input star file: {tomo_names}.')

    if 'rlnOriginXAngst' not in particles_df.columns:
        particles_df['rlnOriginXAngst'] = 0.0
    if 'rlnOriginYAngst' not in particles_df.columns:
        particles_df['rlnOriginYAngst'] = 0.0
    if 'rlnOriginZAngst' not in particles_df.columns:
        particles_df['rlnOriginZAngst'] = 0.0

    if 'rlnTomoSubtomogramRot' not in particles_df.columns:
        particles_df['rlnTomoSubtomogramRot'] = 0.0
    if 'rlnTomoSubtomogramTilt' not in particles_df.columns:
        particles_df['rlnTomoSubtomogramTilt'] = 0.0
    if 'rlnTomoSubtomogramPsi' not in particles_df.columns:
        particles_df['rlnTomoSubtomogramPsi'] = 0.0

    for tomo_name in tomo_names:
        anno_file_name = f'{tomo_name}_particles.star'
        anno_file = annotations_directory / anno_file_name
    
        if anno_file.exists():
            console.log(f'  {anno_file_name} already exists, moving on.')
            continue 
    
        tomo_df = particles_df.loc[particles_df['rlnTomoName'] == tomo_name]
        tomo_bin = tomo_data.rlnTomoTomogramBinning[
                tomo_data.rlnTomoName == tomo_name
        ]
        tomo_x = tomo_data.rlnTomoSizeX[tomo_data.rlnTomoName == tomo_name]
        tomo_y = tomo_data.rlnTomoSizeY[tomo_data.rlnTomoName == tomo_name]
        tomo_z = tomo_data.rlnTomoSizeZ[tomo_data.rlnTomoName == tomo_name]
        assert(len(tomo_bin) == 1)
        assert(len(tomo_x) == 1)
        assert(len(tomo_y) == 1)
        assert(len(tomo_z) == 1)
        tomo_bin = tomo_bin.iloc[0]
        tomo_x = tomo_x.iloc[0]
        tomo_y = tomo_y.iloc[0]
        tomo_z = tomo_z.iloc[0]
        shift = np.array([tomo_x, tomo_y, tomo_z]) / 2

        tilt_series_pixel_size = tomo_data.rlnTomoTiltSeriesPixelSize[
                tomo_data['rlnTomoName'] == tomo_name
        ]
        assert(len(tilt_series_pixel_size) == 1)
        pixel_size = tilt_series_pixel_size.iloc[0]

        new_coords = tomo_df.apply(
            lambda df_row : rlnOrigin_to_rlnCenteredCoordAngst_row(df_row, pixel_size), 
            axis = 1
        )

        new_coords = (new_coords + shift) / tomo_bin
        new_coords.insert(loc=0, column='rlnTomoName', value=tomo_name)
        new_coords.rename(columns={
            'rlnCenteredCoordinateXAngst' : 'rlnCoordinateX',
            'rlnCenteredCoordinateYAngst' : 'rlnCoordinateY',
            'rlnCenteredCoordinateZAngst' : 'rlnCoordinateZ'
        }, inplace=True)
    
        starfile.write(new_coords, anno_file)
        console.log(f'  Wrote {anno_file_name}')


def rlnOrigin_to_rlnCenteredCoordAngst_row(df_row, pixel_size):
    """Given a row of the particles DataFrame and the pixel size of the 
    corresponding tilt series image, incorporate rlnOriginX/Y/ZAgst
    coordinates into rlnCoordinateX/Y/Z."""

    angles_subtomo = df_row[[
        'rlnTomoSubtomogramRot', 
        'rlnTomoSubtomogramTilt',
        'rlnTomoSubtomogramPsi'
    ]]

    A_subtomo = Rotation.from_euler(
            seq='ZYZ', angles=angles_subtomo, degrees=True
    ).as_matrix()

    coords = df_row[['rlnCenteredCoordinateXAngst',
            'rlnCenteredCoordinateYAngst',
            'rlnCenteredCoordinateZAngst']]

    offset = df_row[['rlnOriginXAngst', 'rlnOriginYAngst', 'rlnOriginZAngst']]

    return coords / pixel_size - A_subtomo.T @ offset / pixel_size
