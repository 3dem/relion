from pathlib import Path

import pandas as pd
import starfile
import typer
from rich.console import Console

from ._cli import cli
from .._utils.relion import relion_pipeline_job

console = Console(record=True)

@cli.command(name='combine-particles', no_args_is_help=True)
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

    df2 = pd.DataFrame({'rlnTomoParticlesFile' : [output_file],
                        'rlnTomoTomogramsFile' : [tilt_series_star_file]})
    opt_file = output_directory / 'optimisation_set.star'
    starfile.write({'optimisation_set': df2}, opt_file, overwrite=True)


@cli.command(name='split-particles', no_args_is_help=True)
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
    optics_df = star_data['optics']
    particles_df = star_data['particles']
    tomo_names = particles_df.rlnTomoName.unique()
    console.log(f'  Tomograms found in the input star file: {tomo_names}.')

    for tomo_name in tomo_names:
        anno_file_name = f'{tomo_name}_particles.star'
        anno_file = annotations_directory / anno_file_name
    
        if anno_file.exists():
            console.log(f'  {anno_file_name} already exists, moving on.')
            continue 
    
        # TODO: take 'rlnOriginXAngst' into account when it exists

        tomo_df = particles_df.loc[particles_df['rlnTomoName'] == tomo_name]
        tomo_bin = tomo_data.rlnTomoTomogramBinning[
                tomo_data.rlnTomoName == tomo_name
        ]
        assert(len(tomo_bin) == 1)
        tomo_bin = tomo_bin.iloc[0]
        console.log(f'    tomo_bin = {tomo_bin}')

        # TODO: first check if rlnOriginXAngst column exists 
        tilt_series_pixel_size = optics_df.rlnTomoTiltSeriesPixelSize[
                optics_df['rlnOpticsGroupName'] == tomo_name
        ]
        assert(len(tilt_series_pixel_size) == 1)
        tilt_series_pixel_size = tilt_series_pixel_size.iloc[0]

        anno_df = pd.DataFrame({
            # TODO: when doing the full thing,
            # the division by tomo_bin comes last
            'rlnTomoName'    : tomo_df['rlnTomoName'],
            'rlnCoordinateX' : tomo_df['rlnCoordinateX'] / tomo_bin,
            'rlnCoordinateY' : tomo_df['rlnCoordinateY'] / tomo_bin,
            'rlnCoordinateZ' : tomo_df['rlnCoordinateZ'] / tomo_bin
        })

        starfile.write(anno_df, anno_file)

        console.log(f'  Wrote {anno_file_name}')
