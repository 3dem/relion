from pathlib import Path

import pandas as pd
import starfile
import typer
from morphosamplers import Sphere, sphere_samplers
from scipy.spatial.transform.rotation import Rotation as R

from ._cli import cli
from .._utils.relion import relion_pipeline_job

COMMAND_NAME = 'spheres'


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def derive_poses_on_spheres(
    tilt_series_star_file: Path = typer.Option(
        ..., help='tilt-series STAR file containing tomogram'
    ),
    annotations_directory: Path = typer.Option(
        ..., help='directory containing annotations in each tomogram'
    ),
    output_directory: Path = typer.Option(
        ..., help="directory into which 'particles.star' will be written."
    ),
    spacing_angstroms: float = typer.Option(
        ..., help="target spacing between particle poses on spheres in angstroms."
    ),
):
    global_df = starfile.read(tilt_series_star_file, parse_as_string=['rlnTomoName'])
    global_df = global_df.set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_spheres.star')
    dfs = []
    for file in annotation_files:
        sphere_df = starfile.read(file, parse_as_string=['rlnTomoName'])
        tilt_series_id = '_'.join(file.name.split('_')[:-1])
        pixel_size = float(
            global_df.loc[tilt_series_id, 'rlnTomoTiltSeriesPixelSize'])

        tomo_size = global_df.loc[tilt_series_id][
            ['rlnTomoSizeX', 'rlnTomoSizeY', 'rlnTomoSizeZ']
        ].to_numpy().astype(float)
        tomo_center = (tomo_size / 2 - 1) * pixel_size

        scale_factor = pixel_size * \
            float(global_df.loc[tilt_series_id, 'rlnTomoTomogramBinning'])
        centers = sphere_df[['rlnCoordinateX',
                             'rlnCoordinateY',
                             'rlnCoordinateZ']]
        centers = centers.to_numpy() * scale_factor - tomo_center
        radii = sphere_df['rlnSphereRadius'].to_numpy() * scale_factor

        for center, radius in zip(centers, radii):
            sphere = Sphere(center=center, radius=radius)
            sampler = sphere_samplers.PoseSampler(spacing=spacing_angstroms)
            poses = sampler.sample(sphere)

            # rot/psi are locked when tilt==0
            # rotate particles 90 degrees around the Y-axis so that tilt ~=90
            # during refinement
            rotated_basis = R.from_euler(
                'y', angles=90, degrees=True).as_matrix()
            rotated_orientations = poses.orientations @ rotated_basis
            eulers = R.from_matrix(rotated_orientations).inv().as_euler(
                seq='ZYZ', degrees=True,
            )

            data = {
                'rlnTomoName': [tilt_series_id] * len(poses),
                'rlnCenteredCoordinateXAngst': poses.positions[:, 0],
                'rlnCenteredCoordinateYAngst': poses.positions[:, 1],
                'rlnCenteredCoordinateZAngst': poses.positions[:, 2],
                'rlnTomoSubtomogramRot': eulers[:, 0],
                'rlnTomoSubtomogramTilt': eulers[:, 1],
                'rlnTomoSubtomogramPsi': eulers[:, 2],
            }
            dfs.append(pd.DataFrame(data))
    df = pd.concat(dfs)
    rot_prior, tilt_prior, psi_prior = R.from_matrix(rotated_basis).as_euler(
        seq='ZYZ', degrees=True
    )
    df['rlnAngleRot'] = [rot_prior] * len(df)
    df['rlnAngleTilt'] = [tilt_prior] * len(df)
    df['rlnAnglePsi'] = [psi_prior] * len(df)
    df['rlnAngleTiltPrior'] = [tilt_prior] * len(df)
    df['rlnAnglePsiPrior'] = [psi_prior] * len(df)
    output_file = output_directory / 'particles.star'
    starfile.write({'particles': df}, output_file, overwrite=True)

    df2 = pd.DataFrame({'rlnTomoParticlesFile': [output_file],
                        'rlnTomoTomogramsFile': [tilt_series_star_file]})
    opt_file = output_directory / 'optimisation_set.star'
    starfile.write({'optimisation_set': df2}, opt_file, overwrite=True)
