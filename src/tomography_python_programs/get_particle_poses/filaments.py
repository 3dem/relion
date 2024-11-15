import numpy as np
import pandas as pd
import starfile
import typer
import pathlib
from morphosamplers import Path, path_samplers
from scipy.spatial.transform.rotation import Rotation as R

from ._cli import cli
from .._utils.relion import relion_pipeline_job

COMMAND_NAME = 'filaments'


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def get_poses_along_filament_backbones(
    tilt_series_star_file: pathlib.Path = typer.Option(
        ..., help='tilt-series STAR file containing tomogram'
    ),
    annotations_directory: pathlib.Path = typer.Option(
        ..., help='directory containing annotations in each tomogram'
    ),
    output_directory: pathlib.Path = typer.Option(
        ..., help="directory into which 'particles.star' will be written."
    ),
    spacing_angstroms: float = typer.Option(
        ..., help="spacing between particles along filaments in angstroms."
    ),
    filament_polarity_known: bool = typer.Option(
        True, help="Whether filament polarity from annotations should be fixed "
                   "during refinement."
    )
):
    global_df = starfile.read(tilt_series_star_file,
                              parse_as_string=['rlnTomoName'])
    global_df = global_df.set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_filaments.star')
    dfs = []
    for file in annotation_files:
        filament_df = starfile.read(file, parse_as_string=['rlnTomoName'])
        tilt_series_id = '_'.join(file.name.split('_')[:-1])
        pixel_size = float(
            global_df.loc[tilt_series_id, 'rlnTomoTiltSeriesPixelSize']
        )
        scale_factor = pixel_size * \
            float(global_df.loc[tilt_series_id, 'rlnTomoTomogramBinning'])
        tomo_size = global_df.loc[tilt_series_id][
            ['rlnTomoSizeX', 'rlnTomoSizeY', 'rlnTomoSizeZ']
        ].to_numpy().astype(float)
        tomo_center = (tomo_size / 2 - 1) * pixel_size
        for filament_id, df in filament_df.groupby('rlnTomoManifoldIndex'):
            xyz = df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']]
            xyz = xyz.to_numpy() * scale_factor - tomo_center

            # If two consecutive points on the path are identical, remove one of them
            # or scipy interpolate.splprep throws error.
            # This can happen due to user or napari error.
            indices_to_delete = []
            for i in range(xyz.shape[0]-1):
                if np.linalg.norm(xyz[i]-xyz[i+1]) < 1e-6:
                    indices_to_delete.append(i+1)
            xyz = np.delete(xyz, indices_to_delete, axis=0)
            # derive equidistant poses along length of filament
            path = Path(control_points=xyz)
            pose_sampler = path_samplers.HelicalPoseSampler(
                spacing=spacing_angstroms, twist=0
            )
            poses = pose_sampler.sample(path)

            # rot/psi are coupled when tilt==0,
            # pre-rotate particles such that they have tilt=90 relative to a reference
            # filament aligned along the z-axis
            rotated_basis = R.from_euler(
                'y', angles=90, degrees=True).as_matrix()
            rotated_orientations = poses.orientations @ rotated_basis
            rotated_eulers = R.from_matrix(rotated_orientations).inv().as_euler(
                seq='ZYZ', degrees=True,
            )

            # calculate total length of filament in pixels at tilt series pixel size
            point_sampler = path_samplers.PointSampler(spacing=0.1)
            points = point_sampler.sample(path)
            differences = np.diff(points, axis=0)
            distances = np.linalg.norm(differences, axis=1)
            total_length = np.sum(distances)

            # how far along the helix is each particle? in angstroms
            total_length = total_length / pixel_size
            distance_along_helix = np.linspace(
                0, 1, num=len(poses)) * total_length

            data = {
                'rlnTomoName': [tilt_series_id] * len(poses),
                'rlnHelicalTubeID': [int(filament_id) + 1] * len(poses),
                'rlnHelicalTrackLengthAngst': distance_along_helix,
                'rlnCenteredCoordinateXAngst': poses.positions[:, 0],
                'rlnCenteredCoordinateYAngst': poses.positions[:, 1],
                'rlnCenteredCoordinateZAngst': poses.positions[:, 2],
                'rlnTomoSubtomogramRot': rotated_eulers[:, 0],
                'rlnTomoSubtomogramTilt': rotated_eulers[:, 1],
                'rlnTomoSubtomogramPsi': rotated_eulers[:, 2],
            }
            dfs.append(pd.DataFrame(data))
    df = pd.concat(dfs)

    # add priors on orientations
    rot_prior, tilt_prior, psi_prior = R.from_matrix(rotated_basis).as_euler(
        seq='ZYZ', degrees=True
    )
    df['rlnAngleRot'] = [rot_prior] * len(df)
    df['rlnAngleTilt'] = [tilt_prior] * len(df)
    df['rlnAnglePsi'] = [psi_prior] * len(df)
    df['rlnAngleTiltPrior'] = [tilt_prior] * len(df)
    df['rlnAnglePsiPrior'] = [psi_prior] * len(df)

    if filament_polarity_known is False:
        df['rlnAnglePsiFlipRatio'] = [0.5] * len(df)

    # write output
    output_file = output_directory / 'particles.star'
    starfile.write({'particles': df}, output_file, overwrite=True)

    df2 = pd.DataFrame({'rlnTomoParticlesFile': [output_file],
                        'rlnTomoTomogramsFile': [tilt_series_star_file]})
    opt_file = output_directory / 'optimisation_set.star'
    starfile.write({'optimisation_set': df2}, opt_file, overwrite=True)
