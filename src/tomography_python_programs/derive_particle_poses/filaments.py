import einops
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
def derive_poses_along_filament_backbones(
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
    twist_degrees: float = typer.Option(
        ..., help="twist between particles in degrees."
    )
):
    global_df = starfile.read(tilt_series_star_file)
    global_df = global_df.set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_filaments.star')
    dfs = []
    for file in annotation_files:
        filament_df = starfile.read(file)
        tilt_series_id = '_'.join(file.name.split('_')[:-1])
        pixel_size = float(global_df.loc[tilt_series_id, 'rlnTomoTiltSeriesPixelSize'])
        scale_factor = float(global_df.loc[tilt_series_id, 'rlnTomoTomogramBinning'])
        for filament_id, df in filament_df.groupby('rlnTomoManifoldIndex'):
            xyz = df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']]
            xyz = xyz.to_numpy() * scale_factor
            path = Path(control_points=xyz)

            # calculate total length of filament in pixels at tilt series pixel size
            point_sampler = path_samplers.PointSampler(spacing=0.1)
            points = point_sampler.sample(path)
            differences = np.diff(points, axis=0)
            distances = np.linalg.norm(differences, axis=1)
            total_length = np.sum(distances)

            # derive equidistant poses along length
            pose_sampler = path_samplers.HelicalPoseSampler(
                spacing=spacing_angstroms / pixel_size, twist=twist_degrees
            )
            poses = pose_sampler.sample(path)
            eulers = R.from_matrix(poses.orientations).inv().as_euler(
                seq='ZYZ', degrees=True,
            )

            # how far along the helix is each particle? in angstroms
            total_length = total_length / pixel_size
            distance_along_helix = np.linspace(0, 1, num=len(poses)) * total_length

            data = {
                'rlnTomoName': [tilt_series_id] * len(poses),
                'rlnHelicalTubeID': [int(filament_id)] * len(poses),
                'rlnHelicalTrackLengthAngst': distance_along_helix,
                'rlnCoordinateX': poses.positions[:, 0],
                'rlnCoordinateY': poses.positions[:, 1],
                'rlnCoordinateZ': poses.positions[:, 2],
                'rlnAngleRot': eulers[:, 0],
                'rlnAngleTilt': eulers[:, 1],
                'rlnAnglePsi': eulers[:, 2],
            }
            dfs.append(pd.DataFrame(data))
    df = pd.concat(dfs)
    output_file = output_directory / 'particles.star'
    starfile.write({'particles': df}, output_file, overwrite=True)
