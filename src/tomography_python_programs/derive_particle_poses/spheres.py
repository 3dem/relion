from pathlib import Path

import pandas as pd
import starfile
import typer
from morphosamplers import Sphere, sphere_samplers
from scipy.spatial.transform.rotation import Rotation as R

from ._cli import cli

COMMAND_NAME = 'spheres'


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
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
    star = starfile.read(tilt_series_star_file)
    global_table = star['global'].set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_spheres.star')
    dfs = []
    for file in annotation_files:
        sphere_df = starfile.read(file)
        tilt_series_id = ''.join(file.name.split('_')[:-1])
        pixel_size = float(global_table[tilt_series_id]['rlnTomoTiltSeriesPixelSize'])
        scale_factor = float(global_table[tilt_series_id]['rlnTomoTomogramBinning'])
        centers = sphere_df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']]
        centers = centers.to_numpy() * scale_factor
        radii = sphere_df['rlnSphereRadius'].to_numpy() * scale_factor
        for center, radius in zip(centers, radii):
            sphere = Sphere(center=center, radius=radius)
            sampler = sphere_samplers.PoseSampler(spacing=spacing_angstroms / pixel_size)
            poses = sampler.sample(sphere)
            eulers = R.from_matrix(poses.orientations).inv().as_euler(
                seq='ZYZ', degrees=True,
            )
            data = {
                'rlnTomoName': [tilt_series_id] * len(poses),
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
