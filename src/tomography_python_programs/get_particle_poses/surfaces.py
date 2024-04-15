from pathlib import Path

import numpy as np
import pandas as pd
import starfile
import typer
from scipy.spatial.transform.rotation import Rotation as R

from ._cli import cli
from .._utils.relion import relion_pipeline_job
from .._utils.cgal import Polyhedron

COMMAND_NAME = 'surfaces'

def vector2matrix(z_vectors: np.ndarray):
    """Turn normal vector into rotation matrix
    Stolen from the PoseSampler class in MorphoSampler.
    """
    z_vectors /= np.linalg.norm(z_vectors, axis=-1, keepdims=True)
    y_vectors = np.cross(z_vectors, [1, 1, 1])
    y_vectors /= np.linalg.norm(y_vectors, axis=-1, keepdims=True)
    x_vectors = np.cross(y_vectors, z_vectors)
    x_vectors /= np.linalg.norm(x_vectors, axis=-1, keepdims=True)
    orientations = np.empty(shape=(len(z_vectors), 3, 3), dtype=np.float32)
    orientations[:, :, 0] = x_vectors
    orientations[:, :, 1] = y_vectors
    orientations[:, :, 2] = z_vectors

    # randomise in plane rotation
    angles = np.random.uniform(0, 360, size=(len(z_vectors)))
    Rz = R.from_euler('z', angles=angles, degrees=True).as_matrix()
    return orientations @ Rz



@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def derive_poses_on_surfaces(
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
    global_df = starfile.read(tilt_series_star_file)
    global_df = global_df.set_index('rlnTomoName')
    annotation_files = annotations_directory.glob('*_surfaces.npz')
    dfs = []
    rotated_basis = R.from_euler('y', angles=90, degrees=True).as_matrix()
    for file in annotation_files:
        surface_file = np.load(file)
        tilt_series_id = '_'.join(file.name.split('_')[:-1])
        pixel_size = float(global_df.loc[tilt_series_id, 'rlnTomoTiltSeriesPixelSize'])
        scale_factor = float(global_df.loc[tilt_series_id, 'rlnTomoTomogramBinning'])

        polyhedron = Polyhedron(
            surface_file['vertices'] * scale_factor,
            surface_file['indices']
        )
        polyhedron.remove_isolated_vertices()
        polyhedron.isotropic_remeshing(spacing_angstroms * pixel_size)
        vertices = polyhedron.vertices_array()
        normals = polyhedron.compute_vertex_normals()
        assert len(vertices) == len(normals)
        rotated_orientations = vector2matrix(normals[:, (2, 1, 0)]) @ rotated_basis
        eulers = R.from_matrix(rotated_orientations).inv().as_euler(
            seq='ZYZ', degrees=True,
        )
        assert len(vertices) == len(eulers)
        data = {
            'rlnTomoName': [tilt_series_id] * len(vertices),
            'rlnCoordinateX': vertices[:, 2],
            'rlnCoordinateY': vertices[:, 1],
            'rlnCoordinateZ': vertices[:, 0],
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

    df2 = pd.DataFrame({'rlnTomoParticlesFile' : [output_file],
                        'rlnTomoTomogramsFile' : [tilt_series_star_file]})
    opt_file = output_directory / 'optimisation_set.star'
    starfile.write({'optimisation_set': df2}, opt_file, overwrite=True)
