from enum import Enum
from pathlib import Path

import napari
import numpy as np
import starfile
import magicgui
from scipy.spatial.transform import Rotation as R

from ..view._cli import cli

COMMAND_NAME = 'particles'

@cli.command(name=COMMAND_NAME, no_args_is_help=True)
def view_particles(
    particle_star_file: Path
):
    star = starfile.read(particle_star_file)

    df = star['particles']

    grouped = df.groupby('rlnTomoName')
    TiltSeriesId = Enum('TiltSeriesIds',
                        {ts_id: ts_id for ts_id in grouped.groups.keys()})
    for ts_id in TiltSeriesId:
        first_tilt_series = ts_id
        break

    widget = magicgui.widgets.create_widget(annotation=TiltSeriesId)

    viewer = napari.Viewer()
    viewer.window.add_dock_widget(widget, area='left', name='tilt-series id')


@widget.changed.connect
def load_tilt_series(tilt_series_id: TiltSeriesId):
    tilt_series_id = tilt_series_id.value

    # positions
    zyx = _get_zyx(tilt_series_id)
    if 'particle positions' not in viewer.layers:
        viewer.add_points(zyx, name='particle positions', size=40)
    else:
        viewer.layers['particle positions'].data = zyx

    # basis vectors
    z_vec, y_vec, x_vec = _get_zyx_vectors(tilt_series_id)
    z_vec = np.stack([zyx, z_vec], axis=1)
    y_vec = np.stack([zyx, y_vec], axis=1)
    x_vec = np.stack([zyx, x_vec], axis=1)
    if 'particle z' not in viewer.layers:
        viewer.add_vectors(z_vec, name='particle z', length=30, edge_color='orange')
    else:
        viewer.layers['particle z'].data = z_vec
    if 'particle y' not in viewer.layers:
        viewer.add_vectors(y_vec, name='particle y', length=10, edge_color='blue')
    else:
        viewer.layers['particle y'].data = y_vec
    if 'particle x' not in viewer.layers:
        viewer.add_vectors(x_vec, name='particle x', length=10, edge_color='purple')
    else:
        viewer.layers['particle x'].data = x_vec


def _get_zyx(tilt_series_id: str) -> np.ndarray:
    df = grouped.get_group(tilt_series_id)
    zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()
    if 'rlnOriginZAngst' in df.columns:
        shifts = df[
            ['rlnOriginZAngst', 'rlnOriginYAngst', 'rlnOriginXAngst']].to_numpy()
        zyx -= shifts
    return zyx


def _get_zyx_vectors(tilt_series_id: str) -> np.ndarray:
    df = grouped.get_group(tilt_series_id)
    eulers = df[['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']].to_numpy()
    rotation = R.from_euler(seq='ZYZ', angles=eulers, degrees=True).inv()
    rotation_matrix = rotation.as_matrix()
    if 'rlnTomoSubtomogramRot' in df.columns:
        eulers = df[['rlnTomoSubtomogramRot', 'rlnTomoSubtomogramTilt',
                     'rlnTomoSubtomogramPsi']]
        rotation = R.from_euler(seq='ZYZ', angles=eulers, degrees=True).inv()
        subtomogram_rotation_matrix = rotation.as_matrix()
        rotation_matrix = rotation_matrix @ subtomogram_rotation_matrix
    z_vec = rotation_matrix[:, :, 2]
    y_vec = rotation_matrix[:, :, 1]
    x_vec = rotation_matrix[:, :, 0]
    return z_vec, y_vec, x_vec


load_tilt_series(first_tilt_series)
napari.run()
