from enum import Enum
from functools import partial
from pathlib import Path
from typing import Optional

import mrcfile
import napari
import numpy as np
import pandas as pd
import starfile
import magicgui
import typer
from pandas.core.groupby import DataFrameGroupBy
from scipy.spatial.transform import Rotation as R

from ..view._cli import cli

COMMAND_NAME = 'particles'


def _get_zyx(
    particle_df: pd.DataFrame,
) -> np.ndarray:
    zyx = particle_df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()

    # get particle shifts if present
    origin_headings = [
        'rlnOriginZAngst',
        'rlnOriginYAngst',
        'rlnOriginXAngst',
    ]
    if all(heading in particle_df.columns for heading in origin_headings):
        origin_zyx = particle_df[origin_headings].to_numpy()
        origin_in_dataframe = True
    else:
        origin_in_dataframe = False

    # get the subtomogram orientation within the tomogram, if present
    subtomogram_orientation_headings = [
        'rlnTomoSubtomogramRot',
        'rlnTomoSubtomogramTilt',
        'rlnTomoSubtomogramPsi'
    ]
    if all(heading in particle_df.columns for heading in subtomogram_orientation_headings):
        eulers = particle_df[subtomogram_orientation_headings]
        subtomogram_rotation = R.from_euler(seq='ZYZ', angles=eulers, degrees=True).inv()
        subtomogram_orientation_in_dataframe = True
    else:
        subtomogram_orientation_in_dataframe = False

    # rotate the shifts if necessary
    if subtomogram_orientation_in_dataframe and origin_in_dataframe:
        origin_xyz = origin_zyx[:, ::-1]
        origin_xyz = subtomogram_rotation.apply(origin_xyz)
        origin_zyx = origin_xyz[:, ::-1]

    if origin_in_dataframe:
        zyx -= origin_zyx / particle_df['rlnImagePixelSize'].to_numpy()
    return zyx


def _get_vectors_zyx(
    particle_df: pd.DataFrame
) -> np.ndarray:
    # get particle orientation relative to subtomograms
    eulers = particle_df[['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']].to_numpy()
    particle_rotation = R.from_euler(seq='ZYZ', angles=eulers, degrees=True).inv()

    # get the subtomogram orientation within the tomogram, if present
    subtomogram_orientation_headings = [
        'rlnTomoSubtomogramRot',
        'rlnTomoSubtomogramTilt',
        'rlnTomoSubtomogramPsi'
    ]
    if all(heading in particle_df.columns for heading in subtomogram_orientation_headings):
        eulers = particle_df[subtomogram_orientation_headings]
        subtomogram_rotation = R.from_euler(seq='ZYZ', angles=eulers, degrees=True).inv()
    else:
        subtomogram_rotation = R.identity()
    subtomogram_matrix = subtomogram_rotation.as_matrix()
    particle_matrix = particle_rotation.as_matrix()

    # compose rotations and extract vectors in zyx coords for visualisation
    rotation_matrix = subtomogram_matrix @ particle_matrix
    z_vec = rotation_matrix[:, ::-1, 2]
    y_vec = rotation_matrix[:, ::-1, 1]
    x_vec = rotation_matrix[:, ::-1, 0]
    return z_vec, y_vec, x_vec


def _load_particles(
    tilt_series_id,  # TiltSeriesId, the type is created dynamically at runtime
    viewer: napari.viewer.Viewer,
    particle_df_grouped: DataFrameGroupBy,
):
    # get particles in current tilt-series
    df = particle_df_grouped.get_group(tilt_series_id.value)

    # get positions and add to viewer
    zyx = _get_zyx(df)
    if 'particle positions' not in viewer.layers:
        viewer.add_points(
            zyx, name='particle positions', size=40, n_dimensional=True, opacity=0.5
        )
    else:
        viewer.layers['particle positions'].data = zyx

    # get basis vectors and add to viewer
    z_vec, y_vec, x_vec = _get_vectors_zyx(df)
    z_vec = np.stack([zyx, z_vec], axis=1)
    y_vec = np.stack([zyx, y_vec], axis=1)
    x_vec = np.stack([zyx, x_vec], axis=1)
    if 'particle z' not in viewer.layers:
        viewer.add_vectors(
            z_vec, name='particle z', length=30, edge_color='orange', edge_width=8
        )
    else:
        viewer.layers['particle z'].data = z_vec
    if 'particle y' not in viewer.layers:
        viewer.add_vectors(
            y_vec, name='particle y', length=10, edge_color='blue', edge_width=8
        )
    else:
        viewer.layers['particle y'].data = y_vec
    if 'particle x' not in viewer.layers:
        viewer.add_vectors(
            x_vec, name='particle x', length=10, edge_color='purple', edge_width=8
        )
    else:
        viewer.layers['particle x'].data = x_vec


def _load_tomograms(
    tilt_series_id,  # TiltSeriesId, the type is created dynamically at runtime
    viewer: napari.viewer.Viewer,
    tomogram_df_grouped: DataFrameGroupBy,
):
    # get data for current tilt-series
    df = tomogram_df_grouped.get_group(tilt_series_id.value)

    # load tomogram
    tomogram_file = df['rlnTomoReconstructedTomogram'].iloc[0]
    tomogram = mrcfile.read(tomogram_file)
    tomogram_binning = float(df['rlnTomoTomogramBinning'].iloc[0])

    # add to napari viewer
    if 'tomogram' in viewer.layers:
        viewer.layers.remove(viewer.layers['tomogram'])
    image_layer = viewer.add_image(
        tomogram,
        name='tomogram',
        rendering='minip',
        blending='translucent',
        depiction='plane',
        plane={
            'position': np.array(tomogram.shape) // 2,
            'thickness': 5,
        },
        scale=[tomogram_binning, tomogram_binning, tomogram_binning]
    )

    # reorder layers and recenter camera
    layer_index = viewer.layers.index(image_layer)
    viewer.layers.insert(0, viewer.layers.pop(layer_index))
    viewer.camera.center = (np.array(tomogram.shape) * tomogram_binning) // 2


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
def view_particles(
    particle_star_file: Path = typer.Option(...),
    tilt_series_star_file: Optional[Path] = typer.Option(
        None,
        help='provide this if you would like to visualise particles in their tomograms'
    )
):
    # read in particle STAR file and get particle dataframe per tomogram
    star = starfile.read(particle_star_file, always_dict=True, parse_as_string=['rlnTomoName'])
    particle_df = star['particles']
    if 'optics' in star:
        particle_df = particle_df.merge(star['optics'], on='rlnOpticsGroup')
    particle_df_grouped = particle_df.groupby('rlnTomoName')

    # if requested, read in tomogram STAR file and get dataframe per tomogram
    if tilt_series_star_file is not None:
        tomogram_df = starfile.read(tilt_series_star_file, parse_as_string=['rlnTomoName'])
        tomogram_df_grouped = tomogram_df.groupby('rlnTomoName')

    # construct a widget for switching between different tomograms
    # we create an Enum type for the different possible IDs, magicgui then
    # creates a dropdown from this type
    TiltSeriesId = Enum(
        'TiltSeriesIds', {
            ts_id: ts_id
            for ts_id in particle_df_grouped.groups.keys()
        }
    )
    container = magicgui.widgets.Container()
    dropdown = magicgui.widgets.create_widget(annotation=TiltSeriesId)
    container.append(dropdown)

    # instantiate a napari viewer and add the widget for switching between tomograms
    viewer = napari.Viewer(ndisplay=3)
    viewer.window.add_dock_widget(container, area='left', name='tilt-series id')

    # when the dropdown changes, call _load_particles to load new particles
    load_particles = partial(
        _load_particles, viewer=viewer, particle_df_grouped=particle_df_grouped
    )
    load_particles = dropdown.changed.connect(load_particles)

    # when the dropdown changes, call _load_tomograms to load the tomogram
    load_tomograms = partial(
        _load_tomograms, viewer=viewer, tomogram_df_grouped=tomogram_df_grouped
    )
    load_tomograms = dropdown.changed.connect(load_tomograms)

    # load particles from first tilt-series on startup
    for ts_id in TiltSeriesId:
        first_tilt_series = ts_id
        break
    load_particles(first_tilt_series)
    load_tomograms(first_tilt_series)

    # start the event loop
    napari.run()
