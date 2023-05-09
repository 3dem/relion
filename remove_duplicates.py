from enum import Enum
from functools import lru_cache
from pathlib import Path

import napari
import numpy as np
import pandas as pd
import starfile
import magicgui
from scipy.spatial import KDTree

from magicgui.experimental import guiclass
from magicgui.widgets import Button

STAR_FILE = 'run_it007_data.star'
star = starfile.read(STAR_FILE)

df = star['particles']

grouped = df.groupby('rlnTomoName')
TiltSeriesId = Enum('TiltSeriesIds', {ts_id: ts_id for ts_id in grouped.groups.keys()})
for ts_id in TiltSeriesId:
    first_tilt_series = ts_id
    break

viewer = napari.Viewer(ndisplay=3)
widget = magicgui.widgets.create_widget(annotation=TiltSeriesId)


@guiclass
class ParameterClass:
    tilt_series_id: TiltSeriesId = first_tilt_series
    max_distance: int = 1
    output: Path = 'deduplicated.star'


parameters = ParameterClass()


def get_zyx(tilt_series_id: TiltSeriesId) -> np.ndarray:
    df = grouped.get_group(tilt_series_id.value)
    zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()
    if 'rlnOriginZAngst' in df.columns:
        shifts = df[
            ['rlnOriginZAngst', 'rlnOriginYAngst', 'rlnOriginXAngst']].to_numpy()
        zyx -= shifts
    return zyx


@parameters.events.max_distance.connect
@parameters.events.tilt_series_id.connect
def remove_duplicates():
    points = get_zyx(parameters.tilt_series_id)
    points = _collapse_knn(
        points=points,
        max_distance=parameters.max_distance,

    )
    if 'collapsed points' not in viewer.layers:
        viewer.add_points(points, size=40, name='collapsed points')
    else:
        viewer.layers['collapsed points'].data = points
    viewer.camera.center = np.mean(points, axis=0)


def _collapse_knn(
    points: np.ndarray,
    max_distance: float,
    k: int = 15,
) -> np.ndarray:
    tree = KDTree(data=points)
    distance, idx = tree.query(points, k=k, distance_upper_bound=max_distance)

    # remove distances to self
    distance, idx = distance[:, 1:], idx[:, 1:]

    # collapse knn up to distance
    idx_removed = []
    collapsed_points = []
    for i, (_distance, _idx) in enumerate(zip(distance, idx)):
        if i in idx_removed:
            continue
        valid_idx = _idx[_idx < len(points)]
        if len(valid_idx) == 0:
            collapsed_points.append(points[i])
            continue
        point_group = points[valid_idx]
        collapsed_points.append(point_group.mean(axis=0))
        idx_removed.extend(valid_idx)
    return np.stack(collapsed_points, axis=0)


def save_star_file():
    path = parameters.output
    dfs = []
    for ts_id in TiltSeriesId:
        zyx = get_zyx(ts_id)
        tree = KDTree(data=zyx)
        zyx_final = _collapse_knn(zyx, max_distance=parameters.max_distance)
        _, idx = tree.query(zyx_final, k=1)
        df = grouped.get_group(ts_id.value)
        df = df.iloc[idx]
        dfs.append(df)
        print(f'deduplicated {ts_id.value}')
    df = pd.concat(dfs)
    new_star = star.copy()
    new_star['particles'] = df
    starfile.write(star, path, overwrite=True)
    print(f'saving particles to {path}')

    pass


button = Button(text='save STAR file')
button.clicked.connect(save_star_file)
parameters.gui.append(button)

viewer.window.add_dock_widget(parameters.gui, area='left', name='collapse kNN')
remove_duplicates()
napari.run()
