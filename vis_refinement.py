from enum import Enum

import napari
import numpy as np
import starfile
import magicgui

STAR_FILE = 'run_it007_data.star'
star = starfile.read(STAR_FILE)

df = star['particles']

grouped = df.groupby('rlnTomoName')
TiltSeriesId = Enum('TiltSeriesIds', {ts_id: ts_id for ts_id in grouped.groups.keys()})
for ts_id in TiltSeriesId:
    first_tilt_series = ts_id
    break

widget = magicgui.widgets.create_widget(annotation=TiltSeriesId)

viewer = napari.Viewer()
viewer.window.add_dock_widget(widget, area='left', name='tilt-series id')


@widget.changed.connect
def load_tilt_series(tilt_series_id: TiltSeriesId):
    tilt_series_id = tilt_series_id.value
    zyx = get_zyx(tilt_series_id)
    if 'particle positions' not in viewer.layers:
        viewer.add_points(zyx, name='particle positions', size=40)
    else:
        viewer.layers['particle positions'].data = zyx


def get_zyx(tilt_series_id: str) -> np.ndarray:
    df = grouped.get_group(tilt_series_id)
    zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()
    if 'rlnOriginZAngst' in df.columns:
        shifts = df[['rlnOriginZAngst', 'rlnOriginYAngst', 'rlnOriginXAngst']].to_numpy()
        zyx -= shifts
    return zyx

load_tilt_series(first_tilt_series)
napari.run()