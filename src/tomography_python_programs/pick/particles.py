from pathlib import Path

import napari
import numpy as np
import pandas as pd
import starfile
from napari_threedee.annotators import PointAnnotator
from napari_threedee.data_models import N3dPoints
from qtpy.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSizePolicy
from typer import Option
from napari.utils.notifications import show_info

from ._cli import cli
from .._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from .._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from .._qt.components.save_dialog import SaveDialog
from .._qt.components.tomogram_browser import TomogramBrowserWidget
from .._utils.relion import relion_pipeline_job

COMMAND_NAME = 'particles'


class PickIsolatedParticlesWidget(QWidget):
    def __init__(
        self,
        viewer: napari.Viewer,
        tilt_series: GuiTiltSeriesSet,
        output_directory: Path,
        cache_size: int = 3,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.output_directory = output_directory

        self.tomogram_browser_widget = TomogramBrowserWidget(
            viewer=viewer,
            tilt_series=tilt_series,
            cache_size=cache_size,
        )
        self.save_button = QPushButton('save particles')
        self.save_button.setToolTip(
            'save particles for the current tomogram to disk.'
        )
        self.annotator = PointAnnotator(
            viewer=viewer, enabled=True
        )

        self.tomogram_browser_widget.changing_tomogram.connect(
            self._open_save_dialog
        )
        self.tomogram_browser_widget.tomogram_changed.connect(
            self.synchronise_annotator
        )
        self.save_button.clicked.connect(self.save_particles)
        self.points_layer = None

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tomogram_browser_widget)
        self.layout().addWidget(self.save_button)
        self.tomogram_browser_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.tomogram_browser_widget.layout().setContentsMargins(0, 0, 0, 0)
        self.synchronise_annotator()

    def synchronise_annotator(self):
        self.annotator.enabled = False
        if self.points_layer is not None:
            self.viewer.layers.remove(self.points_layer)
        try:
            n3d_points = self.load_particles()
            self.points_layer = n3d_points.as_layer()
        except FileNotFoundError:
            self.points_layer = N3dPoints(data=[]).as_layer()
        self.points_layer.current_face_color = 'cornflowerblue'
        self.points_layer.face_color = 'cornflowerblue'
        self.viewer.add_layer(self.points_layer)
        self.annotator = PointAnnotator(
            viewer=self.viewer,
            image_layer=self.tomogram_browser_widget.image_layer,
            points_layer=self.points_layer,
            enabled=True
        )
        self.viewer.layers.selection = [self.tomogram_browser_widget.image_layer]

    @property
    def current_particle_star_file(self) -> Path:
        tilt_series_id = self.tomogram_browser_widget.selected_tilt_series.name
        return self.output_directory / 'annotations' / f'{tilt_series_id}_particles.star'

    def save_particles(self):
        zyx = self.points_layer.data
        if len(zyx) == 0:
            raise ValueError('no particles to save')
        output_directory = self.current_particle_star_file.parent
        output_directory.mkdir(parents=True, exist_ok=True)
        ts_id = self.tomogram_browser_widget.selected_tilt_series.name
        xyz = {
            'rlnTomoName': [ts_id] * len(zyx),
            'rlnCoordinateX': zyx[:, -1],
            'rlnCoordinateY': zyx[:, -2],
            'rlnCoordinateZ': zyx[:, -3],
        }
        df = pd.DataFrame(xyz)
        starfile.write(df, self.current_particle_star_file, overwrite=True)
        show_info(f'{len(df)} particles saved to {self.current_particle_star_file}')

    def load_particles(self) -> N3dPoints:
        if not self.current_particle_star_file.exists():
            raise FileNotFoundError
        df = starfile.read(self.current_particle_star_file, parse_as_string=['rlnTomoName'])
        zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()

        return N3dPoints(data=zyx)

    def _open_save_dialog(self):
        if self.points_layer is None:
            return
        elif len(self.annotator.points_layer.data) > 0:
            dlg = SaveDialog(self)
            dlg.buttonBox.accepted.connect(self.save_particles)
            dlg.buttonBox.accepted.connect(dlg.close)
            dlg.buttonBox.rejected.connect(dlg.close)
            dlg.exec()


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def pick_particles(
    tilt_series_star_file: Path = Option(...),
    output_directory: Path = Option(...),
    cache_size: int = 3
):
    viewer = napari.Viewer(ndisplay=3)
    relion_tilt_series_set = RlnTiltSeriesSet.from_star_file(tilt_series_star_file)
    gui_tilt_series_set = relion_tilt_series_set.as_gui_model()
    dock_widget = PickIsolatedParticlesWidget(
        viewer=viewer,
        tilt_series=gui_tilt_series_set,
        output_directory=output_directory,
        cache_size=cache_size
    )
    viewer.window.add_dock_widget(
        dock_widget, name='RELION: 3D particle picker', area='left'
    )
    viewer.axes.visible = True
    viewer.axes.labels = False
    viewer.camera.angles = (-15, 30, 150)
    viewer.camera.zoom *= 0.7
    viewer.text_overlay.text = """
    keyboard shortcuts
    - '[' and ']' for previous/next tomogram
    - 'x'/'y'/'z'/'o' to reorient the plane (around the cursor)

    mouse controls
    - 'Shift' + click + drag to move the plane along its normal
    - 'Alt' + click to pick a particle
    """
    viewer.text_overlay.visible = True
    napari.run()
