from pathlib import Path

import napari
import pandas as pd
import starfile
from napari_threedee.annotators import SphereAnnotator
from napari_threedee.data_models import N3dSpheres
from qtpy.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSizePolicy
from typer import Option

from ._cli import cli
from .._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from .._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from .._qt.components.tomogram_browser import TomogramBrowserWidget
from .._qt.components.save_dialog import SaveDialog
from .._utils.relion import relion_pipeline_job

COMMAND_NAME = 'spheres'


class PickSpheresWidget(QWidget):
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
        self.save_button = QPushButton('save spheres')
        self.save_button.setToolTip(
            'save spheres for the current tomogram to disk.'
        )
        self.annotator = SphereAnnotator(
            viewer=viewer, enabled=True
        )
        self.points_layer = None

        self.tomogram_browser_widget.changing_tomogram.connect(
            self._open_save_dialog
        )
        self.tomogram_browser_widget.tomogram_changed.connect(
            self.synchronise_annotator
        )
        self.save_button.clicked.connect(self.save_particles)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tomogram_browser_widget)
        self.layout().addWidget(self.save_button)
        self.tomogram_browser_widget.setSizePolicy(QSizePolicy.Expanding,
                                                   QSizePolicy.Expanding)
        self.tomogram_browser_widget.layout().setContentsMargins(0, 0, 0, 0)
        self.synchronise_annotator()

    def synchronise_annotator(self):
        self.annotator.enabled = False
        if self.points_layer is not None:
            self.viewer.layers.remove(self.points_layer)
            self.viewer.layers.remove(self.annotator.surface_layer)
        try:
            n3d_spheres = self.load_particles()
            self.points_layer = n3d_spheres.as_layer()
        except FileNotFoundError:
            self.points_layer = N3dSpheres(centers=[], radii=[]).as_layer()
        self.viewer.add_layer(self.points_layer)
        self.annotator = SphereAnnotator(
            viewer=self.viewer,
            image_layer=self.tomogram_browser_widget.image_layer,
            points_layer=self.points_layer,
            enabled=True,
        )
        self.viewer.layers.selection = [self.tomogram_browser_widget.image_layer]

    @property
    def current_particle_star_file(self) -> Path:
        tilt_series_id = self.tomogram_browser_widget.selected_tilt_series.name
        return self.output_directory / 'annotations' / f'{tilt_series_id}_spheres.star'

    def save_particles(self):
        spheres = N3dSpheres.from_layer(self.points_layer)
        if len(spheres.centers) == 0:
            raise ValueError('no spheres to save')
        output_directory = self.current_particle_star_file.parent
        output_directory.mkdir(parents=True, exist_ok=True)
        ts_id = self.tomogram_browser_widget.selected_tilt_series.name
        sphere_data = {
            'rlnTomoName': [ts_id] * len(spheres.centers),
            'rlnCoordinateX': spheres.centers[:, -1],
            'rlnCoordinateY': spheres.centers[:, -2],
            'rlnCoordinateZ': spheres.centers[:, -3],
            'rlnSphereRadius': spheres.radii,
        }
        df = pd.DataFrame(sphere_data)
        starfile.write(df, self.current_particle_star_file, overwrite=True)

    def load_particles(self) -> N3dSpheres:
        if not self.current_particle_star_file.exists():
            raise FileNotFoundError
        df = starfile.read(self.current_particle_star_file, parse_as_string=['rlnTomoName'])
        zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()
        radii = df['rlnSphereRadius'].to_numpy()
        return N3dSpheres(centers=zyx, radii=radii)

    def _open_save_dialog(self):
        if self.annotator.points_layer is None:
            return
        elif len(self.annotator.points_layer.data) > 0:
            dlg = SaveDialog(self)
            dlg.buttonBox.accepted.connect(self.save_particles)
            dlg.buttonBox.accepted.connect(dlg.close)
            dlg.buttonBox.rejected.connect(dlg.close)
            dlg.exec()


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def pick_spheres(
    tilt_series_star_file: Path = Option(...),
    output_directory: Path = Option(...),
    cache_size: int = 3
):
    viewer = napari.Viewer(ndisplay=3)
    relion_tilt_series_set = RlnTiltSeriesSet.from_star_file(tilt_series_star_file)
    gui_tilt_series_set = relion_tilt_series_set.as_gui_model()
    dock_widget = PickSpheresWidget(
        viewer=viewer,
        tilt_series=gui_tilt_series_set,
        output_directory=output_directory,
        cache_size=cache_size
    )
    viewer.window.add_dock_widget(
        dock_widget, name='RELION: 3D sphere picker', area='left'
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
    
    sphere picking
    - 'Alt' + click to pick/move the center of a sphere
    - 'r' to set the radius of the sphere
    - 'n' to enable adding a new sphere
    """
    viewer.text_overlay.visible = True
    napari.run()
