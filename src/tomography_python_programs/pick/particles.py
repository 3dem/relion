from pathlib import Path

import napari
import pandas as pd
import starfile
from napari_threedee.annotators import PointAnnotator
from qtpy.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSizePolicy
from typer import Option

from ._cli import cli
from .._metadata_models.relion.tilt_series_set import TiltSeriesSet
from .._metadata_models.gui.tilt_series_set import TiltSeriesSet as GuiTiltSeriesSet
from .._gui.components.tomogram_browser import TomogramBrowserWidget
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

        self.tomogram_browser_widget.tomogram_changed.connect(
            self.synchronise_annotator
        )
        self.save_button.clicked.connect(self.save_particles)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tomogram_browser_widget)
        self.layout().addWidget(self.save_button)
        self.tomogram_browser_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.tomogram_browser_widget.layout().setContentsMargins(0, 0, 0, 0)
        self.synchronise_annotator()

    def synchronise_annotator(self):
        if 'particles' not in self.viewer.layers:
            self.points_layer = self.viewer.add_points(
                data=[],
                name='particles',
                ndim=3,
                size=20,
                out_of_slice_display=True,
                face_color='cornflowerblue',
            )
        try:
            self.load_particles()
        except FileNotFoundError:
            self.points_layer.data = []
        self.annotator.image_layer = self.viewer.layers['tomogram']
        self.annotator.points_layer = self.viewer.layers['particles']
        self.annotator.enabled = True
        self.viewer.layers.selection = [self.viewer.layers['tomogram']]

    @property
    def current_particle_star_file(self) -> Path:
        tilt_series_id = self.tomogram_browser_widget.selected_tilt_series.name
        return self.output_directory / tilt_series_id / 'particles.star'

    def save_particles(self):
        data = self.viewer.layers['particles'].data
        if len(data) == 0:
            raise ValueError('no particles to save')
        output_directory = self.current_particle_star_file.parent
        output_directory.mkdir(parents=True, exist_ok=True)
        xyz = {
            'rlnCoordinateX': data[:, -1],
            'rlnCoordinateY': data[:, -2],
            'rlnCoordinateZ': data[:, -3],
        }
        df = pd.DataFrame(xyz)
        starfile.write(df, self.current_particle_star_file, overwrite=True)

    def load_particles(self):
        if not self.current_particle_star_file.exists():
            raise FileNotFoundError
        df = starfile.read(self.current_particle_star_file)
        zyx = df[['rlnCoordinateZ', 'rlnCoordinateY', 'rlnCoordinateX']].to_numpy()
        self.points_layer.data = zyx


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def pick_isolated_particles_cli(
    tilt_series_star_file: Path = Option(...),
    output_directory: Path = Option(...),
    cache_size: int = 3
):
    viewer = napari.Viewer(ndisplay=3)
    relion_tilt_series_set = TiltSeriesSet.from_star_file(tilt_series_star_file)
    gui_tilt_series_set = relion_tilt_series_set.as_gui_model()
    dock_widget = PickIsolatedParticlesWidget(
        viewer=viewer,
        tilt_series=gui_tilt_series_set,
        output_directory=output_directory,
        cache_size=cache_size
    )
    viewer.window.add_dock_widget(
        dock_widget, name='RELION particle picker', area='left'
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
