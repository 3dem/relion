from pathlib import Path
from typing import Optional, Tuple

import napari
import napari.layers
import numpy as np
from qtpy.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSizePolicy
from magicgui.widgets import create_widget
from typer import Option

from ._cli import cli
from .._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from .._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from .._qt.components.tomogram_browser import TomogramBrowserWidget
from .._qt.components.save_dialog import SaveDialog
from .._utils.relion import relion_pipeline_job

COMMAND_NAME = 'surfaces'


class PickSurfacesWidget(QWidget):
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
        self._surface_layer_combo = create_widget(
            value=None,
            annotation=Optional[napari.layers.Surface],
            label="Surface layer",
        )
        self.viewer.layers.events.inserted.connect(self._surface_layer_combo.reset_choices)
        self.viewer.layers.events.removed.connect(self._surface_layer_combo.reset_choices)
        self.viewer.layers.events.reordered.connect(self._surface_layer_combo.reset_choices)

        self.tomogram_browser_widget = TomogramBrowserWidget(
            viewer=viewer,
            tilt_series=tilt_series,
            cache_size=cache_size,
        )
        self.save_button = QPushButton('save surfaces')
        self.save_button.setToolTip(
            'save surfaces for the current tomogram to disk.'
        )

        self.tomogram_browser_widget.changing_tomogram.connect(
            self._open_save_dialog
        )
        self.tomogram_browser_widget.tomogram_changed.connect(
            self.synchronise_annotator
        )
        self.save_button.clicked.connect(self.save_particles)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tomogram_browser_widget)
        self.layout().addWidget(self._surface_layer_combo.native)
        self.layout().addWidget(self.save_button)
        self.tomogram_browser_widget.setSizePolicy(QSizePolicy.Expanding,
                                                   QSizePolicy.Expanding)
        self.tomogram_browser_widget.layout().setContentsMargins(0, 0, 0, 0)
        self.synchronise_annotator()

    def synchronise_annotator(self):
        surface_layer = self._surface_layer_combo.value
        if surface_layer is not None:
            self.viewer.layers.remove(surface_layer)
        try:
            layer_data = self.load_particles()
            self._surface_layer_combo.value = self.viewer.add_surface(layer_data)
        except FileNotFoundError:
            pass
        self.viewer.layers.selection = [self.tomogram_browser_widget.image_layer]

    @property
    def current_particle_star_file(self) -> Path:
        tilt_series_id = self.tomogram_browser_widget.selected_tilt_series.name
        return self.output_directory / 'annotations' / f'{tilt_series_id}_surfaces.star'

    @property
    def current_surfaces_file(self) -> Path:
        return self.current_particle_star_file.with_suffix(".npz")

    def surface_empty(self) -> bool:
        vertices, indices = self.vertices_indices()
        return vertices is None or indices is None or len(vertices) == 0 or len(indices) == 0

    def vertices_indices(self) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        if self._surface_layer_combo.value is None:
            return None, None
        surface_data = self._surface_layer_combo.value.data
        if len(surface_data) == 2:
            vertices, indices = surface_data
        else:
            vertices, indices, _ = surface_data
        return vertices, indices

    def save_particles(self):
        if self.surface_empty():
            raise ValueError('no surfaces to save')
        vertices, indices = self.vertices_indices()
        output_directory = self.current_particle_star_file.parent
        output_directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(self.current_surfaces_file, vertices=vertices, indices=indices)

    def load_particles(self) -> Tuple[np.ndarray, np.ndarray]:
        file = np.load(self.current_surfaces_file)
        return file["vertices"], file["indices"]

    def _open_save_dialog(self):
        if self.surface_empty():
            return
        dlg = SaveDialog(self)
        dlg.buttonBox.accepted.connect(self.save_particles)
        dlg.buttonBox.accepted.connect(dlg.close)
        dlg.buttonBox.rejected.connect(dlg.close)
        dlg.exec()


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
@relion_pipeline_job
def pick_surfaces(
    tilt_series_star_file: Path = Option(...),
    output_directory: Path = Option(...),
    cache_size: int = 3
):
    viewer = napari.Viewer(ndisplay=3)
    relion_tilt_series_set = RlnTiltSeriesSet.from_star_file(tilt_series_star_file)
    gui_tilt_series_set = relion_tilt_series_set.as_gui_model()
    dock_widget = PickSurfacesWidget(
        viewer=viewer,
        tilt_series=gui_tilt_series_set,
        output_directory=output_directory,
        cache_size=cache_size
    )
    viewer.window.add_dock_widget(
        dock_widget, name='RELION: 3D surface picker', area='left'
    )
    viewer.axes.visible = True
    viewer.axes.labels = False
    viewer.camera.angles = (-15, 30, 150)
    viewer.camera.zoom *= 0.7
    viewer.text_overlay.text = """
    Load a surface into napari by e.g. dragging an IMOD model into the viewer.
    """
    viewer.text_overlay.visible = True
    napari.run()
