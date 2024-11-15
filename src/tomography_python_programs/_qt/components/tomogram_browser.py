import napari
import mrcfile
import numpy as np
from psygnal import Signal
from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel, QSizePolicy
from lru import LRU

from ..._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from ..._metadata_models.gui.tilt_series import GuiTiltSeries
from .tilt_series_list import TiltSeriesListWidget

IMAGE_LAYER_NAME = 'tomogram'


class TomogramBrowserWidget(QWidget):
    changing_tomogram: Signal = Signal()
    tomogram_changed: Signal = Signal()
    image_layer: napari.layers.Image

    def __init__(
        self,
        viewer: napari.Viewer,
        tilt_series: GuiTiltSeriesSet,
        cache_size: int,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.tilt_series = tilt_series

        self.cache_size = cache_size
        self._cache = LRU(cache_size)

        self.tilt_series_list_label = QLabel('tomograms:')
        self.tilt_series_list_widget = TiltSeriesListWidget(tilt_series)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tilt_series_list_label)
        self.layout().addWidget(self.tilt_series_list_widget)

        self.setMaximumHeight(200)
        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        )

        self.tilt_series_list_widget.currentItemChanged.connect(
            self._on_tomogram_selection_change
        )
        self.previous_tilt_series = self.viewer.bind_key('[', self.previous_tilt_series)
        self.next_tilt_series = self.viewer.bind_key(']', self.next_tilt_series)
        self._on_tomogram_selection_change()

    @property
    def image_layer(self) -> napari.layers.Image:
        return self.viewer.layers[IMAGE_LAYER_NAME]

    @property
    def selected_tilt_series(self) -> GuiTiltSeries:
        return self._selected_tilt_series

    @selected_tilt_series.setter
    def selected_tilt_series(self, value: GuiTiltSeries):
        self._selected_tilt_series = value

    def _on_tomogram_selection_change(self):
        self.changing_tomogram.emit()
        self.selected_tilt_series = self.tilt_series_list_widget.selected_tilt_series_in_view()
        self._load_tomogram(self.selected_tilt_series.name, add_to_viewer=True)
        self._preload_next_tomograms()

    def _update_tomogram_in_viewer(self, tomogram: np.ndarray):
        if 'tomogram' in self.viewer.layers:
            self.viewer.layers[IMAGE_LAYER_NAME].data = tomogram
        else:
            self.viewer.add_image(
                data=tomogram,
                name=IMAGE_LAYER_NAME,
                depiction='plane',
                plane={'thickness': 1},
                rendering='minip',
                blending='translucent',
            )
        self._on_tomogram_loaded()

    def _load_tomogram(self, tilt_series_id: str, add_to_viewer: bool) -> None:
        if tilt_series_id in self._cache:
            self._load_tomogram_from_cache(tilt_series_id, add_to_viewer)
            return
        denoised_tomogram_file = self.tilt_series[tilt_series_id].denoised_tomogram_file
        tomogram_file = self.tilt_series[tilt_series_id].tomogram_file
        tomogram_file = denoised_tomogram_file if denoised_tomogram_file is not None else tomogram_file
        tomogram = mrcfile.read(tomogram_file)
        self._cache[tilt_series_id] = tomogram
        if add_to_viewer:
            self._update_tomogram_in_viewer(tomogram)

    def _load_tomogram_from_cache(self, tilt_series_id: str, add_to_viewer: bool):
        if add_to_viewer is True:
            self._update_tomogram_in_viewer(self._cache[tilt_series_id])

    def _on_tomogram_loaded(self):
        layer = self.viewer.layers[IMAGE_LAYER_NAME]
        layer.depiction = 'plane'
        layer.plane.position = np.array(layer.data.shape) / 2
        layer.plane.normal = (1, 0, 0)
        self.viewer.layers.selection = [layer]
        self.tomogram_changed.emit()

    def _preload_next_tomograms(self):
        tilt_series_ids = list(self.tilt_series.keys())
        current_idx = tilt_series_ids.index(self.selected_tilt_series.name)
        starting_idx = current_idx + 1
        n_to_load = 1
        for i in range(n_to_load):
            idx = starting_idx + i
            if idx >= len(tilt_series_ids):
                break
            self._load_tomogram(tilt_series_ids[idx], add_to_viewer=False)

    def previous_tilt_series(self, event=None):
        self.tilt_series_list_widget.previous()

    def next_tilt_series(self, event=None):
        self.tilt_series_list_widget.next()
