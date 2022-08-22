from functools import lru_cache
from pathlib import Path
from typing import List

import mrcfile
import napari
import numpy as np
from psygnal import Signal
from napari.qt.threading import thread_worker
from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel, QSizePolicy

from .data_models import TiltSeriesSet, TiltSeries
from .tilt_series_list import TiltSeriesListWidget


class TiltSeriesBrowserWidget(QWidget):
    changing_tilt_series: Signal = Signal()
    tilt_series_changed: Signal = Signal()

    def __init__(
            self,
            viewer: napari.Viewer,
            tilt_series: TiltSeriesSet,
            cache_size,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.tilt_series = tilt_series

        self.cache_size = cache_size

        # decorating here to enable setting cache size dynamically
        self._read_image_data = \
            lru_cache(maxsize=cache_size)(self._read_image_data)

        self.tilt_series_list_label = QLabel('tilt-series:')
        self.tilt_series_list_widget = TiltSeriesListWidget(tilt_series)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tilt_series_list_label)
        self.layout().addWidget(self.tilt_series_list_widget)

        self.setMaximumHeight(200)
        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        )

        self.tilt_series_list_widget.currentItemChanged.connect(
            self._on_tilt_series_change
        )
        self.previous_tilt_series = self.viewer.bind_key(
            '[', self.previous_tilt_series
        )
        self.next_tilt_series = self.viewer.bind_key(']', self.next_tilt_series)
        self._on_tilt_series_change()

    @property
    def selected_tilt_series(self) -> TiltSeries:
        return self._selected_tilt_series

    @selected_tilt_series.setter
    def selected_tilt_series(self, value: TiltSeries):
        self._selected_tilt_series = value

    def _on_tilt_series_change(self):
        self.selected_tilt_series = self.tilt_series_list_widget.selected_tilt_series_in_view()
        self.changing_tilt_series.emit()
        self._hide_tilt_series_in_viewer()
        self._load_tilt_series(self.selected_tilt_series.name, add_to_viewer=True)
        self._preload_next_tilt_series()

    def _read_image_data(self, tilt_image_files: List[Path]) -> np.ndarray:
        return np.stack([mrcfile.read(f) for f in tilt_image_files])

    def _add_tilt_series_to_viewer(self, image_data: np.ndarray):
        if 'tilt-series' in self.viewer.layers:
            self.viewer.layers['tilt-series'].data = image_data
        else:
            self.viewer.add_image(image_data, name='tilt-series',
                                  interpolation='bicubic')
        n_images = image_data.shape[0]
        new_dims = np.asarray(self.viewer.dims.point)
        new_dims[0] = n_images // 2
        self.viewer.dims.current_step = tuple(new_dims)
        self._show_tilt_series_in_viewer()

    def _load_tilt_series(self, tilt_series_id: str, add_to_viewer: bool):
        tilt_images = self.tilt_series[tilt_series_id].tilt_image_files
        delayed_read_image_data = thread_worker(self._read_image_data)
        worker = delayed_read_image_data(tuple(tilt_images))
        if add_to_viewer:
            worker.returned.connect(self._add_tilt_series_to_viewer)
            worker.returned.connect(self.tilt_series_changed.emit)
        worker.start()

    def _preload_next_tilt_series(self):
        tilt_series_ids = list(self.tilt_series.keys())
        current_idx = tilt_series_ids.index(self.selected_tilt_series.name)
        starting_idx = current_idx + 1
        n_to_load = self.cache_size - 2
        for i in range(n_to_load):
            idx = starting_idx + i
            if idx >= len(tilt_series_ids):
                break
            self._load_tilt_series(tilt_series_ids[idx], add_to_viewer=False)

    def _hide_tilt_series_in_viewer(self):
        if 'tilt-series' in self.viewer.layers:
            self.viewer.layers['tilt-series'].visible = False

    def _show_tilt_series_in_viewer(self):
        if 'tilt-series' in self.viewer.layers:
            self.viewer.layers['tilt-series'].visible = True

    def previous_tilt_series(self, event=None):
        self.tilt_series_list_widget.previous()

    def next_tilt_series(self, event=None):
        self.tilt_series_list_widget.next()
