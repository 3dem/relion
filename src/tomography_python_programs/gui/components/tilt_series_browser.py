from pathlib import Path
from typing import List, Tuple, Optional
from random import randint

import mrcfile
import napari
import numpy as np
from napari._qt.qthreading import GeneratorWorker
from psygnal import Signal
from napari.qt.threading import thread_worker
from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel, QSizePolicy
from lru import LRU

from .data_models import TiltSeriesSet, TiltSeries
from .tilt_series_list import TiltSeriesListWidget


class TiltSeriesBrowserWidget(QWidget):
    changing_tilt_series: Signal = Signal()
    tilt_series_changed: Signal = Signal()

    def __init__(
            self,
            viewer: napari.Viewer,
            tilt_series: TiltSeriesSet,
            cache_size: int,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.tilt_series = tilt_series

        self.cache_size = cache_size
        self._cache = LRU(cache_size)
        self._active_worker: Optional[GeneratorWorker] = None

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
        if self._active_worker is not None:
            self._active_worker.quit()
        self.selected_tilt_series = self.tilt_series_list_widget.selected_tilt_series_in_view()
        self.changing_tilt_series.emit()
        self._hide_tilt_series_in_viewer()
        self._load_new_tilt_series(self.selected_tilt_series.name, add_to_viewer=True)
        self._preload_next_tilt_series()

    def _set_dims_to_middle_of_tilt_series(self):
        n_images = self.viewer.layers['tilt-series'].data.shape[0]
        new_dims = np.asarray(self.viewer.dims.point)
        new_dims[0] = n_images // 2
        self.viewer.dims.current_step = tuple(new_dims)

    def _update_tilt_series_in_viewer(
            self, read_tilt_series_output: Tuple[str, int, np.ndarray, int, int]
    ):
        tilt_series_id, image_idx, image_data, n_images_completed, n_images_total = read_tilt_series_output
        if 'tilt-series' in self.viewer.layers:
            self.viewer.layers['tilt-series'].data = image_data
        else:
            self.viewer.add_image(
                data=image_data, name='tilt-series', interpolation='bicubic'
            )
        self.viewer.status = f'loaded image: {n_images_completed:02d} / {n_images_total:02d}'
        if image_idx == len(image_data) // 2:
            self._on_tilt_series_loaded()
        if n_images_completed == n_images_total:
            self.viewer.status = 'tilt-series loaded successfully'
        self._show_tilt_series_in_viewer()

    def _load_new_tilt_series(self, tilt_series_id: str, add_to_viewer: bool) -> None:
        if tilt_series_id in self._cache:
            self._load_new_tilt_series_from_cache(tilt_series_id, add_to_viewer)
            return
        tilt_images = self.tilt_series[tilt_series_id].tilt_image_files
        worker = _read_tilt_series(tilt_series_id, list(tilt_images))
        if add_to_viewer:
            self._active_worker = worker
            worker.yielded.connect(self._update_tilt_series_in_viewer)
            worker.finished.connect(self.tilt_series_changed.emit)
        worker.yielded.connect(self._cache_tilt_series)
        worker.start()

    def _load_new_tilt_series_from_cache(self, tilt_series_id: str, add_to_viewer: bool):
        if add_to_viewer is True:
            tilt_series = self._cache[tilt_series_id]
            n = len(tilt_series)
            self._update_tilt_series_in_viewer((tilt_series_id, n, tilt_series, n, n))
            self._on_tilt_series_loaded()

    def _on_tilt_series_loaded(self):
        self._set_dims_to_middle_of_tilt_series()
        self._show_tilt_series_in_viewer()
        self.viewer.layers['tilt-series'].reset_contrast_limits_range()
        self.viewer.layers['tilt-series'].reset_contrast_limits()

    def _cache_tilt_series(
            self, read_tilt_series_output: Tuple[str, int, np.ndarray, int, int]
    ):
        tilt_series_id, image_idx, image_data, n_images_completed, n_images_total = read_tilt_series_output
        if n_images_completed == n_images_total:
            self._cache[tilt_series_id] = image_data

    def _preload_next_tilt_series(self):
        tilt_series_ids = list(self.tilt_series.keys())
        current_idx = tilt_series_ids.index(self.selected_tilt_series.name)
        starting_idx = current_idx + 1
        n_to_load = 2
        for i in range(n_to_load):
            idx = starting_idx + i
            if idx >= len(tilt_series_ids):
                break
            self._load_new_tilt_series(tilt_series_ids[idx], add_to_viewer=False)

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


def _create_empty_tilt_series(tilt_image_files: List[Path], dtype=np.float32) -> np.ndarray:
    with mrcfile.open(tilt_image_files[0], header_only=True) as mrc:
        tilt_series_shape = (len(tilt_image_files), mrc.header.ny, mrc.header.nx)
    return np.zeros(shape=tilt_series_shape, dtype=dtype)


@thread_worker(progress={'total': 40}, start_thread=False)
def _read_tilt_series(tilt_series_id: str, tilt_image_files: List[Path]) -> Tuple:
    result = _create_empty_tilt_series(tilt_image_files, dtype=np.float16)
    indexed_filenames = list(enumerate(tilt_image_files))
    n_images = len(indexed_filenames)

    # load images from the middle outwards and yield
    indexed_filenames.sort(key=lambda t: abs(t[0] - (n_images // 2)))
    for i, (image_index, filename) in enumerate(indexed_filenames, start=1):
        result[image_index] = mrcfile.read(filename)

        if i == 1 or (i % randint(1, 3) == 0) or i == n_images:
            yield tilt_series_id, image_index, result, i, n_images
