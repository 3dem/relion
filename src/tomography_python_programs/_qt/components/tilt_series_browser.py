from pathlib import Path
from typing import List, Dict, Iterable

import mrcfile
import napari
import numpy as np
from napari._qt.qthreading import GeneratorWorker
from psygnal import Signal
from napari.qt.threading import thread_worker
from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel, QSizePolicy
from superqt.combobox import QSearchableComboBox
from lru import LRU

from ..._metadata_models.gui.lazy_tilt_series_data import LazyTiltSeriesData
from ..._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from ..._metadata_models.gui.tilt_series import GuiTiltSeries
from .tilt_series_list import TiltSeriesListWidget


class TiltSeriesBrowserWidget(QWidget):
    """A widget for browsing sets of tilt-series.

    The widget supports
    - lazy loading/viewing of tilt-series sets
    - eager loading of the next tilt-series
    - caching of tilt-series
    - ability to pause/resume loading when browsing a dataset
    """
    changing_tilt_series: Signal = Signal()
    tilt_series_changed: Signal = Signal()

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
        self._cache = LRU(cache_size)  # {ts_id: LazyTiltSeriesData}
        self._workers: Dict[str, GeneratorWorker] = {}  # {ts_id: worker}

        self.tilt_series_list_label = QLabel('tilt-series:')
        self.tilt_series_combo_box = QSearchableComboBox(parent=self)
        for tilt_series_id in tilt_series:
            self.tilt_series_combo_box.addItem(tilt_series_id)
        self.tilt_series_list_widget = TiltSeriesListWidget(tilt_series)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tilt_series_list_label)
        self.layout().addWidget(self.tilt_series_combo_box)
        self.layout().addWidget(self.tilt_series_list_widget)

        self.setMaximumHeight(200)
        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        )

        self.tilt_series_list_widget.currentItemChanged.connect(
            self.on_tilt_series_change
        )
        self.tilt_series_combo_box.currentIndexChanged.connect(
            self._update_list_from_combobox
        )
        self.tilt_series_list_widget.currentItemChanged.connect(
            self._update_combobox_from_list
        )

        self.previous_tilt_series = self.viewer.bind_key(
            '[', self.previous_tilt_series
        )
        self.next_tilt_series = self.viewer.bind_key(']', self.next_tilt_series)

        self.selected_tilt_series = list(self.tilt_series.values())[0]
        self.on_tilt_series_change()

    @property
    def selected_tilt_series(self) -> GuiTiltSeries:
        return self._selected_tilt_series

    @selected_tilt_series.setter
    def selected_tilt_series(self, value: GuiTiltSeries):
        self._selected_tilt_series = value

    def on_tilt_series_change(self):
        """Called when selected tilt-series in the 'view' has changed.

        model/view - this widget is the model, the tilt-series list widget is the view
        """
        # pause running worker and update selected tilt-series in 'model'
        if self.selected_tilt_series.name in self._workers:
            self._workers[self.selected_tilt_series.name].pause()
        self.selected_tilt_series = self.tilt_series_list_widget.selected_tilt_series_in_view()
        self.changing_tilt_series.emit()

        for worker in self.background_workers:
            worker.pause()
            self._disconnect_worker_safe(worker)

        if self.selected_tilt_series.name not in self._cache:
            self.load_tilt_series_async(self.selected_tilt_series.name,
                                        update_viewer=True)
        else:
            self.load_tilt_series_from_cache(self.selected_tilt_series.name)
        self._remove_workers_for_uncached_tilt_series()

    def load_tilt_series_async(self, tilt_series_id: str, update_viewer: bool = False):
        worker = self._create_thread_worker(tilt_series_id)
        worker = self._connect_events_for_background_loading(worker)
        if update_viewer is True:
            worker = self._connect_events_for_gui_updates(worker)
        self._workers[tilt_series_id] = worker
        worker.start()

    def load_tilt_series_from_cache(self, tilt_series_id: str):
        tilt_series = self._cache[tilt_series_id]
        worker = self._workers[tilt_series_id]
        worker = self._connect_events_for_gui_updates(worker)
        self._update_tilt_series_in_viewer(tilt_series)
        self._setup_viewer_for_tilt_series()
        if worker.is_paused:
            worker.resume()
        elif tilt_series.completed is True:
            self.tilt_series_changed.emit()
            self.preload_next_tilt_series()

    def preload_next_tilt_series(self):
        tilt_series_ids = list(self.tilt_series.keys())
        current_idx = tilt_series_ids.index(self.selected_tilt_series.name)
        starting_idx = current_idx + 1
        n_to_load = 2
        for i in range(n_to_load):
            idx = starting_idx + i
            if idx >= len(tilt_series_ids):
                break
            tilt_series_id = tilt_series_ids[idx]
            if tilt_series_id in self._cache:
                worker = self._workers[tilt_series_id]
                worker = self._connect_events_for_gui_updates(worker)
                if worker.is_paused:
                    worker.resume()
            else:
                self.load_tilt_series_async(tilt_series_id, update_viewer=False)

    def _update_tilt_series_in_viewer(self, tilt_series: LazyTiltSeriesData):
        self._remove_workers_for_uncached_tilt_series()
        if tilt_series.name != self.selected_tilt_series.name:
            return
        if 'tilt-series' in self.viewer.layers:
            self.viewer.layers['tilt-series'].data = tilt_series.data
        else:
            self.viewer.add_image(
                data=tilt_series.data, name='tilt-series', interpolation='bicubic'
            )
        self.viewer.status = _status_from_lazy_tilt_series(tilt_series)
        if tilt_series.n_images_loaded == 1:
            self._setup_viewer_for_tilt_series()

    def _create_thread_worker(self, tilt_series_id: str) -> GeneratorWorker:
        tilt_images = self.tilt_series[tilt_series_id].tilt_image_files
        worker = _read_tilt_series(tilt_series_id, list(tilt_images))
        return worker

    def _connect_events_for_background_loading(self,
                                               worker: GeneratorWorker) -> GeneratorWorker:
        worker.yielded.connect(self._cache_tilt_series)
        return worker

    def _connect_events_for_gui_updates(self,
                                        worker: GeneratorWorker) -> GeneratorWorker:
        worker.yielded.connect(self._update_tilt_series_in_viewer)
        worker.finished.connect(self.tilt_series_changed.emit)
        worker.finished.connect(self.preload_next_tilt_series)
        return worker

    def _disconnect_worker_safe(self, worker: GeneratorWorker):
        try:
            worker.finished.disconnect()
        except TypeError:
            pass
        try:
            worker.yielded.disconnect(self._update_tilt_series_in_viewer)
        except TypeError:
            pass

    def _setup_viewer_for_tilt_series(self):
        self.viewer.dims.ndisplay = 2
        n_images = self.viewer.layers['tilt-series'].data.shape[0]
        new_dims = np.asarray(self.viewer.dims.point)
        new_dims[0] = n_images // 2
        self.viewer.dims.current_step = tuple(new_dims)
        self.viewer.layers['tilt-series'].reset_contrast_limits_range()
        self.viewer.layers['tilt-series'].reset_contrast_limits()

    def _cache_tilt_series(self, tilt_series: LazyTiltSeriesData):
        self._cache[tilt_series.name] = tilt_series

    def _remove_workers_for_uncached_tilt_series(self):
        """Remove workers for tilt-series not in cache."""
        uncached_tilt_series = self._workers.keys() - self._cache.keys()
        for key in uncached_tilt_series:
            if key != self.selected_tilt_series.name:
                worker = self._workers.pop(key)
                worker.quit()

    def previous_tilt_series(self, event=None):
        self.tilt_series_list_widget.previous()

    def next_tilt_series(self, event=None):
        self.tilt_series_list_widget.next()

    @property
    def background_workers(self) -> Iterable[GeneratorWorker]:
        background_worker_keys = self._workers.keys() - {self.selected_tilt_series.name}
        return (self._workers[tilt_series_id] for tilt_series_id in
                background_worker_keys)

    def _update_list_from_combobox(self):
        idx = self.tilt_series_combo_box.currentIndex()
        self.tilt_series_list_widget.setCurrentRow(idx)

    def _update_combobox_from_list(self):
        idx = self.tilt_series_list_widget.currentRow()
        self.tilt_series_combo_box.blockSignals(True)
        self.tilt_series_combo_box.setCurrentIndex(idx)
        self.tilt_series_combo_box.blockSignals(False)


def _create_empty_tilt_series_data(tilt_image_files: List[Path],
                                   dtype=np.float32) -> np.ndarray:
    with mrcfile.open(tilt_image_files[0], header_only=True) as mrc:
        tilt_series_shape = (len(tilt_image_files), mrc.header.ny, mrc.header.nx)
    return np.zeros(shape=tilt_series_shape, dtype=dtype)


@thread_worker(progress={'total': 40}, start_thread=False)
def _read_tilt_series(tilt_series_id: str,
                      tilt_image_files: List[Path]) -> LazyTiltSeriesData:
    tilt_series = _create_empty_tilt_series_data(tilt_image_files, dtype=np.float16)
    lazy_tilt_series_data = LazyTiltSeriesData(
        name=tilt_series_id, data=tilt_series, last_loaded_index=None, n_images_loaded=0
    )
    yield lazy_tilt_series_data

    # load images from the middle outwards and yield
    indexed_filenames = list(enumerate(tilt_image_files))
    n_images = len(indexed_filenames)
    indexed_filenames.sort(key=lambda t: abs(t[0] - (n_images // 2)))

    for i, (image_index, filename) in enumerate(indexed_filenames, start=1):
        tilt_series[image_index] = mrcfile.read(filename)
        lazy_tilt_series_data = LazyTiltSeriesData(
            name=tilt_series_id,
            data=tilt_series,
            last_loaded_index=image_index,
            n_images_loaded=i
        )
        yield lazy_tilt_series_data


def _status_from_lazy_tilt_series(tilt_series: LazyTiltSeriesData) -> str:
    if tilt_series.completed is True:
        status = 'tilt-series loaded successfully'
    else:
        status = f'loaded image: {tilt_series.n_images_loaded:02d} / {tilt_series.n_images:02d}'
    return status
