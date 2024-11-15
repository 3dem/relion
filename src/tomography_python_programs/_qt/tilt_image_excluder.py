from typing import Dict

import napari
import numpy as np
from qtpy.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSizePolicy

from .._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from .components.tilt_image_selection_widget import TiltImageSelectionWidget
from .components.tilt_series_browser import TiltSeriesBrowserWidget


class TiltImageExcluderWidget(QWidget):
    def __init__(
            self,
            viewer: napari.Viewer,
            tilt_series: GuiTiltSeriesSet,
            cache_size: int,
            save_button_label: str,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.tilt_series = tilt_series
        self.selected_tilt_images: Dict[str, np.ndarray] = {}  # {tilt_series_id: idx}
        self._set_default_excluded_images()

        first_tilt_series = self.tilt_series[list(self.tilt_series.keys())[0]]
        self.tilt_image_selection_widget = TiltImageSelectionWidget(
            viewer=self.viewer,
            tilt_series=first_tilt_series
        )
        self.tilt_series_browser_widget = TiltSeriesBrowserWidget(
            viewer=self.viewer,
            tilt_series=self.tilt_series,
            cache_size=cache_size,
        )
        self.save_button = QPushButton(save_button_label)

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tilt_image_selection_widget)
        self.layout().addWidget(self.tilt_series_browser_widget)
        self.layout().addWidget(self.save_button)

        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        )

        self.tilt_series_browser_widget.changing_tilt_series.connect(self._on_tilt_series_change)
        self.tilt_series_browser_widget.tilt_series_changed.connect(self._post_tilt_series_change)
        self.tilt_image_selection_widget.selection_changed.connect(self._on_tilt_image_selection_change)
        self.save_button.clicked.connect(self._pre_save)
        self.save_button.clicked.connect(self.save_output)

        self._on_tilt_series_change()

    def _on_tilt_series_change(self):
        self.tilt_image_selection_widget.setDisabled(True)
        self.tilt_image_selection_widget.tilt_series = self.tilt_series_browser_widget.selected_tilt_series
        self._update_tilt_image_checkboxes()

    def _post_tilt_series_change(self):
        self.tilt_image_selection_widget.setEnabled(True)

    def _update_tilt_image_checkboxes(self):
        selection = self.selected_tilt_images[self.tilt_image_selection_widget.tilt_series.name]
        self.tilt_image_selection_widget.selected_images_in_view = selection

    def _set_default_excluded_images(self):
        for tilt_series_id, tilt_series in self.tilt_series.items():
            n_images = len(tilt_series.tilt_image_files)
            self.selected_tilt_images[tilt_series_id] = np.ones(n_images).astype(bool)

    def _on_tilt_image_selection_change(self):
        current_tilt_series_id = self.tilt_series_browser_widget.selected_tilt_series.name
        current_selection = self.tilt_image_selection_widget.selected_images_in_view
        self.selected_tilt_images[current_tilt_series_id] = current_selection

    def _pre_save(self):
        self.save_button.setDisabled(True)
        self.tilt_series_browser_widget.setDisabled(True)
        self.viewer.status = 'saving tilt-image excluder output...'

    def save_output(self):
        """implement in subclasses to save output data."""
        raise NotImplementedError

    def _post_successful_save(self):
        """Call this after saving output in save_output implementations."""
        self.viewer.status = 'tilt-image excluder output saved successfully!'
        self.save_button.setEnabled(True)
        self.tilt_series_browser_widget.setEnabled(True)


