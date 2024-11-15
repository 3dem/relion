import napari
import numpy as np
from psygnal import Signal
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QWidget, QLabel, QVBoxLayout

from ..._metadata_models.gui.tilt_series import GuiTiltSeries
from .selection_buttons import SelectionButtons
from .checkable_tilt_image_list import TiltImageListWidget


class TiltImageSelectionWidget(QWidget):
    selection_changed: Signal = Signal()

    def __init__(
            self,
            viewer: napari.Viewer,
            tilt_series: GuiTiltSeries,
            *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.viewer = viewer
        self.tilt_image_list_label = QLabel('images:')
        self.tilt_image_list_widget = TiltImageListWidget()
        self.selection_buttons = SelectionButtons()
        self.tilt_series = tilt_series

        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.tilt_image_list_label)
        self.layout().addWidget(self.tilt_image_list_widget)
        self.layout().addWidget(self.selection_buttons)

        self.viewer.dims.events.current_step.connect(self._on_dims_change)
        self.tilt_image_list_widget.currentItemChanged.connect(self._sync_dims_with_tilt_image_list)
        self.tilt_image_list_widget.itemChanged.connect(self.selection_changed.emit)
        self.selection_buttons.select_all_button.clicked.connect(self.tilt_image_list_widget.select_all)
        self.selection_buttons.deselect_all_button.clicked.connect(self.tilt_image_list_widget.deselect_all)

        self.previous_tilt_image = self.viewer.bind_key('Up', self.select_previous_tilt_image)
        self.next_tilt_image = self.viewer.bind_key('Down', self.select_next_tilt_image)
        self.toggle_current_image_check_state = self.viewer.bind_key('Enter', self.toggle_current_image_check_state)
        self.first_tilt_image = self.viewer.bind_key('PageUp', self.select_first_tilt_image)
        self.last_tilt_image = self.viewer.bind_key('PageDown', self.select_last_tilt_image)

    def _on_dims_change(self):
        self.tilt_image_list_widget.blockSignals(True)
        item = self.tilt_image_list_widget.item(int(self.viewer.dims.point[0]))
        self.tilt_image_list_widget.setCurrentItem(item)
        self.tilt_image_list_widget.blockSignals(False)

    def _sync_dims_with_tilt_image_list(self):
        current_index = self.tilt_image_list_widget.currentIndex().row()
        if current_index != -1:  # happens when current item not set
            dims = np.asarray(self.viewer.dims.point)
            dims[0] = current_index
            self.viewer.dims.current_step = tuple(dims)

    @property
    def tilt_series(self) -> GuiTiltSeries:
        return self._tilt_series

    @tilt_series.setter
    def tilt_series(self, value: GuiTiltSeries):
        self._tilt_series = value
        self.tilt_image_list_widget.images = self.tilt_series.tilt_image_files

    @property
    def selected_images_in_view(self) -> np.ndarray:
        """boolean array indicating whether each image is selected or not."""
        n_images = len(self.tilt_series)
        selection = np.zeros(n_images).astype(bool)
        for idx in range(n_images):
            item = self.tilt_image_list_widget.item(idx)
            check_state = True if item.checkState() == Qt.CheckState.Checked else False
            selection[idx] = check_state
        return selection

    @selected_images_in_view.setter
    def selected_images_in_view(self, selection: np.ndarray):
        self.tilt_image_list_widget.blockSignals(True)
        for idx in range(len(self.tilt_series)):
            item = self.tilt_image_list_widget.item(idx)
            check_state = Qt.CheckState.Checked if selection[idx] else Qt.CheckState.Unchecked
            item.setCheckState(check_state)
        self.tilt_image_list_widget.blockSignals(False)

    def select_first_tilt_image(self, event=None):
        self.tilt_image_list_widget._go_to_index(0)

    def select_last_tilt_image(self, event=None):
        self.tilt_image_list_widget._go_to_index(len(self.tilt_image_list_widget) - 1)

    def select_previous_tilt_image(self, event=None):
        self.tilt_image_list_widget.previous()

    def select_next_tilt_image(self, event=None):
        self.tilt_image_list_widget.next()

    def toggle_current_image_check_state(self, event=None):
        current_item = self.tilt_image_list_widget.currentItem()
        if current_item.checkState() == Qt.CheckState.Checked:
            current_item.setCheckState(Qt.CheckState.Unchecked)
        else:
            current_item.setCheckState(Qt.CheckState.Checked)


