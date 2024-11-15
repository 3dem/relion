from PyQt5.QtWidgets import QListWidget, QListWidgetItem

from ..._metadata_models.gui.tilt_series_set import GuiTiltSeriesSet
from ..._metadata_models.gui.tilt_series import GuiTiltSeries


class QTiltSeriesItem(QListWidgetItem):
    def __init__(self, tilt_series: GuiTiltSeries, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tilt_series = tilt_series
        self.setText(self.tilt_series.name)


class TiltSeriesListWidget(QListWidget):
    """QListWidget of QTiltSeriesItem."""
    def __init__(self, tilt_series: GuiTiltSeriesSet):
        super().__init__()
        self.setSortingEnabled(False)
        self.setMaximumHeight(200)

        self.tilt_series = tilt_series
        for tilt_series in self.tilt_series.values():
            self.addItem(QTiltSeriesItem(tilt_series))
        self.setCurrentItem(self.item(0))

    def selected_tilt_series_in_view(self) -> GuiTiltSeries:
        return self.currentItem().tilt_series

    def _go_to_index(self, idx: int):
        if self.isEnabled() is False:
            return
        if 0 <= idx < self.count():
            self.setCurrentItem(self.item(idx))

    def next(self):
        current_index = self.currentIndex().row()
        self._go_to_index(current_index + 1)

    def previous(self):
        current_index = self.currentIndex().row()
        self._go_to_index(current_index - 1)

    def __len__(self):
        return len(self.tilt_series)