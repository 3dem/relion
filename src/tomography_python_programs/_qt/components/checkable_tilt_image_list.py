from pathlib import Path
from typing import Optional, List, Sequence

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QListWidgetItem, QListWidget, QSizePolicy


class QTiltImageItem(QListWidgetItem):
    def __init__(self, tilt_image_file: Path):
        super().__init__()
        self.file = tilt_image_file
        self.setText(self.file.name)
        self.setFlags(self.flags() | Qt.ItemIsUserCheckable)
        self.setCheckState(Qt.CheckState.Checked)


class TiltImageListWidget(QListWidget):
    """QListWidget of QTiltImageItem"""

    def __init__(self):
        super().__init__()
        self.images = []
        self.setSizePolicy(
            QSizePolicy(
                QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding
            )
        )

    @property
    def selected_tilt_image(self) -> Optional[Path]:
        return self.currentItem().file if self.currentItem() else None

    @property
    def images(self) -> List[Path]:
        return self._tilt_images

    @images.setter
    def images(self, tilt_images: Sequence[Path]):
        self._tilt_images = list(tilt_images)
        self._repopulate_list()

    def _repopulate_list(self):
        self.clear()
        for tilt_image in self.images:
            self.addItem(QTiltImageItem(tilt_image))

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

    def select_all(self):
        for i in range(self.count()):
            self.item(i).setCheckState(Qt.CheckState.Checked)

    def deselect_all(self):
        for i in range(self.count()):
            self.item(i).setCheckState(Qt.CheckState.Unchecked)
