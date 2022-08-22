from PyQt5.QtWidgets import QWidget, QHBoxLayout, QPushButton


class SelectionButtons(QWidget):
    def __init__(self):
        super().__init__()
        self.setLayout(QHBoxLayout())

        self.deselect_all_button = QPushButton('deselect all')
        self.select_all_button = QPushButton('select all')

        self.layout().addWidget(self.deselect_all_button)
        self.layout().addWidget(self.select_all_button)
