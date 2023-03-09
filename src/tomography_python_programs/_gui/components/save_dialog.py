from qtpy.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QLabel


class SaveDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent=parent)

        QBtn = QDialogButtonBox.Yes | QDialogButtonBox.No

        self.buttonBox = QDialogButtonBox(QBtn)
        # self.buttonBox.accepted.connect(self.accept)
        # self.buttonBox.rejected.connect(self.reject)

        self.layout = QVBoxLayout()
        message = QLabel("Save particles before moving to next tomogram?")
        self.layout.addWidget(message)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)
