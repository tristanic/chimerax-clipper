
from Qt import QtCore, QtGui, QtWidgets


class ChooserWindowBase(QtWidgets.QDialog):
    def __init__(self, session, title, chooser_class, allow_multiple=False):
        super().__init__()
        self.setWindowTitle(title)
        from Qt.QtWidgets import QVBoxLayout, QDialogButtonBox

        if allow_multiple:
            selection_mode='extended'
        else:
            selection_mode='single'

        vl = QVBoxLayout()
        mw = self.main_widget=chooser_class(session, parent=self, selection_mode=selection_mode)
        vl.addWidget(mw)
        bb = self.button_box = QDialogButtonBox(self)
        bb.setOrientation(QtCore.Qt.Horizontal)
        cb = self.cancel_button = bb.addButton(QDialogButtonBox.Cancel)
        ok = self.ok_button = bb.addButton(QDialogButtonBox.Ok)
        cb.clicked.connect(self.reject)
        ok.clicked.connect(self.accept)
        vl.addWidget(bb)
        self.setLayout(vl)
    
    @property
    def value(self):
        return self.main_widget.value

class StructureChooserWindow(ChooserWindowBase):
    def __init__(self, session, title, allow_multiple=False):
        from chimerax.atomic.widgets import AtomicStructureListWidget
        super().__init__(session, title, AtomicStructureListWidget, allow_multiple=allow_multiple)
        from Qt.QtWidgets import QDialogButtonBox
        ff = self.from_file_button = self.button_box.addButton('Model from file', QDialogButtonBox.ButtonRole.ActionRole)
        ff.clicked.connect(self._choose_from_file_cb)


    def _choose_from_file_cb(self):
        self.done(2)





