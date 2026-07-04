# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Note that this software makes use of modified versions of the Clipper, LibCCP4
# and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
# of these is redistributed under its own license terms.

'''
Optional GUI for choosing exactly what to load from a structure-factor file.

Strictly opt-in: nothing here runs unless the user asks for it (``clipper open
... browse true`` or the Tools menu entry). It is a thin override layer over the
automatic classifier — every default shown is what the automatic path would have
chosen, so clicking OK with no changes reproduces the silent behaviour.
'''

from Qt import QtCore, QtWidgets


def run_mtz_browser(session, container):
    '''
    Show the browser for an already-loaded ReflectionDataContainer and return the
    user's selections as {'fsigf_name': str|None, 'map_columns': [str]|None}, or
    None if the user cancelled. Any change to the free-R flag selection is applied
    to `container` directly before returning.
    '''
    dlg = MTZBrowser(session, container)
    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
        return None
    return dlg.apply_and_get_selection()


def _known_experimental_types():
    from chimerax.clipper import (HKL_data_F_sigF, HKL_data_F_sigF_ano,
        HKL_data_I_sigI, HKL_data_I_sigI_ano)
    # Value is the priority used to pick a sensible default (lower = preferred).
    return {
        HKL_data_F_sigF: 0, HKL_data_I_sigI: 1,
        HKL_data_F_sigF_ano: 2, HKL_data_I_sigI_ano: 3,
    }


class MTZBrowser(QtWidgets.QDialog):
    def __init__(self, session, container):
        super().__init__()
        self.session = session
        self.container = container
        self._separate_file = None
        self._free_candidates = []
        self.setWindowTitle('Choose data and maps to load')

        from Qt.QtWidgets import (QVBoxLayout, QLabel, QComboBox, QListWidget,
            QListWidgetItem, QDialogButtonBox, QHBoxLayout, QPushButton)
        layout = QVBoxLayout()

        import os
        layout.addWidget(QLabel('<b>{}</b>'.format(
            os.path.basename(container.filename or 'structure factors'))))

        # --- experimental data (drives the live, model-coupled maps) ----------
        layout.addWidget(QLabel('Experimental data (for live maps):'))
        self.exp_combo = QComboBox()
        self._populate_experimental()
        layout.addWidget(self.exp_combo)

        # --- map coefficients to open as static maps --------------------------
        layout.addWidget(QLabel('Maps to open:'))
        self.map_list = QListWidget()
        self._populate_maps()
        layout.addWidget(self.map_list)

        # --- free-R flags -----------------------------------------------------
        layout.addWidget(QLabel('Free-R flags:'))
        free_row = QHBoxLayout()
        self.free_combo = QComboBox()
        self._populate_free_flags()
        free_row.addWidget(self.free_combo, stretch=1)
        from_file = QPushButton('From separate file…')
        from_file.clicked.connect(self._choose_free_flags_file)
        free_row.addWidget(from_file)
        layout.addLayout(free_row)
        self.free_file_label = QLabel('')
        self.free_file_label.setWordWrap(True)
        layout.addWidget(self.free_file_label)

        bb = QDialogButtonBox(self)
        bb.setOrientation(QtCore.Qt.Orientation.Horizontal)
        bb.addButton(QDialogButtonBox.StandardButton.Cancel).clicked.connect(self.reject)
        bb.addButton(QDialogButtonBox.StandardButton.Ok).clicked.connect(self.accept)
        layout.addWidget(bb)
        self.setLayout(layout)

    # --- population ----------------------------------------------------------

    def _populate_experimental(self):
        from ..io.reflection_classify import experimental_sort_key
        priority = _known_experimental_types()
        datasets = self.container.experimental_data.datasets
        entries = [(name, d) for name, d in datasets.items()
            if d.dtype in priority]
        entries.sort(key=lambda nd: experimental_sort_key(
            nd[0], priority[nd[1].dtype]))
        for name, d in entries:
            self.exp_combo.addItem(
                '{}  [{}]'.format(name, d.dtype.__name__), name)
        if not entries:
            self.exp_combo.addItem('(none found)', None)

    def _populate_maps(self):
        from Qt.QtWidgets import QListWidgetItem
        from ..io.reflection_classify import map_coeff_info
        for dataset in self.container.calculated_data:
            info = map_coeff_info(dataset.name)
            item = QListWidgetItem('{}   —   {}'.format(
                dataset.name, info.human_label))
            item.setData(QtCore.Qt.ItemDataRole.UserRole, dataset.name)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            # Default: open everything except probable Fcalc maps (matches the
            # automatic path's implicit choice).
            item.setCheckState(QtCore.Qt.CheckState.Unchecked if info.is_fcalc
                else QtCore.Qt.CheckState.Checked)
            self.map_list.addItem(item)

    def _populate_free_flags(self):
        from ..clipper_mtz import list_free_flag_candidates
        from ..io.reflection_classify import rank_free_flags
        try:
            _, candidates = list_free_flag_candidates(
                self.session, self.container.filename)
        except Exception as e:
            self.session.logger.info(
                'Could not enumerate free-flag columns: {}'.format(e))
            candidates = []
        self._free_candidates = candidates
        current = self.container.free_flags
        current_name = current.name if current is not None else None
        # Present recognised free-set names first.
        ranked = rank_free_flags([n for n, _ in candidates])
        order = ranked + [n for n, _ in candidates if n not in ranked]
        default_index = 0
        for i, name in enumerate(order):
            self.free_combo.addItem(name, name)
            if name == current_name:
                default_index = i
        if not candidates:
            self.free_combo.addItem('(auto / none)', None)
        self.free_combo.setCurrentIndex(default_index)

    # --- callbacks -----------------------------------------------------------

    def _choose_free_flags_file(self):
        from Qt.QtWidgets import QFileDialog
        import os
        start = os.path.dirname(self.container.filename or '')
        filename, _ = QFileDialog.getOpenFileName(
            self, 'Choose a file containing free-R flags', start,
            'Reflection files (*.mtz *.cif)')
        if filename:
            self._separate_file = filename
            self.free_file_label.setText(
                'Free flags from: {}'.format(os.path.basename(filename)))

    # --- result --------------------------------------------------------------

    def apply_and_get_selection(self):
        # Free-R flags: a separate file takes precedence over the in-file choice.
        if self._separate_file is not None:
            from ..clipper_mtz import adopt_free_flags_from_file
            adopt_free_flags_from_file(
                self.session, self.container, self._separate_file)
        else:
            chosen = self.free_combo.currentData()
            current = self.container.free_flags
            current_name = current.name if current is not None else None
            if chosen is not None and chosen != current_name:
                self._adopt_in_file_free_column(chosen)

        fsigf_name = self.exp_combo.currentData()
        map_columns = []
        for i in range(self.map_list.count()):
            item = self.map_list.item(i)
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                map_columns.append(item.data(QtCore.Qt.ItemDataRole.UserRole))
        return {'fsigf_name': fsigf_name, 'map_columns': map_columns}

    def _adopt_in_file_free_column(self, name):
        from ..clipper_mtz import _match_flag_array
        for cand_name, array in self._free_candidates:
            if cand_name == name:
                matched = _match_flag_array(array, self.container.hklinfo)
                self.container.set_free_flags(name, matched)
                self.session.logger.info(
                    'Using free-R flags "{}".'.format(name))
                return
