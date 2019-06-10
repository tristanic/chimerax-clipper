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

import numpy

class StructureStepper:
    '''
    Class to aid in stepping through an atomic model some number of secondary
    structure elements (plus flanking loops) at a time. Starts at the beginning
    of the first polymer in the structure.
    '''

    FORWARD = 1
    BACKWARD = -1
    def __init__(self, session, model, step = 2, min_gap = 2,
                num_non_protein_residues_per_step = 15,
                non_protein_overlap = 5):
        '''
        Args:
            * session:
                - the ChimeraX session
            * model:
                - a :py:class:`AtomicStructure` instance
            * step:
                - the number of secondary structure elements to be covered in each
                  step
            * min_gap:
                - If a secondary structure element contains fewer than this number
                  of residues, it will be rolled into the next element
            * num_non_protein_residues_per_step:
                - Non-protein chains do not have secondary structure definitions,
                  so will just step through this number of residues at a time.
            * non_protein_overlap:
                - When stepping through a non-protein chain, keep this many
                  residues from the previous selection as the start of the new
                  one.
        '''
        self.session = session
        self.model = model
        self.step = step
        self.min_gap = min_gap
        self._np_sel_size = num_non_protein_residues_per_step
        self._np_overlap_size = non_protein_overlap
        cp = self.current_polymer = self._get_polymers()[0]
        self._current_polymer_index = 0
        self._current_sel = None
        self._new_chain = False

    def _get_polymers(self):
        return self.model.polymers(self.model.PMS_NEVER_CONNECTS)

    def step_forward(self):
        return self._step(direction=self.FORWARD)

    def step_backward(self):
        return self._step(direction=self.BACKWARD)

    def _step(self, direction):
        from chimerax.atomic import Residue
        if self.current_polymer[1] == Residue.PT_AMINO:
            return self._step_protein(direction)
        else:
            return self._step_generic(direction)

    def _step_generic(self, direction):
        cs = self._current_sel
        cp = self.current_polymer[0]
        nres = len(cp)

        sel_size = self._np_sel_size
        overlap = self._np_overlap_size

        if direction == self.BACKWARD:
            cp = cp[::-1]
            if cs is not None:
                cs = cs[::-1]

        if cs is None or self._new_chain:
            self._new_chain = False
            res = cp[:sel_size]
            if direction == self.BACKWARD:
                res = res[::-1]
            sel = self._current_sel = res.atoms
        else:
            indices = cp.indices(cs)
            first_index = indices[-overlap]
            last_index = first_index+sel_size+1
            if last_index > nres:
                self._go_to_next_chain(direction)
            res = cp[first_index:first_index+sel_size+1]
            if direction == self.BACKWARD:
                res = res[::-1]
            sel = self._current_sel = res.atoms
        return sel




    def _step_protein(self, direction):
        from chimerax.atomic import Residue, Residues
        cs = self._current_sel
        cp = self.current_polymer[0]
        nres = len(cp)

        if direction == self.BACKWARD:
            cp = cp[::-1]
            if cs is not None:
                cs = cs[::-1]

        c_ssids = cp.secondary_structure_ids
        c_sstypes = cp.ss_types


        if cs is None or self._new_chain:
            self._new_chain = False
            first_element = c_ssids[0]
        else:
            cr = cs.unique_residues
            ids = cr.secondary_structure_ids
            defined_ss = ids[cr.ss_types != Residue.SS_COIL]
            # Last secondary structure element of previous selection
            first_element = defined_ss[-1]

        indices = []

        this_element = first_element
        element_count = 0


        while True:
            new_indices = (numpy.argwhere(c_ssids == this_element).ravel())
            if len(new_indices) >= self.min_gap and c_sstypes[new_indices[0]] != Residue.SS_COIL:
                element_count += 1
            indices.extend(new_indices)
            if indices[-1] == nres-1:
                # End of this fragment. Move on to the next on the next iteration
                self._go_to_next_chain(direction)
                break
            next_index = indices[-1]+1
            this_element = c_ssids[next_index]
            this_element_type = c_sstypes[next_index]
            if element_count == self.step:
                # Grab the flanking loop and finish
                if this_element_type == Residue.SS_COIL:
                    indices.extend(numpy.argwhere(c_ssids == this_element).ravel())
                if indices[-1] == nres-1:
                    self._go_to_next_chain(direction)
                break

        indices = numpy.array(indices)
        selres = cp[indices]
        if direction == self.BACKWARD:
            # put the selection back in forward order
            selres = selres[::-1]
        selatoms = self._current_sel = selres.atoms
        return selatoms




    def _go_to_next_chain(self, direction):
            polymers = self._get_polymers()
            np = len(polymers)
            ci = self._current_polymer_index
            if direction == self.FORWARD and ci == np-1:
                self.current_polymer = polymers[0]
                self._current_polymer_index = 0
            elif direction == self.BACKWARD and ci == 0:
                self.current_polymer = polymers[-1]
                self._current_polymer_index = np-1
            else:
                self._current_polymer_index += direction
                self.current_polymer = polymers[self._current_polymer_index]
            self._new_chain = True



def _res_ss_type(residue):
  from chimerax.atomic import Residue
  if residue.is_strand:
    return Residue.SS_STRAND
  if residue.is_helix:
    return Residue.SS_HELIX
  return -1
