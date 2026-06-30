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
Choose the correct X-ray scattering-factor species for each atom in a model.

Clipper tabulates form factors for neutral elements *and* their common ions
(``Na1+``, ``Mg2+``, ``Fe2+``/``Fe3+``, ``Cl1-``, ...). The neutral and ionic
curves differ substantially for metals and halides, but the independent-atom
model used in refinement keeps *bonded* atoms neutral. So we assign an ion only
to a genuine monatomic ion -- a single-atom residue (e.g. an ordered Na+, Ca2+
or Cl-) -- and leave everything else neutral.

The ionic charge is taken (in priority order) from an explicit per-atom override,
then the CCD (via the offline ``chimerax.chemcomp`` store), then a small built-in
table of the most common ions. Every candidate identifier is validated against
Clipper's own table before use, because an untabulated identifier makes Clipper
abort with a fatal error rather than fall back to the neutral element.
'''

# Most common monatomic ions, keyed by (upper-case) residue/component name.
# Used only as a fallback when the CCD store has not been populated; the CCD is
# otherwise authoritative (and disambiguates e.g. FE2 -> Fe2+ from FE -> Fe3+).
_DEFAULT_ION_CHARGE = {
    'LI': 1, 'NA': 1, 'K': 1, 'RB': 1, 'CS': 1,
    'BE': 2, 'MG': 2, 'CA': 2, 'SR': 2, 'BA': 2,
    'MN': 2, 'FE2': 2, 'FE': 3, 'FE3': 3, 'CO': 2, '3CO': 3, 'NI': 2,
    'CU1': 1, 'CU': 2, 'ZN': 2, 'CD': 2, 'HG': 2,
    'AL': 3, 'CR': 3,
    'CL': -1, 'BR': -1, 'IOD': -1, 'F': -1,
}

_valid_species = None

# Assignment signatures already announced to the log, so the per-frame live-map
# recalculation does not spam an identical message every redraw.
_logged_signatures = set()


def valid_scattering_species():
    '''Set of every element/ion identifier Clipper can compute (cached).'''
    global _valid_species
    if _valid_species is None:
        from .clipper_python import scattering_factor_names
        _valid_species = set(scattering_factor_names())
    return _valid_species


def _normalize(identifier):
    '''Reformat an element/ion string the way Clipper's AtomShapeFn does:
    first letter upper, remaining letters lower, charge digits/sign kept,
    whitespace dropped. ('FE2+' -> 'Fe2+').'''
    out = []
    nalpha = 0
    for ch in identifier:
        if ch.isalpha():
            nalpha += 1
            out.append(ch.upper() if nalpha == 1 else ch.lower())
        elif not ch.isspace():
            out.append(ch)
    return ''.join(out)


def clipper_species_from_type_symbol(type_symbol, neutral_fallback):
    '''Map a CIF ``_atom_site_type_symbol`` to a Clipper scattering-factor species.

    Small-molecule CIFs commonly declare an ionic scattering type even for bonded
    atoms (e.g. ``O2-``, ``Ca2+``) -- that is the species the structure was refined
    against. Unlike the macromolecular path (which keeps bonded atoms neutral per
    the IAM convention), this honours the CIF-declared charge directly: it returns
    the charged species when Clipper tabulates it, otherwise ``neutral_fallback``.
    '''
    if not type_symbol:
        return neutral_fallback
    valid = valid_scattering_species()
    candidate = _normalize(type_symbol)
    if candidate in valid:
        return candidate
    # Some CIFs omit the '1' for monovalent ions ("Na+" rather than the table's
    # "Na1+"); retry with it inserted.
    import re
    m = re.fullmatch(r'([A-Za-z]+)([+-])', candidate)
    if m:
        alt = '{}1{}'.format(m.group(1), m.group(2))
        if alt in valid:
            return alt
    return neutral_fallback


def _session_for(atoms):
    try:
        structures = atoms.unique_structures
        if len(structures):
            return structures[0].session
    except Exception:
        pass
    return None


def _ccd_charge(session, resname, cache):
    '''Formal charge of a *monatomic* CCD component, or None.

    Only single-atom components yield a charge -- this is what restricts ionic
    factors to genuine monatomic ions even when a polyatomic ligand happens to
    have a single modelled atom.'''
    if resname in cache:
        return cache[resname]
    charge = None
    if session is not None:
        try:
            from chimerax.chemcomp import lookup
            rec = lookup(session, resname)
            if rec is not None:
                ccd_atoms = rec[0]
                if len(ccd_atoms) == 1:
                    charge = int(ccd_atoms[0][2])
        except Exception:
            charge = None
    cache[resname] = charge
    return charge


def ionic_scattering_names(atoms):
    '''Return a list of scattering-factor identifiers, one per atom.

    Bonded atoms keep their neutral element symbol; genuine monatomic ions get
    the matching ionic identifier (e.g. ``Ca2+``, ``Cl1-``) when one exists in
    Clipper's table.
    '''
    elements = atoms.element_names.tolist()
    n = len(elements)
    if n == 0:
        return elements

    import numpy
    num_atoms = atoms.residues.num_atoms
    candidates = numpy.where(num_atoms == 1)[0]
    if not len(candidates):
        return elements

    valid = valid_scattering_species()
    session = _session_for(atoms)
    res_names = atoms.residues.names
    ccd_cache = {}
    assigned = {}   # (resname, element, species) -> count, for logging

    for i in candidates:
        i = int(i)
        element = elements[i]
        atom = atoms[i]

        # 1. Explicit per-atom override of the full identifier.
        override = getattr(atom, 'clipper_scattering_type', None)
        if override:
            candidate = _normalize(override)
            if candidate in valid and candidate != element:
                elements[i] = candidate
                key = (str(res_names[i]), element, candidate)
                assigned[key] = assigned.get(key, 0) + 1
            continue

        # 2. Charge: explicit attribute, then CCD, then built-in defaults.
        charge = getattr(atom, 'formal_charge', None)
        if charge is None:
            charge = _ccd_charge(session, res_names[i], ccd_cache)
        if charge is None:
            charge = _DEFAULT_ION_CHARGE.get(res_names[i].upper())
        if not charge:
            continue

        charge = int(charge)
        candidate = '{}{}{}'.format(element, abs(charge), '+' if charge > 0 else '-')
        if candidate in valid:
            elements[i] = candidate
            key = (str(res_names[i]), element, candidate)
            assigned[key] = assigned.get(key, 0) + 1

    _log_assignments(session, assigned)
    return elements


def _log_assignments(session, assigned):
    '''Announce, once per distinct assignment pattern, which atoms received
    ionic (rather than neutral) scattering factors. De-duplicated so the live
    map's per-frame recalculation does not repeat the message.'''
    if not assigned or session is None:
        return
    signature = frozenset(assigned.items())
    if signature in _logged_signatures:
        return
    _logged_signatures.add(signature)
    parts = ['{} {}→{} (x{})'.format(rname, el, sp, n)
             for (rname, el, sp), n in sorted(assigned.items())]
    total = sum(assigned.values())
    session.logger.info(
        'ChimeraX-Clipper: using ionic X-ray scattering factors for {} atom(s): {}'
        .format(total, ', '.join(parts)))
