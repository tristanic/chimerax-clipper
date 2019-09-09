
_MIN_VALID_DETERMINANT = 5e-7 # equivalent to about 0.5A**2 B-factor

def remove_invalid_anisou(session, atoms):
    logger = session.logger
    import numpy
    from .util import anisou_determinant, anisou_determinants
    # Check all current altlocs at once
    anisou_atoms = atoms[atoms.has_aniso_u]
    determinants = anisou_determinants(anisou_atoms.aniso_u6)
    bad_aniso_atoms = anisou_atoms[determinants<_MIN_VALID_DETERMINANT]
    bad_aniso_atoms.aniso_u6=None
    # Now the other altlocs. Current API doesn't make it easy to do these in a lump
    from chimerax.atomic import Atoms
    altloc_atoms = Atoms([a for a in atoms if len(a.alt_locs)>1])
    current_altlocs = altloc_atoms.alt_locs
    bad_altloc_aniso_atoms = []
    for a, cal in zip(altloc_atoms, current_altlocs):
        for altloc in a.alt_locs:
            if altloc==cal:
                continue
            with a.suppress_alt_loc_change_notifications():
                # suppress_alt_loc_change_notifications() prevents the change
                # tracker from firing and triggering a redraw on the next frame,
                # and automatically returns the atom to the original altloc when
                # done.
                a.alt_loc = altloc
                anisou = a.aniso_u6
                if anisou is not None:
                    if anisou_determinant(anisou) < _MIN_VALID_DETERMINANT:
                        a.aniso_u6 = None
                        bad_altloc_aniso_atoms.append((a, altloc))
    if len(bad_aniso_atoms) or len(bad_altloc_aniso_atoms):
        warning_string = ("The following atoms were found to have physically "
            "impossible anisotropic B-factor (ANISOU) entries. These have been "
            "removed in favour of their isotropic B-factors.")
        if len(bad_aniso_atoms):
            warning_string += '\n'+'\n'.join(['{}{}: {}'.format(
                a.residue.chain_id, a.residue.number, a.name)
                    for a in bad_aniso_atoms])
        if len(bad_altloc_aniso_atoms):
            warning_string += '\n'+'\n'.join(['{}{}: {} altloc {}'.format(
                a.residue.chain_id, a.residue.number, a.name, altloc)
                    for (a, altloc) in bad_altloc_aniso_atoms])
        session.logger.warning(warning_string)
