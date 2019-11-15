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

from chimerax.core.errors import UserError

_valid_db_commands = {
    'pdb':          'rcsb',
    'pdbe_updated': 'pdbe',
    'pdbe':         'pdbe',
    'pdbj':         'pdbj'
}

_cif_sources = {
    'rcsb': 'http://files.rcsb.org/download/{}-sf.cif',
    'pdbe': 'http://www.ebi.ac.uk/pdbe/entry-files/download/r{}sf.ent',
    # Only the coordinates are updated, not the data
    'pdbj': 'https://pdbj.org/rest/downloadPDBfile?format=sf&id={}'
}

_cif_filenames = {
    'rcsb': '{}-sf.cif',
    'pdbe': 'r{}sf.cif',
    'pdbj': 'r{}sf.cif',
}

_compressed = {
    'rcsb': False,
    'pdbe': False,
    'pdbj': True
}

def fetch_structure_factors(session, pdb_id, fetch_source='rcsb', ignore_cache=False, **kw):
    '''Get a structure factor file in CIF format by PDB identifier via the Internet'''
    if len(pdb_id) != 4:
        raise UserError('PDB identifiers are 4 characters long, got "{}"'.format(pdb_id))
    if fetch_source not in _cif_sources.keys():
        fetch_source = _valid_db_commands.get(fetch_source, None)
    if fetch_source is None:
        raise UserError('Fetching structure factors is not implemented for "fromDatabase {}"! Must be one of the following: {}'.format(
            fetch_source, ', '.join(_valid_db_commands.keys())
        ))
    import os
    pdb_id = pdb_id.lower()
    save_name = _cif_filenames[fetch_source].format(pdb_id)
    url = _cif_sources[fetch_source].format(pdb_id)
    from chimerax.core.fetch import fetch_file
    filename = fetch_file(session, url, '{} structure factors'.format(pdb_id),
        save_name, 'PDB-SF',
        uncompress=_compressed[fetch_source],
        ignore_cache=ignore_cache
        )

    # Double check that a cif file was downloaded instead of an HTML error
    # message saying the ID does not exist
    with open(filename, 'r') as f:
        line = f.readline()
        if not line.startswith(('data_', '#')):
            f.close()
            os.remove(filename)
            raise UserError('Structure factors could not be retrieved! Are you '
                'sure this is an x-ray structure?')

    return filename

def fetch_structure_factors_pdbe(session, pdb_id, **kw):
    return fetch_structure_factors(session, pdb_id, fetch_source='pdbe', **kw)

def fetch_structure_factors_pdbe_updated(session, pdb_id, **kw):
    return fetch_structure_factors(session, pdb_id, fetch_source='pdbe', **kw)

def fetch_structure_factors_pdbj(session, pdb_id, **kw):
    return fetch_structure_factors(session, pdb_id, fetch_source='pdbj', **kw)

def fetch_wrapper(fetch_func):
    def _fetch(session, pdb_id, fetch_source='rcsb', ignore_cache=False,
            structure_factors=False, over_sampling=2.0, **kw):
        models, status = fetch_func(session, pdb_id, ignore_cache=ignore_cache, **kw)
        if structure_factors:
            if len(models) != 1:
                raise UserError('Structure factors can only be used with a single model!')
            m = models[0]
            sf_file = fetch_structure_factors(session, pdb_id, fetch_source=fetch_source,
                ignore_cache=ignore_cache, **kw)
            from chimerax.clipper.symmetry import get_map_mgr
            mmgr = get_map_mgr(m, create=True)
            mmgr.add_xmapset_from_mtz(sf_file, oversampling_rate = over_sampling)
            return [mmgr.crystal_mgr], status
        return models, status
    return _fetch
