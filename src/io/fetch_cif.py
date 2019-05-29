# @Author: Tristan Croll <tic20>
# @Date:   29-May-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

from chimerax.core.errors import UserError

_cif_sources = {
    'rcsb': 'http://files.rcsb.org/download/{}-sf.cif',
    'pdbe': 'http://www.ebi.ac.uk/pdbe/entry-files/download/r{}sf.ent',
    'pdbj': 'https://pdbj.org/rest/downloadPDBfile?format=sf&id={}'
}

_cif_filenames = {
    'rcsb': '{}-sf.cif',
    'pdbe': 'r{}sf.ent',
    'pdbj': 'r{}sf.ent',
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
        raise UserError('Unrecognised database name! Must be one of the following: {}'.format(
            ', '.join(_cif_sources.keys())
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
        models, status = fetch_func(session, pdb_id, fetch_source=fetch_source,
            ignore_cache=ignore_cache, **kw)
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
    return _fetch
