/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-May-2020
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 23-May-2020
 * @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
 * @Copyright: 2016-2019 Tristan Croll
 */
#include "../molc.h"
#include <atomstruct/Atom.h>

using Atom = atomstruct::Atom;

extern "C" EXPORT void is_nonpolar_hydrogen(void *atoms, size_t n, npy_bool* hydrophobic)
{
    Atom **aa = static_cast<Atom **>(atoms);
    try {
        for (size_t i=0; i<n; ++i)
        {
            Atom* a = aa[i];
            hydrophobic[i]=false;
            if(a->element().number()==1)
            {
                for(auto nn: a->neighbors())
                    if (nn->element().number()==6)
                    {
                        hydrophobic[i]=true;
                    }
            }
        }
    } catch (...) {
        molc_error();
    }

}
