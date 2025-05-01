// PyBind11 Python bindings for Clipper
// Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
// Note that this software makes use of modified versions of the Clipper, LibCCP4
// and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
// of these is redistributed under its own license terms.


#include "../molc.h"

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-cns.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-phs.h>

using namespace clipper;

#define DEFAULT_NEW(FNAME, CNAME) \
extern "C" EXPORT void* \
FNAME##_new() \
{ \
    try { \
        auto ptr = new CNAME(); \
        return ptr; \
    } catch (...) { \
        molc_error(); \
        return 0; \
    } \
}

#define DEFAULT_DEL(FNAME, CNAME) \
extern "C" EXPORT void \
FNAME##_delete(void* ptr) \
{ \
    try { \
        auto p = static_cast<CNAME *>(ptr); \
        delete p; \
    } catch (...) { \
        molc_error(); \
    } \
}

DEFAULT_NEW(mtzfile, CCP4MTZfile)
DEFAULT_DEL(mtzfile, CCP4MTZfile)

extern "C" EXPORT void
mtzfile_open_read(void *handler, const char *fname)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::open_read, String(fname));
}

extern "C" EXPORT void
mtzfile_close_read(void *handler)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::close_read);
}

extern "C" EXPORT void
mtzfile_open_write(void *handler, const char *fname)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::open_write, String(fname));
}

extern "C" EXPORT void
mtzfile_close_write(void *handler)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::close_write);
}
