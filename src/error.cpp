/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/



/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#include "src/error.h"
#ifdef __GNUC__
#include <execinfo.h>
#endif


// Object Constructor
RelionError::RelionError(const std::string &what, const std::string &fileArg, const long lineArg)
{
#ifdef __GNUC__
    const int SZ_BUF = 100;
    backtrace_buffer = new void*[SZ_BUF];
    size = backtrace(backtrace_buffer, SZ_BUF);
#endif

    msg = "ERROR: \n" + what;
    file= fileArg;
    line=lineArg;

    std::cerr << "in: " << file << ", line " << line << "\n";
	std::cerr << msg << std::endl;

}

// Show message
std::ostream& operator << (std::ostream& o, RelionError& XE)
{

#ifdef __GNUC__
    o << "=== Backtrace  ===" << std::endl;
    char **bt_symbols = backtrace_symbols(XE.backtrace_buffer, XE.size);
    for (int i = 0; i < XE.size; i++) {
        o << bt_symbols[i] << std::endl;
    }
    o << "==================" << std::endl;
    delete[] XE.backtrace_buffer;
    free(bt_symbols);
#endif

    o << XE.msg << std::endl;

    return o;
}
