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

#ifndef ERROR_H
#define ERROR_H
//ROB
#include <cstdlib>

#include <string>
#include <iostream>

/** Show message and throw exception
 * @ingroup ErrorHandling
 *
 * This macro shows the given message and exits with the error code.
 *
 * @code
 * if (...)
 *     REPORT_ERROR("Error 1");
 * @endcode
 */
#define REPORT_ERROR(ErrormMsg) throw RelionError(ErrormMsg, __FILE__, __LINE__)

/** Exception class
 * @ingroup ErrorHandling
 *
 * This is the class type for the errors thrown by the exceptions
 */

class RelionError
{
public:
    /** Error code */
    int __errno;

    /** Message shown */
    std::string msg;

    /** File produstd::cing the error */
    std::string file;

    /** Line number */
    long line;

    RelionError(const std::string& what, const std::string &fileArg, const long lineArg);
    friend std::ostream& operator<<(std::ostream& o, RelionError& XE);
};

#endif
