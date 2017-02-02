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

#define DEVERR "\n\
This is a developer error message which you cannot fix \n\
through changing the run config. Either your data is broken or\n\
an unforseen combination of options was encountered. Please this\n\
report error, the command used and a brief description to\n\
the relion developers at \n\n github.com/3dem/relion/issues \n\n"

#define ADVERR "\n\
This error is normally only displayed when using advanced \n\
features or build-utilities for code development or benchmarking.\n\
You can ask the relion developers for assistance at \n\n github.com/3dem/relion/issues"

#define ERR_GPUID ("\
There was an issue with the GPU-ids. Either \n \t\
- you have specified a GPU index following the --gpu flag which is too high \n \t\
- relion has detected more GPUs than there is one or more node(s). \n\
Try running without ids following the --gpu flag, or specify different such indices.\n\
Remember that the numbering of GPUs start with 0!\n")

#define ERRGPUKERN ("\n\
A GPU-function failed to execute.\n\n \
If this occured at the start of a run, you might have GPUs which\n\
are incompatible with either the data or your installation of relion.\n\
If you \n\n\
\t-> INSTALLED RELION YOURSELF: if you e.g. specified -DCUDA_ARCH=50\n\
\t   and are trying ot run on a compute 3.5 GPU (-DCUDA_ARCH=3.5), \n\
\t   this may happen.\n\n\
\t-> HAVE MULTIPLE GPUS OF DIFFERNT VERISONS: relion needs GPUS with\n\
\t   at least compute 3.5. You may be trying to use a GPU older than\n\
\t   this. If you have multiple generations, try specifying --gpu <X>\n\
\t   with X=0. Then try X=1 in a new run, and so on. The numbering of\n\
\t   GPUs may not be obvious from the driver or intuition. For a list\n\
\t   of GPU compute generations, see \n\n\
\t   en.wikipedia.org/wiki/CUDA#Version_features_and_specifications\n\n\
\t-> ARE USING DOUBLE-PRECISION GPU CODE: relion was been written so\n\
\t   as to not require this, and may thus have unforeseen requirements\n\
\t   when run in this mode. If you think it is nonetheless necessary,\n\
\t   please consult the developers with this error.\n\n\
If this occurred at the middle or end of a run, it might be that\n\n\
\t-> YOUR DATA OR PARAMETERS WERE UNEXPECTED: execution on GPUs is \n\
\t   subject to many restrictions, and relion is written to work within\n\
\t   common restraints. If you have exotic data or settings, unexpected\n\
\t   configurations may occur. See also above point regarding \n\
\t   double precision.\n\
If none of the above applies, please report the error to the relion\n\
developers at    github.com/3dem/relion/issues\n\n")


#define ERRCUDACAOOM ("\n\
You ran out of memory on the GPU(s).\n\n\
Each MPI-rank running on a GPU increases the use of GPU-memory. Relion\n\
tries to distribute load over multiple GPUs to increase performance,\n\
but doing this in a general and memory-efficient way is difficult.\n\n\
1. Check the device-mapping presented at the beginning of each run,\n\
   and be particularly wary of 'device X is plit between N slaves', which \n\
   will result in a higher memory cost on GPU X. In classifications, GPU-\n\
   sharing between MPI-ranks is typically fine, whereas it will usually \n\
   cause out-of-memory during the last iteration of high-resolution refinement.\n\n\
2. If you are not GPU-sharing across MPI-slave ranks, then you might be using a\n\
   too-big box-size for the GPU memory. Currently, N-pixel particle images\n\
   will require *roughly* \n\n\
\t\t    (1.1e-8)*(N*2)^3  GB  \n\n\
   of memory (per rank) during the final iteration of refinement (using\n\
   single-precision GPU code, which is default). 450-pixel images can therefore\n\
   just about fit into a GPU with 8GB of memory, since 11*(450*2)^3 ~= 8.02\n\
   During classifications, resolution is typically lower and N is suitably\n\
   reduced, which means that memory use is much lower.\n\n\
3. If the above estimation fits onto (all of) your GPU(s), you may have \n\
   a very large number of orientations which are found as possible during\n\
   the expectation step, which results in large arrays beeing needed on the \n\
   GPU. If this is the case, you should find large (>10'000) values of \n\
   '_rlnNrOfSignificantSamples' in your _data.star output files. You can try\n\
   adding the --maxsig <P>, flag, where P is an integer limit, but you \n\
   should probably also consult expertise or re-evaluate your data and/or \n\
   input reference. Seeing large such values means relion is finding nothing\n\
   to align.\n\n\
If none of the above applies, please report the error to the relion\n\
developers at    github.com/3dem/relion/issues\n\n")

#define ERR_CANZ      ("There is an allocation on the GPU left between iterations." DEVERR)
#define ERR_CAMUX     ("A mutex could not be created for a GPU memory allocation." DEVERR)
#define ERR_STAGEMEM  ("A zero-size array was attempted to be made, which should not happen." DEVERR)
#define ERR_MDLDIM    ("The model dimension was not set properly." DEVERR)
#define ERR_MDLSET    ("The model was set twice." DEVERR)

#define ERRCTIC ("You are trying to benchmark a (CPU) section, but started timing it twice."    ADVERR)
#define ERRCTOC ("You are trying to benchmark a (CPU) section, but this section has not begun." ADVERR)
#define ERRGTIC ("You are trying to benchmark a (GPU) section, but started timing it twice."    ADVERR)
#define ERRGTOC ("You are trying to benchmark a (GPU) section, but this section has not begun." ADVERR)
#define ERRTPC  ("You are trying to benchmark a (GPU) section, but there is nothing to print."  ADVERR)

#define ERRCUFFTDIR  ("You are setting the direction of a CUFFT-transform to something other than forward/inverse" DEVERR)
#define ERRCUFFTDIRF ("You are trying to run a forward CUFFT-transform for an inverse transform" DEVERR)
#define ERRCUFFTDIRR ("You are trying to run an inverse CUFFT-transform for a forward transform" DEVERR)





#endif
