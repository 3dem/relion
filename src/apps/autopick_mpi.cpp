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
#include <src/autopicker_mpi.h>
#ifdef _CUDA_ENABLED
#include <src/acc/cuda/cuda_autopicker.h>
#endif

int main(int argc, char *argv[])
{
	AutoPickerMpi prm;

	try
	{
		prm.read(argc, argv);

		prm.initialise(prm.getRank());

#ifdef _CUDA_ENABLED
		std::stringstream didSs;
		if (prm.do_gpu)
		{
			didSs << "APr" << prm.getRank();
			prm.deviceInitialise();
		}

		if (prm.do_gpu && !(prm.do_topaz_train || prm.do_topaz_extract))
		{
			prm.cudaPicker = (void*) new AutoPickerCuda((AutoPickerMpi*)&prm, didSs.str().c_str());
			((AutoPickerCuda*)prm.cudaPicker)->run();
		}
		else
#endif
		{
			if (prm.do_topaz_train)
			{
				// only leader trains!
				if (prm.getRank() == 0) prm.trainTopaz();
				else std::cerr << " WARNNG: rank " << prm.getRank() << " is doing nothing in training as it hasn't been parallelised ..." << std::endl;
			}
			else prm.run();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (prm.getRank() == 0 && !prm.do_topaz_train)
			prm.generatePDFLogfile();
	}

	catch (RelionError XE)
	{
		if (prm.verb > 0)
    			//prm.usage();
		std::cerr << XE;
		MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);

		return RELION_EXIT_FAILURE;
	}

        MPI_Barrier(MPI_COMM_WORLD);
	return RELION_EXIT_SUCCESS;
}

