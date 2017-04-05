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
#include <src/autopicker.h>
#ifdef CUDA
#include <src/gpu_utils/cuda_autopicker.h>
#endif

int main(int argc, char *argv[])
{
	AutoPicker prm;

	try
    {
		prm.read(argc, argv);

		prm.initialise();

		if (prm.todo_anything)
		{
#ifdef CUDA
			if (prm.do_gpu)
			{
				std::stringstream didSs;
				didSs << "AP";
				int dev_id = prm.deviceInitialise();
				prm.cudaPicker = (void*) new AutoPickerCuda((AutoPicker*)&prm, dev_id, didSs.str().c_str() );

				((AutoPickerCuda*)prm.cudaPicker)->run();
			}
			else
#endif
				prm.run();
#ifdef TIMING
		std::cout << "timings:" << std::endl;
		prm.timer.printTimes(false);
#endif
		}
    }

    catch (RelionError XE)
    {
        //prm.usage();
        std::cerr << XE;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}

