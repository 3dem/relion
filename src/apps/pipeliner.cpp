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

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include "src/pipeliner.h"
#include <src/args.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


int main(int argc, char *argv[])
{

	try
	{
		// Fill the window, but don't show it!
		FileName fn_pipe = getParameter(argc, argv, "--pipeline", "default");
		FileName fn_sched = getParameter(argc, argv, "--schedule");
		int nr_repeat = textToInteger(getParameter(argc, argv, "--repeat", "1"));
		long int minutes_wait =  textToInteger(getParameter(argc, argv, "--min_wait", "10"));
		FileName fn_jobids  = getParameter(argc, argv, "--jobids", "");


		PipeLine pipeline;
		pipeline.name = fn_pipe;
		pipeline.read(DO_LOCK);
		pipeline.write(DO_LOCK);
		pipeline.runScheduledJobs(fn_sched, fn_jobids, nr_repeat, minutes_wait);

	}

    catch (RelionError XE)
    {
        std::cerr << XE;
        exit(1);
    }

    return 0;
}
