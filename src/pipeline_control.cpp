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
#include "src/pipeline_control.h"

std::string pipeline_control_outputname = "";

int pipeline_control_relion_exit(int mode)
{

    if (pipeline_control_outputname != "")
    {
		std::ofstream  fh;
		std::string fn = pipeline_control_outputname;
		if (mode==0)
		{
			fn += RELION_JOB_EXIT_SUCCESS;
		}
		else if (mode==1)
		{
			fn += RELION_JOB_EXIT_FAILURE;
			std::cout << std::endl << " RELION version: " << g_RELION_VERSION << std::endl << " exiting with an error ..." << std::endl;
		}
		else if (mode==2)
		{
			fn += RELION_JOB_EXIT_ABORTED;
			std::cout << std::endl << " exiting with an abort ..." << std::endl;
		}
		else
		{
			std::cerr << " ERROR: undefined mode! " << std::endl;
			return 12;
		}
		fh.open(fn.c_str(), std::ios::out);
		if (!fh)
		{
			std::cerr << " ERROR: cannot touch file: " << fn << std::endl;
			return 13;
		}
		fh.close();
	}

	// Still return 0 for success, and non-zero for failure/abort as in stdlib
	return mode;

}

bool is_under_pipeline_control()
{
	return (pipeline_control_outputname != "");
}

bool pipeline_control_check_abort_job()
{
	if (pipeline_control_outputname == "")
    	return false;

	struct stat buffer;
	if (stat((pipeline_control_outputname+RELION_JOB_ABORT_NOW).c_str(), &buffer) == 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}

void pipeline_control_delete_exit_files()
{
	struct stat buffer;
	if (stat((pipeline_control_outputname+RELION_JOB_EXIT_SUCCESS).c_str(), &buffer) == 0)
	{
		remove((pipeline_control_outputname+RELION_JOB_EXIT_SUCCESS).c_str());
	}

	if (stat((pipeline_control_outputname+RELION_JOB_EXIT_FAILURE).c_str(), &buffer) == 0)
	{
		remove((pipeline_control_outputname+RELION_JOB_EXIT_FAILURE).c_str());
	}

	if (stat((pipeline_control_outputname+RELION_JOB_EXIT_ABORTED).c_str(), &buffer) == 0)
	{
		remove((pipeline_control_outputname+RELION_JOB_EXIT_ABORTED).c_str());
	}
}
