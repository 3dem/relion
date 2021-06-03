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

#include <unistd.h>
#include <string.h>
#include "src/gui_mainwindow.h"
#include <src/args.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


int main(int argc, char *argv[])
{
	Fl::scheme("gtk+");

#define _MAX_PATH 200

	char my_dir[_MAX_PATH];
	char short_dir[49];
	char* res = getcwd(my_dir, _MAX_PATH);

	// Get last 45 characters of my_dir to fit in titlebar of window
	if (strlen(my_dir) > 45)
	{
		short_dir[0]=short_dir[1]=short_dir[2]='.';
		int j = 3;
		for (int i = strlen(my_dir)-45; i < strlen(my_dir); i++, j++)
		{
			short_dir[j] = my_dir[i];
		}
		short_dir[j] = '\0';
	}
	else
	{
		int i;
		for (i = 0; i < strlen(my_dir); i++)
			short_dir[i] = my_dir[i];
		short_dir[i] = '\0';
	}

	char titletext[256];
	snprintf(titletext, 256, "RELION-%s", g_RELION_VERSION);
#ifdef PACKAGE_VERSION
	strcat(titletext,PACKAGE_VERSION);
#endif
	strcat(titletext,": ");

	strcat (titletext, short_dir);

	try
	{
		// Fill the window
		if (checkParameter(argc, argv, "--help"))
		{
			std::cerr << " [--refresh 2]  : refresh rate in seconds" << std::endl;
			std::cerr << " [--idle 3600]  : quit GUI after this many second" << std::endl;
			std::cerr << " [--readonly]   : limited version of GUI that does not touch any files" << std::endl;
			std::cerr << " [--tomo]       : show tomography-specific GUI" << std::endl;
			std::cerr << " [--ccpem]      : use the ccpem pipeliner" << std::endl;
			std::cerr << " [--do_projdir] : Don't confirm the creation of a new project directory, just make it if it doesn't exist" << std::endl;
			std::cerr << " [--version]    : show the version of this program" << std::endl;
			exit(0);
		}
		else if (checkParameter(argc, argv, "--version"))
		{
			// Although our parser checks for --version, we do it here. Otherwise GuiMainWindow asks for a new project directory.
			PRINT_VERSION_INFO();
			exit(0);
		}
		FileName fn_pipe = getParameter(argc, argv, "--pipeline", "default");
		int _update_every_sec = textToInteger(getParameter(argc, argv, "--refresh", "2"));
		int _exit_after_sec = textToInteger(getParameter(argc, argv, "--idle", "3600"));
		bool _do_read_only = checkParameter(argc, argv, "--readonly");
		bool _do_tomo = checkParameter(argc, argv, "--tomo");
		bool _use_ccpem_pipeliner = checkParameter(argc, argv, "--ccpem");
		bool _do_projdir = checkParameter(argc, argv, "--do_projdir");
		GuiMainWindow window(GUIWIDTH, GUIHEIGHT_EXT, titletext, fn_pipe, _update_every_sec, _exit_after_sec, _do_read_only, _do_tomo, _use_ccpem_pipeliner, _do_projdir);

		// Show and run the window
		window.show();
		Fl::run();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	//return Fl::run();
	return RELION_EXIT_SUCCESS;
}
