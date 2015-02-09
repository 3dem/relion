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
#include "src/gui_mainwindow.h"
#include <unistd.h>
#include <string.h>


int main(int argc, char *argv[])
{
	Fl::scheme("gtk+");

#define _MAX_PATH 200

	char my_dir[_MAX_PATH];
	char short_dir[49];
	getcwd(my_dir, _MAX_PATH);

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

	char titletext[57];
	strcpy (titletext,"RELION: ");
	strcat (titletext, short_dir);
	RelionMainWindow window(GUIWIDTH, GUIHEIGHT, titletext);

    window.show(argc, argv);

    return Fl::run();

}
