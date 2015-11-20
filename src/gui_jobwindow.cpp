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
#include "src/gui_jobwindow.h"

#define DEBUG
// Construct continuation output filenames always in the same manner
void getContinueOutname(std::string &outputname, FileNameEntry &fn_cont)
{
	int pos_it = fn_cont.getValue().rfind("_it");
	int pos_op = fn_cont.getValue().rfind("_optimiser");
	if (pos_it < 0 || pos_op < 0)
		std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << fn_cont.getValue() << std::endl;
	int it = (int)textToFloat((fn_cont.getValue().substr(pos_it+3, 6)).c_str());
	outputname += "_ct" + floatToString(it);
}

std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies)
{
	std::vector<Node> result;

	if (dim < 2 || dim > 3)
		REPORT_ERROR("getOutputNodesRefine ERROR: invalid dim value");

	FileName fn_out;
	if (iter < 0)
	{
		// 3D auto-refine
		fn_out = outputname;
	}
	else
	{
		// 2D or 3D classification
		fn_out.compose(outputname+"_it", iter, "", 3);
	}

	// Data and model.star files
	Node node1(fn_out + "_data.star", NODE_PART_DATA);
	result.push_back(node1);
	int nodetype = (dim==2) ? NODE_2DREF : NODE_3DREF;
	Node node2(fn_out + "_model.star", nodetype);
	result.push_back(node2);

	//For 3D classification or 3D auto-refine, also use individual 3D maps as outputNodes
	if (dim == 3)
	{
    	FileName fn_tmp;
		for (int iclass = 0; iclass < K; iclass++)
    	{
    		if (nr_bodies > 1)
    			fn_tmp.compose(fn_out+"_body", iclass+1, "mrc", 3);
    		else
    			fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);

    		Node node3(fn_tmp, nodetype);
    		result.push_back(node3);
    	}
	}

	// For auto-refine: also output the half1_class001_unfil.mrc map
	if (iter < 0)
	{
		Node node4(fn_out+"_half1_class001_unfil.mrc", NODE_HALFMAP);
		result.push_back(node4);
	}
	return result;

}

RelionJobWindow::RelionJobWindow(int nr_tabs, bool _has_mpi, bool _has_thread, bool _has_run,
		int x, int y, int w, int h, const char* title) : Fl_Box(x,y,w,h,title)
{

	current_y = y;
	has_mpi = _has_mpi;
	has_thread = _has_thread;

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * my_minimum_dedicated = getenv ("RELION_MINIMUM_DEDICATED");
	minimum_nr_dedicated = (my_minimum_dedicated == NULL) ? DEFAULTMININIMUMDEDICATED : textToInteger(my_minimum_dedicated);

	char * my_allow_change_dedicated = getenv ("RELION_ALLOW_CHANGE_MINIMUM_DEDICATED");
	if (my_allow_change_dedicated == NULL)
		do_allow_change_minimum_dedicated = DEFAULTMININIMUMDEDICATED;
	else
	{
		int check_allow =  textToInteger(my_allow_change_dedicated);
		do_allow_change_minimum_dedicated = (check_allow == 0) ? false : true;
	}

	// Set up tabs
    if (nr_tabs >= 1) // there is always the running tab, which is not counted on the input nr_tabs!
    {
    	tabs = new Fl_Tabs(x, current_y, w, h - MENUHEIGHT);
    	current_y += TABHEIGHT;
    	tabs->begin();
		tab1 = new Fl_Group(x, current_y , w, h - MENUHEIGHT, "");
		tab1->end();
		tab1->color(GUI_BACKGROUND_COLOR);
		tab1->selection_color(GUI_BACKGROUND_COLOR2);
		if (nr_tabs >= 2)
		{

			tab2 = new Fl_Group(x, current_y , w, h - MENUHEIGHT, "");
			tab2->end();
			tab2->color(GUI_BACKGROUND_COLOR);
			tab2->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 3)
		{
			tab3 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab3->end();
			tab3->color(GUI_BACKGROUND_COLOR);
			tab3->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 4)
		{
			tab4 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab4->end();
			tab4->color(GUI_BACKGROUND_COLOR);
			tab4->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 5)
		{
			tab5 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab5->end();
			tab5->color(GUI_BACKGROUND_COLOR);
			tab5->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 6)
		{
			tab6 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab6->end();
			tab6->color(GUI_BACKGROUND_COLOR);
			tab6->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 7)
		{
			std::cerr << "ERROR: only 6 job-specific tabs implemented..." << std::endl;
			exit(1);
		}
		current_y += 15;
	    start_y = current_y;

	    if (_has_run)
	    {
			runtab = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			runtab->label("Running");
			setupRunTab();
			runtab->end();
			runtab->color(GUI_BACKGROUND_COLOR);
			runtab->selection_color(GUI_BACKGROUND_COLOR2);
	    }

	    tabs->end();

    }

}

void RelionJobWindow::resetHeight()
{
	current_y = start_y;
}

void RelionJobWindow::setupRunTab()
{
    resetHeight();

	if (has_mpi)
		nr_mpi.place(current_y, "Number of MPI procs:", 1, 1, 64, 1, "Number of MPI nodes to use in parallel. When set to 1, MPI will not be used.");

	if (has_thread)
	{
		nr_threads.place(current_y, "Number of threads:", 1, 1, 16, 1, "Number of shared-memory (POSIX) threads to use in parallel. \
When set to 1, no multi-threading will be used. Multi-threading is often useful in 3D refinements to have more memory. 2D class averaging often proceeds more efficiently without threads.");

		ram_per_thread.place(current_y, "Available RAM (in Gb) per thread:", 4, 1, 16, 1, "Computer memory in Gigabytes that is avaliable for each thread. This will only affect some of the warnings about required computer memory.");

	}

	// Add a little spacer
	if (has_mpi || has_thread)
		current_y += STEPY/2;

    // Set up queue groups for running tab
    queue_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    queue_group->end();

    do_queue.place(current_y, "Submit to queue?", false, "Is set to Yes, the job will be submit to a queue, otherwise \
the job will be executed locally. Note that only MPI jobs may be sent to a queue.", queue_group);

	queue_group->begin();

	queuename.place(current_y, "Queue name: ", "openmpi", "Name of the queue to which to submit the job.");

	qsub.place(current_y, "Queue submit command:", "qsub", "Name of the command used to submit scripts to the queue, e.g. qsub or bsub.\n\n\
Note that the person who installed RELION should have made a custom script for your cluster/queue setup. Check this is the case \
(or create your own script following the RELION WIKI) if you have trouble submitting jobs.");

	// Two additional options that may be set through environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 (for more flexibility)
	char * extra1_text = getenv ("RELION_QSUB_EXTRA1");
	if (extra1_text != NULL)
	{
		have_extra1 = true;
		char * extra1_default = getenv ("RELION_QSUB_EXTRA1_DEFAULT");
		char emptychar[] = "";
		if (extra1_default == NULL)
			extra1_default=emptychar;
		qsub_extra1.place(current_y, extra1_text, extra1_default, "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra1XXX will be changed by this value.");
	}
	else
		have_extra1 = false;

	char * extra2_text = getenv ("RELION_QSUB_EXTRA2");
	if (extra2_text != NULL)
	{
		have_extra2 = true;
		char * extra2_default = getenv ("RELION_QSUB_EXTRA2_DEFAULT");
		char emptychar[] = "";
		if (extra2_default == NULL)
			extra2_default=emptychar;
		qsub_extra2.place(current_y, extra2_text, extra2_default, "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra2XXX will be changed by this value.");
	}
	else
		have_extra2 = false;


	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_QSUB_TEMPLATE");
	if (default_location==NULL)
	{
		char mydefault[]=DEFAULTQSUBLOCATION;
		default_location=mydefault;
	}

	qsubscript.place(current_y, "Standard submission script:", default_location, "Script Files (*.{csh,sh,bash,script})",
"The template for your standard queue job submission script. \
Its default location may be changed by setting the environment variable RELION_QSUB_TEMPLATE. \
In the template script a number of variables will be replaced: \n \
XXXcommandXXX = relion command + arguments; \n \
XXXqueueXXX = The queue name; \n \
XXXmpinodesXXX = The number of MPI nodes; \n \
XXXthreadsXXX = The number of threads; \n \
XXXcoresXXX = XXXmpinodesXXX * XXXthreadsXXX; \n \
XXXdedicatedXXX = The minimum number of dedicated cores on each node; \n \
XXXnodesXXX = The number of requested nodes = CEIL(XXXcoresXXX / XXXdedicatedXXX); \n \
If these options are not enough for your standard jobs, you may define two extra variables: XXXextra1XXX and XXXextra2XXX \
Their help text is set by the environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 \
For example, setenv RELION_QSUB_EXTRA1 \"Max number of hours in queue\" will result in an additional (text) ein the GUI \
Any variables XXXextra1XXX in the template script will be replaced by the corresponding value.\
Likewise, default values for the extra entries can be set through environment variables RELION_QSUB_EXTRA1_DEFAULT and  RELION_QSUB_EXTRA2_DEFAULT. \
But note that (unlike all other entries in the GUI) the extra values are not remembered from one run to the other.");

	min_dedicated.place(current_y, "Minimum dedicated cores per node:", minimum_nr_dedicated, 1, 64, 1, "Minimum number of dedicated cores that need to be requested on each node. This is useful to force the queue to fill up entire nodes of a given size.");
	if (do_allow_change_minimum_dedicated)
		min_dedicated.deactivate(false);
	else
		min_dedicated.deactivate(true);

	queue_group->end();
	do_queue.cb_menu_i(); // This is to make the default effective



    // Add a little spacer
    current_y += STEPY/2;

    other_args.place(current_y, "Additional arguments:", "", "In this box command-line arguments may be provided that are not generated by the GUI. \
This may be useful for testing developmental options and/or expert use of the program. \
The command 'relion_refine' will print a list of possible options.");

}

void RelionJobWindow::toggle_new_continue(bool is_continue)
{
	std::cerr << "toggle default is_continue=" << is_continue << std::endl;
	return;
}

void RelionJobWindow::openWriteFile(std::string fn, std::ofstream &fh)
{
#ifdef DEBUG
	std::cerr << " opening " << fn << ".job for writing ... ";
#endif

	fh.open((fn+".job").c_str(), std::ios::out);
    if (!fh)
    {
    	std::cerr << "Cannot write to file: "<<fn<<std::endl;
    	exit(1);
    }

    // Write the job type
    fh << "job_type == " << type << std::endl;

    // is_continue flag
    if (is_continue)
    	fh << "is_continue == true" << std::endl;
    else
    	fh << "is_continue == false" << std::endl;
}

bool RelionJobWindow::openReadFile(std::string fn, std::ifstream &fh)
{
#ifdef DEBUG
	std::cerr << " opening " << fn << ".job for reading ... ";
#endif

	fh.open((fn+".job").c_str(), std::ios_base::in);
    if (fh.fail())
    	return false;
    else
    {
		std::string line;
		// Get job type from first line
		getline(fh, line, '\n');
		size_t idx = line.find("==");
		idx++;
		type = (int)textToFloat((line.substr(idx+1,line.length()-idx)).c_str());
		if (!(type > 0 && type <= NR_BROWSE_TABS))
			REPORT_ERROR("RelionJobWindow::openReadFile ERROR: cannot find job type in " + fn + ".job");
    	// Get is_continue from second line
		getline(fh, line, '\n');
		if (line.rfind("is_continue == true") == 0)
			is_continue = true;
		else
			is_continue = false;

		return true;
    }
}

void RelionJobWindow::closeWriteFile(std::ofstream& fh)
{
	if (has_mpi)
		nr_mpi.writeValue(fh);
	if (has_thread)
	{
		nr_threads.writeValue(fh);
		ram_per_thread.writeValue(fh);
	}
	do_queue.writeValue(fh);
	queuename.writeValue(fh);
	qsub.writeValue(fh);
	if (have_extra1)
		qsub_extra1.writeValue(fh);
	if (have_extra2)
		qsub_extra2.writeValue(fh);
	qsubscript.writeValue(fh);
	other_args.writeValue(fh);

	fh.close();
#ifdef DEBUG
	std::cerr << " done! " << std::endl;
#endif
}

void RelionJobWindow::closeReadFile(std::ifstream& fh)
{
	if (has_mpi)
		nr_mpi.readValue(fh);
	if (has_thread)
	{
		nr_threads.readValue(fh);
		ram_per_thread.readValue(fh);
	}
	do_queue.readValue(fh);
	queuename.readValue(fh);
	qsub.readValue(fh);
	if (have_extra1)
		qsub_extra1.readValue(fh);
	if (have_extra2)
		qsub_extra2.readValue(fh);
	qsubscript.readValue(fh);
	other_args.readValue(fh);

#ifdef DEBUG
	std::cerr << " done! " << std::endl;
#endif

}

void RelionJobWindow::saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands)
{
	Fl_Text_Buffer *textbuf = new Fl_Text_Buffer;

	// Open the standard job submission file
	int errno;
	if (errno = textbuf->loadfile(qsubscript.getValue().c_str()))
	    fl_alert("Error reading from file \'%s\':\n%s.", qsubscript.getValue().c_str(), strerror(errno));

	// default to a single thread
	int nmpi = nr_mpi.getValue();
	int nthr = (has_thread) ? nr_threads.getValue() : 1;
	int ncores = nr_mpi.getValue() * nthr;
	int ndedi = min_dedicated.getValue();
	float fnodes = (float)ncores / (float)ndedi;
	int nnodes = CEIL(fnodes);
	if (fmod(fnodes, 1) > 0)
	{
		std:: cout << std::endl;
		std::cout << " Warning! You're using " << nmpi << " MPI processes with " << nthr << " threads each (i.e. " << ncores << " cores), while asking for " << nnodes << " nodes with " << ndedi << " cores." << std::endl;
		std::cout << " It is more efficient to make the number of cores (i.e. mpi*threads) a multiple of the minimum number of dedicated cores per node " << std::endl;
	}

	replaceStringAll(textbuf, "XXXmpinodesXXX", floatToString(nmpi) );
	replaceStringAll(textbuf, "XXXthreadsXXX", floatToString(nthr) );
	replaceStringAll(textbuf, "XXXcoresXXX", floatToString(ncores) );
	replaceStringAll(textbuf, "XXXdedicatedXXX", floatToString(ndedi) );
	replaceStringAll(textbuf, "XXXnodesXXX", floatToString(nnodes) );
	replaceStringAll(textbuf, "XXXnameXXX", outputname);
	replaceStringAll(textbuf, "XXXerrfileXXX", outputname + ".err");
	replaceStringAll(textbuf, "XXXoutfileXXX", outputname + ".out");
	replaceStringAll(textbuf, "XXXqueueXXX", queuename.getValue() );
	if (have_extra1)
		replaceStringAll(textbuf, "XXXextra1XXX", qsub_extra1.getValue() );
	if (have_extra2)
		replaceStringAll(textbuf, "XXXextra2XXX", qsub_extra2.getValue() );

	// Get commands.size() entries with the actual command
	if (commands.size() > 1)
		appendLineString(textbuf, "XXXcommandXXX", commands.size() - 1);

	for (int icom = 0; icom < commands.size(); icom++)
		replaceStringOnce(textbuf, "XXXcommandXXX", commands[icom] );

	// Make sure the file ends with an empty line
	textbuf->append("\n");

	// Save the modified job submission script using a local name
	if (errno = textbuf->savefile(newfilename.c_str()))
	    fl_alert("Error writing to file \'%s\':\n%s.", newfilename.c_str(), strerror(errno));

}

void RelionJobWindow::changeDateNTimeInOutputname(std::string &outputname)
{
	// If the outputname contains DATENTIMENRUN, then replace that now for the current time
	int datentime = outputname.rfind("DATENTIMENRUN");
	std::cerr << " old outputname = " << outputname << std::endl;
	if (datentime < outputname.size())
	{
		time_t now = time(0);
		tm *ltm = localtime(&now);
		std::string replacestr = integerToString(ltm->tm_year%100, 2);
		replacestr+= integerToString(1 + ltm->tm_mon, 2);
		replacestr+= integerToString(1 + ltm->tm_mday, 2);
		replacestr+= integerToString(1 + ltm->tm_hour, 2);
		replacestr+= integerToString(1 + ltm->tm_min, 2);
		replacestr+= integerToString(1 + ltm->tm_sec, 2);
		outputname.replace(datentime, 1+outputname.length(), replacestr + "/run");
	}
	std::cerr << " new outputname = " << outputname << std::endl;

}


void RelionJobWindow::prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command)
{

	// Create output directory if the outname contains a "/"
	int last_slash = outputname.rfind("/");
	if (last_slash < outputname.size())
	{
		std::string dirs = outputname.substr(0, last_slash);
		std::string makedirs = "mkdir -p " + dirs;
		int res = system(makedirs.c_str());
	}

	// Prepare full mpi commands or save jobsubmission script to disc
	if (do_queue.getValue())
	{
		// Make the submission script and write it to disc
		std::string output_script = outputname + "_submit.script";
		saveJobSubmissionScript(output_script, outputname, commands);
		final_command = qsub.getValue() + " " + output_script + " &";
	}
	else
	{
		// If there are multiple commands, then join them all on a single line (final_command)
		// Also add mpirun in front of all commands if no submission via the queue is done
		std::string one_command;
		final_command = "";
		for (int icom = 0; icom < commands.size(); icom++)
		{
			if (has_mpi && nr_mpi.getValue() > 1)
				one_command = "mpirun -n " + floatToString(nr_mpi.getValue()) + " " + commands[icom] ;
			else
				one_command = commands[icom];
			final_command += one_command;
			if (icom == commands.size() - 1)
				final_command += " & "; // end by putting composite job in the background
			else
				final_command += " && "; // execute one command after the other...
		}
	}

}

/*
XXXXJobWindow::XXXXJobWindow() : RelionJobWindow(4, HAS_MPI, HAS_THREAD)
{
	tab1->begin();
	tab1->label("I/O");
	resetHeight();
	//ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	//ctf_group->end();

	tab1->end();

	// read settings if hidden file exists
	read(".gui_general.settings", is_continue);
}
void XXXXJobWindow::write(std::string fn)
{
	std::ofstream fh;
	openWriteFile(fn + ".gui_general.settings", fh);
	closeWriteFile(fh);
}
void XXXXJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void XXXXJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
}
void XXXXJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command)
{
	commands.clear();
	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_XXX_mpi`";
	else
		command="`which relion_XXX`";


	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);
	outputname = "run_ctffind";
	prepareFinalCommand(outputname, commands, final_command);
}
*/

GeneralJobWindow::GeneralJobWindow() : RelionJobWindow(1, HAS_MPI, HAS_NOT_THREAD, HAS_NOT_RUN)
{

	type = PROC_GENERAL;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	angpix.place(current_y, "Magnified pixel size (Angstrom):", 1.0, 0.1, 5.0, 0.1, "Magnified pixel size in Angstroms (preferably calculated using a calibrated magnification). \
However, when an exact calibrated pixel size is not available, one may use a preliminary one throughout the entire RELION processing workflow. This will all be internally consistent. \
Then, at the postprocessing step (when one has a better idea of the actual pixel size, e.g. through the fitting of an atomic model) one may change to a more accurate pixel size at that stage.");

	particle_diameter.place(current_y, "Particle mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");

	tab1->end();

	// read settings if hidden file exists
	read(".gui_general", is_continue);
}

void GeneralJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_general";

	std::ofstream fh;
	openWriteFile(fn, fh);

	angpix.writeValue(fh);
	particle_diameter.writeValue(fh);

	closeWriteFile(fh);

}

void GeneralJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		angpix.readValue(fh);
		particle_diameter.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void GeneralJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

}

void GeneralJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command)
{
	commands.clear();
	outputname = "run_general";
	final_command="";

}

CtffindJobWindow::CtffindJobWindow() : RelionJobWindow(3, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_CTFFIND;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	mic_names.place(current_y, "Input micrographs for CTF:", "micrographs_selected.star", "Input micrographs (*.{star,mrc})", "STAR file with the filenames of all micrographs on which to run CTFFIND, OR a unix-type wildcard to the filenames of the micrograph(s) (e.g. Micrographs/*.mrc).\
Note that the micrographs should be in a subdirectory (e.g. called Micrographs/) of the project directory, i.e. the directory from where you are launching the GUI. \
If this is not the case, then make a symbolic link inside the project directory to the directory where your micrographs are stored.");

	output_star_ctf_mics.place(current_y, "Output STAR file:", "micrographs_ctf.star", "Name of the output STAR file with all CTF information for each micrograph");

	// Add a little spacer
	current_y += STEPY/2;

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_CTFFIND_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTCTFFINDLOCATION;
		default_location=mydefault;
	}

	fn_ctffind_exe.place(current_y, "CTFFIND executable:", default_location, "*.exe", "Location of the CTFFIND executable. You can control the default of this field by setting environment variable RELION_CTFFIND_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	ctf_win.place(current_y, "Estimate CTF on window size (pix) ", -1, -16, 4096, 16, "If a positive value is given, a squared window of this size at the center of the micrograph will be used to estimate the CTF. This may be useful to exclude parts of the micrograph that are unsuitable for CTF estimation, e.g. the labels at the edge of phtographic film. \n \n The original micrograph will be used (i.e. this option will be ignored) if a negative value is given.");

	tab1->end();


	tab2->begin();
	tab2->label("Microscopy");
	resetHeight();

	cs.place(current_y, "Spherical aberration (mm):", 2, 0, 8, 0.1, "Spherical aberration of the microscope used to collect these images (in mm)");

	kv.place(current_y, "Voltage (kV):", 300, 50, 500, 10, "Voltage the microscope was operated on (in kV)");

	q0.place(current_y, "Amplitude contrast:", 0.1, 0, 0.3, 0.01, "Fraction of amplitude contrast. Often values around 10% work better than theoretically more accurate lower values...");

	dstep.place(current_y, "Physical pixel size on detector (um):", 14, 1, 32, 1, "Physical pixel size of the detector (in micrometer), e.g. Falcon is 14 um, K2 is 5 um");

	tab2->end();

	tab3->begin();
	tab3->label("CTFFIND");
	resetHeight();

	box.place(current_y, "FFT box size (pix):", 512, 64, 1024, 8, "CTFFIND's Box parameter");

	resmin.place(current_y, "Minimum resolution (A):", 30, 10, 200, 10, "CTFFIND's ResMin parameter");

	resmax.place(current_y, "Maximum resolution (A):", 5, 1, 20, 1, "CTFFIND's ResMax parameter");

	dfmin.place(current_y, "Minimum defocus value (A):", 5000, 0, 25000, 1000, "CTFFIND's dFMin parameter");

	dfmax.place(current_y, "Maximum defocus value (A):", 50000, 20000, 100000, 1000, "CTFFIND's dFMax parameter");

	dfstep.place(current_y, "Defocus step size (A):", 500, 200, 2000, 100,"CTFFIND's FStep parameter");

	dast.place(current_y, "Amount of astigmatism (A):", 100, 0, 2000, 100,"CTFFIND's dAst parameter");

	tab3->end();

	// read settings if hidden file exists
	read(".gui_ctffind", is_continue);
}

void CtffindJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_ctffind";

	std::ofstream fh;
	openWriteFile(fn, fh);

	mic_names.writeValue(fh);
	output_star_ctf_mics.writeValue(fh);
	cs.writeValue(fh);
	kv.writeValue(fh);
	q0.writeValue(fh);
	dstep.writeValue(fh);
	box.writeValue(fh);
	resmin.writeValue(fh);
	resmax.writeValue(fh);
	dfmin.writeValue(fh);
	dfmax.writeValue(fh);
	dfstep.writeValue(fh);
	dast.writeValue(fh);
	fn_ctffind_exe.writeValue(fh);
	ctf_win.writeValue(fh);

	closeWriteFile(fh);
}

void CtffindJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		mic_names.readValue(fh);
		output_star_ctf_mics.readValue(fh);
		cs.readValue(fh);
		kv.readValue(fh);
		q0.readValue(fh);
		dstep.readValue(fh);
		box.readValue(fh);
		resmin.readValue(fh);
		resmax.readValue(fh);
		dfmin.readValue(fh);
		dfmax.readValue(fh);
		dfstep.readValue(fh);
		dast.readValue(fh);
		fn_ctffind_exe.readValue(fh);
		ctf_win.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void CtffindJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	mic_names.deactivate(is_continue);
	output_star_ctf_mics.deactivate(is_continue);
	cs.deactivate(is_continue);
	kv.deactivate(is_continue);
	q0.deactivate(is_continue);
	dstep.deactivate(is_continue);
	fn_ctffind_exe.deactivate(is_continue);

	// TODO: check which log files do not have Final values and re-run on those
	// for that: modify run_ctffind wrapper program
}

void CtffindJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix)
{

	commands.clear();
	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_run_ctffind_mpi`";
	else
		command="`which relion_run_ctffind`";

	// Calculate magnification from user-specified pixel size in Angstroms
	RFLOAT magn = ROUND((dstep.getValue() * 1e-6) / (angpix * 1e-10));

	command += " --i \"" + mic_names.getValue()+"\"";
	command += " --o \"" + output_star_ctf_mics.getValue()+"\"";
	command += " --ctfWin " + floatToString(ctf_win.getValue());
	command += " --CS " + floatToString(cs.getValue());
	command += " --HT " + floatToString(kv.getValue());
	command += " --AmpCnst " + floatToString(q0.getValue());
	command += " --XMAG " + floatToString(magn);
	command += " --DStep " + floatToString(dstep.getValue());
	command += " --Box " + floatToString(box.getValue());
	command += " --ResMin " + floatToString(resmin.getValue());
	command += " --ResMax " + floatToString(resmax.getValue());
	command += " --dFMin " + floatToString(dfmin.getValue());
	command += " --dFMax " + floatToString(dfmax.getValue());
	command += " --FStep " + floatToString(dfstep.getValue());
	command += " --dAst " + floatToString(dast.getValue());
	command += " --ctffind_exe " + fn_ctffind_exe.getValue();

	if (is_continue)
		command += " --only_do_unfinished ";

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	int last_slash_out = mic_names.getValue().rfind("/");
	if (last_slash_out < mic_names.getValue().size())
	{
		// The output name contains a directory: use that one for output
		outputname = mic_names.getValue().substr(0, last_slash_out + 1) + "run_ctffind";
	}
	else
	{
		outputname = "run_ctffind";
	}
	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);


	prepareFinalCommand(outputname, commands, final_command);

}

ManualpickJobWindow::ManualpickJobWindow() : RelionJobWindow(3, HAS_NOT_MPI, HAS_NOT_THREAD)
{
	type = PROC_MANUALPICK;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input micrographs:", "micrographs_ctf.star", "Input micrographs (*.{star,mrc})", "Input STAR file (with or without CTF information), OR a unix-type wildcard with all micrographs in MRC format (in this case no CTFs can be used).");

	fn_out.place(current_y, "Output STAR file: ", "selected_micrographs_ctf.star", "Output STAR file (*.star)", "The selected micrographs will all be joined in a new STAR file (which has fewer lines than the original input one if micrographs are deselected on the GUI.");

	manualpick_rootname.place(current_y, "Picking rootname: ", "manualpick", "Rootname for the coordinate files of all manually picked particles.");

	tab1->end();
	tab2->begin();
	tab2->label("Display");
	resetHeight();

	micscale.place(current_y, "Scale for micrographs:", 0.2, 0.1, 1, 0.05, "The micrographs will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );
	sigma_contrast.place(current_y, "Sigma contrast:", 3, 0, 10, 0.5, "The micrographs will be displayed with the black value set to the average of all values MINUS this values times the standard deviation of all values in the micrograph, and the white value will be set \
to the average PLUS this value times the standard deviation. Use zero to set the minimum value in the micrograph to black, and the maximum value to white ");
	white_val.place(current_y, "White value:", 0, 0, 512, 16, "Use non-zero values to set the value of the whitest pixel in the micrograph.");
	black_val.place(current_y, "Black value:", 0, 0, 512, 16, "Use non-zero values to set the value of the blackest pixel in the micrograph.");

	current_y += STEPY/2;
	lowpass.place(current_y, "Lowpass filter (A)", 0, 0, 100, 5, "Lowpass filter that will be applied to the micrograph before displaying (zero for no filter).");

	current_y += STEPY/2;
	ctfscale.place(current_y, "Scale for CTF image:", 1, 0.1, 2, 0.1, "CTFFINDs CTF image (with the Thonrings) will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );

	tab2->end();
	tab3->begin();
	tab3->label("Colors");
	color_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	color_group->end();

	resetHeight();
	do_color.place(current_y, "Blue<>red color particles?", false, "If set to true, then the circles for each particles are coloured from red to blue (or the other way around) for a given metadatalabel. If this metadatalabel is not in the picked coordinates STAR file \
(basically only the rlnAutopickFigureOfMerit or rlnClassNumber) would be useful values there, then you may provide an additional STAR file (e.g. after classification/refinement below. Particles with values -999, or that are not in the additional STAR file will be coloured the default color: green", color_group);

	color_group->begin();
	color_label.place(current_y, "MetaDataLabel for color:", "rlnParticleSelectZScore", "The Metadata label of the value to plot from red<>blue. Useful examples might be: \n \
rlnParticleSelectZScore \n rlnClassNumber \n rlnAutopickFigureOfMerit \n rlnAngleTilt \n rlnLogLikeliContribution \n rlnMaxValueProbDistribution \n rlnNrOfSignificantSamples\n");

	fn_color.place(current_y, "STAR file with color label: ", "", "STAR file (*.star)", "The program will figure out which particles in this STAR file are on the current micrograph and color their circles according to the value in the corresponding column. \
Particles that are not in this STAR file, but present in the picked coordinates file will be colored green. If this field is left empty, then the color label (e.g. rlnAutopickFigureOfMerit) should be present in the coordinates STAR file.");

	blue_value.place(current_y, "Blue value: ", 0., 0., 4., 0.1, "The value of this entry will be blue. There will be a linear scale from blue to red, according to this value and the one given below.");
	red_value.place(current_y, "Red value: ", 2., 0., 4., 0.1, "The value of this entry will be red. There will be a linear scale from blue to red, according to this value and the one given above.");
	color_group->end();
	do_color.cb_menu_i(); // make default active

	tab3->end();

	// read settings if hidden file exists
	read(".gui_manualpick", is_continue);
}

void ManualpickJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_manualpick";

	std::ofstream fh;
	openWriteFile(fn, fh);
	fn_in.writeValue(fh);
	fn_out.writeValue(fh);
	manualpick_rootname.writeValue(fh);
	lowpass.writeValue(fh);
	micscale.writeValue(fh);
	ctfscale.writeValue(fh);
	sigma_contrast.writeValue(fh);
	white_val.writeValue(fh);
	black_val.writeValue(fh);
	do_color.writeValue(fh);
	color_label.writeValue(fh);
	fn_color.writeValue(fh);
	blue_value.writeValue(fh);
	red_value.writeValue(fh);


	closeWriteFile(fh);
}

void ManualpickJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		fn_out.readValue(fh);
		manualpick_rootname.readValue(fh);
		lowpass.readValue(fh);
		micscale.readValue(fh);
		ctfscale.readValue(fh);
		sigma_contrast.readValue(fh);
		white_val.readValue(fh);
		black_val.readValue(fh);
		do_color.readValue(fh);
		color_label.readValue(fh);
		fn_color.readValue(fh);
		blue_value.readValue(fh);
		red_value.readValue(fh);
		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void ManualpickJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
	do_queue.deactivate(true);
}

void ManualpickJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
		RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	command="`which relion_manualpick`";

	command += " --i \"" + fn_in.getValue() + "\"";
	command += " --o " + fn_out.getValue();
	command += " --pickname " + manualpick_rootname.getValue();

	command += " --scale " + floatToString(micscale.getValue());
	command += " --sigma_contrast " + floatToString(sigma_contrast.getValue());
	command += " --black " + floatToString(black_val.getValue());
	command += " --white " + floatToString(white_val.getValue());

	command += " --lowpass " + floatToString(lowpass.getValue());
	command += " --angpix " + floatToString(angpix);

	command += " --ctf_scale " + floatToString(ctfscale.getValue());

	command += " --particle_diameter " + floatToString(particle_diameter);

	if (do_color.getValue())
	{
		command += " --color_label " + color_label.getValue();
		command += " --blue " + floatToString(blue_value.getValue());
		command += " --red " + floatToString(red_value.getValue());
		if (fn_color.getValue().length() > 0)
			command += " --color_star " + fn_color.getValue();
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);

	outputname = "run_manualpick";
	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	prepareFinalCommand(outputname, commands, final_command);
}


AutopickJobWindow::AutopickJobWindow() : RelionJobWindow(3, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_AUTOPICK;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_input_autopick.place(current_y, "Input micrographs for autopick:", "micrographs_ctf.star", "Input micrographs (*.{star,mrc})", "Input STAR file (with CTF information), or a unix-type wildcard with all micrographs in MRC format (in this case no CTFs can be used).");
	autopick_rootname.place(current_y, "Autopick rootname:", "autopick", "Output coordinate files will end in rootname.star");

	tab1->end();
	tab2->begin();
	tab2->label("References");
	resetHeight();

	//set up group
	autopick_ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	autopick_ctf_group->end();

	fn_refs_autopick.place(current_y, "References:", "", "Input references (*.{star,mrc,mrcs})", "Input STAR file or MRC (stack of) image(s) with the references to be used for picking. Note that the absolute greyscale needs to be correct, so only use images created by RELION itself, e.g. by 2D class averaging or projecting a RELION reconstruction.");

	lowpass_autopick.place(current_y, "Lowpass filter references (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the references before template matching. Do NOT use very high-resolution templates to search your micrographs. The signal will be too weak at high resolution anyway, and you may find Einstein from noise....");

	psi_sampling_autopick.place(current_y, "Angular sampling (deg)", 5, 1, 30, 1, "Angular sampling in degrees for exhaustive searches of the in-plane rotations for all references.");

	// Add a little spacer
	current_y += STEPY/2;

	do_invert_refs.place(current_y, "References have inverted contrast?", true, "Set to Yes to indicate that the reference have inverted contrast with respect to the particles in the micrographs.");

	do_ctf_autopick.place(current_y, "Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n \n If set to Yes, the input micrographs can only be given as a STAR file, which should contain the CTF information for each micrograph.", autopick_ctf_group);

	autopick_ctf_group->begin();

	do_ignore_first_ctfpeak_autopick.place(current_y, "Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");

	autopick_ctf_group->end();
	do_ctf_autopick.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("autopicking");
	resetHeight();

	threshold_autopick.place(current_y, "Picking threshold:", 0.05, 0, 1., 0.01, "Use lower thresholds to pick more particles (and more junk probably)");

	mindist_autopick.place(current_y, "Minimum inter-particle distance (A):", 100, 0, 1000, 20, "Particles closer together than this distance will be consider to be a single cluster. From each cluster, only one particle will be picked.");

	maxstddevnoise_autopick.place(current_y, "Maximum stddev noise:", 1.1, 0.9, 1.5, 0.02, "This is useful to prevent picking in carbon areas, or areas with big contamination features. Peaks in areas where the background standard deviation in the normalized micrographs is higher than this value will be ignored. Useful values are probably in the range 1.0 to 1.2. Set to -1 to switch off the feature to eliminate peaks due to high background standard deviations.");

	current_y += STEPY/2;

	do_write_fom_maps.place(current_y, "Write FOM maps?", false, "If set to Yes, intermediate probability maps will be written out, which (upon reading them back in) will speed up tremendously the optimization of the threshold and inter-particle distance parameters. However, with this option, one cannot run in parallel, as disc I/O is very heavy with this option set.");

	do_read_fom_maps.place(current_y, "Read FOM maps?", false, "If written out previously, read the FOM maps back in and re-run the picking to quickly find the optimal threshold and inter-particle distance parameters");

	tab3->end();

	// read settings if hidden file exists
	read(".gui_autopick", is_continue);
}

void AutopickJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_autopick";

	std::ofstream fh;
	openWriteFile(fn, fh);

	fn_input_autopick.writeValue(fh);
	fn_refs_autopick.writeValue(fh);
	autopick_rootname.writeValue(fh);
	do_invert_refs.writeValue(fh);
	do_ctf_autopick.writeValue(fh);
	do_ignore_first_ctfpeak_autopick.writeValue(fh);
	lowpass_autopick.writeValue(fh);
	psi_sampling_autopick.writeValue(fh);
	do_write_fom_maps.writeValue(fh);
	do_read_fom_maps.writeValue(fh);
	threshold_autopick.writeValue(fh);
	mindist_autopick.writeValue(fh);
	maxstddevnoise_autopick.writeValue(fh);

	closeWriteFile(fh);
}

void AutopickJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	if (openReadFile(fn, fh))
	{

		fn_input_autopick.readValue(fh);
		fn_refs_autopick.readValue(fh);
		autopick_rootname.readValue(fh);
		do_invert_refs.readValue(fh);
		do_ctf_autopick.readValue(fh);
		do_ignore_first_ctfpeak_autopick.readValue(fh);
		lowpass_autopick.readValue(fh);
		psi_sampling_autopick.readValue(fh);
		do_write_fom_maps.readValue(fh);
		do_read_fom_maps.readValue(fh);
		threshold_autopick.readValue(fh);
		mindist_autopick.readValue(fh);
		maxstddevnoise_autopick.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void AutopickJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
	return;
}

void AutopickJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_autopick_mpi`";
	else
		command="`which relion_autopick`";

	command += " --i " + fn_input_autopick.getValue();
	command += " --o " + autopick_rootname.getValue();
	command += " --particle_diameter " + floatToString(particle_diameter);
	command += " --angpix " + floatToString(angpix);
	command += " --ref " + fn_refs_autopick.getValue();

	if (do_invert_refs.getValue())
		command += " --invert ";

	if (do_ctf_autopick.getValue())
	{
		command += " --ctf ";
		if (do_ignore_first_ctfpeak_autopick.getValue())
			command += " --ctf_intact_first_peak ";
	}
	command += " --ang " + floatToString(psi_sampling_autopick.getValue());
	command += " --lowpass " + floatToString(lowpass_autopick.getValue());

	if (do_write_fom_maps.getValue())
		command += " --write_fom_maps ";

	if (do_read_fom_maps.getValue())
		command += " --read_fom_maps ";

	command += " --threshold " + floatToString(threshold_autopick.getValue());
	command += " --min_distance " + floatToString(mindist_autopick.getValue());
	command += " --max_stddev_noise " + floatToString(maxstddevnoise_autopick.getValue());

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	outputname = autopick_rootname.getValue();

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	prepareFinalCommand(outputname, commands, final_command);

}


ExtractJobWindow::ExtractJobWindow() : RelionJobWindow(3, HAS_MPI, HAS_NOT_THREAD)
{
	type = PROC_EXTRACT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	star_mics.place(current_y,"micrograph STAR file: ", "selected_micrographs_ctf.star", "Input STAR file (*.{star})", "Filename of the STAR file that contains all micrographs from which to extract particles.");
	pick_suffix.place(current_y,"Coordinate-file suffix: ", "_autopick.star", "Suffix for all the particle coordinate files. The micrograph rootnames will be extracted from each line in the input micrograph STAR file by removing the (.mrc) extension. \
Then the coordinate filenames are formed by the micrograph rootname + this suffix. For example, a suffix of _autopick.star yields a coordinate filename of mic001_autopick.star for micrograph mic001.mrc. Likewise, a .box suffix would yield mic001.box. \n \n \
Possible formats for coordinate files are RELION-generated STAR files (.star), EMAN boxer (.box) or ASCII files (with any other extension), where each line has the X and Y coordinates of a particle, possibly preceded by a single-line non-numeric header.");
	extract_rootname.place(current_y, "Extract rootname:", "particles", "Output rootname. All particle stacks will contain this rootname, and the final particles STAR file will be called this rootname plus a .star extension. This rootname should NOT contain a directory structure. \n \n Also, when extracting movie-particles, this rootname should be THE SAME ONE as the one you used to extract the average particles. \
On disc, movie particles will be distinguished from the average particles by the movie rootname on the movie tab. If you change the extract rootname upon extration of the movies, then movie processing will NOT work! ");

	tab1->end();
	tab2->begin();
	tab2->label("extract");
	resetHeight();
	extract_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	extract_group->end();

	do_extract.place(current_y, "Extract particles from micrographs?", true, "If set to Yes, particles will be extracted from the micrographs using all selected coordinate files. \
Niko Grigorieff's program CTFFIND will be used for this.", extract_group);

	extract_group->begin();

	extract_size.place(current_y,"Particle box size :", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
	do_invert.place(current_y, "Invert contrast?", false, "If set to Yes, the contrast in the particles will be inverted.");

	// Add a little spacer
	current_y += STEPY/2;

	rescale_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	rescale_group->end();
	do_rescale.place(current_y, "Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.", rescale_group);
	rescale_group->begin();
	rescale.place(current_y, "Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number");
	rescale_group->end();
	do_rescale.cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	norm_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	norm_group->end();
	do_norm.place(current_y, "Normalize particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it.", norm_group);

	norm_group->begin();
	white_dust.place(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	black_dust.place(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	norm_group->end();
	do_norm.cb_menu_i();

	extract_group->end();
	do_extract.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("movies");
	resetHeight();
	movie_extract_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	movie_extract_group->end();

	do_movie_extract.place(current_y, "Extract from movies?", false, "If set to yes, then particles will be extracted from all frames of the MRC stacks that hold the movies.\n \
The name of the MCR stacks should be the rootname of the micrographs + '_movierootname.mrcs', where the movierootname is given below.", movie_extract_group);

	movie_extract_group->begin();

	movie_rootname.place(current_y, "Rootname of movies files:", "movie", "rootname to relate each movie to the single-frame averaged micropgraph. With a rootname of 'movie', the movie for mic001.mrc should be called mic001_movie.mrcs");

	first_movie_frame.place(current_y, "First movie frame to extract: ", 1, 1, 20, 1, "Extract from this movie frame onwards. The first frame is number 1.");

	last_movie_frame.place(current_y, "Last movie frame to extract: ", 0, 0, 64, 1, "Extract until this movie frame. Zero means: extract all frames in the movie");

	movie_extract_group->end();
	do_movie_extract.cb_menu_i();

	tab3->end();

	// read settings if hidden file exists
	read(".gui_extract", is_continue);
}


void ExtractJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_extract";

	std::ofstream fh;
	openWriteFile(fn, fh);

	// I/O
	star_mics.writeValue(fh);
	pick_suffix.writeValue(fh);
	extract_rootname.writeValue(fh);

	// extract
	do_extract.writeValue(fh);
	extract_size.writeValue(fh);
	do_rescale.writeValue(fh);
	rescale.writeValue(fh);
	do_norm.writeValue(fh);
	white_dust.writeValue(fh);
	black_dust.writeValue(fh);
	do_invert.writeValue(fh);

	// movies
	do_movie_extract.writeValue(fh);
	movie_rootname.writeValue(fh);
	first_movie_frame.writeValue(fh);
	last_movie_frame.writeValue(fh);

	closeWriteFile(fh);
}
void ExtractJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{

		// I/O
		star_mics.readValue(fh);
		pick_suffix.readValue(fh);
		extract_rootname.readValue(fh);

		// extract
		do_extract.readValue(fh);
		extract_size.readValue(fh);
		do_rescale.readValue(fh);
		rescale.readValue(fh);
		do_norm.readValue(fh);
		white_dust.readValue(fh);
		black_dust.readValue(fh);
		do_invert.readValue(fh);

		// movies
		do_movie_extract.readValue(fh);
		movie_rootname.readValue(fh);
		first_movie_frame.readValue(fh);
		last_movie_frame.readValue(fh);


		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void ExtractJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
}

void ExtractJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
		RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";

	command += " --o " + extract_rootname.getValue();
	command += " --mic_star " + star_mics.getValue();
	command += " --coord_suffix " + pick_suffix.getValue();

	if (do_extract.getValue())
	{
		command += " --extract";
		command += " --extract_size " + floatToString(extract_size.getValue());
		if (do_movie_extract.getValue())
		{
			command += " --extract_movies";
			command += " --movie_rootname " + movie_rootname.getValue();
			command += " --first_movie_frame " + floatToString(first_movie_frame.getValue());
			command += " --last_movie_frame " + floatToString(last_movie_frame.getValue());
		}

		// Operate stuff
		RFLOAT bg_radius;
		bg_radius = (particle_diameter / (2. * angpix));
		if (do_rescale.getValue())
		{
			command += " --scale " + floatToString(rescale.getValue());
			bg_radius *= rescale.getValue() / extract_size.getValue();
		}
		// Get an integer number for the bg_radius
		bg_radius = (int)(bg_radius);
		if (do_norm.getValue())
		{
			command += " --norm --bg_radius " + floatToString(bg_radius);
			command += " --white_dust " + floatToString(white_dust.getValue());
			command += " --black_dust " + floatToString(black_dust.getValue());
		}
		if (do_invert.getValue())
			command += " --invert_contrast ";

	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);
	outputname = extract_rootname.getValue();
	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	prepareFinalCommand(outputname, commands, final_command);
}

SortJobWindow::SortJobWindow() : RelionJobWindow(1, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_SORT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();
	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	input_star.place(current_y, "Input particles to be sorted:", "", "Input particles(*.{star})", "This STAR file should contain in-plane rotations, in-plane translations and a class number that were obtained by alignment (class2D/class3D or auto3D) OR auto-picking. A column called rlnParticleSelectZScore will be added to this same STAR file with the sorting result. This column can then be used in the display programs to sort the particles on.");

	// Add a little spacer
	current_y += STEPY/2;

	fn_refs.place(current_y, "References:", "", "Input references (*.{star,mrc,mrcs})", "Input STAR file or MRC (stack of) image(s) with the references to be used for sorting. These should be the same references as used to determine the class number and in-plane orientations as given in the STAR file with the input particles");
	do_ctf.place(current_y, "Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n ", ctf_group);

	ctf_group->begin();
	do_ignore_first_ctfpeak.place(current_y, "Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");
	ctf_group->end();
	do_ctf.cb_menu_i();


	tab1->end();

	// read settings if hidden file exists
	read(".gui_sort", is_continue);
}

void SortJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_sort";

	std::ofstream fh;
	openWriteFile(fn, fh);

	input_star.writeValue(fh);
	fn_refs.writeValue(fh);
	do_ctf.writeValue(fh);
	do_ignore_first_ctfpeak.writeValue(fh);

	closeWriteFile(fh);
}
void SortJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		input_star.readValue(fh);
		fn_refs.readValue(fh);
		do_ctf.readValue(fh);
		do_ignore_first_ctfpeak.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void SortJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
}

void SortJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_particle_sort_mpi`";
	else
		command="`which relion_particle_sort`";

	command += " --i " + input_star.getValue();
	command += " --ref " + fn_refs.getValue();
	command += " --angpix " + floatToString(angpix);
	command += " --particle_diameter " + floatToString(particle_diameter);

	if (do_ctf.getValue())
	{
		command += " --ctf ";
		if (do_ignore_first_ctfpeak.getValue())
			command += " --ctf_intact_first_peak ";
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);


	int last_slash = input_star.getValue().rfind("/");
	if (last_slash < input_star.getValue().size())
	{
		// The output name contains a directory: use that one for output
		outputname = input_star.getValue().substr(0, last_slash + 1) + "run_sort";
	}
	else
	{
		outputname = "run_sort";
	}

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	prepareFinalCommand(outputname, commands, final_command);
}


Class2DJobWindow::Class2DJobWindow() : RelionJobWindow(4, HAS_MPI, HAS_THREAD)
{

	type = PROC_2DCLASS;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_out.place(current_y, "Output rootname:", "Class2D/run1", "Output rootname for all files of this run. \
If this rootname contains a directory structure (e.g. 20110724/run1), the directory (20110724) will be created if it does not exist.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	// Add a little spacer
	current_y += STEPY/2;

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");


	// Add a little spacer
	current_y += STEPY/2;

	nr_classes.place(current_y, "Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");

	// Add a little spacer
	current_y += STEPY/2;




	tab1->end();

	tab2->begin();
	tab2->label("CTF");

	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	resetHeight();
	do_ctf_correction.place(current_y, "Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)", ctf_group);

	ctf_group->begin();

	ctf_phase_flipped.place(current_y, "Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");

	ctf_intact_first_peak.place(current_y, "Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	ctf_group->end();

	do_ctf_correction.cb_menu_i(); // To make default effective

	tab2->end();

	tab3->begin();
	tab3->label("Optimisation");
	resetHeight();

	nr_iter.place(current_y, "Number of iterations:", 25, 1, 50, 1, "Number of iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");

	tau_fudge.place(current_y, "Regularisation parameter T:", 2 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. \
Too small values yield too-low resolution structures; too high values result in over-estimated resolutions, mostly notable by the apparition of high-frequency noise in the references.");

	// Add a little spacer
	current_y += STEPY/2;

	do_zero_mask.place(current_y, "Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");

	// Add a little spacer
	current_y += STEPY/2;

	highres_limit.place(current_y, "Limit resolution E-step to (A): ", -1, -1, 20, 1, "If set to a positive number, then the expectation step (i.e. the alignment) will be done only including the Fourier components up to this resolution (in Angstroms). \
This is useful to prevent overfitting, as the classification runs in RELION are not to be guaranteed to be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). In particular for very difficult data sets, e.g. of very small or featureless particles, this has been shown to give much better class averages. \
In such cases, values in the range of 7-12 Angstroms have proven useful.");

	tab3->end();

	tab4->begin();
	tab4->label("Sampling");

	//set up groups
	dont_skip_align_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	dont_skip_align_group->end();

	resetHeight();

	dont_skip_align.place(current_y, "Perform image alignment?", true, "If set to No, then rather than \
performing both alignment and classification, only classification will be performed. This allows the use of very focused masks.\
This requires that the optimal orientations of all particles are already stored in the input STAR file. ", dont_skip_align_group);
	dont_skip_align_group->begin();

	psi_sampling.place(current_y, "In-plane angular sampling:", 5., 0.5, 20, 0.5, "The sampling rate for the in-plane rotation angle (psi) in degrees. \
Using fine values will slow down the program. Recommended value for most 2D refinements: 5 degrees.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");


	offset_range.place(current_y, "Offset search range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	offset_step.place(current_y, "Offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	dont_skip_align_group->end();
	dont_skip_align.cb_menu_i(); // to make default effective

	tab4->end();

	// read settings if hidden file exists
	read(".gui_class2d", is_continue);

}

void Class2DJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_class2d";

	std::ofstream fh;
	openWriteFile(fn, fh);

	// I/O
	fn_out.writeValue(fh);
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);
	nr_classes.writeValue(fh);
	do_parallel_discio.writeValue(fh);

	// CTF
	do_ctf_correction.writeValue(fh);
	ctf_phase_flipped.writeValue(fh);
	ctf_intact_first_peak.writeValue(fh);

	// Optimisation
	nr_iter.writeValue(fh);
	tau_fudge.writeValue(fh);
	do_zero_mask.writeValue(fh);
	highres_limit.writeValue(fh);

	// Sampling
	dont_skip_align.writeValue(fh);
	psi_sampling.writeValue(fh);
	offset_range.writeValue(fh);
	offset_step.writeValue(fh);

	closeWriteFile(fh);
}

void Class2DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	if (openReadFile(fn, fh))
	{

		// I/O
		fn_out.readValue(fh);
		fn_cont.readValue(fh);
		fn_img.readValue(fh);
		nr_classes.readValue(fh);
		do_parallel_discio.readValue(fh);

		// CTF
		do_ctf_correction.readValue(fh);
		ctf_phase_flipped.readValue(fh);
		ctf_intact_first_peak.readValue(fh);

		// Optimisation
		nr_iter.readValue(fh);
		tau_fudge.readValue(fh);
		do_zero_mask.readValue(fh);
		highres_limit.readValue(fh);

		// Sampling
		dont_skip_align.readValue(fh);
		psi_sampling.readValue(fh);
		offset_range.readValue(fh);
		offset_step.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void Class2DJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_cont.deactivate(!is_continue);
	fn_img.deactivate(is_continue);
	nr_classes.deactivate(is_continue);
	do_zero_mask.deactivate(is_continue);
	do_ctf_correction.deactivate(is_continue);
	ctf_phase_flipped.deactivate(is_continue);
	ctf_intact_first_peak.deactivate(is_continue);

}

void Class2DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

	// I/O
	// Save the real output name (could be with _ctX for continuation)
	// This name will also be used for the stderr and stdout outputs and the submit script and gui settings filenames
	outputname = fn_out.getValue();

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	if (is_continue)
		getContinueOutname(outputname, fn_cont);

	command += " --o " + outputname;
	pipelineOutputName = outputname;
	pipelineOutputNodes = getOutputNodesRefine(outputname, nr_iter.getValue(), nr_classes.getValue(), 2, 1);

	if (is_continue)
	{
		command += " --continue " + fn_cont.getValue();
		Node node(fn_cont.getValue(), NODE_OPTIMISER);
		pipelineInputNodes.push_back(node);
	}
	else
	{
		command += " --i " + fn_img.getValue();
		Node node(fn_img.getValue(), NODE_PART_DATA);
		pipelineInputNodes.push_back(node);
		command += " --particle_diameter " + floatToString(particle_diameter);
		command += " --angpix " + floatToString(angpix);
	}
	// Parallel disc I/O?
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";

	// CTF stuff
	if (!is_continue)
	{

		if (do_ctf_correction.getValue())
		{
			command += " --ctf ";
			if (ctf_phase_flipped.getValue())
				command += " --ctf_phase_flipped ";
			if (ctf_intact_first_peak.getValue())
				command += " --ctf_intact_first_peak ";
		}
	}

	// Optimisation
	command += " --iter " + floatToString(nr_iter.getValue());
	command += " --tau2_fudge " + floatToString(tau_fudge.getValue());
	if (!is_continue)
	{
		command += " --K " + floatToString(nr_classes.getValue());
		// Always flatten the solvent
		command += " --flatten_solvent ";
		if (do_zero_mask.getValue())
			command += " --zero_mask ";
		if (highres_limit.getValue() > 0)
			command += " --strict_highres_exp " + floatToString(highres_limit.getValue());
	}

	// Sampling
	int iover = 1;
	command += " --oversampling " + floatToString((float)iover);

	if (!dont_skip_align.getValue())
	{
		command += " --skip_align ";
	}
	else
	{
		// The sampling given in the GUI will be the oversampled one!
		command += " --psi_step " + floatToString(psi_sampling.getValue() * pow(2., iover));
		// Offset range
		command += " --offset_range " + floatToString(offset_range.getValue());
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));
	}

	// Always do norm and scale correction
	if (!is_continue)
		command += " --norm --scale ";

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());
	if (!is_continue)
		command += " --memory_per_thread " + floatToString(ram_per_thread.getValue());

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	prepareFinalCommand(outputname, commands, final_command);

}

Class3DJobWindow::Class3DJobWindow() : RelionJobWindow(5, HAS_MPI, HAS_THREAD)
{

	type = PROC_3DCLASS;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_out.place(current_y, "Output rootname:", "Class3D/run1", "Output rootname for all files of this run. \
If this rootname contains a directory structure (e.g. 20110724/run1), the directory (20110724) will be created if it does not exist.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	// Add a little spacer
	current_y += STEPY/2;
	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");

	// Add a little spacer
	current_y += STEPY/2;

	nr_classes.place(current_y, "Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");

	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();

	fn_ref.place(current_y, "Reference map:", "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");

	ref_correct_greyscale.place(current_y, "Ref. map is on absolute greyscale?", false, "Probabilities are calculated based on a Gaussian noise model, \
which contains a squared difference term between the reference and the experimental image. This has a consequence that the \
reference needs to be on the same absolute intensity grey-scale as the experimental images. \
RELION and XMIPP reconstruct maps at their absolute intensity grey-scale. \
Other packages may perform internal normalisations of the reference density, which will result in incorrect grey-scales. \
Therefore: if the map was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. \
If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the first iteration, \
and prior to the second iteration the map will be filtered again using the initial low-pass filter. \
This procedure is relatively quick and typically does not negatively affect the outcome of the subsequent MAP refinement. \
Therefore, if in doubt it is recommended to set this option to No.");

	ini_high.place(current_y, "Initial low-pass filter (A):", 60, 0, 200, 5, "It is recommended to strongly low-pass filter your initial reference map. \
If it has not yet been low-pass filtered, it may be done internally using this option. \
If set to 0, no low-pass filter will be applied to the initial reference(s).");

	// Add a little spacer
	current_y += STEPY/2;

	sym_name.place(current_y, "Symmetry:", "C1", "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	tab2->end();
	tab3->begin();
	tab3->label("CTF");

	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	resetHeight();
	do_ctf_correction.place(current_y, "Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)", ctf_group);

	ctf_group->begin();

	ctf_corrected_ref.place(current_y, "Has reference been CTF-corrected?", false, "Set this option to Yes if the reference map \
represents density that is unaffected by CTF phases and amplitudes, e.g. it was created using CTF correction (Wiener filtering) inside RELION or from a PDB. \n\n\
If set to No, then in the first iteration, the Fourier transforms of the reference projections are not multiplied by the CTFs.");

	ctf_phase_flipped.place(current_y, "Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");

	ctf_intact_first_peak.place(current_y, "Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	ctf_group->end();

	do_ctf_correction.cb_menu_i(); // To make default effective

	tab3->end();
	tab4->begin();
	tab4->label("Optimisation");
	resetHeight();

	nr_iter.place(current_y, "Number of iterations:", 25, 1, 50, 1, "Number of iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");

	tau_fudge.place(current_y, "Regularisation parameter T:", 4 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. \
Too small values yield too-low resolution structures; too high values result in over-estimated resolutions, mostly notable by the apparition of high-frequency noise in the references.");
	// Add a little spacer
	current_y += STEPY/2;

	do_zero_mask.place(current_y, "Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");

	fn_mask.place(current_y, "Reference mask (optional):", "", "Image Files (*.{spi,vol,msk,mrc})", "\
If no mask is provided, a soft spherical mask based on the particle diameter will be used.\n\
\n\
Otherwise, provide a Spider/mrc map containing a (soft) mask with the same \
dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. \
The reconstructed reference map will be multiplied by this mask.\n\
\n\
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask \
the corresponding pixels in the reconstructed map are set to the average value of these pixels. \
Thereby, for example, the higher density inside the virion may be set to a constant. \
Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. \
To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Running tab).");

	// Add a little spacer
	current_y += STEPY/2;

	highres_limit.place(current_y, "Limit resolution E-step to (A): ", -1, -1, 20, 1, "If set to a positive number, then the expectation step (i.e. the alignment) will be done only including the Fourier components up to this resolution (in Angstroms). \
This is useful to prevent overfitting, as the classification runs in RELION are not to be guaranteed to be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). In particular for very difficult data sets, e.g. of very small or featureless particles, this has been shown to give much better class averages. \
In such cases, values in the range of 7-12 Angstroms have proven useful.");

	tab4->end();

	tab5->begin();
	tab5->label("Sampling");

	//set up groups
	dont_skip_align_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	dont_skip_align_group->end();

	resetHeight();

	dont_skip_align.place(current_y, "Perform image alignment?", true, "If set to No, then rather than \
performing both alignment and classification, only classification will be performed. This allows the use of very focused masks.\
This requires that the optimal orientations of all particles are already stored in the input STAR file. ", dont_skip_align_group);
	dont_skip_align_group->begin();

	sampling.place(current_y, "Angular sampling interval:", sampling_options, &sampling_options[2], "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	offset_range.place(current_y, "Offset search range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	offset_step.place(current_y, "Offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	localsearch_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	localsearch_group->end();

	do_local_ang_searches.place(current_y, "Perform local angular searches?", false, "If set to Yes, then rather than \
performing exhaustive angular searches, local searches within the range given below will be performed. \
A prior Gaussian distribution centered at the optimal orientation in the previous iteration and \
with a stddev of 1/3 of the range given below will be enforced.", localsearch_group);
	localsearch_group->begin();

	sigma_angles.place(current_y, "Local angular search range:", 5., 0, 15, 0.1, "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.");

	localsearch_group->end();
	do_local_ang_searches.cb_menu_i(); // to make default effective

	dont_skip_align_group->end();
	dont_skip_align.cb_menu_i(); // to make default effective

	tab5->end();

	// read settings if hidden file exists
	read(".gui_class3d", is_continue);

}

void Class3DJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_class3d";

	std::ofstream fh;
	openWriteFile(fn, fh);

	// I/O
	fn_out.writeValue(fh);
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);
	nr_classes.writeValue(fh);
	do_parallel_discio.writeValue(fh);

	// Reference
	fn_ref.writeValue(fh);
	ref_correct_greyscale.writeValue(fh);
	ini_high.writeValue(fh);
	sym_name.writeValue(fh);

	// CTF
	do_ctf_correction.writeValue(fh);
	ctf_corrected_ref.writeValue(fh);
	ctf_phase_flipped.writeValue(fh);
	ctf_intact_first_peak.writeValue(fh);

	// Optimisation
	nr_iter.writeValue(fh);
	tau_fudge.writeValue(fh);
	do_zero_mask.writeValue(fh);
	fn_mask.writeValue(fh);
	highres_limit.writeValue(fh);

	// Sampling
	dont_skip_align.writeValue(fh);
	sampling.writeValue(fh);
	offset_range.writeValue(fh);
	offset_step.writeValue(fh);
	do_local_ang_searches.writeValue(fh);
	sigma_angles.writeValue(fh);

	closeWriteFile(fh);
}

void Class3DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	if (openReadFile(fn, fh))
	{

		// I/O
		fn_out.readValue(fh);
		fn_cont.readValue(fh);
		fn_img.readValue(fh);
		nr_classes.readValue(fh);
		do_parallel_discio.readValue(fh);

		// Reference
		fn_ref.readValue(fh);
		ref_correct_greyscale.readValue(fh);
		ini_high.readValue(fh);
		sym_name.readValue(fh);

		// CTF
		do_ctf_correction.readValue(fh);
		ctf_corrected_ref.readValue(fh);
		ctf_phase_flipped.readValue(fh);
		ctf_intact_first_peak.readValue(fh);

		// Optimisation
		nr_iter.readValue(fh);
		tau_fudge.readValue(fh);
		do_zero_mask.readValue(fh);
		fn_mask.readValue(fh);
		highres_limit.readValue(fh);

		// Sampling
		dont_skip_align.readValue(fh);
		sampling.readValue(fh);
		offset_range.readValue(fh);
		offset_step.readValue(fh);
		do_local_ang_searches.readValue(fh);
		sigma_angles.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void Class3DJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_cont.deactivate(!is_continue);
	fn_img.deactivate(is_continue);
	nr_classes.deactivate(is_continue);

	// Reference
	fn_ref.deactivate(is_continue);
	ref_correct_greyscale.deactivate(is_continue);
	ini_high.deactivate(is_continue);
	sym_name.deactivate(is_continue);

	//CTF
	do_ctf_correction.deactivate(is_continue);
	ctf_corrected_ref.deactivate(is_continue);
	ctf_phase_flipped.deactivate(is_continue);
	ctf_intact_first_peak.deactivate(is_continue);

	//Optimisation
	do_zero_mask.deactivate(is_continue);

	// Sampling


}

void Class3DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

	// I/O
	// Save the real output name (could be with _ctX for continuation)
	// This name will also be used for the stderr and stdout outputs and the submit script and gui settings filenames
	outputname = fn_out.getValue();
	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);
	std::cerr << " class23d outputname= " << outputname << std::endl;
	if (is_continue)
		getContinueOutname(outputname, fn_cont);

	command += " --o " + outputname;
	pipelineOutputName = outputname;
	pipelineOutputNodes = getOutputNodesRefine(outputname, nr_iter.getValue(), nr_classes.getValue(), 3, 1);

	if (is_continue)
	{
		command += " --continue " + fn_cont.getValue();
		Node node(fn_cont.getValue(), NODE_OPTIMISER);
		pipelineInputNodes.push_back(node);
	}
	else
	{
		command += " --i " + fn_img.getValue();
		Node node(fn_img.getValue(), NODE_PART_DATA);
		pipelineInputNodes.push_back(node);
		command += " --particle_diameter " + floatToString(particle_diameter);
		command += " --angpix " + floatToString(angpix);
		if (fn_ref.getValue() != "None")
		{
			command += " --ref " + fn_ref.getValue();
			Node node(fn_ref.getValue(), NODE_3DREF);
			pipelineInputNodes.push_back(node);
		}
		if (!ref_correct_greyscale.getValue() && fn_ref.getValue() != "None") // dont do firstiter_cc when giving None
			command += " --firstiter_cc";
		if (ini_high.getValue() > 0.)
			command += " --ini_high " + floatToString(ini_high.getValue());

	}
	// Parallel disc I/O?
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";

	// CTF stuff
	if (!is_continue)
	{

		if (do_ctf_correction.getValue())
		{
			command += " --ctf";
			if (ctf_corrected_ref.getValue())
				command += " --ctf_corrected_ref";
			if (ctf_phase_flipped.getValue())
				command += " --ctf_phase_flipped";
			if (ctf_intact_first_peak.getValue())
				command += " --ctf_intact_first_peak";
		}
	}

	// Optimisation
	command += " --iter " + floatToString(nr_iter.getValue());
	command += " --tau2_fudge " + floatToString(tau_fudge.getValue());
	if (!is_continue)
	{
		command += " --K " + floatToString(nr_classes.getValue());
		// Always flatten the solvent
		command += " --flatten_solvent";
		if (do_zero_mask.getValue())
			command += " --zero_mask";
		if (highres_limit.getValue() > 0)
			command += " --strict_highres_exp " + floatToString(highres_limit.getValue());
	}
	if (fn_mask.getValue().length() > 0)
	{
		command += " --solvent_mask " + fn_mask.getValue();
		Node node(fn_mask.getValue(), NODE_3DMASK);
		pipelineInputNodes.push_back(node);
	}

	// Sampling
	int iover = 1;
	command += " --oversampling " + floatToString((float)iover);

	if (!dont_skip_align.getValue())
	{
		command += " --skip_align ";
	}
	else
	{
		for (int i = 0; i < 10; i++)
		{
			if (strcmp((sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
			{
				// The sampling given in the GUI will be the oversampled one!
				command += " --healpix_order " + floatToString((float)i + 1 - iover);
				break;
			}
		}
		// Manually input local angular searches
		if (do_local_ang_searches.getValue())
			command += " --sigma_ang " + floatToString(sigma_angles.getValue() / 3.);

		// Offset range
		command += " --offset_range " + floatToString(offset_range.getValue());
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));
	}

	// Provide symmetry, and always do norm and scale correction
	if (!is_continue)
	{
		command += " --sym " + sym_name.getValue();
		command += " --norm --scale ";
	}

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());
	if (!is_continue)
		command += " --memory_per_thread " + floatToString(ram_per_thread.getValue());

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	prepareFinalCommand(outputname, commands, final_command);

}

Auto3DJobWindow::Auto3DJobWindow() : RelionJobWindow(6, HAS_MPI, HAS_THREAD)
{

	type = PROC_3DAUTO;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_out.place(current_y, "Output rootname:", "Refine3D/run1", "Output rootname for all files of this run. \
If this rootname contains a directory structure (e.g. 20110724/run1), the directory (20110724) will be created if it does not exist.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run. \n \
Besides restarting jobs that were somehow stopped before convergence, also use the continue-option after the last iteration to do movie processing.");

	// Add a little spacer
	current_y += STEPY/2;

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");


	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();

	fn_ref.place(current_y, "Reference map:", "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");

	ref_correct_greyscale.place(current_y, "Ref. map is on absolute greyscale?", false, "Probabilities are calculated based on a Gaussian noise model, \
which contains a squared difference term between the reference and the experimental image. This has a consequence that the \
reference needs to be on the same absolute intensity grey-scale as the experimental images. \
RELION and XMIPP reconstruct maps at their absolute intensity grey-scale. \
Other packages may perform internal normalisations of the reference density, which will result in incorrect grey-scales. \
Therefore: if the map was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. \
If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the first iteration, \
and prior to the second iteration the map will be filtered again using the initial low-pass filter. \
This procedure is relatively quick and typically does not negatively affect the outcome of the subsequent MAP refinement. \
Therefore, if in doubt it is recommended to set this option to No.");

	ini_high.place(current_y, "Initial low-pass filter (A):", 60, 0, 200, 5, "It is recommended to strongly low-pass filter your initial reference map. \
If it has not yet been low-pass filtered, it may be done internally using this option. \
If set to 0, no low-pass filter will be applied to the initial reference(s).");

	// Add a little spacer
	current_y += STEPY/2;


	sym_name.place(current_y, "Symmetry:", "C1", "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	tab2->end();
	tab3->begin();
	tab3->label("CTF");

	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	resetHeight();
	do_ctf_correction.place(current_y, "Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)", ctf_group);

	ctf_group->begin();

	ctf_corrected_ref.place(current_y, "Has reference been CTF-corrected?", false, "Set this option to Yes if the reference map \
represents density that is unaffected by CTF phases and amplitudes, e.g. it was created using CTF correction (Wiener filtering) inside RELION or from a PDB. \n\n\
If set to No, then in the first iteration, the Fourier transforms of the reference projections are not multiplied by the CTFs.");

	ctf_phase_flipped.place(current_y, "Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");

	ctf_intact_first_peak.place(current_y, "Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	ctf_group->end();

	do_ctf_correction.cb_menu_i(); // To make default effective

	tab3->end();
	tab4->begin();
	tab4->label("Optimisation");
	resetHeight();

	do_zero_mask.place(current_y, "Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");

	fn_mask.place(current_y, "Reference mask (optional):", "", "Image Files (*.{spi,vol,msk,mrc})", "\
If no mask is provided, a soft spherical mask based on the particle diameter will be used.\n\
\n\
Otherwise, provide a Spider/mrc map containing a (soft) mask with the same \
dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. \
The reconstructed reference map will be multiplied by this mask.\n\
\n\
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask \
the corresponding pixels in the reconstructed map are set to the average value of these pixels. \
Thereby, for example, the higher density inside the virion may be set to a constant. \
Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. \
To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Running tab).");


	tab4->end();
	tab5->begin();
	tab5->label("Auto-sampling");
	resetHeight();

	autosample_text.place(current_y, "Note that initial sampling rates will be auto-incremented!");

	sampling.place(current_y, "Initial angular sampling:", sampling_options, &sampling_options[2], "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");

	offset_range.place(current_y, "Initial offset range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");

	offset_step.place(current_y, "Initial offset step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");

	// Add a little spacer
	current_y += STEPY/2;

	auto_local_sampling.place(current_y, "Local searches from auto-sampling:", sampling_options, &sampling_options[4], "In the automated procedure to \
increase the angular samplings, local angular searches of -6/+6 times the sampling rate will be used from this angular sampling rate onwards. For most \
lower-symmetric particles a value of 1.8 degrees will be sufficient. Perhaps icosahedral symmetries may benefit from a smaller value such as 0.9 degrees.");

	tab5->end();
	tab4->end();
	tab6->begin();
	tab6->label("Movies");
	resetHeight();
	movie_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	movie_group->end();

	do_movies.place(current_y, "Realign movie frames?", false, "If set to Yes, then running averages of the individual frames of recorded movies will be aligned as independent particles.", movie_group);

	movie_group->begin();

	fn_movie_star.place(current_y, "Input movie frames STAR file:", "", "STAR Files (*.{star})", "Select the output STAR file from the preprocessing \
procedure of the movie frames.");

	movie_runavg_window.place(current_y, "Running average window:", 5, 1, 15, 1, "The individual movie frames will be averaged using a running \
average window with the specified width. Use an odd number. The optimal value will depend on the SNR in the individual movie frames. For ribosomes, we used a value of 5, where \
each movie frame integrated approximately 1 electron per squared Angstrom.");

	movie_sigma_offset.place(current_y, "Stddev on the translations (pix):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 2 pixels");

	alsorot_movie_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	alsorot_movie_group->end();

	do_alsorot_movies.place(current_y, "Also include rotational searches?", false, "If set to Yes, then running averages of the individual frames of recorded movies will also be aligned rotationally. \n \
If one wants to perform particle polishing, then rotational alignments of the movie frames is NOT necessary and will only take more computing time.", alsorot_movie_group);

	alsorot_movie_group->begin();

	movie_sigma_angles.place(current_y, "Stddev on the rotations (deg):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 1 degree");

	alsorot_movie_group->end();
	do_alsorot_movies.cb_menu_i(); // to make default effective

	movie_group->end();
	do_movies.cb_menu_i(); // to make default effective

	tab6->end();
	// read settings if hidden file exists
	read(".gui_auto3d", is_continue);

}

void Auto3DJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_auto3d";

	std::ofstream fh;
	openWriteFile(fn, fh);

	// I/O
	fn_out.writeValue(fh);
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);
	do_parallel_discio.writeValue(fh);

	// Reference
	fn_ref.writeValue(fh);
	ref_correct_greyscale.writeValue(fh);
	ini_high.writeValue(fh);
	sym_name.writeValue(fh);

	// CTF
	do_ctf_correction.writeValue(fh);
	ctf_corrected_ref.writeValue(fh);
	ctf_phase_flipped.writeValue(fh);
	ctf_intact_first_peak.writeValue(fh);

	// Optimisation
	do_zero_mask.writeValue(fh);
	fn_mask.writeValue(fh);

	// Sampling
	sampling.writeValue(fh);
	offset_range.writeValue(fh);
	offset_step.writeValue(fh);
	auto_local_sampling.writeValue(fh);

	// Movies
	do_movies.writeValue(fh);
	fn_movie_star.writeValue(fh);
	movie_runavg_window.writeValue(fh);
	movie_sigma_offset.writeValue(fh);
	do_alsorot_movies.writeValue(fh);
	movie_sigma_angles.writeValue(fh);

	closeWriteFile(fh);
}

void Auto3DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;
	if (openReadFile(fn, fh))
	{

		// I/O
		fn_out.readValue(fh);
		fn_cont.readValue(fh);
		fn_img.readValue(fh);
		do_parallel_discio.readValue(fh);

		// Reference
		fn_ref.readValue(fh);
		ref_correct_greyscale.readValue(fh);
		ini_high.readValue(fh);
		sym_name.readValue(fh);

		// CTF
		do_ctf_correction.readValue(fh);
		ctf_corrected_ref.readValue(fh);
		ctf_phase_flipped.readValue(fh);
		ctf_intact_first_peak.readValue(fh);

		// Optimisation
		do_zero_mask.readValue(fh);
		fn_mask.readValue(fh);

		// Sampling
		sampling.readValue(fh);
		offset_range.readValue(fh);
		offset_step.readValue(fh);
		auto_local_sampling.readValue(fh);

		// Movies
		do_movies.readValue(fh);
		fn_movie_star.readValue(fh);
		movie_runavg_window.readValue(fh);
		movie_sigma_offset.readValue(fh);
		do_alsorot_movies.readValue(fh);
		movie_sigma_angles.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void Auto3DJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_cont.deactivate(!is_continue);
	fn_img.deactivate(is_continue);

	// Reference
	fn_ref.deactivate(is_continue);
	ref_correct_greyscale.deactivate(is_continue);
	ini_high.deactivate(is_continue);
	sym_name.deactivate(is_continue);

	//CTF
	do_ctf_correction.deactivate(is_continue);
	ctf_corrected_ref.deactivate(is_continue);
	ctf_phase_flipped.deactivate(is_continue);
	ctf_intact_first_peak.deactivate(is_continue);

	//Optimisation
	do_zero_mask.deactivate(is_continue);

	// Sampling
	sampling.deactivate(is_continue);
	offset_range.deactivate(is_continue);
	offset_step.deactivate(is_continue);
	auto_local_sampling.deactivate(is_continue);

	// Movies
	do_movies.deactivate(!is_continue);
	fn_movie_star.deactivate(!is_continue);
	movie_runavg_window.deactivate(!is_continue);
	movie_sigma_offset.deactivate(!is_continue);
	do_alsorot_movies.deactivate(!is_continue);
	movie_sigma_angles.deactivate(!is_continue);


}

void Auto3DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter)
{

	commands.clear();
	std::string command;
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

	// I/O
	// Save the real output name (could be with _ctX for continuation)
	// This name will also be used for the stderr and stdout outputs and the submit script and gui settings filenames
	outputname = fn_out.getValue();

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	if (is_continue)
		getContinueOutname(outputname, fn_cont);
	command += " --o " + outputname;
	pipelineOutputName = outputname;
	pipelineOutputNodes = getOutputNodesRefine(outputname, -1, 1, 3, 1); // TODO: add nr_bodies....

	if (is_continue)
	{
		command += " --continue " + fn_cont.getValue();
		Node node(fn_cont.getValue(), NODE_OPTIMISER);
		pipelineInputNodes.push_back(node);
	}
	else
	{
		command += " --auto_refine --split_random_halves --i " + fn_img.getValue();
		Node node(fn_img.getValue(), NODE_PART_DATA);
		pipelineInputNodes.push_back(node);
		command += " --particle_diameter " + floatToString(particle_diameter);
		command += " --angpix " + floatToString(angpix);
		if (fn_ref.getValue() != "None")
		{
			command += " --ref " + fn_ref.getValue();
			Node node(fn_ref.getValue(), NODE_3DREF);
			pipelineInputNodes.push_back(node);
		}
		if (!ref_correct_greyscale.getValue() && fn_ref.getValue() != "None") // dont do firstiter_cc when giving None
			command += " --firstiter_cc";
		if (ini_high.getValue() > 0.)
			command += " --ini_high " + floatToString(ini_high.getValue());

	}
	// Parallel disc I/O?
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";

	// CTF stuff
	if (!is_continue)
	{

		if (do_ctf_correction.getValue())
		{
			command += " --ctf";
			if (ctf_corrected_ref.getValue())
				command += " --ctf_corrected_ref";
			if (ctf_phase_flipped.getValue())
				command += " --ctf_phase_flipped";
			if (ctf_intact_first_peak.getValue())
				command += " --ctf_intact_first_peak";
		}
	}

	// Optimisation
	if (!is_continue)
	{
		// Always flatten the solvent
		command += " --flatten_solvent";
		if (do_zero_mask.getValue())
			command += " --zero_mask";
	}
	if (fn_mask.getValue().length() > 0)
	{
		command += " --solvent_mask " + fn_mask.getValue();

		// TODO: what if this is a continuation run: re-put the mask as an input node? Or only if it changes? Also for 3Dclass
		Node node(fn_mask.getValue(), NODE_3DMASK);
		pipelineInputNodes.push_back(node);
	}

	// Sampling
	int iover = 1;
	command += " --oversampling " + floatToString((float)iover);

	for (int i = 0; i < 10; i++)
	{
		if (strcmp((sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
		{
			// The sampling given in the GUI will be the oversampled one!
			command += " --healpix_order " + floatToString((float)i + 1 - iover);
			break;
		}

	}
	// Minimum sampling rate to perform local searches (may be changed upon continuation
	for (int i = 0; i < 10; i++)
	{
		if (strcmp((auto_local_sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
		{
			command += " --auto_local_healpix_order " + floatToString((float)i + 1 - iover);
			break;
		}
	}

	// Offset range
	command += " --offset_range " + floatToString(offset_range.getValue());
	// The sampling given in the GUI will be the oversampled one!
	command += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));

	// Provide symmetry, and always do norm and scale correction
	if (!is_continue)
	{

		command += " --sym " + sym_name.getValue();
                // Always join low-res data, as some D&I point group refinements may fall into different hands!
                command += " --low_resol_join_halves 40";
		command += " --norm --scale ";
	}

	// Movies
	if (is_continue && do_movies.getValue())
	{
		command += " --realign_movie_frames " + fn_movie_star.getValue();
		Node node(fn_movie_star.getValue(), NODE_MOVIE_DATA);
		pipelineInputNodes.push_back(node);
		command += " --movie_frames_running_avg " + floatToString(movie_runavg_window.getValue());
		command += " --sigma_off " + floatToString(movie_sigma_offset.getValue());

		if (do_alsorot_movies.getValue())
		{
			command += " --sigma_ang " + floatToString(movie_sigma_angles.getValue());
		}
		else
		{
			command += " --skip_rotate --skip_maximize ";
		}
	}

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());
	if (!is_continue)
		command += " --memory_per_thread " + floatToString(ram_per_thread.getValue());

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	prepareFinalCommand(outputname, commands, final_command);

}

PostJobWindow::PostJobWindow() : RelionJobWindow(4, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_POST;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();
	fn_in.place(current_y, "One of the 2 unfiltered half-maps:", "", "MRC map files (*half1_class001_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");

	fn_out.place(current_y, "Output rootname", "postprocess", "Output rootname. All output files will be saved in the same directory as the unfiltered maps, unless the output name contains a forward slash. In that case, the corresponding directory will be created .");

	//ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	//ctf_group->end();

	tab1->end();
	tab2->begin();
	tab2->label("Mask");
	resetHeight();
	automask_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	automask_group->end();

	do_automask.place(current_y, "Perform auto-masking?", true, "If set to Yes, the program will automatically calculate a mask around the reconstructed particle. \
 A nicely fitting mask will result in an optimal resolution-estimation, as a maximal amount of noise will be removed from the resolution-estimation process. \
The effect of the masking on the FSC-curve will be measured using the randomised phase-approach as in Shaoxia Chen et al, (2013) Ultramicroscopy", automask_group);

	automask_group->begin();

	inimask_threshold.place(current_y, "Initial binarisation threshold:", 0.02, 0., 0.5, 0.01, "This threshold is used to make an initial binary mask from the average of the two unfiltered half-reconstructions. \
If you don't know what value to use, display one of the unfiltered half-maps in a 3D surface rendering viewer and find the lowest threshold that gives no noise peaks outside the reconstruction.");
	extend_inimask.place(current_y, "Extend binary map this many pixels:", 3, 0, 20, 1, "The initial binary mask is extended this number of pixels in all directions." );
	width_mask_edge.place(current_y, "Add a soft-edge of this many pixels:", 3, 0, 20, 1, "The extended binary mask is further extended with a raised-cosine soft edge of the specified width." );
	automask_group->end();
	do_automask.cb_menu_i();

	usermask_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	usermask_group->end();

	do_usermask.place(current_y, "Provide your own mask?", false, "If set to Yes, the program will mask the unfiltered half-reconstructions using a user-provided mask below. This also allows you to save execution time, by providing a previously determined automask for the same particle.", usermask_group);
	usermask_group->begin();
	fn_mask.place(current_y, "User-provided mask:", "", "Image Files (*.{spi,vol,msk,mrc})", "Use this to skip auto-masking by providing your own mask. You may also save execution time by providing a previously determined automask for the same particle.");
	usermask_group->end();
	do_usermask.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Sharpen");
	resetHeight();
	//ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	//ctf_group->end();

	fn_mtf.place(current_y, "MTF of the detector (STAR file)", "", "STAR Files (*.star)", "The MTF of the detector is used in the (later) post-processing and particle polishing stages of refinement.  \
If you know the MTF of your detector, provide it here. Curves for some well-known detectors may be downloaded from the RELION Wiki. Also see there for the exact format \
\n If you do not know the MTF of your detector and do not want to measure it, then by leaving this entry empty, you include the MTF of your detector in your overall estimated B-factor upon sharpening the map.\
Although that is probably slightly less accurate, the overall quality of your map will probably not suffer very much.");

	current_y += STEPY/2;

	autobfac_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	autobfac_group->end();
	do_auto_bfac.place(current_y,"Estimate B-factor automatically?", true, "If set to Yes, then the program will use the automated procedure described by Rosenthal and Henderson (2003, JMB) to estimate an overall B-factor for your map, and sharpen it accordingly. \
Note that your map must extend well beyond the lowest resolution included in the procedure below, which should not be set to resolutions much lower than 10 Angstroms. ", autobfac_group);

	autobfac_group->begin();
	autob_lowres.place(current_y,"Lowest resolution for auto-B fit (A):", 10, 8, 15, 0.5, "This is the lowest frequency (in Angstroms) that will be included in the linear fit of the Guinier plot as described in Rosenthal and Henderson (2003, JMB). Dont use values much lower or higher than 10 Angstroms. If your map does not extend beyond 10 Angstroms, then instead of the automated procedure use your own B-factor.");
	autobfac_group->end();
	do_auto_bfac.cb_menu_i();

	adhocbfac_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	adhocbfac_group->end();
	do_adhoc_bfac.place(current_y,"Use your own B-factor?", false, "Instead of using the automated B-factor estimation, provide your own value. Use negative values for sharpening the map. \
This option is useful if your map does not extend beyond the 10A needed for the automated procedure, or when the automated procedure does not give a suitable value (e.g. in more disordered parts of the map).",adhocbfac_group);

	adhocbfac_group->begin();
	adhoc_bfac.place(current_y,"User-provided B-factor:", -1000, -2000, 0, -50, "Use negative values for sharpening. Be careful: if you over-sharpen your map, you may end up interpreting noise for signal!");
	adhocbfac_group->end();
	do_adhoc_bfac.cb_menu_i();

	tab3->end();
	tab4->begin();
	tab4->label("Filter");
	resetHeight();
	skipweight_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	skipweight_group->end();

	do_skip_fsc_weighting.place(current_y, "Skip FSC-weighting?", false, "If set to No (the default), then the output map will be low-pass filtered according to the mask-corrected, gold-standard FSC-curve. \
Sometimes, it is also useful to provide an ad-hoc low-pass filter (option below), as due to local resolution variations some parts of the map may be better and other parts may be worse than the overall resolution as measured by the FSC. \
In such cases, set this option to Yes and provide an ad-hoc filter as described below.", skipweight_group);

	skipweight_group->begin();
	low_pass.place(current_y, "Ad-hoc low-pass filter (A):",5,1,40,1,"This option allows one to low-pass filter the map at a user-provided frequency (in Angstroms). When using a resolution that is higher than the gold-standard FSC-reported resolution, take care not to interpret noise in the map for signal...");
	skipweight_group->end();
	do_skip_fsc_weighting.cb_menu_i();

	tab4->end();


	// read settings if hidden file exists
	read(".gui_post", is_continue);
}

void PostJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_post";

	std::ofstream fh;
	openWriteFile(fn, fh);
	fn_in.writeValue(fh);
	fn_out.writeValue(fh);
	do_automask.writeValue(fh);
	inimask_threshold.writeValue(fh);
	extend_inimask.writeValue(fh);
	width_mask_edge.writeValue(fh);
	do_usermask.writeValue(fh);
	fn_mask.writeValue(fh);
	do_auto_bfac.writeValue(fh);
	autob_lowres.writeValue(fh);
	do_adhoc_bfac.writeValue(fh);
	adhoc_bfac.writeValue(fh);
	fn_mtf.writeValue(fh);
	do_skip_fsc_weighting.writeValue(fh);
	low_pass.writeValue(fh);
	closeWriteFile(fh);
}

void PostJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		fn_out.readValue(fh);
		do_automask.readValue(fh);
		inimask_threshold.readValue(fh);
		extend_inimask.readValue(fh);
		width_mask_edge.readValue(fh);
		do_usermask.readValue(fh);
		fn_mask.readValue(fh);
		do_auto_bfac.readValue(fh);
		autob_lowres.readValue(fh);
		do_adhoc_bfac.readValue(fh);
		adhoc_bfac.readValue(fh);
		fn_mtf.readValue(fh);
		do_skip_fsc_weighting.readValue(fh);
		low_pass.readValue(fh);
		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void PostJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
}

void PostJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
		RFLOAT angpix)
{

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	commands.clear();
	std::string command;
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	command="`which relion_postprocess`";

	// Get the input rootname from the half-map name
	// run1_half1_class001_unfil.mrc -> run1
	Node node(fn_in.getValue(), NODE_HALFMAP);
	pipelineInputNodes.push_back(node);
	int pos_half = fn_in.getValue().rfind("_half");
	if (pos_half < fn_in.getValue().size())
		command += " --i " + fn_in.getValue().substr(0, pos_half);
	else
	{
		std::cerr << "PostJobWindow::getCommands ERROR: cannot find _half substring in input filename: " << fn_in.getValue() << std::endl;
		exit(1);
	}

	// Get the output name. If the specified fn_out contains a directory, then make this directory and use it for output
	// If it doesn't contain a directory, then use the directory as in fn_in.
	int last_slash_out = fn_out.getValue().rfind("/");
	int last_slash_in = fn_in.getValue().rfind("/");
	if (last_slash_out < fn_out.getValue().size())
	{
		// The output name contains a directory: use that one for output
		outputname = fn_out.getValue();
	}
	else if  (last_slash_in < fn_in.getValue().size())
	{
		// Otherwise: the input name contains a directory: use that one for output
		std::string dirs = fn_in.getValue().substr(0, last_slash_in + 1);
		outputname = dirs + fn_out.getValue();
	}
	// The output name contains a directory: use it for output
	command += " --o " + outputname;
	command += "  --angpix " + floatToString(angpix);

	pipelineOutputName = outputname;
	Node node1(outputname+".mrc", NODE_FINALMAP);
	pipelineOutputNodes.push_back(node1);
	Node node2(outputname+"_masked.mrc", NODE_FINALMAP);
	pipelineOutputNodes.push_back(node2);

	// Masking
	if (do_automask.getValue())
	{
		command += " --auto_mask ";
		command += " --inimask_threshold " + floatToString(inimask_threshold.getValue());
		command += " --extend_inimask " + floatToString(extend_inimask.getValue());
		command += " --width_mask_edge " + floatToString(width_mask_edge.getValue());

		Node(outputname + "_automask.mrc", NODE_3DMASK);
		pipelineOutputNodes.push_back(node);
	}

	/// TODO: make a CREATE_MASK process and get it out of postprocessing
	if (do_usermask.getValue())
	{
		command += " --mask " + fn_mask.getValue();
		Node node(fn_mask.getValue(), NODE_3DMASK);
		pipelineInputNodes.push_back(node);
	}

	// Sharpening
	if (fn_mtf.getValue().length() > 0)
	{
		command += " --mtf " + fn_mtf.getValue();
	}
	if (do_auto_bfac.getValue())
	{
		command += " --auto_bfac ";
		command += " --autob_lowres "  + floatToString(autob_lowres.getValue());
	}
	if (do_adhoc_bfac.getValue())
	{
		command += " --adhoc_bfac " + floatToString(adhoc_bfac.getValue());
	}

	// Filtering
	if (do_skip_fsc_weighting.getValue())
	{
		command += " --skip_fsc_weighting ";
		command += " --low_pass "  + floatToString(low_pass.getValue());
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);
	prepareFinalCommand(outputname, commands, final_command);
}


PolishJobWindow::PolishJobWindow() : RelionJobWindow(3, HAS_MPI, HAS_THREAD)
{

	type = PROC_POLISH;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input STAR file with aligned movies:", "", "STAR files (*_data.star)",  "Provide the data.star file that was output by the movie-processing option in the auto-refine job.");

	fn_out.place(current_y, "Output rootname", "shiny", "Output rootname. Note that the program will write out ");

	tab1->end();
	tab2->begin();
	tab2->label("Movement");
	resetHeight();

	fit_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	fit_group->end();


	do_fit_movement.place(current_y, "Linear fit particle movements?", true, "If set to Yes, then the program will fit linear tracks (in X-Y and in time) through \
the estimated movie tracks in the input STAR file. For small particles (e.g. < 1MDa) this will lead to more robust beam-induced movement modelling. Because particles that are close to each other on a \
micrograph often move in similar directions, the estimated tracks from neighbouring particles may also be included in the fitting of each particle. Again, in particular for smaller particles \
this may improve the robustness of the fits.", fit_group);

	fit_group->begin();

	movie_runavg_window.place(current_y, "Running average window:", 5, 1, 15, 1, "Provide the same value as the one that was used to estimate the movement tracks in the movie-processing tab of the auto-refine job.");

	sigma_nb.place(current_y, "Stddev on particle distance (pix)", 100, 0, 1000, 50, "This value determines how much neighbouring particles contribute to the fit of the movements of each particle. \
This value is the standard deviation of a Gaussian on the inter-particle distance. Larger values mean that particles that are further away still contribute more. Particles beyond 3 standard deviations are excluded \
from the fit. Very large values will lead to all fitted tracks pointing in the same direction. A value of zero means that each particle is fitted independently.");
	fit_group->end();
	do_fit_movement.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Damage");
	resetHeight();

	weight_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	weight_group->end();

	do_bfactor_weighting.place(current_y, "Perform B-factor weighting?", true, "If set to Yes, then the program will estimate a resolution-dependent weight for each movie frames by calculating independent half-reconstructions for each movie frame separately. \
Gold-standard FSCs between these are then converted into relative Guinier plots, through which straight lines are fitted. Linear fits are often suitable beyond 20A resolution. Small particles may not yield such high resolutions for the individual-frame reconstructions. \
Therefore, in some cases it may be better to skip this step. It is however recommended to always try and perform B-factor weighting, and to inspect the output bfactors.star and guinier.star files, as adequate weighting may significantly improve resolution in the final map.", weight_group);

	weight_group->begin();
	perframe_highres.place(current_y, "Highres-limit per-frame maps (A)", 6, 1, 25, 1, "To estimate the resolution and dose dependency of the radiation damage, the program \
will calculate reconstructions from all first movie frames, second movie frames, etc. These per-frame reconstructions will have lower resolution than the reconstruction from all-frames. \
To speed up the calculations (and reduce memory requirements), the per-frame reconstructions may be limited in resolution using this parameter. One can inspect the output STAR files of the per-frame reconstructions \
to check afterwards that this value was not chosen lower than the actual resolution of these reconstructions");

	perframe_bfac_lowres.place(current_y, "Lowres-limit B-factor estimation (A)", 20 , 1, 40, 1, "This value describes the lowest resolution that is included in the B-factor estimation of the per-frame reconstructions. \
Because the power spectrum of per-frame reconstructions is compared to the power spectrum of the reconstruction from all frames, a much lower value than the 10A described in the Rosenthal and Henderson (2003) paper in JMB can be used. Probably a value around 20A is still OK.");

	weight_group->end();
	do_bfactor_weighting.cb_menu_i();

	current_y += STEPY/2;

	fn_mask.place(current_y, "Mask for the reconstructions", "", "Image Files (*.{spi,vol,msk,mrc})", "A continuous mask with values between 0 (solvent) and 1 (protein). You may provide the same map that was obtained in the post-processing of the corresponding auto-refine jobs before the movie processing.");

	sym_name.place(current_y, "Symmetry:", "C1", "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	tab3->end();
	// read settings if hidden file exists
	read(".gui_polish", is_continue);
}

void PolishJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_polish";

	std::ofstream fh;
	openWriteFile(fn, fh);
	fn_in.writeValue(fh);
	fn_out.writeValue(fh);
	movie_runavg_window.writeValue(fh);
	do_fit_movement.writeValue(fh);
	sigma_nb.writeValue(fh);
	do_bfactor_weighting.writeValue(fh);
	perframe_highres.writeValue(fh);
	perframe_bfac_lowres.writeValue(fh);
	fn_mask.writeValue(fh);
	sym_name.writeValue(fh);

	closeWriteFile(fh);
}
void PolishJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		fn_out.readValue(fh);
		movie_runavg_window.readValue(fh);
		do_fit_movement.readValue(fh);
		sigma_nb.readValue(fh);
		do_bfactor_weighting.readValue(fh);
		perframe_highres.readValue(fh);
		perframe_bfac_lowres.readValue(fh);
		fn_mask.readValue(fh);
		sym_name.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void PolishJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
}

void PolishJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
		RFLOAT angpix, RFLOAT particle_diameter, RFLOAT black_dust, RFLOAT white_dust)
{
	commands.clear();
	std::string command;
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	if (nr_mpi.getValue() > 1)
		command="`which relion_particle_polish_mpi`";
	else
		command="`which relion_particle_polish`";

	// General
	command += " --i " + fn_in.getValue();

	Node node(fn_in.getValue(), NODE_MOVIE_DATA);
	pipelineInputNodes.push_back(node);

	command += " --o " + fn_out.getValue();

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	Node node1(fn_out.getValue() + ".star", NODE_MOVIE_DATA); // TODO: output from Shiny in its own directory
	pipelineOutputNodes.push_back(node1);
	pipelineOutputName = fn_out.getValue();  // TODO: change OUTPUTNAME to its own output directory!

	command += "  --angpix " + floatToString(angpix);
	command += " --movie_frames_running_avg " + floatToString(movie_runavg_window.getValue());
	// If this is not a continue job, then re-start from scratch....
	if (!is_continue)
		command += " --dont_read_old_files ";

	// Beam-induced movement fitting options
	if (do_fit_movement.getValue())
		command += " --sigma_nb " + floatToString(sigma_nb.getValue());
	else
		command += " --no_fit ";

	// Damage
	if (do_bfactor_weighting.getValue())
	{
		command += " --perframe_highres " + floatToString(perframe_highres.getValue());
		command += " --autob_lowres " + floatToString(perframe_bfac_lowres.getValue());
	}
	else
	{
		command += " --skip_bfactor_weighting ";
	}

	if (fn_mask.getValue().length() > 0)
	{
		Node node(fn_mask.getValue(), NODE_3DMASK);
		pipelineInputNodes.push_back(node);
		command += " --mask " + fn_mask.getValue();
	}

	// Symmetry group
	command += " --sym " + sym_name.getValue();

	// Normalisation
	RFLOAT bg_rad = ROUND(particle_diameter / (2. * angpix));
	command += " --bg_radius " + floatToString(bg_rad);
	command += " --white_dust " + floatToString(white_dust);
	command += " --black_dust " + floatToString(black_dust);

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);

	// Place the job submission script and the stdout and stderr files in the same directory as the input data.star

	int last_slash = fn_in.getValue().rfind("/");
	if (last_slash < fn_in.getValue().size())
	{
		std::string dirs = fn_in.getValue().substr(0, last_slash + 1); // +1 includes the slash
		outputname = dirs + fn_out.getValue();
	}
	else
	{
		outputname = fn_out.getValue();
	}
	prepareFinalCommand(outputname, commands, final_command);
}

ResmapJobWindow::ResmapJobWindow() : RelionJobWindow(1, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_RESMAP;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "One of the 2 unfiltered half-maps:", "", "MRC map files (*_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");

	current_y += STEPY /2 ;

	pval.place(current_y, "P-value:", 0.05, 0., 1., 0.01, "This value is typically left at 0.05. If you change it, report the modified value in your paper!");
	minres.place(current_y, "Highest resolution (A): ", 0., 0., 10., 0.1, "ResMaps minRes parameter. By default (0), the program will start at just above 2x the pixel size");
	maxres.place(current_y, "Lowest resolution (A): ", 0., 0., 10., 0.1, "ResMaps maxRes parameter. By default (0), the program will stop at 4x the pixel size");
	stepres.place(current_y, "Resolution step size (A)", 1., 0.1, 3, 0.1, "ResMaps stepSize parameter." );

	current_y += STEPY /2 ;
	fn_mask.place(current_y, "User-provided mask (optional):", "", "Image Files (*.{spi,vol,msk,mrc})", "A binary (!) mask with values 0 for solvent and 1 for protein. \
Note that values larger than zero will be changed to 1 by ResMap, therefore the continuous masks from the postprocessing may be too wide. If left empty (default), then ResMap will determine its own mask");

	current_y += STEPY /2 ;

	// Check for environment variable RELION_RESMAP_TEMPLATE
	char * default_location = getenv ("RELION_RESMAP_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[] = DEFAULTRESMAPLOCATION;
		default_location = mydefault;
	}

	fn_resmap.place(current_y, "ResMap executable:", default_location, "ResMap*", "Location of the ResMap executable. You can control the default of this field by setting environment variable RELION_RESMAP_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	tab1->end();

	// read settings if hidden file exists
	read(".gui_resmap", is_continue);
}

void ResmapJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_resmap";

	std::ofstream fh;
	openWriteFile(fn, fh);
	fn_resmap.writeValue(fh);
	fn_in.writeValue(fh);
	pval.writeValue(fh);
	minres.writeValue(fh);
	maxres.writeValue(fh);
	stepres.writeValue(fh);
	fn_mask.writeValue(fh);
	closeWriteFile(fh);
}

void ResmapJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;
	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_resmap.readValue(fh);
		fn_in.readValue(fh);
		pval.readValue(fh);
		minres.readValue(fh);
		maxres.readValue(fh);
		stepres.readValue(fh);
		fn_mask.readValue(fh);
		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void ResmapJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
	// never submit this to queue, as ResMap needs user interaction
	do_queue.deactivate(true);
}

void ResmapJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
		RFLOAT angpix)
{

	// Change outputname if it contains DATENTIMENRUN
	changeDateNTimeInOutputname(outputname);

	commands.clear();
	std::string command;
	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	if (fn_resmap.getValue().length() == 0)
	{
		std::cerr << "ResmapJobWindow::getCommands ERROR: please provide an executable for the ResMap program." << std::endl;
		exit(1);
	}

	command = fn_resmap.getValue();

	// Get the two half-reconstruction names from the single one
	std::string fn_half1, fn_half2;
	int pos_half = fn_in.getValue().rfind("_half");
	if (pos_half < fn_in.getValue().size())
	{
		fn_half1 = fn_in.getValue().substr(0, pos_half) + "_half1_class001_unfil.mrc";
		fn_half2 = fn_in.getValue().substr(0, pos_half) + "_half2_class001_unfil.mrc";
	}
	else
	{
		std::cerr << "ResMapJobWindow::getCommands ERROR: cannot find _half substring in input filename: " << fn_in.getValue() << std::endl;
		exit(1);
	}
	Node node(fn_in.getValue(), NODE_HALFMAP);
	pipelineInputNodes.push_back(node);

	FileName fn_out = fn_in.getValue();
	fn_out = fn_out.insertBeforeExtension("_resmap");
	Node node1(fn_out, NODE_RESMAP);
	pipelineOutputNodes.push_back(node1);
	pipelineOutputName = fn_in.getValue().substr(0, pos_half) + "_resmap";

	command += " --vis2D --noguiSplit " + fn_half1 + " " + fn_half2;
	command += " --vxSize=" + floatToString(angpix);
	command += " --pVal=" + floatToString(pval.getValue());
	command += " --minRes=" + floatToString(minres.getValue());
	command += " --maxRes=" + floatToString(maxres.getValue());
	command += " --stepRes=" + floatToString(stepres.getValue());
	if (fn_mask.getValue().length() > 0)
	{
		Node node2(fn_mask.getValue(), NODE_3DMASK);
		pipelineInputNodes.push_back(node2);
		command += " --maskVol=" + fn_mask.getValue();
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);

	int last_slash = fn_in.getValue().rfind("/");
	if (last_slash < fn_in.getValue().size())
	{
		// The output name contains a directory: use that one for output
		outputname = fn_in.getValue().substr(0, last_slash + 1) + "run_resmap";
	}
	else
	{
		outputname = "run_resmap";
	}

	prepareFinalCommand(outputname, commands, final_command);
}


PublishJobWindow::PublishJobWindow() : RelionJobWindow(2, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_PUBLISH;

	tab1->begin();
	tab1->label("cite RELION");
	resetHeight();

	cite_text.place(current_y, "\
If RELION is useful in your work, please cite us. Relevant papers are:\n \n \
 * General Bayesian approach (and first mention of RELION): \n \
     Scheres (2012) J. Mol. Biol. (PMID: 22100448)	 \n \n\
 * RELION implementation details and the 3D auto-refine procedure: \n \
     Scheres (2012) J. Struct. Biol. (PMID: 23000701)	 \n \n\
 * Gold-standard FSC and the relevance of the 0.143 criterion: \n \
     Scheres & Chen (2012) Nat. Meth. (PMID: 22842542)	 \n \n\
 * Movie-processing procedure: \n \
     Bai et al. (2013) eLife (PMID: 23427024 )	 \n \n\
 * Correction of mask effects on the FSC curve by randomised phases: \n \
     Chen et al. (2013) Ultramicroscopy (PMID: 23872039)	 \n \n\
 * Particle-polishing: \n \
     Scheres (2014) eLife (PMID: 25122622)	 \n \n\
 * Auto-picking : \n \
     Scheres (2014) J. Struct. Biol. (PMID: 25486611) \n \n \
 * Sub-tomogram averaging : \n \
     Bharat et al. (2015) Structure (submitted). \n \n "
, GUIWIDTH - WCOL0 - 50, GUIHEIGHT_OLD - 150);

	//cite_text.mydisp->textsize(12);

	tab1->end();
	tab2->begin();
	tab2->label("cite others");
	resetHeight();

	cite_external_text.place(current_y, "\
Please also cite the following EXTERNAL programs: \n \n \
* CTFFIND for CTF-estimation: \n \
    Mindell & Grigorieff (2003) J. Mol. Biol. (PMID: 12781660) \n \n\
* ResMap for local-resolution estimation:  \n\
    Kucukelbir et al. (2014) Nat. Meth. (PMID: 24213166)"
, GUIWIDTH - WCOL0 - 50, GUIHEIGHT_OLD - 150);

	tab2->end();

}

void PublishJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
	// never submit this to queue, as ResMap needs user interaction
	do_queue.deactivate(true);
}

void PublishJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command)
{
	commands.clear();
	std::string command = " echo 'Sorry, you still need to write your own paper... ;-)' ";

	commands.push_back(command);
	outputname = "run_publish";
	prepareFinalCommand(outputname, commands, final_command);
}

