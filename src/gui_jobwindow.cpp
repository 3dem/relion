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
std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies, bool do_movies, bool do_also_rot)
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
	int node_type = (do_movies) ? NODE_MOVIE_DATA : NODE_PART_DATA;
	Node node1(fn_out + "_data.star", node_type);
	result.push_back(node1);

	if (!do_movies || do_also_rot)
	{
		if (iter > 0)
		{
			Node node2(fn_out + "_model.star", NODE_MODEL);
			result.push_back(node2);
		}

		// For 3D classification or 3D auto-refine, also use individual 3D maps as outputNodes
		if (dim == 3)
		{
			FileName fn_tmp;
			for (int iclass = 0; iclass < K; iclass++)
			{
				if (nr_bodies > 1)
					fn_tmp.compose(fn_out+"_body", iclass+1, "mrc", 3);
				else
					fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);

				Node node3(fn_tmp, NODE_3DREF);
				result.push_back(node3);
			}
		}

		// For auto-refine: also output the half1_class001_unfil.mrc map
		if (iter < 0)
		{
			Node node4(fn_out+"_half1_class001_unfil.mrc", NODE_HALFMAP);
			result.push_back(node4);
		}
	}
	return result;

}

RelionJobWindow::RelionJobWindow(int nr_tabs, bool _has_mpi, bool _has_thread,
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
			tab7 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab7->end();
			tab7->color(GUI_BACKGROUND_COLOR);
			tab7->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 8)
		{
			std::cerr << "ERROR: only 7 job-specific tabs implemented..." << std::endl;
			exit(1);
		}
		current_y += 15;
	    start_y = current_y;

		runtab = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
		runtab->label("Running");
		setupRunTab();
		runtab->end();
		runtab->color(GUI_BACKGROUND_COLOR);
		runtab->selection_color(GUI_BACKGROUND_COLOR2);

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
		nr_threads.place(current_y, "Number of threads:", 1, 1, 16, 1, "Number of shared-memory (POSIX) threads to use in parallel. \
When set to 1, no multi-threading will be used. Multi-threading is often useful in 3D refinements to have more memory. 2D class averaging often proceeds more efficiently without threads.");

	// Add a little spacer
	if (has_mpi || has_thread)
		current_y += STEPY/4;

    // Set up queue groups for running tab
    queue_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    queue_group->end();

    do_queue.place(current_y, "Submit to queue?", false, "If set to Yes, the job will be submit to a queue, otherwise \
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

	qsubscript.place(current_y, "Standard submission script:", default_location, "Script Files (*.{csh,sh,bash,script})", NULL,
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
    current_y += STEPY/4;

    other_args.place(current_y, "Additional arguments:", "", "In this box command-line arguments may be provided that are not generated by the GUI. \
This may be useful for testing developmental options and/or expert use of the program. \
The command 'relion_refine' will print a list of possible options.");

}

void RelionJobWindow::toggle_new_continue(bool is_continue)
{
	return;
}

void RelionJobWindow::openWriteFile(std::string fn, std::ofstream &fh)
{

	fh.open((fn+"run.job").c_str(), std::ios::out);
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

	fh.open((fn+"run.job").c_str(), std::ios_base::in);
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
		if (!(type >= 0 && type < NR_BROWSE_TABS))
			REPORT_ERROR("RelionJobWindow::openReadFile ERROR: cannot find job type in " + fn + "run.job");
    	// Get is_continue from second line
		getline(fh, line, '\n');
		if (line.rfind("is_continue == true") == 0)
			is_continue = true;
		else
			is_continue = false;

		return true;
    }
}

void RelionJobWindow::closeWriteFile(std::ofstream& fh, std::string fn)
{
	if (has_mpi)
		nr_mpi.writeValue(fh);
	if (has_thread)
		nr_threads.writeValue(fh);
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
	std::string command = "chmod 664 " + fn + "run.job";
	int res = system(command.c_str());

}

void RelionJobWindow::closeReadFile(std::ifstream& fh)
{
	if (has_mpi)
		nr_mpi.readValue(fh);
	if (has_thread)
		nr_threads.readValue(fh);
	do_queue.readValue(fh);
	queuename.readValue(fh);
	qsub.readValue(fh);
	if (have_extra1)
		qsub_extra1.readValue(fh);
	if (have_extra2)
		qsub_extra2.readValue(fh);
	qsubscript.readValue(fh);
	other_args.readValue(fh);

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
	replaceStringAll(textbuf, "XXXerrfileXXX", outputname + "run.err");
	replaceStringAll(textbuf, "XXXoutfileXXX", outputname + "run.out");
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

void RelionJobWindow::initialisePipeline(std::string &outputname, std::string defaultname, int job_counter)
{

	pipelineOutputNodes.clear();
	pipelineInputNodes.clear();

	if (outputname == "") // for continue jobs, use the same outputname
	{
		if (job_counter < 1000)
			outputname = defaultname + "/job" + integerToString(job_counter, 3) + "/";
		else
			outputname = defaultname + "/job" + integerToString(job_counter) + "/";
	}

	pipelineOutputName = outputname;

}

bool RelionJobWindow::prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir)
{

	// Create output directory if the outname contains a "/"
	if (do_makedir)
	{
		int last_slash = outputname.rfind("/");
		if (last_slash < outputname.size())
		{
			std::string dirs = outputname.substr(0, last_slash);
			std::string makedirs = "mkdir -p " + dirs;
			int res = system(makedirs.c_str());
		}
	}

	// Prepare full mpi commands or save jobsubmission script to disc
	if (do_queue.getValue() && do_makedir)
	{
		// Make the submission script and write it to disc
		std::string output_script = outputname + "run_submit.script";
		saveJobSubmissionScript(output_script, outputname, commands);
		final_command = qsub.getValue() + " " + output_script + " &";
	}
	else
	{
		// If there are multiple commands, then join them all on a single line (final_command)
		// Also add mpirun in front of those commands that have _mpi` in it (if no submission via the queue is done)
		std::string one_command;
		final_command = "";
		for (size_t icom = 0; icom < commands.size(); icom++)
		{
			if (has_mpi && nr_mpi.getValue() > 1 && (commands[icom]).find("_mpi`") != std::string::npos)
				one_command = "mpirun -n " + floatToString(nr_mpi.getValue()) + " " + commands[icom] ;
			else
				one_command = commands[icom];
			// Save stdout and stderr to a .out and .err files
			// But only when a re-direct '>' is NOT already present on the command line!
			if (std::string::npos == commands[icom].find(">"))
				one_command += " >> " + outputname + "run.out 2>> " + outputname + "run.err";
			final_command += one_command;
			if (icom == commands.size() - 1)
				final_command += " & "; // end by putting composite job in the background
			else
				final_command += " && "; // execute one command after the other...
		}
	}

	char * my_warn = getenv ("RELION_WARNING_LOCAL_MPI");
	int my_nr_warn = (my_warn == NULL) ? DEFAULTWARNINGLOCALMPI : textToInteger(my_warn);
	if (has_mpi && nr_mpi.getValue() > my_nr_warn && !do_queue.getValue())
	{
		std::string ask;
		ask = "You're submitting a local job with " + integerToString(my_nr_warn) + " parallel MPI processes. Do you really want to run this?\n";
		return fl_choice(ask.c_str(), "Don't run", "Run", NULL);
	}
	else
		return true;
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
	closeWriteFile(fh, fn);
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
bool XXXXJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command)
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
	return prepareFinalCommand(outputname, commands, final_command);
}
*/

ImportJobWindow::ImportJobWindow() : RelionJobWindow(1, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_IMPORT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input files:", "Micrographs/*.mrcs", "Input file (*.*)", NULL, "Select any file(s), possibly using Linux wildcards for 2D micrographs or 3D tomograms that you want to import into a STAR file, which will also be saved as a data Node. \n \n \
Note that for importing coordinate files, one has to give a Linux wildcard, where the *-symbol is before the coordinate-file suffix, e.g. if the micrographs are called mic1.mrc and the coordinate files mic1.box or mic1_autopick.star, one HAS to give '*.box' or '*_autopick.star', respectively.\n \n \
Also note that micrographs, movies and coordinate files all need to be in the same directory (with the same rootnames, e.g.mic1 in the example above) in order to be imported correctly. 3D masks or references can be imported from anywhere. \n \n \
Note that movie-particle STAR files cannot be imported from a previous version of RELION, as the way movies are handled has changed in RELION-2.0. \n \n \
For the import of a particle, 2D references or micrograph STAR file or of a 3D reference or mask, only a single file can be imported at a time.");

	// Add a little spacer
	current_y += STEPY/2;

	node_type.place(current_y, "Node type:", node_type_options, &node_type_options[0], "Select the type of Node this is.");

	tab1->end();

	// read settings if hidden file exists
	read(".gui_import", is_continue);

}

void ImportJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_import";

	std::ofstream fh;
	openWriteFile(fn, fh);

	fn_in.writeValue(fh);
	node_type.writeValue(fh);

	closeWriteFile(fh, fn);

}

void ImportJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_import";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		node_type.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void ImportJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;
	do_queue.deactivate(true);
}

bool ImportJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Import", job_counter);

	std::string command;
	FileName outputstar;

	if (node_type.getValue() == "2D micrograph movies (*.mrcs)")
	{
		outputstar = outputname+"movies.star";
		command = "relion_star_loopheader rlnMicrographMovieName > " + outputstar;;
		commands.push_back(command);
		command = "ls " + fn_in.getValue() + " >> " + outputstar;
		commands.push_back(command);
		Node node(outputstar, NODE_MOVIES);
		pipelineOutputNodes.push_back(node);
	}
	else if (node_type.getValue() == "2D micrographs/tomograms (*.mrc)")
	{
		outputstar = outputname+"micrographs.star";
		command = "relion_star_loopheader rlnMicrographName > " + outputstar;;
		commands.push_back(command);
		command = "ls " + fn_in.getValue() + " >> " + outputstar;
		commands.push_back(command);
		Node node(outputstar, NODE_MICS);
		pipelineOutputNodes.push_back(node);
	}
	else if (node_type.getValue() == "2D/3D particle coordinates (*.box, *_pick.star)")
	{
		// Make the same directory structure of the coordinates
		// Copy all coordinate files into the same subdirectory in the Import directory
		command = "cp --parents " + fn_in.getValue() + " " + outputname;
		commands.push_back(command);
		// Get the coordinate-file suffix separately
		FileName fn_suffix = fn_in.getValue();
		FileName fn_suffix2 = fn_suffix.beforeLastOf("*");
		fn_suffix = fn_suffix.afterLastOf("*");
		fn_suffix = "coords_suffix" + fn_suffix;
		Node node(outputname + fn_suffix, NODE_MIC_COORDS);
		pipelineOutputNodes.push_back(node);
		// Make a suffix file, which contains the actual suffix as a suffix
		command = " echo \\\"" + fn_suffix2 + "*.mrc\\\" > " + outputname + fn_suffix;
		commands.push_back(command);
	}
	else if (node_type.getValue() == "Particles STAR file (.star)" ||
			 node_type.getValue() == "Movie-particles STAR file (.star)" ||
			 node_type.getValue() == "Micrographs STAR file (.star)" ||
			 node_type.getValue() == "2D references (.star or .mrcs)" ||
			 node_type.getValue() == "3D reference (.mrc)" ||
			 node_type.getValue() == "3D mask (.mrc)" ||
			 node_type.getValue() == "Unfiltered half-map (unfil.mrc)")
	{
		FileName fnt = "/" + fn_in.getValue();
		fnt = fnt.afterLastOf("/");
		command = "cp " + fn_in.getValue() + " " + outputname + fnt;
		commands.push_back(command);

		int mynodetype;
		if (node_type.getValue() == "Particles STAR file (.star)")
			mynodetype = NODE_PART_DATA;
		else if (node_type.getValue() == "Movie-particles STAR file (.star)")
			mynodetype = NODE_MOVIE_DATA;
		else if (node_type.getValue() == "Micrographs STAR file (.star)")
			mynodetype = NODE_MICS;
		else if (node_type.getValue() == "2D references (.star or .mrcs)")
			mynodetype = NODE_2DREFS;
		else if (node_type.getValue() == "3D reference (.mrc)")
			mynodetype = NODE_3DREF;
		else if (node_type.getValue() == "3D mask (.mrc)")
			mynodetype = NODE_MASK;
		else if (node_type.getValue() == "Unfiltered half-map (unfil.mrc)")
			mynodetype = NODE_HALFMAP;

		Node node(outputname + fnt, mynodetype);
		pipelineOutputNodes.push_back(node);

		// Also get the other half-map
		if (mynodetype == NODE_HALFMAP)
		{
			FileName fn_inb = fn_in.getValue();
			size_t pos = fn_inb.find("half1");
			if (pos != std::string::npos)
			{
				fn_inb.replace(pos, 5, "half2");

			}
			else
			{
				pos = fn_inb.find("half2");
				if (pos != std::string::npos)
				{
					fn_inb.replace(pos, 5, "half1");
				}
			}
			fnt = "/" + fn_inb;
			fnt = fnt.afterLastOf("/");
			command = "cp " + fn_inb + " " + outputname + fnt;
			commands.push_back(command);

			Node node2(outputname + fnt, mynodetype);
			pipelineOutputNodes.push_back(node2);
		}

	}
	else
	{
		std::cerr << " node_type.getValue()= " << node_type.getValue() << std::endl;
		REPORT_ERROR("ImportJobWindow::getCommands ERROR: Unrecognized menu option.");
	}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

MotioncorrJobWindow::MotioncorrJobWindow() : RelionJobWindow(2, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_MOTIONCORR;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	input_star_mics.place(current_y, "Input movies STAR file:", NODE_MOVIES, "", "STAR files (*.star)", "A STAR file with all micrographs to run MOTIONCORR on");
	tab1->end();

	tab2->begin();
	tab2->label("Motioncorr");
	resetHeight();

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_MOTIONCORR_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTMOTIONCORRLOCATION;
		default_location=mydefault;
	}

	fn_motioncorr_exe.place(current_y, "MOTIONCORR executable:", default_location, "*.*", NULL, "Location of the MOTIONCORR executable. You can control the default of this field by setting environment variable RELION_MOTIONCORR_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	// Add a little spacer
	current_y += STEPY/2;

	bin_factor.place(current_y, "Binning factor (1 or 2):", 1, 1, 2, 1, "Bin the micrographs this much by a windowing operation in the Fourier Tranform. Binning at this level is hard to un-do later on, but may be useful to down-scale super-resolution images.");
	first_frame_ali.place(current_y, "First frame for alignment:", 1, 1, 32, 1, "First frame to use in alignment and corrected average (starts counting at 1). This will be used for MOTIONCORRs -nst and -nss");
	last_frame_ali.place(current_y, "Last frame for alignment:", 0, 0, 32, 1, "Last frame to use in alignment and corrected average (0 means use all). This will be used for MOTIONCORRs -ned and -nes");
	first_frame_sum.place(current_y, "First frame for corrected sum:", 1, 1, 32, 1, "First frame to use in corrected average (starts counting at 1). This will be used for MOTIONCORRs -nst and -nss");
	last_frame_sum.place(current_y, "Last frame for corrected sum:", 0, 0, 32, 1, "Last frame to use in corrected average (0 means use all). This will be used for MOTIONCORRs -ned and -nes");
	do_save_movies.place(current_y, "Save aligned movie stacks?", true,"Save the aligned movie stacks? Say Yes if you want to perform movie-processing in RELION as well. Say No if you only want to correct motions in MOTIONCOR");
	bfactor.place(current_y, "Bfactor:", 150, 0, 1500, 50, "The B-factor (in pixel^2) that MOTIONCORR will apply to the micrographs. The MOTIONCORR Readme.txt says: A bfactor 150 or 200pix^2 is good for most cryoEM image with 2x binned super-resolution image. For unbined image, a larger bfactor is needed.");

	// Add a little spacer
	current_y += STEPY/2;
	gpu_ids.place(current_y, "Which GPUs to use: ", "0", "Provide a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example, to place two ranks on device 0 and one rank on device 1, provide '0:0:1'");

	// Add a little spacer
	current_y += STEPY/2;
	other_motioncorr_args.place(current_y, "Other MOTIONCORR arguments", "", "Additional arguments that need to be passed to MOTIONCORR.");

	tab2->end();

	// read settings if hidden file exists
	read(".gui_motioncorr", is_continue);
}



void MotioncorrJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_motioncorr";

	std::ofstream fh;
	openWriteFile(fn, fh);

	input_star_mics.writeValue(fh);
	fn_motioncorr_exe.writeValue(fh);
	bin_factor.writeValue(fh);
	first_frame_ali.writeValue(fh);
	last_frame_ali.writeValue(fh);
	first_frame_sum.writeValue(fh);
	last_frame_sum.writeValue(fh);
	do_save_movies.writeValue(fh);
	bfactor.writeValue(fh);
	gpu_ids.writeValue(fh);
	other_motioncorr_args.writeValue(fh);

	closeWriteFile(fh, fn);
}

void MotioncorrJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_motioncorr";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		input_star_mics.readValue(fh);
		fn_motioncorr_exe.readValue(fh);
		bin_factor.readValue(fh);
		first_frame_ali.readValue(fh);
		last_frame_ali.readValue(fh);
		first_frame_sum.readValue(fh);
		last_frame_sum.readValue(fh);
		do_save_movies.readValue(fh);
		bfactor.readValue(fh);
		gpu_ids.readValue(fh);
		other_motioncorr_args.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}


void MotioncorrJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	input_star_mics.deactivate(is_continue);
	bin_factor.deactivate(is_continue);
	fn_motioncorr_exe.deactivate(is_continue);
	first_frame_ali.deactivate(is_continue);
	last_frame_ali.deactivate(is_continue);
	first_frame_sum.deactivate(is_continue);
	last_frame_sum.deactivate(is_continue);
	do_save_movies.deactivate(is_continue);
	bfactor.deactivate(is_continue);
	other_motioncorr_args.deactivate(is_continue);

}

bool MotioncorrJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "MotionCorr", job_counter);

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_run_motioncorr_mpi`";
	else
		command="`which relion_run_motioncorr`";

	// I/O

	command += " --i " + input_star_mics.getValue();
	Node node(input_star_mics.getValue(), input_star_mics.type);
	pipelineInputNodes.push_back(node);

	command += " --o " + outputname;
	pipelineOutputName = outputname;
	Node node2(outputname + "corrected_micrographs.star", NODE_MICS);
	pipelineOutputNodes.push_back(node2);
	Node node3(outputname + "corrected_micrograph_movies.star", NODE_MOVIES);
	pipelineOutputNodes.push_back(node3);

	// Motioncorr-specific stuff
	command += " --bin_factor " + floatToString(bin_factor.getValue());
	command += " --motioncorr_exe " + fn_motioncorr_exe.getValue();
	command += " --first_frame_ali " + floatToString(first_frame_ali.getValue());
	command += " --last_frame_ali " + floatToString(last_frame_ali.getValue());
	command += " --first_frame_sum " + floatToString(first_frame_sum.getValue());
	command += " --last_frame_sum " + floatToString(last_frame_sum.getValue());
	command += " --bfactor " + floatToString(bfactor.getValue());

	if (do_save_movies.getValue())
		command += " --save_movies ";

	if ((other_motioncorr_args.getValue()).length() > 0)
		command += " --other_motioncorr_args \"" + other_motioncorr_args.getValue() + "\"";

	if (is_continue)
		command += " --only_do_unfinished ";

	// Which GPUs to use?
	command += " --gpu " + gpu_ids.getValue();

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}


CtffindJobWindow::CtffindJobWindow() : RelionJobWindow(5, HAS_MPI, HAS_NOT_THREAD)
{
	type = PROC_CTFFIND;

	char *default_location;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();


	input_star_mics.place(current_y, "Input micrographs STAR file:", NODE_MICS, "", "STAR files (*.star)", "A STAR file with all micrographs to run CTFFIND or Gctf on");

	tab1->end();


	tab2->begin();
	tab2->label("Microscopy");
	resetHeight();

	cs.place(current_y, "Spherical aberration (mm):", 2, 0, 8, 0.1, "Spherical aberration of the microscope used to collect these images (in mm)");

	kv.place(current_y, "Voltage (kV):", 300, 50, 500, 10, "Voltage the microscope was operated on (in kV)");

	q0.place(current_y, "Amplitude contrast:", 0.1, 0, 0.3, 0.01, "Fraction of amplitude contrast. Often values around 10% work better than theoretically more accurate lower values...");

	angpix.place(current_y, "Magnified pixel size (Angstrom):", 1.4, 0.5, 3, 0.1, "Pixel size in Angstroms. ");

	tab2->end();

	tab3->begin();
	tab3->label("CTFFIND");
	resetHeight();

	// Check for environment variable RELION_CTFFIND_EXECUTABLE
	default_location = getenv ("RELION_CTFFIND_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTCTFFINDLOCATION;
		default_location=mydefault;
	}
	fn_ctffind_exe.place(current_y, "CTFFIND executable:", default_location, "*.exe", NULL, "Location of the CTFFIND executable. You can control the default of this field by setting environment variable RELION_CTFFIND_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	// Add a little spacer
	current_y += STEPY/2;

	box.place(current_y, "FFT box size (pix):", 512, 64, 1024, 8, "CTFFIND's Box parameter");

	resmin.place(current_y, "Minimum resolution (A):", 30, 10, 200, 10, "CTFFIND's ResMin parameter");

	resmax.place(current_y, "Maximum resolution (A):", 5, 1, 20, 1, "CTFFIND's ResMax parameter");

	dfmin.place(current_y, "Minimum defocus value (A):", 5000, 0, 25000, 1000, "CTFFIND's dFMin parameter");

	dfmax.place(current_y, "Maximum defocus value (A):", 50000, 20000, 100000, 1000, "CTFFIND's dFMax parameter");

	dfstep.place(current_y, "Defocus step size (A):", 500, 200, 2000, 100,"CTFFIND's FStep parameter");

	dast.place(current_y, "Amount of astigmatism (A):", 100, 0, 2000, 100,"CTFFIND's dAst parameter");

	// Add a little spacer
	current_y += STEPY/2;

	ctf_win.place(current_y, "Estimate CTF on window size (pix) ", -1, -16, 4096, 16, "If a positive value is given, a squared window of this size at the center of the micrograph will be used to estimate the CTF. This may be useful to exclude parts of the micrograph that are unsuitable for CTF estimation, e.g. the labels at the edge of phtographic film. \n \n The original micrograph will be used (i.e. this option will be ignored) if a negative value is given.");

	tab3->end();


	tab4->begin();
	tab4->label("CTFFIND4");
	resetHeight();

	ctffind4_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctffind4_group->end();

	is_ctffind4.place(current_y, "Is this a CTFFIND 4.1+ executable?", false, "If set to Yes, the wrapper will use the extended functionaility of CTFFIND4 (version 4.1 or newer). This includes thread-support, calculation of Thon rings from movie frames and phase-shift estimation for phase-plate data.", ctffind4_group);
	ctffind4_group->begin();

	movie_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	movie_group->end();

	do_movie_thon_rings.place(current_y, "Estimate Thon rongs from movies?", false, "If set to Yes, CTFFIND4 will calculate power spectra of averages of several movie frames and then average those power spectra to calculate Thon rings. This may give better rings than calculation the power spectra of averages of all movie frames, although it does come at increased costs of processing and disk access", movie_group);

	movie_group->begin();
	movie_rootname.place(current_y, "Movie rootname plus extension", "_movie.mrcs", "Give the movie rootname and extension for all movie files. Movies are assumed to be next to the average micrographs in the same directory.");
	avg_movie_frames.place(current_y, "Nr of movie frames to average:", 4, 1, 20, 1,"Calculate averages over so many movie frames to calculate power spectra. Often values corresponding to an accumulated dose of ~ 4 electrons per squared Angstrom work well.");
	movie_group->end();
	do_movie_thon_rings.cb_menu_i(); // make default active

	phaseshift_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	phaseshift_group->end();

	do_phaseshift.place(current_y, "Estimate phase shifts?", false, "If set to Yes, CTFFIND4 will estimate the phase shift, e.g. as introduced by a Volta phase-plate", phaseshift_group);

	phaseshift_group->begin();
	phase_min.placeOnSameYPosition(current_y, "Phase shift - Min, Max, Step (deg):", "Phase shift search (deg) - Min:", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	phase_max.placeOnSameYPosition(current_y, "", "Phase shift search (deg) - Max:", "180", NULL, XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	phase_step.placeOnSameYPosition(current_y, "", "Phase shift search (deg) - Step:", "10", "Minimum, maximum and step size (in degrees) for the search of the phase shift", XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	phaseshift_group->end();
	do_phaseshift.cb_menu_i(); // make default active

	ctffind4_group->end();
	is_ctffind4.cb_menu_i(); // make default active

	tab4->end();
	tab5->begin();
	tab5->label("Gctf");
	resetHeight();

	gctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	gctf_group->end();

	use_gctf.place(current_y, "Use Gctf instead of CTFFIND?", false, "If set to Yes, Kai Zhang's Gctf program (which runs on NVIDIA GPUs) will be used instead of Niko Grigorieff's CTFFIND.", gctf_group);

	gctf_group->begin();
	// Check for environment variable RELION_CTFFIND_EXECUTABLE
	default_location = getenv ("RELION_GCTF_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTGCTFLOCATION;
		default_location=mydefault;
	}
	fn_gctf_exe.place(current_y, "Gctf executable:", default_location, "*", NULL, "Location of the Gctf executable. You can control the default of this field by setting environment variable RELION_GCTF_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	do_ignore_ctffind_params.place(current_y, "Ignore CTFFIND parameters?", true, "If set to Yes, all parameters on the CTFFIND tab will be ignored, and Gctf's default parameters will be used (box.size=1024; min.resol=50; max.resol=4; min.defocus=500; max.defocus=90000; step.defocus=500; astigm=1000) \n \
\n If set to No, all parameters on the CTFFIND tab will be passed to Gctf.");

	do_EPA.place(current_y, "Perform equi-phase averaging?", false, "If set to Yes, equi-phase averaging is used in the defocus refinement, otherwise basic rotational averaging will be performed.");

	// Add a little spacer
	current_y += STEPY/2;

	gpu_ids.place(current_y, "Which GPUs to use: ", "", "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. ");

	//	other_gctf_args.place(current_y, "Perform equi-phase averaging?", true, "If set to Yes, equi-phase averaging is used in the defocus refinement, otherwise basic rotational averaging will be performed.");

	gctf_group->end();
	use_gctf.cb_menu_i(); // make default active

	tab5->end();

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

	input_star_mics.writeValue(fh);

	cs.writeValue(fh);
	kv.writeValue(fh);
	q0.writeValue(fh);
	angpix.writeValue(fh);

	box.writeValue(fh);
	resmin.writeValue(fh);
	resmax.writeValue(fh);
	dfmin.writeValue(fh);
	dfmax.writeValue(fh);
	dfstep.writeValue(fh);
	dast.writeValue(fh);
	fn_ctffind_exe.writeValue(fh);
	ctf_win.writeValue(fh);

	is_ctffind4.writeValue(fh);
	do_movie_thon_rings.writeValue(fh);
	movie_rootname.writeValue(fh);
	avg_movie_frames.writeValue(fh);
	do_phaseshift.writeValue(fh);
	phase_min.writeValue(fh);
	phase_max.writeValue(fh);
	phase_step.writeValue(fh);

	use_gctf.writeValue(fh);
	fn_gctf_exe.writeValue(fh);
	do_ignore_ctffind_params.writeValue(fh);
	do_EPA.writeValue(fh);
	gpu_ids.writeValue(fh);

	closeWriteFile(fh, fn);
}

void CtffindJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_ctffind";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		input_star_mics.readValue(fh);

		cs.readValue(fh);
		kv.readValue(fh);
		q0.readValue(fh);
		angpix.readValue(fh);

		box.readValue(fh);
		resmin.readValue(fh);
		resmax.readValue(fh);
		dfmin.readValue(fh);
		dfmax.readValue(fh);
		dfstep.readValue(fh);
		dast.readValue(fh);
		fn_ctffind_exe.readValue(fh);

		is_ctffind4.readValue(fh);
		do_movie_thon_rings.readValue(fh);
		movie_rootname.readValue(fh);
		avg_movie_frames.readValue(fh);
		do_phaseshift.readValue(fh);
		phase_min.readValue(fh);
		phase_max.readValue(fh);
		phase_step.readValue(fh);

		ctf_win.readValue(fh);
		use_gctf.readValue(fh);
		fn_gctf_exe.readValue(fh);
		do_ignore_ctffind_params.readValue(fh);
		do_EPA.readValue(fh);
		gpu_ids.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void CtffindJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	input_star_mics.deactivate(is_continue);

	cs.deactivate(is_continue);
	kv.deactivate(is_continue);
	q0.deactivate(is_continue);
	angpix.deactivate(is_continue);

	box.deactivate(is_continue);
	resmin.deactivate(is_continue);
	resmax.deactivate(is_continue);
	dfmin.deactivate(is_continue);
	dfmax.deactivate(is_continue);
	dfstep.deactivate(is_continue);
	dast.deactivate(is_continue);
	fn_ctffind_exe.deactivate(is_continue);
	ctf_win.deactivate(is_continue);

	is_ctffind4.deactivate(is_continue);
	do_movie_thon_rings.deactivate(is_continue);
	movie_rootname.deactivate(is_continue);
	avg_movie_frames.deactivate(is_continue);
	do_phaseshift.deactivate(is_continue);
	phase_min.deactivate(is_continue);
	phase_max.deactivate(is_continue);
	phase_step.deactivate(is_continue);
	use_gctf.deactivate(is_continue);
	fn_gctf_exe.deactivate(is_continue);
	do_ignore_ctffind_params.deactivate(is_continue);
	do_EPA.deactivate(is_continue);

}

bool CtffindJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "CtfFind", job_counter);
	std::string command;

	FileName fn_outstar = outputname + "micrographs_ctf.star";
	Node node(fn_outstar, NODE_MICS);
	pipelineOutputNodes.push_back(node);

	Node node2(input_star_mics.getValue(), input_star_mics.type);
	pipelineInputNodes.push_back(node2);

	if (nr_mpi.getValue() > 1)
		command="`which relion_run_ctffind_mpi`";
	else
		command="`which relion_run_ctffind`";

	// Just use 10000 for magnification, so dstep==angpix
	RFLOAT magn = 10000.;

	command += " --i " + input_star_mics.getValue();
	command += " --o " + outputname;
	command += " --CS " + floatToString(cs.getValue());
	command += " --HT " + floatToString(kv.getValue());
	command += " --AmpCnst " + floatToString(q0.getValue());
	command += " --XMAG " + floatToString(magn);
	command += " --DStep " + floatToString(angpix.getValue());
	command += " --Box " + floatToString(box.getValue());
	command += " --ResMin " + floatToString(resmin.getValue());
	command += " --ResMax " + floatToString(resmax.getValue());
	command += " --dFMin " + floatToString(dfmin.getValue());
	command += " --dFMax " + floatToString(dfmax.getValue());
	command += " --FStep " + floatToString(dfstep.getValue());
	command += " --dAst " + floatToString(dast.getValue());
	if (use_gctf.getValue())
	{
		command += " --use_gctf --gctf_exe " + fn_gctf_exe.getValue();
		command += " --angpix " + floatToString(angpix.getValue());
		if (do_ignore_ctffind_params.getValue())
			command += " --ignore_ctffind_params";
		if (do_EPA.getValue())
			command += " --EPA";

		// GPU-allocation
		command += " --gpu " + gpu_ids.getValue();
	}
	else
	{
		command += " --ctffind_exe " + fn_ctffind_exe.getValue();
		command += " --ctfWin " + floatToString(ctf_win.getValue());
		if (is_ctffind4.getValue())
		{
			command += " --is_ctffind4 ";
			if (do_movie_thon_rings.getValue())
			{
				command += " --do_movie_thon_rings --avg_movie_frames " + floatToString(avg_movie_frames.getValue());
				command += " --movie_rootname " + movie_rootname.getValue();
			}
			if (do_phaseshift.getValue())
			{
				command += " --do_phaseshift ";
				command += " --phase_min " + floatToString(textToFloat(phase_min.getValue()));
				command += " --phase_max " + floatToString(textToFloat(phase_max.getValue()));
				command += " --phase_step " + floatToString(textToFloat(phase_step.getValue()));
			}
		}
	}


	if (is_continue)
		command += " --only_do_unfinished ";

	// Other arguments
	command += " " + other_args.getValue();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}



ManualpickJobWindow::ManualpickJobWindow() : RelionJobWindow(3, HAS_NOT_MPI, HAS_NOT_THREAD)
{
	type = PROC_MANUALPICK;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input micrographs:", NODE_MICS, "", "Input micrographs (*.{star,mrc})", "Input STAR file (with or without CTF information), OR a unix-type wildcard with all micrographs in MRC format (in this case no CTFs can be used).");
	tab1->end();

	tab2->begin();
	tab2->label("Display");
	resetHeight();

	diameter.place(current_y, "Particle diameter (A):", 100, 0, 500, 50, "The radius of the circle used around picked particles (in original pixels). Only used for display." );
	micscale.place(current_y, "Scale for micrographs:", 0.2, 0.1, 1, 0.05, "The micrographs will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );
	sigma_contrast.place(current_y, "Sigma contrast:", 3, 0, 10, 0.5, "The micrographs will be displayed with the black value set to the average of all values MINUS this values times the standard deviation of all values in the micrograph, and the white value will be set \
to the average PLUS this value times the standard deviation. Use zero to set the minimum value in the micrograph to black, and the maximum value to white ");
	white_val.place(current_y, "White value:", 0, 0, 512, 16, "Use non-zero values to set the value of the whitest pixel in the micrograph.");
	black_val.place(current_y, "Black value:", 0, 0, 512, 16, "Use non-zero values to set the value of the blackest pixel in the micrograph.");

	current_y += STEPY/2;
	lowpass.place(current_y, "Lowpass filter (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the micrographs. Give a negative value to skip the lowpass filter.");
	highpass.place(current_y, "Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter. Useful values are often in the range of 200-400 Angstroms.");
	angpix.place(current_y, "Pixel size (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. This will be used to calculate the filters and the particle diameter in pixels. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");

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

	fn_color.place(current_y, "STAR file with color label: ", "", "STAR file (*.star)", NULL, "The program will figure out which particles in this STAR file are on the current micrograph and color their circles according to the value in the corresponding column. \
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
	lowpass.writeValue(fh);
	highpass.writeValue(fh);
	angpix.writeValue(fh);
	diameter.writeValue(fh);
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


	closeWriteFile(fh, fn);
}

void ManualpickJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_manualpick";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		lowpass.readValue(fh);
		highpass.readValue(fh);
		angpix.readValue(fh);
		diameter.readValue(fh);
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

	fn_in.deactivate(is_continue);
	do_queue.deactivate(true);
}

bool ManualpickJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "ManualPick", job_counter);

	std::string command;
	command="`which relion_manualpick`";

	command += " --i " + fn_in.getValue();
	Node node(fn_in.getValue(), fn_in.type);
	pipelineInputNodes.push_back(node);

	command += " --odir " + outputname;
	command += " --pickname manualpick";

	FileName fn_suffix = outputname + "coords_suffix_manualpick.star";
	Node node2(fn_suffix, NODE_MIC_COORDS);
	pipelineOutputNodes.push_back(node2);

	command += " --scale " + floatToString(micscale.getValue());
	command += " --sigma_contrast " + floatToString(sigma_contrast.getValue());
	command += " --black " + floatToString(black_val.getValue());
	command += " --white " + floatToString(white_val.getValue());

	if (lowpass.getValue() > 0.)
		command += " --lowpass " + floatToString(lowpass.getValue());
	if (highpass.getValue() > 0.)
		command += " --highpass " + floatToString(highpass.getValue());
	if (angpix.getValue() > 0.)
		command += " --angpix " + floatToString(angpix.getValue());

	command += " --ctf_scale " + floatToString(ctfscale.getValue());

	command += " --particle_diameter " + floatToString(diameter.getValue());

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

	// Also make the suffix file (do this after previous command was pushed back!)
	// Inside it, store the name of the micrograph STAR file, so we can display these later
	command = "echo " + fn_in.getValue() + " > " + fn_suffix;
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);
}


AutopickJobWindow::AutopickJobWindow() : RelionJobWindow(4, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_AUTOPICK;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_input_autopick.place(current_y, "Input micrographs for autopick:", NODE_MICS, "", "Input micrographs (*.{star})", "Input STAR file (preferably with CTF information) with all micrographs to pick from.");
	fn_refs_autopick.place(current_y, "References:", NODE_2DREFS, "", "Input references (*.{star,mrcs})", "Input STAR file or MRC stack with the 2D references to be used for picking. Note that the absolute greyscale needs to be correct, so only use images created by RELION itself, e.g. by 2D class averaging or projecting a RELION reconstruction.");

	// Add a little spacer
	current_y += STEPY/2;

	angpix.place(current_y, "Pixel size in micrographs (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");

	tab1->end();
	tab2->begin();
	tab2->label("References");
	resetHeight();

	//set up group
	autopick_ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	autopick_ctf_group->end();


	lowpass.place(current_y, "Lowpass filter references (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the references before template matching. Do NOT use very high-resolution templates to search your micrographs. The signal will be too weak at high resolution anyway, and you may find Einstein from noise.... Give a negative value to skip the lowpass filter.");
	highpass.place(current_y, "Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter.  Useful values are often in the range of 200-400 Angstroms.");
	angpix_ref.place(current_y, "Pixel size in references (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms for the provided reference images. This will be used to calculate the filters and the particle diameter in pixels. If a negative value is given here, the pixel size in the references will be assumed to be the same as the one in the micrographs, i.e. the particles that were used to make the references were not rescaled upon extraction.");
	particle_diameter.place(current_y, "Mask diameter (A)", -1, 0, 2000, 20, "Diameter of the circular mask that will be applied around the templates in Angstroms. When set to a negative value, this value is estimated automatically from the templates themselves.");

	// Add a little spacer
	current_y += STEPY/2;

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


	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	shrink.place(current_y, "Shrink factor:", 1, 0, 1, 0.1, "This is useful to speed up the calculations, and to make them less memory-intensive. The micrographs will be downscaled (shrunk) to calculate the cross-correlations, and peak searching will be done in the downscaled FOM maps. When set to 0, the micrographs will de downscaled to the lowpass filter of the references, a value between 0 and 1 will downscale the micrographs by that factor. Note that the results will not be exactly the same when you shrink micrographs!");
	gpu_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    gpu_group->end();
	use_gpu.place(current_y, "Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.", gpu_group);
	gpu_group->begin();
	gpu_ids.place(current_y, "Which GPUs to use:", "", "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example: 0:1:0:1:0:1");
    gpu_group->end();
	use_gpu.cb_menu_i(); // This is to make the default effective

	tab3->end();
	tab4->begin();
	tab4->label("Helix");
	resetHeight();

	autopick_helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	autopick_helix_group->end();

	do_pick_helical_segments.place(current_y, "Pick 2D helical segments?", false, "Set to Yes if you want to pick 2D helical segments.", autopick_helix_group);

	autopick_helix_group->begin();

	helical_tube_kappa_max.place(current_y, "Maximum curvature (kappa): ", 0.1, 0.05, 0.5, 0.01, "Maximum curvature allowed for picking helical tubes. \
Kappa = 0.3 means that the curvature of the picked helical tubes should not be larger than 30% the curvature of a circle (diameter = particle mask diameter). \
Kappa ~ 0.05 is recommended for long and straight tubes (e.g. TMV, VipA/VipB and AChR tubes) while 0.20 ~ 0.40 seems suitable for flexible ones (e.g. ParM and MAVS-CARD filaments).");

	helical_tube_outer_diameter.place(current_y, "Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of the tubes.");

	helical_tube_length_min.place(current_y, "Minimum length (A): ", 200, 100, 1000, 10, "Minimum length (in Angstroms) of helical tubes for auto-picking. \
Helical tubes with shorter lengths will not be picked. Note that a long helical tube seen by human eye might be treated as short broken pieces due to low FOM values or high picking threshold.");

	autopick_helix_group->end();

	do_pick_helical_segments.cb_menu_i();

	tab4->end();

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
	do_invert_refs.writeValue(fh);
	do_ctf_autopick.writeValue(fh);
	do_ignore_first_ctfpeak_autopick.writeValue(fh);
	lowpass.writeValue(fh);
	highpass.writeValue(fh);
	angpix.writeValue(fh);
	angpix_ref.writeValue(fh);
	particle_diameter.writeValue(fh);
	psi_sampling_autopick.writeValue(fh);
	do_write_fom_maps.writeValue(fh);
	do_read_fom_maps.writeValue(fh);
	threshold_autopick.writeValue(fh);
	mindist_autopick.writeValue(fh);
	maxstddevnoise_autopick.writeValue(fh);
	do_pick_helical_segments.writeValue(fh);
	helical_tube_kappa_max.writeValue(fh);
	helical_tube_outer_diameter.writeValue(fh);
	helical_tube_length_min.writeValue(fh);
	use_gpu.writeValue(fh);
	gpu_ids.writeValue(fh);
	shrink.writeValue(fh);

	closeWriteFile(fh, fn);
}

void AutopickJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_autopick";

	if (openReadFile(fn, fh))
	{

		fn_input_autopick.readValue(fh);
		fn_refs_autopick.readValue(fh);
		do_invert_refs.readValue(fh);
		do_ctf_autopick.readValue(fh);
		do_ignore_first_ctfpeak_autopick.readValue(fh);
		lowpass.readValue(fh);
		highpass.readValue(fh);
		angpix.readValue(fh);
		angpix_ref.readValue(fh);
		particle_diameter.readValue(fh);
		psi_sampling_autopick.readValue(fh);
		do_write_fom_maps.readValue(fh);
		do_read_fom_maps.readValue(fh);
		threshold_autopick.readValue(fh);
		mindist_autopick.readValue(fh);
		maxstddevnoise_autopick.readValue(fh);
		do_pick_helical_segments.readValue(fh);
		helical_tube_kappa_max.readValue(fh);
		helical_tube_outer_diameter.readValue(fh);
		helical_tube_length_min.readValue(fh);
		use_gpu.readValue(fh);
		gpu_ids.readValue(fh);
		shrink.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void AutopickJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_input_autopick.deactivate(is_continue);
	fn_refs_autopick.deactivate(is_continue);
	do_invert_refs.deactivate(is_continue);
	do_ctf_autopick.deactivate(is_continue);
	do_ignore_first_ctfpeak_autopick.deactivate(is_continue);
	lowpass.deactivate(is_continue);
	highpass.deactivate(is_continue);
	angpix.deactivate(is_continue);
	angpix_ref.deactivate(is_continue);
	particle_diameter.deactivate(is_continue);
	psi_sampling_autopick.deactivate(is_continue);
	shrink.deactivate(is_continue);

}

bool AutopickJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{


	commands.clear();
	initialisePipeline(outputname, "AutoPick", job_counter);

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_autopick_mpi`";
	else
		command="`which relion_autopick`";

	// Input
	command += " --i " + fn_input_autopick.getValue();
	Node node(fn_input_autopick.getValue(), fn_input_autopick.type);
	pipelineInputNodes.push_back(node);
	command += " --ref " + fn_refs_autopick.getValue();
	Node node2(fn_refs_autopick.getValue(), fn_refs_autopick.type);
	pipelineInputNodes.push_back(node2);

	// Output
	Node node3(outputname + "coords_suffix_autopick.star", NODE_MIC_COORDS);
	pipelineOutputNodes.push_back(node3);

	command += " --odir " + outputname;
	command += " --pickname autopick";

	if (do_invert_refs.getValue())
		command += " --invert ";

	if (do_ctf_autopick.getValue())
	{
		command += " --ctf ";
		if (do_ignore_first_ctfpeak_autopick.getValue())
			command += " --ctf_intact_first_peak ";
	}
	command += " --ang " + floatToString(psi_sampling_autopick.getValue());

	command += " --shrink " + floatToString(shrink.getValue());
	if (lowpass.getValue() > 0.)
		command += " --lowpass " + floatToString(lowpass.getValue());
	if (highpass.getValue() > 0.)
		command += " --highpass " + floatToString(highpass.getValue());
	if (angpix.getValue() > 0.)
		command += " --angpix " + floatToString(angpix.getValue());
	if (angpix_ref.getValue() > 0.)
		command += " --angpix_ref " + floatToString(angpix_ref.getValue());
	if (particle_diameter.getValue() > 0.)
		command += " --particle_diameter " + floatToString(particle_diameter.getValue());

	if (do_write_fom_maps.getValue())
		command += " --write_fom_maps ";

	if (do_read_fom_maps.getValue())
		command += " --read_fom_maps ";

	command += " --threshold " + floatToString(threshold_autopick.getValue());
	command += " --min_distance " + floatToString(mindist_autopick.getValue());
	command += " --max_stddev_noise " + floatToString(maxstddevnoise_autopick.getValue());

	// Helix
	if (do_pick_helical_segments.getValue())
	{
		command += " --helix";
		command += " --helical_tube_kappa_max " + floatToString(helical_tube_kappa_max.getValue());
		command += " --helical_tube_outer_diameter " + floatToString(helical_tube_outer_diameter.getValue());
		command += " --helical_tube_length_min " + floatToString(helical_tube_length_min.getValue());
	}

	if (is_continue && !(do_read_fom_maps.getValue() || do_write_fom_maps.getValue()))
		command += " --only_do_unfinished ";

	// GPU-stuff
	if (use_gpu.getValue())
	{
		// for the moment always use --shrink 0 with GPUs ...
		command += " --gpu " + gpu_ids.getValue();
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	// Also touch the suffix file. Do this after the first command had completed
	command = "echo " + fn_input_autopick.getValue() + " > " +  outputname + "coords_suffix_autopick.star";
	commands.push_back(command.c_str());

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}


ExtractJobWindow::ExtractJobWindow() : RelionJobWindow(3, HAS_MPI, HAS_NOT_THREAD)
{
	type = PROC_EXTRACT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

    star_mics.place(current_y,"micrograph STAR file: ", NODE_MICS, "", "Input STAR file (*.{star})", "Filename of the STAR file that contains all micrographs from which to extract particles.");

	current_y += STEPY/2;
	coords_suffix.place(current_y,"Input coordinates: ", NODE_MIC_COORDS, "", "Input coords_suffix file ({coords_suffix}*)", "Filename of the coords_suffix file with the directory structure and the suffix of all coordinate files.");

	reextract_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	reextract_group->end();

	do_reextract.place(current_y, "OR re-extract refined particles? ", false, "If set to Yes, the input Coordinates above will be ignored. Instead, one uses a _data.star file from a previous 2D or 3D refinement to re-extract the particles in that refinement, possibly re-centered with their refined origin offsets. This is particularly useful when going from binned to unbinned particles.", reextract_group);

	reextract_group->begin();

	fndata_reextract.place(current_y,"Refined particles STAR file: ", NODE_PART_DATA, "", "Input STAR file (*.{star})", "Filename of the STAR file with the refined particle coordinates, e.g. from a previous 2D or 3D classification or auto-refine run.");
	do_recenter.place(current_y, "Re-center refined coordinates? ", true, "If set to Yes, the input coordinates will be re-centered according to the refined origin offsets in the provided _data.star file .");

	reextract_group->end();
	do_reextract.cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	set_angpix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	set_angpix_group->end();
	do_set_angpix.place(current_y, "Manually set pixel size? ", false, "If set to Yes, the rlnMagnification and rlnDetectorPixelSize will be set in the resulting STAR file. Only use this option if CTF information is NOT coming from the input coordinate STAR file(s). For example, because you decided not to estimate the CTF for your micrographs.", set_angpix_group);
	set_angpix_group->begin();
	angpix.place(current_y, "Pixel size (A)", 1, 0.3, 5, 0.1, "Provide the pixel size in Angstroms in the micrograph (so before any re-scaling).  If you provide input CTF parameters, then leave this value to the default of -1.");
	set_angpix_group->end();
	do_set_angpix.cb_menu_i();

	tab1->end();

	tab2->begin();
	tab2->label("extract");
	resetHeight();

	extract_size.place(current_y,"Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
	do_invert.place(current_y, "Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted.");

	// Add a little spacer
	current_y += STEPY/2;

	norm_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	norm_group->end();
	do_norm.place(current_y, "Normalize particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it.", norm_group);

	norm_group->begin();


	bg_diameter.place(current_y, "Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");

	white_dust.place(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	black_dust.place(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	norm_group->end();
	do_norm.cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	rescale_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	rescale_group->end();
	do_rescale.place(current_y, "Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.", rescale_group);
	rescale_group->begin();
	rescale.place(current_y, "Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number");
	rescale_group->end();
	do_rescale.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Helix");
	resetHeight();

	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	do_extract_helix.place(current_y, "Extract helical segments?", false, "Set to Yes if you want to extract helical segments. RELION (.star), EMAN2 (.box) and XIMDISP (.coords) formats of tube or segment coordinates are supported.", helix_group);

	helix_group->begin();

	helical_tube_outer_diameter.place(current_y, "Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of helical tubes.");

	current_y += STEPY/2;
	helical_bimodal_angular_priors.place(current_y, "Use bimodal angular priors?", true, "Normally it should be set to Yes and bimodal angular priors will be applied in the following classification and refinement jobs. \
Set to No if the 3D helix looks the same when rotated upside down.");

	helical_tubes_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helical_tubes_group->end();

	current_y += STEPY/2;
	do_extract_helical_tubes.place(current_y, "Coordinates are start-end only?", true, "Set to Yes if you want to extract helical segments from manually picked tube coordinates (starting and end points of helical tubes in RELION, EMAN or XIMDISP format). \
Set to No if segment coordinates (RELION auto-picked results or EMAN / XIMDISP segments) are provided.", helical_tubes_group);

	helical_tubes_group->begin();

	do_cut_into_segments.place(current_y, "Cut helical tubes into segments?", true, "Set to Yes if you want to extract multiple helical segments with a fixed inter-box distance. \
If it is set to No, only one box at the center of each helical tube will be extracted.");

	helical_nr_asu.place(current_y, "Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. This integer should not be less than 1. The inter-box distance (pixels) = helical rise (Angstroms) * number of asymmetrical units / pixel size (Angstroms). \
The optimal inter-box distance might also depend on the box size, the helical rise and the flexibility of the structure. In general, an inter-box distance of ~10% * the box size seems appropriate.");

	helical_rise.place(current_y, "Helical rise (A):", 1, 0, 100, 0.01, "Helical rise in Angstroms. (Please click '?' next to the option above for details about how the inter-box distance is calculated.)");

	helical_tubes_group->end();

	do_extract_helical_tubes.cb_menu_i();

	helix_group->end();

	do_extract_helix.cb_menu_i();

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
	coords_suffix.writeValue(fh);
	do_set_angpix.writeValue(fh);
	angpix.writeValue(fh);
	do_reextract.writeValue(fh);
	fndata_reextract.writeValue(fh);
	do_recenter.writeValue(fh);

	// extract
	extract_size.writeValue(fh);
	do_rescale.writeValue(fh);
	rescale.writeValue(fh);
	do_norm.writeValue(fh);
	bg_diameter.writeValue(fh);
	white_dust.writeValue(fh);
	black_dust.writeValue(fh);
	do_invert.writeValue(fh);

	// Helix
	do_extract_helix.writeValue(fh);
	do_extract_helical_tubes.writeValue(fh);
	do_cut_into_segments.writeValue(fh);
	helical_nr_asu.writeValue(fh);
	helical_rise.writeValue(fh);
	helical_tube_outer_diameter.writeValue(fh);
	helical_bimodal_angular_priors.writeValue(fh);

	closeWriteFile(fh, fn);
}
void ExtractJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_extract";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{

		// I/O
		star_mics.readValue(fh);
		coords_suffix.readValue(fh);
		do_set_angpix.readValue(fh);
		angpix.readValue(fh);
		do_reextract.readValue(fh);
		fndata_reextract.readValue(fh);
		do_recenter.readValue(fh);

		// extract
		extract_size.readValue(fh);
		do_rescale.readValue(fh);
		rescale.readValue(fh);
		do_norm.readValue(fh);
		bg_diameter.readValue(fh);
		white_dust.readValue(fh);
		black_dust.readValue(fh);
		do_invert.readValue(fh);

		// Helix
		do_extract_helix.readValue(fh);
		do_extract_helical_tubes.readValue(fh);
		do_cut_into_segments.readValue(fh);
		helical_nr_asu.readValue(fh);
		helical_rise.readValue(fh);
		helical_tube_outer_diameter.readValue(fh);
		helical_bimodal_angular_priors.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void ExtractJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	// For continuation: deactivate everything exceptthe movie stuff
	star_mics.deactivate(is_continue);
	coords_suffix.deactivate(is_continue);
	do_set_angpix.deactivate(is_continue);
	angpix.deactivate(is_continue);
	do_reextract.deactivate(is_continue);
	fndata_reextract.deactivate(is_continue);
	do_recenter.deactivate(is_continue);

	// extract
	extract_size.deactivate(is_continue);
	do_rescale.deactivate(is_continue);
	rescale.deactivate(is_continue);
	do_norm.deactivate(is_continue);
	bg_diameter.deactivate(is_continue);
	white_dust.deactivate(is_continue);
	black_dust.deactivate(is_continue);
	do_invert.deactivate(is_continue);

	// Helix
	do_extract_helix.deactivate(is_continue);
	do_extract_helical_tubes.deactivate(is_continue);
	do_cut_into_segments.deactivate(is_continue);
	helical_nr_asu.deactivate(is_continue);
	helical_rise.deactivate(is_continue);
	helical_tube_outer_diameter.deactivate(is_continue);
	helical_bimodal_angular_priors.deactivate(is_continue);

}

bool ExtractJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{


	commands.clear();
	initialisePipeline(outputname, "Extract", job_counter);

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";

	// Input
	command += " --i " + star_mics.getValue();
	Node node(star_mics.getValue(), star_mics.type);
	pipelineInputNodes.push_back(node);


	if (do_reextract.getValue())
	{
		command += " --reextract_data_star " + fndata_reextract.getValue();
		Node node2(fndata_reextract.getValue(), fndata_reextract.type);
		pipelineInputNodes.push_back(node2);
		if (do_recenter.getValue())
		{
			command += " --recenter";
		}
	}
	else
	{
		FileName mysuffix = coords_suffix.getValue();
		command += " --coord_dir " + mysuffix.beforeLastOf("/") + "/";
		command += " --coord_suffix " + (mysuffix.afterLastOf("/")).without("coords_suffix");
		Node node2(coords_suffix.getValue(), coords_suffix.type);
		pipelineInputNodes.push_back(node2);
	}

	// Output
	FileName fn_ostar = outputname + "particles.star";
	Node node3(fn_ostar, NODE_PART_DATA);
	pipelineOutputNodes.push_back(node3);
	command += " --part_star " + fn_ostar;

	command += " --part_dir " + outputname;
	command += " --extract";
	command += " --extract_size " + floatToString(extract_size.getValue());

	// Operate stuff
	// Get an integer number for the bg_radius
	RFLOAT bg_radius = (bg_diameter.getValue() < 0.) ? 0.75 * extract_size.getValue() : bg_diameter.getValue();
	bg_radius /= 2.; // Go from diameter to radius
	if (do_rescale.getValue())
	{
		command += " --scale " + floatToString(rescale.getValue());
		bg_radius *= rescale.getValue() / extract_size.getValue();
	}
	if (do_norm.getValue())
	{
		// Get an integer number for the bg_radius
		bg_radius = (int)bg_radius;
		command += " --norm --bg_radius " + floatToString(bg_radius);
		command += " --white_dust " + floatToString(white_dust.getValue());
		command += " --black_dust " + floatToString(black_dust.getValue());
	}
	if (do_invert.getValue())
		command += " --invert_contrast ";

	if (do_set_angpix.getValue())
	{
		command += " --set_angpix " + floatToString(angpix.getValue());
	}

	// Helix
	if (do_extract_helix.getValue())
	{
		command += " --helix";
		command += " --helical_outer_diameter " + floatToString(helical_tube_outer_diameter.getValue());
		if (helical_bimodal_angular_priors.getValue())
			command += " --helical_bimodal_angular_priors";
		if (do_extract_helical_tubes.getValue())
		{
			command += " --helical_tubes";
			if (do_cut_into_segments.getValue())
			{
				command += " --helical_cut_into_segments";
				command += " --helical_nr_asu " + integerToString(helical_nr_asu.getValue());
				command += " --helical_rise " + floatToString(helical_rise.getValue());
			}
			else
				command += " --helical_nr_asu 1 --helical_rise 1";
		}
	}

	if (is_continue)
		command += " --only_extract_unfinished ";


	// Other arguments for extraction
	command += " " + other_args.getValue();
	commands.push_back(command);

	if (do_reextract.getValue() || (do_extract_helix.getValue() && do_extract_helical_tubes.getValue()) )
	{
		// Also touch the suffix file. Do this after the first command had completed
		command = "echo " + star_mics.getValue() + " > " +  outputname + "coords_suffix_extract.star";
		commands.push_back(command.c_str());

		Node node(outputname + "coords_suffix_extract.star", NODE_MIC_COORDS);
		pipelineOutputNodes.push_back(node);

	}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

SortJobWindow::SortJobWindow() : RelionJobWindow(2, HAS_MPI, HAS_NOT_THREAD)
{

	type = PROC_SORT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	input_star.place(current_y, "Input particles to be sorted:", NODE_PART_DATA, "", "Input particles(*.{star})", "This STAR file should contain in-plane rotations, in-plane translations and a class number that were obtained by alignment (class2D/class3D or auto3D) OR auto-picking. A column called rlnParticleSelectZScore will be added to this same STAR file with the sorting result. This column can then be used in the display programs to sort the particles on.");

	tab1->end();

	tab2->begin();
	tab2->label("CTF");
	resetHeight();

	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	do_ctf.place(current_y, "Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n ", ctf_group);

	ctf_group->begin();
	do_ignore_first_ctfpeak.place(current_y, "Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");
	ctf_group->end();
	do_ctf.cb_menu_i();

	tab2->end();

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
	do_ctf.writeValue(fh);
	do_ignore_first_ctfpeak.writeValue(fh);

	closeWriteFile(fh, fn);
}

void SortJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_sort";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		input_star.readValue(fh);
		do_ctf.readValue(fh);
		do_ignore_first_ctfpeak.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void SortJobWindow::toggle_new_continue(bool _is_continue)
{

	is_continue = _is_continue;

	input_star.deactivate(is_continue);
	do_ctf.deactivate(is_continue);
	do_ignore_first_ctfpeak.deactivate(is_continue);

}

bool SortJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Sort", job_counter);

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_particle_sort_mpi`";
	else
		command="`which relion_particle_sort`";

	command += " --i " + input_star.getValue();
	Node node(input_star.getValue(), input_star.type);
	pipelineInputNodes.push_back(node);

	// Determine the --ref automatically, from the particle input filename
	FileName fn_ref, fn_in = input_star.getValue();
	int node_type;
	if (fn_in.contains("_data.star") && (fn_in.contains("Class2D/") || fn_in.contains("Class3D/")) )
	{
		fn_ref = fn_in.without("_data.star") + "_model.star";
		node_type= NODE_MODEL;
	}
	else if (fn_in.contains("_data.star") && fn_in.contains("Refine3D/"))
	{
		if (fn_in.contains("_it0"))
			fn_ref = fn_in.without("_data.star") + "_half1_model.star";
		else
			fn_ref = fn_in.without("_data.star") + "_model.star";
		node_type= NODE_MODEL;
	}
	else if (fn_in.contains("Extract/"))
	{
		// TODO!
		REPORT_ERROR("ERROR: for the moment you can only use inputs from Class2D, Class3D or Refine3D runs!");
		node_type= NODE_2DREFS;
	}
	command += " --ref " + fn_ref;
	Node node2(fn_ref, node_type);
	pipelineInputNodes.push_back(node2);

	command += " --o " + outputname + "particles_sort.star";
	Node node3(outputname + "particles_sort.star", NODE_PART_DATA);
	pipelineOutputNodes.push_back(node3);

	if (do_ctf.getValue())
	{
		command += " --ctf ";
		if (do_ignore_first_ctfpeak.getValue())
			command += " --ctf_intact_first_peak ";
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);
}


Class2DJobWindow::Class2DJobWindow() : RelionJobWindow(6, HAS_MPI, HAS_THREAD)
{

	type = PROC_2DCLASS;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "CURRENT_ODIR",  "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

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

	nr_classes.place(current_y, "Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");

	// Add a little spacer
	current_y += STEPY/2;

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

	particle_diameter.place(current_y, "Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");

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

	// Add a little spacer
	current_y += STEPY/2;

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");


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
	tab5->begin();
	tab5->label("Helix");
	resetHeight();

	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	do_helix.place(current_y, "Classify 2D helical segments?", false, "Set to Yes if you want to classify 2D helical segments. Note that the helical segments should come with priors of psi angles", helix_group);

	helix_group->begin();

	helical_tube_outer_diameter.place(current_y, "Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of the tubes. You may want to copy the value from previous particle extraction job. \
If negative value is provided, this option is disabled and ordinary circular masks will be applied. Sometimes '--dont_check_norm' option is useful to prevent errors in normalisation of helical segments.");

	do_bimodal_psi.place(current_y, "Do bimodal angular searches?", true, "Do bimodal search for psi angles? \
Set to Yes if you want to classify 2D helical segments with priors of psi angles. The priors should be bimodal due to unknown polarities of the segments. \
Set to No if the 3D helix looks the same when rotated upside down. If it is set to No, ordinary angular searches will be performed.\n\nThis option will be invalid if you choose not to perform image alignment on 'Sampling' tab.");

	range_psi.place(current_y, "Angular search range - psi (deg):", 6, 3, 30, 1, "Local angular searches will be performed \
within +/- the given amount (in degrees) from the psi priors estimated through helical segment picking. \
A range of 15 degrees is the same as sigma = 5 degrees. Note that the ranges of angular searches should be much larger than the sampling.\
\n\nThis option will be invalid if you choose not to perform image alignment on 'Sampling' tab.");

	helix_group->end();
	do_helix.cb_menu_i(); // to make default effective

	tab5->end();


	tab6->begin();
	tab6->label("Compute");
	resetHeight();

	do_combine_thru_disc.place(current_y, "Combine iterations through disc?", true, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");

	nr_pool.place(current_y, "Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");

	do_preread_images.place(current_y, "Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in double-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 30Gb, or 120 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM.");


	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    gpu_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    gpu_group->end();
	use_gpu.place(current_y, "Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.", gpu_group);
	gpu_group->begin();
	gpu_ids.place(current_y, "Which GPUs to use:", "", "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. For example: '0,0:1,1:0,0:1,1'");
    gpu_group->end();
	use_gpu.cb_menu_i(); // This is to make the default effective

	tab6->end();


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
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);
	nr_classes.writeValue(fh);

	// CTF
	do_ctf_correction.writeValue(fh);
	ctf_phase_flipped.writeValue(fh);
	ctf_intact_first_peak.writeValue(fh);

	// Optimisation
	nr_iter.writeValue(fh);
	tau_fudge.writeValue(fh);
	particle_diameter.writeValue(fh);
	do_zero_mask.writeValue(fh);
	highres_limit.writeValue(fh);

	// Sampling
	dont_skip_align.writeValue(fh);
	psi_sampling.writeValue(fh);
	offset_range.writeValue(fh);
	offset_step.writeValue(fh);

	// Helix
	do_helix.writeValue(fh);
	do_bimodal_psi.writeValue(fh);
	range_psi.writeValue(fh);
	helical_tube_outer_diameter.writeValue(fh);

	// Compute
	do_combine_thru_disc.writeValue(fh);
	do_parallel_discio.writeValue(fh);
	nr_pool.writeValue(fh);
	do_preread_images.writeValue(fh);
	use_gpu.writeValue(fh);
	gpu_ids.writeValue(fh);

	closeWriteFile(fh, fn);
}

void Class2DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_class2d";

	if (openReadFile(fn, fh))
	{

		// I/O
		fn_cont.readValue(fh);
		fn_img.readValue(fh);
		nr_classes.readValue(fh);

		// CTF
		do_ctf_correction.readValue(fh);
		ctf_phase_flipped.readValue(fh);
		ctf_intact_first_peak.readValue(fh);

		// Optimisation
		nr_iter.readValue(fh);
		tau_fudge.readValue(fh);
		particle_diameter.readValue(fh);
		do_zero_mask.readValue(fh);
		highres_limit.readValue(fh);

		// Sampling
		dont_skip_align.readValue(fh);
		psi_sampling.readValue(fh);
		offset_range.readValue(fh);
		offset_step.readValue(fh);

		// Helix
		do_helix.readValue(fh);
		do_bimodal_psi.readValue(fh);
		range_psi.readValue(fh);
		helical_tube_outer_diameter.readValue(fh);

		// Compute
		do_combine_thru_disc.readValue(fh);
		do_parallel_discio.readValue(fh);
		nr_pool.readValue(fh);
		do_preread_images.readValue(fh);
		use_gpu.readValue(fh);
		gpu_ids.readValue(fh);

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

bool Class2DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Class2D", job_counter);

	std::string command;
	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		int pos_it = fn_cont.getValue().rfind("_it");
		int pos_op = fn_cont.getValue().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << fn_cont.getValue() << std::endl;
		int it = (int)textToFloat((fn_cont.getValue().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + fn_cont.getValue();
    }

    command += " --o " + outputname + fn_run;
	pipelineOutputNodes = getOutputNodesRefine(outputname + fn_run, nr_iter.getValue(), nr_classes.getValue(), 2, 1);

	if (!is_continue)
	{
		command += " --i " + fn_img.getValue();
		Node node(fn_img.getValue(), fn_img.type);
		pipelineInputNodes.push_back(node);
	}

	// Always do compute stuff
	if (!do_combine_thru_disc.getValue())
		command += " --dont_combine_weights_via_disc";
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";
	if (do_preread_images.getValue())
		command += " --preread_images --pool 1 " ;
	else
		command += " --pool " + floatToString(nr_pool.getValue());

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
    command += " --particle_diameter " + floatToString(particle_diameter.getValue());
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

	// Helix
	if (do_helix.getValue())
	{
		command += " --helical_outer_diameter " + floatToString(helical_tube_outer_diameter.getValue());

		if (dont_skip_align.getValue())
		{
			if (do_bimodal_psi.getValue())
				command += " --bimodal_psi";

			RFLOAT val = range_psi.getValue();
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
		}
	}

	// Always do norm and scale correction
	if (!is_continue)
		command += " --norm --scale ";

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());

	// GPU-stuff
	if (use_gpu.getValue())
	{
		command += " --gpu " + gpu_ids.getValue();
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

Class3DJobWindow::Class3DJobWindow() : RelionJobWindow(7, HAS_MPI, HAS_THREAD)
{

	type = PROC_3DCLASS;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	fn_ref.place(current_y, "Reference map:", NODE_3DREF, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");

	fn_mask.place(current_y, "Reference mask (optional):", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "\
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

	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();


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

	nr_classes.place(current_y, "Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");

	// Add a little spacer
	current_y += STEPY/2;

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

	particle_diameter.place(current_y, "Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");

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

	tab6->begin();
	tab6->label("Helix");
	resetHeight();
	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	//helix_text.place(current_y, "Nov 21, 2015");
	do_helix.place(current_y, "Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.", helix_group);
	helix_group->begin();
	helical_tube_inner_diameter.placeOnSameYPosition(current_y, "Tube diameter - inner, outer (A):", "Tube diameter - inner (A):", "-1", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	helical_tube_outer_diameter.placeOnSameYPosition(current_y, "", "Tube diameter - outer (A):", "-1", "Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;
	helical_nr_asu.place(current_y, "Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.");
	helical_twist_initial.placeOnSameYPosition(current_y, "Initial twist (deg), rise (A):", "Initial helical twist (deg):", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	helical_rise_initial.placeOnSameYPosition(current_y, "", "Initial helical rise (A):", "0", "Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;

	// Add a little spacer
	current_y += STEPY/2;

	helix_symmetry_search_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_symmetry_search_group->end();
	do_local_search_helical_symmetry.place(current_y, "Do local searches of symmetry?", true, "If set to Yes, then perform local searches of helical twist and rise within given ranges.", helix_symmetry_search_group);
	helix_symmetry_search_group->begin();
	helical_twist_min.placeOnSameYPosition(current_y, "Twist search - Min, Max, Step (deg):", "Helical twist search (deg) - Min:", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_twist_max.placeOnSameYPosition(current_y, "", "Helical twist search (deg) - Max:", "0", NULL, XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_twist_inistep.placeOnSameYPosition(current_y, "", "Helical twist search (deg) - Step:", "0", "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.", XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y += STEPY + 2;
	helical_rise_min.placeOnSameYPosition(current_y, "Rise search - Min, Max, Step (A):", "Helical rise search (A) - Min:", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_rise_max.placeOnSameYPosition(current_y, "", "Helical rise search (A) - Max:", "0", NULL, XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_rise_inistep.placeOnSameYPosition(current_y, "", "Helical rise search (A) - Step:", "0", "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.", XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y += STEPY + 2;
	helix_symmetry_search_group->end();
	do_local_search_helical_symmetry.cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	helical_z_percentage.place(current_y, "Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis, where Fourier artefacts are also present if the \
number of helical asymmetrical units is larger than 1. Therefore, information from the central part of the box is used for searching and imposing \
helical symmetry in real space. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used.");
	range_tilt.placeOnSameYPosition(current_y, "Angular search range - tilt, psi (deg):", "Angular search range - tilt (deg):", "15", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	range_psi.placeOnSameYPosition(current_y, "", "Angular search range - psi (deg):", "15", "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;
	helical_range_distance.place(current_y, "Range factor of local averaging:", -1., 1., 5., 0.1, "Local averaging of orientations and translations will be performed within a range of +/- this value * the box size. Polarities are also set to be the same for segments coming from the same tube during local refinement. \
Values of ~ 2.0 are recommended for flexible structures such as MAVS-CARD filaments, ParM, MamK, etc. This option might not improve the reconstructions of helices formed from curled 2D lattices (TMV and VipA/VipB). Set to negative to disable this option.");
	helix_group->end();
	do_helix.cb_menu_i(); // to make default effective
	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	do_combine_thru_disc.place(current_y, "Combine iterations through disc?", true, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");

	nr_pool.place(current_y, "Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");

	do_preread_images.place(current_y, "Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in double-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 30Gb, or 120 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM.");


	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    gpu_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    gpu_group->end();
	use_gpu.place(current_y, "Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.", gpu_group);
	gpu_group->begin();
	gpu_ids.place(current_y, "Which GPUs to use:", "", "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");
    gpu_group->end();
	use_gpu.cb_menu_i(); // This is to make the default effective

	tab7->end();


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
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);
	nr_classes.writeValue(fh);

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
	particle_diameter.writeValue(fh);
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

	// Helix
	do_helix.writeValue(fh);
	helical_tube_inner_diameter.writeValue(fh);
	helical_tube_outer_diameter.writeValue(fh);
	helical_nr_asu.writeValue(fh);
	helical_twist_initial.writeValue(fh);
	helical_rise_initial.writeValue(fh);
	do_local_search_helical_symmetry.writeValue(fh);
	helical_twist_min.writeValue(fh);
	helical_twist_max.writeValue(fh);
	helical_twist_inistep.writeValue(fh);
	helical_rise_min.writeValue(fh);
	helical_rise_max.writeValue(fh);
	helical_rise_inistep.writeValue(fh);
	helical_z_percentage.writeValue(fh);
	range_tilt.writeValue(fh);
	range_psi.writeValue(fh);
	helical_range_distance.writeValue(fh);

	// Compute
	do_combine_thru_disc.writeValue(fh);
	do_parallel_discio.writeValue(fh);
	nr_pool.writeValue(fh);
	do_preread_images.writeValue(fh);
	use_gpu.writeValue(fh);
	gpu_ids.writeValue(fh);

	closeWriteFile(fh, fn);
}

void Class3DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_class3d";

	if (openReadFile(fn, fh))
	{

		// I/O
		fn_cont.readValue(fh);
		fn_img.readValue(fh);
		nr_classes.readValue(fh);

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
		particle_diameter.readValue(fh);
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

		// Helix
		do_helix.readValue(fh);
		helical_tube_inner_diameter.readValue(fh);
		helical_tube_outer_diameter.readValue(fh);
		helical_nr_asu.readValue(fh);
		helical_twist_initial.readValue(fh);
		helical_rise_initial.readValue(fh);
		do_local_search_helical_symmetry.readValue(fh);
		helical_twist_min.readValue(fh);
		helical_twist_max.readValue(fh);
		helical_twist_inistep.readValue(fh);
		helical_rise_min.readValue(fh);
		helical_rise_max.readValue(fh);
		helical_rise_inistep.readValue(fh);
		helical_z_percentage.readValue(fh);
		range_tilt.readValue(fh);
		range_psi.readValue(fh);
		helical_range_distance.readValue(fh);

		// Compute
		do_combine_thru_disc.readValue(fh);
		do_parallel_discio.readValue(fh);
		nr_pool.readValue(fh);
		do_preread_images.readValue(fh);
		use_gpu.readValue(fh);
		gpu_ids.readValue(fh);

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

	// Helix
	do_helix.deactivate(is_continue);
	helical_tube_inner_diameter.deactivate(is_continue);
	helical_tube_outer_diameter.deactivate(is_continue);
	helical_nr_asu.deactivate(is_continue);
	helical_twist_initial.deactivate(is_continue);
	helical_rise_initial.deactivate(is_continue);
	do_local_search_helical_symmetry.deactivate(is_continue);
	helical_twist_min.deactivate(is_continue);
	helical_twist_max.deactivate(is_continue);
	helical_twist_inistep.deactivate(is_continue);
	helical_rise_min.deactivate(is_continue);
	helical_rise_max.deactivate(is_continue);
	helical_rise_inistep.deactivate(is_continue);
	helical_z_percentage.deactivate(is_continue);
	range_tilt.deactivate(is_continue);
	range_psi.deactivate(is_continue);
	helical_range_distance.deactivate(is_continue);

	// Sampling


}

bool Class3DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	std::string command;

	initialisePipeline(outputname, "Class3D", job_counter);

	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		int pos_it = fn_cont.getValue().rfind("_it");
		int pos_op = fn_cont.getValue().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << fn_cont.getValue() << std::endl;
		int it = (int)textToFloat((fn_cont.getValue().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + fn_cont.getValue();
    }

    command += " --o " + outputname + fn_run;
	pipelineOutputNodes = getOutputNodesRefine(outputname + fn_run, nr_iter.getValue(), nr_classes.getValue(), 3, 1);

	if (!is_continue)
	{
		command += " --i " + fn_img.getValue();
		Node node(fn_img.getValue(), fn_img.type);
		pipelineInputNodes.push_back(node);
		if (fn_ref.getValue() != "None")
		{
			command += " --ref " + fn_ref.getValue();
			Node node(fn_ref.getValue(), fn_ref.type);
			pipelineInputNodes.push_back(node);
		}
		if (!ref_correct_greyscale.getValue() && fn_ref.getValue() != "None") // dont do firstiter_cc when giving None
			command += " --firstiter_cc";
		if (ini_high.getValue() > 0.)
			command += " --ini_high " + floatToString(ini_high.getValue());

	}

	// Always do compute stuff
	if (!do_combine_thru_disc.getValue())
		command += " --dont_combine_weights_via_disc";
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";
	if (do_preread_images.getValue())
		command += " --preread_images --pool 1 " ;
	else
		command += " --pool " + floatToString(nr_pool.getValue());

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
    command += " --particle_diameter " + floatToString(particle_diameter.getValue());
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
		Node node(fn_mask.getValue(), fn_mask.type);
		pipelineInputNodes.push_back(node);
	}

	// Sampling
	if (!dont_skip_align.getValue())
	{
		command += " --skip_align ";
	}
	else
	{
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

	if ( (!is_continue) && (do_helix.getValue()) )
	{
		command += " --helix";
		if (textToFloat(helical_tube_inner_diameter.getValue()) > 0.)
			command += " --helical_inner_diameter " + floatToString(textToFloat(helical_tube_inner_diameter.getValue()));
		command += " --helical_outer_diameter " + floatToString(textToFloat(helical_tube_outer_diameter.getValue()));
		command += " --helical_z_percentage " + floatToString(helical_z_percentage.getValue() / 100.);
		command += " --helical_nr_asu " + integerToString(helical_nr_asu.getValue());
		command += " --helical_twist_initial " + floatToString(textToFloat(helical_twist_initial.getValue()));
		command += " --helical_rise_initial " + floatToString(textToFloat(helical_rise_initial.getValue()));
		if (do_local_search_helical_symmetry.getValue())
		{
			command += " --helical_symmetry_search";
			command += " --helical_twist_min " + floatToString(textToFloat(helical_twist_min.getValue()));
			command += " --helical_twist_max " + floatToString(textToFloat(helical_twist_max.getValue()));
			if (textToFloat(helical_twist_inistep.getValue()) > 0.)
				command += " --helical_twist_inistep " + floatToString(textToFloat(helical_twist_inistep.getValue()));
			command += " --helical_rise_min " + floatToString(textToFloat(helical_rise_min.getValue()));
			command += " --helical_rise_max " + floatToString(textToFloat(helical_rise_max.getValue()));
			if (textToFloat(helical_rise_inistep.getValue()) > 0.)
				command += " --helical_rise_inistep " + floatToString(textToFloat(helical_rise_inistep.getValue()));
		}
		if ( (dont_skip_align.getValue()) && (!do_local_ang_searches.getValue()) )
		{
			RFLOAT val;
			val = textToFloat(range_tilt.getValue());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_tilt " + floatToString(val / 3.);
			val = textToFloat(range_psi.getValue());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
			if (helical_range_distance.getValue() > 0.)
				command += " --helical_sigma_distance " + floatToString(helical_range_distance.getValue() / 3.);
		}
	}

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());

	// GPU-stuff
	if (use_gpu.getValue())
	{
		command += " --gpu " + gpu_ids.getValue();
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

Auto3DJobWindow::Auto3DJobWindow() : RelionJobWindow(7, HAS_MPI, HAS_THREAD)
{

	type = PROC_3DAUTO;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_img.place(current_y, "Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	fn_cont.place(current_y, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	fn_ref.place(current_y, "Reference map:", NODE_3DREF, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");

	fn_mask.place(current_y, "Reference mask (optional):", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "\
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


	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();

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
	do_ctf_correction.place(current_y, "Do CTF-correction?", true, "If set to Yes, CTFs will be applied to the projections of the map. This requires that CTF information is present in the input STAR file.", ctf_group);

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

	particle_diameter.place(current_y, "Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");

	do_zero_mask.place(current_y, "Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");

	// Add a little spacer
	current_y += STEPY/2;

	do_solvent_fsc.place(current_y, "Use solvent-flattened FSCs?", false, "If set to Yes, then instead of using unmasked maps to calculate the gold-standard FSCs during refinement, \
masked half-maps are used and a post-processing-like correction of the FSC curves (with phase-randomisation) is performed every iteration. This only works when a reference mask is provided on the I/O tab. \
This may yield higher-resolution maps, especially when the mask contains only a relatively small volume inside the box.");

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
	tab6->begin();
	tab6->label("Helix");
	resetHeight();
	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	//helix_text.place(current_y, "Nov 21, 2015");
	do_helix.place(current_y, "Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.", helix_group);
	helix_group->begin();
	helical_tube_inner_diameter.placeOnSameYPosition(current_y, "Tube diameter - inner, outer (A):", "Tube diameter - inner (A):", "-1", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	helical_tube_outer_diameter.placeOnSameYPosition(current_y, "", "Tube diameter - outer (A):", "-1", "Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;
	helical_nr_asu.place(current_y, "Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.");
	helical_twist_initial.placeOnSameYPosition(current_y, "Initial twist (deg), rise (A):", "Initial helical twist (deg):", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	helical_rise_initial.placeOnSameYPosition(current_y, "", "Initial helical rise (A):", "0", "Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;

	// Add a little spacer
	current_y += STEPY/2;

	helix_symmetry_search_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_symmetry_search_group->end();
	do_local_search_helical_symmetry.place(current_y, "Do local searches of symmetry?", true, "If set to Yes, then perform local searches of helical twist and rise within given ranges.", helix_symmetry_search_group);
	helix_symmetry_search_group->begin();
	helical_twist_min.placeOnSameYPosition(current_y, "Twist search - Min, Max, Step (deg):", "Helical twist search (deg) - Min:", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_twist_max.placeOnSameYPosition(current_y, "", "Helical twist search (deg) - Max:", "0", NULL, XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_twist_inistep.placeOnSameYPosition(current_y, "", "Helical twist search (deg) - Step:", "0", "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.", XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y += STEPY + 2;
	helical_rise_min.placeOnSameYPosition(current_y, "Rise search - Min, Max, Step (A):", "Helical rise search (A) - Min:", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_rise_max.placeOnSameYPosition(current_y, "", "Helical rise search (A) - Max:", "0", NULL, XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	helical_rise_inistep.placeOnSameYPosition(current_y, "", "Helical rise search (A) - Step:", "0", "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.", XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y += STEPY + 2;
	helix_symmetry_search_group->end();

	// Add a little spacer
	current_y += STEPY/2;

	do_local_search_helical_symmetry.cb_menu_i(); // to make default effective
	helical_z_percentage.place(current_y, "Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis, where Fourier artefacts are also present if the \
number of helical asymmetrical units is larger than 1. Therefore, information from the central part of the box is used for searching and imposing \
helical symmetry in real space. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used.");
	range_tilt.placeOnSameYPosition(current_y, "Angular search range - tilt, psi (deg):", "Angular search range - tilt (deg):", "15", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	range_psi.placeOnSameYPosition(current_y, "", "Angular search range - psi (deg):", "15", "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y += STEPY + 2;
	helical_range_distance.place(current_y, "Range factor of local averaging:", -1., 1., 5., 0.1, "Local averaging of orientations and translations will be performed within a range of +/- this value * the box size. Polarities are also set to be the same for segments coming from the same tube during local refinement. \
Values of ~ 2.0 are recommended for flexible structures such as MAVS-CARD filaments, ParM, MamK, etc. This option might not improve the reconstructions of helices formed from curled 2D lattices (TMV and VipA/VipB). Set to negative to disable this option.");
	helix_group->end();
	do_helix.cb_menu_i(); // to make default effective

	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	do_combine_thru_disc.place(current_y, "Combine iterations through disc?", true, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	do_parallel_discio.place(current_y, "Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");

	nr_pool.place(current_y, "Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");

	do_preread_images.place(current_y, "Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in double-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 30Gb, or 120 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM.");


	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    gpu_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    gpu_group->end();
	use_gpu.place(current_y, "Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.", gpu_group);
	gpu_group->begin();
	gpu_ids.place(current_y, "Which GPUs to use:", "", "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");
    gpu_group->end();
	use_gpu.cb_menu_i(); // This is to make the default effective

	tab7->end();

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
	fn_cont.writeValue(fh);
	fn_img.writeValue(fh);

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
	particle_diameter.writeValue(fh);
	do_zero_mask.writeValue(fh);
	fn_mask.writeValue(fh);
	do_solvent_fsc.writeValue(fh);

	// Sampling
	sampling.writeValue(fh);
	offset_range.writeValue(fh);
	offset_step.writeValue(fh);
	auto_local_sampling.writeValue(fh);

	// Helix
	do_helix.writeValue(fh);
	helical_tube_inner_diameter.writeValue(fh);
	helical_tube_outer_diameter.writeValue(fh);
	helical_nr_asu.writeValue(fh);
	helical_twist_initial.writeValue(fh);
	helical_rise_initial.writeValue(fh);
	do_local_search_helical_symmetry.writeValue(fh);
	helical_twist_min.writeValue(fh);
	helical_twist_max.writeValue(fh);
	helical_twist_inistep.writeValue(fh);
	helical_rise_min.writeValue(fh);
	helical_rise_max.writeValue(fh);
	helical_rise_inistep.writeValue(fh);
	helical_z_percentage.writeValue(fh);
	range_tilt.writeValue(fh);
	range_psi.writeValue(fh);
	helical_range_distance.writeValue(fh);

	// Compute
	do_combine_thru_disc.writeValue(fh);
	do_parallel_discio.writeValue(fh);
	nr_pool.writeValue(fh);
	do_preread_images.writeValue(fh);
	use_gpu.writeValue(fh);
	gpu_ids.writeValue(fh);

	closeWriteFile(fh, fn);
}

void Auto3DJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_auto3d";

	if (openReadFile(fn, fh))
	{

		// I/O
		fn_cont.readValue(fh);
		fn_img.readValue(fh);

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
		particle_diameter.readValue(fh);
		do_zero_mask.readValue(fh);
		fn_mask.readValue(fh);
		do_solvent_fsc.readValue(fh);

		// Sampling
		sampling.readValue(fh);
		offset_range.readValue(fh);
		offset_step.readValue(fh);
		auto_local_sampling.readValue(fh);

		// Helix
		do_helix.readValue(fh);
		helical_tube_inner_diameter.readValue(fh);
		helical_tube_outer_diameter.readValue(fh);
		helical_nr_asu.readValue(fh);
		helical_twist_initial.readValue(fh);
		helical_rise_initial.readValue(fh);
		do_local_search_helical_symmetry.readValue(fh);
		helical_twist_min.readValue(fh);
		helical_twist_max.readValue(fh);
		helical_twist_inistep.readValue(fh);
		helical_rise_min.readValue(fh);
		helical_rise_max.readValue(fh);
		helical_rise_inistep.readValue(fh);
		helical_z_percentage.readValue(fh);
		range_tilt.readValue(fh);
		range_psi.readValue(fh);
		helical_range_distance.readValue(fh);

		// Compute
		do_combine_thru_disc.readValue(fh);
		do_parallel_discio.readValue(fh);
		nr_pool.readValue(fh);
		do_preread_images.readValue(fh);
		use_gpu.readValue(fh);
		gpu_ids.readValue(fh);

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

	// Helix
	do_helix.deactivate(is_continue);
	helical_tube_inner_diameter.deactivate(is_continue);
	helical_tube_outer_diameter.deactivate(is_continue);
	helical_nr_asu.deactivate(is_continue);
	helical_twist_initial.deactivate(is_continue);
	helical_rise_initial.deactivate(is_continue);
	do_local_search_helical_symmetry.deactivate(is_continue);
	helical_twist_min.deactivate(is_continue);
	helical_twist_max.deactivate(is_continue);
	helical_twist_inistep.deactivate(is_continue);
	helical_rise_min.deactivate(is_continue);
	helical_rise_max.deactivate(is_continue);
	helical_rise_inistep.deactivate(is_continue);
	helical_z_percentage.deactivate(is_continue);
	range_tilt.deactivate(is_continue);
	range_psi.deactivate(is_continue);
	helical_range_distance.deactivate(is_continue);
}

bool Auto3DJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	std::string command;

	initialisePipeline(outputname, "Refine3D", job_counter);

	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		int pos_it = fn_cont.getValue().rfind("_it");
		int pos_op = fn_cont.getValue().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << fn_cont.getValue() << std::endl;
		int it = (int)textToFloat((fn_cont.getValue().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + fn_cont.getValue();
    }

    command += " --o " + outputname + fn_run;
	// TODO: add bodies!! (probably in next version)
	pipelineOutputNodes = getOutputNodesRefine(outputname + fn_run, -1, 1, 3, 1, false, false); // false false means dont do movies

	if (!is_continue)
	{
		command += " --auto_refine --split_random_halves --i " + fn_img.getValue();
		Node node(fn_img.getValue(), fn_img.type);
		pipelineInputNodes.push_back(node);
		if (fn_ref.getValue() != "None")
		{
			command += " --ref " + fn_ref.getValue();
			Node node(fn_ref.getValue(), fn_ref.type);
			pipelineInputNodes.push_back(node);
		}
		if (!ref_correct_greyscale.getValue() && fn_ref.getValue() != "None") // dont do firstiter_cc when giving None
			command += " --firstiter_cc";
		if (ini_high.getValue() > 0.)
			command += " --ini_high " + floatToString(ini_high.getValue());

	}

	// Always do compute stuff
	if (!do_combine_thru_disc.getValue())
		command += " --dont_combine_weights_via_disc";
	if (!do_parallel_discio.getValue())
		command += " --no_parallel_disc_io";
	if (do_preread_images.getValue())
		command += " --preread_images --pool 1 " ;
	else
		command += " --pool " + floatToString(nr_pool.getValue());

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
    command += " --particle_diameter " + floatToString(particle_diameter.getValue());
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

		if (do_solvent_fsc.getValue())
			command += " --solvent_correct_fsc ";

		// TODO: what if this is a continuation run: re-put the mask as an input node? Or only if it changes? Also for 3Dclass
		Node node(fn_mask.getValue(), fn_mask.type);
		pipelineInputNodes.push_back(node);
	}

	if (!is_continue)
	{
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

		command += " --sym " + sym_name.getValue();
                // Always join low-res data, as some D&I point group refinements may fall into different hands!
                command += " --low_resol_join_halves 40";
		command += " --norm --scale ";

		// Helix
		if (do_helix.getValue())
		{
			command += " --helix";
			if (textToFloat(helical_tube_inner_diameter.getValue()) > 0.)
				command += " --helical_inner_diameter " + floatToString(textToFloat(helical_tube_inner_diameter.getValue()));
			command += " --helical_outer_diameter " + floatToString(textToFloat(helical_tube_outer_diameter.getValue()));
			command += " --helical_z_percentage " + floatToString(helical_z_percentage.getValue() / 100.);
			command += " --helical_nr_asu " + integerToString(helical_nr_asu.getValue());
			command += " --helical_twist_initial " + floatToString(textToFloat(helical_twist_initial.getValue()));
			command += " --helical_rise_initial " + floatToString(textToFloat(helical_rise_initial.getValue()));
			if (do_local_search_helical_symmetry.getValue())
			{
				command += " --helical_symmetry_search";
				command += " --helical_twist_min " + floatToString(textToFloat(helical_twist_min.getValue()));
				command += " --helical_twist_max " + floatToString(textToFloat(helical_twist_max.getValue()));
				if (textToFloat(helical_twist_inistep.getValue()) > 0.)
					command += " --helical_twist_inistep " + floatToString(textToFloat(helical_twist_inistep.getValue()));
				command += " --helical_rise_min " + floatToString(textToFloat(helical_rise_min.getValue()));
				command += " --helical_rise_max " + floatToString(textToFloat(helical_rise_max.getValue()));
				if (textToFloat(helical_rise_inistep.getValue()) > 0.)
					command += " --helical_rise_inistep " + floatToString(textToFloat(helical_rise_inistep.getValue()));
			}
			RFLOAT val;
			val = textToFloat(range_tilt.getValue());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_tilt " + floatToString(val / 3.);
			val = textToFloat(range_psi.getValue());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
			if (helical_range_distance.getValue() > 0.)
				command += " --helical_sigma_distance " + floatToString(helical_range_distance.getValue() / 3.);
		}
	}

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());

	// GPU-stuff
	if (use_gpu.getValue())
	{
		command += " --gpu " + gpu_ids.getValue();
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

MovieRefineJobWindow::MovieRefineJobWindow() : RelionJobWindow(4, HAS_MPI, HAS_THREAD)
{

	type = PROC_MOVIEREFINE;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();


	fn_cont.place(current_y, "Continue 3D auto-refine from: ", "", "STAR Files (*_optimiser.star)", "Refine3D/", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run. \n \
Besides restarting jobs that were somehow stopped before convergence, also use the continue-option after the last iteration to do movie processing.");
    movie_rootname.place(current_y, "Rootname of movies files:", "movie", "rootname to relate each movie to the single-frame averaged micropgraph. With a rootname of 'movie', the movie for mic001.mrc should be called mic001_movie.mrcs. If you've run the MOTIONCORR wrapper in RELION, then the correct rootname is 'movie'.");

    // Add a little spacer
	current_y += STEPY/2;

	join_nr_mics.place(current_y, "Process micrographs in batches of: ", 50, 10, 200, 10, "All movie-particles will be extracted from each micrograph separately, but the resulting STAR files will be joined together into batches of the size specified here. The movie-refinement will be performed on these batches. \
Very large batches cost more RAM, but the parallelisation in smaller batches is poorer. If a negative number is given, all micrographs are processed in one batch. Note that movie-refinement that INCLUDE rotational searches cannot be performed in batches!");

	tab1->end();

	tab2->begin();
	tab2->label("extract");
	resetHeight();

	first_movie_frame.place(current_y, "First movie frame to extract: ", 1, 1, 20, 1, "Extract from this movie frame onwards. The first frame is number 1.");
	last_movie_frame.place(current_y, "Last movie frame to extract: ", 0, 0, 64, 1, "Extract until this movie frame. Zero means: extract all frames in the movie. You may want to specify the last frame number though, as it will be useful to detect movies which accidentally have fewer frames.");
	avg_movie_frames.place(current_y, "Average every so many frames: ", 1, 1, 8, 1, "Average every so many movie frames together upon the extraction. For example, 32-frame movies may be reduced to 16-frame movie-particles when provding a value of 2 here. This will reduce computational costs in movie-refinement and polishing, but too large values will affect the results. Default is a value of 1, so no averaging");
	max_mpi_nodes.place(current_y, "Maximum number of MPI nodes: ", 8, 2, 24, 1, "The number of MPI nodes used by the relion_preprocess program will be limited to this value, regardless of the number of MPI nodes requested on the Running tab (which is also used for the refinement step). This is useful to protect the file system from too heavy disk I/O.");

	// Add a little spacer
	current_y += STEPY/2;

	extract_size.place(current_y,"Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). Use the same as for the refined particles!");
	do_invert.place(current_y, "Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted. Use the same as for the refined particles!");

	// Add a little spacer
	current_y += STEPY/2;

	rescale_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	rescale_group->end();
	do_rescale.place(current_y, "Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.", rescale_group);
	rescale_group->begin();
	rescale.place(current_y, "Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number. Use the same as for the refined particles! ");
	rescale_group->end();
	do_rescale.cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("normalise");
	resetHeight();

	norm_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	norm_group->end();
	do_norm.place(current_y, "Normalize movie-particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it. Use the same values as for the refined particles!", norm_group);

	norm_group->begin();


	bg_diameter.place(current_y, "Diameter background circle (pix): ", -1, -1, 600, 10, " Use the same as for the refined particles! Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");

	white_dust.place(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, " Use the same as for the refined particles! Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	black_dust.place(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, " Use the same as for the refined particles! Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	norm_group->end();
	do_norm.cb_menu_i();



	tab3->end();
	tab4->begin();
	tab4->label("refine");
	resetHeight();

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

	tab4->end();

	// read settings if hidden file exists
	read(".gui_movierefine", is_continue);

}

void MovieRefineJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_movierefine";

	std::ofstream fh;
	openWriteFile(fn, fh);

	// I/O
	fn_cont.writeValue(fh);
	movie_rootname.writeValue(fh);
	join_nr_mics.writeValue(fh);

	// Extract movie-particles
	first_movie_frame.writeValue(fh);
	last_movie_frame.writeValue(fh);
	avg_movie_frames.writeValue(fh);
	max_mpi_nodes.writeValue(fh);
	extract_size.writeValue(fh);
	do_rescale.writeValue(fh);
	rescale.writeValue(fh);
	do_norm.writeValue(fh);
	bg_diameter.writeValue(fh);
	white_dust.writeValue(fh);
	black_dust.writeValue(fh);
	do_invert.writeValue(fh);

	// Refine Movies
	movie_runavg_window.writeValue(fh);
	movie_sigma_offset.writeValue(fh);
	do_alsorot_movies.writeValue(fh);
	movie_sigma_angles.writeValue(fh);

	closeWriteFile(fh, fn);
}

void MovieRefineJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_movierefine";

	if (openReadFile(fn, fh))
	{

		// I/O
		fn_cont.readValue(fh);
		movie_rootname.readValue(fh);
		join_nr_mics.readValue(fh);

		// Extract movie-particles
		first_movie_frame.readValue(fh);
		last_movie_frame.readValue(fh);
		avg_movie_frames.readValue(fh);
		max_mpi_nodes.readValue(fh);
		extract_size.readValue(fh);
		do_rescale.readValue(fh);
		rescale.readValue(fh);
		do_norm.readValue(fh);
		bg_diameter.readValue(fh);
		white_dust.readValue(fh);
		black_dust.readValue(fh);
		do_invert.readValue(fh);

		// Refine Movies
		movie_runavg_window.readValue(fh);
		movie_sigma_offset.readValue(fh);
		do_alsorot_movies.readValue(fh);
		movie_sigma_angles.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}

void MovieRefineJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_cont.deactivate(is_continue);
	movie_rootname.deactivate(is_continue);
	join_nr_mics.deactivate(is_continue);

	// Extract movie-particles
	first_movie_frame.deactivate(is_continue);
	last_movie_frame.deactivate(is_continue);
	avg_movie_frames.deactivate(is_continue);
	max_mpi_nodes.deactivate(is_continue);
	extract_size.deactivate(is_continue);
	do_rescale.deactivate(is_continue);
	rescale.deactivate(is_continue);
	do_norm.deactivate(is_continue);
	bg_diameter.deactivate(is_continue);
	white_dust.deactivate(is_continue);
	black_dust.deactivate(is_continue);
	do_invert.deactivate(is_continue);

	// Movies: allow changing parameters here!
	// Make specific run-names for different variables of ravg, sigma_offset and sigma_angles,
	// This way you don't have to re-extract all movie particles each time you try a different parameter setting
	//movie_runavg_window.deactivate(is_continue);
	//movie_sigma_offset.deactivate(is_continue);
	//do_alsorot_movies.deactivate(is_continue);
	//movie_sigma_angles.deactivate(is_continue);

}

bool MovieRefineJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	std::string command;

	initialisePipeline(outputname, "MovieRefine", job_counter);

	// A. First get the extract command
	if (nr_mpi.getValue() > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";

	// Input
	// Get the data.star to be used for re-extraction from the optimiser name
	MetaDataTable MDopt;
	MDopt.read(fn_cont.getValue(), "optimiser_general");
	FileName fn_data;
	MDopt.getValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data);
	command += " --reextract_data_star " + fn_data;
	Node node2(fn_cont.getValue(), NODE_OPTIMISER);
	pipelineInputNodes.push_back(node2);

	// Output files of the extraction: to be used in the second step: movie-refinement
	command += " --part_dir " + outputname;

	FileName fn_ostar = outputname + "particles_" + movie_rootname.getValue() + ".star";
	FileName fn_olist = outputname + "micrographs_" + movie_rootname.getValue() + "_list.star";
	command += " --list_star " + fn_olist;
	command += " --join_nr_mics " + floatToString(join_nr_mics.getValue());
	if (join_nr_mics.getValue() <= 0)
	{
		// Only write out a STAR file with all movie-particles if we're NOT processing on a per-micrograph basis
		command += " --part_star " + fn_ostar;
	}

	// Extraction parameters
	command += " --extract --extract_movies";
	command += " --extract_size " + floatToString(extract_size.getValue());
	command += " --movie_name " + movie_rootname.getValue();
	command += " --first_movie_frame " + floatToString(first_movie_frame.getValue());
	command += " --last_movie_frame " + floatToString(last_movie_frame.getValue());
	command += " --avg_movie_frames " + floatToString(avg_movie_frames.getValue());
	// Limit MPI nodes
	if (nr_mpi.getValue() > ROUND(max_mpi_nodes.getValue()))
		command += " --max_mpi_nodes " + floatToString(max_mpi_nodes.getValue());

	// Operate stuff
	// Get an integer number for the bg_radius
	RFLOAT bg_radius = (bg_diameter.getValue() < 0.) ? 0.75 * extract_size.getValue() : bg_diameter.getValue();
	bg_radius /= 2.; // Go from diameter to radius
	if (do_rescale.getValue())
	{
		command += " --scale " + floatToString(rescale.getValue());
		bg_radius *= rescale.getValue() / extract_size.getValue();
	}
	if (do_norm.getValue())
	{
		// Get an integer number for the bg_radius
		bg_radius = (int)bg_radius;
		command += " --norm --bg_radius " + floatToString(bg_radius);
		command += " --white_dust " + floatToString(white_dust.getValue());
		command += " --black_dust " + floatToString(black_dust.getValue());
	}
	if (do_invert.getValue())
		command += " --invert_contrast ";

	if (is_continue)
		command += " --only_extract_unfinished ";

	// Other arguments for extraction
	command += " " + other_args.getValue();
	commands.push_back(command);

	// Also touch the suffix file. Do this after the first command had completed
	command = "echo " + outputname + "all_micrographs.star  > " +  outputname + "coords_suffix_extract.star";
	commands.push_back(command.c_str());
	Node node3(outputname + "coords_suffix_extract.star", NODE_MIC_COORDS);
	pipelineOutputNodes.push_back(node3);

	// B. Then get the actual movie-refinement command
	if (nr_mpi.getValue() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

	// Make specific run-names for different variables of ravg, sigma_offset and sigma_angles,
	// This way you don't have to re-extract all movie particles each time you try a different parameter setting
	std::string runname = "run_ravg" + floatToString(movie_runavg_window.getValue()) + "_off" + floatToString(movie_sigma_offset.getValue());
	if (do_alsorot_movies.getValue())
		runname += "_ang" + floatToString(movie_sigma_angles.getValue());

	command += " --o " + outputname + runname;
	pipelineOutputNodes = getOutputNodesRefine(outputname + runname, -1, 1, 3, 1, true, do_alsorot_movies.getValue() );

	command += " --continue " + fn_cont.getValue();

	if (join_nr_mics.getValue() > 0)
	{
		if (do_alsorot_movies.getValue())
			REPORT_ERROR("MovieRefineJobWindow ERROR: you cannot process micrographs in batches and perform rotational searches!");

		command += " --process_movies_in_batches --realign_movie_frames " + fn_olist;

		if (is_continue)
			command += " --only_do_unfinished_movies ";
	}
	else
	{
		command += " --realign_movie_frames " + fn_ostar;
	}

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

	// Running stuff
	command += " --j " + floatToString(nr_threads.getValue());

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

PolishJobWindow::PolishJobWindow() : RelionJobWindow(5, HAS_MPI, HAS_THREAD)
{

	type = PROC_POLISH;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input STAR file with aligned movies:", NODE_MOVIE_DATA, "", "STAR files (*_data.star)",  "Provide the data.star file that was output by the movie-processing option in the auto-refine job.");

	fn_mask.place(current_y, "Mask for the reconstructions", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "A continuous mask with values between 0 (solvent) and 1 (protein). You may provide the same map that was obtained in the post-processing of the corresponding auto-refine jobs before the movie processing.");

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

	average_frame_bfactor.place(current_y, "Average frames B-factor estimation", 1 , 1, 7, 1, "B-factors for each movie frame will be estimated from reconstructions of all particles for that movie frame. Single-frame reconstructions sometimes give not enough signal to estimate reliable B-factors. \
This option allows one to calculate the B-factors from running averages of movie frames. The value specified should be an odd number. Calculating B-factors from multiple movie frames improves the SNR in the reconstructions, but limits the estimation of sudden changes in B-factors throughout the movie, for example in the first few frames when beam-induced movement is very rapid. \
Therefore, one should not use higher values than strictly necessary.");

	weight_group->end();
	do_bfactor_weighting.cb_menu_i();

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

	tab3->end();


	tab4->begin();
	tab4->label("Normalise");
	resetHeight();

	bg_diameter.place(current_y, "Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");

	white_dust.place(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	black_dust.place(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	tab4->end();

	tab5->begin();
	tab5->label("Helix");
	resetHeight();

	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	do_helix.place(current_y, "Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.", helix_group);

	helix_group->begin();

	helical_nr_asu.place(current_y, "Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the interbox distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.\n\nPlease copy this value from previous 3D refinement.");

	helical_twist.placeOnSameYPosition(current_y, "Helical twist (deg), rise (A):", "Helical twist (deg):", "0", NULL, XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	helical_rise.placeOnSameYPosition(current_y, "", "Helical rise (A):", "0", "Helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms.\n\nPlease copy the refined helical symmetry from previous 3D refinement.", XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);

	helix_group->end();
	do_helix.cb_menu_i(); // to make default effective

	tab5->end();

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
	fn_mask.writeValue(fh);
	do_fit_movement.writeValue(fh);
	sigma_nb.writeValue(fh);
	do_bfactor_weighting.writeValue(fh);
	perframe_highres.writeValue(fh);
	perframe_bfac_lowres.writeValue(fh);
	average_frame_bfactor.writeValue(fh);
	sym_name.writeValue(fh);
	bg_diameter.writeValue(fh);
	white_dust.writeValue(fh);
	black_dust.writeValue(fh);
	do_helix.writeValue(fh);
	helical_nr_asu.writeValue(fh);
	helical_twist.writeValue(fh);
	helical_rise.writeValue(fh);

	closeWriteFile(fh, fn);
}
void PolishJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_polish";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		fn_mask.readValue(fh);
		do_fit_movement.readValue(fh);
		sigma_nb.readValue(fh);
		do_bfactor_weighting.readValue(fh);
		perframe_highres.readValue(fh);
		perframe_bfac_lowres.readValue(fh);
		average_frame_bfactor.readValue(fh);
		sym_name.readValue(fh);
		bg_diameter.readValue(fh);
		white_dust.readValue(fh);
		black_dust.readValue(fh);
		do_helix.readValue(fh);
		helical_nr_asu.readValue(fh);
		helical_twist.readValue(fh);
		helical_rise.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}
void PolishJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_in.deactivate(is_continue);
	fn_mask.deactivate(is_continue);

}

bool PolishJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{
	commands.clear();
	std::string command;
	initialisePipeline(outputname, "Polish", job_counter);

	if (nr_mpi.getValue() > 1)
		command="`which relion_particle_polish_mpi`";
	else
		command="`which relion_particle_polish`";

	// General
	command += " --i " + fn_in.getValue();
	Node node(fn_in.getValue(), fn_in.type);
	pipelineInputNodes.push_back(node);

	command += " --mask " + fn_mask.getValue();
	Node node2(fn_mask.getValue(), fn_mask.type);
	pipelineInputNodes.push_back(node2);

	command += " --o " + outputname;
	Node node3(outputname + "shiny.star", NODE_PART_DATA);
	pipelineOutputNodes.push_back(node3);

	Node node4(outputname + "logfile.pdf", NODE_PDF_LOGFILE);
	pipelineOutputNodes.push_back(node4);

	Node node5(outputname + "shiny_post.mrc", NODE_FINALMAP);
	pipelineOutputNodes.push_back(node5);

	// If this is not a continue job, then re-start from scratch....
	if (is_continue)
		command += " --only_do_unfinished ";

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
		command += " --bfactor_running_avg " + floatToString(average_frame_bfactor.getValue());
	}
	else
	{
		command += " --skip_bfactor_weighting ";
	}

	// Symmetry group
	command += " --sym " + sym_name.getValue();

	// Normalisation
	command += " --bg_radius " + floatToString(FLOOR(bg_diameter.getValue()/2.));
	command += " --white_dust " + floatToString(white_dust.getValue());
	command += " --black_dust " + floatToString(black_dust.getValue());

	// Helix
	if (do_helix.getValue())
	{
		command += " --helical_nr_asu " + integerToString(helical_nr_asu.getValue());
		command += " --helical_twist " + floatToString(textToFloat(helical_twist.getValue()));
		command += " --helical_rise " + floatToString(textToFloat(helical_rise.getValue()));
	}

	// Other arguments for extraction
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}


ClassSelectJobWindow::ClassSelectJobWindow() : RelionJobWindow(2, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_CLASSSELECT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_model.place(current_y, "Select classes from model.star:", NODE_MODEL, "", "STAR files (*.star)", "A _model.star file from a previous 2D or 3D classification run to select classes from.");
	fn_mic.place(current_y, "OR select from micrographs.star:", NODE_MICS, "", "STAR files (*.star)", "A micrographs.star file to select micrographs from.");
	fn_data.place(current_y, "OR select from particles.star:", NODE_PART_DATA, "", "STAR files (*.star)", "A particles.star file to select individual particles from.");
	fn_coords.place(current_y, "OR select from picked coords:", NODE_MIC_COORDS, "", "STAR files (coords_suffix*.star)", "A coordinate suffix .star file to select micrographs while inspecting coordinates (and/or CTFs).");

	tab1->end();

	tab2->begin();
	tab2->label("Class options");
	resetHeight();

	do_recenter.place(current_y, "Re-center the class averages?", true, "This option is only used when selecting particles from 2D classes. The selected class averages will all re-centered on their center-of-mass. This is useful when you plane to use these class averages as templates for auto-picking.");
	regroup_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	regroup_group->end();
	do_regroup.place(current_y, "Regroup the particles?", false, "If set to Yes, then the program will regroup the selected particles in 'more-or-less' the number of groups indicated below. For re-grouping from individual particle _data.star files, a _model.star file with the same prefix should exist, i.e. the particle star file should be generated by relion_refine", regroup_group);
	regroup_group->begin();
	nr_groups.place(current_y, "Approximate nr of groups: ", 1, 50, 20, 1, "It is normal that the actual number of groups may deviate a little from this number. ");
	regroup_group->end();
	do_regroup.cb_menu_i(); // make default active

	tab2->end();

	// read settings if hidden file exists
	read(".gui_select", is_continue);
}



void ClassSelectJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_select";

	std::ofstream fh;
	openWriteFile(fn, fh);

	fn_model.writeValue(fh);
	fn_data.writeValue(fh);
	fn_mic.writeValue(fh);
	fn_coords.writeValue(fh);
	do_recenter.writeValue(fh);
	do_regroup.writeValue(fh);
	nr_groups.writeValue(fh);

	closeWriteFile(fh, fn);
}

void ClassSelectJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_select";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_model.readValue(fh);
		fn_data.readValue(fh);
		fn_mic.readValue(fh);
		fn_coords.readValue(fh);
		do_recenter.readValue(fh);
		do_regroup.readValue(fh);
		nr_groups.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;

		// For re-setting of new jobs
		ori_fn_model = fn_model.getValue();
		ori_fn_data = fn_data.getValue();
		ori_fn_mic = fn_mic.getValue();
		ori_fn_coords = fn_coords.getValue();

	}
}


void ClassSelectJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_model.deactivate(is_continue);
	fn_data.deactivate(is_continue);
	fn_mic.deactivate(is_continue);
	fn_coords.deactivate(is_continue);

	// For new jobs, always reset the input fields to empty
	if (!_is_continue)
	{
		fn_model.setValue("");
		fn_data.setValue("");
		fn_mic.setValue("");
		fn_coords.setValue("");
	}
	else
	{
		fn_model.setValue(ori_fn_model.c_str());
		fn_data.setValue(ori_fn_data.c_str());
		fn_mic.setValue(ori_fn_mic.c_str());
		fn_coords.setValue(ori_fn_coords.c_str());
	}

	do_recenter.deactivate(is_continue);
	do_regroup.deactivate(is_continue);
	nr_groups.deactivate(is_continue);

	do_queue.deactivate(true);

}

bool ClassSelectJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Select", job_counter);

	std::string command;
	command="`which relion_display`";

	// I/O
	if (fn_model.getValue() != "")
	{
		command += " --gui --i " + fn_model.getValue();
		Node node(fn_model.getValue(), fn_model.type);
		pipelineInputNodes.push_back(node);

		FileName fn_parts = outputname+"particles.star";
		command += " --allow_save --fn_parts " + fn_parts;
		Node node2(fn_parts, NODE_PART_DATA);
		pipelineOutputNodes.push_back(node2);

		// Only save the 2D class averages for 2D jobs
		FileName fnt = fn_model.getValue();
		if (fnt.contains("Class2D/"))
		{
			FileName fn_imgs = outputname+"class_averages.star";
			command += " --fn_imgs " + fn_imgs;
			Node node3(fn_imgs, NODE_2DREFS);
			pipelineOutputNodes.push_back(node3);

			if (do_recenter.getValue())
			{
				command += " --recenter";
			}
		}
	}
	else if (fn_mic.getValue() != "")
	{
		command += " --gui --i " + fn_mic.getValue();
		Node node(fn_mic.getValue(), fn_mic.type);
		pipelineInputNodes.push_back(node);

		FileName fn_mics = outputname+"micrographs.star";
		command += " --allow_save --fn_imgs " + fn_mics;
		Node node2(fn_mics, NODE_MICS);
		pipelineOutputNodes.push_back(node2);
	}
	else if (fn_data.getValue() != "")
	{
		command += " --gui --i " + fn_data.getValue();
		Node node(fn_data.getValue(), fn_data.type);
		pipelineInputNodes.push_back(node);

		FileName fn_parts = outputname+"particles.star";
		command += " --allow_save --fn_imgs " + fn_parts;
		Node node2(fn_parts, NODE_PART_DATA);
		pipelineOutputNodes.push_back(node2);
	}
	else if  (fn_coords.getValue() != "")
	{

    	ManualpickJobWindow global_manualpickjob;

    	FileName fn_job = ".gui_manualpickrun.job";
		bool iscont=false;
		if (exists(fn_job))
			global_manualpickjob.read(fn_job.c_str(), iscont);
		else
			REPORT_ERROR("RelionMainWindow::cb_display_io_node_i ERROR: Save a Manual picking job parameters (using the File menu) before displaying coordinate files. ");

		// Get the name of the micrograph STAR file from reading the suffix file
	    FileName fn_suffix = fn_coords.getValue();
		FileName fn_star;
	    if (is_continue)
	    {
	    	fn_star = outputname + "micrographs_selected.star";
	    }
	    else
	    {
	    	std::ifstream in(fn_suffix.data(), std::ios_base::in);
	    	in >> fn_star ;
	    	in.close();
	    }
	    FileName fn_dirs = fn_suffix.beforeLastOf("/")+"/";
		fn_suffix = fn_suffix.afterLastOf("/").without("coords_suffix_");
		fn_suffix = fn_suffix.withoutExtension();

		// Launch the manualpicker...
		command="`which relion_manualpick` --i " + fn_star;
		Node node4(fn_coords.getValue(), fn_coords.type);
		pipelineInputNodes.push_back(node4);

		command += " --odir " + fn_dirs;
		command += " --pickname " + fn_suffix;

		// The output selection
		FileName fn_outstar = outputname + "micrographs_selected.star";
		Node node3(fn_outstar, NODE_MICS);
		pipelineOutputNodes.push_back(node3);
		command += " --selection " + fn_outstar;
		command += " --allow_save  --selection " + fn_outstar;

		// All the stuff from the saved global_manualpickjob
		command += " --scale " + floatToString(global_manualpickjob.micscale.getValue());
		command += " --sigma_contrast " + floatToString(global_manualpickjob.sigma_contrast.getValue());
		command += " --black " + floatToString(global_manualpickjob.black_val.getValue());
		command += " --white " + floatToString(global_manualpickjob.white_val.getValue());

		if (global_manualpickjob.lowpass.getValue() > 0.)
			command += " --lowpass " + floatToString(global_manualpickjob.lowpass.getValue());
		if (global_manualpickjob.highpass.getValue() > 0.)
			command += " --highpass " + floatToString(global_manualpickjob.highpass.getValue());
		if (global_manualpickjob.angpix.getValue() > 0.)
			command += " --angpix " + floatToString(global_manualpickjob.angpix.getValue());

		command += " --ctf_scale " + floatToString(global_manualpickjob.ctfscale.getValue());

		command += " --particle_diameter " + floatToString(global_manualpickjob.diameter.getValue());

		if (global_manualpickjob.do_color.getValue())
		{
			command += " --color_label " + global_manualpickjob.color_label.getValue();
			command += " --blue " + floatToString(global_manualpickjob.blue_value.getValue());
			command += " --red " + floatToString(global_manualpickjob.red_value.getValue());
			if (global_manualpickjob.fn_color.getValue().length() > 0)
				command += " --color_star " + global_manualpickjob.fn_color.getValue();
		}

		// Other arguments for extraction
		command += " " + global_manualpickjob.other_args.getValue() + " &";

	}

	// Re-grouping
	if (do_regroup.getValue() && fn_coords.getValue() == "")
	{
		command += " --regroup " + floatToString(nr_groups.getValue());
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	// For picked coordinates: also write the selected STAR file in the coord_suffix file
	//if  (fn_coords.getValue() != "" &&!is_continue)
	//{
	//	std::string command2 = "echo " + outputname + "micrographs_selected.star > " + fn_coords.getValue();
	//	commands.push_back(command2);
	//}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}


MaskCreateJobWindow::MaskCreateJobWindow() : RelionJobWindow(3, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_MASKCREATE;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "Input 3D map:", NODE_3DREF, "", "MRC map files (*.mrc)", "Provide an input MRC map from which to start binarizing the map.");
	tab1->end();

	tab2->begin();
	tab2->label("Mask");
	resetHeight();

	lowpass_filter.place(current_y, "Lowpass filter map (A)", -1, 10, 100, 5, "Lowpass filter that will be applied to the input map, prior to binarization. To calculate solvent masks, a lowpass filter of 15-20A may work well.");
	angpix.place(current_y, "Pixel size (A)", 1, 0.3, 5, 0.1, "Provide the pixel size in Angstroms of the input map (to calculate the low-pass filter)");

	// Add a little spacer
    current_y += STEPY/2;

	inimask_threshold.place(current_y, "Initial binarisation threshold:", 0.02, 0., 0.5, 0.01, "This threshold is used to make an initial binary mask from the average of the two unfiltered half-reconstructions. \
If you don't know what value to use, display one of the unfiltered half-maps in a 3D surface rendering viewer and find the lowest threshold that gives no noise peaks outside the reconstruction.");
	extend_inimask.place(current_y, "Extend binary map this many pixels:", 3, 0, 20, 1, "The initial binary mask is extended this number of pixels in all directions." );
	width_mask_edge.place(current_y, "Add a soft-edge of this many pixels:", 3, 0, 20, 1, "The extended binary mask is further extended with a raised-cosine soft edge of the specified width." );

	tab2->end();

	tab3->begin();
	tab3->label("Helix");
	resetHeight();

	helix_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	helix_group->end();

	do_helix.place(current_y, "Mask a 3D helix?", false, "Generate a mask for 3D helix which spans across Z axis of the box.", helix_group);

	helix_group->begin();

	helical_z_percentage.place(current_y, "Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used but you may want to try different lengths.");
	helix_group->end();
	do_helix.cb_menu_i(); // to make default effective

	tab3->end();

	// read settings if hidden file exists
	read(".gui_maskcreate", is_continue);
}



void MaskCreateJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_maskcreate";

	std::ofstream fh;
	openWriteFile(fn, fh);

	fn_in.writeValue(fh);
	lowpass_filter.writeValue(fh);
	angpix.writeValue(fh);
	inimask_threshold.writeValue(fh);
	extend_inimask.writeValue(fh);
	width_mask_edge.writeValue(fh);
	do_helix.writeValue(fh);
	helical_z_percentage.writeValue(fh);

	closeWriteFile(fh, fn);
}

void MaskCreateJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_maskcreate";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		lowpass_filter.readValue(fh);
		angpix.readValue(fh);
		inimask_threshold.readValue(fh);
		extend_inimask.readValue(fh);
		width_mask_edge.readValue(fh);
		do_helix.readValue(fh);
		helical_z_percentage.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;
	}
}


void MaskCreateJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	fn_in.deactivate(is_continue);
}

bool MaskCreateJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "MaskCreate", job_counter);

	std::string command;
	command="`which relion_mask_create`";

	// I/O
	command += " --i " + fn_in.getValue();
	Node node(fn_in.getValue(), fn_in.type);
	pipelineInputNodes.push_back(node);

	command += " --o " + outputname + "mask.mrc";
	Node node2(outputname + "mask.mrc", NODE_MASK);
	pipelineOutputNodes.push_back(node2);

	// TODO: implement low-pass filter!
	if (lowpass_filter.getValue() > 0)
	{
		command += " --lowpass " + floatToString(lowpass_filter.getValue());
		command += " --angpix " + floatToString(angpix.getValue());
	}
	command += " --ini_threshold " + floatToString(inimask_threshold.getValue());
	command += " --extend_inimask " + floatToString(extend_inimask.getValue());
	command += " --width_soft_edge " + floatToString(width_mask_edge.getValue());

	if (do_helix.getValue())
	{
		command += " --helix --z_percentage " + floatToString(helical_z_percentage.getValue() / 100.);
	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}

JoinStarJobWindow::JoinStarJobWindow() : RelionJobWindow(1, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_JOINSTAR;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	part_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	part_group->end();
	do_part.place(current_y, "Combine particle STAR files?", true, "", part_group);
	part_group->begin();
	fn_part1.place(current_y, "Particle STAR file 1: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The first of the particle STAR files to be combined.");
	fn_part2.place(current_y, "Particle STAR file 2: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The second of the particle STAR files to be combined.");
	fn_part3.place(current_y, "Particle STAR file 3: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The third of the particle STAR files to be combined. Leave empty if there are only two files to be combined.");
	fn_part4.place(current_y, "Particle STAR file 4: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The fourth of the particle STAR files to be combined. Leave empty if there are only two or three files to be combined.");
	part_group->end();
	do_part.cb_menu_i(); // make default active

	// Add a little spacer
    current_y += STEPY/2;

	mic_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	mic_group->end();
	do_mic.place(current_y, "Combine micrograph STAR files?", false, "", mic_group);
	mic_group->begin();
	fn_mic1.place(current_y, "Micrograph STAR file 1: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The first of the micrograph STAR files to be combined.");
	fn_mic2.place(current_y, "Micrograph STAR file 2: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The second of the micrograph STAR files to be combined.");
	fn_mic3.place(current_y, "Micrograph STAR file 3: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The third of the micrograph STAR files to be combined. Leave empty if there are only two files to be combined.");
	fn_mic4.place(current_y, "Micrograph STAR file 4: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The fourth of the micrograph STAR files to be combined. Leave empty if there are only two or three files to be combined.");
	mic_group->end();
	do_mic.cb_menu_i(); // make default active
	tab1->end();

	// read settings if hidden file exists
	read(".gui_joinstar", is_continue);
}



void JoinStarJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_joinstar";

	std::ofstream fh;
	openWriteFile(fn, fh);

	do_part.writeValue(fh);
	fn_part1.writeValue(fh);
	fn_part2.writeValue(fh);
	fn_part3.writeValue(fh);
	fn_part4.writeValue(fh);

	do_mic.writeValue(fh);
	fn_mic1.writeValue(fh);
	fn_mic2.writeValue(fh);
	fn_mic3.writeValue(fh);
	fn_mic4.writeValue(fh);

	closeWriteFile(fh, fn);
}

void JoinStarJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_joinstar";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		do_part.readValue(fh);
		fn_part1.readValue(fh);
		fn_part2.readValue(fh);
		fn_part3.readValue(fh);
		fn_part4.readValue(fh);

		do_mic.readValue(fh);
		fn_mic1.readValue(fh);
		fn_mic2.readValue(fh);
		fn_mic3.readValue(fh);
		fn_mic4.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;

		// For re-setting of new jobs
		ori_fn_part1 = fn_part1.getValue();
		ori_fn_part2 = fn_part2.getValue();
		ori_fn_part3 = fn_part3.getValue();
		ori_fn_part4 = fn_part4.getValue();
		ori_fn_mic1 = fn_mic1.getValue();
		ori_fn_mic2 = fn_mic2.getValue();
		ori_fn_mic3 = fn_mic3.getValue();
		ori_fn_mic4 = fn_mic4.getValue();

	}
}


void JoinStarJobWindow::toggle_new_continue(bool _is_continue)
{
	is_continue = _is_continue;

	// For new jobs, always reset the input fields to empty
	if (!_is_continue)
	{
		fn_part1.setValue("");
		fn_part2.setValue("");
		fn_part3.setValue("");
		fn_part4.setValue("");
		fn_mic1.setValue("");
		fn_mic2.setValue("");
		fn_mic3.setValue("");
		fn_mic4.setValue("");
	}
	else
	{
		fn_part1.setValue(ori_fn_part1.c_str());
		fn_part2.setValue(ori_fn_part2.c_str());
		fn_part3.setValue(ori_fn_part3.c_str());
		fn_part4.setValue(ori_fn_part4.c_str());
		fn_mic1.setValue(ori_fn_mic1.c_str());
		fn_mic2.setValue(ori_fn_mic2.c_str());
		fn_mic3.setValue(ori_fn_mic3.c_str());
		fn_mic4.setValue(ori_fn_mic4.c_str());
	}

	do_part.deactivate(is_continue);
	fn_part1.deactivate(is_continue);
	fn_part2.deactivate(is_continue);
	fn_part3.deactivate(is_continue);
	fn_part4.deactivate(is_continue);

	do_mic.deactivate(is_continue);
	fn_mic1.deactivate(is_continue);
	fn_mic2.deactivate(is_continue);
	fn_mic3.deactivate(is_continue);
	fn_mic4.deactivate(is_continue);
}

bool JoinStarJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "JoinStar", job_counter);

	std::string command;
	command="`which relion_star_combine`";

	// I/O  TODO!!!
	if (do_part.getValue())
	{
		command += " --i \" " + fn_part1.getValue();
		Node node(fn_part1.getValue(), fn_part1.type);
		pipelineInputNodes.push_back(node);
		command += " " + fn_part2.getValue();
		Node node2(fn_part2.getValue(), fn_part2.type);
		pipelineInputNodes.push_back(node2);
		if (fn_part3.getValue() != "")
		{
			command += " " + fn_part3.getValue();
			Node node3(fn_part3.getValue(), fn_part3.type);
			pipelineInputNodes.push_back(node3);
		}
		if (fn_part4.getValue() != "")
		{
			command += " " + fn_part4.getValue();
			Node node4(fn_part4.getValue(), fn_part4.type);
			pipelineInputNodes.push_back(node4);
		}
		command += " \" ";

		// Check for duplicates
		command += " --check_duplicates rlnImageName ";
		command += " --o " + outputname + "join_particles.star";
		Node node5(outputname + "join_particles.star", fn_part1.type);
		pipelineOutputNodes.push_back(node5);

	}
	else if (do_mic.getValue())
	{
		command += " --i \" " + fn_mic1.getValue();
		Node node(fn_mic1.getValue(), fn_mic1.type);
		pipelineInputNodes.push_back(node);
		command += " " + fn_mic2.getValue();
		Node node2(fn_mic2.getValue(), fn_mic2.type);
		pipelineInputNodes.push_back(node2);
		if (fn_mic3.getValue() != "")
		{
			command += " " + fn_mic3.getValue();
			Node node3(fn_mic3.getValue(), fn_mic3.type);
			pipelineInputNodes.push_back(node3);
		}
		if (fn_mic4.getValue() != "")
		{
			command += " " + fn_mic4.getValue();
			Node node4(fn_mic4.getValue(), fn_mic4.type);
			pipelineInputNodes.push_back(node4);
		}
		command += " \" ";

		// Check for duplicates
		command += " --check_duplicates rlnMicrographName ";
		command += " --o " + outputname + "join_mics.star";
		Node node5(outputname + "join_mics.star", fn_mic1.type);
		pipelineOutputNodes.push_back(node5);

	}

	// Other arguments
	command += " " + other_args.getValue();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}



SubtractJobWindow::SubtractJobWindow() : RelionJobWindow(2, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_SUBTRACT;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_data.place(current_y, "Input particles:", NODE_PART_DATA, "", "STAR files (*.star)", "The input STAR file with the metadata of all particles.");

    current_y += STEPY/2;
	subtract_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	subtract_group->end();

	do_subtract.place(current_y, "Subtract partial signal?", true, "If set to Yes, a copy of the entire set of particle images in this STAR file will be made that contains the subtracted particle images.", subtract_group);

	subtract_group->begin();
	fn_in.place(current_y, "Map to be projected:", NODE_3DREF, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide the map that will be used to calculate projections, which will be subtracted from the experimental particles. Make sure this map was calculated by RELION from the same data set, as it is crucial that the absolute greyscale is the same as in the experimental particles.");
    fn_mask.place(current_y, "Mask to apply to this map:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein density you wish to subtract from the experimental particles is white (1) and the rest of the protein and the solvent is black (0).");
    subtract_group->end();
    do_subtract.cb_menu_i(); // make default active

    do_fliplabel.place(current_y, "OR revert to original particles?", false, "If set to Yes, no signal subtraction is performed. Instead, the labels of rlnImageName and rlnImageOriginalName are flipped in the input STAR file. This will make the STAR file point back to the original, non-subtracted images.");

	tab1->end();

	tab2->begin();
	tab2->label("CTF");
	resetHeight();

	ctf_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	ctf_group->end();

	do_ctf_correction.place(current_y, "Do apply CTFs?", true, "If set to Yes, CTFs will be applied to the projections to subtract from the experimental particles", ctf_group);

	ctf_group->begin();

	ctf_phase_flipped.place(current_y, "Have data been phase-flipped?", false, "Set this to Yes if the experimental particles have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");

	ctf_intact_first_peak.place(current_y, "Ignore CTFs until first peak?", false, "Set this to Yes when you have used the same option in the generation of the input map to be projected.");

	ctf_group->end();
	do_ctf_correction.cb_menu_i(); // To make default effective

	tab2->end();

	// read settings if hidden file exists
	read(".gui_subtract", is_continue);
}



void SubtractJobWindow::write(std::string fn)
{
	// Write hidden file if no name is given
	if (fn=="")
		fn=".gui_subtract";

	std::ofstream fh;
	openWriteFile(fn, fh);

	fn_data.writeValue(fh);
	fn_in.writeValue(fh);
	fn_mask.writeValue(fh);
	do_subtract.writeValue(fh);
	do_fliplabel.writeValue(fh);

	do_ctf_correction.writeValue(fh);
	ctf_phase_flipped.writeValue(fh);
	ctf_intact_first_peak.writeValue(fh);

	closeWriteFile(fh, fn);
}

void SubtractJobWindow::read(std::string fn, bool &_is_continue)
{

	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_subtract";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{

		fn_data.readValue(fh);
		fn_in.readValue(fh);
		fn_mask.readValue(fh);
		do_subtract.readValue(fh);
		do_fliplabel.readValue(fh);

		do_ctf_correction.readValue(fh);
		ctf_phase_flipped.readValue(fh);
		ctf_intact_first_peak.readValue(fh);

		closeReadFile(fh);
		_is_continue = is_continue;

	}
}


void SubtractJobWindow::toggle_new_continue(bool _is_continue)
{
	fn_data.deactivate(is_continue);
	fn_in.deactivate(is_continue);
	fn_mask.deactivate(is_continue);
	do_subtract.deactivate(is_continue);
	do_fliplabel.deactivate(is_continue);

	do_ctf_correction.deactivate(is_continue);
	ctf_phase_flipped.deactivate(is_continue);
	ctf_intact_first_peak.deactivate(is_continue);

	is_continue = _is_continue;
}

bool SubtractJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Subtract", job_counter);

	std::string command;
	if (do_subtract.getValue())
	{
		command="`which relion_project`";

		// I/O
		command += " --subtract_exp --i " + fn_in.getValue();
		Node node(fn_in.getValue(), fn_in.type);
		pipelineInputNodes.push_back(node);
		command += " --mask " + fn_mask.getValue();
		Node node2(fn_mask.getValue(), fn_mask.type);
		pipelineInputNodes.push_back(node2);
		command += " --ang " + fn_data.getValue();
		Node node3(fn_data.getValue(), fn_data.type);
		pipelineInputNodes.push_back(node3);

		command += " --o " + outputname + "subtracted";
		Node node4(outputname + "subtracted.star", NODE_PART_DATA);
		pipelineOutputNodes.push_back(node4);

		if (do_ctf_correction.getValue())
		{
			command += " --ctf --angpix -1";
			if (ctf_phase_flipped.getValue())
				command += " --ctf_phase_flip";
			if (ctf_intact_first_peak.getValue())
				command += " --ctf_intact_first_peak";
		}

		// Other arguments
		command += " " + other_args.getValue();
	}
	else if (do_fliplabel.getValue())
	{
		Node node(fn_data.getValue(), fn_data.type);
		pipelineInputNodes.push_back(node);

		Node node2(outputname + "original.star", NODE_PART_DATA);
		pipelineOutputNodes.push_back(node2);

		command = "awk '{if  ($1==\"_rlnImageName\") {$1=\"_rlnImageOriginalName\"} else if ($1==\"_rlnImageOriginalName\") {$1=\"_rlnImageName\"}; print }' < ";
		command += fn_data.getValue() + " > " + outputname + "original.star";
	}

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);

}





PostJobWindow::PostJobWindow() : RelionJobWindow(3, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_POST;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();
	fn_in.place(current_y, "One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*half1_class001_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	fn_mask.place(current_y, "Solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");

	current_y += STEPY/2;

	angpix.place(current_y, "Calibrated pixel size (A)", 1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	tab1->end();

	tab2->begin();
	tab2->label("Sharpen");
	resetHeight();

	fn_mtf.place(current_y, "MTF of the detector (STAR file)", "", "STAR Files (*.star)", NULL, "The MTF of the detector is used in the (later) post-processing and particle polishing stages of refinement.  \
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

	tab2->end();
	tab3->begin();
	tab3->label("Filter");
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

	tab3->end();


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
	fn_mask.writeValue(fh);
	angpix.writeValue(fh);
	do_auto_bfac.writeValue(fh);
	autob_lowres.writeValue(fh);
	do_adhoc_bfac.writeValue(fh);
	adhoc_bfac.writeValue(fh);
	fn_mtf.writeValue(fh);
	do_skip_fsc_weighting.writeValue(fh);
	low_pass.writeValue(fh);
	closeWriteFile(fh, fn);
}

void PostJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_post";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_in.readValue(fh);
		fn_mask.readValue(fh);
		angpix.readValue(fh);
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

	fn_in.deactivate(is_continue);
	fn_mask.deactivate(is_continue);

}

bool PostJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{


	commands.clear();
	initialisePipeline(outputname, "PostProcess", job_counter);
	std::string command;

	command="`which relion_postprocess`";

	// Input mask
	command += " --mask " + fn_mask.getValue();
	Node node3(fn_mask.getValue(), fn_mask.type);
	pipelineInputNodes.push_back(node3);

	// Get the input rootname from the half-map name
	// run1_half1_class001_unfil.mrc -> run1
	Node node(fn_in.getValue(), fn_in.type);
	pipelineInputNodes.push_back(node);
	int pos_half = fn_in.getValue().rfind("_half");
	if (pos_half < fn_in.getValue().size())
		command += " --i " + fn_in.getValue().substr(0, pos_half);
	else
	{
		std::cerr << "PostJobWindow::getCommands ERROR: cannot find _half substring in input filename: " << fn_in.getValue() << std::endl;
		exit(1);
	}

	// The output name contains a directory: use it for output
	command += " --o " + outputname + "postprocess";
	command += "  --angpix " + floatToString(angpix.getValue());
	Node node1(outputname+"postprocess.mrc", NODE_FINALMAP);
	pipelineOutputNodes.push_back(node1);
	Node node2(outputname+"postprocess_masked.mrc", NODE_FINALMAP);
	pipelineOutputNodes.push_back(node2);

	Node node2b(outputname+"logfile.pdf", NODE_PDF_LOGFILE);
	pipelineOutputNodes.push_back(node2b);


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
	return prepareFinalCommand(outputname, commands, final_command, do_makedir);
}


ResmapJobWindow::ResmapJobWindow() : RelionJobWindow(2, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = PROC_RESMAP;

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	fn_in.place(current_y, "One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");

	current_y += STEPY /2 ;

	fn_mask.place(current_y, "User-provided solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a mask with values between 0 and 1 around all domains of the complex.");

	current_y += STEPY/2;

	angpix.place(current_y, "Calibrated pixel size (A)", 1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	tab1->end();

	tab2->begin();
	tab2->label("ResMap");
	resetHeight();

	// Check for environment variable RELION_RESMAP_TEMPLATE
	char * default_location = getenv ("RELION_RESMAP_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[] = DEFAULTRESMAPLOCATION;
		default_location = mydefault;
	}

	fn_resmap.place(current_y, "ResMap executable:", default_location, "ResMap*", NULL, "Location of the ResMap executable. You can control the default of this field by setting environment variable RELION_RESMAP_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");

	current_y += STEPY /2 ;

	pval.place(current_y, "P-value:", 0.05, 0., 1., 0.01, "This value is typically left at 0.05. If you change it, report the modified value in your paper!");
	minres.place(current_y, "Highest resolution (A): ", 0., 0., 10., 0.1, "ResMaps minRes parameter. By default (0), the program will start at just above 2x the pixel size");
	maxres.place(current_y, "Lowest resolution (A): ", 0., 0., 10., 0.1, "ResMaps maxRes parameter. By default (0), the program will stop at 4x the pixel size");
	stepres.place(current_y, "Resolution step size (A)", 1., 0.1, 3, 0.1, "ResMaps stepSize parameter." );

	tab2->end();

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
	angpix.writeValue(fh);
	pval.writeValue(fh);
	minres.writeValue(fh);
	maxres.writeValue(fh);
	stepres.writeValue(fh);
	fn_mask.writeValue(fh);
	closeWriteFile(fh, fn);
}

void ResmapJobWindow::read(std::string fn, bool &_is_continue)
{
	std::ifstream fh;

	// Read hidden file if no name is given
	if (fn=="")
		fn=".gui_resmap";

	// Only read things if the file exists
	if (openReadFile(fn, fh))
	{
		fn_resmap.readValue(fh);
		fn_in.readValue(fh);
		angpix.readValue(fh);
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

	fn_resmap.deactivate(is_continue);
	fn_in.deactivate(is_continue);
	angpix.deactivate(is_continue);
	pval.deactivate(is_continue);
	minres.deactivate(is_continue);
	maxres.deactivate(is_continue);
	stepres.deactivate(is_continue);
	fn_mask.deactivate(is_continue);

	// never submit this to queue, as ResMap needs user interaction
	do_queue.deactivate(true);
}

bool ResmapJobWindow::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter)
{

	commands.clear();
	initialisePipeline(outputname, "Resmap", job_counter);

	if (fn_resmap.getValue().length() == 0)
	{
		std::cerr << "ERROR: please provide an executable for the ResMap program." << std::endl;
		exit(1);
	}

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
		std::cerr << "ERROR: cannot find _half substring in input filename: " << fn_in.getValue() << std::endl;
		exit(1);
	}

	// Make symbolic links to the half-maps in the output directory
	commands.push_back("ln -s ../../" + fn_half1 + " " + outputname + "half1.mrc");
	commands.push_back("ln -s ../../" + fn_half2 + " " + outputname + "half2.mrc");

	Node node(fn_in.getValue(), fn_in.type);
	pipelineInputNodes.push_back(node);

	Node node2(fn_mask.getValue(), fn_mask.type);
	pipelineInputNodes.push_back(node2);

	Node node3(outputname + "half1_resmap.mrc", NODE_RESMAP);
	pipelineOutputNodes.push_back(node3);

	std::string command = fn_resmap.getValue();
	command += " --maskVol=" + fn_mask.getValue();
	command += " --noguiSplit " + outputname + "half1.mrc " +  outputname + "half2.mrc";
	command += " --vxSize=" + floatToString(angpix.getValue());
	command += " --pVal=" + floatToString(pval.getValue());
	command += " --minRes=" + floatToString(minres.getValue());
	command += " --maxRes=" + floatToString(maxres.getValue());
	command += " --stepRes=" + floatToString(stepres.getValue());

	// Other arguments for extraction
	command += " " + other_args.getValue();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir);
}

/*
PublishJobWindow::PublishJobWindow() : RelionJobWindow(2, HAS_NOT_MPI, HAS_NOT_THREAD)
{

	type = -1;

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

*/

