#!/usr/bin/env python2.7
"""
relion_it.py
============

Simple GUI to set up RELION-3.2 scheduler.

Authors: Sjors H.W. Scheres & Colin M. Palmer

Usage:
    relion_it.py  [extra_options.py [extra_options2.py ....] ] [--nogui] 

To get started, go to the intended location of your RELION project directory and make sure your micrographs are
accessible from within it (e.g. in a subdirectory called `Movies/' - use a symlink if necessary). Then run this
script, providing the names of files containing options if needed. (To call the script, you'll need to enter the full
path to it, put the directory containing it on your PATH environment variable, or put a copy of the script in the
current directory.)

Options files
-------------

relion_it.py uses a large number of options for controlling both the flow of the script and the parameters for
individual jobs. These options can be read from Python script files when relion_it.py is started.

The options files (extra_options.py, extra_options2.py etc) are read and interpreted as Python dictionaries. A simple list of 
"{ 'option_name1' : value1, 'option_name2' : value2, ...}" is all that is needed. 
When this script is executed, all current options are saved in an option file too.

The options are named descriptively so you can probably understand what most of them do quite easily. For more help on
any particular option, look at the comment above its definition in this script, or search the script's code to see
how it is used.

Use double underscores to separate SCHEDULENAME__JOBNAME_JOBOPTION
 E.g. for SCHEDULENAME=preprocess, JOBNAME=importmovies and JOBOPTION=angpix
'preprocess__importmovies__angpix' defines the value for the 'angpix' option in the file 'Schedules/preprocess/importmovies/job.star'

Likewise, use SCHEDULENAME__VARIABLENAME for Schedule variables
 E.g. 'refine__logbatch_size' refers to the variablename 'logbatch_size' in the file 'Schedules/refine/schedule.star'

This python script will set the corresponding values of these joboption and schedule-variable values in the Schedule STAR files. 
The RELION scheduler can then be used to run the scheduler.

Using the GUI
-------------

The GUI provides a simple way to start new projects with relion_it.py. 

The window that appears should be self-explanatory. Fill in the options as needed for your project, and use the check
boxes on the right to control what processing steps will be done. When you're ready, click either "Save options" or
"Save & run". The program will check the values you've entered and then use them to calculate a few extra options for
relion_it.py. The options will then be saved to a file called `relion_it_options.py', and if you clicked "Save & run"
the processing run will start immediately.

If any of the entered values are invalid (for example, if there are letters in a field which should be a number), the
GUI will display a message box with an error when you click one of the buttons. It will also display a warning if any
values appear to be incorrect (but you can choose to ignore the warning by clicking "OK").

Using this script in a non-interactive manner
-------------

When the --nogui option is provided, the options from this script, or any addition option files are used to modify the STAR files 
and the schedulers are launched from this python script. This will allow a non-interactive use of this script.



"""

from __future__ import print_function
from __future__ import division  # always use float division

import argparse
import glob
import inspect
import math
import os
import runpy
import time
import traceback
import ast 
import collections
from shutil import copytree

try:
    import Tkinter as tk
    import tkMessageBox
    import tkFileDialog
except ImportError:
    # The GUI is optional. If the user requests it, it will fail when it tries
    # to open so we can ignore the error for now.
    pass  

OPTIONS_FILE = 'relion_it_options.py'

RelionItOptions = {
    #############################################################################
    #
    # Change the parameters below to reflect your experiment
    # Use double underscores to separate SCHEDULENAME__JOBNAME__JOBOPTION
    #
    # E.g. for SCHEDULENAME=preprocess, JOBNAME=importmovies and JOBOPTION=angpix
    #  'preprocess__importmovies__angpix' defines the value for the 'angpix' option in the file 'Schedules/preprocess/importmovies/job.star'
    #
    # Likewise, use SCHEDULENAME__VARIABLENAME for Schedule variables
    #
    # E.g. 'refine__logbatch_size' refers to the variablename 'logbatch_size' in the file 'Schedules/refine/schedule.star'
    #
    #  This python script will modify the joboption and schedule-variable values in these files
    #
    #############################################################################


    ### General parameters
    # Pixel size in Angstroms in the input movies
    'preprocess__importmovies__angpix' : 1.1,
    # Acceleration voltage (in kV)
    'preprocess__importmovies__kV' : 300,
    # Polara = 2.0; Talos/Krios = 2.7; some Cryo-ARM = 1.4 
    'preprocess__importmovies__Cs' : 2.7,

    ### Import images (Linux wild card; movies as *.mrc, *.mrcs, *.tiff or *.tif; single-frame micrographs as *.mrc)
    'preprocess__importmovies__fn_in_raw' : 'Movies/*.tiff',
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    'preprocess__importmovies__is_multiframe' : True,

    ### MotionCorrection parameters
    # Dose in electrons per squared Angstrom per fraction
    'preprocess__motioncorr__dose_per_frame' : 1.2,
    # Gain-reference image in MRC format (only necessary if input movies are not yet gain-corrected, e.g. compressed TIFFs from K2)
    'preprocess__motioncorr__fn_gain_ref' : 'Movies/gain.mrc',
    # EER fractionation. The dose rate (motioncor_doseperframe) is e/A2/fraction after this fractionation.
    'preprocess__motioncorr__eer_grouping' : 20,
    # For super-resolution movies, bin by factor of 2 straight away, as we don't believe in super-resolution
    'preprocess__motioncorr__bin_factor' : 1,
    # Number of threads to use for Relioncor2
    'preprocess__motioncorr__nr_threads' : 12,

    ### CTF estimation parameters
    # Most cases won't need changes here...
    'preprocess__ctffind__do_phaseshift' : False,
    'preprocess__do_until_ctf' : False,

    ### Autopick parameters
    # Minimum and maximum diameter in Angstrom for the LoG filter
    'preprocess__logpicker__log_diam_min' : 150,
    'preprocess__logpicker__log_diam_max' : 180,
    # Use positive values (0-1) to pick fewer particles; use negative values (-1-0) to pick more particles
    'preprocess__logpicker__log_adjust_thr' : 0.0,

    #
    # For both LoG and refs:
    #
    # Use this to remove false positives from carbon edges (useful range: 1.0-1.2, -1 to switch off)
    'preprocess__logpicker__maxstddevnoise_autopick' : -1,
    # Use this to remove false positives from carbon edges (useful range: -0.5-0.0; -999 to switch off)
    'preprocess__logpicker__minavgnoise_autopick' : -999,


    ### Extract parameters for Logpicker job
    # Box size of particles in the averaged micrographs (in pixels)
    'preprocess__extract_logpick__extract_size' : 256,
    # Down-scale the particles upon extraction?
    'preprocess__extract_logpick__do_rescale' : True,
    # Box size of the down-scaled particles (in pixels)
    'preprocess__extract_logpick__rescale' : 64,
    ### Split parameters after logpick (this will be the maximum number of particles in the first batch)
    'preprocess__split_logpick__split_size' : 5000,

    ### Extract parameters for Topaz job
    # Box size of particles in the averaged micrographs (in pixels)
    'preprocess__extract_topazpick__extract_size' : 256,
    # Down-scale the particles upon extraction?
    'preprocess__extract_topazpick__do_rescale' : True,
    # Box size of the down-scaled particles (in pixels)
    'preprocess__extract_topazpick__rescale' : 64,
    # Minimum FOM for topaz extraction
    'preprocess__extract_topazpick__minimum_pick_fom' : -3,
    
    ### Parameters for Topaz picking
    # Expected number of particles per micrograph
    'refine__train_topaz__topaz_nr_particles' : 300,
    # Which (single) GPU to run training on 
    'refine__train_topaz__gpu_ids' : 0,
    # How many MPI processes to use in parallel for picking
    'preprocess__topazpicker__nr_mpi' : 2,

    ### Parameters for automated 2D class selection
    # Minimum rank score for particles after LoG picking
    'refine__select_logbatch__rank_threshold' : 0.35,
    # Minimum rank score for particles after Topaz picking
    'refine__select_rest__rank_threshold' : 0.35,

    ### Parameters for 2D classification (logbatch and rest)
    # Which (single) GPU to run on for logbatch and rest
    'refine__class2d_rest__gpu_ids' : '0,1',
    'refine__class2d_logbatch__gpu_ids' : '0,1',

    # Minimum number of particles in the first batch of logpicked particles to perform 2D classification on (this should be <= 'preprocess__split_logpick__split_size' above)
    'refine__logbatch_size' : 5000,
    # Diameter of the mask used for 2D classification (in Angstrom)
    'refine__class2d_logbatch__particle_diameter' : 200,
    # Stop after class2d, so skip inimodel3d and refine3d
    'refine__do_until_class2d' : False,

    # Options for inimodel3d and refine3d
    # Symmetry
    'refine__inimodel3d__sym_name' : 'C1',
    'refine__refine3d__sym_name' : 'C1',
    # Diameter of the mask used for inimodel3d and refine3d (in Angstrom)
    'refine__inimodel3d__particle_diameter' : 200,
    'refine__refine3d__particle_diameter' : 200,


    ### End of options
    }

class RelionItGui(object):

    def __init__(self, main_window, options):
        self.main_window = main_window
        self.options = options

        # Convenience function for making file browser buttons
        def new_browse_button(master, var_to_set, filetypes=(('MRC file', '*.mrc'), ('All files', '*'))):
            def browse_command():
                chosen_file = tkFileDialog.askopenfilename(filetypes=filetypes)
                if chosen_file is not None:
                    # Make path relative if it's in the current directory
                    if chosen_file.startswith(os.getcwd()):
                        chosen_file = os.path.relpath(chosen_file)
                    var_to_set.set(chosen_file)
            return tk.Button(master, text="Browse...", command=browse_command)

        ### Create GUI

        main_frame = tk.Frame(main_window)
        main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame = tk.Frame(main_frame)
        left_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        right_frame = tk.Frame(main_frame)
        right_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        ###

        project_frame = tk.LabelFrame(left_frame, text="Movie details", padx=5, pady=5)
        project_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(project_frame, 1, weight=1)

        row = 0

        tk.Label(project_frame, text="Pattern for movies:").grid(row=row, sticky=tk.W)
        self.import_images_var = tk.StringVar()  # for data binding
        self.import_images_entry = tk.Entry(project_frame, textvariable=self.import_images_var)
        self.import_images_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.import_images_entry.insert(0, self.options['preprocess__importmovies__fn_in_raw'])

        import_button = new_browse_button(project_frame, self.import_images_var,
                                          filetypes=(('Image file', '{*.mrc, *.mrcs, *.tif, *.tiff}'), ('All files', '*')))
        import_button.grid(row=row, column=2)

        row += 1
        
        tk.Label(project_frame, text="Gain reference (optional):").grid(row=row, sticky=tk.W)
        self.gainref_var = tk.StringVar()  # for data binding
        self.gainref_entry = tk.Entry(project_frame, textvariable=self.gainref_var)
        self.gainref_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.gainref_entry.insert(0, self.options['preprocess__motioncorr__fn_gain_ref'])

        new_browse_button(project_frame, self.gainref_var).grid(row=row, column=2)
        
        row += 1
         
        tk.Label(project_frame, text="Super-resolution?").grid(row=row, sticky=tk.W)
        self.superres_var = tk.IntVar()
        superres_button = tk.Checkbutton(project_frame, var=self.superres_var)
        superres_button.grid(row=row, column=1, sticky=tk.W)
        if options['preprocess__motioncorr__bin_factor'] == '2':
            superres_button.select()

        ###
        
        expt_frame = tk.LabelFrame(left_frame, text="Experimental details", padx=5, pady=5)
        expt_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        row = 0

        tk.Label(expt_frame, text="Voltage (kV):").grid(row=row, sticky=tk.W)
        self.voltage_entry = tk.Entry(expt_frame)
        self.voltage_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.voltage_entry.insert(0, str(options['preprocess__importmovies__kV']))

        row += 1
        
        tk.Label(expt_frame, text="Cs (mm):").grid(row=row, sticky=tk.W)
        self.cs_entry = tk.Entry(expt_frame)
        self.cs_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.cs_entry.insert(0, str(options['preprocess__importmovies__Cs']))

        row += 1
        
        tk.Label(expt_frame, text="Phase plate?").grid(row=row, sticky=tk.W)
        self.phaseplate_var = tk.IntVar()
        phaseplate_button = tk.Checkbutton(expt_frame, var=self.phaseplate_var)
        phaseplate_button.grid(row=row, column=1, sticky=tk.W)
        if options['preprocess__ctffind__do_phaseshift'] == 'True':
            phaseplate_button.select()

        row += 1

        tk.Label(expt_frame, text=u"Pixel size (\u212B):").grid(row=row, sticky=tk.W)
        self.angpix_var = tk.StringVar()  # for data binding
        self.angpix_entry = tk.Entry(expt_frame, textvariable=self.angpix_var)
        self.angpix_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.angpix_entry.insert(0, str(options['preprocess__importmovies__angpix']))

        row += 1
        
        tk.Label(expt_frame, text=u"Exposure rate (e\u207B / \u212B\u00B2 / frame):").grid(row=row, sticky=tk.W)
        self.exposure_entry = tk.Entry(expt_frame)
        self.exposure_entry.grid(row=row, column=1, sticky=tk.W + tk.E)
        self.exposure_entry.insert(0, str(options['preprocess__motioncorr__dose_per_frame']))

        ###

        particle_frame = tk.LabelFrame(left_frame, text="Particle details", padx=5, pady=5)
        particle_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(particle_frame, 1, weight=1)

        row = 0

        tk.Label(particle_frame, text="Symmetry:").grid(row=row, sticky=tk.W)
        self.symmetry_var = tk.StringVar()  # for data binding
        self.symmetry_entry = tk.Entry(particle_frame, textvariable=self.symmetry_var)
        self.symmetry_entry.grid(row=row, column=1, sticky=tk.W)
        self.symmetry_entry.insert(0, str(options['refine__inimodel3d__sym_name']))

        row += 1

        tk.Label(particle_frame, text="Nr particles per micrograph:").grid(row=row, sticky=tk.W)
        self.partspermic_var = tk.StringVar()  # for data binding
        self.partspermic_entry = tk.Entry(particle_frame, textvariable=self.partspermic_var)
        self.partspermic_entry.grid(row=row, column=1, sticky=tk.W)
        self.partspermic_entry.insert(0, str(options['refine__train_topaz__topaz_nr_particles']))

        row += 1

        tk.Label(particle_frame, text=u"Longest diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.particle_max_diam_var = tk.StringVar()  # for data binding
        self.particle_max_diam_entry = tk.Entry(particle_frame, textvariable=self.particle_max_diam_var)
        self.particle_max_diam_entry.grid(row=row, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_max_diam_entry.insert(0, str(options['preprocess__logpicker__log_diam_max']))

        row += 1

        tk.Label(particle_frame, text=u"Shortest diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.particle_min_diam_entry = tk.Entry(particle_frame)
        self.particle_min_diam_entry.grid(row=row, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_min_diam_entry.insert(0, str(options['preprocess__logpicker__log_diam_min']))

        row += 1
        
        tk.Label(particle_frame, text=u"Mask diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.mask_diameter_var = tk.StringVar()  # for data binding
        self.mask_diameter_entry = tk.Entry(particle_frame, textvariable=self.mask_diameter_var)
        self.mask_diameter_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.mask_diameter_entry.insert(0, str(options['refine__class2d_logbatch__particle_diameter']))
        self.mask_diameter_px = tk.Label(particle_frame, text="= NNN px")
        self.mask_diameter_px.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(particle_frame, text="Box size (px):").grid(row=row, sticky=tk.W)
        self.box_size_var = tk.StringVar()  # for data binding
        self.box_size_entry = tk.Entry(particle_frame, textvariable=self.box_size_var)
        self.box_size_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.box_size_entry.insert(0, str(options['preprocess__extract_logpick__extract_size']))
        self.box_size_in_angstrom = tk.Label(particle_frame, text=u"= NNN \u212B")
        self.box_size_in_angstrom.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(particle_frame, text="Down-sample to (px):").grid(row=row, sticky=tk.W)
        self.extract_small_boxsize_var = tk.StringVar()  # for data binding
        self.extract_small_boxsize_entry = tk.Entry(particle_frame, textvariable=self.extract_small_boxsize_var)
        self.extract_small_boxsize_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.extract_small_boxsize_entry.insert(0, str(options['preprocess__extract_logpick__rescale']))
        self.extract_angpix = tk.Label(particle_frame, text=u"= NNN \u212B/px")
        self.extract_angpix.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(particle_frame, text="Calculate for me:").grid(row=row, sticky=tk.W)
        self.auto_boxsize_var = tk.IntVar()
        auto_boxsize_button = tk.Checkbutton(particle_frame, var=self.auto_boxsize_var)
        auto_boxsize_button.grid(row=row, column=1, sticky=tk.W)
        auto_boxsize_button.select()

        ###

        picking_frame = tk.LabelFrame(right_frame, text="Picking details", padx=5, pady=5)
        picking_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(picking_frame, 1, weight=1)

        row = 0

        tk.Label(picking_frame, text="Nr particles for LoG picking:").grid(row=row, sticky=tk.W)
        self.logbatch_var = tk.StringVar()  # for data binding
        self.logbatch_entry = tk.Entry(picking_frame, textvariable=self.logbatch_var)
        self.logbatch_entry.grid(row=row, column=1, sticky=tk.W)
        self.logbatch_entry.insert(0, str(options['preprocess__split_logpick__split_size']))

        row += 1

        tk.Label(picking_frame, text="LoG picking threshold:").grid(row=row, sticky=tk.W)
        self.log_thresh_var = tk.StringVar()  # for data binding
        self.log_thresh_entry = tk.Entry(picking_frame, textvariable=self.log_thresh_var)
        self.log_thresh_entry.grid(row=row, column=1, sticky=tk.W)
        self.log_thresh_entry.insert(0, str(options['preprocess__logpicker__log_adjust_thr']))

        row += 1

        tk.Label(picking_frame, text="LoG class2d score:").grid(row=row, sticky=tk.W)
        self.log_classscore_var = tk.StringVar()  # for data binding
        self.log_classscore_entry = tk.Entry(picking_frame, textvariable=self.log_classscore_var)
        self.log_classscore_entry.grid(row=row, column=1, sticky=tk.W)
        self.log_classscore_entry.insert(0, str(options['refine__select_logbatch__rank_threshold']))

        row += 1

        tk.Label(picking_frame, text="Topaz picking threshold:").grid(row=row, sticky=tk.W)
        self.topaz_thresh_var = tk.StringVar()  # for data binding
        self.topaz_thresh_entry = tk.Entry(picking_frame, textvariable=self.topaz_thresh_var)
        self.topaz_thresh_entry.grid(row=row, column=1, sticky=tk.W)
        self.topaz_thresh_entry.insert(0, str(options['preprocess__extract_topazpick__minimum_pick_fom']))

        row += 1

        tk.Label(picking_frame, text="Topaz class2d score:").grid(row=row, sticky=tk.W)
        self.topaz_classscore_var = tk.StringVar()  # for data binding
        self.topaz_classscore_entry = tk.Entry(picking_frame, textvariable=self.topaz_classscore_var)
        self.topaz_classscore_entry.grid(row=row, column=1, sticky=tk.W)
        self.topaz_classscore_entry.insert(0, str(options['refine__select_rest__rank_threshold']))

        ###

        compute_frame = tk.LabelFrame(right_frame, text="Computation details", padx=5, pady=5)
        compute_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        row = 0

        tk.Label(compute_frame, text="Do Autopick & Class2D?").grid(row=row, sticky=tk.W)
        self.do_class2d_var = tk.IntVar()
        do_class2d_button = tk.Checkbutton(compute_frame, var=self.do_class2d_var)
        do_class2d_button.grid(row=row, column=1, sticky=tk.W)
        if not options['preprocess__do_until_ctf']:
            do_class2d_button.select()

        row += 1
        
        tk.Label(compute_frame, text="Do Inimodel3D & Refine3D?").grid(row=row, sticky=tk.W)
        self.do_refine3d_var = tk.IntVar()
        do_refine3d_button = tk.Checkbutton(compute_frame, var=self.do_refine3d_var)
        do_refine3d_button.grid(row=row, column=1, sticky=tk.W)
        if not options['refine__do_until_class2d']:
            do_refine3d_button.select()

        row += 1
        
        tk.Label(compute_frame, text="Nr CPU cores:").grid(row=row, sticky=tk.W)
        self.nthreads_var = tk.StringVar()  # for data binding
        self.nthreads_entry = tk.Entry(compute_frame, textvariable=self.nthreads_var)
        self.nthreads_entry.grid(row=row, column=1, sticky=tk.W)
        self.nthreads_entry.insert(0, str(options['preprocess__motioncorr__nr_threads']))

        row += 1

        tk.Label(compute_frame, text="GPUs (comma-separated):").grid(row=row, sticky=tk.W)
        self.gpu_var = tk.StringVar()  # for data binding
        self.gpu_entry = tk.Entry(compute_frame, textvariable=self.gpu_var)
        self.gpu_entry.grid(row=row, column=1, sticky=tk.W)
        self.gpu_entry.insert(0, str(options['refine__class2d_rest__gpu_ids']))

         ### Add logic to the box size boxes

        def calculate_box_size(particle_size_pixels):
            # Use box 50% larger than largest particle diameter and ensure size is even
            box_size_exact = 1.5 * particle_size_pixels
            box_size_int = int(math.ceil(box_size_exact))
            return box_size_int + box_size_int % 2

        def calculate_downscaled_box_size(box_size_pix, angpix):
            for small_box_pix in (48, 64, 96, 128, 160, 192, 256, 288, 300, 320, 360,
                                  384, 400, 420, 450, 480, 512, 640, 768, 896, 1024):
                # Don't go larger than the original box
                if small_box_pix > box_size_pix:
                    return box_size_pix
                # If Nyquist freq. is better than 8.5 A, use this downscaled box, otherwise continue to next size up
                small_box_angpix = angpix * box_size_pix / small_box_pix
                if small_box_angpix < 4.25:
                    return small_box_pix
            # Fall back to a warning message
            return "Box size is too large!"

        def update_box_size_labels(*args_ignored, **kwargs_ignored):
            try:
                angpix = float(self.angpix_entry.get())
                if self.get_var_as_bool(self.superres_var):
                    angpix = angpix * 2
            except ValueError:
                # Can't update any of the labels without angpix
                self.mask_diameter_px.config(text="= NNN px")
                self.box_size_in_angstrom.config(text=u"= NNN \u212B")
                self.extract_angpix.config(text=u"= NNN \u212B/px")
                return
            try:
                mask_diameter = float(self.mask_diameter_entry.get())
                mask_diameter_px = mask_diameter / angpix
                self.mask_diameter_px.config(text="= {:.1f} px".format(mask_diameter_px))
            except (ValueError, ZeroDivisionError):
                self.mask_diameter_px.config(text="= NNN px")
                # Don't return - an error here doesn't stop us calculating the other labels
            try:
                box_size = float(self.box_size_entry.get())
                box_angpix = angpix * box_size
                self.box_size_in_angstrom.config(text=u"= {:.1f} \u212B".format(box_angpix))
            except ValueError:
                # Can't update these without the box size
                self.box_size_in_angstrom.config(text=u"= NNN \u212B")
                self.extract_angpix.config(text=u"= NNN \u212B/px")
                return
            try:
                extract_small_boxsize = float(self.extract_small_boxsize_entry.get())
                small_box_angpix = box_angpix / extract_small_boxsize
                self.extract_angpix.config(text=u"= {:.3f} \u212B/px".format(small_box_angpix))
            except (ValueError, ZeroDivisionError):
                # Can't update the downscaled pixel size unless the downscaled box size is valid
                self.extract_angpix.config(text=u"= NNN \u212B/px")

        def update_box_sizes(*args_ignored, **kwargs_ignored):
            # Always activate entry boxes - either we're activating them anyway, or we need to edit the text.
            # For text editing we need to activate the box first then deactivate again afterwards.
            self.mask_diameter_entry.config(state=tk.NORMAL)
            self.box_size_entry.config(state=tk.NORMAL)
            self.extract_small_boxsize_entry.config(state=tk.NORMAL)
            if self.get_var_as_bool(self.auto_boxsize_var):
                try:
                    particle_size_angstroms = float(self.particle_max_diam_entry.get())
                    mask_diameter = 1.1 * particle_size_angstroms
                    self.mask_diameter_entry.delete(0, tk.END)
                    self.mask_diameter_entry.insert(0, str(mask_diameter))
                    angpix = float(self.angpix_entry.get())
                    if self.get_var_as_bool(self.superres_var):
                        angpix = angpix * 2
                    particle_size_pixels = particle_size_angstroms / angpix
                    box_size = calculate_box_size(particle_size_pixels)
                    self.box_size_entry.delete(0, tk.END)
                    self.box_size_entry.insert(0, str(box_size))
                    small_boxsize = calculate_downscaled_box_size(int(box_size), angpix)
                    self.extract_small_boxsize_entry.delete(0, tk.END)
                    self.extract_small_boxsize_entry.insert(0, str(small_boxsize))
                except:
                    # Ignore errors - they will be picked up if the user tries to save the options
                    pass
                self.mask_diameter_entry.config(state=tk.DISABLED)
                self.box_size_entry.config(state=tk.DISABLED)
                self.extract_small_boxsize_entry.config(state=tk.DISABLED)
            update_box_size_labels()

        self.box_size_var.trace('w', update_box_size_labels)
        self.extract_small_boxsize_var.trace('w', update_box_size_labels)

        self.angpix_var.trace('w', update_box_sizes)
        self.particle_max_diam_var.trace('w', update_box_sizes)
        auto_boxsize_button.config(command=update_box_sizes)

        button_frame = tk.Frame(right_frame)
        button_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)

        self.save_button = tk.Button(button_frame, text="Save options", command=self.save_options)
        self.save_button.pack(padx=5, pady=5, side=tk.RIGHT)

        self.run_button = tk.Button(button_frame, text="Save & run", command=self.run_pipeline)
        self.run_button.pack(padx=5, pady=5, side=tk.RIGHT)

        # Show initial pixel sizes
        update_box_sizes()

    def get_var_as_bool(self, var):
        """Helper function to convert a Tk IntVar (linked to a checkbox) to a boolean value"""
        return True if var.get() == 1 else False

    def fetch_options_from_gui(self):
        """
        Fetch the current values from the GUI widgets and store them in the options object.

        Returns:
            A list of warning messages about possible incorrect option values.

        Raises:
            ValueError: If an option value is invalid.
        """

        opts = self.options
        warnings = []

        opts['preprocess__do_until_ctf'] = not self.get_var_as_bool(self.do_class2d_var)

        opts['refine__do_until_class2d'] = not self.get_var_as_bool(self.do_refine3d_var)

        try:
            opts['preprocess__importmovies__kV'] = float(self.voltage_entry.get())
        except ValueError:
            raise ValueError("Voltage must be a number")
        if opts['preprocess__importmovies__kV'] <= 0.0:
            warnings.append("- Voltage should be a positive number")

        try:
            opts['preprocess__importmovies__Cs'] = float(self.cs_entry.get())
        except ValueError:
            raise ValueError("Cs must be a number")

        opts['preprocess__ctffind__do_phaseshift'] = self.get_var_as_bool(self.phaseplate_var)

        try:
            opts['preprocess__importmovies__angpix'] = float(self.angpix_entry.get())
        except ValueError:
            raise ValueError("Pixel size must be a number")
        if opts['preprocess__importmovies__angpix'] <= 0.0:
            warnings.append("- Pixel size should be a positive number")

        try:
            opts['preprocess__motioncorr__nr_threads'] = int(self.nthreads_entry.get())
        except ValueError:
            raise ValueError("Number of threads for RelionCor2 must be a number")

        try:
            opts['preprocess__motioncorr__dose_per_frame'] = float(self.exposure_entry.get())
        except ValueError:
            raise ValueError("Exposure rate must be a number")
        if opts['preprocess__motioncorr__dose_per_frame'] <= 0.0:
            warnings.append("- Exposure rate should be a positive number")

        if self.get_var_as_bool(self.superres_var):
            opts['preprocess__motioncorr__bin_factor'] = 2;
        else:
            opts['preprocess__motioncorr__bin_factor'] = 1;

        try:
            opts['preprocess__logpicker__log_diam_max'] = float(self.particle_max_diam_entry.get())
            opts['refine__train_topaz__topaz_particle_diameter'] = float(self.particle_max_diam_entry.get())
            opts['preprocess__topazpicker__topaz_particle_diameter'] = float(self.particle_max_diam_entry.get())
        except ValueError:
            if len(self.particle_max_diam_entry.get()) == 0 and opts['preprocess__do_until_ctf'] == 'True':
                # This was left blank and won't be used, set to zero to avoid errors in calculations later
                opts['preprocess__logpicker__log_diam_max'] = 0.0
            else:
                raise ValueError("Particle longest diameter must be a number")

        try:
            opts['preprocess__logpicker__log_diam_min'] = float(self.particle_min_diam_entry.get())
        except ValueError:
            if len(self.particle_min_diam_entry.get()) == 0 and opts['preprocess__do_until_ctf'] == 'True':
                # This was left blank and won't be used, set to zero to avoid errors in calculations later
                opts['preprocess__logpicker__log_diam_min'] = 0.0
            else:
                raise ValueError("Particle shortest diameter must be a number")

        try:
            opts['refine__class2d_logbatch__particle_diameter'] = float(self.mask_diameter_entry.get())
            opts['refine__class2d_rest__particle_diameter'] = float(self.mask_diameter_entry.get())
            opts['refine__inimodel3d__particle_diameter'] = float(self.mask_diameter_entry.get())
            opts['refine__refine3d__particle_diameter'] = float(self.mask_diameter_entry.get())
        except ValueError:
            raise ValueError("Mask diameter must be a number")
        if opts['refine__class2d_logbatch__particle_diameter'] <= 0:
            warnings.append("- Mask diameter should be a positive number")

        try:
            opts['preprocess__extract_logpick__extract_size'] = int(self.box_size_entry.get())
            opts['preprocess__extract_topazpick__extract_size'] = int(self.box_size_entry.get())
        except ValueError:
            raise ValueError("Box size must be a number")
        if opts['preprocess__extract_logpick__extract_size'] <= 0:
            warnings.append("- Box size should be a positive number")

        try:
            opts['preprocess__extract_logpick__rescale'] = int(self.extract_small_boxsize_entry.get())
            opts['preprocess__extract_logpick__do_rescale'] = True
            opts['preprocess__extract_topazpick__rescale'] = int(self.extract_small_boxsize_entry.get())
            opts['preprocess__extract_topazpick__do_rescale'] = True
        except ValueError:
            raise ValueError("Down-sampled box size must be a number")
        if opts['preprocess__extract_logpick__rescale'] <= 0:
            warnings.append("- Down-sampled box size should be a positive number")

        opts['preprocess__importmovies__fn_in_raw'] = self.import_images_entry.get()
        if opts['preprocess__importmovies__fn_in_raw'].startswith(('/', '..')):
            warnings.append("- Movies should be located inside the project directory")
        if '*' not in opts['preprocess__importmovies__fn_in_raw']:
            warnings.append("- Pattern for input movies should normally contain a '*' to select more than one file")

        opts['preprocess__motioncorr__fn_gain_ref'] = self.gainref_entry.get()
        if len(opts['preprocess__motioncorr__fn_gain_ref']) > 0 and not os.path.isfile(opts['preprocess__motioncorr__fn_gain_ref']):
            warnings.append("- Gain reference file '{}' does not exist".format(opts['preprocess__motioncorr__fn_gain_ref']))

        try:
            opts['preprocess__split_logpick__split_size'] = int(self.logbatch_entry.get())
            opts['refine__logbatch_size'] = int(self.logbatch_entry.get())
        except ValueError:
            raise ValueError("Nr particles for topaz training must be a number")
        if opts['preprocess__split_logpick__split_size'] <= 0:
            warnings.append("- Nr particles for topaz training should be a positive number")

        try:
            opts['preprocess__logpicker__log_adjust_thr'] = float(self.log_thresh_var.get())
        except ValueError:
            raise ValueError("LoG picking threshold must be a number")

        try:
            opts['refine__select_logbatch__rank_threshold'] = float(self.log_classscore_var.get())
        except ValueError:
            raise ValueError("LoG class2d score must be a number")

        try:
            opts['preprocess__extract_topazpick__minimum_pick_fom'] = float(self.topaz_thresh_var.get())
        except ValueError:
            raise ValueError("Topaz picking threshold must be a number")

        try:
            opts['refine__select_rest__rank_threshold'] = float(self.topaz_classscore_var.get())
        except ValueError:
            raise ValueError("Topaz class2d score must be a number")

        try:
            opts['refine__train_topaz__topaz_nr_particles'] = int(self.partspermic_entry.get())
            opts['preprocess__topazpicker__topaz_nr_particles'] = int(self.partspermic_entry.get())
        except ValueError:
            raise ValueError("Nr particles per micrograph must be a number")
        if opts['refine__train_topaz__topaz_nr_particles'] <= 0:
            warnings.append("- Nr particles per micrograph should be a positive number")

        opts['refine__inimodel3d__sym_name'] = self.symmetry_entry.get()
        opts['refine__refine3d__sym_name'] = self.symmetry_entry.get()

        try:
            opts['preprocess__topazpicker__nr_mpi'] = int(self.mpi_topaz_entry.get())
        except ValueError:
            raise ValueError("Nr MPI for topaz picking must be a single number")
        
        opts['refine__train_topaz__gpu_ids'] = (self.gpu_entry.get()).split(',')[0]
        opts['refine__class2d_logbatch__gpu_ids'] = self.gpu_entry.get()
        opts['refine__class2d_rest__gpu_ids'] = self.gpu_entry.get()
        opts['refine__inimodel3d__gpu_ids'] = self.gpu_entry.get()
        opts['refine__refine3d__gpu_ids'] = (self.gpu_entry.get()).replace(',',':')


        return warnings

    def save_options(self):
        """
        Update the full set of options from the values in the GUI, and save them to a file.

        Returns:
            True if the options were valid and saved successfully, otherwise False.
        """
        try:
            warnings = self.fetch_options_from_gui()
            if len(warnings) == 0 or tkMessageBox.askokcancel("Warning", "\n".join(warnings), icon='warning',
                                                              default=tkMessageBox.CANCEL):

                # Write the current options to a .py file
                with open(OPTIONS_FILE, 'w') as file:
                    file.write("{\n") 
                    for k,v in self.options.items():
                        file.write("'%s' : '%s', \n" % (k, v))
                    file.write("}\n")                  

                print(" RELION_IT: Written all options to {}".format(OPTIONS_FILE))
                
                # loop over all options and change the scheduler STAR files
                for option, value in self.options.items():
                    
                    if (value == ''):
                        value = '\\"\\"'
                    
                    # Set variables in schedule.star
                    if option.count('__') == 1:
                        splits = option.split('__')
                        schedulename = splits[0]
                        varname = splits[1]
                        schedulestar = 'Schedules/' + schedulename + '/schedule.star'
                        if not os.path.isfile(schedulestar):
                            message = 'Error: ' + schedulestar + ' does not exist'
                            print(message)
                            tkMessageBox.showerror(message)
                            return False
                    
                        command = 'relion_scheduler --schedule ' + schedulename + ' --set_var ' + varname + ' --value \"' + str(value) + '\"' + ' --original_value \"' + str(value) + '\"'
                        print(' RELION_IT: excuting: ', command)
                        os.system(command)
                        
                    # Set joboptions in job.star
                    elif option.count('__') == 2:
                        splits = option.split('__')
                        schedulename = splits[0]
                        jobname = splits[1]
                        joboption = splits[2]
                        jobstar = 'Schedules/' + schedulename + '/' + jobname + '/job.star'
                        if not os.path.isfile(jobstar):
                            message = 'Error: ' + jobstar + 'does not exist'
                            print(message)
                            tkMessageBox.showerror(message)
                        
                        command = 'relion_pipeliner --editJob ' + jobstar + ' --editOption ' + joboption + ' --editValue \"' + str(value) + '\"'
                        print(' RELION_IT: excuting: ', command)
                        os.system(command)

                return True

        except Exception as ex:
            tkMessageBox.showerror("Error", ex.message)
            traceback.print_exc()
        return False

    def run_pipeline(self):
        """
        Update the full set of options from the values in the GUI, close the GUI and run the pipeline.
        """
        if self.save_options():
            self.main_window.destroy()
            run_scheduler(self.options, True) #True means launch the RELION GUI
 
def run_scheduler(options, do_gui):

    command = 'relion_scheduler --schedule preprocess --reset &'
    print(' RELION_IT: excuting: ', command)
    os.system(command)

    command = 'relion_scheduler --schedule preprocess --run --pipeline_control Schedules/preprocess/ >> Schedules/preprocess/run.out 2>> Schedules/preprocess/run.err &'
    print(' RELION_IT: excuting: ', command)
    os.system(command)

    if not options['preprocess__do_until_ctf'] == 'True':

        command = 'relion_scheduler --schedule refine --reset &'
        print(' RELION_IT: excuting: ', command)
        os.system(command)

        command = 'relion_scheduler --schedule refine --run  --pipeline_control Schedules/refine/ >> Schedules/refine/run.out 2>> Schedules/refine/run.err  &'
        print(' RELION_IT: excuting: ', command)
        os.system(command)

    print(' RELION_IT: Now monitor the preprocess and refine Schedules from the RELION GUI ...')

    if do_gui:
        command = 'relion --do_projdir &'
        os.system(command)


def copy_schedule(schedulename):

    ## Only copy the Schedule directory structure from the RELION installation directory if it doesn't exist yet
    if not os.path.isdir('Schedules/'+schedulename):
        mydir = os.path.dirname(os.path.realpath(__file__))
        print(' RELION_IT: copying Schedules/' + schedulename + ' from: ' + mydir)
        copytree(mydir+'/Schedules/'+schedulename, 'Schedules/'+schedulename)


def main():
    """
    Run the RELION 3.2 Scheduler.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """
     # Start by parsing arguments
    # (If --help is given, the program will print a usage message and exit)
    parser = argparse.ArgumentParser()
    parser.add_argument("extra_options", nargs="*", metavar="extra_options.py",
                        help="Python files containing options for relion_it.py")
    parser.add_argument("--nogui", action="store_true", help="don't launch GUI to set options, execute non-interactively")
    args = parser.parse_args()

    print(' RELION_IT: -------------------------------------------------------------------------------------------------------------------')
    print(' RELION_IT: script for automated, on-the-fly single-particle analysis in RELION (>= 3.2)')
    print(' RELION_IT: authors: Sjors H.W. Scheres, Takanori Nakane & Colin M. Palmer')
    print(' RELION_IT: ')
    print(' RELION_IT: usage: ./relion_it.py [extra_options.py [extra_options2.py ....] ] [--nogui]')
    print(' RELION_IT: ')
    print(' RELION_IT: -------------------------------------------------------------------------------------------------------------------')
    print(' RELION_IT: ')
    

    opts = collections.OrderedDict(RelionItOptions)

    """
    for user_opt_file in args.extra_options:
        print(' RELION_IT: reading options from {}'.format(user_opt_file))
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)
    """    
    for user_opt_file in args.extra_options:
        print(' RELION_IT: reading options from {}'.format(user_opt_file))
        with open(user_opt_file) as file: 
            user_opts = collections.OrderedDict(ast.literal_eval(file.read()))
            opts.update(user_opts)

    # Copy Schedules over from RELION directory if they dont exit
    copy_schedule('preprocess')
    copy_schedule('refine')

    if args.nogui:
        run_scheduler(opts, False)
    else:
        print(' RELION_IT: launching GUI...')
        tk_root = tk.Tk()
        tk_root.title("relion_it.py setup")
        RelionItGui(tk_root, opts)
        tk_root.mainloop()

if __name__ == "__main__":
    main()




