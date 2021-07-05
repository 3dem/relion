#!/usr/bin/env python2.7
"""
relion_it.py
============

Simple GUI to set up RELION-4.0 schemer.

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

Use double underscores to separate SCHEMENAME__JOBNAME__JOBOPTION
 E.g. for SCHEMENAME=prep, JOBNAME=importmovies and JOBOPTION=angpix
'prep__importmovies__angpix' defines the value for the 'angpix' option in the file 'Schemes/prep/importmovies/job.star'

Likewise, use SCHEMENAME__VARIABLENAME for Scheme variables
 E.g. 'proc__inibatch_size' refers to the variablename 'inibatch_size' in the file 'Schemes/proc/scheme.star'

This python script will set the corresponding values of these joboption and scheme-variable values in the Scheme STAR files. 
The RELION schemer can then be used to run the schemer.

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
and the schemers are launched from this python script. This will allow a non-interactive use of this script.



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
    # Use double underscores to separate SCHEMENAME__JOBNAME__JOBOPTION
    #
    # E.g. for SCHEMENAME=prep, JOBNAME=importmovies and JOBOPTION=angpix
    #  'prep__importmovies__angpix' defines the value for the 'angpix' option in the file 'Schemes/prep/importmovies/job.star'
    #
    # Likewise, use SCHEMENAME__VARIABLENAME for Scheme variables
    #
    # E.g. 'proc__do_log' refers to the variablename 'do_log' in the file 'Schemes/proc/scheme.star'
    #
    #  This python script will modify the joboption and scheme-variable values in these files
    #
    #############################################################################

    # Perform Import, Motion correction and CTF estimation?
    'do_prep' : True,
    # If the above option is false, then use this to provide a specific star file to the micrographs instead
    'proc__ctffind_mics' : 'Schemes/prep/ctffind/micrographs_ctf.star',
    # Perform autopicking and processing of particles?
    'do_proc' : True,

    ### How many micrograph movies to do in each iteration of the prep Scheme
    'prep__do_at_most' : 50,

    ### General parameters
    # Pixel size in Angstroms in the input movies
    'prep__importmovies__angpix' : 0.885,
    # Acceleration voltage (in kV)
    'prep__importmovies__kV' : 200,
    # Polara = 2.0; Talos/Krios = 2.7; some Cryo-ARM = 1.4 
    'prep__importmovies__Cs' : 1.4,

    ### Import images (Linux wild card; movies as *.mrc, *.mrcs, *.tiff or *.tif; single-frame micrographs as *.mrc)
    'prep__importmovies__fn_in_raw' : 'Movies/*.tiff',
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    'prep__importmovies__is_multiframe' : True,

    ### MotionCorrection parameters
    # Dose in electrons per squared Angstrom per fraction
    'prep__motioncorr__dose_per_frame' : 1.277,
    # Gain-reference image in MRC format (only necessary if input movies are not yet gain-corrected, e.g. compressed TIFFs from K2)
    'prep__motioncorr__fn_gain_ref' : 'Movies/gain.mrc',
    # EER fractionation. The dose rate (motioncor_doseperframe) is e/A2/fraction after this fractionation.
    'prep__motioncorr__eer_grouping' : 20,
    # For super-resolution movies, bin by factor of 2 straight away, as we don't believe in super-resolution
    'prep__motioncorr__bin_factor' : 1,
    # Number of threads to use for Relioncor2
    'prep__motioncorr__nr_threads' : 16,

    ### CTF estimation parameters
    # Name of the ctffind executable
    'prep__ctffind__fn_ctffind_exe' : '/public/EM/ctffind/ctffind.exe',
    # Most cases won't need changes here...
    'prep__ctffind__do_phaseshift' : False,

    ### Selection of micrographs based on CTF resolution 
    'proc__select_mics__select_maxval' : 6,

    # Calculate boxsize automatically for me?
    'do_auto_boxsize' : True,

    ### Parameters for autopicking
    ## Shortest diameter for LoG picking
    'proc__autopick__log_diam_min' : 150.0,
    ## Longest diameter for LoG picking
    'proc__autopick__log_diam_max' : 180.0,
    # Adjust threshold for LoG picking
    'proc__autopick__log_adjust_thr' : 0.0,
    # Model for topaz picking (leave empty for default model)
    'proc__topaz_model' : '',
    # Name of the topaz executable
    'proc__autopick__fn_topaz_exec' : '/public/EM/RELION/topaz',
    # Additional arguments to be passed onto topaz
    'proc__autopick__topaz_other_args' : '',
    # Use GPUs for autopicking? (True for topaz, False for LoG)
    'proc__autopick__use_gpu' : False,
    # Which (single) GPU to run topaz on 
    'proc__autopick__gpu_ids' : 0,
    # How many MPI processes to use in parallel for picking
    'proc__autopick__nr_mpi' : 4,

    ### Extract particles
    # Box size of particles in the averaged micrographs (in pixels)
    'proc__extract__extract_size' : 256,
    # Down-scale the particles upon extraction?
    'proc__extract__do_rescale' : True,
    # Box size of the down-scaled particles (in pixels)
    'proc__extract__rescale' : 64,
    # Use FOM threshold for extraction (will be set for True for topaz, False for LoG picking)
    'proc__extract__do_fom_threshold' : False,
    # Minimum FOM for topaz extraction
    'proc__extract__minimum_pick_fom' : 0,

    ### Perform LoG picking instead of topaz 
    'proc__do_log' : True,


    ### Parameters for automated 2D class selection
    # Minimum rank score for selected classes
    'proc__select_parts__rank_threshold' : 0.25,
    # Minimum number of selected classes (not used on GUI)
    'proc__select_parts__select_nr_classes' : 0,

    ### Parameters for 2D classification
    # Which (comma-separated) GPUs to run on 
    'proc__class2d__gpu_ids' : '0,1',
    # Diameter of the mask used for 2D classification (in Angstrom)
    'proc__class2d__particle_diameter' : 200,

    # Minimum number of particles to continue with Inimodel3D and Refine3D after Class2D?
    'proc__min_nr_parts_3d' : 5000,

    # Name for initial 3D reference for Refine3D. If this file does not exist, a 3D Initial Model job will be run
    'proc__iniref' : 'None',

    # Options for inimodel3d and refine3d
    # Symmetry
    'proc__inimodel3d__sym_name' : 'D2',
    'proc__refine3d__sym_name' : 'D2',
    # Diameter of the mask used for inimodel3d and refine3d (in Angstrom)
    'proc__inimodel3d__particle_diameter' : 200,
    'proc__refine3d__particle_diameter' : 200,


    ### End of options
    }

class RelionItGui(object):

    def __init__(self, main_window, options):
        self.main_window = main_window
        self.options = options

        ### Create GUI

        # Colour definitions
        # Yellowish background for entries
        entry_bg = '#ffffe6'
        # reddish colour for Browse buttons
        button_bg = '#c8506e'
        # Derker red for run buttons
        runbutton_bg = '#a01e3c'

        # Convenience function for making file browser buttons
        def new_browse_button(master, var_to_set, filetypes=(('MRC file', '*.mrc'), ('All files', '*'))):
            def browse_command():
                chosen_file = tkFileDialog.askopenfilename(filetypes=filetypes)
                if chosen_file is not None:
                    # Make path relative if it's in the current directory
                    if chosen_file.startswith(os.getcwd()):
                        chosen_file = os.path.relpath(chosen_file)
                    var_to_set.set(chosen_file)
            return tk.Button(master, text="Browse", command=browse_command, bg=button_bg)

        main_frame = tk.Frame(main_window)
        main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame = tk.Frame(main_frame)
        left_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        right_frame = tk.Frame(main_frame)
        right_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        ###

        compute_frame = tk.LabelFrame(left_frame, text="Computation settings", padx=5, pady=5)
        compute_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(compute_frame, 1, weight=1)

        row = 0

        tk.Label(compute_frame, text="Do micrograph preprocessing?").grid(row=row, sticky=tk.W)
        self.do_prep_var = tk.IntVar()
        self.do_prep_button = tk.Checkbutton(compute_frame, var=self.do_prep_var)
        self.do_prep_button.grid(row=row, column=1, sticky=tk.W)
        if options['do_prep']:
            self.do_prep_button.select()

        row += 1
        
        tk.Label(compute_frame, text="micrographs_ctf.star:").grid(row=row, sticky=tk.W)
        self.mics_var = tk.StringVar()  # for data binding
        self.mics_entry = tk.Entry(compute_frame, textvariable=self.mics_var, bg=entry_bg)
        self.mics_entry.grid(row=row, column=1, sticky=tk.W)
        self.mics_entry.insert(0, str(options['proc__ctffind_mics']))

        mics_button = new_browse_button(compute_frame, self.mics_var, filetypes=(('STAR file', '{*.star}'), ('All files', '*')))
        mics_button.grid(row=row, column=2)

        row += 1
 
        tk.Label(compute_frame, text="Do particle processsing?").grid(row=row, sticky=tk.W)
        self.do_proc_var = tk.IntVar()
        self.do_proc_button = tk.Checkbutton(compute_frame, var=self.do_proc_var)
        self.do_proc_button.grid(row=row, column=1, sticky=tk.W)
        if options['do_proc']:
            self.do_proc_button.select()

        row += 1
        
        tk.Label(compute_frame, text="GPUs (comma-separated):").grid(row=row, sticky=tk.W)
        self.gpu_var = tk.StringVar()  # for data binding
        self.gpu_entry = tk.Entry(compute_frame, textvariable=self.gpu_var, bg=entry_bg)
        self.gpu_entry.grid(row=row, column=1, sticky=tk.W)
        self.gpu_entry.insert(0, str(options['proc__class2d__gpu_ids']))

        ###

        self.prep_frame = tk.LabelFrame(left_frame, text="Preprocessing settings", padx=5, pady=5)
        self.prep_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.prep_frame, 1, weight=1)

        row = 0

        tk.Label(self.prep_frame, text="Pattern for movies:").grid(row=row, sticky=tk.W)
        self.import_images_var = tk.StringVar()  # for data binding
        self.import_images_entry = tk.Entry(self.prep_frame, textvariable=self.import_images_var, bg=entry_bg)
        self.import_images_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.import_images_entry.insert(0, self.options['prep__importmovies__fn_in_raw'])

        import_button = new_browse_button(self.prep_frame, self.import_images_var,
                                          filetypes=(('Image file', '{*.mrc, *.mrcs, *.tif, *.tiff}'), ('All files', '*')))
        import_button.grid(row=row, column=2)

        row += 1
        
        tk.Label(self.prep_frame, text="Gain reference (optional):").grid(row=row, sticky=tk.W)
        self.gainref_var = tk.StringVar()  # for data binding
        self.gainref_entry = tk.Entry(self.prep_frame, textvariable=self.gainref_var, bg=entry_bg)
        self.gainref_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.gainref_entry.insert(0, self.options['prep__motioncorr__fn_gain_ref'])

        new_browse_button(self.prep_frame, self.gainref_var).grid(row=row, column=2)
        
        row += 1
         
        tk.Label(self.prep_frame, text="Super-resolution?").grid(row=row, sticky=tk.W)
        self.superres_var = tk.IntVar()
        superres_button = tk.Checkbutton(self.prep_frame, var=self.superres_var)
        superres_button.grid(row=row, column=1, sticky=tk.W)
        if options['prep__motioncorr__bin_factor'] == '2':
            superres_button.select()

        row += 1

        tk.Label(self.prep_frame, text="Voltage (kV):").grid(row=row, sticky=tk.W)
        self.voltage_entry = tk.Entry(self.prep_frame, bg=entry_bg)
        self.voltage_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.voltage_entry.insert(0, str(options['prep__importmovies__kV']))

        row += 1
        
        tk.Label(self.prep_frame, text="Cs (mm):").grid(row=row, sticky=tk.W)
        self.cs_entry = tk.Entry(self.prep_frame, bg=entry_bg)
        self.cs_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.cs_entry.insert(0, str(options['prep__importmovies__Cs']))

        row += 1
        
        tk.Label(self.prep_frame, text="Phase plate?").grid(row=row, sticky=tk.W)
        self.phaseplate_var = tk.IntVar()
        phaseplate_button = tk.Checkbutton(self.prep_frame, var=self.phaseplate_var)
        phaseplate_button.grid(row=row, column=1, sticky=tk.W)
        if options['prep__ctffind__do_phaseshift']:
            phaseplate_button.select()

        row += 1

        tk.Label(self.prep_frame, text=u"(Super-res) pixel size (\u212B):").grid(row=row, sticky=tk.W)
        self.angpix_var = tk.StringVar()  # for data binding
        self.angpix_entry = tk.Entry(self.prep_frame, textvariable=self.angpix_var, bg=entry_bg)
        self.angpix_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.angpix_entry.insert(0, str(options['prep__importmovies__angpix']))

        row += 1
        
        tk.Label(self.prep_frame, text=u"Exposure rate (e\u207B / \u212B\u00B2 / frame):").grid(row=row, sticky=tk.W)
        self.exposure_entry = tk.Entry(self.prep_frame, bg=entry_bg)
        self.exposure_entry.grid(row=row, column=1, sticky=tk.W + tk.E)
        self.exposure_entry.insert(0, str(options['prep__motioncorr__dose_per_frame']))

        ###

        self.particle_frame = tk.LabelFrame(right_frame, text="Particle settings", padx=5, pady=5)
        self.particle_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.particle_frame, 1, weight=1)

        row = 0

        tk.Label(self.particle_frame, text="Symmetry:").grid(row=row, sticky=tk.W)
        self.symmetry_var = tk.StringVar()  # for data binding
        self.symmetry_entry = tk.Entry(self.particle_frame, textvariable=self.symmetry_var, bg=entry_bg)
        self.symmetry_entry.grid(row=row, column=1, sticky=tk.W)
        self.symmetry_entry.insert(0, str(options['proc__inimodel3d__sym_name']))

        row += 1

        tk.Label(self.particle_frame, text=u"Particle diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.particle_max_diam_var = tk.StringVar()  # for data binding
        self.particle_max_diam_entry = tk.Entry(self.particle_frame, textvariable=self.particle_max_diam_var, bg=entry_bg)
        self.particle_max_diam_entry.grid(row=row, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_max_diam_entry.insert(0, str(options['proc__autopick__log_diam_max']))

        row += 1
        
        tk.Label(self.particle_frame, text=u"Mask diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.mask_diameter_var = tk.StringVar()  # for data binding
        self.mask_diameter_entry = tk.Entry(self.particle_frame, textvariable=self.mask_diameter_var, bg=entry_bg)
        self.mask_diameter_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.mask_diameter_entry.insert(0, str(options['proc__class2d__particle_diameter']))
        self.mask_diameter_px = tk.Label(self.particle_frame, text="= NNN px")
        self.mask_diameter_px.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(self.particle_frame, text="Box size (px):").grid(row=row, sticky=tk.W)
        self.box_size_var = tk.StringVar()  # for data binding
        self.box_size_entry = tk.Entry(self.particle_frame, textvariable=self.box_size_var, bg=entry_bg)
        self.box_size_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.box_size_entry.insert(0, str(options['proc__extract__extract_size']))
        self.box_size_in_angstrom = tk.Label(self.particle_frame, text=u"= NNN \u212B")
        self.box_size_in_angstrom.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(self.particle_frame, text="Down-sample to (px):").grid(row=row, sticky=tk.W)
        self.extract_small_boxsize_var = tk.StringVar()  # for data binding
        self.extract_small_boxsize_entry = tk.Entry(self.particle_frame, textvariable=self.extract_small_boxsize_var, bg=entry_bg)
        self.extract_small_boxsize_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.extract_small_boxsize_entry.insert(0, str(options['proc__extract__rescale']))
        self.extract_angpix = tk.Label(self.particle_frame, text=u"= NNN \u212B/px")
        self.extract_angpix.grid(row=row, column=2,sticky=tk.W)

        row += 1

        tk.Label(self.particle_frame, text="Calculate for me:").grid(row=row, sticky=tk.W)
        self.auto_boxsize_var = tk.IntVar()
        self.auto_boxsize_button = tk.Checkbutton(self.particle_frame, var=self.auto_boxsize_var)
        self.auto_boxsize_button.grid(row=row, column=1, sticky=tk.W)
        if options['do_auto_boxsize']:
            self.auto_boxsize_button.select()

        ###

        self.proc_frame = tk.LabelFrame(right_frame, text="Processing settings", padx=5, pady=5)
        self.proc_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.proc_frame, 1, weight=1)

        row = 0

        tk.Label(self.proc_frame, text="Min. resol. micrographs (A):").grid(row=row, sticky=tk.W)
        self.minres_var = tk.StringVar()  # for data binding
        self.minres_entry = tk.Entry(self.proc_frame, textvariable=self.minres_var, bg=entry_bg)
        self.minres_entry.grid(row=row, column=1, sticky=tk.W)
        self.minres_entry.insert(0, str(options['proc__select_mics__select_maxval']))

        row += 1

        tk.Label(self.proc_frame, text="Topaz model:").grid(row=row, sticky=tk.W)
        self.topaz_model_var = tk.StringVar()  # for data binding
        self.topaz_model_entry = tk.Entry(self.proc_frame, textvariable=self.topaz_model_var, bg=entry_bg)
        self.topaz_model_entry.grid(row=row, column=1, sticky=tk.W)
        self.topaz_model_entry.insert(0, str(options['proc__topaz_model']))

        self.topaz_model_button = new_browse_button(self.proc_frame, self.topaz_model_var, filetypes=(('Topaz model file', '{*.sav}'), ('All files', '*')))
        self.topaz_model_button.grid(row=row, column=2)

        row += 1

        tk.Label(self.proc_frame, text="Min. FOM for topaz extract:").grid(row=row, sticky=tk.W)
        self.extract_topaz_thresh_var = tk.StringVar()  # for data binding
        self.extract_topaz_thresh_entry = tk.Entry(self.proc_frame, textvariable=self.extract_topaz_thresh_var, bg=entry_bg)
        self.extract_topaz_thresh_entry.grid(row=row, column=1, sticky=tk.W)
        self.extract_topaz_thresh_entry.insert(0, str(options['proc__extract__minimum_pick_fom']))

        row += 1

        tk.Label(self.proc_frame, text="LoG-pick instead of topaz?").grid(row=row, sticky=tk.W)
        self.do_log_var = tk.IntVar()
        self.do_log_button = tk.Checkbutton(self.proc_frame, var=self.do_log_var)
        self.do_log_button.grid(row=row, column=1, sticky=tk.W)
        if options['proc__do_log']:
            self.do_log_button.select()

        row += 1

        tk.Label(self.proc_frame, text=u"Ratio short/long diameter:").grid(row=row, sticky=tk.W)
        self.particle_circularity_entry = tk.Entry(self.proc_frame, bg=entry_bg)
        self.particle_circularity_entry.grid(row=row, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_circularity_entry.insert(0, str(float(options['proc__autopick__log_diam_min'])/float(options['proc__autopick__log_diam_max'])))

        row += 1

        tk.Label(self.proc_frame, text="Adjust LoG threshold:").grid(row=row, sticky=tk.W)
        self.autopick_log_thresh_var = tk.StringVar()  # for data binding
        self.autopick_log_thresh_entry = tk.Entry(self.proc_frame, textvariable=self.autopick_log_thresh_var, bg=entry_bg)
        self.autopick_log_thresh_entry.grid(row=row, column=1, sticky=tk.W)
        self.autopick_log_thresh_entry.insert(0, str(options['proc__autopick__log_adjust_thr']))

        row += 1

        tk.Label(self.proc_frame, text="Min. score for class2d select:").grid(row=row, sticky=tk.W)
        self.min_class_score_var = tk.StringVar()  # for data binding
        self.min_class_score_entry = tk.Entry(self.proc_frame, textvariable=self.min_class_score_var, bg=entry_bg)
        self.min_class_score_entry.grid(row=row, column=1, sticky=tk.W)
        self.min_class_score_entry.insert(0, str(options['proc__select_parts__rank_threshold']))

        row += 1
        
        tk.Label(self.proc_frame, text="Min. nr. parts. for 3D?").grid(row=row, sticky=tk.W)
        self.min_parts_3d_var = tk.StringVar()
        self.min_parts_3d_entry = tk.Entry(self.proc_frame, textvariable=self.min_parts_3d_var, bg=entry_bg)
        self.min_parts_3d_entry.grid(row=row, column=1, sticky=tk.W)
        self.min_parts_3d_entry.insert(0, str(options['proc__min_nr_parts_3d']))

        row += 1

        tk.Label(self.proc_frame, text="3D reference:").grid(row=row, sticky=tk.W)
        self.iniref_var = tk.StringVar()  # for data binding
        self.iniref_entry = tk.Entry(self.proc_frame, textvariable=self.iniref_var, bg=entry_bg)
        self.iniref_entry.grid(row=row, column=1, sticky=tk.W)
        self.iniref_entry.insert(0, str(options['proc__iniref']))

        ref_button = new_browse_button(self.proc_frame, self.iniref_var, filetypes=(('MRC file', '{*.mrc}'), ('All files', '*')))
        ref_button.grid(row=row, column=2)
        ###


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
                # If Nyquist freq. is better than 7.5 A, use this downscaled box, otherwise continue to next size up
                small_box_angpix = angpix * box_size_pix / small_box_pix
                if small_box_angpix < 3.75:
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

        def update_prep_status(*args_ignored, **kwargs_ignored):
            if self.get_var_as_bool(self.do_prep_var):
                for child in self.prep_frame.winfo_children():
                    child.configure(state=tk.NORMAL)
                self.mics_entry.delete(0,tk.END)
                self.mics_entry.insert(0, 'Schemes/prep/ctffind/micrographs_ctf.star')
                self.mics_entry.configure(state=tk.DISABLED)
            else:
                for child in self.prep_frame.winfo_children():
                    child.configure(state=tk.DISABLED)
                self.mics_entry.configure(state=tk.NORMAL)

        def update_proc_status(*args_ignored, **kwargs_ignored):
            if self.get_var_as_bool(self.do_proc_var):
                for child in self.particle_frame.winfo_children():
                    child.configure(state=tk.NORMAL)
                for child in self.proc_frame.winfo_children():
                    child.configure(state=tk.NORMAL)
                self.iniref_entry.configure(state=tk.NORMAL)
                self.gpu_entry.configure(state=tk.NORMAL)
                update_box_sizes()
            else:
                for child in self.particle_frame.winfo_children():
                    child.configure(state=tk.DISABLED)
                for child in self.proc_frame.winfo_children():
                    child.configure(state=tk.DISABLED)
                self.iniref_entry.configure(state=tk.DISABLED)
                self.gpu_entry.configure(state=tk.DISABLED)

        def update_autopick_status(*args_ignored, **kwargs_ignored):
            if self.get_var_as_bool(self.do_log_var):
                self.particle_circularity_entry.configure(state=tk.NORMAL)
                self.autopick_log_thresh_entry.configure(state=tk.NORMAL)
                self.extract_topaz_thresh_entry.configure(state=tk.DISABLED)
                self.topaz_model_entry.configure(state=tk.DISABLED)
            else:
                self.particle_circularity_entry.configure(state=tk.DISABLED)
                self.autopick_log_thresh_entry.configure(state=tk.DISABLED)
                self.extract_topaz_thresh_entry.configure(state=tk.NORMAL)
                self.topaz_model_entry.configure(state=tk.NORMAL)

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
        self.auto_boxsize_button.config(command=update_box_sizes)
      
        self.do_prep_button.config(command=update_prep_status)
        self.do_proc_button.config(command=update_proc_status)
        self.do_log_button.config(command=update_autopick_status)

        button_frame = tk.Frame(left_frame)
        button_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)

        self.save_button = tk.Button(button_frame, text="Save options", command=self.save_options_from_gui, bg=runbutton_bg)
        self.save_button.pack(padx=5, pady=5, side=tk.LEFT)

        self.run_button = tk.Button(button_frame, text="Save & run", command=self.run_pipeline, bg=runbutton_bg)
        self.run_button.pack(padx=5, pady=5, side=tk.LEFT)

        # Show initial pixel sizes
        update_box_sizes()
        # Show initial prep status
        update_prep_status()
        # Show initial proc status
        update_proc_status()
        # Show initial autopick status
        update_autopick_status()

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

        opts['do_prep'] = self.get_var_as_bool(self.do_prep_var)
        opts['proc__ctffind_mics'] = self.mics_entry.get()
        opts['do_proc'] = self.get_var_as_bool(self.do_proc_var)
        opts['do_auto_boxsize'] =  self.get_var_as_bool(self.auto_boxsize_var)
        opts['proc__topaz_model'] = self.topaz_model_entry.get()
        opts['prep__ctffind__do_phaseshift'] = self.get_var_as_bool(self.phaseplate_var)
        opts['proc__do_log'] = self.get_var_as_bool(self.do_log_var)

        try:
            opts['prep__importmovies__kV'] = float(self.voltage_entry.get())
        except ValueError:
            raise ValueError("Voltage must be a number")
        if opts['prep__importmovies__kV'] <= 0.0:
            warnings.append("- Voltage should be a positive number")

        try:
            opts['prep__importmovies__Cs'] = float(self.cs_entry.get())
        except ValueError:
            raise ValueError("Cs must be a number")

        try:
            opts['prep__importmovies__angpix'] = float(self.angpix_entry.get())
        except ValueError:
            raise ValueError("Pixel size must be a number")
        if opts['prep__importmovies__angpix'] <= 0.0:
            warnings.append("- Pixel size should be a positive number")

        try:
            opts['prep__motioncorr__dose_per_frame'] = float(self.exposure_entry.get())
        except ValueError:
            raise ValueError("Exposure rate must be a number")
        if opts['prep__motioncorr__dose_per_frame'] <= 0.0:
            warnings.append("- Exposure rate should be a positive number")

        if self.get_var_as_bool(self.superres_var):
            opts['prep__motioncorr__bin_factor'] = 2;
        else:
            opts['prep__motioncorr__bin_factor'] = 1;

        try:
            opts['proc__autopick__log_diam_max'] = float(self.particle_max_diam_entry.get())
            opts['proc__autopick__topaz_particle_diameter'] = float(self.particle_max_diam_entry.get())
        except ValueError:
            raise ValueError("Particle longest diameter must be a number")

        try:
            opts['proc__autopick__log_diam_min'] = float(self.particle_circularity_entry.get()) * float(self.particle_max_diam_entry.get())
        except ValueError:
            raise ValueError("Ratio short/long diameter must be a number from 0 to 1")
        if (opts['proc__autopick__log_diam_min']  < 0. or opts['proc__autopick__log_diam_min'] > opts['proc__autopick__log_diam_max']):
            warnings.append("- Ratio short/long diameter should be between 0 and 1 (thin -> circular particles)")

        try:
            opts['proc__class2d__particle_diameter'] = float(self.mask_diameter_entry.get())
            opts['proc__inimodel3d__particle_diameter'] = float(self.mask_diameter_entry.get())
            opts['proc__refine3d__particle_diameter'] = float(self.mask_diameter_entry.get())
        except ValueError:
            raise ValueError("Mask diameter must be a number")
        if opts['proc__class2d__particle_diameter'] <= 0:
            warnings.append("- Mask diameter should be a positive number")

        try:
            opts['proc__extract__extract_size'] = int(self.box_size_entry.get())
        except ValueError:
            raise ValueError("Box size must be a number")
        if opts['proc__extract__extract_size'] <= 0:
            warnings.append("- Box size should be a positive number")

        try:
            opts['proc__extract__rescale'] = int(self.extract_small_boxsize_entry.get())
            opts['proc__extract__do_rescale'] = True
        except ValueError:
            raise ValueError("Down-sampled box size must be a number")
        if opts['proc__extract__rescale'] <= 0:
            warnings.append("- Down-sampled box size should be a positive number")

        opts['prep__importmovies__fn_in_raw'] = self.import_images_entry.get()
        if opts['prep__importmovies__fn_in_raw'].startswith(('/', '..')):
            warnings.append("- Movies should be located inside the project directory")
        if '*' not in opts['prep__importmovies__fn_in_raw']:
            warnings.append("- Pattern for input movies should normally contain a '*' to select more than one file")

        opts['prep__motioncorr__fn_gain_ref'] = self.gainref_entry.get()
        if len(opts['prep__motioncorr__fn_gain_ref']) > 0 and not os.path.isfile(opts['prep__motioncorr__fn_gain_ref']) and opts['do_prep']:
            warnings.append("- Gain reference file '{}' does not exist".format(opts['prep__motioncorr__fn_gain_ref']))

        try:
            opts['proc__select_mics__select_maxval'] = int(self.minres_entry.get())
        except ValueError:
            raise ValueError("Minimum resolution for micrographs must be a number")
        if opts['proc__select_mics__select_maxval'] <= 0:
            warnings.append("- Minimum resolution for micrographs should be a positive number")

        try:
            opts['proc__autopick__log_adjust_thr'] = float(self.autopick_log_thresh_var.get())
        except ValueError:
            raise ValueError("Adjust LoG threshold must be a number")

        try:
            opts['proc__extract__minimum_pick_fom'] = float(self.extract_topaz_thresh_var.get())
        except ValueError:
            raise ValueError("Minimum FOM for topaz must be a number")

        if opts['proc__do_log']:
            opts['proc__autopick__use_gpu'] = False
            opts['proc__extract__do_fom_threshold'] = False
        else:
            opts['proc__autopick__use_gpu'] = True
            opts['proc__extract__do_fom_threshold'] = True 
            # dont want multiple topaz runs bumping into each other!
            opts['proc__autopick__nr_mpi'] = 1

        try:
            opts['proc__select_parts__rank_threshold'] = float(self.min_class_score_var.get())
        except ValueError:
            raise ValueError("Min. class2d score must be a number")

        try:
            opts['proc__min_nr_parts_3d'] = float(self.min_parts_3d_var.get())
        except ValueError:
            raise ValueError("Min. nr. parts for 3D must be a number")

        opts['proc__inimodel3d__sym_name'] = self.symmetry_entry.get()
        opts['proc__iniref'] = self.iniref_entry.get()
        opts['proc__refine3d__sym_name'] = self.symmetry_entry.get()

        # Set GPU IDs in all jobs
        opts['proc__autopick__gpu_ids'] = (self.gpu_entry.get()).split(',')[0]
        opts['proc__class2d__gpu_ids'] = self.gpu_entry.get()
        opts['proc__inimodel3d__gpu_ids'] = self.gpu_entry.get()
        opts['proc__refine3d__gpu_ids'] = (self.gpu_entry.get()).replace(',',':')

        return warnings

    def save_options_from_gui(self):
        """
        Update the full set of options from the values in the GUI, and save them to a file.

        Returns:
            True if the options were valid and saved successfully, otherwise False.
        """
        try:
            warnings = self.fetch_options_from_gui()
            if len(warnings) == 0 or tkMessageBox.askokcancel("Warning", "\n".join(warnings), icon='warning',
                                                              default=tkMessageBox.CANCEL):

                save_options(self.options)

                return True

        except Exception as ex:
            tkMessageBox.showerror("Error", ex.message)
            traceback.print_exc()
        return False

    def run_pipeline(self):
        """
        Update the full set of options from the values in the GUI, close the GUI and run the pipeline.
        """
        if self.save_options_from_gui():
            self.main_window.destroy()
            run_schemer(self.options, True) #True means launch the RELION GUI and the schemegui.py GUIs
 
def save_options(options):

    # Write the current options to a .py file
    with open(OPTIONS_FILE, 'w') as file:
        file.write("{\n") 
        for k,v in options.items():
            file.write("'%s' : '%s', \n" % (k, v))
        file.write("}\n")                  

    print(" RELION_IT: Written all options to {}".format(OPTIONS_FILE))
                
    # loop over all options and change the schemer STAR files
    for option, value in options.items():
                    
        if (value == ''):
            value = '\\"\\"'
                    
        # Set variables in scheme.star
        if option.count('__') == 1:
            splits = option.split('__')
            schemename = splits[0]
            varname = splits[1]
            schemestar = 'Schemes/' + schemename + '/scheme.star'
            if not os.path.isfile(schemestar):
                message = 'Error: ' + schemestar + ' does not exist'
                print(message)
                tkMessageBox.showerror(message)
                return False
                    
            command = 'relion_schemer --scheme ' + schemename + ' --set_var ' + varname + ' --value \"' + str(value) + '\"' + ' --original_value \"' + str(value) + '\"'
            print(' RELION_IT: executing: ', command)
            os.system(command)
                        
        # Set joboptions in job.star
        elif option.count('__') == 2:
            splits = option.split('__')
            schemename = splits[0]
            jobname = splits[1]
            joboption = splits[2]
            jobstar = 'Schemes/' + schemename + '/' + jobname + '/job.star'
            if not os.path.isfile(jobstar):
                message = 'Error: ' + jobstar + 'does not exist'
                print(message)
                tkMessageBox.showerror(message)
                        
            command = 'relion_pipeliner --editJob ' + jobstar + ' --editOption ' + joboption + ' --editValue \"' + str(value) + '\"'
            print(' RELION_IT: executing: ', command)
            os.system(command)
            
    print(' RELION_IT: done saving all options in the Schemes.') 

def run_schemer(options, do_gui):

    command = 'relion_schemer --scheme prep --reset &'
    print(' RELION_IT: executing: ', command)
    os.system(command)

    if options['do_prep']:
        
        command = 'relion_schemer --scheme prep --run --pipeline_control Schemes/prep/ >> Schemes/prep/run.out 2>> Schemes/prep/run.err &'
        print(' RELION_IT: executing: ', command)
        os.system(command)
        if do_gui:
            command = 'relion_schemegui.py prep &'
            print(' RELION_IT: executing: ', command)
            os.system(command)
          
    if options['do_proc']:

        command = 'relion_schemer --scheme proc --reset &'
        print(' RELION_IT: executing: ', command)
        os.system(command)

        command = 'relion_schemer --scheme proc --run  --pipeline_control Schemes/proc/ >> Schemes/proc/run.out 2>> Schemes/proc/run.err  &'
        print(' RELION_IT: executing: ', command)
        os.system(command)
        if do_gui:
            command = 'relion_schemegui.py proc &'
            print(' RELION_IT: executing: ', command)
            os.system(command)


    print(' RELION_IT: Now monitor the prep (and proc) Schemes from the RELION GUI ...')

    if do_gui:
        command = 'relion --do_projdir &'
        os.system(command)


def copy_scheme(schemename):
  
    ## Only copy the Scheme directory structure from the RELION installation directory if it doesn't exist yet
    if not os.path.isdir('Schemes/'+schemename):
        try:
            mydir = os.environ['RELION_SCRIPT_DIRECTORY']
        except KeyError:
            raise KeyError("Environment variable $RELION_SCRIPT_DIRECTORY has not been set. This is required to copy the prep and proc schemes from.")
        print(' RELION_IT: copying Schemes/' + schemename + ' from: ' + mydir)
        copytree(mydir+'/Schemes/'+schemename, 'Schemes/'+schemename)


def main():
    """
    Run the RELION 3.2 Schemer.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """
     # Start by parsing arguments
    # (If --help is given, the program will print a usage message and exit)
    parser = argparse.ArgumentParser()
    parser.add_argument("extra_options", nargs="*", metavar="extra_options.py",
                        help="Python files containing options for relion_it.py")
    parser.add_argument("--nogui", action="store_true", help="don't launch GUI to set options, execute non-interactively")
    parser.add_argument("--onlysave", action="store_true", help="don't launch GUI, nor execute Schemes, only save options")
    args = parser.parse_args()

    print(' RELION_IT: -------------------------------------------------------------------------------------------------------------------')
    print(' RELION_IT: script for automated, on-the-fly single-particle analysis in RELION (>= 4.0)')
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
            for k,v in user_opts.items():
                if v == 'True':
                    user_opts[k] = True
                elif v == 'False':
                    user_opts[k] = False
            opts.update(user_opts)

    # Copy Schemes over from RELION directory if they dont exit
    copy_scheme('prep')
    copy_scheme('proc')

    if args.onlysave:
        save_options(opts)
    elif args.nogui:
        save_options(opts)
        run_schemer(opts, False) #False means don't launch RELION GUI after launching the Schemes
    else:
        print(' RELION_IT: launching GUI...')
        tk_root = tk.Tk()
        tk_root.title("relion_it.py setup")
        RelionItGui(tk_root, opts)
        tk_root.mainloop()

if __name__ == "__main__":
    main()




