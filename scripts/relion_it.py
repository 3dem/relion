#!/usr/bin/env python2.7
"""
relion_it.py
============

Script for automated, on-the-fly single-particle analysis in RELION 3

Authors: Sjors H.W. Scheres, Takanori Nakane & Colin M. Palmer

Usage:
    relion_it.py  [extra_options.py [extra_options2.py ....] ] [--gui] [--continue]

To get started, go to the intended location of your RELION project directory and make sure your micrographs are
accessible from within it (e.g. in a subdirectory called `Movies/' - use a symlink if necessary). Then run this
script, providing the names of files containing options if needed. (To call the script, you'll need to enter the full
path to it, put the directory containing it on your PATH environment variable, or put a copy of the script in the
current directory.)

Run with the `--gui' option to launch a simple GUI which will set up a run from a few basic options. (The GUI can
also be used to save a complete options file that you can then edit as required.)

Once the script is running, open a normal RELION GUI to see what's happening and visualise the results.

See below for full instructions including how to handle errors. If you have any problems, please edit the script as
needed, call on your local Python expert or email the CCP-EM mailing list (https://www.jiscmail.ac.uk/ccpem).


Overview
--------

relion_it.py creates a number of RELION jobs and then runs one or more `relion_pipeliner' processes to schedule them
(exactly like using the "Schedule" button in the RELION GUI). Instructions and information are printed to the terminal
by relion_it.py as it runs.

relion_it.py uses a large number of options to control how the jobs are run. It's designed to be very flexible and so
these options can be changed in a number of ways:

- The easiest way is to use the simple GUI (enabled by passing the `--gui' argument), which allows you to set a few
  simple options. These are then used to calculate appropriate values for the complete set of options. (See "Using the
  GUI" below for more information on this.)

- For more control, options can be put into one or more Python files (with a simple "option_name = value" format or
  with more complicated calculations - see "Options files" below for more information). The names of these options
  files can passed as command line arguments to relion_it.py.

- For maximum control, you can make your own copy of this script and change the option values and the code itself
  however you want.

Before running relion_it.py, you need to make sure you're in your intended RELION project directory, and that your
movie files are accessible by relative paths within that directory (as usual for a RELION project). You could do this
by moving the files from the microscope straight into the project directory, using a symlink from your project
directory to the real location of the data, or running a script to create a new symlink to each micrograph as it is
collected.


Options files
-------------

relion_it.py uses a large number of options for controlling both the flow of the script and the parameters for
individual jobs. These options can be read from Python script files when relion_it.py is started.

The options are all listed the body of the script below, with a comment to explain each option. One way to use this
script is to copy it in its entirety into your project directory, edit the options directly in the script and then
run it (with no command line arguments). However, it's often better to keep the script in the RELION source directory
(where it can be updated easily) and use options files to configure it.

An example of a simple options file is:

angpix = 1.06

This would override the default pixel size value, but leave all other options at their defaults.

The options files are read and interpreted as Python scripts. A simple list of "option_name = value" lines is all
that is needed, though you can also use any Python commands you like to do more complex calculations. To generate
an example file containing all of the options, run "relion_it.py --gui" and then click the "Save options" button,
which will save all the current options to a file called `relion_it_options.py' in the working directory.

The options are named descriptively so you can probably understand what most of them do quite easily. For more help on
any particular option, look at the comment above its definition in this script, or search the script's code to see
how it is used.

Options files can be useful as templates. As an example, at Diamond Light Source's eBIC facility, we have a template
file called `dls_cluster_options.py' that contains the necessary settings to make relion_it.py submit most of its jobs
to run on the DLS GPU cluster. You could also set up standard templates for a particular microscope (say, voltage and
Cs settings) or for a particular project or computer configuration.

When relion_it.py starts, it reads all options files in the order they are given on the command line. Subsequent files
will override earlier ones, so the last value given for any particular option will be the value that is used.

If you start relion_it.py with the `--continue' argument, it will automatically add `relion_it_options.py' to the end
of the list of options files. This means that if you are in a project directory where the relion_it.py GUI has
previously been used, all options will be defined in the relion_it_options.py file and they will override any other
options files given on the command line. (This is very useful for restarting the script after a problem, but it would
be pointless to combine `--continue' with any options template files.)

Note that if relion_it.py finds option names that it doesn't recognise while it's reading an options file, it will
print a warning (but continue anyway). If you've been editing options files by hand, you should check the output from
relion_it.py when it starts to make sure there are no typos in the options you wanted to set. (If you're using local
variables for intermediate Python calculations in an options file, it's a good idea to use names starting with a
leading underscore so you can immediately tell them apart from warnings about genuine spelling mistakes.)


Using the GUI
-------------

The GUI provides a simple way to start new projects with relion_it.py. If you want to use it, prepare your project
directory as described above, then start the GUI with "relion_it.py --gui". (If you're using any template options
files, you can give those too, for example "relion_it.py /path/to/site/options.py --gui".)

The window that appears should be self-explanatory. Fill in the options as needed for your project, and use the check
boxes on the right to control what processing steps will be done. When you're ready, click either "Save options" or
"Save & run". The program will check the values you've entered and then use them to calculate the full set of options
for relion_it.py. The options will then be saved to a file called `relion_it_options.py', and if you clicked "Save &
run" the processing run will start immediately.

If any of the entered values are invalid (for example, if there are letters in a field which should be a number), the
GUI will display a message box with an error when you click one of the buttons. It will also display a warning if any
values appear to be incorrect (but you can choose to ignore the warning by clicking "OK").

The GUI will try to calculate a sensible set of options from the values you enter using the following steps:

1. Set the mask diameter to the next multiple of 10 above (1.1 * particle longest diameter in pixels) - for example,
   if the particle size is 180 pixels, the mask diameter will be 190 pixels.

2. Set the box size to (1.8 * particle longest diameter in pixels), rounded up to the next value in the sequence (32,
   48, 64,..., 768, 1024) - for example, if the particle size is 180 pixels, the box size will be 384 pixels.

3. If the box size is larger than 256 pixels, downscale the particles to 256-pixel boxes when they are extracted. (Or
   128 pixels in the first pass of a reference-free run.)

4a. If a 3D reference is given, use a single pass with reference-based autopicking and a batch size of 100,000
    particles.

4b. If no 3D reference is given, run a first pass with reference-free LoG autopicking and a batch size of 10,000, and
    then a second pass with reference-based autopicking and a batch size of 100,000.

These options should be sensible in many cases, but if you'd like to change them, save the options from the GUI using
the "Save options" button, close the GUI, and edit the `relion_it_options.py' file to change the option values as
needed. You can then start the processing run with "relion_it.py --continue".


Running the pipelines
---------------------

relion_it.py uses several different scheduling pipelines to run its jobs. While each one is running, a file is
created in the project directory called `RUNNING_PIPELINER_<name>'. A log of the jobs run by that pipeline is stored
in `pipeline_<name>.log'.

If you want to stop one of the pipelines for any reason, delete its `RUNNING_' file and within a minute or two the
pipeliner will notice that the file has been removed and stop.

relion_it.py itself uses a similar file called `RUNNING_RELION_IT', and you can delete this to stop the script (which
will not affect any pipelines that are already running). It keeps a list of all of the jobs it has submitted in a
file called `RELION_IT_SUBMITTED_JOBS'. This file can be edited manually if necessary (but not while the script is
running!) Most of the jobs are run by the `preprocessing' pipeline. This will do the following:

  1. Import movies
  2. Motion correction
  3. CTF estimation
  4. Particle auto-picking
  5. Particle extraction
  6. Batch selection

After a number of particles have been extracted (1,000 by default), a 2D classification job will be run to provide
feedback on the quality of the data collection and particle picking.

Particles are split into batches of a fixed size (default 10,000 for the first pass with no reference, or 100,000
otherwise). The first batch is special: as it grows, the 2D classification job is re-run repeatedly to provide early
feedback on the quality of the data. For subsequent batches, the script waits for each batch to be complete before
running 2D classification on it.

You can provide reference structures for auto-picking and 3D classification. (If you provide a 3D reference in the
GUI it will automatically be used for both tasks.)

If you do not provide a reference for auto-picking, reference-free LoG picking will be used. If you do not provide a
reference for classification, relion_it.py will run the preprocessing pipeline twice. In the first pass, an initial
model will be generated, and then a second pass of preprocessing will be done using the initial model as a reference
for auto-picking and classification.

relion_it.py makes an effort to try to identify a suitable reference to use from the classes produced by the
InitialModel job, but if it selects an inappropriate reference, you can change it by stopping the pipelines and
script ("rm RUNNING_*"), updating the reference filename stored in the file named `RELION_IT_2NDPASS_3DREF', deleting
the relevant jobs (`autopick2_job' and those following) from the `RELION_IT_SUBMITTED_JOBS' file, then restarting the
pipeline with "relion_it.py --continue".


Fixing problems
---------------

One-off job failure
```````````````````
Occasionally, a single job can fail with an isolated error, for example if there are temporary network problems while
working on a remote filesystem. If this happens, RELION will wait forever for the files to appear that would indicate
the job has finished. In the meantime, no new jobs will be run, which can cause a backlog of micrographs to build up.

To fix this (for a preprocessing job), you can just try to re-run the job from the RELION GUI. Select the job in the
"Running jobs" list, then click "Job actions" -> "Mark as finished". Select the job again in the "Finished jobs"
list, then click "Continue!" to re-start the job.

That approach should work for preprocessing jobs, but probably won't work for classification or inital model
generation jobs, since those cannot be continued and must instead be restarted from the beginning. The best way to do
that is to restart the job manually, outside the RELION GUI, and then when the job finishes RELION should continue as
if the job had never failed.

For example, with a failed local job:

    ps -e | grep relion                            # to check if the job is still active
    kill <process_id>                              # to stop the job
    # now re-run the commands from the job's `note.txt' file

or with a job that was submitted to an SGE cluster queue:

    qstat                                          # to check if the job is still active in the queue
    qdel <job_id>                                  # to remove the job from the queue
    qsub job_type/job_directory/run_submit.script  # to re-submit the job

The other option is to just run a new job from the RELION GUI in the normal way (select the job you want to "copy" in
the jobs list, make a "new" job by clicking on the job type in the list in the top-left of the GUI, then click
"Run!"). However, if you do this, relion_it.py will not know about the new job and will not run any further
downstream processing based on it. In this situation, you can either continue to process your data manually in RELION,
or you could edit the `RELION_IT_SUBMITTED_JOBS' file to replace the failed job with the manual one, and delete the
jobs that followed the original one. After that, if you re-run the script it should continue as normal from that
job onwards.


Repeated job failure
````````````````````
If a job fails repeatedly, it usually indicates that there is some problem with the job parameters or the files that
the job needs to access.

In favourable cases, it's possible you could fix the problem by selecting the job in the RELION GUI, changing one of
the parameters that is not greyed out, then clicking "Continue!". Often, though, the problem will be with one of the
parameters that can't be changed for a job that already exists, so the job will need to be deleted and recreated with
a different set of parameters.

To handle this situation, stop all of the pipelines and the relion_it.py script ("rm RUNNING_*"), then identify and
fix the problem. Often, the problem will be an error in one of the job parameters, which can usually be fixed by
changing one of the script options (for example by changing the settings in `relion_it_options.py', if you originally
used the GUI to start the run).

If the problem is caused by missing files from an upstream job, you might need to check the output of previous jobs
and look in the job directories to figure out what the problem is. Again, if it's an error in the parameters for a
job, you can probably fix it by editing `relion_it_options.py'.

After changing any script options, you'll need to use the RELION GUI to delete the affected job and all jobs
downstream of it, and also remove them from the list in the `RELION_IT_SUBMITTED_JOBS' file. Then you should be able
to restart the pipelines by running "relion_it.py --continue".

If you still can't get a particular job to run without errors, you can at least continue to run the upstream jobs
that are working properly. You can do this either by changing the options for relion_it.py (there are options to
switch off 2D or 3D classification, or to stop after CTF estimation), or by manually scheduling the jobs you want
using the RELION GUI. Remember that after running relion_it.py, you have a normal RELION project, so if the script
can't do what you want, you can simply stop it and then use all of RELION's normal job management and scheduling
abilities.


Advanced usage
--------------

It's possible to customise many aspects of the way relion_it.py works, but the details go beyond the scope of this
introduction. Simple customisation can be done by setting appropriate option values (see "Option files" above). For
more substantial changes, you might need to edit the script's Python code to get the behaviour you want. Most of the
important logic is in the `run_pipeline()' function so that's a good place to start. Good luck!

"""

import argparse
import glob
import math
import os
import runpy
import time

try:
    import Tkinter as tk
    import tkMessageBox
except ImportError:
    # The GUI is optional. If the user requests it, it will fail when it tries
    # to open so we can ignore the error for now.
    pass  

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING_RELION_IT'
SECONDPASS_REF3D_FILE = 'RELION_IT_2NDPASS_3DREF'
SETUP_CHECK_FILE = 'RELION_IT_SUBMITTED_JOBS'
PREPROCESS_SCHEDULE_PASS1 = 'PREPROCESS'
PREPROCESS_SCHEDULE_PASS2 = 'PREPROCESS_PASS2'
OPTIONS_FILE = 'relion_it_options.py'

class RelionItOptions(object):
    """
    Options for the relion_it pipeline setup script.
    
    When initialised, this contains default values for all options. Call
    ``update_from()`` to override the defaults with a dictionary of new values.
    """
    #############################################################################
    # Change the parameters below to reflect your experiment                    #
    # Current defaults reflect cryo-ARM betagal data set of RELION-3.0 tutorial #
    #############################################################################

    ### General parameters
    # Pixel size in Angstroms in the input movies
    angpix = 0.885
    # Acceleration voltage (in kV)
    voltage = 200
    # Polara = 2.0; Talos/Krios = 2.7; some Cryo-ARM = 1.4 
    Cs = 1.4


    ### Import images (Linux wild card; movies as *.mrc, *.mrcs, *.tiff or *.tif; single-frame micrographs as *.mrc)
    import_images = 'Movies/*.tiff'
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    images_are_movies = True

    
    ### MotionCorrection parameters
    # Dose in electrons per squared Angstrom per frame
    motioncor_doseperframe = 1.277
    # Gain-reference image in MRC format (only necessary if input movies are not yet gain-corrected, e.g. compressed TIFFs from K2)
    motioncor_gainreference = 'Movies/gain.mrc'


    ### CTF estimation parameters
    # Most cases won't need changes here...


    ### Autopick parameters
    # Use reference-free Laplacian-of-Gaussian picking (otherwise use reference-based template matching instead)
    autopick_do_LoG = True
    # Minimum and maximum diameter in Angstrom for the LoG filter
    autopick_LoG_diam_min = 150
    autopick_LoG_diam_max = 180
    # Use positive values (0-1) to pick fewer particles; use negative values (-1-0) to pick more particles
    autopick_LoG_adjust_threshold = 0.0
    #
    # OR:
    #
    # References for reference-based picking (when autopick_do_LoG = False)
    autopick_2dreferences = ''
    # OR: provide a 3D references for reference-based picking (when autopick_do_LoG = False)
    autopick_3dreference = ''

    # Threshold for reference-based autopicking (threshold 0 will pick too many particles. Default of 0.4 is hopefully better. Ultimately, just hope classification will sort it all out...)
    autopick_refs_threshold = 0.4
    # Minimum inter-particle distance for reference-based picking (~70% of particle diameter often works well)
    autopick_refs_min_distance = 120
    #
    # For both LoG and refs:
    #
    # Use this to remove false positives from carbon edges (useful range: 1.0-1.2, -1 to switch off)
    autopick_stddev_noise = -1
    # Use this to remove false positives from carbon edges (useful range: -0.5-0.0; -999 to switch off)
    autopick_avg_noise = -999


    ### Extract parameters
    # Box size of particles in the averaged micrographs (in pixels)
    extract_boxsize = 256
    # Down-scale the particles upon extraction?
    extract_downscale = False
    # Box size of the down-scaled particles (in pixels)
    extract_small_boxsize = 64
    # In second pass, down-scale the particles upon extraction?
    extract2_downscale = False
    # In second pass, box size of the down-scaled particles (in pixels)
    extract2_small_boxsize = 128

    
    ### Now perform 2D and/or 3D classification with the extracted particles?
    do_class2d = True
    # And/or perform 3D classification?
    do_class3d = True
    # Repeat 2D and/or 3D-classification for batches of this many particles
    batch_size = 10000
    # Number of 2D classes to use
    class2d_nr_classes  = 50
    # Diameter of the mask used for 2D/3D classification (in Angstrom)
    mask_diameter = 190
    # Symmetry group (when using SGD for initial model generation, C1 may work best)
    symmetry = 'C1'
    #
    ### 3D-classification parameters
    # Number of 3D classes to use
    class3d_nr_classes  = 4
    # Have initial 3D model? If not, calculate one using SGD initial model generation
    have_3d_reference = False
    # Initial reference model
    class3d_reference = ''
    # Is reference on correct greyscale?
    class3d_ref_is_correct_greyscale = False
    # Has the initial reference been CTF-corrected?
    class3d_ref_is_ctf_corrected = True
    # Initial lowpass filter on reference
    class3d_ini_lowpass = 40 


    ### Use the largest 3D class from the first batch as a 3D reference for a second pass of autopicking? (only when do_class3d is True)
    do_second_pass = True
    # Only move on to template-based autopicking if the 3D references achieves this resolution (in A)
    minimum_resolution_3dref_2ndpass = 20
    # In the second pass, perform 2D classification?
    do_class2d_pass2 = True
    # In the second pass, perform 3D classification?
    do_class3d_pass2 = False
    # Batch size in the second pass
    batch_size_pass2 = 100000
    

    ###################################################################################
    ############ Often the parameters below can be kept the same for a given set-up
    ###################################################################################

    ### Repeat settings for entire pipeline
    # Repeat the pre-processing runs this many times (or until RUNNING_PIPELINER_default_PREPROCESS file is deleted)
    preprocess_repeat_times = 999
    # Wait at least this many minutes between each repeat cycle
    preprocess_repeat_wait = 1
    ### Stop after CTF estimation? I.e., skip autopicking, extraction, 2D/3D classification, etc?
    stop_after_ctf_estimation = False
    # Check every this many minutes if enough particles have been extracted for a new batch of 2D-classification
    batch_repeat_time = 1


    ### MotionCorrection parameters
    # Use RELION's own implementation of motion-correction (CPU-only) instead of the UCSF implementation?
    motioncor_do_own = False
    # The number of threads (only for RELION's own implementation) is optimal when nr_movie_frames/nr_threads = integer
    motioncor_threads = 12
    # Exectutable of UCSF MotionCor2
    motioncor_exe = '/public/EM/MOTIONCOR2/MotionCor2'
    # On which GPU(s) to execute UCSF MotionCor2
    motioncor_gpu = '0'
    # How many MPI processes to use for running motion correction?
    motioncor_mpi = 1
    # Local motion-estimation patches for MotionCor2
    motioncor_patches_x = 5
    motioncor_patches_y = 5
    # B-factor in A^2 for downweighting of high-spatial frequencies
    motioncor_bfactor = 150
    # Use binning=2 for super-resolution K2 movies
    motioncor_binning = 1
    # Provide a defect file for your camera if you have one
    motioncor_defectfile = ''
    # orientation of the gain-reference w.r.t your movies (if input movies are not yet gain-corrected, e.g. TIFFs)
    motioncor_gainflip = 'No flipping (0)'
    motioncor_gainrot = 'No rotation (0)'
    # Other arguments for MotionCor2
    motioncor_other_args = ''
    # Submit motion correction job to the cluster?
    motioncor_submit_to_queue = False
    

    ### CTF estimation parameters
    # Amplitude contrast (Q0)
    ampl_contrast = 0.1
    # CTFFIND-defined parameters
    ctffind_boxsize = 512
    ctffind_astigmatism = 100
    ctffind_maxres = 5
    ctffind_minres = 30
    ctffind_defocus_max = 50000
    ctffind_defocus_min = 5000
    ctffind_defocus_step = 500
    # For Gctf: ignore parameters on the 'Searches' tab?
    ctffind_do_ignore_search_params = True
    # For Gctf: perform equi-phase averaging?
    ctffind_do_EPA = True
    # Also estimate phase shifts (for VPP data)
    ctffind_do_phaseshift = False
    # Executable to Kai Zhang's Gctf
    gctf_exe = '/public/EM/Gctf/bin/Gctf'
    # On which GPU(s) to execute Gctf
    gctf_gpu = '0'
    # Use Alexis Rohou's CTFFIND4 (CPU-only) instead?
    use_ctffind_instead = False
    # Executable for Alexis Rohou's CTFFIND4
    ctffind4_exe = '/public/EM/ctffind/ctffind.exe'
    # How many MPI processes to use for running CTF estimation?
    ctffind_mpi = 1
    # Submit CTF estimation job to the cluster?
    ctffind_submit_to_queue = False


    ### Autopick parameters
    # Use GPU-acceleration for autopicking?
    autopick_do_gpu = True
    # Which GPU(s) to use for autopicking
    autopick_gpu = '0'
    # Low-pass filter for auto-picking the micrographs
    autopick_lowpass = 20
    # Shrink factor for faster picking (0 = fastest; 1 = slowest)
    autopick_shrink_factor = 0
    # How many MPI processes to use for running auto-picking?
    autopick_mpi = 1
     # Additional arguments for autopicking
    autopick_other_args = ''
    # Submit Autopick job to the cluster?
    autopick_submit_to_queue = False
    # Are the references CTF-corrected?
    autopick_refs_are_ctf_corrected = True
    # Do the references have inverted contrast wrt the micrographs?
    autopick_refs_have_inverted_contrast = True
    # Ignore CTFs until the first peak
    autopick_refs_ignore_ctf1stpeak = False
    # Diameter of mask for the references (in A; negative value for automated detection of mask diameter)
    autopick_refs_mask_diam = -1
    # In-plane angular sampling interval
    autopick_inplane_sampling = 10
    # Symmetry of the 3D reference for autopicking
    autopick_3dref_symmetry = 'C1'
    # 3D angular sampling for generating projections of the 3D reference for autopicking (30 degrees is usually enough)
    autopick_3dref_sampling = '30 degrees'
    # Pixel size in the provided 2D/3D references (negative for same as in motion-corrected movies)
    autopick_ref_angpix = -1

    ### Extract parameters
    # Diameter for background normalisation (in pixels; negative value: default is 75% box size)
    extract_bg_diameter = -1
    # How many MPI processes to use for running particle extraction?
    extract_mpi = 1
    # Submit Extract job to the cluster?
    extract_submit_to_queue = False


    ## Discard particles based on average/stddev values? (this may be important for SGD initial model generation)
    do_discard_on_image_statistics = False
    # Discard images that have average/stddev values that are more than this many sigma away from the ensemble average
    discard_sigma = 4
    # Submit discard job to the cluster?
    discard_submit_to_queue = False


    #### Common relion_refine paremeters used for 2D/3D classification and initial model generation 
    # Read all particles in one batch into memory?
    refine_preread_images = False
    # Or copy particles to scratch disk?
    refine_scratch_disk = ''
    # Number of pooled particles?
    refine_nr_pool = 10
    # Use GPU-acceleration?
    refine_do_gpu = True
    # Which GPU to use (different from GPU used for pre-processing?)
    refine_gpu = '1'
    # How many MPI processes to use
    refine_mpi = 1
    # How many threads to use
    refine_threads = 6
    # Skip padding?
    refine_skip_padding = False
    # Submit jobs to the cluster?
    refine_submit_to_queue = False
    # Use fast subsets in 2D/3D classification when batch_size is bigger than this
    refine_batchsize_for_fast_subsets = 100000


    ### 2D classification parameters
    # Wait with the first 2D classification batch until at least this many particles are extracted
    minimum_batch_size = 1000
    # Number of iterations to perform in 2D classification
    class2d_nr_iter = 20
    # Rotational search step (in degrees)
    class2d_angle_step = 6
    # Offset search range (in pixels)
    class2d_offset_range = 5
    # Offset search step (in pixels)
    class2d_offset_step = 1
    # Option to ignore the CTFs until their first peak (try this if all particles go into very few classes) 
    class2d_ctf_ign1stpeak = False
    # Additional arguments to pass to relion-refine
    class2d_other_args = ''
    

    ### 3D classification parameters
    # Number of iterations to perform in 3D classification
    class3d_nr_iter = 20
    # Reference mask
    class3d_reference_mask = ''
    # Option to ignore the CTFs until their first peak (try this if all particles go into very few classes) 
    class3d_ctf_ign1stpeak = False
    # Regularisation parameter (T)
    class3d_T_value = 4
    # Angular sampling step
    class3d_angle_step = '7.5 degrees'
    # Offset search range (in pixels)
    class3d_offset_range = 5
    # Offset search step (in pixels)
    class3d_offset_step = 1
    # Additional arguments to pass to relion-refine
    class3d_other_args = ''


    ## SGD initial model generation
    # Number of models to generate simulatenously (K>1 may be useful for getting rid of outliers in the particle images)
    inimodel_nr_classes = 4
    # Ignore CTFs until first peak?
    inimodel_ctf_ign1stpeak = False
    # Enforce non-negative solvent?
    inimodel_solvent_flatten = True
    # Initial angular sampling
    inimodel_angle_step = '15 degrees'
    # Initial search range (in pixels)
    inimodel_offset_range = 6 
    # Initial offset search step (in pixels)
    inimodel_offset_step = 2
    # Number of initial iterations
    inimodel_nr_iter_initial = 50
    # Number of in-between iterations
    inimodel_nr_iter_inbetween = 200
    # Number of final iterations
    inimodel_nr_iter_final = 50
    # Frequency to write out information
    inimodel_freq_writeout = 10
    # Initial resolution (in A)
    inimodel_resol_ini = 35
    # Final resolution (in A)
    inimodel_resol_final = 15
    # Initial mini-batch size
    inimodel_batchsize_ini = 100
    # Final mini-batch size
    inimodel_batchsize_final = 500
    # Increased noise variance half-life (off, i.e. -1, by default; values of ~1000 have been observed to be useful in difficult cases)
    inimodel_sigmafudge_halflife = -1
    # Additional arguments to pass to relion_refine (skip annealing to get rid of outlier particles)
    inimodel_other_args = ' --sgd_skip_anneal '


    ### Cluster submission settings
    # Name of the queue to which to submit the job
    queue_name = 'openmpi'
    # Name of the command used to submit scripts to the queue
    queue_submit_command = 'qsub'
    # The template for your standard queue job submission script
    queue_submission_template = '/public/EM/RELION/relion/bin/qsub.csh'
    # Minimum number of dedicated cores that need to be requested on each node
    queue_minimum_dedicated = 1
    
    #######################################################################
    ############ typically no need to change anything below this line
    #######################################################################

    def update_from(self, other):
        """
        Update this RelionItOptions object from a dictionary.
        
        Special values (with names like '__xxx__') are removed, allowing this
        method to be given a dictionary containing the namespace from a script
        run with ``runpy``.
        """
        while len(other) > 0:
            key, value = other.popitem()
            if not (key.startswith('__') and key.endswith('__')): # exclude __name__, __builtins__ etc.
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print 'Unrecognised option {}'.format(key)
    
    def print_options(self, out_file=None):
        """
        Print the current options.
        
        This method prints the options in the same format as they are read,
        allowing options to be written to a file and re-used.

        Args:
            out_file: A file object (optional). If supplied, options will be
                written to this file, otherwise they will be printed to
                sys.stdout.
        """
        for key in sorted(dir(self)):
            if (not (key.startswith('__') and key.endswith('__'))
                and not callable(getattr(self, key))):
                print >>out_file, '{} = {}'.format(key, repr(getattr(self, key)))


class RelionItGui(object):

    def __init__(self, main_window, options):
        self.main_window = main_window
        self.options = options

        ### Create GUI

        main_frame = tk.Frame(main_window)
        main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame = tk.Frame(main_frame)
        left_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        right_frame = tk.Frame(main_frame)
        right_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        ###
        
        expt_frame = tk.LabelFrame(left_frame, text="Experimental details", padx=5, pady=5)
        expt_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        row = 0

        tk.Label(expt_frame, text="Voltage (kV):").grid(row=row, sticky=tk.W)
        self.voltage_entry = tk.Entry(expt_frame)
        self.voltage_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.voltage_entry.insert(0, str(options.voltage))

        row += 1
        
        tk.Label(expt_frame, text="Cs (mm):").grid(row=row, sticky=tk.W)
        self.cs_entry = tk.Entry(expt_frame)
        self.cs_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.cs_entry.insert(0, str(options.Cs))

        row += 1
        
        tk.Label(expt_frame, text="Phase plate?").grid(row=row, sticky=tk.W)
        self.phaseplate_var = tk.IntVar()
        phaseplate_button = tk.Checkbutton(expt_frame, var=self.phaseplate_var)
        phaseplate_button.grid(row=row, column=1, sticky=tk.W)
        if options.ctffind_do_phaseshift:
            phaseplate_button.select()

        row += 1

        tk.Label(expt_frame, text=u"Pixel size (\u212B):").grid(row=row, sticky=tk.W)
        self.angpix_entry = tk.Entry(expt_frame)
        self.angpix_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.angpix_entry.insert(0, str(options.angpix))

        row += 1
        
        tk.Label(expt_frame, text=u"Exposure rate (e\u207B / \u212B\u00B2 / frame):").grid(row=row, sticky=tk.W)
        self.exposure_entry = tk.Entry(expt_frame)
        self.exposure_entry.grid(row=row, column=1, sticky=tk.W + tk.E)
        self.exposure_entry.insert(0, str(options.motioncor_doseperframe))

        ###

        particle_frame = tk.LabelFrame(left_frame, text="Particle details", padx=5, pady=5)
        particle_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(particle_frame, 1, weight=1)

        row = 0

        tk.Label(particle_frame, text=u"Longest diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.particle_max_diam_entry = tk.Entry(particle_frame)
        self.particle_max_diam_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.particle_max_diam_entry.insert(0, str(options.autopick_LoG_diam_max))

        row += 1

        tk.Label(particle_frame, text=u"Shortest diameter (\u212B):").grid(row=row, sticky=tk.W)
        self.particle_min_diam_entry = tk.Entry(particle_frame)
        self.particle_min_diam_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.particle_min_diam_entry.insert(0, str(options.autopick_LoG_diam_min))

        row += 1
        
        tk.Label(particle_frame, text="3D reference (optional):").grid(row=row, sticky=tk.W)
        ref_3d_var = tk.StringVar()  # for data binding
        self.ref_3d_entry = tk.Entry(particle_frame, textvariable=ref_3d_var)
        self.ref_3d_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.ref_3d_entry.insert(0, str(options.autopick_3dreference))

        ###

        project_frame = tk.LabelFrame(left_frame, text="Project details", padx=5, pady=5)
        project_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(project_frame, 1, weight=1)

        row = 0

        tk.Label(project_frame, text="Project directory:").grid(row=row, sticky=tk.W)
        tk.Label(project_frame, text=os.getcwd(), anchor=tk.W).grid(row=row, column=1, sticky=tk.W)

        row += 1

        tk.Label(project_frame, text="Path to movies:").grid(row=row, sticky=tk.W)
        self.import_images_entry = tk.Entry(project_frame)
        self.import_images_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.import_images_entry.insert(0, self.options.import_images)

        row += 1
        
        tk.Label(project_frame, text="Gain reference (optional):").grid(row=row, sticky=tk.W)
        self.gainref_entry = tk.Entry(project_frame)
        self.gainref_entry.grid(row=row, column=1, sticky=tk.W+tk.E)
        self.gainref_entry.insert(0, self.options.motioncor_gainreference)

        ###

        pipeline_frame = tk.LabelFrame(right_frame, text="Pipeline control", padx=5, pady=5)
        pipeline_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        row = 0

        tk.Label(pipeline_frame, text="Stop after CTF estimation?").grid(row=row, sticky=tk.W)
        self.stop_after_ctf_var = tk.IntVar()
        stop_after_ctf_button = tk.Checkbutton(pipeline_frame, var=self.stop_after_ctf_var)
        stop_after_ctf_button.grid(row=row, column=1, sticky=tk.W)
        if options.stop_after_ctf_estimation:
            stop_after_ctf_button.select()

        row += 1

        tk.Label(pipeline_frame, text="Do 2D classification?").grid(row=row, sticky=tk.W)
        self.class2d_var = tk.IntVar()
        class2d_button = tk.Checkbutton(pipeline_frame, var=self.class2d_var)
        class2d_button.grid(row=row, column=1, sticky=tk.W)
        if options.do_class2d:
            class2d_button.select()

        row += 1

        tk.Label(pipeline_frame, text="Do 3D classification?").grid(row=row, sticky=tk.W)
        self.class3d_var = tk.IntVar()
        class3d_button = tk.Checkbutton(pipeline_frame, var=self.class3d_var)
        class3d_button.grid(row=row, column=1, sticky=tk.W)
        if options.do_class3d:
            class3d_button.select()

        row += 1

        tk.Label(pipeline_frame, text="Do second pass? (only if no 3D ref)").grid(row=row, sticky=tk.W)
        self.second_pass_var = tk.IntVar()
        second_pass_button = tk.Checkbutton(pipeline_frame, var=self.second_pass_var)
        second_pass_button.grid(row=row, column=1, sticky=tk.W)
        if options.do_second_pass:
            second_pass_button.select()

        row += 1

        tk.Label(pipeline_frame, text="Do 2D classification (2nd pass)?").grid(row=row, sticky=tk.W)
        self.class2d_pass2_var = tk.IntVar()
        class2d_pass2_button = tk.Checkbutton(pipeline_frame, var=self.class2d_pass2_var)
        class2d_pass2_button.grid(row=row, column=1, sticky=tk.W)
        class2d_pass2_button.select()
        if options.do_class2d_pass2:
            class2d_pass2_button.select()

        row += 1

        tk.Label(pipeline_frame, text="Do 3D classification (2nd pass)?").grid(row=row, sticky=tk.W)
        self.class3d_pass2_var = tk.IntVar()
        class3d_pass2_button = tk.Checkbutton(pipeline_frame, var=self.class3d_pass2_var)
        class3d_pass2_button.grid(row=row, column=1, sticky=tk.W)
        if options.do_class3d_pass2:
            class3d_pass2_button.select()

        ### Now add logic to the check boxes

        def update_gui_state(*args_ignored, **kwargs_ignored):
            new_state = tk.DISABLED if self.stop_after_ctf_var.get() else tk.NORMAL
            class2d_button.config(state=new_state)
            class3d_button.config(state=new_state)
            self.particle_max_diam_entry.config(state=new_state)
            self.particle_min_diam_entry.config(state=new_state)
            self.ref_3d_entry.config(state=new_state)
            can_do_second_pass = (self.class3d_var.get()
                                  and len(ref_3d_var.get()) == 0
                                  and not self.stop_after_ctf_var.get())
            second_pass_button.config(state=tk.NORMAL if can_do_second_pass else tk.DISABLED)
            will_do_second_pass = can_do_second_pass and self.second_pass_var.get()
            class2d_pass2_button.config(state=tk.NORMAL if will_do_second_pass else tk.DISABLED)
            class3d_pass2_button.config(state=tk.NORMAL if will_do_second_pass else tk.DISABLED)

        stop_after_ctf_button.config(command=update_gui_state)
        class3d_button.config(command=update_gui_state)
        second_pass_button.config(command=update_gui_state)
        ref_3d_var.trace('w', update_gui_state)

        ###

        button_frame = tk.Frame(right_frame)
        button_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)

        self.run_button = tk.Button(button_frame, text="Save & run", command=self.run_pipeline)
        self.run_button.pack(padx=5, pady=5, side=tk.RIGHT)

        self.save_button = tk.Button(button_frame, text="Save options", command=self.save_options)
        self.save_button.pack(padx=5, pady=5, side=tk.RIGHT)

    def fetch_options_from_gui(self):
        """
        Fetch the current values from the GUI widgets and store them in the options object.

        Returns:
            A list of warning messages about possible incorrect option values.

        Raises:
            ValueError: If an option value is invalid.
        """
        # Helper function to convert a Tk IntVar (linked to a checkbox) to a boolean value
        def get_var_as_bool(var):
            return True if var.get() == 1 else False

        opts = self.options
        warnings = []

        opts.stop_after_ctf_estimation = get_var_as_bool(self.stop_after_ctf_var)
        opts.do_class2d = get_var_as_bool(self.class2d_var)
        opts.do_class3d = get_var_as_bool(self.class3d_var)
        opts.do_second_pass = get_var_as_bool(self.second_pass_var)
        opts.do_class2d_pass2 = get_var_as_bool(self.class2d_pass2_var)
        opts.do_class3d_pass2 = get_var_as_bool(self.class3d_pass2_var)

        try:
            opts.voltage = float(self.voltage_entry.get())
        except ValueError:
            raise ValueError("Voltage must be a number")
        if opts.voltage <= 0.0:
            warnings.append("- Voltage should be a positive number")

        try:
            opts.Cs = float(self.cs_entry.get())
        except ValueError:
            raise ValueError("Cs must be a number")

        opts.ctffind_do_phaseshift = get_var_as_bool(self.phaseplate_var)

        try:
            opts.angpix = float(self.angpix_entry.get())
        except ValueError:
            raise ValueError("Pixel size must be a number")
        if opts.angpix <= 0.0:
            warnings.append("- Pixel size should be a positive number")

        try:
            opts.motioncor_doseperframe = float(self.exposure_entry.get())
        except ValueError:
            raise ValueError("Exposure rate must be a number")
        if opts.motioncor_doseperframe <= 0.0:
            warnings.append("- Exposure rate should be a positive number")

        try:
            opts.autopick_LoG_diam_max = float(self.particle_max_diam_entry.get())
        except ValueError:
            if len(self.particle_max_diam_entry.get()) == 0 and opts.stop_after_ctf_estimation:
                # This was left blank and won't be used, set to zero to avoid errors in calculations later
                opts.autopick_LoG_diam_max = 0.0
            else:
                raise ValueError("Particle longest diameter must be a number")

        try:
            opts.autopick_LoG_diam_min = float(self.particle_min_diam_entry.get())
        except ValueError:
            if len(self.particle_min_diam_entry.get()) == 0 and opts.stop_after_ctf_estimation:
                # This was left blank and won't be used, set to zero to avoid errors in calculations later
                opts.autopick_LoG_diam_min = 0.0
            else:
                raise ValueError("Particle shortest diameter must be a number")

        opts.autopick_3dreference = self.ref_3d_entry.get()
        if len(opts.autopick_3dreference) > 0 and not os.path.isfile(opts.autopick_3dreference):
            warnings.append("- 3D reference file '{}' does not exist".format(opts.autopick_3dreference))

        opts.import_images = self.import_images_entry.get()
        if opts.import_images.startswith(('/', '..')):
            warnings.append("- Movies should be located inside the project directory")

        opts.motioncor_gainreference = self.gainref_entry.get()
        if len(opts.motioncor_gainreference) > 0 and not os.path.isfile(opts.motioncor_gainreference):
            warnings.append("- Gain reference file '{}' does not exist".format(opts.motioncor_gainreference))

        return warnings

    def calculate_full_options(self):
        """
        Update the options from the values that have been fetched from the GUI.

        This method uses the values that the user has set in the GUI to calculate a number of other options for the
        script.
        """
        opts = self.options

        particle_size_pixels = opts.autopick_LoG_diam_max / opts.angpix

        # Set mask diameter to next multiple of ten above (1.1 * particle size)
        opts.mask_diameter = 10 * math.ceil(1.1 * particle_size_pixels / 10.0)

        # Calculate required box size and round up to next size in sequence
        minimum_box_size = 1.8 * particle_size_pixels
        for box in (32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024):
            if box > minimum_box_size:
                opts.extract_boxsize = box
                break
        else:
            raise ValueError("Required box size {} is too large!".format(minimum_box_size))

        # If we have a 3D reference, do a single pass with a large batch size
        if len(opts.autopick_3dreference) > 0:
            opts.autopick_do_LoG = False
            opts.autopick_refs_min_distance = opts.autopick_LoG_diam_max * 0.7
            opts.class3d_reference = opts.autopick_3dreference
            opts.do_second_pass = False
        else:
            # No 3D reference - do LoG autopicking in the first pass
            opts.autopick_do_LoG = True
            opts.class3d_reference = ''

        # Now set sensible batch size and downscaling options (leaving batch_size_pass2 at its default 100,000)
        if opts.do_second_pass:
            opts.batch_size = 10000
            if opts.extract_boxsize > 128:
                opts.extract_downscale = True
                opts.extract_small_boxsize = 128
            if opts.extract_boxsize > 256:
                opts.extract2_downscale = True
                opts.extract2_small_boxsize = 256
        else:
            opts.batch_size = 100000
            if opts.extract_boxsize > 256:
                opts.extract_downscale = True
                opts.extract_small_boxsize = 256

    def save_options(self):
        """
        Update the full set of options from the values in the GUI, and save them to a file.

        Returns:
            True if the options were valid and saved successfully, otherwise False.
        """
        try:
            warnings = self.fetch_options_from_gui()
            if len(warnings) == 0 or tkMessageBox.askokcancel("Warning", "\n".join(warnings), icon='warning'):
                self.calculate_full_options()
                print " RELION_IT: Writing all options to {}".format(OPTIONS_FILE)
                if os.path.isfile(OPTIONS_FILE):
                    print " RELION_IT: File {0} already exists; renaming old copy to {0}~".format(OPTIONS_FILE)
                    os.rename(OPTIONS_FILE, OPTIONS_FILE + '~')
                with open(OPTIONS_FILE, 'w') as optfile:
                    self.options.print_options(optfile)
                return True
        except Exception as ex:
            tkMessageBox.showerror("Error", ex.message)
        return False

    def run_pipeline(self):
        """
        Update the full set of options from the values in the GUI, close the GUI and run the pipeline.
        """
        if self.save_options():
            self.main_window.destroy()
            run_pipeline(self.options)


def safe_load_star(filename, max_try=5, wait=10, expected=[]):
    for _ in xrange(max_try):
        try:
            star = load_star(filename)
            entry = star

            # make sure the expected key is present
            for key in expected:
               entry = entry[key]
            return star
        except:
            print "safe_load_star is retrying to read: ", filename, ", expected key:", expected
            import time
            time.sleep(wait)
    assert False, "Failed to read a star file: " + filename

def load_star(filename):
    from collections import OrderedDict
    
    datasets = OrderedDict()
    current_data = None
    current_colnames = None
    
    in_loop = 0 # 0: outside 1: reading colnames 2: reading data

    for line in open(filename):
        line = line.strip()
        
        # remove comments
        comment_pos = line.find('#')
        if comment_pos > 0:
            line = line[:comment_pos]

        if line == "":
            continue

        if line.startswith("data_"):
            in_loop = 0

            data_name = line[5:]
            current_data = OrderedDict()
            datasets[data_name] = current_data

        elif line.startswith("loop_"):
            current_colnames = []
            in_loop = 1

        elif line.startswith("_"):
            if in_loop == 2:
                in_loop = 0

            elems = line[1:].split()
            if in_loop == 1:
                current_colnames.append(elems[0])
                current_data[elems[0]] = []
            else:
                current_data[elems[0]] = elems[1]

        elif in_loop > 0:
            in_loop = 2
            elems = line.split()
            assert len(elems) == len(current_colnames), ("Error in STAR file {}, number of elements in {} does not match number of column names {}"
                                                         .format(filename, elems, current_colnames))
            for idx, e in enumerate(elems):
                current_data[current_colnames[idx]].append(e)        
        
    return datasets

# Don't get stuck in infinite while True loops....
def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print " RELION_IT:", RUNNING_FILE, "file no longer exists, exiting now ..."
        exit(0)

# Allow direct progressing to the second pass
def getSecondPassReference():
    if os.path.isfile(SECONDPASS_REF3D_FILE):
        with open(SECONDPASS_REF3D_FILE, 'r') as myfile:
            filename, angpix = myfile.readlines()
    else:
        filename = ''
        angpix = '0'
    return filename.replace('\n',''), angpix.replace('\n','')

def getJobName(name_in_script, done_file):
    jobname = None
    # See if we've done this job before, i.e. whether it is in the done_file
    if (os.path.isfile(done_file)):
        f = open(done_file,'r')
        for line in f:
            elems = line.split()
            if len(elems) < 3: continue 
            if elems[0] == name_in_script:
                jobname = elems[2]
                break
        f.close()

    return jobname

def addJob(jobtype, name_in_script, done_file, options, alias=None):
    jobname = getJobName(name_in_script, done_file)

    # If we hadn't done it before, add it now
    if (jobname is not None):
        already_had_it = True 
    else:
        already_had_it = False
        optionstring = ''
        for opt in options[:]:
            optionstring += opt + ';'

        command = 'relion_pipeliner --addJob ' + jobtype + ' --addJobOptions "' + optionstring + '"'
        if alias is not None:
                command += ' --setJobAlias "' + alias + '"'
        os.system(command)

        pipeline = safe_load_star(PIPELINE_STAR, expected=['pipeline_processes', 'rlnPipeLineProcessName'])
        jobname = pipeline['pipeline_processes']['rlnPipeLineProcessName'][-1]
        
        # Now add the jobname to the done_file
        f = open(done_file,'a')
        f.write(name_in_script + ' = ' + jobname + '\n')
        f.close()

    # return the name of the job in the RELION pipeline, e.g. 'Import/job001/'
    return jobname, already_had_it

def RunJobs(jobs, repeat, wait, schedulename):

    runjobsstring = ''
    for job in jobs[:]:
        runjobsstring += job + ' '

    command = 'relion_pipeliner --schedule ' + schedulename + ' --repeat ' + str(repeat) + ' --min_wait ' + str(wait) + ' --RunJobs "' + runjobsstring + '" &' 

    os.system(command)

def WaitForJob(wait_for_this_job, seconds_wait):
    time.sleep(seconds_wait)
    print " RELION_IT: waiting for job to finish in", wait_for_this_job
    while True:
        pipeline = safe_load_star(PIPELINE_STAR, expected=['pipeline_processes', 'rlnPipeLineProcessName'])
        myjobnr = -1
        for jobnr in range(0,len(pipeline['pipeline_processes']['rlnPipeLineProcessName'])):
            jobname = pipeline['pipeline_processes']['rlnPipeLineProcessName'][jobnr]
            if jobname == wait_for_this_job:
                myjobnr = jobnr
        if myjobnr < 0:
            print " ERROR: cannot find ", wait_for_this_job, " in ", PIPELINE_STAR
            exit(1)

        status = int(pipeline['pipeline_processes']['rlnPipeLineProcessStatus'][myjobnr])
        if status == 2:
            print " RELION_IT: job in", wait_for_this_job, "has finished now"
            return
        else:
            CheckForExit()
            time.sleep(seconds_wait)

def writeManualPickingGuiFile(my_part_diam):
    if not os.path.isfile('.gui_manualpickrun.job'):
        with open('.gui_manualpickrun.job', 'w') as g:
            g.write("""job_type == 3
Pixel size (A) == -1
Black value: == 0
Blue value:  == 0
MetaDataLabel for color: == rlnParticleSelectZScore
Scale for CTF image: == 1
Particle diameter (A): == {}
Blue<>red color particles? == No
Highpass filter (A) == -1
Lowpass filter (A) == 20
Scale for micrographs: == 0.2
Red value:  == 2
Sigma contrast: == 3
White value: == 0
""".format(my_part_diam))

    return


def findBestClass(model_star_file, use_resol=True):
    model_star = safe_load_star(model_star_file)
    best_resol = 999
    best_size = 0 
    best_class = 0
    for iclass in range(0, len(model_star['model_classes']['rlnReferenceImage'])):
        mysize = float(model_star['model_classes']['rlnClassDistribution'][iclass])
        myresol = float(model_star['model_classes']['rlnEstimatedResolution'][iclass])
        if (not use_resol and (mysize > best_size or (mysize == best_size and myresol < best_resol))) \
        or (use_resol and (myresol < best_resol or (myresol == best_resol and mysize > best_size))):
            best_size = mysize
            best_class = model_star['model_classes']['rlnReferenceImage'][iclass]
            best_resol = myresol

    print " RELION_IT: found best class:",best_class,"with class size of",best_size,"and resolution of",best_resol
    return best_class, best_resol, model_star['model_general']['rlnPixelSize']

def findOutputModelStar(job_dir):
    found = None
    try:
        job_star = safe_load_star(job_dir + "job_pipeline.star", expected=['pipeline_output_edges', 'rlnPipeLineEdgeToNode'])
        for output_file in job_star["pipeline_output_edges"]['rlnPipeLineEdgeToNode']:
            if output_file.endswith("_model.star"):
                found = output_file
                break
    except:
        pass

    return found

def run_pipeline(opts):
    """
    Configure and run the RELION 3 pipeline with the given options.
    
    Args:
        opts: options for the pipeline, as a RelionItOptions object.
    """

    # if this really necessary? dont think so...
    if (os.path.isfile(PIPELINE_STAR) == False): 
        g = open(PIPELINE_STAR,'w')
        g.write('data_pipeline_general\n')
        g.write('_rlnPipeLineJobCounter 1\n')
        g.close()

    # Write RUNNING_RELION_IT file, when deleted, this script will stop
    with open(RUNNING_FILE, 'w'):
        pass

    # Write mainGUI project file, so GUI won't ask to set up a project
    with open('.gui_projectdir', 'w'):
        pass

    #### Set up GUI file for Manualpick job to allow easy viewing of autopick results
    if opts.autopick_do_LoG:
        my_part_diam = opts.autopick_LoG_diam_min
    else:
        my_part_diam = opts.autopick_refs_min_distance
    writeManualPickingGuiFile(my_part_diam)

    ### Prepare the list of queue arguments for later use
    queue_options = ['Submit to queue? == Yes',
                     'Queue name:  == {}'.format(opts.queue_name),
                     'Queue submit command: == {}'.format(opts.queue_submit_command),
                     'Standard submission script: == {}'.format(opts.queue_submission_template),
                     'Minimum dedicated cores per node: == {}'.format(opts.queue_minimum_dedicated)]

    # If we're only doing motioncorr and ctf estimation, then forget about the second pass and the batch processing
    if opts.stop_after_ctf_estimation:
        opts.do_class2d = False
        opts.do_class3d = False
        opts.do_second_pass = False

    if opts.do_second_pass:
        nr_passes = 2
    else:
        nr_passes = 1

    # if SECONDPASS_REF3D_FILE exists, go straight into the second pass
    first_pass = 0
    if opts.do_second_pass:
        secondpass_ref3d, secondpass_ref3d_angpix = getSecondPassReference()
        if not secondpass_ref3d == '':
            print ' RELION_IT: found', secondpass_ref3d,'with angpix=',secondpass_ref3d_angpix,'as a 3D reference for second pass in file',SECONDPASS_REF3D_FILE
            print ' RELION_IT: if the automatic selection of the reference turned out to be unsatisfactory,'
            print ' RELION_IT: you can re-run the second pass with another reference by:'
            print ' RELION_IT:  stopping the pipeline by deleting RUNNING_*'
            print ' RELION_IT:  updating the reference filename in',SECONDPASS_REF3D_FILE
            print ' RELION_IT:  deleting relevant jobs (autopick2_job and followings) in',SETUP_CHECK_FILE
            print ' RELION_IT:  and restarting the pipeline.'
            first_pass = 1
            opts.autopick_3dreference = secondpass_ref3d
            opts.autopick_ref_angpix = secondpass_ref3d_angpix
            opts.autopick_2dreferences = ''
            opts.autopick_do_LoG = False
            opts.class3d_reference = secondpass_ref3d
            opts.have_3d_reference = True

    # Allow to perform two passes through the entire pipeline (PREPROCESS and CLASS2D/3D batches)
    # The second pass, a 3D reference generated in the first pass will be used for template-based autopicking
    for ipass in range(first_pass, nr_passes):


        #### Set up the Import job
        import_options = ['Input files: == {}'.format(opts.import_images)]

        if opts.images_are_movies:
            import_options.append('Node type: == 2D micrograph movies (*.mrcs)')
        else:
            import_options.append('Node type: == 2D micrographs/tomograms (*.mrc)')

        import_job, already_had_it  = addJob('Import','import_job', SETUP_CHECK_FILE, import_options)

        if opts.images_are_movies:
            #### Set up the MotionCor job
            motioncorr_options = ['Input movies STAR file: == {}movies.star'.format(import_job),
                                  'MOTIONCOR2 executable: == {}'.format(opts.motioncor_exe),
                                  'Defect file: == {}'.format(opts.motioncor_defectfile),
                                  'Gain-reference image: == {}'.format(opts.motioncor_gainreference),
                                  'Gain flip: == {}'.format(opts.motioncor_gainflip),
                                  'Gain rotation: == {}'.format(opts.motioncor_gainrot),
                                  'Do dose-weighting? == Yes',
                                  'Voltage (kV): == {}'.format(opts.voltage),
                                  'Dose per frame (e/A2): == {}'.format(opts.motioncor_doseperframe),
                                  'Pixel size (A): == {}'.format(opts.angpix),
                                  'Number of patches X: ==  {}'.format(opts.motioncor_patches_x),
                                  'Number of patches Y: == {}'.format(opts.motioncor_patches_y),
                                  'Bfactor: ==  {}'.format(opts.motioncor_bfactor),
                                  'Binning factor: == {}'.format(opts.motioncor_binning),
                                  'Which GPUs to use: == {}'.format(opts.motioncor_gpu),
                                  'Other MOTIONCOR2 arguments == {}'.format(opts.motioncor_other_args),
                                  'Number of threads: == {}'.format(opts.motioncor_threads),
                                  'Number of MPI procs: == {}'.format(opts.motioncor_mpi)]

            if (opts.motioncor_do_own):
                motioncorr_options.append('Use RELION\'s own implementation? == Yes')
            else:
                motioncorr_options.append('Use RELION\'s own implementation? == No')

            if opts.motioncor_submit_to_queue:
                motioncorr_options.extend(queue_options)

            motioncorr_job, already_had_it  = addJob('MotionCorr', 'motioncorr_job', SETUP_CHECK_FILE, motioncorr_options)


        #### Set up the CtfFind job
        star_name = 'corrected_micrographs.star' if opts.images_are_movies else 'micrographs.star'
        ctffind_options = ['Input micrographs STAR file: == {}{}'.format(motioncorr_job, star_name),
                           'Voltage (kV): == {}'.format(opts.voltage),
                           'Spherical aberration (mm): == {}'.format(opts.Cs),
                           'Amplitude contrast: == {}'.format(opts.ampl_contrast),
                           'Amount of astigmatism (A): == {}'.format(opts.ctffind_astigmatism),
                           'FFT box size (pix): == {}'.format(opts.ctffind_boxsize),
                           'Maximum defocus value (A): == {}'.format(opts.ctffind_defocus_max),
                           'Minimum defocus value (A): == {}'.format(opts.ctffind_defocus_min),
                           'Defocus step size (A): == {}'.format(opts.ctffind_defocus_step),
                           'Magnified pixel size (Angstrom): == {}'.format(opts.angpix * opts.motioncor_binning),
                           'Maximum resolution (A): == {}'.format(opts.ctffind_maxres),
                           'Minimum resolution (A): == {}'.format(opts.ctffind_minres),
                           'Gctf executable: == {}'.format(opts.gctf_exe),
                           'Which GPUs to use: == {}'.format(opts.gctf_gpu),
                           'CTFFIND-4.1 executable: == {}'.format(opts.ctffind4_exe),
                           'Number of MPI procs: == {}'.format(opts.ctffind_mpi)]

        if opts.use_ctffind_instead:
            ctffind_options.append('Use CTFFIND-4.1? == Yes')
            ctffind_options.append('Use Gctf instead? == No')
        else:
            ctffind_options.append('Use CTFFIND-4.1? == No')
            ctffind_options.append('Use Gctf instead? == Yes')
            if (opts.ctffind_do_ignore_search_params):
                ctffind_options.append('Ignore \'Searches\' parameters? == Yes')
            else:
                ctffind_options.append('Ignore \'Searches\' parameters? == No')
            if (opts.ctffind_do_EPA):
                 ctffind_options.append('Perform equi-phase averaging? == Yes')
            else:
                 ctffind_options.append('Perform equi-phase averaging? == No')

        if opts.ctffind_do_phaseshift:
            ctffind_options.append('Estimate phase shifts? == Yes')
        else:
            ctffind_options.append('Estimate phase shifts? == No')

        if opts.ctffind_submit_to_queue:
            ctffind_options.extend(queue_options)

        ctffind_job, already_had_it  = addJob('CtfFind', 'ctffind_job', SETUP_CHECK_FILE, ctffind_options)

        runjobs = [import_job, motioncorr_job, ctffind_job]

        # There is an option to stop on-the-fly processing after CTF estimation
        if not opts.stop_after_ctf_estimation:
            autopick_options = ['Input micrographs for autopick: == {}micrographs_ctf.star'.format(ctffind_job),
                                'Min. diameter for LoG filter (A) == {}'.format(opts.autopick_LoG_diam_min),
                                'Max. diameter for LoG filter (A) == {}'.format(opts.autopick_LoG_diam_max),
                                'Maximum resolution to consider (A) == {}'.format(opts.autopick_lowpass),
                                'Adjust default threshold == {}'.format(opts.autopick_LoG_adjust_threshold),
                                '2D references: == {}'.format(opts.autopick_2dreferences),
                                '3D reference: == {}'.format(opts.autopick_3dreference),
                                'Symmetry: == {}'.format(opts.autopick_3dref_symmetry),
                                'Pixel size in references (A) == {}'.format(opts.autopick_ref_angpix),
                                '3D angular sampling: == {}'.format(opts.autopick_3dref_sampling),
                                'In-plane angular sampling (deg) == {}'.format(opts.autopick_inplane_sampling),
                                'Picking threshold: == {}'.format(opts.autopick_refs_threshold),
                                'Minimum inter-particle distance (A): == {}'.format(opts.autopick_refs_min_distance),
                                'Mask diameter (A) == {}'.format(opts.autopick_refs_mask_diam),
                                'Maximum stddev noise: == {}'.format(opts.autopick_stddev_noise),
                                'Minimum avg noise: == {}'.format(opts.autopick_avg_noise),
                                'Shrink factor: == {}'.format(opts.autopick_shrink_factor),
                                'Which GPUs to use: == {}'.format(opts.autopick_gpu),
                                'Additional arguments: == {}'.format(opts.autopick_other_args),
                                'Number of MPI procs: == {}'.format(opts.autopick_mpi)]

            if not opts.autopick_3dreference == '':
                autopick_options.append('OR: provide a 3D reference? == Yes')
            else:
                autopick_options.append('OR: provide a 3D reference? == No')

            if opts.autopick_do_LoG:
                autopick_options.append('OR: use Laplacian-of-Gaussian? == Yes')
            else:
                autopick_options.append('OR: use Laplacian-of-Gaussian? == No')

            if opts.autopick_refs_are_ctf_corrected:
                autopick_options.append('Are References CTF corrected? == Yes')
            else:
                autopick_options.append('Are References CTF corrected? == No')

            if opts.autopick_refs_have_inverted_contrast:
                autopick_options.append('References have inverted contrast? == Yes')
            else:
                autopick_options.append('References have inverted contrast? == No')

            if opts.autopick_refs_ignore_ctf1stpeak:
                autopick_options.append('Ignore CTFs until first peak? == Yes')
            else:
                autopick_options.append('Ignore CTFs until first peak? == No')

            if opts.autopick_do_gpu and (not opts.autopick_do_LoG):
                autopick_options.append('Use GPU acceleration? == Yes')
            else:
                autopick_options.append('Use GPU acceleration? == No')

            if opts.autopick_submit_to_queue:
                autopick_options.extend(queue_options)
            
            if ipass == 0:
                autopick_job_name = 'autopick_job'
                autopick_alias = 'pass 1'
            else:
                autopick_job_name = 'autopick2_job'
                autopick_alias = 'pass 2'

            autopick_job, already_had_it  = addJob('AutoPick', autopick_job_name, SETUP_CHECK_FILE, autopick_options, alias=autopick_alias)
            runjobs.append(autopick_job)

            #### Set up the Extract job
            extract_options = ['Input coordinates:  == {}coords_suffix_autopick.star'.format(autopick_job),
                               'micrograph STAR file:  == {}micrographs_ctf.star'.format(ctffind_job),
                               'Diameter background circle (pix):  == {}'.format(opts.extract_bg_diameter),
                               'Particle box size (pix): == {}'.format(opts.extract_boxsize),
                               'Number of MPI procs: == {}'.format(opts.extract_mpi)]

            if ipass == 0:
                if opts.extract_downscale:
                    extract_options.append('Rescale particles? == Yes')
                    extract_options.append('Re-scaled size (pixels):  == {}'.format(opts.extract_small_boxsize))
            else:
                if opts.extract2_downscale:
                    extract_options.append('Rescale particles? == Yes')
                    extract_options.append('Re-scaled size (pixels):  == {}'.format(opts.extract2_small_boxsize))

            if opts.extract_submit_to_queue:
                extract_options.extend(queue_options)

            if ipass == 0:
                extract_job_name = 'extract_job'
                extract_alias = 'pass 1'
            else:
                extract_job_name = 'extract2_job'
                extract_alias = 'pass 2'

            extract_job, already_had_it  = addJob('Extract', extract_job_name, SETUP_CHECK_FILE, extract_options, alias=extract_alias)
            runjobs.append(extract_job)

            if (ipass == 0 and (opts.do_class2d or opts.do_class3d)) or (ipass == 1 and (opts.do_class2d_pass2 or opts.do_class3d_pass2)):
                #### Set up the Select job to split the particle STAR file into batches
                split_options = ['OR select from particles.star: == {}particles.star'.format(extract_job),
                                 'OR: split into subsets? == Yes',
                                 'OR: number of subsets:  == -1']

                if ipass == 0:
                    split_job_name = 'split_job'
                    split_options.append('Subset size:  == {}'.format(opts.batch_size))
                    split_alias = 'into {}'.format(opts.batch_size) 
                else:
                    split_job_name = 'split2_job'
                    split_options.append('Subset size:  == {}'.format(opts.batch_size_pass2))
                    split_alias = 'into {}'.format(opts.batch_size_pass2)

                split_job, already_had_it = addJob('Select', split_job_name, SETUP_CHECK_FILE, split_options, alias=split_alias)

                # Now start running stuff
                runjobs.append(split_job)

        # Now execute the entire preprocessing pipeliner
        if ipass == 0:
            preprocess_schedule_name = PREPROCESS_SCHEDULE_PASS1
        else:
            preprocess_schedule_name = PREPROCESS_SCHEDULE_PASS2
        RunJobs(runjobs, opts.preprocess_repeat_times, opts.preprocess_repeat_wait, preprocess_schedule_name)
        print ' RELION_IT: submitted',preprocess_schedule_name,'pipeliner with', opts.preprocess_repeat_times,'repeats of the preprocessing jobs'
        print ' RELION_IT: this pipeliner will run in the background of your shell. You can stop it by deleting the file RUNNING_PIPELINER_'+preprocess_schedule_name


        ########## From now on, process extracted particles in batches for 2D or 3D classification, only perform SGD inimodel for first batch and if no 3D reference is available

        # There is again an option to stop here...
        if (ipass == 0 and (opts.do_class2d or opts.do_class3d)) or (ipass == 1 and (opts.do_class2d_pass2 or opts.do_class3d_pass2)):

            ### If necessary, rescale the 3D reference in the second pass!
            if ipass == 1 and (opts.extract_downscale or opts.extract2_downscale):
                particles_angpix = opts.angpix
                if opts.images_are_movies:
                    particles_angpix = particles_angpix * opts.motioncor_binning
                if opts.extract2_downscale:
                    particles_angpix = particles_angpix * opts.extract_boxsize / opts.extract2_small_boxsize
                    particles_boxsize = opts.extract2_small_boxsize
                else:
                    particles_boxsize = opts.extract_boxsize
                if abs(float(particles_angpix) - float(opts.autopick_ref_angpix)) > 0.01:
                    # Now rescale the reference for 3D classification
                    opts.class3d_reference = opts.autopick_3dreference.replace('.mrc','_rescaled.mrc')
                    print ' RELION_IT: rescaling the 3D reference from pixel size',opts.autopick_ref_angpix,'to',particles_angpix,'and saving the new reference as',opts.class3d_reference
                    command = 'relion_image_handler --i ' + opts.autopick_3dreference + ' --o ' + opts.class3d_reference + ' --angpix ' + str(opts.autopick_ref_angpix) + ' --rescale_angpix ' + str(particles_angpix) + ' --new_box ' + str(particles_boxsize) 
                    os.system(command)


            print ' RELION_IT: now entering an infinite loop for batch-processing of particles. You can stop this loop by deleting the file', RUNNING_FILE
            
            # It could be that this is a restart, so check previous_batch1_size in the output directory.
            # Also check the presence of class2d_job_batch_001 in case the first job was not submitted yet.
            if getJobName("class2d_job_batch_001", SETUP_CHECK_FILE) is not None and \
               os.path.isfile(split_job + 'particles_split001.star'):
                batch1 = safe_load_star(split_job + 'particles_split001.star', expected=['', 'rlnMicrographName'])
                previous_batch1_size = len(batch1['']['rlnMicrographName'])
            else:
                previous_batch1_size = 0

            continue_this_pass = True
            while continue_this_pass:

                have_new_batch = False
                nr_batches = len(glob.glob(split_job + "particles_split*.star"))
                for ibatch  in range(0, nr_batches):
                    iibatch = ibatch + 1
                    batch_name = split_job + 'particles_split%03d.star' % iibatch

                    batch = safe_load_star(batch_name, expected=['', 'rlnMicrographName'])
                    batch_size = len(batch['']['rlnMicrographName'])
                    rerun_batch1 = False
                    if ( iibatch == 1 and batch_size > previous_batch1_size and batch_size > opts.minimum_batch_size ):
                        previous_batch1_size = batch_size
                        rerun_batch1 = True

                    particles_star_file = batch_name

                    # The first batch is special: perform 2D classification with smaller batch size (but at least minimum_batch_size) and keep overwriting in the same output directory
                    if ( rerun_batch1 or batch_size == opts.batch_size):


                        # Discard particles with odd average/stddev values
                        if opts.do_discard_on_image_statistics:

                            #### Run a Select job to get rid of particles with outlier average/stddev values...
                            discard_options = ['OR select from particles.star: == {}'.format(batch_name),
                                               'OR: select on image statistics? == Yes',
                                               'Sigma-value for discarding images: == {}'.format(opts.discard_sigma),
                                               'Metadata label for images: == rlnImageName']

                            if ipass == 0:
                                discard_job_name = 'discard_job'
                            else:
                                discard_job_name = 'discard2_job'

                            if opts.discard_submit_to_queue:
                                discard_options.extend(queue_options)

                            discard_job, already_had_it = addJob('Select', discard_job_name, SETUP_CHECK_FILE, discard_options)

                            if ((not already_had_it) or rerun_batch1):
                                have_new_batch = True
                                RunJobs([discard_job], 1, 1, 'DISCARD')
                                print " RELION_IT: submitted job to discard based on image statistics for", batch_size ,"particles in", batch_name

                                # Wait here until this Discard job is finished. Check every thirty seconds
                                WaitForJob(discard_job, 30)

                            particles_star_file = discard_job + 'particles.star'


                        # 2D classification
                        if (ipass == 0 and opts.do_class2d) or (ipass == 1 and opts.do_class2d_pass2):

                            class2d_options = ['Input images STAR file: == {}'.format(particles_star_file),
                                               'Number of classes: == {}'.format(opts.class2d_nr_classes),
                                               'Mask diameter (A): == {}'.format(opts.mask_diameter),
                                               'Number of iterations: == {}'.format(opts.class2d_nr_iter),
                                               'Angular search range - psi (deg): == {}'.format(opts.class2d_angle_step), 
                                               'Offset search range (pix): == {}'.format(opts.class2d_offset_range),
                                               'Offset search step (pix): == {}'.format(opts.class2d_offset_step),
                                               'Number of pooled particles: == {}'.format(opts.refine_nr_pool),
                                               'Which GPUs to use: == {}'.format(opts.refine_gpu),
                                               'Number of MPI procs: == {}'.format(opts.refine_mpi),
                                               'Number of threads: == {}'.format(opts.refine_threads),
                                               'Copy particles to scratch directory: == {}'.format(opts.refine_scratch_disk),
                                               'Additional arguments: == {}'.format(opts.class2d_other_args)]

                            if batch_size > opts.refine_batchsize_for_fast_subsets:
                                class2d_options.append('Use fast subsets (for large data sets)? == Yes')
                            else:
                                class2d_options.append('Use fast subsets (for large data sets)? == No')

                            if opts.refine_do_gpu:
                                class2d_options.append('Use GPU acceleration? == Yes')
                            else:
                                class2d_options.append('Use GPU acceleration? == No')

                            if opts.class2d_ctf_ign1stpeak:
                                class2d_options.append('Ignore CTFs until first peak? == Yes')
                            else:
                                class2d_options.append('Ignore CTFs until first peak? == No')

                            if opts.refine_preread_images:
                                class2d_options.append('Pre-read all particles into RAM? == Yes')
                            else:
                                class2d_options.append('Pre-read all particles into RAM? == No')

                            if opts.refine_submit_to_queue:
                                class2d_options.extend(queue_options)

                            if ipass == 0:
                                jobname = 'class2d_job_batch_{:03d}'.format(iibatch)
                                alias = 'pass1_batch_{:03d}'.format(iibatch)
                            else:
                                jobname = 'class2d_pass2_job_batch_{:03d}'.format(iibatch)
                                alias = 'pass2_batch_{:03d}'.format(iibatch)

                            class2d_job, already_had_it = addJob('Class2D', jobname, SETUP_CHECK_FILE, class2d_options, alias=alias)
 
                            if ((not already_had_it) or rerun_batch1):
                                have_new_batch = True
                                RunJobs([class2d_job], 1, 1, 'CLASS2D')
                                print " RELION_IT: submitted 2D classification with", batch_size ,"particles in", class2d_job

                                # Wait here until this Class2D job is finished. Check every thirty seconds
                                WaitForJob(class2d_job, 30)

                    # Perform 3D classification
                    if (ipass == 0 and opts.do_class3d) or (ipass == 1 and opts.do_class3d_pass2):

                        # Do SGD initial model generation only in the first pass, when no reference is provided AND only for the first (complete) batch, for subsequent batches use that model
                        if (not opts.have_3d_reference) and ipass == 0 and iibatch == 1 and batch_size == opts.batch_size:

                            inimodel_options = ['Input images STAR file: == {}'.format(particles_star_file),
                                                'Symmetry: == {}'.format(opts.symmetry),
                                                'Mask diameter (A): == {}'.format(opts.mask_diameter),
                                                'Number of classes: == {}'.format(opts.inimodel_nr_classes),
                                                'Initial angular sampling: == {}'.format(opts.inimodel_angle_step), 
                                                'Offset search range (pix): == {}'.format(opts.inimodel_offset_range),
                                                'Offset search step (pix): == {}'.format(opts.inimodel_offset_step),
                                                'Number of initial iterations: == {}'.format(opts.inimodel_nr_iter_initial), 
                                                'Number of in-between iterations: == {}'.format(opts.inimodel_nr_iter_inbetween), 
                                                'Number of final iterations: == {}'.format(opts.inimodel_nr_iter_final),
                                                'Write-out frequency (iter): == {}'.format(opts.inimodel_freq_writeout), 
                                                'Initial resolution (A): == {}'.format(opts.inimodel_resol_ini), 
                                                'Final resolution (A): == {}'.format(opts.inimodel_resol_final), 
                                                'Initial mini-batch size: == {}'.format(opts.inimodel_batchsize_ini), 
                                                'Final mini-batch size: == {}'.format(opts.inimodel_batchsize_final), 
                                                'SGD increased noise variance half-life: == {}'.format(opts.inimodel_sigmafudge_halflife), 
                                                'Number of pooled particles: == 1',
                                                'Which GPUs to use: == {}'.format(opts.refine_gpu),
                                                'Number of MPI procs: == {}'.format(opts.refine_mpi),
                                                'Number of threads: == {}'.format(opts.refine_threads),
                                                'Copy particles to scratch directory: == {}'.format(opts.refine_scratch_disk),
                                                'Additional arguments: == {}'.format(opts.inimodel_other_args)]

                            if opts.inimodel_solvent_flatten:
                                inimodel_options.append('Flatten and enforce non-negative solvent? == Yes')
                            else:    
                                inimodel_options.append('Flatten and enforce non-negative solvent? == No')

                            if opts.refine_skip_padding:
                                inimodel_options.append('Skip padding? == Yes')
                            else:    
                                inimodel_options.append('Skip padding? == No')

                            if opts.refine_do_gpu:
                                inimodel_options.append('Use GPU acceleration? == Yes')
                            else:
                                inimodel_options.append('Use GPU acceleration? == No')

                            if opts.inimodel_ctf_ign1stpeak:
                                inimodel_options.append('Ignore CTFs until first peak? == Yes')
                            else:
                                inimodel_options.append('Ignore CTFs until first peak? == No')

                            if opts.refine_preread_images:
                                inimodel_options.append('Pre-read all particles into RAM? == Yes')
                            else:
                                inimodel_options.append('Pre-read all particles into RAM? == No')

                            if opts.refine_submit_to_queue:
                                inimodel_options.extend(queue_options)

                            inimodel_job, already_had_it = addJob('InitialModel', 'inimodel', SETUP_CHECK_FILE, inimodel_options)              

                            if (not already_had_it):
                                have_new_batch = True
                                RunJobs([inimodel_job], 1, 1, 'INIMODEL')
                                print " RELION_IT: submitted initial model generation with", batch_size ,"particles in", inimodel_job

                                # Wait here until this inimodel job is finished. Check every thirty seconds
                                WaitForJob(inimodel_job, 30)

                            sgd_model_star = findOutputModelStar(inimodel_job)
                            if sgd_model_star is None:
                                print " RELION_IT: Initial model generation " + inimodel_job + " does not contain expected output maps."
                                print " RELION_IT: This job should have finished, but you may continue it from the GUI. "
                                raise Exception("ERROR!! quitting the pipeline.") # TODO: MAKE MORE ROBUST

                            # Use the model of the largest class for the 3D classification below
                            total_iter = opts.inimodel_nr_iter_initial + opts.inimodel_nr_iter_inbetween + opts.inimodel_nr_iter_final
                            best_inimodel_class, best_inimodel_resol, best_inimodel_angpix = findBestClass(sgd_model_star, use_resol=True)
                            opts.class3d_reference = best_inimodel_class
                            opts.class3d_ref_is_correct_greyscale = True
                            opts.class3d_ref_is_ctf_corrected = True
                            opts.have_3d_reference = True


                        if opts.have_3d_reference:
                            # Now perform the actual 3D classification
                            class3d_options = ['Input images STAR file: == {}'.format(particles_star_file),
                                               'Reference map: == {}'.format(opts.class3d_reference),
                                               'Initial low-pass filter (A): == {}'.format(opts.class3d_ini_lowpass),
                                               'Symmetry: == {}'.format(opts.symmetry),
                                               'Regularisation parameter T: == {}'.format(opts.class3d_T_value),
                                               'Reference mask (optional): == {}'.format(opts.class3d_reference_mask),
                                               'Number of classes: == {}'.format(opts.class3d_nr_classes),
                                               'Mask diameter (A): == {}'.format(opts.mask_diameter),
                                               'Number of iterations: == {}'.format(opts.class3d_nr_iter),
                                               'Angular sampling interval: == {}'.format(opts.class3d_angle_step), 
                                               'Offset search range (pix): == {}'.format(opts.class3d_offset_range),
                                               'Offset search step (pix): == {}'.format(opts.class3d_offset_step),
                                               'Number of pooled particles: == {}'.format(opts.refine_nr_pool),
                                               'Which GPUs to use: == {}'.format(opts.refine_gpu),
                                               'Number of MPI procs: == {}'.format(opts.refine_mpi),
                                               'Number of threads: == {}'.format(opts.refine_threads),
                                               'Copy particles to scratch directory: == {}'.format(opts.refine_scratch_disk),
                                               'Additional arguments: == {}'.format(opts.class3d_other_args)]

                            if batch_size > opts.refine_batchsize_for_fast_subsets:
                                class3d_options.append('Use fast subsets (for large data sets)? == Yes')
                            else:
                                class3d_options.append('Use fast subsets (for large data sets)? == No')
                            
                            if opts.class3d_ref_is_correct_greyscale:
                                class3d_options.append('Ref. map is on absolute greyscale? == Yes')
                            else:
                                class3d_options.append('Ref. map is on absolute greyscale? == No')

                            if opts.class3d_ref_is_ctf_corrected:
                                class3d_options.append('Has reference been CTF-corrected? == Yes')
                            else:
                                class3d_options.append('Has reference been CTF-corrected? == No')

                            if opts.refine_skip_padding:
                                class3d_options.append('Skip padding? == Yes')
                            else:    
                                class3d_options.append('Skip padding? == No')

                            if opts.refine_do_gpu:
                                class3d_options.append('Use GPU acceleration? == Yes')
                            else:
                                class3d_options.append('Use GPU acceleration? == No')

                            if opts.class3d_ctf_ign1stpeak:
                                class3d_options.append('Ignore CTFs until first peak? == Yes')
                            else:
                                class3d_options.append('Ignore CTFs until first peak? == No')

                            if opts.refine_preread_images:
                                class3d_options.append('Pre-read all particles into RAM? == Yes')
                            else:
                                class3d_options.append('Pre-read all particles into RAM? == No')

                            if opts.refine_submit_to_queue:
                                class3d_options.extend(queue_options)

                            if ipass == 0:
                                jobname = 'class3d_job_batch_{:03d}'.format(iibatch)
                                alias = 'pass1_batch_{:03d}'.format(iibatch)
                            else:
                                jobname = 'class3d2_job_batch_{:03d}'.format(iibatch)
                                alias = 'pass2_batch_{:03d}'.format(iibatch)

                            class3d_job, already_had_it = addJob('Class3D', jobname, SETUP_CHECK_FILE, class3d_options, alias=alias) 

                            if ((not already_had_it) or rerun_batch1):
                                have_new_batch = True
                                RunJobs([class3d_job], 1, 1, 'CLASS3D')
                                print ' RELION_IT: submitted 3D classification with', batch_size ,'particles in', class3d_job

                                # Wait here until this Class2D job is finished. Check every thirty seconds
                                WaitForJob(class3d_job, 30)

                            class3d_model_star = findOutputModelStar(class3d_job)
                            if class3d_model_star is None:
                                print " RELION_IT: 3D Classification " + class3d_job + " does not contain expected output maps."
                                print " RELION_IT: This job should have finished, but you may continue it from the GUI."
                                raise Exception("ERROR!! quitting the pipeline.") # TODO: MAKE MORE ROBUST

                            best_class3d_class, best_class3d_resol, best_class3d_angpix = findBestClass(class3d_model_star, use_resol=True)

                            # Once the first batch in the first pass is completed: move on to the second pass
                            if (ipass == 0 and opts.do_second_pass and iibatch == 1 and best_class3d_resol < opts.minimum_resolution_3dref_2ndpass):
                                opts.autopick_3dreference = best_class3d_class
                                opts.autopick_ref_angpix = best_class3d_angpix
                                opts.autopick_2dreferences = ''
                                opts.autopick_do_LoG = False
                                opts.class3d_reference = best_class3d_class
                                opts.have_3d_reference = True
                                opts.autopick_3dref_symmetry = opts.symmetry

                                # Stop the PREPROCESS pipeliner of the first pass by removing its RUNNING file
                                filename_to_remove = 'RUNNING_PIPELINER_'+preprocess_schedule_name
                                if os.path.isfile(filename_to_remove):
                                    print ' RELION_IT: removing file',filename_to_remove,'to stop the pipeliner from the first pass'
                                    os.remove(filename_to_remove)

                                # Generate a file to indicate we're in the second pass, so that restarts of the python script will be smooth
                                g = open(SECONDPASS_REF3D_FILE,'w')
                                g.write(str(best_class3d_class)+'\n'+str(best_class3d_angpix)+'\n')
                                g.close()

                                # Move out of this ipass of the passes loop....
                                ibatch = nr_batches+1
                                continue_this_pass = False
                                print ' RELION_IT: moving on to the second pass using',opts.autopick_3dreference,'for template-based autopicking'
                                # break out of the for-loop over the batches
                                break


                if not have_new_batch:
                    CheckForExit()
                    # The following prevents checking the particles.star file too often
                    time.sleep(60*opts.batch_repeat_time)


def main():
    """
    Run the RELION 3 pipeline.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """
    # Start by parsing arguments
    # (If --help is given, the program will print a usage message and exit)
    parser = argparse.ArgumentParser()
    parser.add_argument("extra_options", nargs="*", metavar="extra_options.py",
                        help="Python files containing options for relion_it.py")
    parser.add_argument("--gui", action="store_true", help="launch a simple GUI to set options")
    parser.add_argument("--continue", action="store_true", dest="continue_",
                        help="continue a previous run by loading options from ./relion_it_options.py")
    args = parser.parse_args()

    print ' RELION_IT: -------------------------------------------------------------------------------------------------------------------'
    print ' RELION_IT: script for automated, on-the-fly single-particle analysis in RELION (>= 3.0-alpha-5)'
    print ' RELION_IT: authors: Sjors H.W. Scheres, Takanori Nakane & Colin M. Palmer'
    print ' RELION_IT: '
    print ' RELION_IT: usage: ./relion_it.py [extra_options.py [extra_options2.py ....] ] [--gui] [--continue]'
    print ' RELION_IT: '
    print ' RELION_IT: this script will check whether processes are still running using files with names starting with RUNNING' 
    print ' RELION_IT:   you can restart this script after stopping previous processes by deleting all RUNNING files'
    print ' RELION_IT: this script keeps track of already submitted jobs in a filed called',SETUP_CHECK_FILE
    print ' RELION_IT:   upon a restart, jobs present in this file will be continued (for preprocessing), or ignored when already finished'
    print ' RELION_IT: if you would like to re-do a specific job from scratch (e.g. because you changed its parameters)' 
    print ' RELION_IT:   remove that job, and those that depend on it, from the',SETUP_CHECK_FILE
    print ' RELION_IT: -------------------------------------------------------------------------------------------------------------------'
    print ' RELION_IT: '
    
    # Make sure no other version of this script are running...
    if os.path.isfile(RUNNING_FILE):
        print " RELION_IT: ERROR:", RUNNING_FILE, "is already present: delete this file and make sure no other copy of this script is running. Exiting now ..."
        exit(0)

    # Also make sure the preprocessing pipeliners are stopped before re-starting this script
    for checkfile in ('RUNNING_PIPELINER_'+PREPROCESS_SCHEDULE_PASS1, 'RUNNING_PIPELINER_'+PREPROCESS_SCHEDULE_PASS2):
        if os.path.isfile(checkfile):
            print " RELION_IT: ERROR:", checkfile, "is already present: delete this file and make sure no relion_pipeliner job is still running. Exiting now ..."
            exit(0)

    if args.continue_:
        print ' RELION_IT: continuing a previous run. Options will be loaded from ./relion_it_options.py'
        args.extra_options.append(OPTIONS_FILE)

    opts = RelionItOptions()
    for user_opt_file in args.extra_options:
        print ' RELION_IT: reading options from {}'.format(user_opt_file)
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)

    if args.gui:
        print ' RELION_IT: launching GUI...'
        tk_root = tk.Tk()
        tk_root.title("relion_it.py setup")
        RelionItGui(tk_root, opts)
        tk_root.mainloop()
    else:
        run_pipeline(opts)


if __name__ == "__main__":
    main()




