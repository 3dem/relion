#!/usr/bin/env python2.7
"""
relion_it
---------

Pipeline setup script for live processing with RELION 3.

Authors: Sjors H.W. Scheres, Takanori Nakane & Colin Palmer

Call this from the intended location of the RELION project directory, and provide
the name of a file containing options if needed. See the relion_it_options.py
file for an example.

Usage:
    /path/to/relion_it.py [options_file...]
"""

import collections
import os
import runpy
import sys
import time
import glob

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING_RELION_IT'

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
    voltage = 300
    # Polara = 2.0; Talos/Krios = 2.7; Cryo-ARM = 1.4 
    Cs = 1.4


    ### Import images (Linux wild card; movies as *.mrcs or *.tif; single-frame micrographs as *.mrc)
    import_images = 'Movies/*tiff'
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    images_are_movies = True

    
    ### MotionCorrection parameters
    # The number of threads (only for RELION's own implementation) is optimal when nr_movie_frames/nr_threads = integer
    motioncor_threads = 12
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
    autopick_references = ''
    # Threshold for reference-based autopicking
    autopick_refs_threshold = 0.4
    # Minimum inter-particle distance for reference-based picking
    autopick_refs_min_distance = 100
    #
    # For both LoG and refs:
    #
    # Use this to remove false positives from carbon edges (useful range: 1.0-1.2, -1 to switch off)
    autopick_stddev_noise = -1
    # Use this to remove false positives from carbon edges (useful range: -0.5-0.0; -999 to switch off)
    autopick_avg_noise = -999


    ### Extract parameters
    # Box size of particles in the averaged micrographs (in pixels)
    extract_boxsize = 260
    # Down-scale the particles upon extraction?
    extract_downscale = False
    # Box size of the down-scaled particles (in pixels)
    extract_small_boxsize = 128
 
    
    ### Now perform 2D and/or 3D classification with the extracted particles?
    do_class2d = False
    do_class3d = True
    # Repeat 2D and/or 3D-classification for batches of this many particles
    class_batch_size = 10000
    # Wait with the first 2D and/or 3D-classification batch until at least this many particles are extracted
    class_min_batch_size = 500
    # Diameter of the mask used for 2D/3D classification (in Angstrom)
    class_mask_diameter = 190
    #
    ### 2D-classification parameters
    # Number of 2D classes to use
    class2d_nr_classes  = 50
    #
    ### 3D-classification parameters
    # Number of 3D classes to use
    class3d_nr_classes  = 4
    # Initial reference model
    class3d_reference = 'myref.mrc'
    # Is reference on correct greyscale?
    class3d_ref_is_correct_greyscale = False
    # Has the initial reference been CTF-corrected?
    class3d_ref_is_ctf_corrected = True
    # Initial lowpass filter on reference
    class3d_ini_lowpass = 40
    # Symmetry group
    class3d_symmetry = 'D2'
    
    
    ###################################################################################
    ############ Often the parameters below can be kept the same for a given set-up
    ###################################################################################

    ### Repeat settings for entire pipeline
    # Repeat the pre-processing runs this many times (or until RUNNING_PIPELINER_default_PREPROCESS file is deleted)
    preprocess_repeat_times = 999
    # Wait at least this many minutes between each repeat cycle
    preprocess_repeat_wait = 1
    ### Stop after CTF estimation? I.e., skip autopicking, extraction, etc?
    stop_after_ctf_estimation = False
    # Check every this many minutes if enough particles have been extracted for a new batch of 2D-classification
    class2d_repeat_time = 1


    ### MotionCorrection parameters
    # Use RELION's own implementation of motion-correction (CPU-only) instead of the UCSF implementation?
    motioncor_do_own = True
    # Exectutable of UCSF MotionCor2
    motioncor_exe = '/public/EM/MOTIONCOR2/MotionCor2-1.0.0'
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


    ### Extract parameters
    # Diameter for background normalisation (in pixels; negative value: default is 75% box size)
    extract_bg_diameter = -1
    # How many MPI processes to use for running particle extraction?
    extract_mpi = 1
    # Submit Extract job to the cluster?
    extract_submit_to_queue = False


    #### 2D/3D classification parameters
    # Read all particles in one batch into memory?
    class_preread_images = False
    # Or copy particles to scratch disk?
    class_scratch_disk = ''
    # Number of iterations to perform in 2D/3D classification
    class_nr_iter = 20
    #
    ### 2D classification parameters
    # Rotational search step (in degrees)
    class2d_angle_step = 6
    # Offset search range (in pixels)
    class2d_offset_range = 5
    # Offset search step (in pixels)
    class2d_offset_step = 1
    # Option to ignore the CTFs until their first peak (try this if all particles go into very few classes) 
    class2d_ctf_ign1stpeak = False
    # Use GPU-acceleration for 2D classification?
    class2d_do_gpu = True
    # Which GPU to use for 2D-classification (different from GPU used for pre-processing?)
    class2d_gpu = '1'
    # How many MPI processes to use for 2D classification
    class2d_mpi = 1
    # How many threads to use for 2D classification
    class2d_threads = 6
    # Additional arguments to pass to relion-refine
    class2d_other_args = ''
    # Submit 2D classification job to the cluster?
    class2d_submit_to_queue = False
    #
    ### 3D classification parameters
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
    # Skip padding in 3D classification?
    class3d_skip_padding = False
    # Use GPU-acceleration for 3D classification?
    class3d_do_gpu = True
    # Which GPU to use for 3D-classification (different from GPU used for pre-processing?)
    class3d_gpu = '1'
    # How many MPI processes to use for 2D classification
    class3d_mpi = 1
    # How many threads to use for 2D classification
    class3d_threads = 6
    # Additional arguments to pass to relion-refine
    class3d_other_args = ''
    # Submit 2D classification job to the cluster?
    class3d_submit_to_queue = False


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

    setup_check_file = 'RELION_IT_DONE_PREPROCESS_SETUP'
    
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
            assert len(elems) == len(current_colnames)
            for idx, e in enumerate(elems):
                current_data[current_colnames[idx]].append(e)        
        
    return datasets

# Don't get stuck in infinite while True loops....
def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print " RELION_IT:", RUNNING_FILE, "file no longer exists, exiting now ..."
        exit(0)

def addJob(jobtype, name_in_script, done_file, options):

    jobname = ""
    # See if we've done this job before, i.e. whether it is in the done_file
    if (os.path.isfile(done_file)):
        f = open(done_file,'r')
        for line in f: 
            if name_in_script in line:
                jobname = line.split()[2]
        f.close()
        
    # If we hadn't done it before, add it now
    if (jobname != ""):
        already_had_it = True 
    else:
        already_had_it = False
        optionstring = ''
        for opt in options[:]:
            optionstring += opt + ';'

        command = 'relion_pipeliner --addJob ' + jobtype + ' --addJobOptions "' + optionstring + '"'

        os.system(command)

        pipeline = load_star(PIPELINE_STAR)
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
        pipeline = load_star(PIPELINE_STAR)
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
        

def run_pipeline(opts):
    """
    Configure and run the RELION 3 pipeline with the given options.
    
    Args:
        opts: options for the pipeline, as a RelionItOptions object.
    """
    # Only setup the jobs in the pipeline the first time

    if (os.path.isfile(PIPELINE_STAR) == False): 
        g = open(PIPELINE_STAR,'w')
        g.write('data_pipeline_general')
        g.write('_rlnPipeLineJobCounter 1')
        g.close()

    # Write RUNNING_RELION_IT file, when deleted, this script will stop
    with open(RUNNING_FILE, 'w'):
        pass

    ### Prepare the list of queue arguments for later use
    queue_options = ['Submit to queue? == Yes',
                     'Queue name:  == {}'.format(opts.queue_name),
                     'Queue submit command: == {}'.format(opts.queue_submit_command),
                     'Standard submission script: == {}'.format(opts.queue_submission_template),
                     'Minimum dedicated cores per node: == {}'.format(opts.queue_minimum_dedicated)]

    #### Set up the Import job
    import_options = ['Input files: == {}'.format(opts.import_images)]

    if opts.images_are_movies:
        import_options.append('Node type: == 2D micrograph movies (*.mrcs)')
    else:
        import_options.append('Node type: == 2D micrographs/tomograms (*.mrc)')
                
    import_job, already_had_it  = addJob('Import','import_job', opts.setup_check_file, import_options)

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
                              'Number of threads: == {}'.format(opts.motioncor_threads),
                              'Number of MPI procs: == {}'.format(opts.motioncor_mpi)]

        if (opts.motioncor_do_own):
            motioncorr_options.append('Use RELION\'s own implementation? == Yes')
        else:
            motioncorr_options.append('Use RELION\'s own implementation? == No')

        if opts.motioncor_submit_to_queue:
            motioncorr_options.extend(queue_options)
        
	motioncorr_job, already_had_it  = addJob('MotionCorr', 'motioncorr_job', opts.setup_check_file, motioncorr_options)


    #### Set up the CtfFind job
    star_name = 'corrected_micrographs.star' if opts.images_are_movies else 'micrographs.star'
    ctffind_options = ['Input micrographs STAR file: == {}{}'.format(motioncorr_job, star_name),
                       'Voltage (kV): == {}'.format(opts.voltage),
                       'Spherical aberration (mm): == {}'.format(opts.Cs),
                       'Amplitude contrast: == {}'.format(opts.ampl_contrast),
                       'Amount of astigmatism (A): == {}'.format(opts.ctffind_astigmatism),
                       'FFT box size (pix): == {}'.format(opts.ctffind_boxsize),
                       'Perform equi-phase averaging? == Yes', 
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

    if opts.ctffind_do_phaseshift:
        ctffind_options.append('Estimate phase shifts? == Yes')
    else:
        ctffind_options.append('Estimate phase shifts? == No')


    ctffind_job, already_had_it  = addJob('CtfFind', 'ctffind_job', opts.setup_check_file, ctffind_options)

    runjobs = [import_job, motioncorr_job, ctffind_job]

    if opts.stop_after_ctf_estimation:
        opts.do_class2d = False
        opts.do_class3d = False
    else:
        #### Set up the Autopick job
        autopick_options = ['Input micrographs for autopick: == {}micrographs_ctf.star'.format(ctffind_job),
                            'Min. diameter for LoG filter (A) == {}'.format(opts.autopick_LoG_diam_min),
                            'Max. diameter for LoG filter (A) == {}'.format(opts.autopick_LoG_diam_max),
                            'Maximum resolution to consider (A) == {}'.format(opts.autopick_lowpass),
                            'Adjust default threshold == {}'.format(opts.autopick_LoG_adjust_threshold),
                            'References: == {}'.format(opts.autopick_references),
                            'Picking threshold: == {}'.format(opts.autopick_refs_threshold),
                            'Minimum inter-particle distance (A): == {}'.format(opts.autopick_refs_min_distance),
                            'Mask diameter (A) == {}'.format(opts.autopick_refs_mask_diam),
                            'Maximum stddev noise: == {}'.format(opts.autopick_stddev_noise),
                            'Minimum avg noise: == {}'.format(opts.autopick_avg_noise),
                            'Shrink factor: == {}'.format(opts.autopick_shrink_factor),
                            'Which GPUs to use: == {}'.format(opts.autopick_gpu),
                            'Additional arguments: == {}'.format(opts.autopick_other_args),
                            'Number of MPI procs: == {}'.format(opts.autopick_mpi)]

        if opts.autopick_do_LoG:
            autopick_options.append('Or use Laplacian-of-Gaussian? == Yes')
        else:
            autopick_options.append('Or use Laplacian-of-Gaussian? == No')
        
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

        if opts.autopick_do_gpu:
            autopick_options.append('Use GPU acceleration? == Yes')
        else:
            autopick_options.append('Use GPU acceleration? == No')


        if opts.autopick_submit_to_queue:
            autopick_options.extend(queue_options)
        
        autopick_job, already_had_it  = addJob('AutoPick', 'autopick_job', opts.setup_check_file, autopick_options)
        runjobs.append(autopick_job)

        #### Set up GUI file for Manualpick job to allow easy viewing of autopick results
        if opts.autopick_do_LoG:
            my_part_diam = opts.autopick_LoG_diam_min
        else:
            my_part_diam = opts.autopick_refs_min_distance

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

        #### Set up the Extract job
        extract_options = ['Input coordinates:  == {}coords_suffix_autopick.star'.format(autopick_job),
                           'micrograph STAR file:  == {}micrographs_ctf.star'.format(ctffind_job),
                           'Particle box size (pix): == {}'.format(opts.extract_boxsize),
                           'Number of MPI procs: == {}'.format(opts.extract_mpi)]

        if opts.extract_downscale:
            extract_options.append('Rescale particles? == Yes')
            extract_options.append('Re-scaled size (pixels):  == {}'.format(opts.extract_small_boxsize))
        
        if opts.extract_submit_to_queue:
            extract_options.extend(queue_options)

        extract_job, already_had_it  = addJob('Extract', 'extract_job', opts.setup_check_file, extract_options)
        runjobs.append(extract_job)
    
        if opts.do_class2d or opts.do_class3d:
            #### Set up the Select job to split the particle STAR file into batches
            select_options = ['OR select from particles.star: == {}particles.star'.format(extract_job),
                              'Split into subsets? == Yes',
                              'Subset size:  == {}'.format(opts.class_batch_size),
                              'OR: number of subsets:  == -1']

            select_job, already_had_it = addJob('Select', 'select_job', opts.setup_check_file, select_options)

            # Now start running stuff
            runjobs.append(select_job)

    # Write mainGUI project file, so GUI won't ask to set up a project
    with open('.gui_projectdir', 'w'):
        pass

    # Now execute the entire preprocessing pipeliner
    RunJobs(runjobs, opts.preprocess_repeat_times, opts.preprocess_repeat_wait, 'PREPROCESS')
    print " RELION_IT: submitted PREPROCESS pipeliner with", opts.preprocess_repeat_times, "repeats of the preprocessing jobs"

    # Possibility to stop here...
    if opts.do_class2d or opts.do_class3d:

        previous_batch1_size = 0
        while True:
		    
            have_new_batch = False
            nr_batches = len(glob.glob(select_job + "particles_split*.star"))
            for ibatch  in range(0, nr_batches):
                iibatch = ibatch + 1
                batch_name = select_job + 'particles_split%03d.star' % iibatch

                batch = load_star(batch_name)
                batch_size = len(batch['']['rlnMicrographName'])
                
                # The first batch is special: process with smaller size
                if ( ( iibatch == 1 and batch_size > opts.class_min_batch_size ) or batch_size == opts.class_batch_size):

                    rerun_batch1 = False
                    if ( iibatch == 1 and batch_size > previous_batch1_size ):
                        previous_batch1_size = batch_size
                        rerun_batch1 = True

                    if opts.do_class2d:

                        class2d_options = ['Input images STAR file: == {}'.format(batch_name),
                                           'Number of classes: == {}'.format(opts.class2d_nr_classes),
                                           'Mask diameter (A): == {}'.format(opts.class_mask_diameter),
                                           'Number of iterations: == {}'.format(opts.class_nr_iter),
                                           'Angular search range - psi (deg): == {}'.format(opts.class2d_angle_step), 
                                           'Offset search range (pix): == {}'.format(opts.class2d_offset_range),
                                           'Offset search step (pix): == {}'.format(opts.class2d_offset_step),
                                           'Number of pooled particles: == 30',
                                           'Which GPUs to use: == {}'.format(opts.class2d_gpu),
                                           'Number of MPI procs: == {}'.format(opts.class2d_mpi),
                                           'Number of threads: == {}'.format(opts.class2d_threads),
                                           'Copy particles to scratch directory: == {}'.format(opts.class_scratch_disk),
                                           'Additional arguments: == {}'.format(opts.class2d_other_args)]
                
                        if opts.class2d_do_gpu:
                            class2d_options.append('Use GPU acceleration? == Yes')
                        else:
                            class2d_options.append('Use GPU acceleration? == No')
                
                        if opts.class2d_ctf_ign1stpeak:
                            class2d_options.append('Ignore CTFs until first peak? == Yes')
                        else:
                            class2d_options.append('Ignore CTFs until first peak? == No')
                
                        if opts.class_preread_images:
                            class2d_options.append('Pre-read all particles into RAM? == Yes')
                        else:
                            class2d_options.append('Pre-read all particles into RAM? == No')
                
                        if opts.class2d_submit_to_queue:
                            class2d_options.extend(queue_options)
                
                        jobname = 'class2d_job_batch_{:03d}'.format(iibatch)
                        class2d_job, already_had_it = addJob('Class2D', jobname, opts.setup_check_file, class2d_options)              


                        if ((not already_had_it) or rerun_batch1):
                            have_new_batch = True
                            RunJobs([class2d_job], 1, 1, 'CLASS2D')
                            print " RELION_IT: submitted 2D classification with", batch_size ,"particles in", class2d_job

                            # Wait here until this Class2D job is finished. Check every thirty seconds
                            WaitForJob(class2d_job, 30)


                    if opts.do_class3d:

                        class3d_options = ['Input images STAR file: == {}'.format(batch_name),
                                           'Reference map: == {}'.format(opts.class3d_reference),
                                           'Initial low-pass filter (A): == {}'.format(opts.class3d_ini_lowpass),
                                           'Symmetry: == {}'.format(opts.class3d_symmetry),
                                           'Regularisation parameter T: == {}'.format(opts.class3d_T_value),
                                           'Reference mask (optional): == {}'.format(opts.class3d_reference_mask),
                                           'Number of classes: == {}'.format(opts.class3d_nr_classes),
                                           'Mask diameter (A): == {}'.format(opts.class_mask_diameter),
                                           'Number of iterations: == {}'.format(opts.class_nr_iter),
                                           'Angular sampling interval: == {}'.format(opts.class3d_angle_step), 
                                           'Offset search range (pix): == {}'.format(opts.class3d_offset_range),
                                           'Offset search step (pix): == {}'.format(opts.class3d_offset_step),
                                           'Number of pooled particles: == 30',
                                           'Which GPUs to use: == {}'.format(opts.class3d_gpu),
                                           'Number of MPI procs: == {}'.format(opts.class3d_mpi),
                                           'Number of threads: == {}'.format(opts.class3d_threads),
                                           'Copy particles to scratch directory: == {}'.format(opts.class_scratch_disk),
                                           'Additional arguments: == {}'.format(opts.class3d_other_args)]
                
                        if opts.class3d_ref_is_correct_greyscale:
                            class3d_options.append('Ref. map is on absolute greyscale? == Yes')
                        else:
                            class3d_options.append('Ref. map is on absolute greyscale? == No')

                        if opts.class3d_ref_is_ctf_corrected:
                            class3d_options.append('Has reference been CTF-corrected? == Yes')
                        else:
                            class3d_options.append('Has reference been CTF-corrected? == No')

                        if opts.class3d_skip_padding:
                            class3d_options.append('Skip padding? == Yes')
                        else:    
                            class3d_options.append('Skip padding? == No')

                        if opts.class3d_do_gpu:
                            class3d_options.append('Use GPU acceleration? == Yes')
                        else:
                            class3d_options.append('Use GPU acceleration? == No')
                
                        if opts.class3d_ctf_ign1stpeak:
                            class3d_options.append('Ignore CTFs until first peak? == Yes')
                        else:
                            class3d_options.append('Ignore CTFs until first peak? == No')
                
                        if opts.class_preread_images:
                            class3d_options.append('Pre-read all particles into RAM? == Yes')
                        else:
                            class3d_options.append('Pre-read all particles into RAM? == No')
                
                        if opts.class3d_submit_to_queue:
                            class3d_options.extend(queue_options)
                
                        jobname = 'class3d_job_batch_{:03d}'.format(iibatch)
                        class3d_job, already_had_it = addJob('Class3D', jobname, opts.setup_check_file, class3d_options)              

                        if ((not already_had_it) or rerun_batch1):
                            have_new_batch = True
                            RunJobs([class3d_job], 1, 1, 'CLASS3D')
                            print " RELION_IT: submitted 3D classification with", batch_size ,"particles in", class3d_job

                            # Wait here until this Class2D job is finished. Check every thirty seconds
                            WaitForJob(class3d_job, 30)


            if not have_new_batch:
                CheckForExit()
                # The following prevents checking the particles.star file too often
                time.sleep(60*opts.class2d_repeat_time)


def main():
    """
    Run the RELION 3 pipeline.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """
    opts = RelionItOptions()
    for user_opt_file in sys.argv[1:]:
        print ' RELION_IT: reading options from {}'.format(user_opt_file)
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)
    run_pipeline(opts)


if __name__ == "__main__":
    main()

