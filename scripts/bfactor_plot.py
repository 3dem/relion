#!/usr/bin/env python
"""
bfactor_plot
---------

Pipeline setup script for automated processing with RELION 4.

Authors: Sjors H.W. Scheres, Takanori Nakane & Colin Palmer

Call this from the intended location of the RELION project directory, and provide
the name of a file containing options if needed. See the relion_it_options.py
file for an example.

Usage:
    /path/to/relion_it.py [options_file...]
"""

from __future__ import print_function

import collections
import os
import runpy
import sys
import time
import glob
from math import log, sqrt

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING' # prefix is appended in main()
SETUP_CHECK_FILE = 'SUBMITTED_JOBS' # prefix is appended in main()

class RelionItOptions(object):
    """
    Options for the relion_it pipeline setup script.
    
    When initialised, this contains default values for all options. Call
    ``update_from()`` to override the defaults with a dictionary of new values.
    """
    #############################################################################
    # Change the parameters below to reflect your experiment                    #
    #############################################################################

    # job prefix
    prefix = 'BFACTOR_PLOT_'

    # If program crahses saying "'utf-8' codec can't decode byte 0xXX in position YY",
    # most likely run.job file in the job directory contains garbage bytes.

    # Refine3D job with all particles
    # This must be a job from RELION 3.1, not 3.0.
    input_refine3d_job = 'Refine3D/job040/'
    # PostProcess job for resolution assessment
    input_postprocess_job = 'PostProcess/job083/'
    # Minimum number of particles
    minimum_nr_particles = 100
    # Maximum number of particles
    maximum_nr_particles = 9999999

    #### relion_refine paremeters 
    # Initial low-pass filter for the refinements
    refine_ini_lowpass = 40
    # Read all particles in one batch into memory?
    refine_preread_images = False
    # Or copy particles to scratch disk?
    refine_scratch_disk = ''
    # Number of pooled particles?
    refine_nr_pool = 10
    # Use GPU-acceleration?
    refine_do_gpu = True
    # Which GPU to use (different from GPU used for pre-processing?)
    refine_gpu = ''
    # How many MPI processes to use
    refine_mpi = 5
    # How many threads to use
    refine_threads = 6
    # Skip padding?
    refine_skip_padding = False
    # Submit jobs to the cluster?
    refine_submit_to_queue = False

    ### Cluster submission settings
    # Name of the queue to which to submit the job
    queue_name = 'openmpi'
    # Name of the command used to submit scripts to the queue
    queue_submit_command = 'qsub -l gpu=4'
    # The template for your standard queue job submission script
    queue_submission_template = '/public/EM/RELION/relion/bin/qsub.csh'
    # Minimum number of dedicated cores that need to be requested on each node
    queue_minimum_dedicated = 32
    
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
                    print('Unrecognised option {}'.format(key))

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
            if in_loop == 2:
                in_loop = 0
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

def addJob(jobtype, name_in_script, done_file, options, template=None, alias=None):
    jobname = getJobName(name_in_script, done_file)

    # If we hadn't done it before, add it now
    if (jobname is not None):
        already_had_it = True 
    else:
        already_had_it = False
        optionstring = ''
        for opt in options[:]:
            optionstring += opt + ';'

        command = 'relion_pipeliner'

        if template is None:
            command += ' --addJob ' + jobtype
        else:
            command += ' --addJobFromStar ' + template

        command += ' --addJobOptions "' + optionstring + '"'
        if alias is not None:
            command += ' --setJobAlias "' + alias + '"'

        #print("Debug: addJob executes " + command)
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

    #print("Debug: RunJobs executes " + command)
    os.system(command)

def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print(" RELION_IT:", RUNNING_FILE, "file no longer exists, exiting now ...")
        exit(0)

def WaitForJob(wait_for_this_job, seconds_wait):
    time.sleep(seconds_wait)
    print(" RELION_IT: waiting for job to finish in", wait_for_this_job)
    while True:
        pipeline = load_star(PIPELINE_STAR)
        myjobnr = -1
        for jobnr in range(0,len(pipeline['pipeline_processes']['rlnPipeLineProcessName'])):
            jobname = pipeline['pipeline_processes']['rlnPipeLineProcessName'][jobnr]
            if jobname == wait_for_this_job:
                myjobnr = jobnr
        if myjobnr < 0:
            print(" ERROR: cannot find ", wait_for_this_job, " in ", PIPELINE_STAR)
            exit(1)

        status = pipeline['pipeline_processes']['rlnPipeLineProcessStatusLabel'][myjobnr]
        if status == "Succeeded":
            print(" RELION_IT: job in", wait_for_this_job, "has finished now")
            return
        else:
            CheckForExit()
            time.sleep(seconds_wait)

def find_split_job_output(prefix, n, max_digits=6):
    import os.path
    for i in range(max_digits):
        filename = prefix + str(n).rjust(i, '0') + '.star'
        if os.path.isfile(filename):
            return filename
    return None

def line_fit(xs, ys):
    n = len(xs)
    assert n == len(ys)

    mean_x = 0.0
    mean_y = 0.0
    for x, y in zip(xs, ys):
        mean_x += x
        mean_y += y

    mean_x /= n
    mean_y /= n

    var_x = 0.0
    cov_xy = 0.0
    for x, y in zip(xs, ys):
        var_x += (x - mean_x) ** 2
        cov_xy += (x - mean_x) * (y - mean_y)

    slope = cov_xy / var_x
    intercept = mean_y - slope * mean_x

    return slope, intercept

def get_postprocess_result(post_star):
    result = load_star(post_star)['general']
    resolution = float(result['rlnFinalResolution'])
    pp_bfactor = float(result['rlnBfactorUsedForSharpening'])
    return resolution, pp_bfactor

def run_pipeline(opts):
    """
    Configure and run the RELION 3 pipeline with the given options.
    
    Args:
        opts: options for the pipeline, as a RelionItOptions object.
    """

    # Write RUNNING_RELION_IT file, when deleted, this script will stop
    with open(RUNNING_FILE, 'w'):
        pass

    ### Prepare the list of queue arguments for later use
    queue_options = ['Submit to queue? == Yes',
                     'Queue name:  == {}'.format(opts.queue_name),
                     'Queue submit command: == {}'.format(opts.queue_submit_command),
                     'Standard submission script: == {}'.format(opts.queue_submission_template),
                     'Minimum dedicated cores per node: == {}'.format(opts.queue_minimum_dedicated)]

    # Get the original STAR file
    refine3d_run_file = opts.input_refine3d_job+'job.star'
    all_particles_star_file = None
    if os.path.exists(refine3d_run_file):
        for line in open(refine3d_run_file,'r'):
            if 'fn_img' in line:
                all_particles_star_file = line.split()[1].replace('\n','')
                break
    else:
        refine3d_run_file = opts.input_refine3d_job+'run.job' # old style
        for line in open(refine3d_run_file,'r'):
            if 'Input images STAR file' in line:
                all_particles_star_file = line.split(' == ')[1].replace('\n','')
                break
    if all_particles_star_file is None:
        print(' ERROR: cannot find input STAR file in', refine3d_run_file)
        exit(1)

    all_particles = load_star(all_particles_star_file)
    all_nr_particles = len(all_particles['particles']['rlnImageName'])
    all_particles_resolution, all_particles_bfactor = get_postprocess_result(opts.input_postprocess_job + 'postprocess.star') 

    nr_particles = []
    resolutions = []
    pp_bfactors = []

    current_nr_particles = opts.minimum_nr_particles
    while current_nr_particles <= opts.maximum_nr_particles and current_nr_particles < all_nr_particles:

        schedule_name = 'batch_' + str(current_nr_particles)

        # A. Split the STAR file
        split_options = ['OR select from particles.star: == {}'.format(all_particles_star_file),
                         'OR: split into subsets? == Yes',
                         'Subset size:  == {}'.format(current_nr_particles),
                         'Randomise order before making subsets?: == Yes',
                         'OR: number of subsets:  == 1']

        split_job_name = 'split_job_' + str(current_nr_particles)
        split_alias = opts.prefix + 'split_' + str(current_nr_particles)
        split_job, already_had_it = addJob('Select', split_job_name, SETUP_CHECK_FILE, split_options, None, split_alias)
        if not already_had_it:
            RunJobs([split_job], 1, 0, schedule_name)
            WaitForJob(split_job, 30)

        # B. Run Refine3D
        split_filename = find_split_job_output('{}particles_split'.format(split_job), 1)
        assert split_filename is not None
        refine_options = ['Input images STAR file: == {}'.format(split_filename),
                          'Number of pooled particles: == {}'.format(opts.refine_nr_pool),
                          'Which GPUs to use: == {}'.format(opts.refine_gpu),
                          'Number of MPI procs: == {}'.format(opts.refine_mpi),
                          'Initial low-pass filter (A): == {}'.format(opts.refine_ini_lowpass),
                          'Number of threads: == {}'.format(opts.refine_threads)]

        if opts.refine_skip_padding:
            refine_options.append('Skip padding? == Yes')
        else:    
            refine_options.append('Skip padding? == No')

        if opts.refine_do_gpu:
            refine_options.append('Use GPU acceleration? == Yes')
        else:
            refine_options.append('Use GPU acceleration? == No')

        if opts.refine_preread_images:
            refine_options.append('Pre-read all particles into RAM? == Yes')
            refine_options.append('Copy particles to scratch directory: == ')
        else:
            refine_options.append('Pre-read all particles into RAM? == No')
            refine_options.append('Copy particles to scratch directory: == {}'.format(opts.refine_scratch_disk))
        
        if opts.refine_submit_to_queue:
            refine_options.extend(queue_options)
        else:
            refine_options.append('Submit to queue? == No')

        refine_job_name = 'refine_job_' + str(current_nr_particles)
        refine_alias = opts.prefix + str(current_nr_particles)
        refine_job, already_had_it = addJob('Refine3D', refine_job_name, SETUP_CHECK_FILE, refine_options, refine3d_run_file, refine_alias)
        if not already_had_it:
            RunJobs([refine_job], 1, 0, schedule_name)
            WaitForJob(refine_job, 30)

        halfmap_filename = None
        try:
            job_star = load_star(refine_job + "job_pipeline.star")
            for output_file in job_star["pipeline_output_edges"]['rlnPipeLineEdgeToNode']:
                if output_file.endswith("half1_class001_unfil.mrc"):
                    halfmap_filename = output_file
                    break
            assert halfmap_filename != None
        except:
            print(" RELION_IT: Refinement job " + refine_job + " does not contain expected output maps.")
            print(" RELION_IT: This job should have finished, but you may continue it from the GUI.")
            print(" RELION_IT: For now, making the plot without this job.")

        if halfmap_filename is not None:
            # C. Run PostProcess
            postprocess_run_file = opts.input_postprocess_job+'job.star'
            if not os.path.exists(postprocess_run_file):
                postprocess_run_file = opts.input_postprocess_job+'run.job'
            post_options = ['One of the 2 unfiltered half-maps: == {}'.format(halfmap_filename)]
            post_job_name = 'post_job_' + str(current_nr_particles)
            post_alias = opts.prefix + str(current_nr_particles)
            post_job, already_had_it = addJob('PostProcess', post_job_name, SETUP_CHECK_FILE, post_options, postprocess_run_file, post_alias)
            if not already_had_it:
                RunJobs([post_job], 1, 0, schedule_name)
                WaitForJob(post_job, 30)
        
            # Get resolution from
            post_star = post_job + 'postprocess.star'
            try:
                resolution, pp_bfactor = get_postprocess_result(post_star)
                nr_particles.append(current_nr_particles)
                resolutions.append(resolution)
                pp_bfactors.append(pp_bfactor)
            except:
                print(' RELION_IT: WARNING: Failed to get post-processed resolution for {} particles'.format(current_nr_particles))

        # Update the current number of particles
        current_nr_particles = 2 * current_nr_particles

    # Also include the result from the original PostProcessing job
    if all_nr_particles <= opts.maximum_nr_particles:
        nr_particles.append(all_nr_particles)
        resolutions.append(all_particles_resolution)
        pp_bfactors.append(all_particles_bfactor)

    # Now already make preliminary plots here, e.g
    print()
    print('NrParticles Ln(NrParticles) Resolution(A) 1/Resolution^2 PostProcessBfactor')
    xs = []
    ys = []
    for n_particles, resolution, pp_bfactor in zip(nr_particles, resolutions, pp_bfactors):
        log_n_particles = log(n_particles)
        inv_d2 = 1.0 / (resolution * resolution)
        print('{0:11d} {1:15.3f} {2:13.2f} {3:14.4f} {4:18.2f}'.format(n_particles,log_n_particles, resolution, inv_d2, -pp_bfactor))
        xs.append(log_n_particles)
        ys.append(inv_d2)
    slope, intercept = line_fit(xs, ys)
    b_factor = 2.0 / slope
    print()
    print(" RELION_IT: ESTIMATED B-FACTOR from {0:d} points is {1:.2f}".format(len(xs), b_factor))
    print(" RELION_IT: The fitted line is: Resolution = 1 / Sqrt(2 / {0:.3f} * Log_e(#Particles) + {1:.3f})".format(b_factor, intercept))
    print(" RELION_IT: IF this trend holds, you will get:")
    for x in (1.5, 2, 4, 8):
        current_nr_particles = int(all_nr_particles * x)
        resolution = 1 / sqrt(slope * log(current_nr_particles) + intercept)
        print(" RELION_IT:   {0:.2f} A from {1:d} particles ({2:d} % of the current number of particles)".format(resolution, current_nr_particles, int(x * 100)))
    if True:#try: # Try plotting
        import matplotlib as mpl
        mpl.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
        fitted = []
        for x in xs:
            fitted.append(x * slope + intercept)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(xs, ys, '.')
        ax1.plot(xs, fitted)
        ax1.set_xlabel("ln(#particles)")
        ax1.set_ylabel("1/Resolution$^2$ in 1/$\AA^2$")
        ax1.set_title("Rosenthal & Henderson plot: B = 2.0 / slope = {:.1f}".format(b_factor));

        ax2 = ax1.twiny()
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")
        ax2.set_xlim(ax1.get_xlim())
        ax2.spines["bottom"].set_position(("axes", -0.15)) # In matplotlib 1.2, the order seems to matter
        ax2.set_xlabel("#particles")
        ax2.set_xticklabels(np.exp(ax1.get_xticks()).astype(np.int))

        ax3 = ax1.twinx()
        ax3.set_ylabel("Resolution in $\AA$")
        ax3.set_ylim(ax1.get_ylim())
        ax3.yaxis.set_ticks_position("right")
        ax3.yaxis.set_label_position("right")
        yticks = ax1.get_yticks()
        yticks[yticks <= 0] = 1.0 / (999 * 999) # to avoid zero division and negative sqrt
        ndigits = 1
        if np.max(yticks) > 0.25:
            ndigits = 2
        ax3.set_yticklabels(np.sqrt(1 / yticks).round(ndigits))

        output_name = opts.prefix + "rosenthal-henderson-plot.pdf"
        plt.savefig(output_name, bbox_inches='tight')
        print(" RELION_IT: Plot written to " + output_name)
    else:#except:
        print('WARNING: Failed to plot. Probably matplotlib and/or numpy is missing.')

    if os.path.isfile(RUNNING_FILE):
        os.remove(RUNNING_FILE)

    print(' RELION_IT: exiting now... ')

def main():
    """
    Run the RELION 3 pipeline.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """

    global RUNNING_FILE
    global SETUP_CHECK_FILE

    opts = RelionItOptions()
    for user_opt_file in sys.argv[1:]:
        print(' RELION_IT: reading options from {}'.format(user_opt_file))
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)

    SETUP_CHECK_FILE = opts.prefix + SETUP_CHECK_FILE
    RUNNING_FILE = opts.prefix + RUNNING_FILE

    # Make sure no other version of this script are running...
    if os.path.isfile(RUNNING_FILE):
        print(" RELION_IT: ERROR:", RUNNING_FILE, "is already present: delete this file and make sure no other copy of this script is running. Exiting now ...")
        exit(0)

    print(' RELION_IT: -------------------------------------------------------------------------------------------------------------------')
    print(' RELION_IT: Script for automated Bfactor-plot generation in RELION (>= 3.1)')
    print(' RELION_IT: Authors: Sjors H.W. Scheres & Takanori Nakane')
    print(' RELION_IT: ')
    print(' RELION_IT: Usage: ./bfactor_plot.py [extra_options.py ...]')
    print(' RELION_IT: ')
    print(' RELION_IT: This script keeps track of already submitted jobs in a filed called', SETUP_CHECK_FILE)
    print(' RELION_IT:   upon a restart, jobs present in this file will be ignored.')
    print(' RELION_IT: If you would like to re-do a specific job from scratch (e.g. because you changed its parameters)')
    print(' RELION_IT:   remove that job, and those that depend on it, from the', SETUP_CHECK_FILE)
    print(' RELION_IT: -------------------------------------------------------------------------------------------------------------------')
    print(' RELION_IT: ')
    
    run_pipeline(opts)

if __name__ == "__main__":
    main()
