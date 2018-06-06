#!/usr/bin/env python2.7
"""
bfactor_plot
---------

Pipeline setup script for automated processing with RELION 3.

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
from math import log

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING_BFACTOR_PLOT'
SETUP_CHECK_FILE = 'BFACTOR_PLOT_SUBMITTED_JOBS'

class RelionItOptions(object):
    """
    Options for the relion_it pipeline setup script.
    
    When initialised, this contains default values for all options. Call
    ``update_from()`` to override the defaults with a dictionary of new values.
    """
    #############################################################################
    # Change the parameters below to reflect your experiment                    #
    #############################################################################

    # Refine3D job with all particles
    input_refine3d_job = 'Refine3D/job040/' # 24
    # PostProcess job for resolution assessment
    input_postprocess_job = 'PostProcess/job083/' # 26
    # Minimum number of particles
    minimum_nr_particles = 100
    # Maximum number of particles
    maximum_nr_particles = 9999999

    #### relion_refine paremeters 
    # Initial low-pass filter for the refinements
    refine_ini_lowpass = 40
    # Read all particles in one batch into memory?
    refine_preread_images = True
    # Or copy particles to scratch disk?
    refine_scratch_disk = '/ssd'
    # Number of pooled particles?
    refine_nr_pool = 10
    # Use GPU-acceleration?
    refine_do_gpu = True
    # Which GPU to use (different from GPU used for pre-processing?)
    refine_gpu = ''
    # How many MPI processes to use
    refine_mpi = 3
    # How many threads to use
    refine_threads = 12
    # Skip padding?
    refine_skip_padding = False
    # Submit jobs to the cluster?
    refine_submit_to_queue = False


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

def appendJobOptionsFromRunJobFile(filename, job_options):

    f = open(filename,'r')
    for line in f: 
        label = line.split(' == ')[0]
        already_have_label = False
        if not (('job_type' in label) or ('is_continue' in label)):
            # search for this label in the input job_options
            for option in job_options:
                if ( label in option ):
                    already_have_label = True
            if not already_have_label:
                job_options.append(line.replace('\n',''))
    f.close()
    return


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

def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print " RELION_IT:", RUNNING_FILE, "file no longer exists, exiting now ..."
        exit(0)

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
    refine3d_run_file = opts.input_refine3d_job+'run.job'
    f = open(refine3d_run_file,'r')
    all_particles_star_file = ''
    for line in f: 
        if 'Input images STAR file' in line:
            all_particles_star_file = line.split(' == ')[1].replace('\n','')
            break
    if (all_particles_star_file == ''):
        print ' ERROR: cannot find input STAR file in', refine3d_run_file
        exit(1)

    all_particles = load_star(all_particles_star_file)
    all_nr_particles = len(all_particles['']['rlnImageName'])
    all_particles_resolution = float(load_star(opts.input_postprocess_job + 'postprocess.star')['general']['rlnFinalResolution'])

    nr_particles = []
    resolutions = []

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
        split_job, already_had_it = addJob('Select', split_job_name, SETUP_CHECK_FILE, split_options)
        if not already_had_it:
            RunJobs([split_job], 1, 0, schedule_name)
            WaitForJob(split_job, 30)

        # B. Run Refine3D
        refine_options = ['Input images STAR file: == {}particles_split001.star'.format(split_job),
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

        appendJobOptionsFromRunJobFile(refine3d_run_file, refine_options)
        refine_job_name = 'refine_job_' + str(current_nr_particles)
        refine_job, already_had_it = addJob('Refine3D', refine_job_name, SETUP_CHECK_FILE, refine_options)
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
		print " RELION_IT: Refinement job " + refine_job + " does not contain expected output maps."
		print " RELION_IT: This job should have finished, but you may continue it from the GUI. "
		print " RELION_IT: For now, making the plot without this job."

	if halfmap_filename is not None:
	        # C. Run PostProcess
	        postprocess_run_file = opts.input_postprocess_job+'run.job'
	        post_options = ['One of the 2 unfiltered half-maps: == {}'.format(halfmap_filename)]
	        appendJobOptionsFromRunJobFile(postprocess_run_file, post_options)
	        post_job_name = 'post_job_' + str(current_nr_particles)
	        post_job, already_had_it = addJob('PostProcess', post_job_name, SETUP_CHECK_FILE, post_options)
	        if not already_had_it:
	            RunJobs([post_job], 1, 0, schedule_name)
	            WaitForJob(post_job, 30)
        
	        # Get resolution from
	        post_star = post_job + 'postprocess.star'
        	try:
	            resolution = float(load_star(post_star)['general']['rlnFinalResolution'])
	            nr_particles.append(current_nr_particles)
	            resolutions.append(resolution)
	        except:
	            print ' RELION_IT: WARNING: Failed to get post-processed resolution for {} particles'.format(current_nr_particles)

        # Update the current number of particles
        current_nr_particles = 2 * current_nr_particles

    # Also include the result from the original PostProcessing job
    if all_nr_particles <= opts.maximum_nr_particles:
        nr_particles.append(all_nr_particles)
        resolutions.append(all_particles_resolution)

    # Now already make preliminary plots here, e.g
    print
    print 'NrParticles Ln(NrParticles) Resolution(A) 1/Resolution^2'
    xs = []
    ys = []
    for n_particles, resolution in zip(nr_particles, resolutions):
        log_n_particles = log(n_particles)
        inv_d2 = 1.0 / (resolution * resolution)
        print '{0:11d} {1:15.3f} {2:13.2f} {3:14.4f}'.format(n_particles,log_n_particles, resolution, inv_d2)
        xs.append(log_n_particles)
        ys.append(inv_d2)
    slope, intercept = line_fit(xs, ys)
    b_factor = 2.0 / slope
    print
    print " RELION_IT: ESTIMATED B-FACTOR from {0:d} points is {1:.2f}".format(len(xs), b_factor)
    try: # Try plotting
        import matplotlib.pyplot as plt
        fitted = []
        for x in xs:
            fitted.append(x * slope + intercept)
        plt.plot(xs, ys, '.')
        plt.plot(xs, fitted)
        plt.xlabel("ln(#particles)")
        plt.ylabel("1/Resolution$^2$ in 1/$\AA^2$")
        plt.title("Rosenthal & Henderson plot: B = 2.0 / slope = {:.1f}".format(b_factor));
        plt.savefig("rosenthal-henderson-plot.pdf", bbox_inches='tight')
        print "Plot written to rosenthal-henderson-plot.pdf."
    except:
        print 'WARNING: Failed to plot. Probably matplotlib is missing.'

    if os.path.isfile(RUNNING_FILE):
        os.remove(RUNNING_FILE)

    print ' RELION_IT: Finished all refinements, the plot was written to rosenthal-henderson-plot.pdf.'
    print ' RELION_IT: exiting now... '


def main():
    """
    Run the RELION 3 pipeline.
    
    Options files given as command line arguments will be opened in order and
    used to update the default options.
    """
    print ' RELION_IT: -------------------------------------------------------------------------------------------------------------------'
    print ' RELION_IT: script for automated bfactor-plot generation in RELION (>= 3.0-alpha-5)'
    print ' RELION_IT: authors: Sjors H.W. Scheres & Takanori Nakane'
    print ' RELION_IT: '
    print ' RELION_IT: usage: ./bfactor_plot.py [extra_options.py ...]'
    print ' RELION_IT: '
    print ' RELION_IT: this script keeps track of already submitted jobs in a filed called',SETUP_CHECK_FILE
    print ' RELION_IT:   upon a restart, jobs present in this file will be ignored'
    print ' RELION_IT: if you would like to re-do a specific job from scratch (e.g. because you changed its parameters)' 
    print ' RELION_IT:   remove that job, and those that depend on it, from the',SETUP_CHECK_FILE
    print ' RELION_IT: -------------------------------------------------------------------------------------------------------------------'
    print ' RELION_IT: '
    
    # Make sure no other version of this script are running...
    if os.path.isfile(RUNNING_FILE):
        print " RELION_IT: ERROR:", RUNNING_FILE, "is already present: delete this file and make sure no other copy of this script is running. Exiting now ..."
        exit(0)

    opts = RelionItOptions()
    for user_opt_file in sys.argv[1:]:
        print ' RELION_IT: reading options from {}'.format(user_opt_file)
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)
    run_pipeline(opts)


if __name__ == "__main__":
    main()






