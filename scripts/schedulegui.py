#!/usr/bin/env python2.7
"""
relion_it.py
============

Simple GUI to abort and restart RELION-4.0 Schedules.
"""

from __future__ import print_function
from __future__ import division  # always use float division

import argparse
import sys
import glob
import inspect
import math
from os import path
import os
import re
import runpy
import time
import traceback
import ast 
import collections
from shutil import copytree
import Tkinter as tk
#from Tkinter import scrolledtext 
import tkMessageBox

OPTIONS_FILE = 'relion_it_options.py'
IS_RUNNING = False

# Colour definitions
# Yellowish background for entries
entry_bg = '#ffffe6'
# reddish colour for Browse buttons
button_bg = '#c8506e'
# Darker red for run buttons
runbutton_bg = '#a01e3c'

def load_star(filename):
    from collections import OrderedDict
    import shlex  
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
            elems = shlex.split(line)
            assert len(elems) == len(current_colnames)
            for idx, e in enumerate(elems):
                current_data[current_colnames[idx]].append(e)

    return datasets

class ScheduleGui(object):

    def __init__(self, main_window, schedulename):
        self.main_window = main_window
        self.schedulename = schedulename

        self.main_frame = tk.Frame(self.main_window)
        self.main_frame.pack(fill=tk.BOTH, expand=1)
        self.runoutsize = 0
        self.runerrsize = 0

        ### Status frame

        self.status_frame = tk.LabelFrame(self.main_frame, text="Current status", padx=3, pady=3)
        self.status_frame.pack(padx=3, pady=3, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        tk.Label(self.status_frame, text="Current:").grid(row=0, column=0, sticky=tk.W)
        self.current_node_var = tk.StringVar()
        self.current_node_entry = tk.Entry(self.status_frame, textvariable=self.current_node_var, bg=entry_bg)
        self.current_node_entry.grid(row=0, column=1, sticky=tk.W)
        self.set_current_node()

        self.abort_button = tk.Button(self.status_frame, text="Abort", command=self.abort_schedule, bg=runbutton_bg)
        self.abort_button.grid(row=0, column=2, sticky=tk.W+tk.E)

        self.unlock_schedule_button = tk.Button(self.status_frame, text="Unlock", command=self.unlock_schedule, bg=runbutton_bg)
        self.unlock_schedule_button.grid(row=0, column=3, sticky=tk.W)

        self.reset_schedule_button = tk.Button(self.status_frame, text="Reset", command=self.reset_schedule, bg=runbutton_bg)
        self.reset_schedule_button.grid(row=0, column=4, sticky=tk.W)

        self.restart_schedule_button = tk.Button(self.status_frame, text="Restart", command=self.restart_schedule, bg=runbutton_bg)
        self.restart_schedule_button.grid(row=0, column=5)
        tk.Label(self.status_frame, text="at:").grid(row=0, column=6, sticky=tk.W)

        schedule = load_star('Schedules/'+self.schedulename+'/schedule.star')
        jobnames = schedule['schedule_jobs']['rlnScheduleJobNameOriginal']
        operators = schedule['schedule_operators']['rlnScheduleOperatorName']
        restartvars = jobnames + operators
        self.restart_var = tk.StringVar()
        self.restart_entry = tk.OptionMenu(self.status_frame, self.restart_var, *restartvars)
        self.restart_entry.grid(row=0, column=7, sticky=tk.W)

        left_frame = tk.Frame(self.main_frame)
        left_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        right_frame = tk.Frame(self.main_frame)
        right_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)
        

        ### Joboption frame
        
        self.joboption_frame = tk.LabelFrame(left_frame, text="Set Job option", padx=3, pady=3)
        self.joboption_frame.pack(padx=3, pady=3, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        row = 0

        ### Fill the jobnames and schedule variables pull-down menus
        
        tk.Label(self.joboption_frame, text="Job name:").grid(row=row, sticky=tk.W)
        self.jobname_var = tk.StringVar()
        self.jobname_entry = tk.OptionMenu(self.joboption_frame, self.jobname_var, *jobnames)
        self.jobname_entry.grid(row=row, column=1, sticky=tk.W)
        self.jobname_var.trace('w', self.change_jobname)

        row += 1

        tk.Label(self.joboption_frame, text="Job option:").grid(row=row, sticky=tk.W)
        self.joboption_var = tk.StringVar()
        self.joboptions = ["undefined"] 
        self.joboption_entry = tk.OptionMenu(self.joboption_frame, self.joboption_var, *self.joboptions)
        self.joboption_entry.grid(row=row, column=1, sticky=tk.W)
        self.joboption_var.trace('w', self.change_joboption)
        
        row += 1

        tk.Label(self.joboption_frame, text="New value:").grid(row=row, sticky=tk.W)
        self.joboption_value_var = tk.StringVar()
        self.joboption_value_entry = tk.Entry(self.joboption_frame, textvariable=self.joboption_value_var, bg=entry_bg)
        self.joboption_value_entry.grid(row=row, column=1, sticky=tk.W)

        self.set_joboption_button = tk.Button(self.joboption_frame, text="Set", command=self.set_joboption, bg=button_bg)
        self.set_joboption_button.grid(row=row, column=2)

        ### Schedulevar frame

        self.schedulevar_frame = tk.LabelFrame(left_frame, text="Set Schedule variable", padx=3, pady=3)
        self.schedulevar_frame.pack(padx=3, pady=3, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        row = 0

        self.schedulevars = []
        try:
            self.schedulevars += schedule['schedule_floats']['rlnScheduleFloatVariableName']
            self.schedulevars += schedule['schedule_bools']['rlnScheduleBooleanVariableName']
            self.schedulevars += schedule['schedule_strings']['rlnScheduleStringVariableName']
        except KeyError:
            pass
            
        tk.Label(self.schedulevar_frame, text="Variable:").grid(row=row, sticky=tk.W)
        self.schedulevar_var = tk.StringVar()
        self.schedulevar_entry = tk.OptionMenu(self.schedulevar_frame, self.schedulevar_var, *self.schedulevars)
        self.schedulevar_entry.grid(row=row, column=1, sticky=tk.W)
        self.schedulevar_var.trace('w', self.change_schedulevar)

        row += 1

        tk.Label(self.schedulevar_frame, text="New value:").grid(row=row, sticky=tk.W)
        self.schedulevar_value_var = tk.StringVar()
        self.schedulevar_value_entry = tk.Entry(self.schedulevar_frame, textvariable=self.schedulevar_value_var, bg=entry_bg)
        self.schedulevar_value_entry.grid(row=row, column=1, sticky=tk.W)

        set_schedulevar_button = tk.Button(self.schedulevar_frame, text="Set", command=self.set_schedulevar, bg=button_bg)
        set_schedulevar_button.grid(row=row, column=2)

        #### Output frame
        
        self.output_frame = tk.LabelFrame(right_frame, text='Schedules/'+self.schedulename+'/run.out', padx=3, pady=1)
        self.output_frame.pack(padx=3, pady=1, fill=tk.X, expand=0)
        
        scroll = tk.Scrollbar(self.output_frame)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.stdout_entry = tk.Text(self.output_frame, wrap=tk.WORD, yscrollcommand=scroll.set, bg="white", height=8, font=("TkDefaultFont", 9))
        self.stdout_entry.pack()
        self.stdout_entry.yview_moveto('1.0')

        #### Error frame
        
        self.error_frame = tk.LabelFrame(right_frame, text='Schedules/'+self.schedulename+'/run.err', padx=3, pady=1)
        self.error_frame.pack(padx=3, pady=1, fill=tk.X, expand=0)
        
        scroll2 = tk.Scrollbar(self.error_frame)
        scroll2.pack(side=tk.RIGHT, fill=tk.Y)
        self.stderr_entry = tk.Text(self.error_frame, wrap=tk.WORD, yscrollcommand=scroll2.set, fg="red", bg="white", height=3, font=("TkDefaultFont", 9))
        self.stderr_entry.pack()
        self.stderr_entry.yview_moveto('1.0')
        
        # Check whether the Schedule is running every 3 seconds
        self.main_frame.after(0, self.check_running())
        
    def change_jobname(self, *args):
        job = load_star('Schedules/'+self.schedulename+'/'+self.jobname_var.get()+'/job.star')
        self.joboptions = job['joboptions_values']['rlnJobOptionVariable']
        self.joboption_entry['menu'].delete(0, 'end')
        for opt in self.joboptions:
            self.joboption_entry['menu'].add_command(label=opt, command=tk._setit(self.joboption_var, opt))
        self.joboption_value_var.set("")

    def change_joboption(self, *args):
        job = load_star('Schedules/'+self.schedulename+'/'+self.jobname_var.get()+'/job.star')
        myopt = self.joboption_var.get()
        myval = ""
        for i in range(0,len(job['joboptions_values']['rlnJobOptionVariable'])):
            if job['joboptions_values']['rlnJobOptionVariable'][i] == myopt:
                myval += job['joboptions_values']['rlnJobOptionValue'][i]
        self.joboption_value_var.set(myval)

    def change_schedulevar(self, *args):
        schedule = load_star('Schedules/'+self.schedulename+'/schedule.star')
        myvar = self.schedulevar_var.get()
        myval = ""
        try:
            ifloat = ibool = istring = -1
            for i in range(0,len(schedule['schedule_floats']['rlnScheduleFloatVariableName'])):
                if schedule['schedule_floats']['rlnScheduleFloatVariableName'][i] == myvar:
                    myval += schedule['schedule_floats']['rlnScheduleFloatVariableValue'][i]
            for i in range(0,len(schedule['schedule_bools']['rlnScheduleBooleanVariableName'])):
                if schedule['schedule_bools']['rlnScheduleBooleanVariableName'][i] == myvar:
                    myval += schedule['schedule_bools']['rlnScheduleBooleanVariableValue'][i]
            for i in range(0,len(schedule['schedule_strings']['rlnScheduleStringVariableName'])):
                if schedule['schedule_strings']['rlnScheduleStringVariableName'][i] == myvar:
                    myval += schedule['schedule_strings']['rlnScheduleStringVariableValue'][i]
        except KeyError:
            pass
        self.schedulevar_value_var.set(myval)

    def update_relion_it_options(self, myoption):
        with open(OPTIONS_FILE) as file: 
            user_opts = collections.OrderedDict(ast.literal_eval(file.read()))
        user_opts.update(myoption)
        with open(OPTIONS_FILE, 'w') as file:
            file.write("{\n") 
            for k,v in user_opts.items():
                file.write("'%s' : '%s', \n" % (k, v))
            file.write("}\n")                  
        print(' RELION_IT: updated ', OPTIONS_FILE)

    def set_joboption(self, *args_ignored, **kwargs_ignored):
        jobstar = 'Schedules/' + self.schedulename + '/' + self.jobname_var.get() + '/job.star'
        command = 'relion_pipeliner --editJob ' + jobstar + ' --editOption ' + self.joboption_var.get() + ' --editValue \"' + str(self.joboption_value_entry.get()) + '\"'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        # Also set has_started to false for this job in the scheduler
        command2 = 'relion_scheduler --schedule ' + self.schedulename + ' --set_has_started ' + self.jobname_var.get() + ' --value False'
        print(' RELION_IT: excuting: ', command2)
        os.system(command2)
        myoption = self.schedulename + '__' + self.jobname_var.get() + '__' + self.joboption_var.get()
        mydic = {myoption : self.joboption_value_entry.get()}
        self.update_relion_it_options(mydic)
            
    def set_schedulevar(self, *args_ignored, **kwargs_ignored):
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --set_var ' + self.schedulevar_var.get() + ' --value \"' + self.schedulevar_value_var.get() + '\"' + ' --original_value \"' + self.schedulevar_value_var.get() + '\"'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        myoption = self.schedulename + '__' + self.schedulevar_var.get()
        mydic = {myoption : self.schedulevar_value_var.get()}
        self.update_relion_it_options(mydic)

    def abort_schedule(self, *args_ignored, **kwargs_ignored):
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --abort'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        self.set_current_node()
        self.current_node_entry.configure(state=tk.DISABLED)

    def reset_schedule(self, *args_ignored, **kwargs_ignored):
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --reset'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        self.set_current_node()

    def unlock_schedule(self, *args_ignored, **kwargs_ignored):
        if tkMessageBox.askokcancel("Confirm unlock", 'Use:\n\nps -ef | grep relion_scheduler \n\nto confirm this scheduler is no longer running. \n\nOK to unlock?', 
                                    icon='warning', default=tkMessageBox.CANCEL):
            lockname = ".relion_lock_schedule_" + self.schedulename
            command = 'rm -rf ' + lockname
            print(' RELION_IT: excuting: ', command)
            os.system(command)

    def restart_schedule(self, *args_ignored, **kwargs_ignored):
        # Set the current node, as per the GUI
        myrestartpoint = self.restart_var.get()
        if myrestartpoint == "":
            myrestartpoint = self.current_node_var.get() 
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --set_current_node ' + myrestartpoint
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        # And then run the Schedule again
        command2 = 'relion_scheduler --schedule ' + self.schedulename + ' --run --pipeline_control Schedules/' + self.schedulename + '/ >> Schedules/' +self.schedulename + '/run.out 2>> Schedules/' + self.schedulename + '/run.err &'
        print(' RELION_IT: excuting: ', command2)
        os.system(command2)

    def set_current_node(self):
        current_node = "unknown"
        schedulestar = 'Schedules/' + self.schedulename + '/schedule.star'
        if path.exists(schedulestar):
            for line in open(schedulestar, 'r'):
                if re.search('rlnScheduleCurrentNodeName', line):
                    current_node = line.split()[1]
        self.current_node_entry.configure(state=tk.NORMAL)
        self.current_node_entry.delete(0,tk.END)
        self.current_node_entry.insert(0, current_node)

    def check_running(self):
        lockname = ".relion_lock_schedule_" + self.schedulename
        # read run,out and run.err
        fn = 'Schedules/' + self.schedulename + '/run.'
        if path.exists(fn+'out'):
            mysize = os.path.getsize(fn+'out')
            if mysize > self.runoutsize:
                self.runoutsize = mysize
                with open(fn+'out', 'r') as f:
                    self.stdout_entry.configure(state=tk.NORMAL)
                    self.stdout_entry.delete('1.0', tk.END)
                    self.stdout_entry.insert(tk.END, f.read())
                    self.stdout_entry.yview_moveto('1.0') 
                    self.stdout_entry.configure(state=tk.DISABLED)
        if path.exists(fn+'err'):
            mysize = os.path.getsize(fn+'err')
            if mysize > self.runerrsize:
                self.runerrsize = mysize
                with open(fn+'err', 'r') as f:
                    self.stderr_entry.configure(state=tk.NORMAL)
                    self.stderr_entry.delete('1.0', tk.END)
                    self.stderr_entry.insert(tk.END, f.read())
                    self.stderr_entry.yview_moveto('1.0')
                    self.stderr_entry.configure(state=tk.DISABLED)
        if path.exists(lockname):
            self.set_current_node()
            self.status_frame.configure(text="Currently running", fg="dark green")
            self.current_node_entry.configure(state=tk.DISABLED)
            self.abort_button.configure(state=tk.NORMAL)
            for child in self.joboption_frame.winfo_children():
                child.configure(state=tk.DISABLED)
            for child in self.schedulevar_frame.winfo_children():
                child.configure(state=tk.DISABLED)
            self.unlock_schedule_button.configure(state=tk.NORMAL)
            self.restart_schedule_button.configure(state=tk.DISABLED)
            self.restart_entry.configure(state=tk.DISABLED)
            self.reset_schedule_button.configure(state=tk.DISABLED)
        else:
            self.status_frame.configure(text="Currently stopped", fg="dark red")
            self.current_node_entry.configure(state=tk.NORMAL)
            self.abort_button.configure(state=tk.DISABLED)
            for child in self.joboption_frame.winfo_children():
                child.configure(state=tk.NORMAL)
            for child in self.schedulevar_frame.winfo_children():
                child.configure(state=tk.NORMAL)
            self.unlock_schedule_button.configure(state=tk.DISABLED)
            self.restart_schedule_button.configure(state=tk.NORMAL)
            self.restart_entry.configure(state=tk.NORMAL)
            self.reset_schedule_button.configure(state=tk.NORMAL)
        self.main_frame.after(3000, self.check_running)
        
def main():

    if len(sys.argv) < 2:
        print(' RELION_IT: ERROR you need to provide the name of the Schedule as the first argument')
        exit(1)
    else:
       schedulename = sys.argv[1] 

    print(' RELION_IT: launching GUI...')
    tk_root = tk.Tk()
    titlegui = "Schedule: " + schedulename
    tk_root.title(titlegui)
    ScheduleGui(tk_root, schedulename)
    tk_root.mainloop()

if __name__ == "__main__":
    main()

