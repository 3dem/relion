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

try:
    import Tkinter as tk
    from tkinter import scrolledtext 
    import tkMessageBox
    import tkFileDialog
except ImportError:
    # The GUI is optional. If the user requests it, it will fail when it tries
    # to open so we can ignore the error for now.
    pass  

OPTIONS_FILE = 'relion_it_options.py'
IS_RUNNING = False

# Colour definitions
# Yellowish background for entries
entry_bg = '#ffffe6'
# reddish colour for Browse buttons
button_bg = '#c8506e'
# Darker red for run buttons
runbutton_bg = '#a01e3c'

class ScheduleGui(object):

    def __init__(self, main_window, schedulename):
        self.main_window = main_window
        self.schedulename = schedulename

        self.main_frame = tk.Frame(self.main_window)
        self.main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame = tk.Frame(self.main_frame)
        left_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)

        right_frame = tk.Frame(self.main_frame)
        right_frame.pack(side=tk.LEFT, anchor=tk.N, fill=tk.X, expand=1)
        
        ### Status frame

        self.status_frame = tk.LabelFrame(left_frame, text="Current status", padx=5, pady=5)
        self.status_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        row = 0

        tk.Label(self.status_frame, text="Current:").grid(row=row, sticky=tk.W)
        self.current_node_var = tk.StringVar()
        self.current_node_entry = tk.Entry(self.status_frame, textvariable=self.current_node_var, bg=entry_bg)
        self.current_node_entry.grid(row=row, column=1, sticky=tk.W)
        self.set_current_node()

        self.abort_button = tk.Button(self.status_frame, text="Abort", command=self.abort_schedule, bg=button_bg)
        self.abort_button.grid(row=row, column=2, sticky=tk.W)

        row += 1

        ### Joboption frame
        
        self.joboption_frame = tk.LabelFrame(left_frame, text="Set Job option", padx=5, pady=5)
        self.joboption_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        row = 0

        tk.Label(self.joboption_frame, text="Job name:").grid(row=row, sticky=tk.W)
        self.jobname_var = tk.StringVar()
        self.jobname_entry = tk.Entry(self.joboption_frame, textvariable=self.jobname_var, bg=entry_bg)
        self.jobname_entry.grid(row=row, column=1, sticky=tk.W)

        row += 1

        tk.Label(self.joboption_frame, text="Job option:").grid(row=row, sticky=tk.W)
        self.joboption_var = tk.StringVar()
        self.joboption_entry = tk.Entry(self.joboption_frame, textvariable=self.joboption_var, bg=entry_bg)
        self.joboption_entry.grid(row=row, column=1, sticky=tk.W)

        row += 1

        tk.Label(self.joboption_frame, text="New value:").grid(row=row, sticky=tk.W)
        self.joboption_value_var = tk.StringVar()
        self.joboption_value_entry = tk.Entry(self.joboption_frame, textvariable=self.joboption_value_var, bg=entry_bg)
        self.joboption_value_entry.grid(row=row, column=1, sticky=tk.W)

        self.set_joboption_button = tk.Button(self.joboption_frame, text="Set", command=self.set_joboption, bg=button_bg)
        self.set_joboption_button.grid(row=row, column=2)

        ### Schedulevar frame

        self.schedulevar_frame = tk.LabelFrame(left_frame, text="Set Schedule variable", padx=5, pady=5)
        self.schedulevar_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        row = 0

        tk.Label(self.schedulevar_frame, text="Variable:").grid(row=row, sticky=tk.W)
        self.schedulevar_var = tk.StringVar()
        self.schedulevar_entry = tk.Entry(self.schedulevar_frame, textvariable=self.schedulevar_var, bg=entry_bg)
        self.schedulevar_entry.grid(row=row, column=1, sticky=tk.W)

        row += 1

        tk.Label(self.schedulevar_frame, text="New value:").grid(row=row, sticky=tk.W)
        self.schedulevar_value_var = tk.StringVar()
        self.schedulevar_value_entry = tk.Entry(self.schedulevar_frame, textvariable=self.schedulevar_value_var, bg=entry_bg)
        self.schedulevar_value_entry.grid(row=row, column=1, sticky=tk.W)

        set_schedulevar_button = tk.Button(self.schedulevar_frame, text="Set", command=self.set_schedulevar, bg=button_bg)
        set_schedulevar_button.grid(row=row, column=2)

        #### Restart frame
        
        self.restart_frame = tk.LabelFrame(left_frame, text="Restart Schedule", padx=5, pady=5)
        self.restart_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.main_frame, 1, weight=1)

        self.restart_schedule_button = tk.Button(self.restart_frame, text="Restart", command=self.restart_schedule, bg=runbutton_bg)
        self.restart_schedule_button.grid(row=0, column=0)

        self.reset_schedule_button = tk.Button(self.restart_frame, text="Reset", command=self.reset_schedule, bg=runbutton_bg)
        self.reset_schedule_button.grid(row=0, column=1)

        self.unlock_schedule_button = tk.Button(self.restart_frame, text="Unlock", command=self.unlock_schedule, bg=runbutton_bg)
        self.unlock_schedule_button.grid(row=0, column=2)

        #### Output frame
        
        self.output_frame = tk.LabelFrame(right_frame, text='Schedules/'+self.schedulename+'/run.out')
        self.output_frame.pack(padx=5, pady=5, fill=tk.X, expand=0)
        
        scroll = tk.Scrollbar(self.output_frame)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.stdout_entry = tk.Text(self.output_frame, wrap=tk.WORD, yscrollcommand=scroll.set, bg="white", height=13)
        self.stdout_entry.pack()
        self.stdout_entry.yview_moveto('1.0')

        #### Error frame
        
        self.error_frame = tk.LabelFrame(right_frame, text='Schedules/'+self.schedulename+'/run.err', height=30)
        self.error_frame.pack(padx=5, pady=2, fill=tk.X, expand=0)
        
        scroll2 = tk.Scrollbar(self.error_frame)
        scroll2.pack(side=tk.RIGHT, fill=tk.Y)
        self.stderr_entry = tk.Text(self.error_frame, wrap=tk.WORD, yscrollcommand=scroll2.set, fg="red", bg="white", height=5)
        self.stderr_entry.pack()
        self.stderr_entry.yview_moveto('1.0')
        
        # Check whether the Schedule is running every 3 seconds
        self.main_frame.after(0, self.check_running())
        

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
        jobstar = 'Schedules/' + self.schedulename + '/' + self.jobname_entry.get() + '/job.star'
        command = 'relion_pipeliner --editJob ' + jobstar + ' --editOption ' + self.joboption_entry.get() + ' --editValue \"' + str(self.joboption_value_entry.get()) + '\"'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        # Also set has_started to false for this job in the scheduler
        command2 = 'relion_scheduler --schedule ' + self.schedulename + ' --set_has_started ' + self.jobname_entry.get() + ' --value False'
        print(' RELION_IT: excuting: ', command2)
        os.system(command2)
        myoption = self.schedulename + '__' + self.jobname_entry.get() + '__' + self.joboption_entry.get()
        mydic = {myoption : self.joboption_value_entry.get()}
        self.update_relion_it_options(mydic)
            
    def set_schedulevar(self, *args_ignored, **kwargs_ignored):
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --set_var ' + self.schedulevar_entry.get() + ' --value \"' + self.schedulevar_value_entry.get() + '\"' + ' --original_value \"' + self.schedulevar_value_entry.get() + '\"'
        print(' RELION_IT: excuting: ', command)
        os.system(command)
        myoption = self.schedulename + '__' + self.schedulevar_entry.get()
        mydic = {myoption : self.schedulevar_value_entry.get()}
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

    def unlock_schedule(self, *args_ignored, **kwargs_ignored):
        print(' RELION_IT: use the following command to ensure the relion_scheduler process is no longer running:') 
        print('   ps -ef | grep relion_scheduler')
        print(' RELION_IT: if you see the relion_scheduler process, use the Abort button instead!')
        print(' RELION_IT: otherwise, you can unlock the Schedule by typing:')
        lockname = ".relion_lock_schedule_" + self.schedulename
        print('   rm -rf ', lockname)

    def restart_schedule(self, *args_ignored, **kwargs_ignored):
        # Set the current node, as per the GUI
        command = 'relion_scheduler --schedule ' + self.schedulename + ' --set_current_node ' + self.current_node_entry.get()
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
            with open(fn+'out', 'r') as f:
                self.stdout_entry.configure(state=tk.NORMAL)
                self.stdout_entry.delete('1.0', tk.END)
                self.stdout_entry.insert(tk.END, f.read())
                self.stdout_entry.yview_moveto('1.0') 
                self.stdout_entry.configure(state=tk.DISABLED)
        if path.exists(fn+'err'):
            with open(fn+'err', 'r') as f:
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

