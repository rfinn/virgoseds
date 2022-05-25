#!/usr/bin/env python

'''


* run this from the directory where you want all the scripts created, e.g. ~/scripts/

'''


import os
import subprocess
import sys
import glob


HOME = os.getenv("HOME")


def write_output(filename,dirname,submit=False):
    ''' copying from Matt Bellis '''
    output = ""
    output += "#!/bin/bash"
    output += "\n"
    output += "# Set Job Name"
    output += "#SBATCH -J job"
    output += "\n"
    output += "# Set file to capture standard out and standard error and append the jobID (%j)"
    output += "#SBATCH -o job.out.%j"
    output += "\n"
    output += "#SBATCH --partition=normal"
    
    output += "#Set the number of nodes"
    output += "#SBATCH -N 1"
    output += "#SBATCH --ntasks=1"
    output += "\n"
    output += "#Set the time limit for the job - 10.5 hour is specified here"
    output += "#SBATCH --time=10:30:00"
    output += "\n"
    output += "#SBATCH --cpus-per-task=1"
    output += "\n"
    output += "# Load any environmental modules needed"
    output += "module load Python3"
    output += "module load gnu9"
    
    output += "# Move to the directory needed - defaults to the submission directory"
    output += "cd $HOME/research/Virgo/magphysParallel/output/"
    output += "\n"
    output += "# perform calculation"
    output += "\n"
    output += f"python {HOME}/github/virgoseds/python/run1magphys.py {dirname}"

    outfname = f"JOB_{filename}.sh"
    outfile = open(outfname,'w')
    
    outfile.write(output)
    outfile.close()

    cmds = ['sbatch', outfname]
    if submit:
        print(f"Submitting job to process {filename}")
        process = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()
        print(stdout.decode())
        print(stderr.decode())
        pass
        


data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
dirlist = glob.glob(f"{data_dir}VFID????")

dirlist.sort()

# write out files and submit jobs
for d in dirlist:
    # remove full path to directory so just VFID???? is passed in
    gname = os.path.basename(d)
    write_output(gname,gname,submit=False)
