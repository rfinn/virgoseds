#!/usr/bin/env python

'''

GOAL:
* write an sbatch array file to compute the grid of libary models for Virgo
* this should save time b/c we need 134 models, if we let dz=.0001
* this is def faster than making library for each of 6k galaxies
* will need to make two sets, one for the N and one for S filters
USAGE:
* min and max recession velocities and dz are hardcoded
  * need to adjust those manually (or add argparse)

* to launch, set submit=True

DECIDED NOT TO GO THIS ROUTE...


'''


import os
import subprocess
import sys
import glob


HOME = os.getenv("HOME")


def write_output(script_id,input_file,narray=1000,submit=False):
    ''' copying from Matt Bellis '''
    output = ""
    output += "#!/bin/bash\n"
    output += "\n"
    output += "# Set Job Name\n"
    output += "#SBATCH -J job\n"
    output += "\n"
    output += "# Set file to capture standard out and standard error and append the jobID (%j)\n"
    output += "#SBATCH -o job.out.%j\n"
    output += "\n"
    output += "#SBATCH --partition=normal\n"
    output += "\n"    
    output += "# for testing\n"
    output += f"#SBATCH --array=1-{narray}\n"
    output += "\n"
    output += "#Set the number of nodes\n"
    output += "#SBATCH -N 1\n"
    output += "#SBATCH --ntasks=1\n"
    output += "\n"
    output += "#Set the time limit for the job - 10.5 hour is specified here\n"
    output += "#SBATCH --time=10:30:00\n"
    output += "\n"
    output += "#SBATCH --cpus-per-task=1\n"
    output += "\n"
    output += "# Load any environmental modules needed\n"
    output += "module load Python3\n"
    output += "module load gnu9\n"
    output += "\n"
    output += "# Move to the directory needed - defaults to the submission directory\n"
    output += "\n"
    output += "# perform calculation\n"
    output += "#\n"

    # set redshift using SLURM_ARRAY_TASK_ID
    output += "z=0.0001*$
    echo $z "70.,0.30,0.70" | xargs -n 1 | $magphys/get_optic_colors
    echo $z "70.,0.30,0.70" | xargs -n 1 | $magphys/get_infrared_colors


    
    output += f'LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $HOME/research/Virgo/magphysParallel/output/{input_file})\n'
    output += "#\n"    
    output += f"python {HOME}/github/virgoseds/python/run1magphys.py $LINE\n"

    outfname = f"JOB_{script_id}.sh"
    
    outfile = open(outfname,'w')
    outfile.write(output)
    outfile.close()

    cmds = ['sbatch', outfname]
    if submit:
        print(f"Submitting job to process {input_file}")
        process = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()
        print(stdout.decode())
        print(stderr.decode())
        pass
        

###########################################################
cwd = os.getcwd()

data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
os.chdir(data_dir)
nfiles = []
for i in range(7):
    outfile = f"Dirs{i}.txt"
    os.system(f"ls -d VFID{i}??? > {outfile}")
    # count lines
    infile = open(outfile,'r')
    nfiles.append(len(infile.readlines()))
    infile.close()
os.chdir(cwd)


# write out files and submit jobs
#for d in dirlist:
for i in range(7):
    # remove full path to directory so just VFID???? is passed in
    script_id = f"VFID{i}000"
    input_file = f"Dirs{i}.txt"
    write_output(script_id,input_file,narray=nfiles[i],submit=False)
