#!/usr/bin/env python

'''

* this will split sample into groups of 1000 to meet max limit of arrays
* run this from the directory where you want all the scripts created, e.g. ~/scripts/

'''


import os
import subprocess
import sys
import glob


HOME = os.getenv("HOME")


def write_output(script_id,input_file,submit=False):
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
    output += "#SBATCH --array=1-1000\n"
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
    output += "cd $HOME/research/Virgo/magphysParallel/output/\n"
    output += "\n"
    output += "# perform calculation\n"
    output += "\n"
    output += f'LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {input_file})'    
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
for i in range(7):
    os.system(f"ls -d VFID{i}??? > Dirs{i}.txt")
os.chdir(cwd)


# write out files and submit jobs
#for d in dirlist:
for i in range(7):
    # remove full path to directory so just VFID???? is passed in
    script_id = f"VFID{i}000"
    input_file = f"Dirs{i}.txt"
    write_output(script_id,input_file,submit=False)
