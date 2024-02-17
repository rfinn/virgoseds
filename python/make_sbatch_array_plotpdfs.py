#!/usr/bin/env python

'''
GOAL:
* This will create a bash script to run multiple serial jobs in parallel using slurms array option

USAGE:
* Run this in the output directory that has a subfolder for each galaxy.
* Try running this once and check the output of the JOB_{}.sh script.
* If it looks good, run again using the --submit flag to submit the script to slurm.

'''


import os
import subprocess
import sys
import glob


HOME = os.getenv("HOME")


def write_output(script_id,input_file,narray=1000,data_dir=None,submit=False):
    ''' copying python from Matt Bellis, commands from Ryan Decker '''
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
    # running in array mode, rather than spawning narray independent processes
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
    output += "# Purge modules \n"
    output += "module purge\n"
    output += "\n"    
    output += "# Load any environmental modules needed\n"
    output += "module load Python3\n"
    output += "module load gnu9\n"
    output += "\n"
    output += "# Move to the directory needed - defaults to the submission directory\n"
    output += "\n"
    output += "# perform calculation\n"
    output += "#\n"
    if data_dir is not None:
        s = f'LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {data_dir}/{input_file})\n'
    else:
        print('please provide a valid data directory')
        return
    output += s
    output += "#\n"    
    output += f"python {HOME}/github/virgoseds/python/plot1magphys.py $LINE\n"

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




###########################
##### SET UP ARGPARSE
###########################
import argparse
parser = argparse.ArgumentParser(description ='Program to create bash script to run magphys in parallel.  This uses the array option in slurm so that the processes are associated with the same process id, rather than submitting a bunch of individual jobs.  The script can be submitted by setting the --submit flag.')
parser.add_argument('--ext', dest = 'ext', default = 0, help = 'extinction correction to apply.  0=None; 1=Legacy Survey; 2=Salim/Leroy. Default is zero. The main difference between 1 and 2 is how the GALEX fluxes are handled.  See Leroy+2019 and Salim+2016 for more details.')
parser.add_argument('--nozband', dest = 'nozband', default = False,action='store_true', help = 'do not use z-band in sed fits.  Default is false. usually you will not adjust this.  adding option for testing to see if this is the root of the systematic difference in magphys results between N and S samples.')
parser.add_argument('--submit', dest = 'submit', default = False, action='store_true',help = 'Set this option to submit the script to slurm. Default is false. You can check output JOB_{scriptid}.sh.  If everything looks good, then add submit flag.')
args = parser.parse_args()


###########################################################
cwd = os.getcwd()




if not args.nozband:
    data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
    script_id = "VFIDall"
    if int(args.ext) == 1:
        data_dir = f"{HOME}/research/Virgo/magphysParallel/output-legacyExt/"
        script_id = "VFIDall-legacyExt"
    if int(args.ext) == 2:
        data_dir = f"{HOME}/research/Virgo/magphysParallel/output-salimExt/"
        script_id = "VFIDall-salimExt"    
elif args.nozband:
    data_dir = f"{HOME}/research/Virgo/magphysParallel/output-nozband/"
    script_id = "VFIDall-nozband"
    if int(args.ext) == 1:
        data_dir = f"{HOME}/research/Virgo/magphysParallel/output-nozband-legacyExt/"
        script_id = "VFIDall-nozband-legacyExt"
    if int(args.ext) == 2:
        data_dir = f"{HOME}/research/Virgo/magphysParallel/output-nozband-salimExt/"
        script_id = "VFIDall-nozband-salimExt"    
print('data_dir = ',data_dir)
print()
print('script_id = ',script_id)
print()
os.chdir(data_dir)

'''
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
'''

outfile = "Dirs.txt"
os.system(f"ls -d ???? > {outfile}")
# count lines
infile = open(outfile,'r')
nfiles = (len(infile.readlines()))
infile.close()



# write out files and submit jobs
#for d in dirlist:

input_file = "Dirs.txt"
write_output(script_id,input_file,narray=nfiles,data_dir=data_dir,submit=args.submit)

os.chdir(cwd)
