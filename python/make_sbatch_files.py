#!/usr/bin/env python

'''


* run this from the top level directory, where all the VFIDs are in subdirectories
* the log files will go in log/
* the script files will go in the top level directory

'''


import os
import subprocess
import sys
import glob
import time


HOME = os.getenv("HOME")


def write_output(filename,dirname,submit=False):
    ''' copying from Matt Bellis '''
    output = ""
    output += "#!/bin/bash\n"
    output += "\n"
    output += "# Set Job Name\n"
    output += "#SBATCH -J job\n"
    output += "\n"
    output += "# Set file to capture standard out and standard error and append the jobID (%j)\n"
    output += "#SBATCH -o logs/job.out.%j\n"
    output += "\n"
    output += "#SBATCH --partition=high\n"
    
    output += "#Set the number of nodes\n"
    output += "#SBATCH -N 1\n"
    output += "#SBATCH --ntasks=1\n"
    output += "\n"
    output += "#Set the time limit for the job - 10 hour is specified here\n"
    output += "#SBATCH --time=5:00:00\n"
    output += "\n"
    output += "#SBATCH --cpus-per-task=1\n"
    output += "\n"
    output += "# Load any environmental modules needed\n"
    output += "module load Python3\n"
    output += "module load gnu9\n"
    output += "\n"
    #output += "# Move to the directory needed - defaults to the submission directory\n"
    #output += "cd $HOME/research/Virgo/magphysParallel/output/\n"
    output += "\n"
    output += "# perform calculation\n"
    output += "\n"
    output += f"python {HOME}/github/virgoseds/python/run1magphys.py {dirname}\n"
    # adding test to see if it makes the full number of jobs
    # this is basically a "hello world!" program
    #output += f"python {HOME}/github/virgoseds/python/grawpTest.py {dirname}\n"

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
        #time.sleep(0.1)
        

###########################################################
cwd = os.getcwd()
#data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
#data_dir = f"{HOME}/research/Virgo/magphysParallel/damien_input_files/output/"
data_dir = os.getcwd()
if not os.path.exists(data_dir+'/logs'):
    os.mkdir(data_dir+'/logs')
dirlist = glob.glob(f"{data_dir}/????")
# reversing list b/c later galaxies are not getting run
dirlist.sort(reverse=True)


nrun=0
#for d in dirlist[0:5]:    
for d in dirlist:

    # remove full path to directory so just VFID???? is passed in
    gname = os.path.basename(d)
    
    os.chdir(data_dir)
    # check to see if magphys results exist    
    if os.path.exists(f"{gname}/{gname}.fit") & os.path.exists(f"{gname}/{gname}.sed"):
            #os.rename(f"{vfid}", f"{done_dir}/{vfid}")
            print(f"output exists for {gname} - moving to next galaxy")
            continue

    write_output(gname,gname,submit=True)
    nrun += 1
print('number of jobs to run = {}'.format(nrun))
# write out files and submit jobs
