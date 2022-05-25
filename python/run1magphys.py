#!/usr/bin/env python

'''
This should be run from a directory that contains

filters.dat
observations.dat

I will put each galaxy in its own folder, and the parallel program 
will move to each folder and run this script

time parallel -j+0 --eta ‘cd {} && python $magphys/run1magphys.py' ::: */ 

'''
import os
import sys
import glob

# galaxy to work with
galid = sys.argv[1]
HOME = os.getenv("HOME")
# directory where galaxy folders are
data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
data_dir = os.getcwd()+'/'
code_dir = f"{HOME}/software/magphys/"
#print(data_dir)
#print(code_dir)

# move to data directory
os.chdir(data_dir)

# check to see if the 
if os.path.exists(f"{galid}/{galid}.fit"):
    # exist this script
    print(f"output exists for {galid}")
    print("moving to the next galaxy")
    # make sure lbr files have been removed
    #os.system('rm *.lbr')
    exit()


os.chdir(data_dir+'/'+galid)
print(os.getcwd())
print('available files: ',glob.glob('*.dat'))
os.system(f"source {code_dir}/.magphys_bashrc")
os.system(f"{code_dir}/make_zgrid")

# could check to see if lbr files exist
# and skip this if they do
os.system(f"source {code_dir}/get_libs_bash")

print('running fit_sed')
os.system(f"{code_dir}/fit_sed")
#os.system('rm *.lbr')
    
os.chdir(data_dir)
