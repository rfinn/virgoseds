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
#import glob

# galaxy to work with
galid = sys.argv[1]
HOME = os.getenv("HOME")
# directory where galaxy folders are
data_dir = '{}/research/Virgo/magphysParallel/output/'.format(HOME)
code_dir = '{}/software/magphys/'.format(HOME)
#print(data_dir)
#print(code_dir)

# move to data directory
os.chdir(data_dir)

# check to see if the 
if os.path.exists('{galid}/{galid}.fit'):
    # exist this script
    print('output exists for {galid}')
    print('moving to the next galaxy')
    # make sure lbr files have been removed
    os.system('rm *.lbr')
    exit()


os.chdir(data_dir+'/'+galid)
#print(os.getcwd())
#print('available files: ',glob.glob('*.dat'))
s = "source {}/.magphys_bashrc".format(code_dir)
os.system(s)
s = '{}/make_zgrid'.format(code_dir)
os.system(s)

# could check to see if lbr files exist
# and skip this if they do
os.system('source {}/get_libs_bash'.format(code_dir))

os.system('{}/fit_sed'.format(code_dir))
os.system('rm *.lbr')
    
os.chdir(data_dir)
