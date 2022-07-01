#!/usr/bin/env python
'''
GOAL:
this program takes the output from legacy2magphys.py 
and sets up for parallel processing.  

This program:
* creates a directory for each galaxy with non-nan phot values, 
* puts photometry into observations.dat.
* copies the appropriate filter file.

The parallel processing can then be run by calling make_sbatch_array_files.py



'''
from astropy.table import Table
import numpy as np
import os




###########################
##### SET UP ARGPARSE
###########################
import argparse
parser = argparse.ArgumentParser(description ='Program to set up directories for individual galaxies.  Can then run magphys in parallel using make_sbatch_array_files.py.')
#parser.add_argument('--sbmag', dest = 'sbmag', default = 24, help = 'sb to fit.  default is 24, as this has the highest fraction of non-zero entries.  The options are [22,26] in increments of 0.5.')
parser.add_argument('--ext', dest = 'ext', default = 0, help = 'extinction correction to apply.  0=None; 1=Legacy Survey; 2=Salim/Leroy.  Default is zero. The main difference between 1 and 2 is how the GALEX fluxes are handled.  See Leroy+2019 and Salim+2016 for more details.')
parser.add_argument('--nozband', dest = 'nozband', default = False,action='store_true', help = 'do not use z-band in sed fits.  Default is false. usually you will not adjust this.  adding option for testing to see if this is the root of the systematic difference in magphys results between N and S samples.')


args = parser.parse_args()
############################


###########################
##### SET UP INPUT FILES
##### ACCORDING TO EXTINCTION
##### AND ZBAND OPTIONS
###########################

legdir = homedir+'/research/Virgo/legacy-phot/'

legdir = os.getcwd()
file_suffix=''

# these are the default values for running
# without the extinction correction and
# including the z-band data
outdir = os.path.join(os.getcwd(),'output/')
Nfilters = os.path.join(legdir,'legacyFiltersN.dat')
Sfilters = os.path.join(legdir,'legacyFiltersS.dat')

if int(args.ext) == 1:
    outdir = outdir.replace('output/',"output-legacyExt/")
    file_suffix='-legacyExt'
if int(args.ext) == 2:
    outdir = outdir.replace('output/',"output-salimExt/")
    file_suffix='-salimExt'    
if args.nozband:
    outdir = outdir.replace('output/',"output-nozband/")
    Nfilters = os.path.join(legdir,'legacyFiltersN-nozband.dat')
    Sfilters = os.path.join(legdir,'legacyFiltersS-nozband.dat')
    

Nphot = os.path.join(legdir,f'magphysInputN{file_suffix}.dat')
Sphot = os.path.join(legdir,f'magphysInputS{file_suffix}.dat')





homedir = os.getenv("HOME")
# input files


#homedir+'/research/Virgo/magphysParallel/'
#outdir = legdir+'/maphysParallel/'
if not(os.path.exists(outdir)):
    os.mkdir(outdir)

# loop through North file
infile = open(Nphot,'r')
i = 0
nskip = 0
for line in infile:
    if i == 0:
        header = line
        i += 1
        continue
    else:
        t = line.split()
        if np.isnan(float(t[8])):
            #print('skipping ',t[0])
            nskip += 1
            continue
        else:
            # adding this for compatibility with the files that I gave Damien
            # these just have the numeric ID, and also use the V1 catalog names...
            if t[0].startswith('VFID'):
                pdir = outdir+"{}".format(t[0].replace("VFID",""))
            else:
                pdir = outdir+t[0]
            if not(os.path.exists(pdir)):
                os.mkdir(pdir)
            outfile = open(pdir+'/observations.dat','w')
            outfile.write(header)
            outfile.write(line.replace("VFID",""))
            outfile.close()
            i += 1

        os.system('cp '+Nfilters+' '+pdir+'/filters.dat')
infile.close()

# loop through South file
infile = open(Sphot,'r')
i = 0
for line in infile:
    if i == 0:
        header = line
        i += 1
        continue
    else:
        t = line.split()
        if np.isnan(float(t[8])):
            #print('skipping ',t[0])
            nskip += 1
            continue
        else:
            if t[0].startswith('VFID'):
                pdir = outdir+"{}".format(t[0].replace("VFID",""))
            else:
                pdir = outdir+t[0]
            if not(os.path.exists(pdir)):
                os.mkdir(pdir)
            
            outfile = open(pdir+'/observations.dat','w')
            outfile.write(header)
            outfile.write(line.replace("VFID",""))
            outfile.close()
            i += 1

        os.system('cp '+Sfilters+' '+pdir+'/filters.dat')
infile.close()


print('Skipped {} galaxies that had Rmag = nan'.format(nskip))
