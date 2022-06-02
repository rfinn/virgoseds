#!/usr/bin/env python
'''
GOAL:
this program takes the output from legacy2magphys.py 
and sets up for parallel processing.  This program creates
a directory for each galaxy with non-nan phot values, 
puts photometry into observations.dat
and creates a symbolic link to the appropriate filter file.

The parallel processing can then be run from 

/home/rfinn/Virgo/magphysParallel

tcsh
source $magphys/runparallel


'''
from astropy.table import Table
import numpy as np
import os

homedir = os.getenv("HOME")
# input files
legdir = homedir+'/research/Virgo/legacy-phot/'

legdir = os.getcwd()
Nphot = os.path.join(legdir,'magphysInputN.dat')
Sphot = os.path.join(legdir,'magphysInputS.dat')
Nfilters = os.path.join(legdir,'legacyFiltersN.dat')
Sfilters = os.path.join(legdir,'legacyFiltersS.dat')


outdir = os.path.join(os.getcwd(),'output/')
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
                pdir = outdir+"{}".format(t[0].replace("VFID","")
            else:
                pdir = outdir+t[0]
            if not(os.path.exists(pdir)):
                os.mkdir(pdir)
            outfile = open(pdir+'/observations.dat','w')
            outfile.write(header)
            outfile.write(line)
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
            if not t[0].startswith('VFID'):
                pdir = outdir+"VFID{:04d}".format(int((t[0])))
            else:
                pdir = outdir+t[0]
            if not(os.path.exists(pdir)):
                os.mkdir(pdir)
            
            outfile = open(pdir+'/observations.dat','w')
            outfile.write(header)
            outfile.write(line)
            outfile.close()
            i += 1

        os.system('cp '+Sfilters+' '+pdir+'/filters.dat')
infile.close()


print('Skipped {} galaxies that had Rmag = nan'.format(nskip))
