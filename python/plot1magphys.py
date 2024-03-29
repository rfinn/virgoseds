#!/usr/bin/env python

'''
This should be run from a directory that contains subfolders for each galaxy

like


filters.dat
observations.dat

I will put each galaxy in its own folder, and the parallel program 
will move to each folder and run this script
 


time parallel -j+0 --eta ‘cd {} && python $magphys/run1magphys.py' ::: */ 

'''
import os
import sys
import glob
import numpy as np
# galaxy to work with
galid = sys.argv[1]
HOME = os.getenv("HOME")
# directory where galaxy folders are
#data_dir = f"{HOME}/research/Virgo/magphysParallel/output/"
data_dir = os.getcwd()+'/'
code_dir = f"{HOME}/github/virgoseds/python/"
#print(data_dir)
sys.path.append(code_dir)

import sedFunctions
# move to data directory
#os.chdir(data_dir)

effective_wavelengths = np.array([ 0.1516,0.2267,0.48623,0.64606,0.91993,3.40025,4.65201,12.81034,22.37528],'d')

overwrite = True
os.chdir(data_dir+'/'+galid)

sed_file = '{}.sed'.format(galid)
fit_file = '{}.fit'.format(galid)

if galid.startswith('VFID'):
    sedplot = '{}-magphys-sed.png'.format(galid)
    histplot = '{}-magphys-pdfs.png'.format(galid)
    sed_pdfs_plot = '{}-magphys-sed-pdfs.png'.format(galid)                
else:
    sedplot = 'VFID{}-magphys-sed.png'.format(galid)
    histplot = 'VFID{}-magphys-pdfs.png'.format(galid)
    sed_pdfs_plot = 'VFID{}-magphys-sed-pdfs.png'.format(galid)                                


files = [sedplot,histplot,sed_pdfs_plot]
if os.path.exists(sed_pdfs_plot):
    if overwrite:
        os.remove(sed_pdfs_plot)
    else:
        print("output found - exiting")
        sys.exit()
else:
    print("making plot")


s = sedFunctions.magphys_sed(sed_file,fit_file,effective_wavelengths,galid=galid)
#s.plot_sed()
#s.plot_histograms()
s.plot_sed_pdfs()
                
#os.rename(sedplot,os.path.join(plotdir,sedplot))
#os.rename(histplot,os.path.join(plotdir,histplot))
#os.rename(sed_pdfs_plot,sed_pdfs_plot) 
#plt.close('all')

    
os.chdir(data_dir)
