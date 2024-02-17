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

os.chdir(data_dir+'/'+galid)

sed_file = '{}.sed'.format(galid)
fit_file = '{}.fit'.format(galid)

if galid.startswith('VFID'):
    sedplot = '{}-magphys-sed.png'.format(d)
    histplot = '{}-magphys-pdfs.png'.format(d)
    sed_pdfs_plot = '{}-magphys-sed-pdfs.png'.format(d)                
else:
    sedplot = 'VFID{}-magphys-sed.png'.format(d)
    histplot = 'VFID{}-magphys-pdfs.png'.format(d)
    sed_pdfs_plot = 'VFID{}-magphys-sed-pdfs.png'.format(d)                                


files = [sedplot,histplot,sed_pdfs_plot]
if os.path.exists(sed_pdfs_plot) & overwrite:
    os.remove(f)
else:
    print("output image exists.  exiting")
    sys.exit()

s = sedFunctions.magphys_sed(sed_file,fit_file,effective_wavelengths)
#s.plot_sed()
#s.plot_histograms()
s.plot_sed_pdfs()
                
#os.rename(sedplot,os.path.join(plotdir,sedplot))
#os.rename(histplot,os.path.join(plotdir,histplot))
#os.rename(sed_pdfs_plot,sed_pdfs_plot) 
plt.close('all')

    
os.chdir(data_dir)
