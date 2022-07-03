#!/usr/bin/env python


'''
GOAL: 
* grab all the output from magphys
* this creates 
  * the SED and PDF plots
  * a file with one line per vf v2 galaxy, with Mstar and SFR for each galaxy 

'''
import os
import glob
import sys
import numpy as np
from astropy.table import Table
HOME = os.getenv("HOME")
sys.path.append(HOME+'/github/Virgo/programs/')
import sedFunctions
from matplotlib import pyplot as plt
from datetime import datetime
import argparse

testSample = False
makeplots=False


###################################################################
#### SET UP ARGPARSE
###################################################################

parser = argparse.ArgumentParser(description ='Gather output from magphys.')
parser.add_argument('--plot',dest = 'plot', default=False,action='store_true',help='make plots of SED and pdf histograms.  This increases the execution time A LOT!!!  Default is false.')

parser.add_argument('--magdir',dest = 'magdir', default='research/Virgo/magphys/magphysParallelGrawp/output/',help='directory to grab the magphys results from.  the default is HOMEDIR+research/Virgo/maphys/magphysParallelGrawp/output/')
parser.add_argument('--id',dest = 'id', default=None, help='identifier for output table, like nozband or legacyExt or salimExt')

    
args = parser.parse_args()


# the magphys_output directory contains subdirectory for all the galaxy folders

dateTimeObj = datetime.now()
myDate = dateTimeObj.strftime("%d-%b-%Y")


if testSample:
    magphys_output = HOME+'/research/Virgo/magphys/magphysParallelGrawp/output-testsample/'
    output_table_dir = HOME+'/research/Virgo/tables-north/v1/'
    output_table = output_table_dir+'/vf_v1_magphys_testsample_'+myDate+'.fits'
    plotdir = HOME+'/research/Virgo/magphys/magphysParallelGrawp/plots-testsample/'    
else:

    output_table_dir = HOME+'/research/Virgo/tables-north/v2/'
    if args.id is not None:
        output_table = output_table_dir+'/vf_v2_magphys_'+args.id+'_'+myDate+'.fits'
        plotdir = HOME+'/research/Virgo/magphys/magphysParallelGrawp/plots_'+args.id+'/'
    else:
        output_table = output_table_dir+'/vf_v2_magphys_'+myDate+'.fits'
        plotdir = HOME+'/research/Virgo/magphys/magphysParallelGrawp/plots/'

magphys_output = os.path.join(HOME,args.magdir,'')

# check to make sure the plot directory exists
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    
effective_wavelengths = np.array([ 0.1516,0.2267,0.48623,0.64606,0.91993,3.40025,4.65201,12.81034,22.37528],'d')
os.chdir(magphys_output)

dirlist = glob.glob('*')
#print(dirlist)
dirlist.sort()
# get ready to gather sfr and mstar results
vfids = []
vfid_numb = []
sfrs  = []
mstars = []
chisq = []
nmagphys = 0
for d in dirlist:
    if d.startswith('job'):
        continue
    #print('checking directory ',d)
    sedfile = d+'/'+d+'.sed'
    fitfile = d+'/'+d+'.fit'    
    if os.path.exists(sedfile):
        #print('\t found results for ',d)
        lastdir = d
        os.chdir(d)
        nmagphys += 1
        if args.plot:
            s = sedFunctions.magphys_sed(d,effective_wavelengths)
            #s.plot_sed()
            #s.plot_histograms()
            s.plot_sed_pdfs()
            if d.startswith('VFID'):
                sedplot = '{}-magphys-sed.png'.format(d)
                histplot = '{}-magphys-pdfs.png'.format(d)
                sed_pdfs_plot = '{}-magphys-sed-pdfs.png'.format(d)                
            else:
                sedplot = 'VFID{}-magphys-sed.png'.format(d)
                histplot = 'VFID{}-magphys-pdfs.png'.format(d)
                sed_pdfs_plot = 'VFID{}-magphys-sed-pdfs.png'.format(d)                                
                
            #os.rename(sedplot,os.path.join(plotdir,sedplot))
            #os.rename(histplot,os.path.join(plotdir,histplot))
            os.rename(sed_pdfs_plot,os.path.join(plotdir,sed_pdfs_plot)) 
            plt.close('all')
            #sfrs.append(np.log10(s.sfr))
            #mstars.append(np.log10(s.mstar))
            
        # gather outputs only
        # this is much faster
        fit_file = d+'.fit'
        in1 = open(fit_file,'r')
        fit_lines = in1.readlines()
        in1.close()
        # flux is given as L_nu in units of L_sun/Hz
        # I am keeping syntax/variable names similar to plot_sed
        # but this is really a luminosity, or luminosity density?
        flux = np.array(fit_lines[2].split(),'d') # observed L
        e_flux = np.array(fit_lines[3].split(),'d') # error
        # not sure exactly what these are yet, except for z of course...
        i_sfh,i_ir,chi2,z = np.array(fit_lines[8].split(),'d')
        fmu_sfh,fmu_ir,mu,tauv,ssfr,mstar,Ldust,T_W,T_C_ISM,\
            xi_Ctot,xi_PAHtot,xi_MIRtot,xi_Wtot,\
            tvism,Mdust,SFR = np.array(fit_lines[10].split(),'d')
            
        # get sfr and mstars
        sfrs.append(np.log10(SFR))
        mstars.append(np.log10(mstar))
        chisq.append(chi2)        

        vfids.append(d)
        vfid_numb.append(int(d.replace('VFID','')))

        os.chdir(magphys_output)
        

print('max directory = ',lastdir)
print('number processed = ',nmagphys)
sfrs = np.array(sfrs,'d')
mstars = np.array(mstars,'d')
# make table row-lined to
VFID = ['VFID{:04d}'.format(i) for i in range(6780)]
vfid_numb = np.array(vfid_numb,'i')
vfid_sfr = np.zeros(len(VFID),'d')
vfid_mstar = np.zeros(len(VFID),'d')
vfid_chisq = np.zeros(len(VFID),'d')
vfid_sfr[vfid_numb]=sfrs
vfid_mstar[vfid_numb]=mstars
vfid_chisq[vfid_numb]=chisq
tab = Table([VFID,vfid_sfr,vfid_mstar,vfid_chisq],names=['VFID','logSFR','logMstar','chisq'])





tab.write(output_table,format='fits',overwrite=True)

    
