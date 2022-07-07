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
#### FUNCTIONS
###################################################################
def parse_percentile(line):
    '''split percentile line, return med and 68% conf interval  ''' 
    t = line.split()
    #print(t)
    return t[3],t[1],t[4]


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
ssfrs  = []
mstars = []
sfrs_percent  = []
ssfrs_percent  = []
mstars_percent = []
sfrs_med  = []
ssfrs_med  = []
mstars_med = []
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

        
        ######################################################        
        # Get the percentiles of the pdfs for SFR and Mstar
        # percentiles are given in the line after the pdf
        # percentiles: 2.5, 14, 50, 68, 97.5
        ######################################################

        # this code is copied from sedFunctions.py
        # leaving all parameters here in case we want to add additional output
        
        #fmu_SFR = parse_pdf(fit_lines[16:36])
        #fmu_IR = parse_pdf(fit_lines[39:59])
        #mu = parse_pdf(fit_lines[62:82])
        #tau_V = parse_pdf(fit_lines[85:133])

        #sSFR = parse_pdf(fit_lines[136:206])
        t = parse_percentile(fit_lines[207])
        ssfrs_med.append(float(t[0]))        
        ssfrs_percent.append(np.array(t[1:]))

        #Mstar = parse_pdf(fit_lines[209:269])
        t = parse_percentile(fit_lines[270])        
        mstars_med.append(float(t[0]))
        mstars_percent.append(np.array(t[1:]))


        #Ld_tot = parse_pdf(fit_lines[272:332])
        #Tc_ISM = parse_pdf(fit_lines[335:345])
        #Tw_BC = parse_pdf(fit_lines[348:378])
        #xi_C_tot = parse_pdf(fit_lines[381:401])
        #xi_W_tot = parse_pdf(fit_lines[450:470])
        #tau_V_ISM = parse_pdf(fit_lines[473:553])
        #Mdust = parse_pdf(fit_lines[556:616])                

        #SFR = parse_pdf(fit_lines[619:679])
            
        t = parse_percentile(fit_lines[680])
        sfrs_med.append(float(t[0]))                    
        sfrs_percent.append(t[1:])
            

            
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
ssfrs = np.array(sfrs,'d')

sfrs_med = np.array(sfrs_med,'d')
mstars_med = np.array(mstars_med,'d')
ssfrs_med = np.array(sfrs_med,'d')

sfrs_perc = np.array(sfrs_percent,'d')
mstars_perc = np.array(mstars_percent,'d')
ssfrs_perc = np.array(sfrs_percent,'d')

# make table row-lined to
VFID = ['VFID{:04d}'.format(i) for i in range(6780)]
vfid_numb = np.array(vfid_numb,'i')
vfid_ssfr = np.zeros(len(VFID),'d')
vfid_sfr = np.zeros(len(VFID),'d')
vfid_mstar = np.zeros(len(VFID),'d')

magphys_flag = np.zeros(len(VFID),'bool')

vfid_ssfr_percent = np.zeros((len(VFID),2),'d')
vfid_sfr_percent = np.zeros((len(VFID),2),'d')
vfid_mstar_percent = np.zeros((len(VFID),2),'d')

vfid_ssfr_med = np.zeros(len(VFID),'d')
vfid_sfr_med = np.zeros(len(VFID),'d')
vfid_mstar_med = np.zeros(len(VFID),'d')

vfid_chisq = np.zeros(len(VFID),'d')

vfid_sfr[vfid_numb]=sfrs
vfid_ssfr[vfid_numb]=ssfrs
vfid_mstar[vfid_numb]=mstars
vfid_chisq[vfid_numb]=chisq

magphys_flag[vfid_numb] = np.ones(len(vfid_numb),'bool')

vfid_ssfr_med[vfid_numb]=ssfrs_med
vfid_sfr_med[vfid_numb]=sfrs_med
vfid_mstar_med[vfid_numb]=mstars_med


vfid_ssfr_percent[vfid_numb]=ssfrs_perc
vfid_sfr_percent[vfid_numb]=sfrs_perc
vfid_mstar_percent[vfid_numb]=mstars_perc



data_columns = [VFID,\
                vfid_sfr,vfid_sfr_med, vfid_sfr_percent,\
                vfid_mstar,vfid_mstar_med,vfid_mstar_percent,\
                vfid_ssfr,vfid_ssfr_med, vfid_ssfr_percent,\
                vfid_chisq,magphys_flag]

names=['VFID',
       'logSFR','logSFR-med','logSFR-68conf',\
       'logMstar','logMstar-med','logMstar-68conf',\
       'logsSFR','logsSFR-med','logsSFR-68conf',\
       'chisq','magphysFlag']

tab = Table(data=data_columns,names=names)





tab.write(output_table,format='fits',overwrite=True)

    
