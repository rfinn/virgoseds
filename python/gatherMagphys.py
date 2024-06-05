#!/usr/bin/env python


'''
GOAL: 
* grab all the output from magphys
* this creates 
  * the SED and PDF plots
  * a file with one line per vf v2 galaxy, with Mstar and SFR for each galaxy 

* the output table is written to ~/research/Virgo/tables-north/v2/




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
    #return float(t[3]),float(t[1]),float(t[4])
    return float(t[2]),float(t[1]),float(t[3])


def calc_percentiles(lines,plotflag=False):
    '''sum total percentile, calc and return med and 68% conf interval  '''
    from scipy.interpolate import interp1d
    
    val = []
    prob = []
    for l in lines:
        t = l.split()
        val.append(float(t[0]))
        prob.append(float(t[1]))        
        #try:
        #    val.append(float(t[0]))
        #    prob.append(float(t[0]))
        #except ValueError:
        #    print("error splitting line: ",l)
        #    sys.exit()

    val = np.array(val,'d')
    prob = np.array(prob,'d')
    
    # calculate the total probability
    ptotal = np.sum(prob)

    # calculate cumulative sum and create inter function
    integral2 = np.cumsum(prob)
    interp_profile = interp1d(integral2,val)
    
    # 0.025,0.16,0.50,0.84,0.975/
    thresholds = [0.16,0.5,0.84]
    percentiles = []
    for p in thresholds:
        try:
            t = interp_profile(p*integral2[-1])
            percentiles.append(float(f"{t:.3f}"))
        except ValueError:
            percentiles.append(np.nan)

    if plotflag:
        plt.figure()
        plt.subplot(1,2,1)
        plt.plot(val,prob,'b.')
        plt.subplot(1,2,2)
        plt.plot(val,integral2)
        for p in percentiles:
            plt.axvline(x=p,ls='--')
    return percentiles[1], percentiles[0], percentiles[2]
    


###################################################################
#### SET UP ARGPARSE
###################################################################

parser = argparse.ArgumentParser(description ='Gather output from magphys. Write output table to ~/research/Virgo/tables-north/v2/')
parser.add_argument('--plot',dest = 'plot', default=False,action='store_true',help='make plots of SED and pdf histograms.  This increases the execution time A LOT!!!  Default is false.')
parser.add_argument('--verbose',dest = 'verbose', default=False,action='store_true',help='set this to print extra statements')

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


# save the parameters from the best-fit models
sfr_best  = []
ssfr_best  = []
mstar_best = []
fmuSFR_best = []
fmuIR_best = []
mu_best = []
tau_V_best = []
Ld_tot_best = []
Tc_ISM_best = []
Tw_BC_best = []
xi_C_tot_best = []
xi_PAH_tot_best = []
xi_MIR_tot_best = []
xi_W_tot_best = []
tau_V_ISM_best = []
Mdust_best = []
chisq_best = []


# parameters from pdf
sfr_med  = []
ssfr_med  = []
mstar_med = []

# additional parameters that we haven't been tracking
fmuSFR_med = []
fmuIR_med = []
mu_med = []
tau_V_med = []
Ld_tot_med = []
Tc_ISM_med = []
Tw_BC_med = []
xi_C_tot_med = []
xi_PAH_tot_med = []
xi_MIR_tot_med = []
xi_W_tot_med = []
tau_V_ISM_med = []
Mdust_med = []


# store the 68% conf interval from the pdf
sfr_percent  = []
ssfr_percent  = []
mstar_percent = []
            
fmuSFR_percent = []
fmuIR_percent = []
mu_percent = []
tau_V_percent = []
Ld_tot_percent = []
Tc_ISM_percent = []
Tw_BC_percent = []
xi_C_tot_percent = []
xi_PAH_tot_percent = []
xi_MIR_tot_percent = []
xi_W_tot_percent = []
tau_V_ISM_percent = []
Mdust_percent = []

            
# save the observed fluxes and errors
obs_flux = []
obs_flux_err = []
best_flux = []



nmagphys = 0
for d in dirlist:
    
    if d.startswith('job'):
        continue
    if not os.path.isdir(d):
        continue
    if args.verbose:
        print('checking directory ',d)
        print(os.getcwd())
    
    sedfile = d+'/'+d+'.sed'
    fitfile = d+'/'+d+'.fit'
    #print(d)
    if os.path.exists(fitfile):
        #print('\t found results for ',d)
        lastdir = d
        os.chdir(d)
        nmagphys += 1
        if args.plot:
            sed_file = '{}.sed'.format(d)
            fit_file = '{}.fit'.format(d)

            s = sedFunctions.magphys_sed(sed_file,fit_file,effective_wavelengths)
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

        # keep all fluxes
        obs_flux.append(flux)
        obs_flux_err.append(e_flux)

        # this is the output from the best-fit model
        # not sure exactly what these are yet, except for z of course...
        i_sfh,i_ir,chi2,z = np.array(fit_lines[8].split(),'d')
        fmu_sfh,fmu_ir,mu,tauv,ssfr,mstar,Ldust,T_W,T_C_ISM,\
            xi_Ctot,xi_PAHtot,xi_MIRtot,xi_Wtot,\
            tvism,Mdust,SFR = np.array(fit_lines[10].split(),'d')

        # best flux
        best_flux.append(np.array(fit_lines[12].split(),'d')) # best L

        
        # get sfr and mstars
        sfr_best.append(np.log10(SFR))
        mstar_best.append(np.log10(mstar))

        # not sure if this is already logged or not
        # it is not logged, so need to log it
        ssfr_best.append(np.log10(ssfr))

        fmuSFR_best.append(fmu_sfh)
        fmuIR_best.append(fmu_ir)
        mu_best.append(mu)
        tau_V_best.append(tauv)
        Ld_tot_best.append(Ldust)
        Tc_ISM_best.append(T_C_ISM)
        Tw_BC_best.append(T_W)
        xi_C_tot_best.append(xi_Ctot)
        xi_PAH_tot_best.append(xi_PAHtot)
        xi_MIR_tot_best.append(xi_MIRtot)        
        xi_W_tot_best.append(xi_Wtot)
        tau_V_ISM_best.append(tvism)
        Mdust_best.append(Mdust)
        
        chisq_best.append(chi2)        
        
        
        ######################################################        
        # Get the percentiles of the pdfs for SFR and Mstar
        # percentiles are given in the line after the pdf
        # percentiles: 2.5, 14, 50, 86, 97.5
        ######################################################

        # this code is copied from sedFunctions.py
        # leaving all parameters here in case we want to add additional output
        
        #fmu_SFR = parse_pdf(fit_lines[16:36])
        #t = parse_percentile(fit_lines[37])
        t = calc_percentiles(fit_lines[16:36])
        fmuSFR_med.append(float(t[0]))        
        fmuSFR_percent.append(np.array(t[1:]))
        
        #fmu_IR = parse_pdf(fit_lines[39:59])
        #t = parse_percentile(fit_lines[60])
        t = calc_percentiles(fit_lines[39:59])        
        fmuIR_med.append(float(t[0]))        
        fmuIR_percent.append(np.array(t[1:]))

        #mu = parse_pdf(fit_lines[62:82])
        t = parse_percentile(fit_lines[83])
        t = calc_percentiles(fit_lines[62:82])        
        mu_med.append(float(t[0]))        
        mu_percent.append(np.array(t[1:]))
        
        #tau_V = parse_pdf(fit_lines[85:133])
        #t = parse_percentile(fit_lines[134])
        t = calc_percentiles(fit_lines[85:133])        
        tau_V_med.append(float(t[0]))        
        tau_V_percent.append(np.array(t[1:]))

        #sSFR = parse_pdf(fit_lines[136:206])
        #t = parse_percentile(fit_lines[207])
        t = calc_percentiles(fit_lines[136:206])        
        ssfr_med.append(float(t[0]))
        ssfr_percent.append(np.array(t[1:]))

        #Mstar = parse_pdf(fit_lines[209:269])
        #t = parse_percentile(fit_lines[270])
        #t = parse_percentile(fit_lines[310])
        t = calc_percentiles(fit_lines[209:309])
        mstar_med.append(float(t[0]))
        mstar_percent.append(np.array(t[1:]))

        #Ld_tot = parse_pdf(fit_lines[272:332])
        #t = parse_percentile(fit_lines[333])
        #t = parse_percentile(fit_lines[413])
        t = calc_percentiles(fit_lines[312:412])
        Ld_tot_med.append(float(t[0]))
        Ld_tot_percent.append(np.array(t[1:]))
        
        #Tc_ISM = parse_pdf(fit_lines[335:345])
        #t = parse_percentile(fit_lines[346])
        #t = parse_percentile(fit_lines[426])
        t = calc_percentiles(fit_lines[415:425])        
        Tc_ISM_med.append(float(t[0]))
        Tc_ISM_percent.append(np.array(t[1:]))

        
        #Tw_BC = parse_pdf(fit_lines[348:378])
        #t = parse_percentile(fit_lines[379])
        #t = parse_percentile(fit_lines[459])
        t = calc_percentiles(fit_lines[428:458])        
        Tw_BC_med.append(float(t[0]))
        Tw_BC_percent.append(np.array(t[1:]))
        
        #xi_C_tot = parse_pdf(fit_lines[381:401])
        #t = parse_percentile(fit_lines[402])
        #t = parse_percentile(fit_lines[482])
        t = calc_percentiles(fit_lines[461:481])
        xi_C_tot_med.append(float(t[0]))
        xi_C_tot_percent.append(np.array(t[1:]))

        # PAH
        #t = parse_percentile(fit_lines[425])
        #t = parse_percentile(fit_lines[505])
        t = calc_percentiles(fit_lines[484:504])
        xi_PAH_tot_med.append(float(t[0]))
        xi_PAH_tot_percent.append(np.array(t[1:]))

        # MIR - hadn't been tracking these...
        #t = parse_percentile(fit_lines[448])
        #t = parse_percentile(fit_lines[528])
        t = calc_percentiles(fit_lines[507:527])
        xi_MIR_tot_med.append(float(t[0]))
        xi_MIR_tot_percent.append(np.array(t[1:]))
        
        #xi_W_tot = parse_pdf(fit_lines[450:470])
        #t = parse_percentile(fit_lines[471])
        #t = parse_percentile(fit_lines[551])
        t = calc_percentiles(fit_lines[530:550])        
        xi_W_tot_med.append(float(t[0]))
        xi_W_tot_percent.append(np.array(t[1:]))
        
        #tau_V_ISM = parse_pdf(fit_lines[473:553])
        #t = parse_percentile(fit_lines[554])
        #t = parse_percentile(fit_lines[634])
        t = calc_percentiles(fit_lines[553:633])        
        tau_V_ISM_med.append(float(t[0]))
        tau_V_ISM_percent.append(np.array(t[1:]))
        
        #Mdust = parse_pdf(fit_lines[556:616])                
        #t = parse_percentile(fit_lines[617])
        #t = parse_percentile(fit_lines[697])
        t = calc_percentiles(fit_lines[636:696])
        Mdust_med.append(float(t[0]))
        Mdust_percent.append(np.array(t[1:]))

        #SFR = parse_pdf(fit_lines[619:679]) 
        #t = parse_percentile(fit_lines[680])
        t = parse_percentile(fit_lines[810])
        if args.verbose:
            print("\nSFR parsing percentile line: ",t)
        t = calc_percentiles(fit_lines[699:809],plotflag=True)
        if args.verbose:
            print("SFR calculating percentile line: ",t)        
        sfr_med.append(float(t[0]))
        sfr_percent.append(np.array(t[1:]))

        vfids.append(d)
        vfid_numb.append(int(d.replace('VFID','')))

        os.chdir(magphys_output)
        

print('max directory = ',lastdir)
print('number processed = ',nmagphys)

## ugh, this is terrible.  I'm sure there is a smarter way to write out a million columns...

# make table row-matched v2 tables
VFID = ['VFID{:04d}'.format(i) for i in range(6780)]
magphys_flag = np.zeros(len(VFID),'bool')

vfid_numb = np.array(vfid_numb,'i')

# best-fit parameters
vfid_ssfr_best = np.zeros(len(VFID),'d')
vfid_sfr_best = np.zeros(len(VFID),'d')
vfid_mstar_best = np.zeros(len(VFID),'d')
vfid_fmuSFR_best = np.zeros(len(VFID),'d')
vfid_fmuIR_best = np.zeros(len(VFID),'d')
vfid_mu_best = np.zeros(len(VFID),'d')
vfid_tau_V_best = np.zeros(len(VFID),'d')
vfid_Ld_tot_best = np.zeros(len(VFID),'d')
vfid_Tc_ISM_best = np.zeros(len(VFID),'d')
vfid_Tw_BC_best = np.zeros(len(VFID),'d')
vfid_xi_C_tot_best = np.zeros(len(VFID),'d')
vfid_xi_PAH_tot_best = np.zeros(len(VFID),'d')
vfid_xi_MIR_tot_best = np.zeros(len(VFID),'d')
vfid_xi_W_tot_best = np.zeros(len(VFID),'d')
vfid_tau_V_ISM_best = np.zeros(len(VFID),'d')
vfid_Mdust_best = np.zeros(len(VFID),'d')
vfid_chisq_best = np.zeros(len(VFID),'d')


vfid_sfr_best[vfid_numb] = np.array(sfr_best,'d')
vfid_ssfr_best[vfid_numb] = np.array(ssfr_best,'d')
vfid_mstar_best[vfid_numb] = np.array(mstar_best,'d')
vfid_chisq_best[vfid_numb] = np.array(chisq_best,'d')
vfid_fmuSFR_best[vfid_numb] = np.array(fmuSFR_best,'d')
vfid_fmuIR_best[vfid_numb] = np.array(fmuIR_best,'d')
vfid_mu_best[vfid_numb] = np.array(mu_best,'d')
vfid_tau_V_best[vfid_numb] = np.array(tau_V_best,'d')
vfid_Ld_tot_best[vfid_numb] = np.array(Ld_tot_best,'d')
vfid_Tc_ISM_best[vfid_numb] = np.array(Tc_ISM_best,'d')
vfid_Tw_BC_best[vfid_numb] = np.array(Tw_BC_best,'d')
vfid_xi_C_tot_best[vfid_numb] = np.array(xi_C_tot_best,'d')
vfid_xi_PAH_tot_best[vfid_numb] = np.array(xi_PAH_tot_best,'d')
vfid_xi_MIR_tot_best[vfid_numb] = np.array(xi_MIR_tot_best,'d')
vfid_xi_W_tot_best[vfid_numb] = np.array(xi_W_tot_best,'d')
vfid_tau_V_ISM_best[vfid_numb] = np.array(tau_V_ISM_best,'d')
vfid_Mdust_best[vfid_numb] = np.array(Mdust_best,'d')


# median parameters

vfid_ssfr_med = np.zeros(len(VFID),'d')
vfid_sfr_med = np.zeros(len(VFID),'d')
vfid_mstar_med = np.zeros(len(VFID),'d')
vfid_fmuSFR_med = np.zeros(len(VFID),'d')
vfid_fmuIR_med = np.zeros(len(VFID),'d')
vfid_mu_med = np.zeros(len(VFID),'d')
vfid_tau_V_med = np.zeros(len(VFID),'d')
vfid_Ld_tot_med = np.zeros(len(VFID),'d')
vfid_Tc_ISM_med = np.zeros(len(VFID),'d')
vfid_Tw_BC_med = np.zeros(len(VFID),'d')
vfid_xi_C_tot_med = np.zeros(len(VFID),'d')
vfid_xi_PAH_tot_med = np.zeros(len(VFID),'d')
vfid_xi_MIR_tot_med = np.zeros(len(VFID),'d')
vfid_xi_W_tot_med = np.zeros(len(VFID),'d')
vfid_tau_V_ISM_med = np.zeros(len(VFID),'d')
vfid_Mdust_med = np.zeros(len(VFID),'d')



vfid_sfr_med[vfid_numb] = np.array(sfr_med,'d')
vfid_ssfr_med[vfid_numb] = np.array(ssfr_med,'d')
vfid_mstar_med[vfid_numb] = np.array(mstar_med,'d')
vfid_fmuSFR_med[vfid_numb] = np.array(fmuSFR_med,'d')
vfid_fmuIR_med[vfid_numb] = np.array(fmuIR_med,'d')
vfid_mu_med[vfid_numb] = np.array(mu_med,'d')
vfid_tau_V_med[vfid_numb] = np.array(tau_V_med,'d')
vfid_Ld_tot_med[vfid_numb] = np.array(Ld_tot_med,'d')
vfid_Tc_ISM_med[vfid_numb] = np.array(Tc_ISM_med,'d')
vfid_Tw_BC_med[vfid_numb] = np.array(Tw_BC_med,'d')
vfid_xi_C_tot_med[vfid_numb] = np.array(xi_C_tot_med,'d')
vfid_xi_PAH_tot_med[vfid_numb] = np.array(xi_PAH_tot_med,'d')
vfid_xi_MIR_tot_med[vfid_numb] = np.array(xi_MIR_tot_med,'d')
vfid_xi_W_tot_med[vfid_numb] = np.array(xi_W_tot_med,'d')
vfid_tau_V_ISM_med[vfid_numb] = np.array(tau_V_ISM_med,'d')
vfid_Mdust_med[vfid_numb] = np.array(Mdust_med,'d')



# median percent

vfid_ssfr_percent = np.zeros((len(VFID),2),'d')
vfid_sfr_percent = np.zeros((len(VFID),2),'d')
vfid_mstar_percent = np.zeros((len(VFID),2),'d')
vfid_fmuSFR_percent = np.zeros((len(VFID),2),'d')
vfid_fmuIR_percent = np.zeros((len(VFID),2),'d')
vfid_mu_percent = np.zeros((len(VFID),2),'d')
vfid_tau_V_percent = np.zeros((len(VFID),2),'d')
vfid_Ld_tot_percent = np.zeros((len(VFID),2),'d')
vfid_Tc_ISM_percent = np.zeros((len(VFID),2),'d')
vfid_Tw_BC_percent = np.zeros((len(VFID),2),'d')
vfid_xi_C_tot_percent = np.zeros((len(VFID),2),'d')
vfid_xi_PAH_tot_percent = np.zeros((len(VFID),2),'d')
vfid_xi_MIR_tot_percent = np.zeros((len(VFID),2),'d')
vfid_xi_W_tot_percent = np.zeros((len(VFID),2),'d')
vfid_tau_V_ISM_percent = np.zeros((len(VFID),2),'d')
vfid_Mdust_percent = np.zeros((len(VFID),2),'d')



vfid_sfr_percent[vfid_numb] = np.array(sfr_percent,'d')
vfid_ssfr_percent[vfid_numb] = np.array(ssfr_percent,'d')
vfid_mstar_percent[vfid_numb] = np.array(mstar_percent,'d')
vfid_fmuSFR_percent[vfid_numb] = np.array(fmuSFR_percent,'d')
vfid_fmuIR_percent[vfid_numb] = np.array(fmuIR_percent,'d')
vfid_mu_percent[vfid_numb] = np.array(mu_percent,'d')
vfid_tau_V_percent[vfid_numb] = np.array(tau_V_percent,'d')
vfid_Ld_tot_percent[vfid_numb] = np.array(Ld_tot_percent,'d')
vfid_Tc_ISM_percent[vfid_numb] = np.array(Tc_ISM_percent,'d')
vfid_Tw_BC_percent[vfid_numb] = np.array(Tw_BC_percent,'d')
vfid_xi_C_tot_percent[vfid_numb] = np.array(xi_C_tot_percent,'d')
vfid_xi_PAH_tot_percent[vfid_numb] = np.array(xi_PAH_tot_percent,'d')
vfid_xi_MIR_tot_percent[vfid_numb] = np.array(xi_MIR_tot_percent,'d')
vfid_xi_W_tot_percent[vfid_numb] = np.array(xi_W_tot_percent,'d')
vfid_tau_V_ISM_percent[vfid_numb] = np.array(tau_V_ISM_percent,'d')
vfid_Mdust_percent[vfid_numb] = np.array(Mdust_percent,'d')












magphys_flag[vfid_numb] = np.ones(len(vfid_numb),'bool')




data_columns = [VFID,\
                vfid_mstar_med,vfid_sfr_med,vfid_ssfr_med,\
                vfid_fmuSFR_med, vfid_fmuIR_med, vfid_mu_med,\
                vfid_tau_V_med, vfid_Ld_tot_med,\
                vfid_Tc_ISM_med, vfid_Tw_BC_med,\
                vfid_xi_C_tot_med, vfid_xi_PAH_tot_med,\
                vfid_xi_MIR_tot_med, vfid_xi_W_tot_med,\
                vfid_tau_V_ISM_med,vfid_Mdust_med,\
                
                vfid_mstar_percent,vfid_sfr_percent,vfid_ssfr_percent,\
                vfid_fmuSFR_percent, vfid_fmuIR_percent, vfid_mu_percent,\
                vfid_tau_V_percent, vfid_Ld_tot_percent,\
                vfid_Tc_ISM_percent, vfid_Tw_BC_percent,\
                vfid_xi_C_tot_percent, vfid_xi_PAH_tot_percent,\
                vfid_xi_MIR_tot_percent, vfid_xi_W_tot_percent,\
                vfid_tau_V_ISM_percent,vfid_Mdust_percent,\
                
                vfid_mstar_best,vfid_sfr_best,vfid_ssfr_best,\
                vfid_fmuSFR_best, vfid_fmuIR_best, vfid_mu_best,\
                vfid_tau_V_best, vfid_Ld_tot_best,\
                vfid_Tc_ISM_best, vfid_Tw_BC_best,\
                vfid_xi_C_tot_best, vfid_xi_PAH_tot_best,\
                vfid_xi_MIR_tot_best, vfid_xi_W_tot_best,\
                vfid_tau_V_ISM_best,vfid_Mdust_best,\
                vfid_chisq_best,magphys_flag]


names=['VFID',
       'logMstar_med','logSFR_med','logsSFR_med',\
       'fmu_SFR_med','fmu_IR_med','mu_med',\
       'tau_V_med','Ldust_tot_med',\
       'Tc_ISM_med','Tw_BC_med',\
       'xi_C_tot_med','xi_PAH_tot_med','xi_MIR_tot_med','xi_W_tot_med',\
       'tau_V_ISM_med','Mdust_med',\
       
       'logMstar_percent','logSFR_percent','logsSFR_percent',\
       'fmu_SFR_percent','fmu_IR_percent','mu_percent',\
       'tau_V_percent','Ldust_tot_percent',\
       'Tc_ISM_percent','Tw_BC_percent',\
       'xi_C_tot_percent','xi_PAH_tot_percent','xi_MIR_tot_percent','xi_W_tot_percent',\
       'tau_V_ISM_percent','Mdust_percent',\
       
       'logMstar_best','logSFR_best','logsSFR_best',\
       'fmu_SFR_best','fmu_IR_best','mu_best',\
       'tau_V_best','Ldust_tot_best',\
       'Tc_ISM_best', 'Tw_BC_best',\
       'xi_C_tot_best','xi_PAH_tot_best','xi_MIR_tot_best','xi_W_tot_best',\
       'tau_V_ISM_best','Mdust_best',\
       
       'chisq_best','magphysFlag']

print('length of columns and names = ',len(data_columns),len(names))
tab = Table(data=data_columns,names=names)





tab.write(output_table,format='fits',overwrite=True)

    
