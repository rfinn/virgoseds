#!/usr/bin/env python

'''
GOAL:
* convert data from John Moustakas's legacy photometry file into magphys input
* by default, we use the fluxes from AP06
* we will worry about the conversion to total mag after running magphys
* extinction correction can be none, Legacy, or Salim

NOTES:
* according to magphys manual, 

  observations.dat
  ID : galaxy name
  redshift : galaxy redshift
  flux in Jy: flux_0 flux_0_err flux_1 flux_1_err

* document with notes on how we decided on aperture and other issues

  https://docs.google.com/document/d/1vdBfzbZlWhm8_edVWIReyqIPKoT_DQPlCIUDs5KuuZA/edit?usp=sharing
  

* notebook with more details on the full magphys process

  https://github.com/rfinn/virgoseds/blob/main/notebooks/creating-magphys-input-files.ipynb
  

* repository with this code and more info on the wiki

  https://github.com/rfinn/virgoseds
  
'''

from astropy.table import Table
from astropy.io import fits
import os
import numpy as np
import speclite.filters
import argparse

###########################
##### SET UP ARGPARSE
###########################

parser = argparse.ArgumentParser(description ='Program to convert phot files from JM to magphys input')
#parser.add_argument('--sbmag', dest = 'sbmag', default = 24, help = 'sb to fit.  default is 24, as this has the highest fraction of non-zero entries.  The options are [22,26] in increments of 0.5.')
parser.add_argument('--ext', dest = 'ext', default = 0, help = 'extinction correction to apply.  0=None; 1=Legacy Survey; 2=Salim/Leroy.  The main difference between 1 and 2 is how the GALEX fluxes are handled.  See Leroy+2019 and Salim+2016 for more details.')
args = parser.parse_args()

#################################################
##### SET UP FILTERS AND EXTINCTION COEFFICIENTS
#################################################


filters = ['FUV','NUV','G','R','Z','W1','W2','W3','W4']
# filter IDS from magphys filters.log
filter_ids_S = ['123','124','426','427','429','280','281','282','414']
filter_ids_N = ['123','124','551','552','553','280','281','282','414']



writeFilterFiles = False
if writeFilterFiles:
    ##################################
    # FILTER EFFECTIVE WAVELENGTHS 
    ##################################
    grz_north = ['BASS-g','BASS-r','MzLS-z']
    grz_south = ['decamDR1-g','decamDR1-r','decamDR1-z']
    South=True
    if South:
        grz = grz_south
    else:
        grz = grz_north
    all_effective_wavelengths = []
    # filters need to be ordered in terms of increasing effective wavelength
    FUV = 1516e-4 # in um https://asd.gsfc.nasa.gov/archive/galex/Documents/ERO_data_description_2.htm
    NUV = 2267e-4
    all_effective_wavelengths.append(str(FUV))
    all_effective_wavelengths.append(str(NUV))

    #grz_filters = speclite.filters.load_filters(grz)
    for f in grz:
        specfilter = speclite.filters.load_filters(f)
        all_effective_wavelengths.append('{:.5f}'.format(specfilter[0].effective_wavelength.to('micron').value))
    
        
    wise = speclite.filters.load_filters('wise2010-*') 
    W1 = '{:.5f}'.format(wise[0].effective_wavelength.to('micron').value)
    W2 = '{:.5f}'.format(wise[1].effective_wavelength.to('micron').value)
    W3 = '{:.5f}'.format(wise[2].effective_wavelength.to('micron').value)
    W4 = '{:.5f}'.format(wise[3].effective_wavelength.to('micron').value)

    wise_filters = [W1,W2,W3,W4]
    all_effective_wavelengths += wise_filters

    # WRITE OUT FILTERS.DAT
    outfile = '/home/rfinn/research/Virgo/legacy-phot/legacyFiltersN.dat'
    outf = open(outfile,'w')
    outf.write('#name  lambda_eff   filter_id   fit?\n')
    for i,f in enumerate(filters):
        s = '{}   {}    {}    {}\n'.format(filters[i],all_effective_wavelengths[i],filter_ids_N[i],1)
        outf.write(s)
    outf.close()
    # WRITE OUT FILTERS.DAT
    outfile = '/home/rfinn/research/Virgo/legacy-phot/legacyFiltersS.dat'
    outf = open(outfile,'w')
    outf.write('#name  lambda_eff   filter_id   fit?\n')
    for i,f in enumerate(filters):
        s = '{}   {}    {}    {}\n'.format(filters[i],all_effective_wavelengths[i],filter_ids_S[i],1)
        outf.write(s)
    outf.close()


##########################################################################################
###  John suggested using surface brightness magnitudes
###  this is what we are doing
###  but we scale the sb fluxes by a factor
###  based on how the total magnitudes from COG compare with sb mag
##########################################################################################

#sb2fit = (args.sbmag) # could be 23, 24, 25, 26
fluxes = ['FLUX_AP06_{}'.format(f) for f in filters]
ivars = ['FLUX_IVAR_AP06_{}'.format(f) for f in filters]

# for determining scaling factor to convert aperture fluxes to total
# we are actually not using this - we can do the scaling after running magphys
# therefore, we are just using AP06 fluxes for running magphys
mtots = ['COG_MTOT_{}'.format(f) for f in filters]

# read in legacy phot
#photfile = '/home/rfinn/research/Virgo/legacy-phot/virgofilaments-legacyphot.fits'


# updating for v2, which is almost the entire catalog
photfile = '/home/rfinn/research/Virgo/legacy-phot/virgofilaments-v2-legacyphot.fits'
mef_table = fits.open(photfile)
# extensions are ['PARENT','ELLIPSE','TRACTOR']
ephot = Table(mef_table['ELLIPSE'].data) # first hdu is the elliptical photometry
ephot1 = Table(mef_table['PARENT'].data) # this one has RA and DEC, which we need to separate N and S
mef_table.close()



# convert fluxes from nmgy to Jy
flux_Jy = np.zeros((len(ephot),len(filters)),'d')
err_Jy = np.zeros((len(ephot),len(filters)),'d')


for i,f in enumerate(fluxes):
    flux_Jy[:,i] = ephot[f]*3.631e-6 # convert from nanomaggy to Jy
    #print(mtots[i])
    #print(ephot[mtots[i]])


for i,f in enumerate(ivars):
    flag =  ephot[f] > 0
    err_Jy[:,i][flag] = 1./np.sqrt(ephot[f][flag])*3.631e-6 # convert inverse variance to error, also convert from nanomaggy to Jy

    # set error for galaxies with ivar=0 to large value
    err_Jy[:,i][~flag] = 1e6*np.ones(sum(~flag),'d')




    
################################################
## FIND THE BEST APERTURE FLUX TO USE
##
## use the r-band data for this
## the default is set by the input parameter
## but need to make adjustments if this is not available
## also need to track which aperture flux we use
##
##
################################################

# create a column to track the input SB used
input_sb_colname = []
input_sb = []
ApFluxInput = np.zeros(len(ephot),'bool')
for i in range(len(ephot)):
    input_sb_colname.append('FLUX_AP06_R')
    input_sb.append(6)    


# find galaxies whose r-band flux is zero
zeroRfluxFlag =  flux_Jy[:,3] <= 0
zeroIndices = np.arange(len(ephot))[zeroRfluxFlag]
print('number of galaxies with AP06 == 0 is {}'.format(len(zeroIndices)))

# search through alternate fluxes to find a valid one
# ordering in terms of SB mags with highest completeness
possible_aps = [5,4,7,8,3,2,1]    
nbadR = 0
for i in zeroIndices:

    foundAltInputFlux = False
    for m in possible_aps:
        colname = 'FLUX_AP{:02d}_R'.format(m)

        if ephot[colname][i] > 0:
            #print(f'found alternate input flux for {i} us SB={m}')
            # exit mag loop once non-zero option is found
            foundAltInputFlux = True


            # track the SB of the valid one if an alternative is found
            input_sb_colname[i] = colname
            input_sb[i] = m            

            # set all fluxes to this SB for this galaxy
            for j,f in enumerate(filters):

                colname = 'FLUX_AP{:02d}_{}'.format(m,f)
                ivarname = 'FLUX_IVAR_AP{:02d}_{}'.format(m,f)        
                flux_Jy[i,j] = ephot[colname][i]*3.631e-6
                err_Jy[i,j] = 1./np.sqrt(ephot[ivarname][i])*3.631e-6 # convert inverse variance to error, also convert from nanomaggy to Jy

            break

    if not foundAltInputFlux:
        print('WARNING: no alternate r-band input flux:{}, VFID{:04d} ({})'.format(i,ephot['VF_ID'][i],ephot['GALAXY'][i]))
        # could check SB fluxes here to see if any are valid
        possible_mags = [26,25.5,25,24.5,24,23.5,23,22.5,22]
        foundAltInputSB = False
        for m in possible_mags:
            colname = f'FLUX_SB{m}_R'
            if ephot[colname][i] >= 0:
                input_sb_colname[i] = colname
                input_sb[i] = m                            
                # set all fluxes to this SB for this galaxy
                for j,f in enumerate(filters):

                    colname = f'FLUX_SB{m}_{f}'
                    ivarname = f'FLUX_IVAR_SB{m}_{f}'
                    flux_Jy[i,j] = ephot[colname][i]*3.631e-6
                    err_Jy[i,j] = 1./np.sqrt(ephot[ivarname][i])*3.631e-6 # convert inverse variance to error, also convert from nanomaggy to Jy
                print('\tGood news: found SB flux :{}, VFID{:04d} ({})'.format(i,ephot['VF_ID'][i],colname))
                foundAltInputSB=True
                break

        if not foundAltInputSB:
            nbadR += 1


        
print()
print('total number with no valid R-band flux as input = {}'.format(nbadR))
print()
print()
print()
print()



scaleFluxes = False

if scaleFluxes:

    ################################################
    ## SCALE FLUXES BY RATIO OF FLUX_MTOT/FLUX_R24
    ## use the r-band data for this
    ##
    ## THIS PART IS OBSOLETE B/C WE ARE USING AP06 FLUXES
    ## AND ALL GALAXIES HAVE A MEASURE OF THIS
    ##
    ## WE ARE NOT SCALING FLUXES TO TOTAL IN THIS PROGRAM!!!
    ##
    ################################################
    rmag_tot = ephot['COG_MTOT_R']
    # check to see if rmag_tot = -1; this happens if the John's cog fails

    total_flux_column = []
    for i in range(len(ephot)):
        total_flux_column.append('COG_MTOT_R')

    possible_mags = [26,25.5,25,24.5,24,23.5,23,22.5,22]
    bad_scale_factor = np.zeros(len(ephot),'bool')
    nbadTotal=0
    badCOG = 0
    for i,r in enumerate(rmag_tot):

        if r == -1:
            foundAltTotal = False
            # just search for fluxes in fainter SB
            print('{}, COG_MTOT_R = {}'.format(i,r))
            badCOG += 1
            if not foundAltTotal:
                foundAltApTotal=False
                print('WARNING: no COG r-band magnitude to use as total for:{}, VFID{:04d} ({}) - checking AP FLUXES'.format(i,ephot['VF_ID'][i],ephot['GALAXY'][i]))
                for nap in np.arange(8,input_sb[i],-1):
                    colname = 'FLUX_AP{:02d}_R'.format(nap)
                    if ephot[colname][i] >= 0:
                        scaling_col = colname
                        total_flux_column[i] = colname
                        foundAltApTotal = True
                        #print('{}: COG_MTOT_R = -1 for {} but found {} {}'.format(i,ephot['GALAXY'][i],scaling_col, ephot[colname][i]))
                        rmag_tot[i] = 22.5 - 2.5*np.log10(ephot[scaling_col][i])
                        break
                if not foundAltApTotal:
                    print("\t could not find an ap flux to use as a total either :(")
                    nbadTotal += 1


    print()
    print('total number with no valid total R-band flux = {}'.format(nbadTotal))
    print()
    print()
    print()
    print()

        
    flux_tot = 10.**((22.5-rmag_tot)/2.5) # flux in nanomaggies
    flux_tot_Jy = flux_tot*3.631e-6
    scale_factor = flux_tot_Jy/flux_Jy[:,3]

    # galaxies with no valid total flux
    # set their scale factor to 1
    noTotalFlux = rmag_tot == -1
    if np.sum(noTotalFlux) > 0:
        scale_factor[noTotalFlux] = np.ones(np.sum(noTotalFlux))

    # scale the flux and error in all wavebands
    scaled_flux_Jy = np.zeros(flux_Jy.shape,'d')
    scaled_err_Jy = np.zeros(flux_Jy.shape,'d')
    for i in range(len(scale_factor)):
        # check for nans and inf
        if (not np.isfinite(scale_factor[i])) | (scale_factor[i] < 1):
            scale_factor[i] = 1
        scaled_flux_Jy[i] = flux_Jy[i]*scale_factor[i]
        scaled_err_Jy[i] = err_Jy[i]*scale_factor[i]

else:
    scaled_flux_Jy = flux_Jy
    scaled_err_Jy = err_Jy
        

###########################################################
## read in vf_main so we can get redshift for each galaxy
##
## updating to v2 catalogs in May 2022
###########################################################
#vffile = '/home/rfinn/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
vffile = '/home/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits'
vfmain = Table.read(vffile)

# use this file to get the E(B-V) values
vffile = '/home/rfinn/research/Virgo/tables-north/v2/vf_v2_extinction.fits'
vfext = Table.read(vffile) 

# store redshift of each galaxy
redshift = np.zeros(len(ephot),'d')
redshift_flag = np.ones(len(ephot),'bool')
for i in range(len(ephot)):
    #print(ephot['VF_ID'][i])
    vfindex = np.arange(len(vfmain))[vfmain['VFID'] == 'VFID{:04d}'.format(ephot['VF_ID'][i])]
    #print(vfindex)
    #print(ephot['VF_ID'][i],vfindex)
    #if len(vfindex) == 0:
    redshift[i] = vfmain['vr'][vfindex]/3.e5 # divide by speed of light to convert recession velocity to redshift





###########################################################
# create output table
# ID redshift flux_0 flux_0_err flux_1 flux_1_err ...
#
# 2022-Jun-03 add extinction correction using
# legacy survey prescription
#
# NEED TO UPDATE TO ALLOW FOR SALIM TREATMENT OF EXTINCTION
# this will allow us to compare with Leroy+2019 values
###########################################################
output_columns = [ephot1['VFID'],redshift]
for i in range(len(filters)):

    # check to see if extinction correction should be applied
    if float(args.ext) > 0:
        if int(args.ext) == 1:
            ecolname = f'A({filters[i]})_SFD'
        elif int(args.ext) == 2:
            # can change suffice to _S16 to use extinction law use in Salim+2016
            ecolname = f'A({filters[i]})_S16'
        else:
            print('Extinction selection is not valid.  Needs to be [0,1,2] for [none,legacy,salim]')
            sys.exit()
        # convert magnitudes of extinction into a linear scale factor
        extcorr = 10.**(vfext[ecolname]/2.5)
        output_columns.append(scaled_flux_Jy[:,i]*extcorr[vfindex])
        output_columns.append(scaled_err_Jy[:,i]*extcorr[vfindex])
    else:
        output_columns.append(scaled_flux_Jy[:,i])
        output_columns.append(scaled_err_Jy[:,i])

#print(len(output_columns))
colnames = ['#VFID','redshift']
for i in range(len(filters)):
    colnames.append(filters[i])
    colnames.append('{}_err'.format(filters[i]))

out_table = Table(output_columns)#,names=colnames)

# define north flag to separate grz filters from DES vs BASS+MzLS
north_flag = ephot1['DEC'] >= 32.375

# write out north phot table
if int(args.ext) == 0:
    outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputN.dat'
elif int(args.ext) == 1:
    outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputN_legacyExt.dat'
elif int(args.ext) == 2:
    outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputN_salimExt.dat'

out_table[north_flag].write(outfile,format='ascii',overwrite=True,names=colnames)

# write out south phot table
outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputS.dat'
if int(args.ext) == 1:
    outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputS_legacyExt.dat'
elif int(args.ext) == 2:
    outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInputS_salimExt.dat'
out_table[~north_flag].write(outfile,format='ascii',overwrite=True,names=colnames)


nsf_flag = ephot1['VFID'] == 'VFID1778'
outfile = '/home/rfinn/research/Virgo/legacy-phot/NSF2021VFID1778.dat'
out_table[nsf_flag].write(outfile,format='ascii',overwrite=True,names=colnames)
# writing out the first galaxy, bc this is the only one with all non-zero fluxes at sb=23
#outfile = '/home/rfinn/research/Virgo/legacy-phot/magphysInput-gal1.dat'
#ptab = Table(out_table[0])
#ptab.write(outfile,format='ascii.fast_no_header',overwrite=True,comment=False)
# 


# save table with input_sb, total_flux_column, scale_factor
