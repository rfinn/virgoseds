#!/usr/bin/env python


'''
GOAL: 
* run gatherMagphys.py to create a table for each magphys run (nozband, legacyExt, etc)

* combine outputs from each run into one table

* include 68% conf interval and chisq


USAGE:

python ~/github/virgoseds/python/mergeMagphysTables.py 

'''
import os
import sys
import numpy as np
from astropy.table import Table,hstack
HOME = os.getenv("HOME")
sys.path.append(HOME+'/github/Virgo/programs/')
from datetime import datetime



ids = ['-nozband','-nozband-legacyExt','-nozband-salimExt']


# get the current date so we know what the output files will be named
dateTimeObj = datetime.now()
myDate = dateTimeObj.strftime("%d-%b-%Y")

output_table_dir = HOME+'/research/Virgo/tables-north/v2/'

# create a list to contain output from 3 versions
output_tables = []
for i,id in enumerate(ids):
    magdir = f'{HOME}/research/Virgo/magphys/magphysParallelGrawp/output{id}/'

    # run gatherMagphys on output-nozband
    myid = id.replace('-nozband','nozband')
    programString = f'python ~/github/virgoseds/python/gatherMagphys.py --magdir {magdir} --id {myid}'
    #print(programString)
    print()    
    print('running gatherMagphys.py on {}'.format(id))
    print()
    os.system(programString)
    
    output_tables.append(Table.read(output_table_dir+'/vf_v2_magphys_'+id.replace('-nozband','nozband')+'_'+myDate+'.fits'))


# add id to column names to avoid duplicate names

for i,t in enumerate(output_tables):
    if i == 0:
        # this will keep all the columns as is for the first table in the list
        continue
    else:
        for c in t.colnames:
            
            if c == 'VFID':
                # remove VFID as this is a duplicate in each table
                t.remove_column(c)
            elif c == 'magphysFlag':
                # remove maphysFlag as this is a duplicate in each table
                t.remove_column(c)
            else:
                t.rename_column(c,c+ids[i].replace('-nozband',''))


# merge tables into one

combined_table = hstack(output_tables)



# write out the combined table

combined_name = output_table_dir+'/vf_v2_magphys_'+myDate+'.fits'

combined_table.write(combined_name,format='fits',overwrite=True)


