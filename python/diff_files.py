#!/usr/bin/env python

import glob
import os

damien_code = "/home/rfinn/software/magphys-damien/"

cb07 = "/home/rfinn/software/magphys/"

filelist = glob.glob(cb07+'*')

for f in filelist:
    fname = os.path.basename(f)
    if fname[-2:] == '.o':
        continue
    if fname[-2:] == '.a':
        continue
    if fname[-1] == '~':
        continue    
    print("FILE : ",fname)
    diff_string = f'diff {cb07}/{fname} {damien_code}/{fname}'
    os.system(diff_string)
    print()
    print()
