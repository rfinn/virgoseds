#! /usr/bin/env bash

#setenv FILTERS            $magphys/FILTERBIN.RES
#setenv OPTILIB		  $magphys/OptiLIB_cb07.bin
#setenv OPTILIBIS	  $magphys/OptiLIBis_cb07.bin
#setenv IRLIB	          $magphys/InfraredLIB.bin

# user input files [this can be changed by user]
#
#setenv USER_FILTERS	$magphys/eg_user_files/filters.dat
#setenv USER_OBS	$magphys/eg_user_files/observations.dat
#setenv USER_FILTERS	filters.dat
#setenv USER_OBS		observations.dat

export magphys=/home/siena.edu/rfinn/software/magphys/
# already defining environment variables above
echo "HEY!!!"
echo $1
if [ -f $1/$1.fit ];
then
    # exist this script
    echo 'output exists for '$1
    echo 'moving to the next galaxy'
    rm *.lbr
    exit 0
fi
cd $1
source $magphys/.magphys_bashrc
printf '\n starting make_zgrid\n\n'
$magphys/make_zgrid
printf '\n starting get_libs_bash\n'
source $magphys/get_libs_bash
printf '\n starting fit_sed\n'
$magphys/fit_sed
rm *.lbr

#cd ..

