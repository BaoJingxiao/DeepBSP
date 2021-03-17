#!/bin/bash
#Setting
pdb="1qf1"

#Get complex names and descriptor files list
rm -rf ComplexList.dat
PPWD1=`pwd`
cd "CASF2016_Docking_"$pdb || exit
LIGS=`ls -1 -v | grep "^"$pdb"_"`
for Lig in $LIGS;do
    cd $Lig
    PPWD2=`pwd`
    echo $Lig" "$PPWD2"/Descriptor.dat.nc" >> $PPWD1"/ComplexList.dat"
    cd ..
done
cd ..

