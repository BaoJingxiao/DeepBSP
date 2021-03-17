#!/bin/bash
#Setting
pdb="1qf1"

#Clean files
cd "CASF2016_Docking_"$pdb || exit
ls -1 -v | grep "^"$pdb"_" | xargs rm -rf
cd Data0 || exit
rm -rf Log_GetPDBQT2PDB.log Rec_noNH.pdb
cd ../../
