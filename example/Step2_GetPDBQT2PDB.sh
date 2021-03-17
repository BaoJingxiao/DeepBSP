#!/bin/bash
#Setting
pdb="1qf1"
AutoDockToolsPath="/home/paul/AutoDockVina/MGLToolsPckgs/AutoDockTools/Utilities24"
AutoDockToolsPythonSH="/home/paul/AutoDockVina/bin/pythonsh"

#Prepare receptor PDB file with prepare_receptor4.py and change PDBQT format to PDB format with pdbqt_to_pdb.py
cd "CASF2016_Docking_"$pdb"/Data0" || exit
$AutoDockToolsPythonSH $AutoDockToolsPath/prepare_receptor4.py -r $pdb"_protein.pdb" -o Rec.pdbqt -A none -U nphs_lps_waters -v > Log_GetPDBQT2PDB.log 2>&1
$AutoDockToolsPythonSH $AutoDockToolsPath/pdbqt_to_pdb.py -f Rec.pdbqt -o Rec_noNH.pdb >> Log_GetPDBQT2PDB.log 2>&1
rm -f Rec.pdbqt
if [ ! -s "Rec_noNH.pdb" ];then
    echo "Failed get Rec_noNH.pdb for "$pdb
    exit
fi
cd ..

#Prepare ligand Mol2 file with prepare_ligand4.py and change PDBQT format to PDB format with pdbqt_to_pdb.py
LIGS=`ls -1 -v | grep "^"$pdb"_"`
for Lig in $LIGS;do
    cd $Lig || exit
    ln -s -f ../Data0/Rec_noNH.pdb .
    $AutoDockToolsPythonSH $AutoDockToolsPath/prepare_ligand4.py -l Lig_Data0.mol2 -o Lig.pdbqt -v > Log_GetPDBQT2PDB.log 2>&1
    $AutoDockToolsPythonSH $AutoDockToolsPath/pdbqt_to_pdb.py -f Lig.pdbqt -o Lig_noNH.pdb >> Log_GetPDBQT2PDB.log 2>&1
    rm -rf Lig.pdbqt
    if [ ! -s "Lig_noNH.pdb" ];then
        echo "Failed get Lig_noNH.pdb for "$Lig
    fi
    cd ..
done
cd ..

