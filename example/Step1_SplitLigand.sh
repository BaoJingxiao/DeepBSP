#!/bin/bash
#Setting
pdb="1qf1"

#Split ligand Mol2 files from decoys.mol2
cd "CASF2016_Docking_"$pdb || exit
rm -rf $pdb"_ligand"
mkdir -p $pdb"_ligand"
cd $pdb"_ligand" || exit
LineID1=`cat "../Data0/"$pdb"_ligand.mol2" | grep -n "^@<TRIPOS>MOLECULE" | awk -F ":" '{print($1)}'`
LineID2=`cat "../Data0/"$pdb"_ligand.mol2" | grep -n "^@<TRIPOS>SUBSTRUCTURE" | awk -F ":" '{print($1)}'`
cat "../Data0/"$pdb"_ligand.mol2" | head -n $[$LineID2+1] | tail -n $[$LineID2+2-$LineID1] > Lig_Data0.mol2
echo " 0.00" > Lig_RMSD.dat
cd ..
TotLines=`cat "Data0/"$pdb"_decoys.mol2" | wc -l`
LineIDS=`cat "Data0/"$pdb"_decoys.mol2" | grep -n "^@<TRIPOS>MOLECULE" | awk -F ":" '{print($1)}'`
LineIDS=`echo $LineIDS" "$TotLines`
Flag=0
for StopLineID in $LineIDS;do
    if [ "$Flag" -eq 0 ];then
        Flag=1
        StartLineID=$StopLineID
    else
        rm -rf $pdb"_Temp"
        mkdir -p $pdb"_Temp"
        cd $pdb"_Temp" || exit
        cat "../Data0/"$pdb"_decoys.mol2" | head -n $[$StopLineID-1] | tail -n $[$StopLineID-$StartLineID] > Temp.mol2
        LineID1=`cat Temp.mol2 | grep -n "^@<TRIPOS>MOLECULE" | awk -F ":" '{print($1)}'`
        LineID2=`cat Temp.mol2 | grep -n "^@<TRIPOS>SUBSTRUCTURE" | awk -F ":" '{print($1)}'`
        cat Temp.mol2 | head -n $[$LineID2+1] | tail -n $[$LineID2+2-$LineID1] > Lig_Data0.mol2
        rm -f Temp.mol2
        LigPoseName=`cat Lig_Data0.mol2 | head -n 2 | tail -n 1`
        cat "../Data0/"$pdb"_rmsd.dat" | grep -w "^$LigPoseName " | awk -F " " '{printf("%5.2f\n",$2)}' > Lig_RMSD.dat
        cd ..
        mv $pdb"_Temp" $LigPoseName
        StartLineID=$StopLineID
    fi
done
cd ..

