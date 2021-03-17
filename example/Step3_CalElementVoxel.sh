#!/bin/bash
#Setting
pdb="1qf1"
ElementVoxelProgram="../../../bin/Cal_ElementVoxel_PDB"
ElementVoxelDataFile="../../../dat/Cal_ElementVoxel_PDB_ElementVDWRadii.dat"

#Calculate ElementVoxel descriptor for each complexes
cd "CASF2016_Docking_"$pdb || exit
LIGS=`ls -1 -v | grep "^"$pdb"_"`
for Lig in $LIGS;do
    cd $Lig || exit
    PPWD=`pwd`
    #Create the config file for descriptor calculation
    echo "#configures for calculate the descriptor of complex"  > Descriptor.config
    echo "#file path (should end with /):"                     >> Descriptor.config
    echo $PPWD"/"                                              >> Descriptor.config
    echo "#receptor pdb file name:"                            >> Descriptor.config
    echo "Rec_noNH.pdb"                                        >> Descriptor.config
    echo "#ligand pdb file name:"                              >> Descriptor.config
    echo "Lig_noNH.pdb"                                        >> Descriptor.config
    echo "#output data file name:"                             >> Descriptor.config
    echo "Descriptor.dat.nc"                                   >> Descriptor.config
    echo "#element VDW radii file:"                            >> Descriptor.config
    echo $ElementVoxelDataFile                                 >> Descriptor.config
    echo "#grid size of X:"                                    >> Descriptor.config
    echo "24.0"                                                >> Descriptor.config
    echo "#grid size of Y:"                                    >> Descriptor.config
    echo "24.0"                                                >> Descriptor.config
    echo "#grid size of Z:"                                    >> Descriptor.config
    echo "24.0"                                                >> Descriptor.config
    echo "#grid offset:"                                       >> Descriptor.config
    echo "1.0"                                                 >> Descriptor.config
    echo "#VDW distance cutoff:"                               >> Descriptor.config
    echo "10.0"                                                >> Descriptor.config
    echo "#end"                                                >> Descriptor.config
    #Run descriptor calculation
    $ElementVoxelProgram Descriptor.config
    rm -f Descriptor.config
    cd ..
done
cd ..

