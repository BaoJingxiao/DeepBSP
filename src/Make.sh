#!/bin/bash
#Settings
NetCDFPath="../packages/NetCDF"
Fortran="gfortran"
Flags="-march=native -O3 -s"

#Make
if [ ! -d ../bin ];then
    mkdir ../bin
fi
$Fortran Cal_ElementVoxel_PDB.f95 $Flags -I $NetCDFPath/include $NetCDFPath/lib/libnetcdff.a $NetCDFPath/lib/libnetcdf.a $NetCDFPath/lib/libhdf5_hl.a $NetCDFPath/lib/libhdf5.a $NetCDFPath/lib/libz.a -ldl -lrt -o ../bin/Cal_ElementVoxel_PDB
rm -f *.mod
