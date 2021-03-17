#!/bin/bash
#Setting
Version_zlib="1.2.11"
Version_hdf5="1.10.7"
Version_netcdf_c="4.7.4"
Version_netcdf_fortran="4.5.3"
Prefix=`pwd`

#zlib
tar -xzf "zlib-"$Version_zlib".tar.gz"
cd "zlib-"$Version_zlib
./configure --prefix=$Prefix
make && make install
cd ..
rm -rf "zlib-"$Version_zlib

#hdf5
tar -xjf "hdf5-"$Version_hdf5".tar.bz2"
cd "hdf5-"$Version_hdf5
./configure --prefix=$Prefix --with-zlib=$Prefix --enable-hl
make && make install
cd ..
rm -rf "hdf5-"$Version_hdf5

#netcdf-c
tar -xzf "netcdf-c-"$Version_netcdf_c".tar.gz"
cd "netcdf-c-"$Version_netcdf_c
./configure CPPFLAGS="-I$Prefix/include" LDFLAGS="-L$Prefix/lib" --prefix=$Prefix --disable-dap
make && make install
cd ..
rm -rf "netcdf-c-"$Version_netcdf_c

#netcdf-fortran
tar -xzf "netcdf-fortran-"$Version_netcdf_fortran".tar.gz"
cd "netcdf-fortran-"$Version_netcdf_fortran
export LD_LIBRARY_PATH=$Prefix"/lib":$LD_LIBRARY_PATH
./configure CPPFLAGS="-I$Prefix/include" LDFLAGS="-L$Prefix/lib" --prefix=$Prefix --enable-shared
make && make install
cd ..
rm -rf "netcdf-fortran-"$Version_netcdf_fortran
