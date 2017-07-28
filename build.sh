#!/bin/bash

./configure

rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

mkdir -p lib
cd lib

echo "Download ZLIB"
if [ -f "zlib-1.2.11.tar.gz" ]
then
	echo "zlib-1.2.11.tar.gz found."
    rm -rf zlib-1.2.11
else
	wget http://zlib.net/zlib-1.2.11.tar.gz
fi
tar -xf zlib-1.2.11.tar.gz

echo "Build ZLIB"
mv zlib-1.2.11 zlib
cd zlib
zlibdir=`pwd`
./configure --prefix=$zlibdir/build --static --64
make all
make all install
cd ../

echo "Download JEMALLOC"
if [ -f "jemalloc-4.5.0.tar.bz2" ]
then
	echo "jemalloc-4.5.0.tar.bz2 found."
    rm -rf jemalloc-4.5.0
else
	wget https://github.com/jemalloc/jemalloc/releases/download/4.5.0/jemalloc-4.5.0.tar.bz2
fi
tar -xf jemalloc-4.5.0.tar.bz2

echo "Build JEMALLOC"
mv jemalloc-4.5.0 jemalloc
cd jemalloc
autoconf
buildir=`pwd`
./configure --prefix=$buildir/build --with-jemalloc-prefix="je_"
make all
make all install
cd ../

echo "Download HTSLIB"
if [ -f "htslib-1.4.tar.bz2" ]
then
	echo "htslib-1.4.tar.bz2 found."
    rm -rf htslib-1.4
else
	wget https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2
fi
tar -xf htslib-1.4.tar.bz2

echo "Build HTSLIB"
mv htslib-1.4 htslib
cd htslib
autoheader
autoconf
buildir=`pwd`
./configure --prefix=$buildir/build
make all
make all install
cd ../

echo "Download HDF5"
if [ -f "hdf5-1.8.18.tar.bz2" ]
then
	echo "hdf5-1.8.18.tar.bz2 found."
    rm -rf hdf5-1.8.18
else
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.18/src/hdf5-1.8.18.tar.bz2
fi
tar -xf hdf5-1.8.18.tar.bz2

echo "Build HDF5"
mv hdf5-1.8.18 hdf5
cd hdf5
mkdir -p objs
buildir=`pwd`
cd objs
cmake -DCMAKE_INSTALL_PREFIX:PATH=$buildir -DBUILD_SHARED_LIBS:BOOL=OFF -DZLIB_DIR:PATH=$zlibdir/build/include -DZLIB_LIBRARY:FILEPATH=$zlibdir/build/lib -DCMAKE_BUILD_TYPE:STRING=ON ../
make all
make all install
cd ../
rm -rf objs
cd ../

echo "Download libdivsufsort"
if [ -d "libdivsufsort" ]
then
	echo "libdivsort found."
else
	git clone https://github.com/lh3/libdivsufsort
fi

echo "Build libdivsufsort"
cd libdivsufsort
mkdir -p build
buildir=`pwd`
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$buildir -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON ../
make all
make all install
cd ../
rm -rf build

echo "Build lib finish"

cd ../../
buildir=`pwd`
echo "Build HERA at $buildir"

UNAME=`uname`

Makefile_mac="Makefile_mac"
Makefile_linux="Makefile_linux"

if [ $UNAME == "Darwin" ]
then
	echo "Mac"
	make -f $Makefile_mac
else
	echo "Linux"
	make -f $Makefile_linux
fi
