#!/bin/bash

# Search for sqlite tgz file
SQLITE_FOUND=`find ext_libs/ -type f -name sqlite.tar.gz`
TOP_DIR=`pwd`

# If the tgz file isn't there, download it
if [ -z $SQLITE_FOUND ]
then
    echo "sqlite3 source code not found"
    # Get the sqlite source code
    curl -o ext_libs/sqlite.tar.gz https://www.sqlite.org/2021/sqlite-autoconf-3350500.tar.gz
fi

# Create Target Directory
cd ext_libs/
mkdir sqlite3

# Extract Contents
tar -xf sqlite.tar.gz -C sqlite3 --strip-components=1

# Build sqlite3
cd sqlite3
mkdir $TOP_DIR/build/ext_libs/sqlite3
./configure --prefix=$TOP_DIR/build/ext_libs/sqlite3
make
make install
