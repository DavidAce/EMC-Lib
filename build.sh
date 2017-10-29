#!/bin/bash

#rm -rf *
buildtype="Release"

if [[ "$@" == *"ebug"* ]]
then
	buildtype="Debug"
fi

if [[ "$@" == *"lean"* ]]
then
    echo "Cleaning build"
	rm -rf build
    exit 0
fi


mkdir build
cd build
mkdir ${buildtype}
cd ${buildtype}
echo "Starting Build"
cmake -DCMAKE_BUILD_TYPE=${buildtype}  ../../
make