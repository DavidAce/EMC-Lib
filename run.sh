#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
        echo "  **  Trapped CTRL-C"
}

echo "Starting job"

if [[ "$@" == *"valgrind"* ]]
then
    valgrind --tool=memcheck --leak-check=full -v ./build/Debug/EMC
elif [[ "$@" == *"gdb"* ]]
then
    gdb ./build/Debug/EMC
elif [[ "$@" == *"ebug"* ]]
then
    ulimit -c unlimited
   ./build/Debug/EMC
elif [[ "$@" == *"elease"* ]]
then
    build/Release/EMC
else
    build/Release/EMC
fi

echo "Finished Job"