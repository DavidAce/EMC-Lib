#!/bin/bash
echo "Running"
# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

if [[ "$@" == *"valgrind"* ]]
then
    valgrind --tool=memcheck --leak-check=full -v mpirun -n 4 ./build/Debug/EMC
elif [[ "$@" == *"gdb"* ]]
then
    mpirun -n 4 xterm -e gdb ./build/Debug/WL
#    mpirun -n 4 gdb build/Debug/WL

elif [[ "$@" == *"ebug"* ]]
then
    ulimit -c unlimited
    mpirun -n 4  -bind-to core:overload-allowed ./build/Debug/EMC
elif [[ "$@" == *"elease"* ]]
then
    ulimit -c unlimited
    mpirun -n 4  -bind-to core:overload-allowed ./build/Release/EMC
else
    ulimit -c unlimited
    mpirun -n 4  -bind-to core:overload-allowed ./build/Release/EMC
fi
function ctrl_c() {
        echo "** Trapped CTRL-C"
}



