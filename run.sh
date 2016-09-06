#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
        echo "  **  Trapped CTRL-C"
}

./build/Release/EMC
echo "Finished Job"
