#!/bin/bash

n=200000
k=10
reg=1
r=.1
verbose=0

export OMP_NUM_THREADS=1

echo ./test.exe_arm $n $k $reg $r $verbose
./test.exe_arm $n $k $reg $r $verbose