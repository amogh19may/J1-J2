#!/bin/bash
echo "hello"
read code
cd scratch
mkdir $code
cd ..
cp $code scratch/$code
cd scratch/$code
mpifort $code -o op.exe
mpirun -n 4 ./op.exe > op.txt
cd ..
mv  $code ../am
echo "run over"
