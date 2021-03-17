#!/bin/bash
#Setting
ComplexList="./ComplexList.dat"
ResultFile="./PoseRMSD.dat"
DeepBSP="../bin/DeepBSP_GPU.py"
NetDataPath="../dat"
BenchSize=32
BenchNumWorkers=2
TotalThreads=2

#Run calculation
rm -rf $ResultFile
python $DeepBSP $ComplexList $NetDataPath $ResultFile $BenchSize $BenchNumWorkers $TotalThreads
