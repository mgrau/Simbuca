#!/bin/bash
echo compiling
make -j 8
echo testing Simbuca with benchmark 2
../simulations/simbuca ../simulations/bench2.sim
echo removing any output files
rm out_bench2_logfile.txt
rm out_bench2_pAll.txt
