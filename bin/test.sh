#!/bin/bash
echo compiling
make -j 8
echo moving to simulation directoy
cd ../simulations
echo testing Simbuca with benchmark 2
# ../simulations/simbuca ../simulation/rf.ini
./simbuca rf.ini
echo removing any output files
# rm out_rf_*
