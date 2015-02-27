1) Simbuca usage
---------------
all the information is stored and up to date on the Simbuca wiki page
https://sourceforge.net/p/simbuca/wiki/Home/ 


2) installing Simbuca
------------------=---
cd bin
make all


3)Working with the GPU or CPU
-----------------------------
nano Makefile (or use another editor)
change the PU variable to 'cpu' or 'gpu'
type make clear to remove all object files
type make all to make all object files again but now fore the other processing unit


4)Compiling with icpc or g++
----------------------------
pico Makefile (or use another editor)
change the COMPILER variable to 'icpc' or 'g++' 


5)Doing a trap simulation
-------------------------
By default the simbuca executable is stored in the simulations directory with some example inputfiles.

You can execute Simbuca by feeding it with an input file. For example
./simbuca example1.sim


It is advised to adapt the example.sim file and use this for your simulations. 
More information about how to set up the example.sim file can be found in the FAQ, see 1)



6)SVN help
----------
good site: http://wiki.apertium.org/wiki/Using_SVN


6)More information
------------------
Check the project website: https://sourceforge.net/projects/simbuca/
