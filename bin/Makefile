############################################ 
# 					   #
#        				   #
#     	Makefile for Simbuca v3.2 	   #	
#					   #
#					   #
############################################ 
VERSION=3.4

#---------------------------------------------#
#------------USER INPUT PARAMETERS------------#
#---------------------------------------------#

#compiler:  mpic++ / g++ / icpc
# COMPILER=g++
COMPILER=icpc
# COMPILER=mpic++

#processing unit: cpu / gpu
PU=gpu

#NBODY algorithm for the calculation: CUNBODY / NBODY
#after changing $PU, type 'make clean' to remove all *.o files.
NBODY_ALGO=CUNBODY
# NBODY_ALGO=NBODY

#simulation directory
SIMULATION_DIR=../simulations

#name of the executable
EXECUTABLE=$(SIMULATION_DIR)/simbuca

#---------------------------------------------#
#--------- SET NBODY ALGO --------------------#
#---------------------------------------------#
ifeq ($(PU),gpu)
  ifeq ($(NBODY_ALGO),NBODY)
  NBODY_ALGO_FLAGS = -D__NBODY_ON__ -I$(CUDASYS)/include 
  else
  NBODY_ALGO_FLAGS = -D__CUNBODY_ON__ -I$(CUDASYS)/include
  endif
endif
ifeq ($(PU),cpu)
  NBODY_ALGO_FLAGS = -D__CPUNBODY_ON__
  NBODY_ALGO = CPU
endif
#---------------------------------------------#


#---------------------------------------------#
#---------COMPILE AND LINK OPTIONS------------#
#---------------------------------------------#
#g++ Compiler
CC_gpp=g++ -O3 -g -Wno-write-strings $(NBODY_ALGO_FLAGS) -c
LINK_gpp=g++ -lm -O3 -o

#icpc Compiler
CC_icpc=icpc -fast $(NBODY_ALGO_FLAGS) -c
LINK_icpc=icpc -o  

# CUDA Compiler
CUDASYS = /usr/local/cuda/7.0
NVCC = $(CUDASYS)/bin/nvcc	
NVCCFLAGS = -ccbin gcc --compiler-options -fno-strict-aliasing,-Wall,-Wextra,-Wno-unused-parameter -O3
#MPI Compiler
MPI_PATH = /usr/local/mpich2/bin
CC_mpicpp=$(MPI_PATH)/mpicxx -fast -D__MPI_ON__ -finline-functions -Wno-write-strings $(NBODY_ALGO_FLAGS) -c 
LINK_mpicpp=$(MPI_PATH)/mpicxx -o 
#libraries
cunbodylib=../libraries/cunbody/libcunbody1.a_TI20140519_CUDA60
#--------------------------------------------#



#---------------------------------------------#
#---------------SET COMPILER------------------#
#---------------------------------------------#
ifneq ($(COMPILER),g++)
  ifneq ($(COMPILER),icpc)
   ifneq ($(COMPILER),mpic++)
    CC=UNSUPPORTED_COMPILER
    LINK=UNSUPPORTED_COMPILER
    $(error WRONG INPUT: use g++ or icpc as compiler. See 'make help'.)
   endif
  endif
endif

ifeq ($(COMPILER),g++)
  CC = $(CC_gpp) 
  LINK = $(LINK_gpp)
endif

ifeq ($(COMPILER),icpc)
  CC= $(CC_icpc) 
  LINK = $(LINK_icpc)
endif

ifeq ($(COMPILER),mpic++)
  CC = $(CC_mpicpp) 
  LINK = $(LINK_mpicpp)
endif


#------------CHECK PROCESSING UNIT------------#
ifneq ($(PU),cpu)
  ifneq ($(PU),gpu)
    CC=UNSUPPORTED_PROCESING_UNIT
    LINK=UNSUPPORTED_PROCESING_UNIT
    $(error WRONG INPUT: use gpu or cpu as processing unit. See 'make help'.)
  endif
endif



#------------BUILDING DEPENDENCIES------------#
all:	$(EXECUTABLE)

help:
	@echo "############################################ "
	@echo "Makefile v3.0 for Simbuca"
	@echo "Default compiler g++ . By default the program runs on CPU instead of the GPU"
	@echo "select cpu/gpu and icpc/g++/mpic++ in Makefile"
	@echo "Type 'make all' to build the executable $(EXECUTABLE)"
	@echo "############################################ "

clean: 
	rm -f *.o $(EXECUTABLE)

delink:
	rm $(EXECUTABLE)

parameters:
	@echo COMPILER=$(COMPILER)
	@echo PU=$(PU)
	@echo EXECUTABLE=$(EXECUTABLE)
	@echo NBODY_ALGO=$(NBODY_ALGO)
export:
	mkdir -p ../backup
	tar --exclude=.svn -c -v -p -s -f  ../backup/Simbuca_v$(VERSION).`date +"%Y%m%d-%H%M%S"`.tar ../bin/Makefile ../bin/*.cpp ../bin/*.h ../bin/*.a \
		../doc/readme.txt \
		../CoulombSnippet/* \
		../libraries/* \
		../tools/ioncloud_parser/main.cpp \
		../tools/ioncloud_parser/Makefile \
		$(SIMULATION_DIR)/example*.sim \
		$(SIMULATION_DIR)/trap_config.dat \
		../tools/fft/*/fft.cpp \
		../tools/fft/*/Makefile \
		../tools/column_selector/Makefile \
		../tools/column_selector/columnselect.cpp \
		../tools/movie/create_movie.txt \
		../tools/movie/draw_movie.C

ifeq ($(PU)$(NBODY_ALGO),gpuNBODY)
$(EXECUTABLE) : main.o parser.o ion_fly.o ode.o force.o matrix.o ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o operation.o nbody_kernel.o trapparameters.h initialization.o mpi_funcs.o nbody.o trapparameters.o numrec.o SLogger.o SLogWriter.o inireader.o ini.o
else
$(EXECUTABLE) : main.o parser.o ion_fly.o ode.o force.o matrix.o ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o operation.o trapparameters.h initialization.o mpi_funcs.o cunbody.o trapparameters.o numrec.o SLogger.o SLogWriter.o inireader.o ini.o
endif
	mkdir -p $(SIMULATION_DIR)
ifeq ($(PU),cpu)
	$(LINK) $(EXECUTABLE) main.o parser.o ion_fly.o ode.o force.o matrix.o ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o operation.o initialization.o mpi_funcs.o trapparameters.o numrec.o SLogger.o SLogWriter.o inireader.o ini.o
endif
ifeq ($(PU),gpu)
  ifeq ($(NBODY_ALGO),NBODY)
    $(LINK) $(EXECUTABLE)  main.o parser.o ion_fly.o ode.o force.o -L $(CUDASYS)/lib64 -lcuda -lcudart -lGL -lGLU -L $(CUDASYS)/lib -lcuda -lcudart -lGL -lGLU matrix.o  ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o nbody_kernel.o operation.o initialization.o mpi_funcs.o nbody.o trapparameters.o numrec.o SLogger.o SLogWriter.o inireader.o ini.o
  else
	$(LINK) $(EXECUTABLE)  main.o parser.o ion_fly.o ode.o force.o matrix.o  ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o operation.o initialization.o mpi_funcs.o cunbody.o trapparameters.o numrec.o  SLogger.o SLogWriter.o inireader.o ini.o $(cunbodylib) -L$(CUDASYS)/lib64 -lcudart -lcuda
  endif
endif	
@echo "compiled: $(EXECUTABLE)" 

#------------OBJECT FILES------------#
ifeq ($(NBODY_ALGO),NBODY)
%.o: %.cu 
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif

main.o : main.cpp parser.h ion_fly.h ode.h force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h MPI_simbuca.h IonTable.h PDGTable.h operation.h  initialization.h SMsgType.h SLogger.h SLogWriter.h
	$(CC) main.cpp
parser.o : parser.cpp parser.h ion_fly.h ode.h force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h SMsgType.h SLogger.h SLogWriter.h
	$(CC) parser.cpp
ion_fly.o : ion_fly.cpp ion_fly.h ode.h force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h SMsgType.h SLogger.h SLogWriter.h
	$(CC) ion_fly.cpp
ode.o : ode.cpp ode.h force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h MPI_simbuca.h mpi_funcs.h SMsgType.h SLogger.h SLogWriter.h
	$(CC) ode.cpp
force.o : force.cpp force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h  MPI_simbuca.h mpi_funcs.h SMsgType.h SLogger.h SLogWriter.h
ifeq ($(PU),cpu)
	$(CC) force.cpp
endif
ifeq ($(PU),gpu)
	$(CC) -I/usr/local/cuda/include  force.cpp
endif	
matrix.o : matrix.cpp matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h
	$(CC) matrix.cpp
ioncloud.o : ioncloud.cpp ioncloud.h particle.h ion.h mtrand.h logfile.h MPI_simbuca.h 
	$(CC) ioncloud.cpp
particle.o : particle.cpp particle.h ion.h mtrand.h logfile.h
	$(CC) particle.cpp
ion.o : ion.cpp ion.h mtrand.h logfile.h
	$(CC) ion.cpp
IonTable.o: IonTable.cpp IonTable.h
	$(CC) IonTable.cpp	
PDGTable.o: PDGTable.cpp PDGTable.h
	$(CC) PDGTable.cpp	
mtrand.o : mtrand.cpp mtrand.h logfile.h
	$(CC) mtrand.cpp
logfile.o: logfile.cpp logfile.h
	$(CC) logfile.cpp
operation.o: operation.cpp operation.h ion_fly.cpp ion_fly.h mpi_funcs.h
	$(CC) operation.cpp
mpi_funcs.o: mpi_funcs.cpp mpi_funcs.h
	$(CC) mpi_funcs.cpp
nbody.o: nbody.cpp nbody.h
	$(CC) nbody.cpp      
cunbody.o: cunbody.cpp cunbody.h
	$(CC) cunbody.cpp
initialization.o: initialization.cpp parser.h ion_fly.h ode.h force.h matrix.h ioncloud.h particle.h ion.h mtrand.h logfile.h MPI_simbuca.h IonTable.h PDGTable.h operation.h initialization.h trapparameters.h numrec.h inireader.h ini.h
	$(CC) initialization.cpp
trapparameters.o: trapparameters.cpp trapparameters.h
	$(CC) trapparameters.cpp   
numrec.o: numrec.cpp numrec.h
	$(CC) numrec.cpp   
SLogger.o: SLogger.cpp SLogger.h
	$(CC) SLogger.cpp
SLogWriter.o: SLogWriter.cpp SLogWriter.h
	$(CC) SLogWriter.cpp
ini.o: ini.c ini.h
	$(CC) ini.c
inireader.o: inireader.cpp inireader.h
	$(CC) inireader.cpp

#****************#
#compilation is done with the intell compiler, for more info see:
# 	http://www.intel.com/cd/software/products/asmo-na/eng/346158.htm
#site with tutorials and stuff:
# 	http://www.jaist.ac.jp/iscenter-new/mpc/altix/altixdata/opt/intel/cc/10.0.026/doc/main_cls/mergedProjects/optaps_cls/common/optaps_samples.htm
#how to compile MKL together with INTEL compilers, can be architecture different
#	http://docs.notur.no/uit/stallo_documentation/faq/test-faq-3
#to find stuff easily  find /opt -iname *.c    *.c or something else that you want
#****************#


#-----------------------
#icpc:  compiler
#-c     just compile file, don`t link stuff... 
#-xT    compile for Intel? Core? 2 Quad 6XXX series generate SSE4 Vectorizing   *see bottom page
#-o3    code/memory management more agressive                                   *see bottom page
#-ipo   the compiler performs inline function expansion for calls to functions defined in separate files.
#-no-prec-div   it enables optimizations for floating point divisions that give slightly less precise results than full IEEE division.
#-static        no real clue, some say to prevent crashes, I think doing this on runtime would be a waste off time. (to be investigated)
#-----------------------

#to link the stuff
#------------------
#-ssp           Software-based Speculative Pre-computation. Use it othwerise you get errors, like undefined references to __kmpc_end and __kmpc_begin
#		this kmpc errors because kmp_stacksize is too small if you use -parallell and -O3 simulatious, check man pages for -parallel
#------------------
