#---------------------------------------------#
#------------USER INPUT PARAMETERS------------#
#---------------------------------------------#

#compiler:  mpic++ / g++ / icpc
COMPILER=g++

IDIR=../../bin

EXECUTABLE=exe_parsefile_average

#--------------------------------------------#



#---------COMPILE AND LINK OPTIONS-----------#
CC_gpp=g++ -O3 -Wno-write-strings -D__CPUNBODY_ON__ -c
LINK_gpp=g++ -lm -O3 -o 

CC_icpc=icpc -O3 -ipo -no-prec-div -D__CPUNBODY_ON__ -c
LINK_icpc=icpc -O3 -ipo -no-prec-div -o

#MPICHDIR= /Users/pdupre/mpich2-1.0.8/bin
CC_mpicpp=mpicxx  -O3 -D__MPI_ON__ -finline-functions -Wno-write-strings  -DMPICH_IGNORE_CXX_SEEK  $(NBODY_ALGO_FLAGS) -c 
#LINK_mpicpp=g++ -O3  -L/Users/pdupre/mpich2/lib -lmpichcxx -lpmpich -lmpich -Wwrite-strings -std -lm -o 
LINK_mpicpp=mpicxx -O3   -std -lm -o 
#--------------------------------------------#


#---------------SET COMPILER-----------------#
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





#------------BUILDING DEPENDENCIES------------#
all: $(EXECUTABLE) 
    

clean: 
	rm -f *.o $(EXECUTABLE)

delink:
	rm $(EXECUTABLE)
	
parameters:
	@echo COMPILER=$(COMPILER)
	@echo PU=$(PU)
	@echo EXECUTABLE=$(EXECUTABLE)

export:
	mkdir -p ../backup
	tar cvfps ../backup/Simbuca_v$(VERSION).`date +"%Y%m%d-%H%M%S"`.tar ../bin/Makefile ../bin/*.cpp ../bin/*.h ../bin/*.a ../bin/Changelog.txt \
	../doc/* \
	../CoulombSnippet/* \
	../cunbody_source/* \
	../libraries/*

copy:
	cp main.cpp main_bckup.cpp



$(EXECUTABLE) : main.o parser.o ionFly.o ode.o force.o coll.o ImageCharges.o CreateRecoilDistr.o fieldcalc.o fieldmap.o matrix.o ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o trapparameters.o numrec.o Potmap.o Potmap.o SLogger.o SLogWriter.o
	$(LINK) $(EXECUTABLE) main.o parser.o ionFly.o ode.o force.o coll.o ImageCharges.o CreateRecoilDistr.o fieldcalc.o fieldmap.o matrix.o ioncloud.o particle.o ion.o mtrand.o logfile.o IonTable.o PDGTable.o trapparameters.o numrec.o Potmap.o SLogger.o SLogWriter.o


main.o : main.cpp $(IDIR)/parser.h $(IDIR)/ionFly.h $(IDIR)/ode.h $(IDIR)/ImageCharges.h $(IDIR)/force.h $(IDIR)/coll.h $(IDIR)/CreateRecoilDistr.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h $(IDIR)/Potmap.h $(IDIR)/SMsgType.h $(IDIR)/SLogger.h $(IDIR)/SLogWriter.h
	$(CC) main.cpp 
parser.o : $(IDIR)/parser.cpp $(IDIR)/parser.h $(IDIR)/ionFly.h $(IDIR)/ode.h $(IDIR)/ImageCharges.h $(IDIR)/force.h $(IDIR)/coll.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h $(IDIR)/SMsgType.h $(IDIR)/SLogger.h $(IDIR)/SLogWriter.h
	$(CC) $(IDIR)/parser.cpp 
ionFly.o : $(IDIR)/ionFly.cpp $(IDIR)/ionFly.h $(IDIR)/ode.h $(IDIR)/ImageCharges.h $(IDIR)/force.h $(IDIR)/coll.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h $(IDIR)/SMsgType.h $(IDIR)/SLogger.h $(IDIR)/SLogWriter.h
	$(CC) $(IDIR)/ionFly.cpp
ode.o : $(IDIR)/ode.cpp $(IDIR)/ode.h $(IDIR)/force.h $(IDIR)/coll.h $(IDIR)/ImageCharges.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h $(IDIR)/Potmap.h $(IDIR)/SMsgType.h $(IDIR)/SLogger.h $(IDIR)/SLogWriter.h
	$(CC) $(IDIR)/ode.cpp
force.o : $(IDIR)/force.cpp $(IDIR)/force.h $(IDIR)/coll.h $(IDIR)/ImageCharges.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h $(IDIR) $(IDIR)/Potmap.h $(IDIR)/SMsgType.h $(IDIR)/SLogger.h $(IDIR)/SLogWriter.h
	$(CC) $(IDIR)/force.cpp
coll.o : $(IDIR)/coll.cpp $(IDIR)/coll.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/coll.cpp
ImageCharges.o: $(IDIR)/ImageCharges.cpp $(IDIR)/ImageCharges.h $(IDIR)/fieldcalc.h $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h 
	$(CC) $(IDIR)/ImageCharges.cpp
CreateRecoilDistr.o : $(IDIR)/CreateRecoilDistr.cpp $(IDIR)/CreateRecoilDistr.h
	$(CC) $(IDIR)/CreateRecoilDistr.cpp
fieldcalc.o : $(IDIR)/fieldcalc.cpp $(IDIR)/fieldcalc.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/fieldcalc.cpp
fieldmap.o : $(IDIR)/fieldmap.cpp $(IDIR)/fieldmap.h $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/fieldmap.cpp
matrix.o : $(IDIR)/matrix.cpp $(IDIR)/matrix.h $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/matrix.cpp
ioncloud.o : $(IDIR)/ioncloud.cpp $(IDIR)/ioncloud.h $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/ioncloud.cpp
particle.o : $(IDIR)/particle.cpp $(IDIR)/particle.h $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/particle.cpp
ion.o : $(IDIR)/ion.cpp $(IDIR)/ion.h $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/ion.cpp
IonTable.o: $(IDIR)/IonTable.cpp $(IDIR)/IonTable.h
	$(CC) $(IDIR)/IonTable.cpp	
PDGTable.o: $(IDIR)/PDGTable.cpp $(IDIR)/PDGTable.h
	$(CC) $(IDIR)/PDGTable.cpp
mtrand.o : $(IDIR)/mtrand.cpp $(IDIR)/mtrand.h $(IDIR)/logfile.h
	$(CC) $(IDIR)/mtrand.cpp
logfile.o: $(IDIR)/logfile.cpp $(IDIR)/logfile.h
	$(CC) $(IDIR)/logfile.cpp
trapparameters.o: $(IDIR)/trapparameters.cpp $(IDIR)/trapparameters.h
	$(CC) $(IDIR)/trapparameters.cpp   
numrec.o: $(IDIR)/numrec.cpp $(IDIR)/numrec.h
	$(CC) $(IDIR)/numrec.cpp   
Potmap.o: $(IDIR)/Potmap.cpp $(IDIR)/Potmap.h
	$(CC) $(IDIR)/Potmap.cpp
SLogger.o: $(IDIR)/SLogger.cpp $(IDIR)/SLogger.h	
	$(CC) $(IDIR)/SLogger.cpp
SLogWriter.o: $(IDIR)/SLogWriter.cpp $(IDIR)/SLogWriter.h
	$(CC) $(IDIR)/SLogWriter.cpp   

