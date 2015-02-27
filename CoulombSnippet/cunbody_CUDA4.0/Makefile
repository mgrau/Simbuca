## Makefile for the snippet exampler for cunbody 
#
#icpc -g -debug -O3 -ipo -no-prec-div -o runsim main.o parser.o ionFly.o ode.o force.o libcunbody1.a -L $CUDASYS/lib -L $CUDASYS/NVIDIA_CUDA_SDK/lib -lcuda -lcudart -lGL -lGLU -lcutil coll.o fieldmap.o matrix.o  ioncloud.o particle.o ion.o mtrand.o logfile.o
##

#cunbodylib=libcunbody1.a_TI20111214f_CUDA4.0
cunbodylib=libcunbody1.a

CUDASYS=/usr/local/cuda

all: aout

clean: 
	rm -f *.o aout

aout : sample.o libcunbody1.a
	g++ -g -O3 -o aout sample.o $(cunbodylib) -L $(CUDASYS)/lib64 -lcuda -lcudart -lGL -lGLU
sample.o : sample.c libcunbody1.a
	g++ -O3 -c -g sample.c
