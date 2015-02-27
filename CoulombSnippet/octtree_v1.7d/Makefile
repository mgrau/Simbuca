# 
#icpc -g -debug -O3 -ipo -no-prec-div -o runsim main.o parser.o ionFly.o ode.o force.o 
#libcunbody1.a -L $CUDASYS/lib -L $CUDASYS/NVIDIA_CUDA_SDK/lib -lcuda -lcudart -lGL -lGLU -lcutil coll.o 
#fieldmap.o matrix.o  ioncloud.o particle.o ion.o mtrand.o logfile.o
# 
#
#	First example of a makefile: made by ME!
#libcunbody1.a
#
all: aout

clean: 
	rm -f *.o aout

aout : sample.o
	g++ -O3 -o aout sample.o liboctgrav.a -L $(CUDASYS)/lib64 -lcuda -lcudart -lGL -lGLU
sample.o : sample.c liboctgrav.a
	 g++ -O3 -c -I/usr/local/cuda/include -I$(CUDASYS)/C/common/inc sample.c -I../../libraries/octgrav_v1.7d

