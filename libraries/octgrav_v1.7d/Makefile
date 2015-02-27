CXX = g++
CC  = gcc
LD  = g++
F90  = ifort

.SUFFIXES: .o .f90 .cpp .cu

#CUDA_SDK = $(CUDASDK)
CUDA_SDK = /usr/local/cuda/NVIDIA_GPU_Computing_SDK
CUDA_TK  = /usr/local/cuda

OFLAGS = -O0 -g -Wall 
OFLAGS = -O3 -g -mpreferred-stack-boundary=4 -funroll-loops -fforce-addr
CFLAGS = $(OFLAGS)  -I$(CUDA_TK)/include -I$(CUDA_SDK)/C/common/inc
CXXFLAGS = $(CFLAGS)

NVCC      = $(CUDA_TK)/bin/nvcc
NVCCFLAGS = -D_DEBUG \
	 --maxrregcount=32  -O0 -g -I$(CUDA_SDK)/C/common/inc

PROG = nbint
OCTLIB = octgrav

LIBOBJS = load_bodies.o \
	evaluate_gravity.o \
	reorder_bodies.o \
	memory.o \
	device_data.o \
	CUDA_evaluate_gravity.o

OBJS = $(PROG).o 

CUOBJS = host_evaluate_gravity.cu_o

#LIBS = -L$(CUDA_TK)/lib64 -lcuda -lcudart
LIBS = -L$(CUDASYS/lib64 -lcuda -lcudart

all: $(OCTLIB) $(PROG) 

$(PROG): $(OBJS)
	$(LD) $(CXXFLAGS) $(LIBS) $^ -o $@ -L. -l$(OCTLIB)

$(OCTLIB): $(LIBOBJS) $(CUOBJS)
	ar qv lib$@.a $^
	ranlib lib$@.a

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

.f90.o:
	$(F90) $(F90FLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

convert:
	ifort -O3 -o convert convert.f90

clean:
	/bin/rm -rf *.o *.cu_o
	/bin/rm -rf $(PROG) lib$(OCTLIB).a

$(PROG).o:  octgrav.h
load_data.o:  octgrav.h
build_balanced_tree.o: octgrav.h
compute_multipole_moments.o: octgrav.h
sort_leaves.o: octgrav.h
reorder_bodies.o: octgrav.h
evaluate_gravity.o: octgrav.h
memory.o: octgrav.h
prepare_data_for_device.o: octgrav.h
copy.o: octgrav.h
CUDA_evaluate_gravity.o: octgrav.h

host_evaluate_gravity.cu_o: dev_evaluate_gravity.cu dev_octgrav_tex.cuh 
