# Makefile
#MPICHDIR= /Users/pdupre/mpich2-1.0.8
CXXFLAGS = -O3 -I $(CUDA_HOME)/include -Wall
cunbodylib=../libcunbody1.a_TI20130122_CUDA5.0
LDFLAGS  = $(cunbodylib) -lm  -L $(CUDA_HOME)/lib64 -lcudart -lcuda 
EXEC= test
CXX= mpicxx
SRC =  mpi_funcs.cpp main.cpp cloud.cpp MultiGPU_cunbody.cpp
OBJ= $(SRC:.cpp=.o)


$(EXEC): $(OBJ)
	$(CXX) -o $@  $^ $(LDFLAGS)
%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<







.PHONY: clean mrproper

clean::
	-rm -f *.o *.exe

# end
