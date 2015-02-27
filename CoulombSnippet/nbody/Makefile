# Makefile 

NVCC = nvcc 	
NVCCFLAGS = -ccbin gcc --compiler-options -fno-strict-aliasing,-Wall,-Wextra,-Wno-unused-parameter -O3  
LIBS= -lcudart -lcuda
CXX = nvcc
CXXFLAGS = -O3  -I/usr/local/cuda/include


OBJS = test_nbody.o nbody_kernel.o mtrand.o
TARGET = test_nbody

LINKLINE  = nvcc $(NVCCFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


%.o: %.cu 
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
	
%.o: %.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
	
$(TARGET): $(OBJS) Makefile
	$(LINKLINE)
		
clean:
	find . -name \*.o -exec rm {} \;
	-rm -f $(EXECUTABLES)

.PHONY: clean
