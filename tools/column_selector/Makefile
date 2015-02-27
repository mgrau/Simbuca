EXECUTABLE=columnselect
CC=icpc
CFLAGS= -pthread -m64
LIBS= 

all: $(EXECUTABLE)

$(EXECUTABLE): 
	$(CC) $(CFLAGS) $(LIBS) -o $(EXECUTABLE) columnselect.cpp

clean:  
	rm $(EXECUTABLE)
