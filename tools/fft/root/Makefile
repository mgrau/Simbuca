EXECUTABLE=fft
CC=icpc
CFLAGS= -pthread -m64 -I/usr/local/root/include
LIBS= -L/usr/local/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic


all: $(EXECUTABLE)

$(EXECUTABLE): 
	$(CC) $(CFLAGS) $(LIBS) -o $(EXECUTABLE) fft.cpp

clean:  
	rm fft
