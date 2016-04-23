MPICC	= mpic++
CC		= g++
CCFLAGS	= -g -ansi -Wall
LIBS	= -lX11 -lpthread

all: $(OBJECTS)
	$(MPICC) $(CCFLAGS) -o seam_carver seam_carver.cpp -lm
	$(CC) $(CCFLAGS) -o imgproc imgproc.cpp $(LIBS)

clean:
	rm -f *.o seam_carver imgproc *~ *.gch