MPICC	= mpic++
CC		= g++
CCFLAGS	= -g -ansi -Wall
LIBS	= -lX11 -lpthread

all: $(OBJECTS)
	$(MPICC) $(CCFLAGS) -o seam_carver seam_carver.cpp -lm

imgproc:
	$(CC) $(CCFLAGS) -o imgproc imgproc.cpp $(LIBS)	

energy:
	$(MPICC) $(CCFLAGS) -o parallel_energy_seam_carver parallel_energy_seam_carver.cpp -lm

seam:
	$(MPICC) $(CCFLAGS) -o parallel_seam_seam_carver parallel_seam_seam_carver.cpp -lm

clean:
	rm -f *.o seam_carver imgproc parallel_energy_seam_carver parallel_seam_seam_carver *~ *.gch