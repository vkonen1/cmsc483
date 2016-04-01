CC		= g++
CCFLAGS	= -g -ansi -Wall
OBJECTS = seam_carver.o

all: $(OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) -o seam_carver -lX11 -lpthread
	rm *.o

seam_carver.o: seam_carver.cpp
	$(CC) $(CCFLAGS) -c seam_carver.cpp

clean:
	rm -f *.o seam_carver *~ *.gch