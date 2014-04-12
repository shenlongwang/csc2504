CC = g++
CFLAGS = #-g -O2 
LIBS = -lm #-fopenmp

all:raytracer 

meshObject.o: meshObject.cpp meshObject.h
	$(CC) -c -o $@ $< $(CFLAGS)

raytracer:  raytracer.o util.o photon.o light_source.o scene_object.o bmp_io.o meshObject.o
	$(CC) $(CFLAGS) -o raytracer \
	raytracer.o util.o photon.o light_source.o scene_object.o bmp_io.o meshObject.o $(LIBS)

clean:
	-rm -f core *.o
	-rm raytracer
	



