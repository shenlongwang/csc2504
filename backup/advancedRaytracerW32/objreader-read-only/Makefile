
CC = g++
FLAGS = -Wall -g
OBJECTS = meshObject.o

all: readObj

readObj: $(OBJECTS)
	$(CC) -o $@ readObj.cpp $(OBJECTS) $(FLAGS)

meshObject.o: meshObject.cpp meshObject.h
	$(CC) -c -o $@ $< $(FLAGS)

