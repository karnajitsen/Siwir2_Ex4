CC=g++
CFLAGS= -Wall -std=c++11 -pedantic -g
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=mdsim.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mdsim
COMMON=

all: clean test

test:
	$(CC) $(CFLAGS) $(SOURCES) -o mdsim
	
mdsim:
	rm -rf ./output/*.*
	./mdsim params.dat

ref0:
	./waveguide 0.01 0.0000000001 0
	rm -rf *.pdf
	gnuplot ./plot.p
test: 
	$(CC) $(CFLAGS) $(SOURCES) -o waveguide
	./waveguide 0.01 0.0000000001

clean:
	rm -rf lbm
	
.PHONY : all clean
