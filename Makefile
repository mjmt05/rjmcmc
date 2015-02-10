CXX=g++
CXXFLAGS=-Wall -Wno-long-long -pedantic -march=native -O3
INCLUDES=-I/usr/include/gsl 
LDLIBS=-lgsl -lgslcblas -lm
HEADERS=$(wildcard *.h)

.PHONY: clean

all: mainRJ_example mainSMC_example

mainRJ_example: mainRJ_example.cpp $(HEADERS)

mainSMC_example: mainSMC_example.cpp $(HEADERS)

clean:
	rm -f mainRJ_example mainSMC_example
