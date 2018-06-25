CXX=g++-6
RM=rm -f
CXXFLAGS=-std=c++11 -Wall -O -I/usr/include -pg -fopenmp
GSLFLAGS=-lgsl -lgslcblas -lm

SRCS=ljmain.cpp ljtps.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: ljmain

ljmain: ljmain.o ljtps.o
	$(CXX) $(CXXFLAGS) -o ljmain ljmain.o ljtps.o $(GSLFLAGS)

ljmain.o: ljmain.cpp ljtps.h
	$(CXX) $(CXXFLAGS) -c ljmain.cpp

ljtps.o: ljtps.h
	$(CXX) $(CXXFLAGS) -c ljtps.cpp

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) ljmain

anew: distclean all
