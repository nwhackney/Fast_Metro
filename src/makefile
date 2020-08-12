CXX = g++
CXXFLAGS = -std=c++14 -Wall  -g

GSLLIBS = -L/usr/local/lib -lgsl -lgslcblas
GSLINCS = -I/usr/local/include

main: main.o lattice.o
	$(CXX) $(CXXFLAGS) -o main main.o lattice.o $(GSLINCS) $(GSLLIBS)

main.o:
	$(CXX) $(CXXFLAGS) -c main.cpp

lattice.o:
	$(CXX) $(CXXFLAGS) -c lattice.cpp

clean:
	rm *.o
	
