CXX = g++
CXXFLAGS = -std=c++14 -Wall

main: main.o lattice.o
	$(CXX) $(CXXFLAGS) -o main main.o lattice.o

main.o:
	$(CXX) $(CXXFLAGS) -c main.cpp

lattice.o:
	$(CXX) $(CXXFLAGS) -c lattice.cpp

