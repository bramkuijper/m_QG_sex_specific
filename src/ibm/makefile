# makefile to build an executable

# name of the executable we will build
EXE=m_matpat.exe

# name of the C++ file
CPP=m_quantgen_matpat.cpp

# compiler command
CXX=g++

# some flags to compile the thing
CXXFLAGS=-Wall -O3 -std=c++17

# then build the actual thing
$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)


clean :
	rm -rf $(EXE)
