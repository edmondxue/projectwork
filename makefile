CXX := g++
CXXFLAGS :=-std=c++11  -Wall -Wconversion -Wextra -Wpedantic 

TARGET:= main
OBJS := main.o COO2CSR.o CGSolver.o matvecops.o
INCS := COO2CSR.hpp CGSolver.hpp matvecops.hpp

$(TARGET): $(OBJS)
	 $(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
	 

%.o: %.cpp $(INCS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) 
 
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET) *~
