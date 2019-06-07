#-----------------------------------------------------------
# Macro definitions
CXXFLAGS = -O1 -g -Wall -std=c++17
OBJ = main.o tools.o CPM.o Group.o Subject.o predictory.o
TARGET = cpm 
#-----------------------------------------------------------
# Rules
all: $(TARGET)
$(TARGET): $(OBJ)
	$(CXX) -o $@ $(OBJ)
clean:
	rm -f $(OBJ) $(TARGET)
#-----------------------------------------------------------
# Dependencies
main.o: main.cpp tools.hpp CPM.hpp Group.hpp
tools.o: tools.cpp tools.hpp 
CPM.o: CPM.cpp CPM.hpp Group.hpp predictory.hpp
predictory.o: predictory.cpp predictory.hpp Group.hpp
RCPM.o: RCPM.cpp RCPM.hpp Group.hpp predictory.hpp
Group.o: Group.hpp Subject.hpp Group.cpp
Subject.o: Subject.hpp Subject.cpp

#Game.o: Game.hpp 
#FlipFlop.o: FlipFlop.hpp 
#ThinkADot.o: ThinkADot.hpp 
