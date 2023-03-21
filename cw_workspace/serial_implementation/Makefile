CXX = g++
CXXFLAGS = -Wall -O0 -g
HDRS = ShallowWater.h
LIBS = -lblas
TARGET = main
OBJS = main.o ShallowWater.o

default: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)


.PHONY: clean

clean:
	-rm -f *.o $(TARGET)
