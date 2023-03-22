CXX = g++
CXXFLAGS = -Wall -O0 -g
HDRS = ShallowWater.h
LIBS = -lblas -lboost_program_options -fopenmp
TARGET = main
OBJS = main.o ShallowWater.o

default: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

profilerFor: $(TARGET)
	rm -r testing.er
	collect -o testing.er ./$(TARGET) --T 80 --ic 4 --intType 1
	analyzer testing.er

profilerBlas: $(TARGET)
	rm -r testing.er
	collect -o testing.er ./$(TARGET) --T 80 --ic 4 --intType 2
	analyzer testing.er

.PHONY: clean

clean:
	-rm -f *.o $(TARGET)
