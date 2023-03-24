CXX = g++
CXXFLAGS = -Wall -O3 -g
HDRS = ShallowWater.h
LIBS = -lblas -lboost_program_options -fopenmp
TARGET = main
OBJS = main.o ShallowWater.o

default: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

test1: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --intType 2

test2: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --intType 2

test3: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --intType 2

test4: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --intType 2


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
