CXXC :=g++
CXXFLAGS := -O3 -fopenmp

CPP_SOURCES := $(wildcard src/*.cpp)
HPP_SOURCES := $(wildcard src/*.hpp)
CPP_OBJECTS := $(patsubst src/%.cpp,build/%.o,$(CPP_SOURCES))

.PHONY: build

mandelbrot: $(CPP_OBJECTS) | build
	$(CXXC) $(CXXFLAGS) -o $@ $(CPP_OBJECTS) 

build/%.o: src/%.cpp $(HPP_SOURCES) | build
	$(CXXC) $(CXXFLAGS) -c -o $@ $<

build:
	- mkdir build

clean:
	- rm build/*
	- rm mandelbrot
