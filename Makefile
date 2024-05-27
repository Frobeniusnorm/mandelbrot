CXXC := acpp
CXXFLAGS := -Og -g 

CPP_SOURCES := $(wildcard src/*.cpp)
HPP_SOURCES := $(wildcard src/*.hpp)
CPP_OBJECTS := $(patsubst src/%.cpp,build/%.o,$(CPP_SOURCES))

.PHONY: build

mandelbrot: $(CPP_OBJECTS) | build
	$(CXXC) -o $@ $(CPP_OBJECTS) 

build/%.o: src/%.cpp $(HPP_SOURCES) | build
	$(CXXC) $(CXXFLAGS) -c -o $@ $<

build:
	- mkdir build

clean:
	- rm build/*
	- rm mandelbrot
