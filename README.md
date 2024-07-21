# Mandelbrot
Mandelbrot Zoom Application with SyCL

## Building
You need [AdaptiveCPP](https://github.com/AdaptiveCpp/AdaptiveCpp/blob/develop/doc/installing.md) to compile the application.
After installation the application `acpp` should be available.
The application can then be build by running the `Makefile` which generates the `mandelbrot` application in the
root of the project. It needs an up-to-date version of whatever API you want to be used for the calculations
(i.e. HIP, CUDA or OpenMP). OpenCL should work for systems with OpenCL 3.0 support.

## Usage
The application has a help option that explains all arguments and flags (`mandelbrot -h`).
In general it is either possible to generate a zoom or a single image. When generating
an image the fixed argument is the filename that should be generated, else it is the folder in which
all the images of the zoom are to be generated. `ffmpeg` can be used to generate a video from the
single images.

This SyCL implementation uses an approximation formula to estimate the needed number of iterations for each
image. Check out the `cpu_only` branch that contains an OpenMP implementation that uses an approximation
by measuring how many pixel reach the upper quartile of the iteration count to either increase or decrease the
needed number of iterations.

## Used libraries
- SyCL
- [CmdParser](https://github.com/FlorianRappl/CmdParser/blob/master/cmdparser.hpp) (MIT License)
- [stb_image_write](https://github.com/nothings/stb/blob/master/stb_image_write.h) (MIT License)
