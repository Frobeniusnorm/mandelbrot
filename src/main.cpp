/*
 *  Mandelbrot - a mandelbrot zoom application with SyCL
 *  Copyright (C) 2024  Sven Vollmar & David Schwarzbeck

 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "cmdparser.hpp"
#include "mandelbrot.hpp"
#include <string>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.hpp"

template <size_t bytes>
static void deep_zoom(FixedFloat<bytes> x1, FixedFloat<bytes> x2,
                 FixedFloat<bytes> y1, FixedFloat<bytes> y2, int i,
                 const int iterations, const FixedFloat<bytes> &zoom,
                 MandelbrotRenderer &renderer, std::string name, int width,
                 int height) {

  std::vector<unsigned char> &img_data =
      renderer.generate_image(x1, x2, y1, y2);
  stbi_write_jpg((name + "/" + std::to_string(i) + ".jpg").c_str(), width,
                 height, 3, img_data.data(), 100);
  std::cout << "generated " << name + "/" + std::to_string(i) + ".jpg"
            << std::endl;
  // reduce the size
  FixedFloat<bytes> x_space = (x2 - x1);
  FixedFloat<bytes> y_space = (y2 - y1);
  FixedFloat<bytes> x_space_new = (x2 - x1) * zoom;
  FixedFloat<bytes> y_space_new = (y2 - y1) * zoom;
  // adapt the start s.t. we zoom towards re and im
  // TODO if i come to implement division it would be helpful here
  auto xdiff = (x_space - x_space_new) * 0.5;
  x1 += xdiff;
  x2 -= xdiff;
  auto ydiff = (y_space - y_space_new) * 0.5;
  y1 += ydiff;
  y2 -= ydiff;
  if (i < iterations)
    deep_zoom(x1, x2, y1, y2, i + 1, iterations, zoom, renderer, name, width,
         height);
}

std::string generate_help() {
  using namespace std;
  const char *COLORS[] = {"\u001b[38;5;34m", "\u001b[38;5;35m",
                          "\u001b[38;5;36m", "\u001b[38;5;37m",
                          "\u001b[38;5;38m", "\u001b[38;5;39m"};
  const char *TEXT[] = {"\t>> ", "Man", "del", "brot ", "Zoom", " <<"};
  const string RESET = "\u001b[0m";
  string help = "\u001b[1m";
  for (int i = 0; i < sizeof(TEXT) / sizeof(char *); i++)
    help += string(COLORS[i]) + TEXT[i];
  help += RESET + "\n" + R"(
Mandelbrot - a mandelbrot zoom application with SyCL
Copyright (C) 2024  Sven Vollmar & David Schwarzbeck)";
  return help;
}

int main(int argc, char **argv) {
  cli::Parser parser(argc, argv, generate_help());
  parser.set_default(false, "output file or folder",
                     std::string("mandelbrot.jpg"));
  parser.set_optional("r", "real", -0.7746806106269039,
                      "Minimum of draw area on the x (real) axis");
  parser.set_optional("i", "imaginary", -0.1374168856037867,
                      "Minimum of draw area on the y (imaginary) axis");
  parser.set_optional("x", "width", 1920, "Image width");
  parser.set_optional("y", "height", 1080, "Image height");
  parser.set_optional(
      "n", "number", 1,
      "If set to 1 generates a single image, else it sets the number "
      "of images that should be generated");
  parser.set_optional("s", "speed", 0.95,
                      "Zoom factor (zoom per frame). Only relevant if the "
                      "number of iterations is set");
  parser.run();
  const int width = parser.get<int>("x");
  const int height = parser.get<int>("y");
  const double re = parser.get<double>("r");
  const double im = parser.get<double>("i");
  const double speed = parser.get<double>("s");
  const int iterations = parser.get<int>("n");
  const std::string name = parser.get_default<std::string>();
  std::cout << iterations << std::endl;
  double x1 = -2;
  double x2 = 1;
  double y1 = -1;
  double y2 = 1;
  MandelbrotRenderer renderer(width, height);
  if (iterations == 1) {
    std::vector<unsigned char> &img_data =
        renderer.generate_image(x1, x2, y1, y2);
    stbi_write_jpg(name.c_str(), width, height, 3, img_data.data(), 100);
  } else {
    for (int i = 0; i < iterations; i++) {
      std::vector<unsigned char> &img_data =
          renderer.generate_image(x1, x2, y1, y2);
      stbi_write_jpg((name + "/" + std::to_string(i) + ".jpg").c_str(), width,
                     height, 3, img_data.data(), 100);
      std::cout << "generated " << name + "/" + std::to_string(i) + ".jpg"
                << std::endl;
      // reduce the size
      double x_space = (x2 - x1);
      double y_space = (y2 - y1);
      double x_space_new = (x2 - x1) * speed;
      double y_space_new = (y2 - y1) * speed;
      // adapt the start s.t. we zoom towards re and im
      double dr = (re - x1) / (x2 - x1);
      double di = (im - y1) / (y2 - y1);
      x1 += (x_space - x_space_new) * dr;
      x2 -= (x_space - x_space_new) * (1 - dr);
      y1 += (y_space - y_space_new) * di;
      y2 -= (y_space - y_space_new) * (1 - di);
      if (x_space_new < 0.00001) {
        std::cout << "switching to BigFloat 16 bytes" << std::endl;
        auto ffspeed = FixedFloat<12>(speed);
        deep_zoom(FixedFloat<12>(x1), FixedFloat<12>(x2), FixedFloat<12>(y1),
             FixedFloat<12>(y2), i, iterations, ffspeed, renderer, name, width,
             height);
      }
    }
  }
}
