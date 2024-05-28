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
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.hpp"

std::string generate_help() {
  using namespace std;
  const char *COLORS[] = {"\u001b[38;5;34m", "\u001b[38;5;35m", "\u001b[38;5;36m",
                          "\u001b[38;5;37m", "\u001b[38;5;38m", "\u001b[38;5;39m"};
  const char *TEXT[] = {"\t>> ", "Man", "del", "brot ", "Zoom", " <<"};
  const string RESET = "\u001b[0m";
  string help = "\u001b[1m";
  for (int i = 0; i < sizeof(TEXT) / sizeof(char*); i++)
    help += string(COLORS[i]) + TEXT[i];
  help += RESET + "\n" + R"(
Mandelbrot - a mandelbrot zoom application with SyCL
Copyright (C) 2024  Sven Vollmar & David Schwarzbeck)";
  return help;
}

int main(int argc, char **argv) {
  cli::Parser parser(argc, argv, generate_help());
  parser.set_default(false, "output file", std::string("mandelbrot.jpg"));
  parser.set_optional("x1", "xmin", -2., "Minimum of draw area on the x (real) axis");
  parser.set_optional("x2", "xmax", 1., "Maximum of draw area on the x (real) axis");
  parser.set_optional("y1", "ymin", -1., "Minimum of draw area on the y (imaginary) axis");
  parser.set_optional("y2", "ymax", 1., "Maximum of draw area on the y (imaginary) axis");
  parser.set_optional("x", "width", 3000, "Image width");
  parser.set_optional("y", "height", 2000, "Image height");
  parser.run();
  const int width = parser.get<int>("x");
  const int height = parser.get<int>("y");
  const double x1 = parser.get<double>("x1");
  const double x2 = parser.get<double>("x2");
  const double y1 = parser.get<double>("y1");
  const double y2 = parser.get<double>("y2");
  const std::string name = parser.get_default<std::string>();
  MandelbrotRenderer renderer(width, height);
  std::vector<unsigned char>& img_data = renderer.generate_image(x1, x2, y1, y2);
  stbi_write_jpg(name.c_str(), width, height, 3, img_data.data(), 100);
}
