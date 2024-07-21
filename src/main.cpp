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
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.hpp"
static struct timeval timestamp_start, timestamp_stop;

inline bool exists(const std::string &name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}
template <size_t bytes>
static void
deep_zoom(FixedFloat<bytes> x1, FixedFloat<bytes> x2, FixedFloat<bytes> y1,
          FixedFloat<bytes> y2, int i, const int iterations,
          const FixedFloat<bytes> &zoom, const FixedFloat<bytes> &ffre,
          const FixedFloat<bytes> &ffim, MandelbrotRenderer &renderer,
          std::string name, int width, int height) {
  // constants
  const FixedFloat<bytes> ffone(1.0);
  const FixedFloat<bytes> coeff(2 << (bytes / 2));
  std::string path = name + "/";
  {
    std::stringstream iss;
    iss << std::setw(log(iterations) / log(10.0) + 1) << std::setfill('0') << i;
    std::string is = iss.str();
    path += is + ".jpg";
  }
  // only generate images that don't have been calculated already
  if (!exists(path)) {
    gettimeofday(&timestamp_start, nullptr);
    std::vector<unsigned char> &img_data =
        renderer.generate_image(x1, x2, y1, y2);
    stbi_write_jpg(path.c_str(), width, height, 3, img_data.data(), 100);
    gettimeofday(&timestamp_stop, nullptr);
    long millis = timestamp_stop.tv_usec - timestamp_start.tv_usec;
    std::cout << "generated " << path << " in "
              << (millis > 10000 ? std::to_string(millis / 1000) + "s"
                                 : std::to_string(millis) + "ms")
              << std::endl;
  } else
    std::cout << "skipping " << path << std::endl;
  // reduce the size
  FixedFloat<bytes> x_space = (x2 - x1);
  FixedFloat<bytes> y_space = (y2 - y1);
  FixedFloat<bytes> x_space_new = (x2 - x1) * zoom;
  FixedFloat<bytes> y_space_new = (y2 - y1) * zoom;
  // adapt the start s.t. we zoom towards re and im
  // TODO if i come to implement division it would be helpful here
  auto dr = FixedFloat<bytes>(*(coeff * (ffre - x1)) / *(coeff * (x2 - x1)));
  auto di = FixedFloat<bytes>(*(coeff * (ffim - y1)) / *(coeff * (y2 - y1)));
  auto xdiff = x_space - x_space_new;
  auto ydiff = y_space - y_space_new;

  x1 += (x_space - x_space_new) * dr;
  x2 -= (x_space - x_space_new) * (ffone - dr);
  y1 += (y_space - y_space_new) * di;
  y2 -= (y_space - y_space_new) * (ffone - di);
  if (i < iterations) {
    deep_zoom(x1, x2, y1, y2, i + 1, iterations, zoom, ffre, ffim, renderer,
              name, width, height);
  }
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
typedef FixedFloat<16> FF16;
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
  FF16 x1 = -2;
  FF16 x2 = 1;
  FF16 y1 = -1;
  FF16 y2 = 1;
  const FF16 ffre(re);
  const FF16 ffim(im);
  const FF16 ffone(1.0);
  const FF16 coeff(2 << 8);
  const FF16 ffspeed(speed);
  MandelbrotRenderer renderer(width, height);
  if (iterations == 1) {
    // generate one image
    std::vector<unsigned char> &img_data =
        renderer.generate_image(*x1, *x2, *y1, *y2);
    stbi_write_jpg(name.c_str(), width, height, 3, img_data.data(), 100);
  } else {
    // generate iterations images, when the precision is too low for the zoom,
    // switch to 16-byte bigfloats.
    for (int i = 0; i < iterations; i++) {
      std::string path = name + "/";
      {
        std::stringstream iss;
        iss << std::setw(log(iterations) / log(10.0) + 1) << std::setfill('0')
            << i;
        std::string is = iss.str();
        path += is + ".jpg";
      }
      if (!exists(path)) {
        gettimeofday(&timestamp_start, nullptr);
        std::vector<unsigned char> &img_data =
            renderer.generate_image(*x1, *x2, *y1, *y2);
        gettimeofday(&timestamp_stop, nullptr);
        stbi_write_jpg(path.c_str(), width, height, 3, img_data.data(), 100);
        long starttime =
            timestamp_start.tv_sec * 1000 + timestamp_start.tv_usec / 1000;
        long stoptime =
            timestamp_stop.tv_sec * 1000 + timestamp_stop.tv_usec / 1000;
        long millis = stoptime - starttime;
        std::cout << "generated " << path << " in "
                  << (millis > 10000000l
                          ? std::to_string(millis / 1000000l) + "s"
                          : std::to_string(millis / 1000) + "ms")
                  << " (x-size: " << (*x2 - *x1) << ")" << std::endl;
      } else
        std::cout << "skipping " << path << std::endl;
      // reduce the size
      FF16 x_space = (x2 - x1);
      FF16 y_space = (y2 - y1);
      FF16 x_space_new = (x2 - x1) * ffspeed;
      FF16 y_space_new = (y2 - y1) * ffspeed;
      // adapt the start s.t. we zoom towards re and im
      auto dr = FF16(*(coeff * (ffre - x1)) / *(coeff * (x2 - x1)));
      auto di = FF16(*(coeff * (ffim - y1)) / *(coeff * (y2 - y1)));
      x1 += (x_space - x_space_new) * dr;
      x2 -= (x_space - x_space_new) * (ffone - dr);
      y1 += (y_space - y_space_new) * di;
      y2 -= (y_space - y_space_new) * (ffone - di);
      if (x_space_new < 3.5e-12) {
        std::cout << "switching to BigFloat 16 bytes" << std::endl;
        deep_zoom(FixedFloat<16>(*x1), FixedFloat<16>(*x2), FixedFloat<16>(*y1),
                  FixedFloat<16>(*y2), i + 1, iterations, ffspeed, ffre, ffim,
                  renderer, name, width, height);
      }
    }
  }
}
