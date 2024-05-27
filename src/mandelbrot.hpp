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
#ifndef MANDELBROT_HPP
#define MANDELBROT_HPP
#include <vector>
struct MandelbrotRenderer {
  unsigned int res_width;
  unsigned int res_height;

  MandelbrotRenderer(unsigned int res_width, unsigned int res_height)
      : res_width(res_width), res_height(res_height),
        working_image(res_width * res_height * 3, 0) {}

  std::vector<unsigned char>& generate_image(double x_min, double x_max, double y_min,
                                    double y_max);

private:
  std::vector<unsigned char> working_image;
};
#endif
