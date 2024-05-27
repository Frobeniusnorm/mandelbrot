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

#include "mandelbrot.hpp"
#include <sycl/sycl.hpp>

using namespace sycl;
const unsigned char COLOR_SET1[][3] = {
	{255, 255, 255},
	{240, 255, 230},
	{230, 255, 220},
	{200, 255, 230},
	{150, 250, 200},
	{100, 250, 150},
	{80, 240, 150},
	{60, 230, 180},
	{50, 210, 190},
	{50, 200, 200},
	{50, 190, 210},
	{50, 190, 230},
	{50, 190, 240},
	{50, 150, 240},
	{50, 130, 230}
};
static size_t max_iterations(double xlims) {
  // no mathematical proof, just approximation
  return 50 + (long)sycl::pow(sycl::log10(((4. / xlims))), 5);
}

void mandelbrot(queue &Q, double x_min, double x_max, double y_min,
                double y_max, std::vector<unsigned char> &res, size_t res_width,
                size_t res_height) {
  buffer image_buffer(res.data(), range<1>{res.size()});
  Q.submit([&](auto &h) {
    accessor img(image_buffer, h, write_only, no_init);
    h.parallel_for(res.size() / 3, [=](item<1> i) {
      const int yi = i / res_width; // screen space in [0, res_width]
      const int xi = i % res_width; // screen space in [0, res_height]
      const double x = xi / (double)(res_width - 1);  // screen space in [0, 1]
      const double y = yi / (double)(res_height - 1); // screen space in [0, 1]
      const double x0 = x * sycl::fabs(x_max - x_min) +
                        x_min; // complex plane in [x_min, x_max]
      const double y0 = y * sycl::fabs(y_max - y_min) +
                        y_min; // complex plane in [y_min, y_max]
      double c_x = 0;
      double c_y = 0;
      const size_t max_iter = max_iterations(sycl::fabs(x_max - x_min));
      // simulate complex conjecture
      size_t iter = 0;
      for (; iter < max_iter && c_x * c_x + c_y * c_y <= 2 * 2; iter++) {
        const double new_x = c_x * c_x - c_y * c_y + x0;
        const double new_y = 2 * c_x * c_y + y0;
        c_x = new_x;
        c_y = new_y;
      }
      // map to color map, higher iterations -> higher index
      int color_idx = (int)((iter / (double)max_iter) *
                            (sizeof(COLOR_SET1) / (sizeof(char) * 3)));
      for (int j = 0; j < 3; j++) {
        img[i * 3 + j] = COLOR_SET1[color_idx][j];
      }
    });
  });
  Q.wait();
  host_accessor h_acc(image_buffer);
  for (int i = 0; i < res.size(); i++)
    res[i] = h_acc[i];
}

std::vector<unsigned char> &MandelbrotRenderer::generate_image(double x_min,
                                                               double x_max,
                                                               double y_min,
                                                               double y_max) {
  default_selector device_selector;
  queue Q(device_selector);
  mandelbrot(Q, x_min, x_max, y_min, y_max, working_image, res_width,
             res_height);
  return working_image;
}
