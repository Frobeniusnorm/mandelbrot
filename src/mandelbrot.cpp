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

static size_t max_iterations(double xlims) {
	// no mathematical proof, just approximation
	return 50 + pow(log10(((4. / xlims))), 5);
}
static float COLOR_SET1[][3] = {
    {0.8, 0.6, 0.4},  {0.85, 0.7, 0.5}, {0.6, 0.85, 0.5}, {0.5, 0.85, 0.6},
    {0.3, 0.85, 0.8}, {0.3, 0.7, 0.85}, {0.3, 0.6, 1}};

std::vector<float> MandelbrotRenderer::generate_image(double x_min,
                                                      double x_max,
                                                      double y_min,
                                                      double y_max) {
  using namespace cl::sycl;
  Q.submit([&](auto &h) {
    accessor img(image_buffer, h, write_only, no_init);
    h.parallel_for(working_image.size(), [=](auto i) {
      const int xi = (i / 3) / res_height; // screen space in [0, res_width]
      const int yi = (i / 3) % res_height; // screen space in [0, res_height]
      const double x = xi / (double)(res_width - 1);  // screen space in [0, 1]
      const double y = yi / (double)(res_height - 1); // screen space in [0, 1]
      const double x0 =
          x * std::abs(x_max - x_min) + x_min; // complex plane in [x_min, x_max]
      const double y0 =
          y * std::abs(y_max - y_min) + y_min; // complex plane in [y_min, y_max]
      double c_x = 0;
      double c_y = 0;
      const size_t max_iter = max_iterations(std::abs(x_max - x_min));
      // simulate complex conjecture
      size_t iter = 0;
      for (; iter < max_iter && c_x * c_x + c_y * c_y <= 2 * 2; iter++) {
        const double new_x = c_x * c_x + c_y * c_y + x0;
        const double new_y = 2 * c_x * c_y + y0;
        c_x = new_x;
        c_y = new_y;
      }
      // map to color map, higher iterations -> higher index
      int color_idx = (int)((iter / (double)max_iter) * sizeof(COLOR_SET1) /
                            (sizeof(float) * 3));
      for (int j = 0; j < 3; j++)
        img[i + j] = COLOR_SET1[color_idx][j];
    });
  });
  host_accessor h_acc(image_buffer);
  for (int i = 0; i < working_image.size(); i++)
    working_image[i] = h_acc[i];
  return working_image;
}
