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
#include "bigfloat.hpp"
#include <sycl/sycl.hpp>

using namespace sycl;
static size_t max_iterations(double xlims) {
  // no mathematical proof, just approximation
  return 50 + (long)sycl::pow(sycl::log10(((4. / xlims))), 5);
}
template <size_t bytes>
void mandelbrot(queue &Q, double x_min, double x_max, double y_min,
                double y_max, std::vector<unsigned char> &res, size_t res_width,
                size_t res_height) {
  buffer image_buffer(res.data(), range<1>{res.size()});
  Q.submit([&](auto &h) {
    accessor img(image_buffer, h, write_only, no_init);
    h.parallel_for(res.size() / 3, [=](item<1> i) {
      // somehow it breaks if a fixedfloat is used more than once?
      const FixedFloat<bytes> two = FixedFloat<bytes>(2.);
      const FixedFloat<bytes> barrier = FixedFloat<bytes>((double)(1 << 16));
      const int yi = i / res_width; // screen space in [0, res_width]
      const int xi = i % res_width; // screen space in [0, res_height]
      const double x = xi / (double)(res_width - 1);  // screen space in [0, 1]
      const double y = yi / (double)(res_height - 1); // screen space in [0, 1]
      const FixedFloat<bytes> x0(x * sycl::fabs(x_max - x_min) +
                                 x_min); // complex plane in [x_min, x_max]
      const FixedFloat<bytes> y0 = (y * sycl::fabs(y_max - y_min) +
                                    y_min); // complex plane in [y_min, y_max]

      FixedFloat<bytes> c_x(0.0);
      FixedFloat<bytes> c_y(0.0);
      const size_t max_iter = max_iterations(sycl::fabs(x_max - x_min));
      // simulate complex conjecture
      FixedFloat<bytes> cx_squared;
      FixedFloat<bytes> cy_squared;
      FixedFloat<bytes> newx;
      FixedFloat<bytes> newy;
      size_t iter = 0;
      for (; iter < max_iter; iter++) {
        cx_squared = c_x * c_x;
        cy_squared = c_y * c_y;
        if (cx_squared + cy_squared > barrier)
          break;
        newx = cx_squared - cy_squared + x0;
        newy = two * c_x * c_y + y0;
        c_x = newx;
        c_y = newy;
      }
      // apply smoothing
      if (iter < max_iter) {
        double rcx = *c_x;
        double rcy = *c_y;
        const double log_zn = sycl::log(rcx * rcx + rcy * rcy) / 2;
        const double nu = sycl::log(log_zn / sycl::log(2.)) / sycl::log(2.);
        const double iters = iter + 1 - nu;
        // map to color map, higher iterations -> higher index
        const double progress = iters / max_iter;
        img[i * 3] =
            (char)(sycl::clamp(sycl::clamp(1.0 - 1.8 * progress, 0., 1.) +
                                   2 * progress * progress,
                               0., 1.) *
                   255);
        img[i * 3 + 1] = (char)((1.0 - 0.4 * progress) * 256);
        img[i * 3 + 2] =
            (char)(sycl::clamp(1.0 - sycl::sqrt(progress) + progress, 0., 1.) *
                   255);
      } else {
        // not diverging
        for (int j = 0; j < 3; j++) {
          img[i * 3 + j] = 0;
        }
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
  mandelbrot<16>(Q, x_min, x_max, y_min, y_max, working_image, res_width,
                 res_height);
  return working_image;
}
