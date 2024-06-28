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
#include "bigfloat.hpp"
#include <cmath>
#include <vector>
struct MandelbrotRenderer {
  unsigned int res_width;
  unsigned int res_height;

  MandelbrotRenderer(unsigned int res_width, unsigned int res_height)
      : res_width(res_width), res_height(res_height),
        working_image(res_width * res_height * 3, 0) {}

  std::vector<unsigned char> &generate_image(double x_min, double x_max,
                                             double y_min, double y_max) {
    mandelbrot_double(x_min, x_max, y_min, y_max, working_image, res_width,
                      res_height);
    return working_image;
  }

  template <size_t bytes>
  std::vector<unsigned char> &
  generate_image(FixedFloat<bytes> x_min, FixedFloat<bytes> x_max,
                 FixedFloat<bytes> y_min, FixedFloat<bytes> y_max) {
    // TODO choose between openmp and mandelbrot_new_gpu
    mandelbrot_openmp(x_min, x_max, y_min, y_max, working_image, res_width,
                      res_height);
    return working_image;
  }

private:
  static size_t max_iterations(double xlims) {
    // no mathematical proof, just approximation
    // was: return 50 + (int)sycl::pow(sycl::log10(((256. / xlims))), 3);
    return 50 + (int)std::pow(std::log10(((256. / xlims))), 5);
  }
  static void mandelbrot_double(double x_min, double x_max, double y_min,
                                double y_max, std::vector<unsigned char> &res,
                                size_t res_width, size_t res_height) {
    const size_t max_iter = max_iterations(std::fabs(x_max - x_min));
    std::cout << "image generation with " << max_iter << " iterations... "
              << std::flush;
#pragma omp parallel for
    for (int i = 0; i < res.size() / 3; i++) {
      const int yi = i / res_width; // screen space in [0, res_width]
      const int xi = i % res_width; // screen space in [0, res_height]
      const double x = xi / (double)(res_width - 1);  // screen space in [0, 1]
      const double y = yi / (double)(res_height - 1); // screen space in [0, 1]
      const double x0(x * std::fabs(x_max - x_min) +
                      x_min); // complex plane in [x_min, x_max]
      const double y0(y * std::fabs(y_max - y_min) +
                      y_min); // complex plane in [y_min, y_max]
      double c_x = 0;
      double c_y = 0;
      double cx_squared(0), cy_squared(0), newx(0), newy(0);
      const double two(2.), barrier(1 << 16);
      // simulate complex conjecture
      size_t iter = 0;
      for (; iter < max_iter; iter++) {
        cx_squared = c_x * c_x;
        cy_squared = c_y * c_y;
        if (cx_squared + cy_squared > barrier)
          break;
        newx = c_x * c_x - c_y * c_y + (x0);
        newy = two * c_x * c_y + (y0);
        c_x = newx;
        c_y = newy;
      }
      // apply smoothing
      if (iter < max_iter) {
        const double log_zn = std::log((c_x * c_x + c_y * c_y)) / 2;
        const double nu = std::log(log_zn / std::log(2.)) / std::log(2.);
        const double iters = iter + 1 - nu;
        // map to color map, higher iterations -> higher index
        const double progress = iters / max_iter;

        res[i * 3] = 0;
        res[i * 3 + 1] = progress * 255;
        res[i * 3 + 2] = progress * 255;
      } else {
        // not diverging
        for (int j = 0; j < 3; j++) {
          res[i * 3 + j] = 0;
        }
      }
    }
    std::cout << "... finished" << std::endl;
  }
  template <size_t bytes>
  static void
  mandelbrot_openmp(FixedFloat<bytes> x_min, FixedFloat<bytes> x_max,
                    FixedFloat<bytes> y_min, FixedFloat<bytes> y_max,
                    std::vector<unsigned char> &res, size_t res_width,
                    size_t res_height) {
    const size_t max_iter = max_iterations(*(x_max - x_min));
    std::cout << "image generation with " << max_iter << " iterations... "
              << std::flush;
#pragma omp parallel for
    for (size_t i = 0; i < res_width * res_height; i++) {
      const int yi = i / res_width; // screen space in [0, res_width]
      const int xi = i % res_width; // screen space in [0, res_height]
      const double x = xi / (double)(res_width - 1);  // screen space in [0, 1]
      const double y = yi / (double)(res_height - 1); // screen space in [0, 1]
      const FixedFloat<bytes> x0 =
          FixedFloat<bytes>(x) * (x_max - x_min).abs() +
          x_min; // complex plane in [x_min, x_max]
      const FixedFloat<bytes> y0 =
          FixedFloat<bytes>(y) * (y_max - y_min).abs() +
          y_min; // complex plane in [y_min, y_max]
      FixedFloat<bytes> c_x = 0;
      FixedFloat<bytes> c_y = 0;
      FixedFloat<bytes> cx_squared(0), cy_squared(0), newx(0), newy(0);
      const FixedFloat<bytes> two(2.), barrier(1 << 16);
      // simulate complex conjecture
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
        const double log_zn = std::log(*(c_x * c_x + c_y * c_y)) / 2;
        const double nu = std::log(log_zn / std::log(2.)) / std::log(2.);
        const double iters = iter + 1 - nu;
        // map to color map, higher iterations -> higher index
        const double progress = iters / max_iter;
        res[i * 3] =
            (char)(std::clamp(std::clamp(1.0 - 1.8 * progress, 0., 1.) +
                                  2 * progress * progress,
                              0., 1.) *
                   255);
        res[i * 3 + 1] = (char)((1.0 - 0.4 * progress) * 256);
        res[i * 3 + 2] =
            (char)(std::clamp(1.0 - std::sqrt(progress) + progress, 0., 1.) *
                   255);
      } else {
        // not diverging
        for (int j = 0; j < 3; j++) {
          res[i * 3 + j] = 0;
        }
      }
    };
    std::cout << "... finished" << std::endl;
  }
  std::vector<unsigned char> working_image;
};
#endif
