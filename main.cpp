#include "problem.h"
#include <boost/align/aligned_allocator.hpp>
#include <chrono>
#include <iostream>
#include <ufc.h>
#include <vector>

int main(int argc, char *argv[]) {
  int ncells = 100000;

  ufc_form a = *form_problem_a;
  int ndofs_cell = a.finite_elements[1]->space_dimension;
  auto &&kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

  const std::vector<double, boost::alignment::aligned_allocator<double, 64>>
      coordinate_dofs = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1,
                         0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::vector<double, boost::alignment::aligned_allocator<double, 64>> Ae(
      ndofs_cell * ndofs_cell);

  std::vector<double, boost::alignment::aligned_allocator<double, 64>>
      coefficients(ndofs_cell);

  for (int i = 0; i < ndofs_cell; i++)
    coefficients[i] = 1.0;

  auto start = std::chrono::steady_clock::now();
  for (int c = 0; c < ncells; c++) {
    kernel(Ae.data(), coefficients.data(), nullptr, coordinate_dofs.data(), 0,
           0, 0);
  }
  auto end = std::chrono::steady_clock::now();
  double duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      1.e6;
  std::cout << ncells << ", " << duration;
  return 0;
}
