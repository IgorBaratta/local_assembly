#include "problem.h"
#include <iostream>
#include <ufc.h>
#include <chrono>

int main(int argc, char *argv[])
{
    int ncells = 100'000;

    ufc_form a = *form_problem_a;
    int ndofs_cell = a.finite_elements[1]->space_dimension;
    auto &&kernel_a = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;
    ufc_form L = *form_problem_L;
    if (ndofs_cell != L.finite_elements[0]->space_dimension)
        throw std::runtime_error("L/a");
    auto &&kernel_L = L.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

    const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0};

    double Ae[ndofs_cell * ndofs_cell];
    double be[ndofs_cell];
    double coefficients[ndofs_cell];

    for (int i = 0; i < ndofs_cell; i ++)
        coefficients[i] = 1.0;

    auto a_start = std::chrono::steady_clock::now();
    for (int c = 0; c < ncells; c++)
    {
        kernel_a(Ae, coefficients, nullptr, coordinate_dofs, 0, 0, 0);
    }
    auto a_end = std::chrono::steady_clock::now();
    auto L_start = std::chrono::steady_clock::now();
    for (int c = 0; c < ncells; c++)
    {
        kernel_L(be, coefficients, nullptr, coordinate_dofs, 0, 0, 0);
    }
    auto L_end = std::chrono::steady_clock::now();
    double duration_a = std::chrono::duration_cast<std::chrono::microseconds>(a_end - a_start).count() / 1.e6;
    double duration_L = std::chrono::duration_cast<std::chrono::microseconds>(L_end - L_start).count() / 1.e6;
    std::cout << ncells << ", " << duration_a << ", " << duration_L;
    return 0;
}