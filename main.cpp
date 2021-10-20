
#include <ufc.h>
#include "problem.h"
#include <iostream>
#include <chrono>

int main(int argc, char *argv[])
{
    int ncells = 100000;

    problem_cell_integral_0_otherwise a;
    
    int ndofs_cell = 10000;

    const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0};

    double Ae[ndofs_cell * ndofs_cell];
    double coefficients[ndofs_cell];
    const double *coeff_ref = &coefficients[0];

    for (int i = 0; i < ndofs_cell; i ++)
        coefficients[i] = 1.0;

    auto start = std::chrono::steady_clock::now();
    for (int c = 0; c < ncells; c++)
    {
      a.tabulate_tensor(Ae, &coeff_ref,  coordinate_dofs, 0);
    }
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
    std::cout << ncells << ", " << duration;
    return 0;
}
