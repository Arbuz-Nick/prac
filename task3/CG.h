#pragma once

#include <unistd.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "omp.h"

double now();
void CG(const int size,
        Matrix& A,
        std::vector<double>& b,
        int run_iters,
        int max_solver_iters);