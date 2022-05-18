#pragma once

#include <unistd.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "mpi.h"

void reflection_method(const int n,
                       Matrix& A,
                       double** cols,
                       double* b,
                       int nrows,
                       int ncols,
                       int proc_num,
                       int rank,
                       int div_size);