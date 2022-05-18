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

void reflection_method(const int n, Matrix& A, std::vector<double>& b);