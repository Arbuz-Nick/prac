#include "matrix.h"

void Matrix::generate() {
  data = new double[nrows * ncols];
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      data[i * nrows + j] = (i + 1) / (j + 1) + 10;
}

void Matrix::print() {
  std::cout << "Matrix size: " << nrows << "x" << ncols << std::endl;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++)
      std::cout << data[i * nrows + j] << " ";
    std::cout << std::endl;
  }
}

void Matrix::T() {
  double* res = new double[nrows * ncols];
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      res[i * nrows + j] = data[j * ncols + i];
  delete[] data;
  data = res;
}