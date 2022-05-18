#include "matrix.h"

void Matrix::generate() {

  data = new double[size * size];
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      data[i * size + j] = (i+1)/(j+1)+10;//(double)(i + 1) / (j + 1) + 10;
}

void Matrix::print() {
  std::cout << "Matrix size: " << size << "x" << size << std::endl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      std::cout << data[i * size + j] << " ";
    std::cout << std::endl;
  }
}

void Matrix::T() {
  double* res = new double[size * size];
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      res[i * size + j] = data[j * size + i];
  delete(data);
  data = res;
}