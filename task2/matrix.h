#pragma once

#include <stdlib.h>
#include <iostream>

struct Matrix {
 public:
  Matrix(int _nrows, int _ncols) : ncols(_ncols), nrows(_nrows) {}

  ~Matrix() { delete[] data; }
  void generate();
  void print();
  void T();
  double* data = nullptr;
  int nrows;
  int ncols;
};
