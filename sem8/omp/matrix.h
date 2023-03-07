#pragma once

#include <stdlib.h>
#include <iostream>

struct Matrix {
 
 public:

  Matrix(int _size) : size(_size) {}
  
  ~Matrix() { delete(data); }
  void generate();
  void print();
  void T();
  double* data = nullptr;
  int size;
};
