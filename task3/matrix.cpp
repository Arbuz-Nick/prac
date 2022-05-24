#include "matrix.h"

#include <math.h>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

static inline double get_element(int i, int j) {
  return cos(i * j + M_PI);
}

void Matrix::print() const {
  if (is_empty)
    return;

  //    std::cout << std::scientific << std::setprecision(3);

  std::cout << "Matrix shape: " << nrows << " x " << ncols << std::endl;
  std::cout << "Maximum nonzeros per rows: " << row_size << std::endl;
  std::cout << "Values: " << std::endl;
  for (int i = 0; i < nrows; i++) {
    std::cout << "row " << std::setw(3) << i << " : ";
    for (int j = 0; j < row_size; j++)
      std::cout << std::setw(5) << (uint64_t)col[i * row_size + j] << ":"
                << std::setw(10) << val[i * row_size + j] << " ";
    std::cout << std::endl;
  }
}

void Matrix::generate(int _size) {
  std::cout << "Generation started\n";

  //    val.resize(nrows * row_size);
  //    col.resize(nrows * row_size);

  size = _size;
  // row_size = 7;

  for (int K = 0; K < size; K++) {
    for (int J = 0; J < size; J++) {
      for (int I = 0; I < size; I++) {
        int i = K * size * size + J * size + I;
        int k;
        int j = 0;
        int diag_position = 0;
        if (K > 0) {
          k = (K - 1) * size * size + J * size + I;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        if (J > 0) {
          k = K * size * size + (J - 1) * size + I;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        if (I > 0) {
          k = K * size * size + J * size + I - 1;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        val.push_back(0.0);
        col.push_back(K * size * size + J * size + I);
        diag_position = j;
        j++;
        if (I + 1 != size) {
          k = K * size * size + J * size + I + 1;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        if (J + 1 != size) {
          k = K * size * size + (J + 1) * size + I;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        if (K + 1 != size) {
          k = (K + 1) * size * size + J * size + I;
          val.push_back(get_element(i, k));
          col.push_back(k);
          j++;
        }
        for (int jj = j; jj < row_size; jj++) {
          val.push_back(0.0);
          col.push_back(0);
        }
        double sum = 0.0;
        for (int jj = 0; jj < row_size; jj++)
          sum += abs(val[i * row_size + jj]);
        val[i * row_size + diag_position] = 1.5 * sum;
      }
    }
  }
  is_empty = false;
}

std::vector<double> Matrix::get_diag() const {
  std::vector<double> result(nrows, 0.0);

  for (int i = 0; i < nrows; i++) {
    int j = 0;
    while (col[i * row_size + j] != i) {
      assert(j < row_size);
      j++;
    }
    result[i] = val[i * row_size + j];
  }

  return result;
}
