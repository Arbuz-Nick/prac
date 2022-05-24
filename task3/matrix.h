#pragma once

#include <vector>
#include <cstdint>

struct Matrix {
public:

    std::vector<double> val;
    std::vector<int> col;

    bool is_empty = true;
    int size = 0;
    int nrows, ncols = 0;
    int row_size = 7;
    
    Matrix(int _size) {
        nrows = _size;
        ncols = _size;
    }

    void print() const;
    void generate(int);
    std::vector<double> get_diag() const;
};