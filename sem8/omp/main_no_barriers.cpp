#include "omp.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.h"

static inline double now()
{
    return omp_get_wtime();
}

Matrix matmul(Matrix &A, Matrix &B)
{
    int size = A.size;
    Matrix C(size);
    C.generate();
#pragma omp parallel
    for (int i = 0; i < size; i++)
    {
        for (int k = 0; k < size; k++)
        {
#pragma omp for nowait
            for (int j = 0; j < size; j++)
            {
                C.data[i * size + j] += A.data[i * size + k] * B.data[k * size + j];
            }
        }
    }
    //C.print();
    return C;
}

int main(int argc, char const *argv[])
{
    int n;

    if (argc > 1)
    {
        sscanf(argv[1], "%d", &n);
    }
    else
    {
        n = 10;
    }

    Matrix A(n);
    A.generate();
    Matrix B(n);
    B.generate();
    // Matrix C;

    double start_time = now();
    matmul(A, B);
    double end_time = now();
    std::ofstream result;
    result.open("result_omp_polus.csv", std::ios_base::app);
    result << end_time - start_time << ";" << omp_get_max_threads() << ";" << n << ";no_barriers" << std::endl;
    result.close();
    // C.print();
    return 0;
}
