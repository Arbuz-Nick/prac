#include <vector>
#include "CG.h"
#include "matrix.h"

std::vector<double> generate_vector(const Matrix& A) {
  std::vector<double> b(A.nrows, 0.0);

  for (unsigned int i = 0; i < A.nrows; i++) {
    double sum = 0.0;
    for (unsigned int j = 0; j < A.row_size; j++)
      sum += A.val[i * A.row_size + j];
    b[i] = sum;
  }

  return b;
}

int main(int argc, char const* argv[]) {
  int n = 10, max_iter = 100, run_iter = 50;

  if (argc > 1)
    sscanf(argv[1], "%d", &n);
  if (argc > 2)
    sscanf(argv[2], "%d", &max_iter);
  if (argc > 3)
    sscanf(argv[3], "%d", &run_iter);

  int size = n * n * n;
  Matrix A(size);

  double start_generation = now();
  A.generate(n);
  double end_generation = now();
  std::cout << "Generation time: " << end_generation - start_generation
            << std::endl;

  std::vector<double> b = generate_vector(A);

  CG(size, A, b, run_iter, max_iter);
  return 0;
}