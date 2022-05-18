#include <vector>
#include "matrix.h"
#include "reflection_method.h"


static void generate_vector(int n, std::vector<double>& b) {
  for (int i = 0; i < n; i++)
    b[i] = (double)((2 * i) % n) / n;
}


int main(int argc, char const* argv[]) {
  int n;
  
  if (argc > 1) {
    sscanf(argv[1], "%d", &n);
  } else {
    n = 10;
  }

  Matrix A(n);

  A.generate();

  std::vector<double> b(n);

  generate_vector(n, b);

  reflection_method(n, A, b);

  return 0;
}