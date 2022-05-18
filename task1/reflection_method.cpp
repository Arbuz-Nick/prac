#include "reflection_method.h"

static inline double now() {
  return omp_get_wtime();
}

double norm(const std::vector<double>& vec, int size) {
  double res = 0;
  for (int i = 0; i < size; i++) {
    res += vec[i] * vec[i];
  }
  return std::sqrt(res);
}

double norm(const std::vector<double>& vec) {
  double res = 0.0;
  for (const auto& it : vec)
    res += it * it;
  return std::sqrt(res);
}

double sgn(double num) {
  if (num < 0)
    return (-1.0);
  else if (num > 0)
    return (1.0);
  return 0.0;
}

double residual(const Matrix& A,
                const std::vector<double>& x,
                const std::vector<double>& b) {
  std::vector<double> Ax(A.size, 0.0);

  for (int i = 0; i < A.size; i++)
    for (int j = 0; j < A.size; j++)
      Ax[i] += A.data[i * A.size + j] * x[j];

  std::vector<double> res(A.size, 0.0);

  for (int i = 0; i < A.size; i++) {
    res[i] = Ax[i] - b[i];
  }

  return norm(res);
}

void reflection_method(const int n, Matrix& A, std::vector<double>& b) {
  std::vector<double> a(n);

  std::vector<double> x(n, 0.0);

  A.T();
  double start_time = now();

  for (int i = 0; i < n - 1; i++) {
    int len = n - i;

#pragma omp parallel for
    for (int j = i; j < n; j++)
      a[j - i] = A.data[i * n + j];

    double a_norm = norm(a, len);
    a[0] -= sgn(a[0]) * a_norm;
    a_norm = norm(a, len);

#pragma omp parallel for
    for (int j = 0; j < len; j++)
      a[j] /= a_norm;

#pragma omp parallel for
    for (int col = i; col < n; col++) {
      double sum = 0.0;
      for (int k = 0; k < len; k++)
        sum += 2.0 * a[k] * A.data[col * n + i + k];
      for (int k = 0; k < len; k++)
        A.data[col * n + i + k] -= sum * a[k];
    }

    double sum = 0.0;
    for (int j = 0; j < len; j++)
      sum += 2.0 * a[j] * b[j + i];
    for (int j = 0; j < len; j++)
      b[j + i] -= sum * a[j];
  }

  double to_r_time = now();

  std::vector<double> start_b(n, 0.0);
  for (int i = 0; i < n; i++)
    start_b[i] = b[i];

  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++)
      b[i] -= A.data[j * n + i] * x[j];

    x[i] = b[i] / A.data[i * n + i];
  }

  double end_time = now();

  A.T();

  std::cout << "To R Matrix time: " << to_r_time - start_time << std::endl;
  std::cout << "Gauss time: " << end_time - to_r_time << std::endl;
  std::cout << "Full time: " << end_time - start_time << std::endl;
  std::cout << "Residual: " << residual(A, x, start_b) << std::endl;
  std::ofstream result;
  result.open("result_omp_polus.csv", std::ios_base::app);
  result << end_time - start_time << ";" << to_r_time - start_time << ";"
         << end_time - to_r_time << ";" << omp_get_max_threads() << ";" << n
         << ";" << residual(A, x, start_b) << std::endl;
  result.close();
}
