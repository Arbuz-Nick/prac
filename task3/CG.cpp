#include "CG.h"

double now() {
  return omp_get_wtime();
}
/*
double norm(const std::vector<double>& vec, int size) {
  double res = 0;
  for (int i = 0; i < size; i++) {
    res += vec[i] * vec[i];
  }
  return std::sqrt(res);
}
*/
double norm(const std::vector<double>& vec) {
  double res = 0.0;
  for (const auto& it : vec)
    res += it * it;
  return std::sqrt(res);
}
/*
double sgn(double num) {
  if (num < 0)
    return (-1.0);
  else if (num > 0)
    return (1.0);
  return 0.0;
}
*/
void spmv(const Matrix& A,
          const std::vector<double>& x,
          std::vector<double>& b) {
  auto nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);

  uint32_t nrows = A.nrows;
  uint32_t ncols = A.ncols;
  uint32_t width = A.row_size;
#pragma omp parallel for
  for (int i = 0; i < nrows; i++) {
    b[i] = 0.0;
    for (int j = 0; j < A.row_size; j++)
      b[i] += x[A.col[i * width + j]] * A.val[i * width + j];
  }
}

double dot(const std::vector<double>& x, const std::vector<double>& y) {
  auto nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);

  double res = 0.0;

#pragma omp parallel for reduction(+ : res)
  for (int i = 0; i < x.size(); i++)
    res += x[i] * y[i];

  return res;
}

void lin_comb(const double a,
              std::vector<double>& x,
              const double b,
              const std::vector<double>& y) {
  auto nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for (int i = 0; i < y.size(); i++)
    x[i] = a * x[i] + b * y[i];
}

double residual(const Matrix& A,
                const std::vector<double>& x,
                const std::vector<double>& b) {
  std::vector<double> Ax(x.size(), 0.0);
  spmv(A, x, Ax);

  for (int i = 0; i < Ax.size(); i++) {
    Ax[i] -= b[i];
  }

  return norm(Ax);
}

double error(std::vector<double> x) {
  for (auto& it : x) {
    it -= 1;
  }
  return norm(x);
}

void set(std::vector<double>& x, double value) {
#pragma omp parallel for
  for (uint32_t i = 0; i < x.size(); i++)
    x[i] = value;
}

void CG(const int size,
        Matrix& A,
        std::vector<double>& b,
        int run_iters,
        int max_solver_iters = 100) {
  std::vector<double> x(size, 0.0);
  std::vector<double> r(size, 0.0);
  std::vector<double> p(size, 0.0);
  std::vector<double> q(size, 0.0);
  std::vector<double> z(size, 0.0);

  double ro = 0.0, ro_prev = 1.0;
  double betta = 1.0, alpha = 1.0;

  int k;
  double local_time, full_time = 0.0;

  bool convergence;

  for (int iter = 0; iter < run_iters; iter++) {
    convergence = false;

    double start = now();
    r = b;
    set(x, 0.0);
    k = 1;
    do {
      ro = dot(r, r);

      if (k == 1) {
        p = r;
      } else {
        betta = ro / ro_prev;
        lin_comb(betta, p, 1.0, r);
      }

      spmv(A, p, q);
      alpha = ro / dot(p, q);

      lin_comb(1.0, x, alpha, p);
      lin_comb(1.0, r, -alpha, q);
      if (norm(r) < 1e-12 || k >= max_solver_iters)
        convergence = true;
      else {
        k++;
        ro_prev = ro;
      }
    } while (!convergence);

    double end = now();

    local_time = end - start;
    full_time += local_time;
  }
  double res = residual(A, x, b);
  std::cout << "Iter Num: " << k << std::endl;
  double err = error(x);
  std::cout << "Full time: " << full_time << std::endl;
  std::cout << "Avg time: " << full_time / run_iters << std::endl;
  std::cout << "Residual: " << norm(r) << std::endl;
  std::cout << "||Ax-b|| = " << res << std::endl;
  std::cout << "Error: " << err << std::endl;
  std::ofstream result;
  result.open("result_omp_polus.csv", std::ios_base::app);
  result << size << ";" << omp_get_max_threads() << ";" << full_time / run_iters
         << ";" << res << ";" << err << std::endl;
  result.close();
}
