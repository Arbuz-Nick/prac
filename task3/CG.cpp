#include "CG.h"

double now() {
  return omp_get_wtime();
}

double norm(const std::vector<double>& vec) {
  double res = 0.0;
  for (const auto& it : vec)
    res += it * it;
  return std::sqrt(res);
}

void spmv(const Matrix& A,
          const std::vector<double>& x,
          std::vector<double>& b) {
  unsigned int nrows = A.nrows;
  unsigned int width = A.row_size;
#pragma omp parallel for
  for (int i = 0; i < nrows; i++) {
    b[i] = 0.0;
    for (int j = 0; j < A.row_size; j++)
      b[i] += x[A.col[i * width + j]] * A.val[i * width + j];
  }
}

double dot(const std::vector<double>& x, const std::vector<double>& y) {
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
  for (unsigned int i = 0; i < x.size(); i++)
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
  double local_time, full_time = 0.0, dot_start, spmv_start, lin_start,
                     full_dot = 0.0, full_spmv = 0.0, full_lin = 0.0;

  bool convergence;

  for (int iter = 0; iter < run_iters; iter++) {
    convergence = false;

    double start = now();
    r = b;
    set(x, 0.0);
    k = 1;
    do {
      dot_start = now();
      ro = dot(r, r);
      full_dot += now() - dot_start;
      if (k == 1) {
        p = r;
      } else {
        betta = ro / ro_prev;
        lin_start = now();
        lin_comb(betta, p, 1.0, r);
      }

      spmv_start = now();
      spmv(A, p, q);
      full_spmv += now() - spmv_start;

      dot_start = now();
      double p2q = dot(p, q);
      full_dot += now() - dot_start;

      alpha = ro / p2q;

      lin_start = now();
      lin_comb(1.0, x, alpha, p);
      full_lin += now() - lin_start;

      lin_start = now();
      lin_comb(1.0, r, -alpha, q);
      full_lin += now() - lin_start;

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
  double err = error(x);
  std::cout << "Iter Num: " << k << std::endl;
  std::cout << "Full time: " << full_time << std::endl;
  std::cout << "Avg time: " << full_time / run_iters << std::endl;
  std::cout << "Avg spmv time: " << full_spmv / run_iters << std::endl;
  std::cout << "Avg dot time: " << full_dot / (2 * run_iters) << std::endl;
  std::cout << "Avg lin comb time: " << full_lin / (2 * run_iters) << std::endl;

  std::cout << "Residual: " << norm(r) << std::endl;
  std::cout << "||Ax-b|| = " << res << std::endl;
  std::cout << "Error: " << err << std::endl;
  std::ofstream result;
  result.open("result_omp_polus.csv", std::ios_base::app);
  result << size << ";" << omp_get_max_threads() << ";" << full_time / run_iters
         << ";" << full_spmv / run_iters << ";" << full_dot / (2 * run_iters) << ";" << full_lin / (2 * run_iters) << ";" << res << ";" << err << std::endl;
  result.close();
}
