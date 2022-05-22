#include "reflection_method.h"

static inline double now() {
  return MPI_Wtime();
}

double norm(double* vec, int size) {
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

double residual(const Matrix& A, double* x, double* b) {
  std::vector<double> Ax(A.nrows, 0.0);
  for (int i = 0; i < A.nrows; i++) {
    for (int j = 0; j < A.ncols; j++) {
      Ax[i] += A.data[i * A.ncols + j] * x[j];
    }
  }

  std::vector<double> res(A.nrows, 0.0);
  for (int i = 0; i < A.nrows; i++) {
    res[i] = Ax[i] - b[i];
  }
  return norm(res);
}

double error(double* x, int size) {
  for (int i = 0; i < size; i++)
    x[i] -= 1;
  return norm(x, size);
}

void reflection_method(const int n,
                       Matrix& A,
                       double** cols,
                       double* b,
                       int nrows,
                       int ncols,
                       int proc_num,
                       int rank,
                       int div_size) {
  double* a = new double[n];
  for (int i = 0; i < n; ++i) {
    a[i] = 0;
  }
  double* start_b = nullptr;
  if (rank == 0) {
    start_b = new double[n];
    for (int i = 0; i < n; ++i)
      start_b[i] = 0;
  }

  double* x = new double[n];
  double* x_global = new double[n];
  for (int i = 0; i < n; ++i) {
    x[i] = 0;
    x_global[i] = 0;
  }

  double start_time = now();

  for (int i = 0; i < nrows - 1; i++) {
    int len = nrows - i;
    if (rank == i % proc_num) {
      for (int j = i; j < nrows; j++)
        a[j - i] = cols[i / proc_num][j];

      double a_norm = norm(a, len);
      a[0] -= sgn(a[0]) * a_norm;
      a_norm = norm(a, len);

      for (int j = 0; j < len; j++)
        a[j] /= a_norm;
    }

    MPI_Bcast(a, len, MPI_DOUBLE, i % proc_num, MPI_COMM_WORLD);

    for (int col = i / proc_num; col < div_size; col++) {
      if (col * proc_num + rank >= i) {
        double sum = 0.0;
        for (int k = 0; k < len; k++)
          sum += 2.0 * a[k] * cols[col][i + k];
        for (int k = 0; k < len; k++)
          cols[col][i + k] -= sum * a[k];
      }
    }

    if (rank == n % proc_num) {
      double sum = 0.0;
      for (int j = 0; j < len; j++)
        sum += 2.0 * a[j] * b[j + i];
      for (int j = 0; j < len; j++)
        b[j + i] -= sum * a[j];
    }
  }

  double to_r_time = now();

  if (n % proc_num != 0) {
    if (rank == n % proc_num) {
      MPI_Send(b, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else if (rank == 0) {
      MPI_Recv(start_b, n, MPI_DOUBLE, n % proc_num, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  } else if (rank == 0) {
    for (int i = 0; i < n; i++)
      start_b[i] = b[i];
  }
  for (int i = n - 1; i >= 0; i--) {
    if (rank == i % proc_num) {
      if (proc_num > 1)
        MPI_Recv(b, n, MPI_DOUBLE, (i + 1) % proc_num, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      x[i] = b[i] / cols[i / proc_num][i];
      for (int j = 0; j <= i; j++) {
        b[j] -= x[i] * cols[i / proc_num][j];
      }
    } else if (rank == (i + 1) % proc_num) {
      if (proc_num > 1) {
        MPI_Send(b, n, MPI_DOUBLE, i % proc_num, 0, MPI_COMM_WORLD);
      }
    }
  }

  MPI_Allreduce(x, x_global, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double end_time = now();

  if (proc_num > 1)
    for (int i = 0; i < n / proc_num; i++)
      MPI_Gather(cols[i], n, MPI_DOUBLE, A.data + i * proc_num * n, n,
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int i = 1; i < n % proc_num; i++) {
    if (rank == i)
      MPI_Send(cols[div_size - 1], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    else if (rank == 0)
      MPI_Recv(A.data + proc_num * n * (n / proc_num) + i * n, n, MPI_DOUBLE, i,
               0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (rank == 0) {
    A.T();
    double resid = residual(A, x_global, start_b);
    double err = error(x_global, n);

    std::cout << "To R Matrix time: " << to_r_time - start_time << std::endl;
    std::cout << "Gauss time: " << end_time - to_r_time << std::endl;
    std::cout << "Full time: " << end_time - start_time << std::endl;
    std::cout << "Residual: " << resid << std::endl;
    std::cout << "Error: " << err << std::endl;

    std::ofstream result;
    result.open("result_mpi_polus.csv", std::ios_base::app);

    result << end_time - start_time << ";" << to_r_time - start_time << ";"
           << end_time - to_r_time << ";" << proc_num << ";" << n << ";"
           << resid << ";" << err << std::endl;
    result.close();
    delete[](start_b);
  }
  delete[](x);
  delete[](x_global);
  delete[](a);
}
