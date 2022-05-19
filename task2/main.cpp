#include <vector>
#include "matrix.h"
#include "reflection_method.h"

int rank, proc_num;

static void generate_vector(const Matrix& A, double* b) {
  for (int i = 0; i < A.nrows; i++) {
    int sum = 0;
    for (int j = 0; j < A.ncols; j++)
      sum += A.data[i * A.nrows + j];
    b[i] = sum;
  }
}

int main(int argc, char* argv[]) {
  std::ofstream stat;
  stat.open("status.txt", std::ios_base::app);
  int r = 0;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Number of proccessors: " << proc_num << std::endl;
  }
  int n;

  if (argc > 1) {
    sscanf(argv[1], "%d", &n);
  } else {
    n = 10;
  }

  int div_size = n / proc_num;
  int rem = n % proc_num;
  for (int i = 0; i < rem; i++) {
    if (rank == i) {
      div_size++;
    }
  }

  Matrix A(n, n);
  double* b = new double[n];
  double** cols = new double*[div_size];
  for (int i = 0; i < n; i++) {
    b[i] = 0;
  }
  if (rank == 0) {
    A.generate();
    generate_vector(A, b);

    A.T();
    for (int i = 0; i < div_size; i++)
      cols[i] = A.data + i * proc_num * n;
  } else {
    for (int i = 0; i < div_size; i++) {
      cols[i] = new double[n];
      for (int j = 0; j < n; j++) {
        cols[i][j] = 0;
      }
    }
  }
  if (rank == 0) {
    stat << "START" << std::endl;
  }
  for (int i = 0; i < div_size; i++)
    MPI_Scatter(A.data + i * proc_num * n, n, MPI_DOUBLE, cols[i], n,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    stat << "SCATTER" << std::endl;
  }
  for (int i = 1; i < rem; i++) {
    if (rank == 0)
      MPI_Send(A.data + n * proc_num * (n / proc_num) + i * n, n, MPI_DOUBLE, i,
               0, MPI_COMM_WORLD);
    else if (rank == i)
      MPI_Recv(cols[div_size - 1], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
  }

  if (rank == 0) {
    stat << "SND/RCV1" << std::endl;
  }
  if (n % proc_num != 0) {
    if (rank == 0)
      MPI_Send(b, n, MPI_DOUBLE, rem, 0, MPI_COMM_WORLD);
    else if (rank == n % proc_num)
      MPI_Recv(b, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (rank == 0) {
    stat << "SND/RCV2" << std::endl;
  }
  if (rank == 0) {
    stat << "END" << std::endl;
  }
  int nrows = A.nrows, ncols = A.ncols;
  MPI_Barrier(MPI_COMM_WORLD);
  reflection_method(n, A, cols, b, nrows, ncols, proc_num, rank, div_size);
  delete[](b);
  if (rank != 0)
    for (int i = 0; i < div_size; i++) {
      delete[](cols[i]);
    }
  delete[](cols);
  MPI_Finalize();
  return 0;
}
