#include <iostream>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>

#define MSG_TAG 0
#define MSG_COUNT 10
double bench_t_start, bench_t_end;

static double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday(&Tp, NULL);
    if (stat != 0)
        printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
    bench_t_start = rtclock();
}

void bench_timer_stop()
{
    bench_t_end = rtclock();
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    int msg_size;
    sscanf(argv[2], "%d", &msg_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *send_buf = new int[msg_size];
    int *recv_buf = new int[msg_size * msg_size];

    for (int i = 0; i < msg_size; i++)
    {
        send_buf[i] = rank;
    }
    MPI_Request *send_reqs = new MPI_Request[msg_size];
    MPI_Request *recv_reqs = new MPI_Request[msg_size];
    MPI_Status *statuses = new MPI_Status[msg_size];
    if (rank == 0)
    {
        bench_timer_start();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int k = -1;
    for (int i = 0; i < msg_size; i++)
    {

        MPI_Send_init(send_buf, msg_size, MPI_INT, (rank + size - 1 - (i % (size - 1))) % size, MSG_TAG, MPI_COMM_WORLD, &send_reqs[i]);
        MPI_Recv_init(&recv_buf[msg_size * i], msg_size, MPI_INT, (rank + 1 + (i % (size - 1))) % size, MSG_TAG, MPI_COMM_WORLD, &recv_reqs[i]);
    }

    MPI_Startall(msg_size, send_reqs);
    MPI_Startall(msg_size, recv_reqs);

    MPI_Waitall(msg_size, send_reqs, statuses);
    MPI_Waitall(msg_size, recv_reqs, statuses);

    if (rank == 0)
    {
        bench_timer_stop();
        double search_time = bench_t_end - bench_t_start;
        std::ofstream out_file;
        std::string path = "result_mpi_";
        path += argv[1];
        path += ".csv";
        out_file.open(path, std::ios_base::app);
        out_file << search_time << ";" << size << ";"
                 << "normal"
                 << ";" << argv[2] << std::endl;
        out_file.close();
    }
    /*
    if (rank == 1)
    {
        for (int i = 0; i < msg_size; i++)
        {
            for (int j = 0; j < msg_size; j++)
            {
                std::cerr << recv_buf[msg_size * i + j] << ' ';
            }
            std::cerr << std::endl;
        }
    }
    */
    delete[] recv_buf;
    delete[] send_reqs;
    delete[] recv_reqs;
    delete[] statuses;
    MPI_Finalize();
    return 0;
}