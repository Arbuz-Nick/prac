#include <iostream>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>

#define MSG_TAG 0
#define MSG_COUNT 5
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

    MPI_Request req_send, req_recv;

    if (rank == 0)
    {
        bench_timer_start();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int k = -1;
    for (int i = 0; i < msg_size; i++)
    {

        /*
        0 1 2 3 5 0 1 2 3
           0 1 2 3 4 5 6 7 8
           0 1 2 3 0 1 2 3 0
        0: 1 2 3 1 2 3 1 2 3
        1: 2 3 0 2 3 0 2 3 0
        2: 3 0 1 3 0 1 3 0 1
        3: 0 1 2 0 1 2 0 1 2
        */
        
        MPI_Isend(send_buf, msg_size, MPI_INT, (rank + size - 1 - (i % (size - 1))) % size, MSG_TAG, MPI_COMM_WORLD, &req_send);
        MPI_Irecv(&recv_buf[msg_size * i], msg_size, MPI_INT, (rank + 1 + (i % (size - 1))) % size, MSG_TAG, MPI_COMM_WORLD, &req_recv);

        MPI_Wait(&req_send, MPI_STATUS_IGNORE);
        MPI_Wait(&req_recv, MPI_STATUS_IGNORE);
        for (int j = 0; j < msg_size; j++)
        {
            if (recv_buf[msg_size * i + j] != (rank + 1 + (i % (size - 1))) % size)
            {
                std::cerr << "Error: received value " << recv_buf[msg_size * (i < rank ? i : i - 1) + j] << " from process " << i << std::endl;
                break;
            }
        }
    }
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
                 << "problem"
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

    MPI_Finalize();
    return 0;
}
