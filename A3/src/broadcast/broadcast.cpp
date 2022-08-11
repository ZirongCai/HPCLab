#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <mpi.h>

#define N 10
#define LOOP 50

using namespace std::chrono;
//sends an array containing n doubles from root process 0
//to all other processes within the ring
void my_bcast(void *data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator)
{
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

    if (world_rank == root)
    {
        int i;
        for (i = 0; i < world_size; i++)
        {
            if (i != world_rank)
            {
                MPI_Send(data, count, datatype, i, 0, communicator);
            }
        }
    }
    else
    {
        MPI_Recv(data, count, datatype, root, 0, communicator,
                 MPI_STATUS_IGNORE);
    }
}

void broadcast_tree(void *data, int count, MPI_Datatype datatype, int root_proc, int end_proc, MPI_Comm communicator)
{
    if (root_proc == end_proc)
    {
        return;
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int recieve_proc = root_proc + (end_proc - root_proc + 1) / 2;
    if (rank == root_proc)
    {
        MPI_Send(data, count, datatype, recieve_proc, 0, communicator);
    }
    else if (rank == recieve_proc)
    {
        MPI_Recv(data, count, datatype, root_proc, 0, communicator, MPI_STATUS_IGNORE);
    }
    broadcast_tree(data, count, datatype, root_proc, recieve_proc - 1, communicator);
    broadcast_tree(data, count, datatype, recieve_proc, end_proc, communicator);
}

void broadcast_buffer(void *data, int count, MPI_Datatype datatype, int root,
                      MPI_Comm communicator)
{
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

    int buffer_size = 4096;
    int packet_num = (count / buffer_size); // it should be count/buffer_size, but i don't know the buffer size
    int remainder = count % buffer_size;
    //MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == root)
    {
        for (size_t j = 0; j < packet_num; j++)
        {
            MPI_Send((double *)data + (j * buffer_size), buffer_size, datatype, (world_rank + 1) % world_size, 0, communicator);
        }
        MPI_Send((double *)data + (packet_num*buffer_size), remainder, datatype, (world_rank + 1) % world_size, 0, communicator);
    }
    else
    {
        for (size_t j = 0; j < packet_num; j++)
        {
            MPI_Recv((double *)data + (j * buffer_size), buffer_size, datatype, (world_rank +world_size - 1)% world_size, 0, communicator, MPI_STATUS_IGNORE);
            MPI_Send((double *)data + (j * buffer_size), buffer_size, datatype, (world_rank + 1) % world_size, 0, communicator);
        }
        MPI_Recv((double *)data + (packet_num * buffer_size), remainder, datatype, (world_rank +world_size - 1)% world_size, 0, communicator, MPI_STATUS_IGNORE);
        MPI_Send((double *)data + (packet_num * buffer_size), remainder, datatype, (world_rank + 1) % world_size, 0, communicator);
    }
}

void my_bcast_tree(void *data, int count, MPI_Datatype datatype, int root,
                   MPI_Comm communicator)
{
    int world_size;
    MPI_Comm_size(communicator, &world_size);
    broadcast_tree(data, count, datatype, root, world_size - 1, communicator);
}

bool check(double *data, double *checksum)
{
    for (int i = 0; i < N; i++)
    {
        if (data[i] != checksum[i])
        {
            return false;
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size;
    int rank;
    size_t length = N;
    if (argc > 1)
    {
        length = atoi(argv[1]);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *data = (double *)malloc(length * sizeof(double));
    double *checksum = (double *)malloc(length * sizeof(double));

    if (0 == data)
    {
        printf("memory allocation failed");
        return 0;
    }

    srand(0);
    for (int i = 0; i < length; i++)
    {
        checksum[i] = (double)rand() / (double)RAND_MAX;
    }
    if (rank == 0)
    {
        for (int i = 0; i < length; i++)
        {
            data[i] = checksum[i];
        }
    }

    // Start time measuring
    //trivial implementation
    double aver_time = 0;
    auto start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        my_bcast(data, length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!(check(data, checksum)))
        {
            printf("FUNCTION NOT WELL IMPLEMENTED\n");
        }
    }
    auto end = steady_clock::now();
    duration<double> elapse_time = end - start;
    aver_time = elapse_time.count() / LOOP;
    if (rank == 0)
    {
        std::cout << "TRIVIAL "
                  << "Data Size: " << length << ",time spent: " << aver_time << std::endl;
    }

    //TREE
    start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        my_bcast_tree(data, length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!(check(data, checksum)))
        {
            printf("FUNCTION NOT WELL IMPLEMENTED\n");
        }
    }
    end = steady_clock::now();
    elapse_time = end - start;
    aver_time = elapse_time.count() / LOOP;
    if (rank == 0)
    {
        std::cout << "TREE "
                  << "Data Size: " << length << ",time spent: " << aver_time << std::endl;
    }

    // Buffer
   start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        broadcast_buffer(data, length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!(check(data, checksum)))
        {
            printf("FUNCTION Buffer NOT WELL IMPLEMENTED\n");
        }
    }
    end = steady_clock::now();
    elapse_time = end - start;
    aver_time = elapse_time.count() / LOOP;
    if (rank == 0)
    {
        std::cout << "Buffer: "
                  << "Data Size: " << length << ",time spent: " << aver_time << std::endl;
    }

    // Bcast
    start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        MPI_Bcast(data, length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!(check(data, checksum)))
        {
            printf("FUNCTION NOT WELL IMPLEMENTED\n");
        }
    }
    end = steady_clock::now();
    elapse_time = end - start;
    aver_time = elapse_time.count() / LOOP;
    if (rank == 0)
    {
        std::cout << "MPI_Bcast "
                  << "Data Size: " << length << ",time spent: " << aver_time << std::endl;
    }
    free(data);
    free(checksum);
    MPI_Finalize();
}