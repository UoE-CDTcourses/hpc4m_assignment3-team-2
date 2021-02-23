//
//
//

#include <iostream>
#include <mpi.h>
using namespace std;

//#define N 100000 // Number of iterations
#define m 300000 // Size of array
int array[m] = {0}; // all elements 0 in C++

int main(){
    int rank,           // Task identifier
        size,           // Number of tasks
        tag,            // Interaction tag
        n;              //
    double start, finish, timed;

    /* Settings for MPI */
    MPI_Comm comm;
    MPI_Status status;
    comm = MPI_COMM_WORLD;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (size < 2){
        cout << "There is only one task.\n";
        cout << "We need at least one worker task!";
        MPI_Abort(comm, 1);
        return 1;
    }

    for (unsigned int k = 1; k <= 20; k++){
        int N = 10000 * k;
        //int m = 125000 + 100000*k;
        //int array[m] = {0};

        MPI_Barrier(comm);
        // Send first message
        if (rank == 0){
            cout << "MPI started with " << size << " tasks" << endl << endl;
            start = MPI_Wtime();     // Start measuring time
            MPI_Send(array, m, MPI_INT, 1, 0, comm);
        }


        long int i = 1;
        while ( i < N){
            /* Process 1 */
            if (rank == 1){
                MPI_Recv(array, m, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
                MPI_Send(array, m, MPI_INT, 0, 0, comm);
            }

            /* Process 0 */
            if (rank == 0){
                MPI_Recv(array, m, MPI_INT, 1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Send(array, m, MPI_INT, 1, 0, comm);
            }
            i += 2;
        }

        if (rank == 1){
            MPI_Recv(array, m, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
        }


        // Present finish time

        MPI_Barrier(comm);
        if (rank == 0){
            finish = MPI_Wtime();
            timed  = finish - start;
            cout << N << " messages ";
            printf("done in %f seconds.\n", timed);
            cout << "Short: " << endl;
            printf("%d %d %e %4.2f %e", m*4, N, timed,
                                (double)N/timed, (double)m*4.0/(timed * 2024));
            cout << endl;
        }

    }

    // Wrap up
    MPI_Finalize();

    return 0;
}
