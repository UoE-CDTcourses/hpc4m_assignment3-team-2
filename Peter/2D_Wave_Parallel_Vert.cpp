#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <string>

using namespace std;

int main(int argc, char* argv[]){


int rank, ierr, size;
MPI_Comm comm;
const int root = 0;
comm = MPI_COMM_WORLD;
MPI_Init(NULL,NULL);
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);

int M = 3601;  // M length intervals.
int J = ((M+1)-2)/size +2;
double T = 1;  // final time.
double dt = 0.2/M;
int N = T/dt;
int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
double dy = 2./M;
double dy2 = dy*dy;
double dtdy2 = dt*dt/dy2;

// double U[3][M+1][M+1];  // 3D array, the first dimension for three different times.
								  // this does not work for large M, however... test it!

double** UFinal = new double* [M+1]; // this creates a 1D array which is ready to hold 1D arrays in each entry.
for (int i = 0; i < M+1; ++i){ // put another array on each entry of the 1D array to get a 2D array.
      UFinal[i] = new double[M+1];
    }


// dynamically allocate enough memory to create the array.
double*** U = new double** [3]; // this creates an array of size 3 which is ready do hold 2D arrays in each entry.
for (int i = 0; i < 3; ++i) {   // this loop puts an array of size M+1 on each of the three entries.
  U[i] = new double*[M+1];		  // we get a 2D array.
  for (int j = 0; j < M+1; ++j){ // put another array on each entry of the 2D array to get a 3D array.
    U[i][j] = new double[J];
  }
}
// apart from array deletion at the end, everything else stays as usual..

// initialize numerical array with given conditions.
for (int i=1; i<M; ++i){
	for (int j=0; j<J; ++j){
		U[0][i][j] = U[1][i][j] = exp( -40 * ( (i*dy-1-0.4)*(i*dy-1-0.4) + ((j+rank*(J-2))*dy-1)*((j+rank*(J-2))*dy-1) ) );
	}
}

if (rank == size - 1){
  for (int i=0; i<=M; ++i){
   U[0][i][J-1] = U[1][i][J-1] = U[2][i][J-1] = 0;
  }
}
if (rank == root){
  for (int i=0; i<=M; ++i){
    U[0][i][0]  = U[1][i][0]  = U[2][i][0] = 0;
  }

}
for (int i=0; i<J; ++i){
  U[0][0][i] = U[0][M][i] = U[1][0][i] = U[1][M][i] = U[2][0][i] = U[2][M][i] = 0;
}


// print initial U values to file, row by row.
//Gathering initial U
if (rank == root){
for (int i=0;i<M+1;i++){
  for (int j=0;j<J;j++){
UFinal[i][j] = U[0][i][j];
}}}

for (int i=0;i<J;i++){
  for (int j=0; j<=M;j++){
    if (rank!=root){

  MPI_Send(&U[0][j][i],1,MPI_DOUBLE,root,rank,MPI_COMM_WORLD);
}
  if (rank == root){

  for (int k=1; k<size; k++){
  MPI_Recv(&UFinal[j][(J-2)*k+i],1,MPI_DOUBLE,k,k,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }}

}}


//printing initial U
if (rank == root){
ofstream out {"U_t0.csv"};
out<<fixed<<setprecision(4);
	for(int i=0; i<=M; ++i){
		for(int j=0; j<=M; ++j){
			out<<UFinal[i][j]<<" ";
		}
		out<<endl;
	}
out.close();}

// use numerical scheme to obtain the future values of U.
for (int t=1; t<=N; ++t){
	for (int i=1; i<M; ++i){
		for (int j=1; j<J-1; ++j){
			U[2][i][j] = 2*U[1][i][j] - U[0][i][j]
						 	 + dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j]);
		}

	}

//Halo swapping
for (int k = 1; k<M;k++){
  if (rank < size-1){
	MPI_Send(&U[2][k][J-2],1,MPI_DOUBLE,rank+1,2,comm);
}
	if (rank > 0){
	MPI_Send(&U[2][k][1],1,MPI_DOUBLE,rank-1,2,comm);
	MPI_Recv(&U[2][k][0],1,MPI_DOUBLE,rank-1,2,comm,MPI_STATUS_IGNORE);
}
	if (rank < size-1){
	MPI_Recv(&U[2][k][J-1],1,MPI_DOUBLE,rank+1,2,comm,MPI_STATUS_IGNORE);
}
}
// update the previous times.
for (int i=1; i<M; ++i){
  for (int j=0; j<J; ++j){
    U[0][i][j] = U[1][i][j];
    U[1][i][j] = U[2][i][j];
    }
}
   // print out files for fixed times

	if(t==t1 || t==t2 || t==t3){


    //Gathering




    if (rank == root){
    for (int i=0;i<M+1;i++){
      for (int j=0;j<J;j++){
    UFinal[i][j] = U[2][i][j];
    }}}

    int MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<J;i++){
      for (int j=0; j<=M;j++){
        if (rank!=root){

      MPI_Send(&U[2][j][i],1,MPI_DOUBLE,root,rank,MPI_COMM_WORLD);
    }
      if (rank == root){

      for (int k=1; k<size; k++){
      MPI_Recv(&UFinal[j][(J-2)*k+i],1,MPI_DOUBLE,k,k,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }}

    }

  }


    if (rank == root){
		stringstream ss;
		ss << fixed << setprecision(2) << t*dt; // this ensures that the double value gets converted
		string time = ss.str();						 // to string with only 2 trailing digits.

		ofstream out {"U_t"+ss.str()+".csv"};
		out<<fixed<<setprecision(4);
			for(int i=0; i<=M; ++i){
				for(int j=0; j<=M; ++j){
					out<<UFinal[i][j]<<" ";
				}
				out<<endl;
			}
		out.close();}
  }
  if (rank == root){
	cout<<"iteration "<<t<<" done"<<endl;
}

}


// now we need to delete the arrays.
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < M+1; j++){
        delete[] U[i][j];
	 }
    delete[] U[i];
}

delete[] U;

for (int i = 0; i < M+1; i++)
{
    delete[] UFinal[i];
}

delete[] UFinal;

//return 0;
MPI_Finalize();
}
