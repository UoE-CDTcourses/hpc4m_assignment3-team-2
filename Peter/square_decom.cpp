#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]){


int rank, ierr, size;
MPI_Comm comm;
const int root = 0;
comm = MPI_COMM_WORLD;
MPI_Init(NULL,NULL);
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);
int X = sqrt(size);

int M = 2307;  // M length intervals.
int J = ((M+1)-2)/X +2; //change this for squares
double T = 1;  // final time.
double dt = 0.2/M;
int N = T/dt;
int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
double dy = 2./M;
double dy2 = dy*dy;
double dtdy2 = dt*dt/dy2;


// now lets create our annoyingly large array...
// double U[3][M+1][M+1];  // 3D array, the first dimension for three different times.
								  // this does not work for large M, however... test it!

double** UFinal = new double* [M+1]; // this creates an 1D array which is ready to hold 1D arrays in each entry.
for (int i = 0; i < M+1; ++i){ // put another array on each entry of the 1D array to get a 2D array.
      UFinal[i] = new double[M+1];
    }


// dynamically allocate enough memory to create the array.
double*** U = new double** [3]; // this creates an array of size 3 which is ready do hold 2D arrays in each entry.
for (int i = 0; i < 3; ++i) {   // this loop puts an array of size M+1 on each of the three entries.
  U[i] = new double*[J];		  // we get a 2D array.
  for (int j = 0; j < J; ++j){ // put another array on each entry of the 2D array to get a 3D array.
    U[i][j] = new double[J];
  }
}

//coordinate system for squares
int column = rank%X; //X is the square root of the number of boxes in each direction.
int row = (rank-column)/X;



// initialize numerical array with given conditions.
for (int i=0; i<J; ++i){
	for (int j=0; j<J; ++j){
		U[0][i][j] = U[1][i][j] = exp( -40 * ( ((i+row*(J-2))*dy-1-0.4)*((i+row*(J-2))*dy-1-0.4) + ((j+column*(J-2))*dy-1)*((j+column*(J-2))*dy-1) ) );
	}
}
if (column == X-1){
  for (int i=0; i<J; ++i){
   U[0][i][J-1] = U[1][i][J-1] = U[2][i][J-1] = 0;
  }
}
if (column == 0){
  for (int i=0; i<J; ++i){
    U[0][i][0]  = U[1][i][0]  = U[2][i][0] = 0;
  }

}
if (row == X-1){
  for (int i=0; i<J; ++i){
  U[0][J-1][i] = U[1][J-1][i] = U[2][J-1][i] = 0;

}
}

if (row == 0){
  for (int i=0; i<J; ++i){
  U[0][0][i] = U[1][0][i] = U[2][0][i] = 0;
}

}
if (rank == root){
for (int i=0;i<J;i++){
  for (int j=0;j<J;j++){
UFinal[i][j] = U[0][i][j];
}}}

for (int i=0;i<J;i++){
  for (int j=0; j<J;j++){
    if (rank!=root){

  MPI_Send(&U[0][i][j],1,MPI_DOUBLE,root,rank,MPI_COMM_WORLD);
}
  if (rank == root){

  for (int k=1; k<size; k++){
  int K_col = k%X; //column
  int K_row = (k - K_col)/X; //row
  MPI_Recv(&UFinal[(J-2)*K_row+i][(J-2)*K_col+j],1,MPI_DOUBLE,k,k,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	for (int i=1; i<J-1; ++i){
		for (int j=1; j<J-1; ++j){
			U[2][i][j] = 2*U[1][i][j] - U[0][i][j]
						 	 + dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j]);
		}

	}


  //Halo swapping rows columns and corners diagonally.
  if (row<X-1){
    MPI_Send(&U[2][J-2][1],J-2,MPI_DOUBLE,column+X*(row+1),2,comm);


  }

  if (row>0){
    MPI_Send(&U[2][1][1],J-2,MPI_DOUBLE,column+X*(row-1),2,comm);
    MPI_Recv(&U[2][0][1],J-2,MPI_DOUBLE,column+X*(row-1),2,comm,MPI_STATUS_IGNORE);



  }


  if (row < X-1){

    MPI_Recv(&U[2][J-1][1],J-2,MPI_DOUBLE,column+X*(row+1),2,comm,MPI_STATUS_IGNORE);


  }
  /////////////////////////////////////////
  for (int k = 1; k<J-1;k++){

  if (column < X-1){
    MPI_Send(&U[2][J-2][k],1,MPI_DOUBLE,rank+1,1,comm);
    }

  if (column >0){
    MPI_Send(&U[2][1][k],1,MPI_DOUBLE,rank-1,1,comm);
    MPI_Recv(&U[2][0][k],1,MPI_DOUBLE,rank-1,1,comm,MPI_STATUS_IGNORE);

    }

  if (column <X-1){
    MPI_Recv(&U[2][J-1][k],1,MPI_DOUBLE,rank+1,1,comm,MPI_STATUS_IGNORE);


    }
  }
  ////Corners/////////
  if (column>0){
    if (row>0){
      MPI_Send(&U[2][1][1],1,MPI_DOUBLE,column-1 + X*(row-1),3,comm);
    }
  }
  if (column < X-1){
    if (row< X-1){
      MPI_Recv(&U[2][J-1][J-1],1,MPI_DOUBLE,column+1 + X*(row+1),3,comm,MPI_STATUS_IGNORE);
      MPI_Send(&U[2][J-2][J-2],1,MPI_DOUBLE,column+1 + X*(row+1),3,comm);
    }
  }
  if (column >0){
    if (row>0){

      MPI_Recv(&U[2][0][0],1,MPI_DOUBLE,column-1 + X*(row-1),3,comm,MPI_STATUS_IGNORE);

    }

  }
  if (column>0){
    if (row < X-1){
      MPI_Send(&U[2][1][J-2],1,MPI_DOUBLE,column-1 + X*(row+1),4,comm);


    }


  }
  if (column < X-1){
    if (row >0){
      MPI_Recv(&U[2][J-1][0],1,MPI_DOUBLE,column+1 + X*(row-1),4,comm,MPI_STATUS_IGNORE);

    }

  }


  if (column<X-1){
    if (row > 0){
      MPI_Send(&U[2][J-2][1],1,MPI_DOUBLE,column+1 + X*(row-1),4,comm);


    }


  }
  if (column > 0){
    if (row <X-1){
      MPI_Recv(&U[2][0][J-1],1,MPI_DOUBLE,column-1 + X*(row+1),4,comm,MPI_STATUS_IGNORE);

    }

  }


// update the previous times.
for (int i=0; i<J; ++i){
  for (int j=0; j<J; ++j){
    U[0][i][j] = U[1][i][j];
    U[1][i][j] = U[2][i][j];
    }
}




   // print out files for fcolumned times




	if(t==t1 || t==t2 || t==t3){



    if (rank == root){
    for (int i=0;i<J;i++){
      for (int j=0;j<J;j++){
    UFinal[i][j] = U[0][i][j];
    }}}

    for (int i=0;i<J;i++){
      for (int j=0; j<J;j++){
        if (rank!=root){

      MPI_Send(&U[2][i][j],1,MPI_DOUBLE,root,rank,MPI_COMM_WORLD);
    }
      if (rank == root){

      for (int k=1; k<size; k++){
      int K_col = k%X;
      int K_row = (k - K_col)/X;
      MPI_Recv(&UFinal[(J-2)*K_row+i][(J-2)*K_col+j],1,MPI_DOUBLE,k,k,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }}

    }}

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
    for (int j = 0; j < J; j++){
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
