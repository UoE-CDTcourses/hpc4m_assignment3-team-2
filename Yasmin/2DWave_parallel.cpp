#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>

using namespace std;


int main(int argc, char* argv[]){

int rank;
int size =4;
//cout << size << endl;
//initialise MPI
MPI_Init(NULL,NULL);
MPI_Comm_rank( MPI_COMM_WORLD, &rank);
MPI_Comm_size( MPI_COMM_WORLD, &size);


int M = 2305;  // M length intervals.
double T = 1;  // final time.
double dt = 0.2/M;
int N = T/dt;
int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
double dy = 2./M;
double dy2 = dy*dy;
double dtdy2 = dt*dt/dy2;
int no = (M-1)/size+2;

// dynamically allocate enough memory to create the array.
double*** U = new double** [3]; // this creates an array of size 3 which is ready do hold 2D arrays in each entry.
for (int i = 0; i < 3; ++i) {   // this loop puts an array of size M+1 on each of the three entries.
  U[i] = new double*[M+1];		  // we get a 2D array.
  for (int j = 0; j < M+1; ++j){ // put another array on each entry of the 2D array to get a 3D array.
    U[i][j] = new double[M+1];
  }
}



//time stepping
MPI_Barrier(MPI_COMM_WORLD);
double t_comp1=MPI_Wtime();

for (int t=1; t<=N; ++t){
if (rank ==0 ){
if (t==1){
for (int i=1; i<no; ++i){
        for (int j=1; j<M; ++j){
                U[0][i][j] = U[1][i][j] = exp( -40 * ( (i*dy-1-0.4)*(i*dy-1-0.4) + (j*dy-1)*(j*dy-1) ) );
        }
}

for (int k =1 ; k<size; k++){
for (int j =1; j<no-1; j++){
	MPI_Recv(U[0][k*(no-2)+j+1],M+1,MPI_DOUBLE,k,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
}
for (int i=0; i<=M; ++i){
        U[0][0][i] = U[0][M][i] = U[1][0][i] = U[1][M][i] = U[2][0][i] = U[2][M][i] = 0;
        U[0][i][0] = U[0][i][M] = U[1][i][0] = U[1][i][M] = U[2][i][0] = U[2][i][M] = 0;
}

// print initial U values to file, row by row.
double t_sav1=MPI_Wtime();
ofstream out {"U_t0.csv"};
out<<fixed<<setprecision(4);
	for(int i=0; i<=M; ++i){
		for(int j=0; j<=M; ++j){					
			out<<U[0][i][j]<<" ";
		}
		out<<endl;
	}	
out.close();
double t_sav2=MPI_Wtime();
cout << "time for saving t=0 " << t_sav2- t_sav1 << endl;
}
else {
for (int j=1; j < no; j++){
for (int i = 1; i<M; i++){
           U[2][j][i] = 2*U[1][j][i] - U[0][j][i] + dtdy2*( U[1][j+1][i] + U[1][j-1][i] + U[1][j][i+1] + U[1][j][i-1] - 4*U[1][j][i] );
}
}
  MPI_Send(U[2][no-2], M+1, MPI_DOUBLE, rank+1, 0,MPI_COMM_WORLD);
  MPI_Recv(U[2][no-1],M+1, MPI_DOUBLE, rank+1 , 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);

// update the previous times.
      for (int j=1; j<no; ++j){
              for (int i=1; i<M; ++i){
                      U[0][j][i] = U[1][j][i];
                      U[1][j][i] = U[2][j][i];
                      }
      }

   // print out files for fixed times
      if(t==t1 || t==t2 || t==t3){
	for (int k =1 ; k<size; k++){
	for (int j =1; j<no-1; j++){
        MPI_Recv(U[2][k*(no-2)+j+1],M+1,MPI_DOUBLE,k,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	}
double t_sav1 = MPI_Wtime();
              stringstream ss;
              ss << fixed << setprecision(2) << t*dt; // this ensures that the double value gets converted
              string time = ss.str();                                          // to string with only 2 trailing digits.

              ofstream out {"U_t"+ss.str()+".csv"};
             out<<fixed<<setprecision(4);
                      for(int i=0; i<=M; ++i){
                              for(int j=0; j<=M; ++j){
                                      out<<U[2][i][j]<<" ";
                              }
                              out<<endl;
                      }
              out.close();
	      double t_sav2 = MPI_Wtime();
	      cout << "time to save time " << t << " : " << t_sav2-t_sav1 << endl;
      }
}
}

 else if (rank == (size-1)){
	if (t==1){

	 for (int j=1; j<no-1; ++j){
        for (int i=1; i<M; ++i){
                U[0][j][i]=U[1][j][i]=exp( -40 * ( ((rank*no+j)*dy-1-0.4)*((rank*no+j)*dy-1-0.4) + (i*dy-1)*(i*dy-1) ) );

        }
        MPI_Send(U[0][j],M+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD);
}
 }
else{
 for (int j=1; j < no-1; j++){
    for (int i = 1; i<M; i++){
    U[2][j][i]= 2*U[1][j][i]-U[0][j][i] + dtdy2*(U[1][j+1][i] + U[1][j-1][i] + U[1][j][i+1] + U[1][j][i-1] - 4*U[1][j][i])  ;
    }
    }
   MPI_Send(U[2][1],M+1,MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
   MPI_Recv(U[2][0],M+1,MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   //update old values
  for (int j = 0 ; j<no; j++){
          for (int i=0; i<=M; i++){
  U[0][j][i]=U[1][j][i];
  U[1][j][i]=U[2][j][i];
  }
  }
if(t==t1 || t==t2 || t==t3){
         for (int j=1; j<no-1; ++j){
		        MPI_Send(U[2][j],M+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD);
				
}
}
}
}

else{
if (t==1){
for (int j=1; j<no-1; ++j){
        for (int i=1; i<M; ++i){
		U[0][j][i]=U[1][j][i]=exp( -40 * ( ((rank*no+j)*dy-1-0.4)*((rank*no+j)*dy-1-0.4) + (i*dy-1)*(i*dy-1) ) );
        }
	MPI_Send(U[0][j],M+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD);
}
}
else{
	  for (int j=1; j < no; j++){
    for (int i = 1; i<M; i++){
    U[2][j][i]= 2*U[1][j][i]-U[0][j][i] + dtdy2*(U[1][j+1][i] + U[1][j-1][i] + U[1][j][i+1] + U[1][j][i-1] - 4*U[1][j][i])  ;
    }
    }
   MPI_Send(U[2][1],M+1,MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
   MPI_Send(U[2][no-2],M+1,MPI_DOUBLE, rank+1 , 0,MPI_COMM_WORLD);
   MPI_Recv(U[2][0],M+1,MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(U[2][no-1],M+1,MPI_DOUBLE, rank+1 , 1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  

   //update old values
  for (int j = 0 ; j<no; j++){
          for (int i=0; i<=M; i++){
  U[0][j][i]=U[1][j][i];
  U[1][j][i]=U[2][j][i];
  }
  }
if(t==t1 || t==t2 || t==t3){
         for (int j=1; j<no-1; ++j){
                        MPI_Send(U[2][j],M+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD);
}
}

}
}

}
//end for loop

// now we need to delete the array.
cout << "final barrier" << endl;
MPI_Barrier(MPI_COMM_WORLD);
double t_comp2= MPI_Wtime();
cout << "time for everything is " << t_comp2-t_comp1 << endl;
MPI_Finalize();

for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < M+1; j++){
        delete[] U[i][j];
	 }    
    delete[] U[i];
}

delete[] U;


return 0;

}
