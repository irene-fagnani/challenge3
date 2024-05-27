#include "utils.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
using json=nlohmann::json;

int main(int argc, char** argv){
int provided;
MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE,&provided);
int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
std::ifstream file("data.json");
json data = json::parse(file);

MuparserFun f(data.value("f",""));
int n=data.value("n",10);
int niter=data.value("niter",1000);
double tol=data.value("tol",1e-4);

if(n==1){
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);
}


double h=1/(n-1);

std::vector<int> recv_counts(size,0), recv_start_idx(size,0);   
int start_idx=0;
for(int i=0;i<size;i++){
    recv_counts[i]=(n%size>i) ? n/size+1:n/size;
    recv_start_idx[i]=start_idx;
    start_idx+=recv_counts[i];
}

int local_n=(n%size>rank) ? n/size+1:n/size;

// usare una forma per definire le matrici che sia contigua nella memoria
std::vector<std::vector<double>> local_U(recv_counts[rank], std::vector<double>(n, 0.0));
std::vector<std::vector<double>> local_U_old(recv_counts[rank], std::vector<double>(n, 0.0));

initialize_grid(local_U_old, 0.0);


std::vector<double> prev_row(n, 0.0);
std::vector<double> next_row(n, 0.0);

double global_error = tol+1;
int iter = 0;
while (global_error > tol && iter < niter) {
    // Exchange boundary rows with neighboring processes
    if (rank > 0){
        std::cout<<"\nInvio prima riga al processo precedente "<<rank<<std::endl;
        MPI_Send(local_U[0].data(), n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        for(int i=0;i<n;i++)
            std::cout<<local_U[0][i]<<" ";
            std::cout<<std::endl;
        //std::cout<<"\nRicevo prev_row dal processo precedente "<<rank<<std::endl;
        //MPI_Recv(prev_row.data(), n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
   
    if (rank < size - 1){
        //std::cout<<"\nInvio ultima riga al processo successivo "<<rank<<std::endl;
        //MPI_Send(local_U[local_n - 1].data(), n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        std::cout<<"\nRicevo next_row dal processo successivo "<<rank<<std::endl;
        MPI_Recv(next_row.data(), n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i=0;i<n;i++)
            std::cout<<next_row[i]<<" ";
            std::cout<<std::endl;
    }
    
    // Perform Jacobi iteration
    int start_i=0, end_i=local_n;

    if(rank==0)
        start_i=1;

    if(rank==size-1)
        end_i=local_n-1;

    local_U_old = local_U;
     for(std::size_t i=start_i;i<end_i;++i){
        for(std::size_t j=1;j<n-1;++j){
            local_U[i][j]=0.25*(local_U_old[i-1][j]+local_U_old[i+1][j]+local_U_old[i][j-1]+local_U_old[i][j+1]+h*h*f(i,j));
        }
    }

    // Compute local error
    double local_error = compute_error(local_U, local_U_old, local_n);
        
    // Compute global error
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global_error = sqrt(h*global_error);

     iter++;
    }

    // Gather results to the master process
    std::vector<std::vector<double>> U_global;
    if (rank == 0) {
        U_global.resize(n, std::vector<double>(n, 0.0));
    }

    for (int i = 1; i <= local_n; ++i) {
        MPI_Gather(U[i].data(), n, MPI_DOUBLE, U_global[i-1].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    // Save the solution to a file
    if (rank == 0) {
        generateVTKFile("solution.vtk",U_global, n, h);
    }



MPI_Finalize();
}