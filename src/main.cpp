#include "utils.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
using json=nlohmann::json;


int main(int argc, char** argv){

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // read data from json file
    std::ifstream file("data.json");
    json data = json::parse(file);
    MuparserFun f(data.value("f","x*y"));
    int n=data.value("n",11);
    int niter=data.value("niter",1000);
    double tol=data.value("tol",1e-4);
    
    // compute the step size
    double h=1.0/(n-1);

    if(n==1){
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);
    }

    // Initialize recv_counts, recv_start_idx vectors. 
    // recv_counts contains the number of rows each process will receive, 
    // recv_start_idx contains the starting index of each process.
    std::vector<int> recv_counts(size,0), recv_start_idx(size,0);   
    initialize_recv_vectors(recv_counts, recv_start_idx, n);
    
    // Initialize local_U and local_U_old vectors, with zero.
    std::vector<std::vector<double>> local_U(recv_counts[rank], std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> local_U_old(recv_counts[rank], std::vector<double>(n, 0.0));

    // Initialize the grid with boundary conditions
    set_boundary_conditions(local_U, 0.0);
    
    // Initialize the total grid, which will contains the final solution.
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    
    // Initialize the previous and next row vectors, which will contain the boundary rows of the neighboring processes.
    std::vector<double> prev_row(n, 0.0);
    std::vector<double> next_row(n, 0.0);
    
    // Initialize error, global error and the number of iterations.
    double global_error = tol+1;
    double error=0;
    int iter = 0;

    while (global_error > tol && iter < niter){

        // Exchange boundary rows with neighboring processes
        if (rank > 0){
            MPI_Send(local_U[0].data(), n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            MPI_Recv(prev_row.data(), n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
   
        if (rank < size - 1){
            MPI_Send(local_U[recv_counts[rank] - 1].data(), n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(next_row.data(), n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      
    // Perform Jacobi iteration  
        for(std::size_t i=1;i<recv_counts[rank]-1;++i){

            for(std::size_t j=1;j<n-1;++j){
                local_U[i][j]=0.25*(local_U_old[i-1][j]+local_U_old[i+1][j]+local_U_old[i][j-1]+
                local_U_old[i][j+1]+h*h*f(h*(recv_start_idx[rank]+i),j*h));
            }
        }
        
        
        
        if(rank<size-1){
            for(std::size_t j=1;j<n-1;++j){
                local_U[recv_counts[rank]-1][j]=0.25*(local_U_old[recv_counts[rank]-2][j]+next_row[j]+local_U_old[recv_counts[rank]-1][j-1]+local_U_old[recv_counts[rank]-1][j+1]+h*h*f(((recv_counts[rank]-1)+recv_start_idx[rank])*h,j*h));
            }
        }

        if(rank>0){
           for(std::size_t j=1;j<n-1;++j){
                local_U[0][j]=0.25*(prev_row[j]+local_U_old[0][j-1]+local_U_old[0][j+1]+local_U_old[1][j]+h*h*f(recv_start_idx[rank]*h,j*h));
            }
        }

     
        // Compute local error
        error=0;
         for(std::size_t i=0;i<recv_counts[rank];++i){
           for(std::size_t j=0;j<n;++j){
            error+=std::abs(local_U[i][j]-local_U_old[i][j])*std::abs(local_U[i][j]-local_U_old[i][j]);
        }
        }
        
        // Compute global error
         MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
         global_error = sqrt(h*global_error);
         local_U_old = local_U; 
     
        iter++;
    }
   


    if(rank == 0){
        for(int p=1;p<size;p++){
            for(std::size_t i=0;i<recv_counts[p];++i){
                MPI_Recv(U[recv_start_idx[p]+i].data(), n, MPI_DOUBLE, p, p+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for(std::size_t i=0;i<recv_counts[rank];++i){
            U[i]=local_U[i];
        }
    }else{
        for(std::size_t i=0;i<recv_counts[rank];++i){
            MPI_Send(local_U[i].data(), n, MPI_DOUBLE, 0, rank+i, MPI_COMM_WORLD);
        }
    
    }

   
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { 
        print_matrix(U);
        generateVTKFile("solution.vtk",U, n, h);
    }
     
    MPI_Barrier(MPI_COMM_WORLD);



MPI_Finalize();
return 0;
}

