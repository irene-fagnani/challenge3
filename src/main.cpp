#include "utils.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
using json=nlohmann::json;

int main(int argc, char** argv){

std::ifstream file("data.json");
json data = json::parse(file);

MuparserFun f(data.value("f",""));
int n=data.value("n",100);
int niter=data.value("niter",1000);
double tol=data.value("tol",1e-4);

MPI_Init(&argc, &argv);
int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);



int local_n=(n%size>rank)?n/size+1:n/size;

int start_row=(n%size<rank)?rank*local_n:rank*(local_n-1);

// Initialize the grid
std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
std::vector<std::vector<double>> U_old(n, std::vector<double>(n, 0.0));
initialize_grid(U, 1.0);

std::vector<std::vector<double>> local_U(local_n, std::vector<double>(n, 0.0));
std::vector<std::vector<double>> local_U_old(local_n, std::vector<double>(n, 0.0));

MPI_Scatter(U.data(), local_n*n, MPI_DOUBLE, local_U.data(), local_n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// Main Jacobi iteration loop
double global_error = tol + 1;
int iter = 0;

while (global_error > tol && iter < niter) {
    // Exchange boundary rows with neighboring processes
    if (rank > 0) {
        MPI_Send(U[1].data(), n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        MPI_Recv(U[0].data(), n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    if (rank < size-1) {
        MPI_Send(U[local_n].data(), n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        MPI_Recv(U[local_n+1].data(), n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Perform Jacobi iteration
    U_old = U;
    jacobi_iteration(U, U_old,local_n, n, f);

    // Compute local error
    double local_error = compute_error(U, U_old, local_n);
        
    // Compute global error
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global_error = sqrt(global_error);

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
    
    double h=1/(n-1);
    // Save the solution to a file
    if (rank == 0) {
        generateVTKFile("solution.vtk",U_global, n, h);
    }


MPI_Finalize();
}