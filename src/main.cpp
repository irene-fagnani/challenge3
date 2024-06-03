#include "utils_parallel.hpp"
#include "utils_sequential.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
using json=nlohmann::json;


int main(int argc, char** argv){

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // read data from json file
    std::ifstream file("data.json");
    json data = json::parse(file);
    MuparserFun f(data.value("f",""));
    MuparserFun u_exact(data.value("u_exact",""));
    int n=data.value("n",11);
    int niter=data.value("niter",1000);
    double tol=data.value("tol",1e-4);
    
    // compute the step size
    double h=1.0/(n-1);

    if(n==1){
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);
    }
    
    long double L2_error=0;
    
    if(size==1){
        // Solve the problem in sequential, and print the computation timing
    auto t0_s=std::chrono::high_resolution_clock::now();
    L2_error=sequential::solve(tol,n,niter,f,u_exact,h);
    auto t1_s=std::chrono::high_resolution_clock::now();
    auto delta_t_s=std::chrono::duration_cast<std::chrono::microseconds>(t1_s-t0_s);
    std::cout<<"Time for the multiplication, in the sequential case: "<<delta_t_s.count()<<" microseconds\n";
    }
    
    if(size>1){
      // Solve the problem in parallel, and print the computation timing
    auto t0_p=std::chrono::high_resolution_clock::now();
    L2_error=parallel::solve(tol,n,niter,f,u_exact,h);
    auto t1_p=std::chrono::high_resolution_clock::now();
    auto delta_t_p=std::chrono::duration_cast<std::chrono::microseconds>(t1_p-t0_p);
        if(rank==0){
            std::cout<<"Time for the computation, in the parallel case: "<<delta_t_p.count()<<" microseconds\n";
        }
    }
    
    if(rank==0){
        std::cout<<"L2 error: "<<L2_error<<std::endl;
    }
    
    

    MPI_Finalize();
    return 0;
}

