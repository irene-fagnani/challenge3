#include "utils_parallel.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <chrono>
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
    int n;
    bool test=false;
    int nmax=0;
    
    if(argc>1){
        n=std::stoi(argv[1]);
    }else{
        n=data.value("n",11);
    }

    if(argc>2){
        test=true;
        nmax=std::stoi(argv[3]);
    }

    int niter=data.value("niter",1000);
    double tol=data.value("tol",1e-4);

    if(n==1){
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);
    }
    // compute the step size
    double h=1.0/(n-1);
    
    long double L2_error=0;
    
    double time=0;

    auto t0=std::chrono::high_resolution_clock::now();
    L2_error=parallel::solve(tol,n,niter,f,u_exact,h);
    auto t1=std::chrono::high_resolution_clock::now();
    auto delta_t=(std::chrono::duration_cast<std::chrono::microseconds>(t1-t0));
    
    if(rank==0){
        std::cout<<"Time for the computation, in the parallel case: "<<delta_t.count()<<" microseconds\n";
        time=delta_t.count();
        std::cout<<"L2 error: "<<L2_error<<std::endl;
        if(test==true){
            save_data(size,n,nmax,L2_error,"test/data/errors.dat");
            save_data(size,n,nmax,time,"test/data/times.dat");
        }
    }
    
    

    MPI_Finalize();
    return 0;
}

