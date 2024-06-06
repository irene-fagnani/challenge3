#include "utils.hpp"
/**
 * @brief Namespace containing the parallel implementation of the Jacobi method.
 * 
 */
namespace parallel{
    /**
 * @brief Initialize the grid with the boundary condition provided as input.
 * @tparam T value type of the entries of the grid.
 * @param u The grid to inizialize.
 * @param boundary_condition The value to assign to the boundary of the grid.
 */

template <typename T>
void set_boundary_conditions(std::vector<std::vector<T>> & u, T boundary_condition){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(std::size_t i=0;i<u.size();++i){
        u[i][0]=u[i][u.size()-1]=boundary_condition;
    }

    if(rank==0){
        for(std::size_t j=0;j<u[0].size();++j){
            u[0][j]=boundary_condition;
        }
    }

    if(rank==size-1){
        for(std::size_t j=0;j<u[u.size()-1].size();++j){
            u[u.size()-1][j]=boundary_condition;
        }
    }
}

/**
 * @brief Inizialize the vectors recv_counts and recv_start_idx.
 * @param recv_counts Vector to inizialize with the number of elements to receive from each process.
 * @param recv_start_idx Vector to inizialize with the starting index of the elements to receive from each process.
 * @param n Number of rows and columns of the grid.
 */
void initialize_recv_vectors(std::vector<int> & recv_counts, std::vector<int> & recv_start_idx, int n){
    int size=0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int start_idx=0;
    for(int i=0;i<recv_counts.size();++i){
        recv_counts[i]=(n%size>i) ? n/size+1:n/size;
        recv_start_idx[i]=start_idx;
        start_idx+=recv_counts[i];
    }
}

/**
 * @brief Exchange boundary rows with neighboring processes
 * 
 * @tparam T Type of the elements in the grid.
 * @param local_U Vector containing the grid of the current process.
 * @param prev_row Vector containing the previous row of the current process.
 * @param next_row Vector containing the next row of the current process.
 * @param n Number of columns in the grid.
 */
template<typename T>
void send_and_receive_neighbors(std::vector<std::vector<T>> & local_U, std::vector<double> & prev_row, std::vector<double> & next_row, int n, int local_n){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank > 0){
            MPI_Send(local_U[0].data(), n, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            MPI_Recv(prev_row.data(), n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
   
        if (rank < size - 1){
            MPI_Send(local_U[local_n - 1].data(), n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(next_row.data(), n, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
}
/**
 * @brief Compute the error at each iteration
 * @tparam T Type of the elements in the grid.
 * @param local_U Local grid of the current process and iteration.
 * @param local_U_old Local grid of the current process and of the last iteration.
 * @param h Lenght of each sub-interval.
 * @param local_n Number of rows of the current process.
 * @param n Number of columns of the grid
 * @return The error at each iteration
 */

template <typename T>
double compute_local_error(std::vector<T> local_U,std::vector<T> local_U_old,double h, int local_n,int n){
    double error=0;
    for(std::size_t i=0;i<local_n;++i){
        for(std::size_t j=0;j<n;++j){
        error+=std::abs(local_U[i][j]-local_U_old[i][j])*std::abs(local_U[i][j]-local_U_old[i][j]);
    }
 }
 return error;
        
}


/**
 * @brief Collect the solution from all the processes.
 * 
 * @param U Total grid containing the solution.
 * @param local_U Local grid of the current process.
 * @param recv_counts Vector containing the number of rows of each process.
 * @param recv_start_idx Vector containing the starting index of the rows of each process.
 * @param n Dimension of the total grid.
 */
void collect_solution(std::vector<std::vector<double>> & U, std::vector<std::vector<double>> & local_U, std::vector<int> & recv_counts, std::vector<int> & recv_start_idx, int n){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
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

}
/**
 * @brief Run a Jacobi method iteration.
 * @tparam T Type of the elements in the grid.
 * @param local_U Local grid of the current process and iteration.
 * @param local_U_old Local grid of the current process and of the last iteration.
 * @param f Force function.
 * @param h Lenght of each sub-interval.
 * @param recv_counts Number of rows of the current process.
 * @param recv_start_idx Starting index of the rows of the current process.
 * @param n Number of columns of the grid.
 */
template<typename T>
void run_jacobi(std::vector<std::vector<T>> & local_U,std::vector<std::vector<T>> & local_U_old,std::vector<T> prev_row,std::vector<T> next_row,MuparserFun f,double h,int local_n,int start_idx,int n){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    # pragma omp shared for(local_U)
    for(std::size_t i=1;i<local_n-1;++i){

        for(std::size_t j=1;j<n-1;++j){

            local_U[i][j]=0.25*(local_U_old[i-1][j]+local_U_old[i+1][j]+local_U_old[i][j-1]+local_U_old[i][j+1]+h*h*f(h*(start_idx+i),j*h));

            }
        }
    MPI_Barrier(MPI_COMM_WORLD);    
        
        
        
    if(rank<size-1){
        
        # pragma omp shared for(local_U)
        for(std::size_t j=1;j<n-1;++j){

            local_U[local_n-1][j]=0.25*(local_U_old[local_n-2][j]+next_row[j]+local_U_old[local_n-1][j-1]+local_U_old[local_n-1][j+1]+h*h*f(((local_n-1)+start_idx)*h,j*h));

        }
    }

    if(rank>0){
        
        # pragma omp shared for(local_U)
        for(std::size_t j=1;j<n-1;++j){

            local_U[0][j]=0.25*(prev_row[j]+local_U_old[0][j-1]+local_U_old[0][j+1]+local_U_old[1][j]+h*h*f(start_idx*h,j*h));
            
            }
        }
}

/**
 * @brief Solve the problem in parallel using the Jacobi method.
 * 
 * @param tol Tolerance of the solution.
 * @param n Dimension of the grid.
 * @param niter Maximum number of iterations.
 * @param f Force function.
 * @param h Length of each sub-interval.
 */
double solve(double tol, int n, int niter, MuparserFun f, MuparserFun u_exact, double h){
    double L2_error=0;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Initialize recv_counts, recv_start_idx vectors. 
    // recv_counts contains the number of rows each process will receive, 
    // recv_start_idx contains the starting index of each process.
    std::vector<int> recv_counts(size,0), recv_start_idx(size,0);   
    parallel::initialize_recv_vectors(recv_counts, recv_start_idx, n);
    
    // Initialize local_U and local_U_old vectors, with zero.
    std::vector<std::vector<double>> local_U(recv_counts[rank], std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> local_U_old(recv_counts[rank], std::vector<double>(n, 0.0));

    // Initialize the grid with boundary conditions
    parallel::set_boundary_conditions(local_U, 0.0);
    
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
    parallel::send_and_receive_neighbors(local_U, prev_row, next_row, n, recv_counts[rank]);

    // Perform Jacobi iteration  
    parallel::run_jacobi(local_U, local_U_old,prev_row,next_row, f, h, recv_counts[rank], recv_start_idx[rank], n);

     
    // Compute local error
    error=parallel::compute_local_error(local_U, local_U_old, h, recv_counts[rank], n);

    // Compute global error
    MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    global_error = sqrt(h*global_error);
    local_U_old = local_U; 
     
    iter++;
    
    }
   

    parallel::collect_solution(U, local_U, recv_counts, recv_start_idx, n);

   
    MPI_Barrier(MPI_COMM_WORLD);

    // print the solution
    if (rank == 0) { 
        generateVTKFile("output/ParallelSol.vtk",U, n, h);
        L2_error=compute_L2_error(U,u_exact,h);
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    return L2_error;

   
}
}