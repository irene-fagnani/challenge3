
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <mpi.h>
#include "muparser_fun.hpp"

#ifndef UTILS_HPP
#define UTILS_HPP
/**
 * @brief Function to print a vector.
 * @tparam T Type of the elements in the vector.
 * @param v Vector to print.
 */

template <typename T>
void print_vector( std::vector<T> & v){
    for(std::size_t i=0;i<v.size();++i){
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}


/**
 * @brief Print a matrix.
 * @tparam T Type of the elements in the matrix.
 * @param matrix Matrix to print.
 */

template <typename T>
void print_matrix( std::vector<std::vector<T>> & matrix){
    for(std::size_t i=0;i<matrix.size();++i){
        print_vector(matrix[i]);
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

/**
 * @brief Initialize the grid with the boundary condition provided as input.
 * @tparam T value type of the entries of the grid.
 * @param u The grid to inizialize.
 * @param n The number of rows and cols of the grid.
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



void collect_solution(std::vector<std::vector<double>> & U, std::vector<std::vector<double>> & local_U, std::vector<int> & recv_counts, std::vector<int> & recv_start_idx, int n){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
     if(rank == 0){
        for(int p=1;p<size;p++){
            for(std::size_t i=0;i<recv_counts[p];++i){
                MPI_Recv(U[recv_start_idx[p]+i].data(), n, MPI_DOUBLE, p, p+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::cout<<"ricevuto"<<std::endl;
            }
        }

        for(std::size_t i=0;i<recv_counts[rank];++i){
            U[i]=local_U[i];
        }

    }else{
        for(std::size_t i=0;i<recv_counts[rank];++i){
            MPI_Send(local_U[i].data(), n, MPI_DOUBLE, 0, rank+i, MPI_COMM_WORLD);
            std::cout<<"inviato"<<std::endl;
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
 * @brief Generate a vtk structured grid file with the scalar field provided as input.
 * @tparam T Elements type in the grid.
 * @param filename Name of the file to generate.
 * @param scalarField The square scalar field to write in the file.
 * @param n Number of rows and columns of the scalar field.
 * @param h Spacing between points in the scalar field, equal to 1/(n-1).
 */
template <typename T>
void generateVTKFile(const std::string & filename, 
                     const std::vector<std::vector<T>> & scalarField, 
                     int n, double h) {

    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }
    
    vtkFile <<  "# vtk DataFile Version 3.0\n";
    vtkFile << "Scalar Field Data\n";
    vtkFile << "ASCII\n";                                
    

  
    vtkFile << "DATASET STRUCTURED_POINTS\n";                             
    vtkFile << "DIMENSIONS " << n << " " << n << " " << 1 << "\n";  
    vtkFile << "ORIGIN 0 0 0\n";                                         
    vtkFile << "SPACING" << " " << h << " " << h << " " << 1 << "\n";   
    vtkFile << "POINT_DATA " << (n) * (n) << "\n";                  
                                                                
    
    
    vtkFile << "SCALARS scalars double\n";               
    vtkFile << "LOOKUP_TABLE default\n";                 

    
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            vtkFile <<  scalarField[i][j] << "\n";
        }
    }

    return;
}

#endif


