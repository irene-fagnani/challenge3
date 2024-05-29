
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
 * @brief Compute the error at each iteration.
 * @tparam T Elements type in the grid.
 * @param u The grid at the current iteration.
 * @param u_old The grid at the previous iteration.
 * @param n Number of rows and columns of the grid.
 * @return double The error computed as the L2 norm of the difference between the two grids.
 */

template <typename T>
double compute_error(const std::vector<std::vector<T>> & u, const std::vector<std::vector<T>> & u_old, int n){

    if(n==1)
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);

    double error=0.0;

    #pragma omp parallel for reduction(+:error)
    for(std::size_t i=1;i<u.size();++i){
        for(std::size_t j=1;j<n;++j){
            error+=std::abs(u[i][j]-u_old[i][j])*std::abs(u[i][j]-u_old[i][j]);
        }
    }

    return error;
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


