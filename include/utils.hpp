
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
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ostringstream oss;
    oss << filename << "_" << size << "_" << n << ".vtk";
    std::string completeFilename = oss.str();
    std::ofstream vtkFile(completeFilename);

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

/**
 * @brief Compute the L2 error between the numerical solution and the exact solution.
 * 
 * @tparam T Type of the elements in the grid.
 * @param U Computed solution.
 * @param U_exact Exact solution.
 * @param h Dimension of each sub-interval.
 * @return double L2 error between the computed and the exact solution.
 */
template<typename T>
double compute_L2_error(std::vector<std::vector<T>> & U,MuparserFun U_exact,double h){
    double error=0;
    for(std::size_t i=0;i<U.size();++i){
        for(std::size_t j=0;j<U.size();++j){
        error+=std::abs(U[i][j]-U_exact(i*h,j*h))*std::abs(U[i][j]-U_exact(i*h,j*h));
        }
    }
    return sqrt(h*error);
}

/**
 * @brief Save the data in a file.
 * 
 * @param cores Number of cores used.
 * @param matrix_size Dimension of the matrix.
 * @param error Value  to save in the file.
 * @param filename Name of the file to save the data.
 */
void save_data(int cores, int matrix_size, int matrix_size_max, auto value,const std::string & filename){

    std::ofstream outfile(filename, std::ios::app);
    if (!outfile) {
        std::cerr << "Error: could not open the file." << std::endl;
        exit(1);
    }
    
    outfile << cores << " " << matrix_size << " " << value << "\n";
    if(matrix_size==matrix_size_max){
        outfile << "\n";
    }
    


    outfile.close();
}
#endif


