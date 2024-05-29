
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


