#include "utils.hpp"

namespace sequential{
    /**
 * @brief Initialize the grid with the boundary condition provided as input.
 * @tparam T value type of the entries of the grid.
 * @param u The grid to inizialize.
 * @param boundary_condition The value to assign to the boundary of the grid.
 */

template <typename T>
void set_boundary_conditions(std::vector<std::vector<T>> & u, T boundary_condition){
    
        for(std::size_t i=0;i<u.size();++i){
            u[i][0]=u[i][u.size()-1]=boundary_condition;
        }
    
        for(std::size_t j=1;j<u[0].size()-1;++j){
            u[0][j]=u[u.size()-1][j]=boundary_condition;
        }
}



/**
 * @brief Compute the total error.
 * @tparam T Type of the elements in the grid.
 * @param U Local grid of the current iteration.
 * @param U_old Local grid of the last iteration.
 * @param h Lenght of each sub-interval.
 * @return The error.
 */

template <typename T>
double compute_error(std::vector<T> U,std::vector<T> U_old,double h){
    double error=0;
    for(std::size_t i=0;i<U.size();++i){
        for(std::size_t j=0;j<U.size();++j){
        error+=std::abs(U[i][j]-U_old[i][j])*std::abs(U[i][j]-U_old[i][j]);
        }
    }
    return sqrt(error);
        
}


/**
 * @brief Run a Jacobi method iteration.
 * @tparam T Type of the elements in the grid.
 * @param U Grid of the current  iteration.
 * @param U_old Grid of the last iteration.
 * @param f Force function.
 * @param h Lenght of each sub-interval.
 */

template<typename T>
void run_jacobi(std::vector<std::vector<T>> & U,std::vector<std::vector<T>> & U_old,MuparserFun f,double h){
  
    for(std::size_t i=1;i<U.size()-1;++i){

        for(std::size_t j=1;j<U.size()-1;++j){

            U[i][j]=0.25*(U_old[i-1][j]+U_old[i+1][j]+U_old[i][j-1]+U_old[i][j+1]+h*h*f(i*h,j*h));

            }
        }

    }


void solve(double tol, int n, int niter, MuparserFun f, double h){

// Initialize local_U and local_U_old vectors, with zero.
std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
std::vector<std::vector<double>> U_old(n, std::vector<double>(n, 0.0));

// Initialize the grid with boundary conditions
sequential::set_boundary_conditions(U, 0.0);

    
// Initialize error, global error and the number of iterations.
double error = tol+1;
int iter = 0;

while (error > tol && iter < niter){


// Perform Jacobi iteration  
sequential::run_jacobi(U, U_old, f, h);

     
// Compute local error
error=sequential::compute_error(U, U_old, h);

U_old = U; 
     
iter++;
    
}
   
generateVTKFile("output/SequentialSol.vtk",U, n, h);
   

}
}