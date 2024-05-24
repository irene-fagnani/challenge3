#include <vector>
#include <iostream>

/**
 * @brief Initialize the grid with the boundary condition provided as input.
 * @tparam T value type of the entries of the grid.
 * @param u The grid to inizialize.
 * @param n The number of rows and cols of the grid.
 * @param boundary_condition The value to assign to the boundary of the grid.
 */

template <typename T>
void intialize_grid(std::vector<std::vector<T>> & u, int n, T boundary_condition){
    for(std::size_t i=0;i<n;++i){
        u[i][0]=u[i][n-1]=u[0][i]=u[n-1][i]=boundary_condition;
    }
}

template <typename T>
double compute_error(const std::vector<std::vector<T>> & u, const std::vector<std::vector<T>> & u_old, int n){
    if(n==1){
        std::cerr<<"Error: n must be greater than 1"<<std::endl;
        exit(1);
    }
    double h=1/(n-1);
    double error=0.0;

    #pragma omp parallel for reduction(+:error)
    for(std::size_t i=1;i<u.size()-1;++i){
        for(std::size_t j=1;j<u.size()-1;++j){
            error+=std::abs(u[i][j]-u_old[i][j])*std::abs(u[i][j]-u_old[i][j]);
        }
    }
    return sqrt(error);
}




