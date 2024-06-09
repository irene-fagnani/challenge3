# Challenge 3: A matrix–free parallel solver for the Laplace equation

### Introduction

Consider the Laplace equation for modeling heat diffusion over a square domain with a prescribed temperature (Dirichlet conditions) on the entire boundary. 

A possible approach to solve this problem is the Jacobi iteration method. Given a uniform Cartesian grid of points along each coordinate direction, the goal is to find the discrete solution $u_{ij}=u(x_i,y_j), \quad i,j=1,...,n$ at each point of this grid.

We aim to represent the solution as a dense matrix U of size n x n.

 The matrix is initialized with zeroes, except for the first and last rows and columns, which contain the boundary condition values defined in the equation above.

The algorithm consists of the following iterative procedure:

1. For k=1,... until convergence:
+ Update each internal entry as the average of the values of a four-point stencil:

$U^{(k+1)}(i,j)=\frac{1}{4}(U^{(k)}(i-1,j)+U^{(k)}(i+1,j)+U^{(k)}(i,j-1)+U^{(k)}(i,j+1))$ for all i,j=2,...,n-1
 
​2. Compute the convergence criterion as the norm of the increment between $U^{(k+1)}$ and $U^{(k)}$ as follows: $e=\sqrt{h \cdot \sum_{i,j}(U^{(k+1)}(i,j)-U^{(k)}(i,j))^2}$. If e is smaller than the fixed tolerance, stop. h represents the mesh spacing.

## Features

+ The user can provide the data needed to solve the problem, by modify the `data.json` file. In particular, the user can provide the force function f and the exact solution u_exact by modify the correspondent fields in the `data.json` file, thanks to the muParser and json library.
+ Thanks to the chrono library, it is possible to compute the time of computing of the problem.
+ It is possible to run a test, that involves different nu,bers of core and different dimensions of the mesh grid. Read the `README.md` file in the `test` folder for further details.
+ It is possible to both run sequentially and in parallel the program and test both cases.


## Content of the repository
+ [test](https://github.com/irene-fagnani/challenge3/tree/main/test) folder, where it is possible to run a test (also, my results in tems of graphs are collected).
+ [include](https://github.com/irene-fagnani/challenge3/tree/main/include) folder, which contains all the `.hpp` files. [muparser_fun.hpp](https://github.com/irene-fagnani/challenge3/blob/main/include/muparser_fun.hpp), the file that allows to read a function as a string and transform it in a `MuParser` type and allows to use the `MuParser.h` library. [utils_parallel.hpp](https://github.com/irene-fagnani/challenge3/blob/main/include/utils_parallel.hpp), which contains all the functions needed to solve the problem. [utils.hpp](https://github.com/irene-fagnani/challenge3/blob/main/include/utils.hpp), which contains functions needed to save datas in files and print structures.
+ [src](https://github.com/irene-fagnani/challenge3/blob/main/src/), the folder that contains `.cpp` files (in this case, only the main.cpp).
+ [output](https://github.com/irene-fagnani/challenge3/blob/main/output/) folder, which contains the `.vtk` files create by the execution of the program (to create if there isn't).
+ [data.json](https://github.com/irene-fagnani/challenge3/blob/main/data.json), can be modified by the user to pass the parameters needed by the problem (force function, the exact solution, the dimension of the grid, the maximum number of iterations allowed and the tolerance).
+ [Doxyfile](https://github.com/irene-fagnani/challenge3/blob/main/Doxyfile), that allows to produce the documentation, thanks to the DoxyComments.
+ [hw.info](https://github.com/irene-fagnani/challenge3/blob/main/hw.info), that contains the informations about processors in my computer.
+ [Makefile](https://github.com/irene-fagnani/challenge3/blob/main/Makefile), that allows the compilation of the code, cleaning of objects and executables and production of documentation.

## Run the code
To clone the project:

```bash
git clone git@github.com:irene-fagnani/challenge3.git
```
To compile the code:
```bash
make 
```
To delete the object and executables:
```bash
make clean
``` 
To produce the documentation:
```bash
make doc
```
**NB**: make sure to have your `PACS_ROOT` defined as the folder where `pacs-examples/Examples` resides. This repository use files in that path. Moreover, also the files in Examples/lib are used, so make sure also to have `LD_LIBRARY_PATH` set as the path in which you have `pacs-examples/Examples/lib`

## Visualize the outputs
To visualize the `.vtk` files, open Paraview and select the path where the files arwe located.
