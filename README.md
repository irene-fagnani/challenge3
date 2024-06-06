# Challenge 3: A matrix–free parallel solver for the Laplace equation

### Introduction

Consider the Laplace equation for modeling heat diffusion over a square domain with a prescribed temperature (Dirichlet conditions) on the entire boundary:

$\left\{ \begin{array}{cl}
-\Delta u =f(x) & in \quad \Omega=(0,1)^2  \\
u=0 & on\quad  \left\{x=0\right\} \\
u=0 & on\quad  \left\{x=1\right\} \\
u=0 & on\quad  \left\{y=0\right\} \\
u=0 & on\quad  \left\{y=1\right\} \\
\end{array} \right.$

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
