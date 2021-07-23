# Structured Grid Generator for Flow Around a Bump

Language used: Fortran

Description: 
The project objective is to create a mesh generator to initialize a structured grid for flow passes through a bump. The mesh obtained will be used further in another project to solve the Euler equation.

The bump has the shape of a sin function with an amplitude of 0.2. 

To generate an orthogonal grid, one solves an elliptic PDE which taking the form of a Laplacian equation. The code first initializes an algebraic grid for the defined domain around a bump. Then, that algebraic grid is used as an initial condition to the Laplacian equation. To numerically solve the Laplacian equation, the Gauss-Seidel method has been chosen as the iterative technique along with a central point discretization scheme. Adequate boundary conditions have been defined.

A sample of algebraic grid and Laplacian grid with a dimension of 90x15 nodes is shown in the figure below.

![alt text](https://raw.githubusercontent.com/AhmedHeshamSaad/StrGridBUMP/main/algebric_and_laplace_grid_90x15.jpeg?raw=true)