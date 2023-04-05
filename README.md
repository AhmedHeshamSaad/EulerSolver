# EulerSolver using JSK scheme for Inviscid flows

Language used: Fortran

Description: 
The project's goal is to develop a CFD solver for inviscid flows. It is an excellent opportunity to learn advanced computational methods, PDEs, and computer science.

The grid is loaded from a vtk file which can be generated using the [StrGridBump code](https://github.com/AhmedHeshamSaad/StrGridBUMP).
The results will be also written in a new vtk file including all the state variables (&rho;, u, v, p)

It is a finite volume solver that use the Jameson-Schmidt-Turkel scheme. It models the dissipation source terms (2nd and 4th order viscosity terms) to solve the flow variable discontinuities across shock waves.

### Inlet BC:              
The inlet boundary state vector was computed using Riem1 at the far stream and Riem2 at the interior's first cells.
### Outlet BC:
The static pressure at the outlet was set to be equal the static pressure of the far stream, while density and velocity was computed assuming zero gradient.
### Top and bottom BCs:
The top and bottom boundary was defined as a wall where ghost cells were used that have mirror image of the interior velocity.

### The Temporal Discretization:
using Runge-Kutta algorithm with 4 steps.

### Computational Domain:
This project makes use of a mesh developed by previous in-house code ([StrGridBump](https://github.com/AhmedHeshamSaad/StrGridBUMP)).
The bump is shaped like a sin function and has an amplitude of 0.1 of its length.
The images below depict Mach contours for mach numbers of 0.3, 0.7 and 0.9.
Transparent grid lines are also shwon.

M = 0.3
![M = 0.3](contours/M0p3.jpg?raw=true)

M = 0.7
![M = 0.7](contours/0p7.png?raw=true)

M = 0.9
![M = 0.9](contours/0p9.png?raw=true)

M = 1.2
![M = 1.2](contours/1p2.png?raw=true)