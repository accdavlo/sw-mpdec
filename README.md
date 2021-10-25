# Modified Patankar Deferred Correction Application for Shallow Water equations

This is the public repository refering to the article [arXiv prepint to be put](https://arxiv.org) where we include the Fortran code which we used to develop the method.

This code is based on a finite volume WENO5 space discretization that is available at [FV-solver-SW](https://github.com/jbnunezd/fv-solver-sw.git) developed by Jonatan Núñez-de la Rosa.

## Compile the code

In order to compile the code the user can either go to [fv-solver-sw](fv-solver-sw) or [fv-solver-sw/src](fv-solver-sw/src)
from a terminal and type **make**. This will compile the code and generate the executable called **main** in [fv-solver-sw/bin](fv-solver-sw/bin).

## Run the test cases presented in the paper

In the paper, six test cases are presented:
1. Unsteady smooth vortex
2. Lake at rest (well-balanced)
3. Perturbation analysis (well-balanced, wet-dry)
4. Circular dam break #1 (wet-dry)
5. Circular dam break #2
6. Wave over dry island

In order to run them, the user can follow two simple procedures:

* go to [fv-solver-sw/bin](tree/main/fv-solver-sw/bin) and type 
**./make TestCase NElems** 
where **TestCase** is the number of the testcase and **NElems** will generate a mesh *NElems x NElems*;

* manually open the script [parameters.f90](fv-solver-sw/src/bin/parameters.f90) and modify the setup considering that

*InitialCondition* - the number of the test case
*TEnd* - the final time of the simulation
*nElemsX* (*nElemsY*) - the number of elements along x-direction (y-direction)
*MESH_X0* and *MESH_X1* allow to impose the computational domain (either a square or a rectangle)
*BoundaryConditionsType* - the color of the boundary condition to apply on the four boundaries (so far mPDeC only work with periodic BC, *BoundaryConditionsType=1*)
*BathymetryFlag* - the parameter used to apply different bathymetry (more bathymetries can be implemented in [equation.f90](fv-solver-sw/src/bin/equation.f90))
*CFL* - the Courant-Friedrich-Levy number
*Reconstruction* - space discretization (1. first order, 2. second order TVD, 3. WENO3, 4. WENO5)
*ReconstructionFix* - space discretization where shocks occur, if necessary (elements are flagged with a shock indicator implemented in [shocksindicator.f90](fv-solver-sw/src/bin/shocksindicator.f90)) 
*timescheme* - time discretization (1. explicit Euler, 2. SSPRK64, 3. RK65, 4. DeC5, 5. PatankarDeC5, 6. PatankarEuler, 7. mPDeC2)
*nOutputFiles* - number of output files printed
