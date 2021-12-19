# Modified Patankar Deferred Correction WENO Code for Shallow Water Equations

This is the public repository refering to the article [arXiv:2110.13509](https://arxiv.org/abs/2110.13509#) where we include the Fortran code which we used to develop the method.

This code is based on a finite volume WENO5 space discretization that is available at [FV-solver-SW](https://github.com/jbnunezd/fv-solver-sw.git) developed by Jonatan Núñez-de la Rosa.

![Island](pictures/island.gif)

## Compile the code

In order to compile the code the user can either go to [fv-solver-sw](fv-solver-sw) or [fv-solver-sw/src](fv-solver-sw/src)
from a terminal and type **make**. This will compile the code and generate the executable called **main** in [fv-solver-sw/bin](fv-solver-sw/bin). Up to this point, the makefile is generated for Linux systems using GNU Debug. For Mac systems, the path to the libraries must be set correctly in BLAS_LAPACK.

Depending on the **FCFLAGS** set in the [makefile](fv-solver-sw/src/makefile), the code can compiled in different ways. Two **FCFLAGS** are already available:
* *WELLBALANCED* - subtract fluxes and source term corresponding to the exact solution one wants to preserve (new solutions can be set in the subroutine ExactFunctionWB in [equation.f90](fv-solver-sw/src/equation.f90))
* *PATANKAR*     - recast the hyperbolic system as a production-destruction system (it is mandatory to set a Patankar scheme for time discretization)
## Run the test cases presented in the paper 

In the paper, six test cases are presented:
1. Unsteady smooth vortex (convergence)
1. Lake at rest (well-balanced, convergence)
1. Perturbation analysis (well-balanced, wet-dry)
1. Circular dry dam break (wet-dry)
1. Circular wet dam break
1. Wave over dry island (wet-dry)

In order to run them, the user can follow two simple procedures:

* go to [fv-solver-sw/bin](fv-solver-sw/bin) and type 
**./make TestCase NElems** 
where **TestCase** is the number of the testcase and **NElems** will generate a mesh *NElems x NElems*;

* manually open the script [parameters.f90](fv-solver-sw/src/parameters.f90) and modify the setup considering that  

  * **InitialCondition**       - the number of the test case     
  * **TEnd**                   - the final time of the simulation   
  * **nElemsX (nElemsY)**      - the number of elements along x-direction (y-direction)    
  * **MESH_X0** and **MESH_X1**- impose the computational domain (either a square or a rectangle)   
  * **BoundaryConditionsType** - the color of the boundary condition to apply on the four boundaries (so far mPDeC only work with periodic BC, *BoundaryConditionsType=1*)   
  * **BathymetryFlag**         - the parameter used to apply different bathymetries (new ones can be implemented in [equation.f90](fv-solver-sw/src/equation.f90))   
  * **CFL**                    - the Courant-Friedrich-Levy condition  
  * **Reconstruction**         - space discretization (1. first order, 2. second order TVD, 3. WENO3, 4. WENO5)    
  * **ReconstructionFix**      - space discretization where shocks occur, if necessary (elements are flagged with a shock indicator implemented in [shocksindicator.f90](fv-solver-sw/src/shocksindicator.f90))     
  * **timescheme**             - [time discretization](fv-solver-sw/src/timediscretization.f90) (1. explicit Euler, 2. SSPRK64, 3. RK65, 4. DeC5, 5. mPDeC5, 6. mPEuler, 7. mPDeC2), here mP stands for modified Patankar
  * **nOutputFiles**           - number of output files printed    

By appropriately modifying these parameter and adding new codes, new test cases can be easily added to those already implemented.

## Reconstruction of the primitive variable

In the directory [weno-weights](weno-weights) the user can find two matlab scripts ([WENO_sym.m](weno-weights/WENO_sym.m) and 
[WENO_double.m](weno-weights/WENO_double.m)) used by the authors to compute all the ingredients needed for a WENO reconstruction. 
In particular, the two scripts provide the same results but the first one prints everything with symbolic notation whereas the second one only prints
what it is needed for the reconstruction as a double. In both scripts, the user can set the number of cell averages, **J**, (odd number) and the quadrature 
points, **x_quads**, where the reconstruction has to be computed (already few gaussian quadrature rules are specified, the code is setup with the 
four-point rule). The codes print out
* smoothness indicators
* coefficients of lower order polynomials
* ideal weights for the linear high order reconstruction.



