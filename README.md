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

In order to run them, the user can follow two simple procedures

go to [fv-solver-sw/bin](tree/main/fv-solver-sw/bin) and type 
**./make** *TestCase* *NElems* 

manually open the script and set the parameters.
