!===============================================================================!
MODULE MOD_JacobiIteration
!===============================================================================!

IMPLICIT NONE

PRIVATE

INTERFACE jacobi
  MODULE PROCEDURE jacobi
END INTERFACE

INTERFACE sparse_matmul
  MODULE PROCEDURE sparse_matmul
END INTERFACE

PUBLIC :: jacobi
PUBLIC :: sparse_matmul

CONTAINS

 !*------------------------------------------------------------------------------------------------
 !* Compute the multiplication matrix * vector where matrix is given as CRS
 !*------------------------------------------------------------------------------------------------
 SUBROUTINE sparse_matmul(rowptr,colptr,matrixvalues, rhs, product)
   IMPLICIT NONE
   REAL, DIMENSION(:), INTENT(in)    :: matrixvalues, rhs
   INTEGER, DIMENSION(:), INTENT(in) :: rowptr, colptr
   REAL, DIMENSION(:), INTENT(inout) :: product

   INTEGER :: n, nnz, row, spIdx
 
   n=SIZE(rhs,DIM=1)
   nnz=SIZE(matrixvalues,DIM=1)
   product = 0.
   DO row=1,n
      DO spIdx=rowptr(row),rowptr(row+1)-1
         product(row)=product(row)+matrixvalues(spIdx)*rhs(colptr(spIdx))
      ENDDO
   ENDDO


  END SUBROUTINE sparse_matmul


  
 !*------------------------------------------------------------------------------------------------
 !* Jacobi iteration algorithm to find solution of (D+L)x=b
 !* The iterations are defined as x^(k+1)=D^-1(b-L*x^k)
 !* Input: L as a CRS (row, col, values), D as vector, b vector, sol vector output
 !*------------------------------------------------------------------------------------------------
 SUBROUTINE jacobi(rowptr,colptr,matrixvalues, diag, rhs, sol)
   IMPLICIT NONE
   REAL, DIMENSION(:), INTENT(in)    :: matrixvalues, rhs, diag
   INTEGER, DIMENSION(:), INTENT(in) :: rowptr, colptr
   REAL, DIMENSION(:), INTENT(inout) :: sol
   REAL, DIMENSION(SIZE(rhs,DIM=1))  :: solp
   REAL                              :: residual, tolerance

   INTEGER :: n, nnz, k, maxIter, i
   n=SIZE(rhs,DIM=1)

   tolerance = 1.e-16
   maxIter = 1000
   
   k=0
   solp=rhs
   sol=0.
   residual = SUM(ABS(sol-solp))/n

   DO WHILE ( (tolerance<residual).AND. (k<maxIter))
      k=k+1
      CALL sparse_matmul(rowptr,colptr,matrixvalues,solp,sol)  !sol = L x^k
      DO  i=1, n
         sol(i)=(rhs(i)-sol(i))/diag(i)
      ENDDO 
      residual = SUM(ABS(sol-solp))/n
      !PRINT*, "res = ", residual

      solp=sol
   ENDDO

   IF (tolerance<residual) THEN
      print*, "Jacobi not converged. residual = ", residual
   ENDIF

  END SUBROUTINE jacobi


END MODULE MOD_JacobiIteration
!-------------------------------------------------------------------------------!
