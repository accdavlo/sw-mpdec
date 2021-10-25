!===============================================================================!
MODULE MOD_PostProcessing
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE PostProcessing  
  MODULE PROCEDURE PostProcessing 
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: PostProcessing
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PostProcessing
USE MOD_Equation, ONLY :  ExactFunction
USE MOD_FiniteVolume2D_vars
IMPLICIT NONE
INTEGER              :: ii, jj, iGP, jGP
REAL, ALLOCATABLE    :: Uexact(:,:,:), L1error(:)
CHARACTER(LEN=255)   :: FNAME
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp

ALLOCATE( L1error(1:nVar) )
ALLOCATE( Uexact(1:nVar,1:nElemsX,1:nElemsY))
 
Uexact(:,:,:) = 0.
L1error(:) = 0. 
Utemp = 0.
DO jj = 1,nElemsY  
  DO ii = 1,nElemsX  
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tEnd, MeshGP(:,ii,jj,iGP,jGP), Utemp(1:nVar,iGP,jGP))
        Uexact(1:nVar,ii,jj) = Uexact(1:nVar,ii,jj) + WeightsGP(iGP,jGP)* Utemp(1:nVar,iGP,jGP)
      END DO
    END DO
    L1error(:) = L1error(:) + ABS( U(:,ii,jj) - Uexact(:,ii,jj) ) 
  ENDDO
ENDDO

L1error(:) = L1error(:) * MESH_DX(1) * MESH_DX(2)

FNAME = 'ErrorL1_XXXX_XXXX.dat'
WRITE(FNAME(9:12),FMT="(I4.4)") nElemsX
WRITE(FNAME(14:17),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) L1error(1:nVar)
CLOSE(666)

DEALLOCATE( L1error )
DEALLOCATE( Uexact )

END SUBROUTINE

END MODULE MOD_PostProcessing
