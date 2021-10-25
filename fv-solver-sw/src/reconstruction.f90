!===============================================================================!
MODULE MOD_Reconstruction
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ReconstructionX
  MODULE PROCEDURE ReconstructionX
END INTERFACE

INTERFACE ReconstructionY
  MODULE PROCEDURE ReconstructionY
END INTERFACE

INTERFACE ReconstructionFixX
  MODULE PROCEDURE ReconstructionFixX
END INTERFACE

INTERFACE ReconstructionFixY
  MODULE PROCEDURE ReconstructionFixY
END INTERFACE

INTERFACE MUSCL 
  MODULE PROCEDURE MUSCL 
END INTERFACE

INTERFACE WENO3_SecondSweep 
  MODULE PROCEDURE WENO3_SecondSweep
END INTERFACE

INTERFACE WENO5_SecondSweep 
  MODULE PROCEDURE WENO5_SecondSweep
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ReconstructionX
PUBLIC :: ReconstructionY
PUBLIC :: ReconstructionFixX
PUBLIC :: ReconstructionFixY
PUBLIC :: MUSCL
PUBLIC :: WENO3_SecondSweep 
PUBLIC :: WENO5_SecondSweep 
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
SUBROUTINE ReconstructionX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        WM(1:nVar,nGPs,ii,jj) = V(1:nVar,ii,jj)
        WP(1:nVar,nGPs,ii,jj) = V(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        CALL MUSCL(&
                  V(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
      END DO
    END DO
  CASE(3,4)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (.NOT. Ind(1,ii,jj)) THEN
          CALL WENO_XDIR(&
                    V(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1),&
                    Reconstruction)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        WM(1:nVar,nGPs,ii,jj) = V(1:nVar,ii,jj)
        WP(1:nVar,nGPs,ii,jj) = V(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        CALL MUSCL(&
                   V(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                   WM(1:nVar,1:nGPs,ii,jj),&
                   WP(1:nVar,1:nGPs,ii,jj),&
                   MESH_DX(2))
      END DO
    END DO
  CASE(3,4)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (.NOT. Ind(2,ii,jj)) THEN
          CALL WENO_YDIR(&
                    V(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2),&
                    Reconstruction)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

IF ((Reconstruction .EQ. 1) .OR. (Reconstruction .EQ. 2)) THEN
  RETURN
END IF

SELECT CASE (ReconstructionFix)
  CASE(1)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj) = V(1:nVar,ii,jj)
              WP(1:nVar,iGP,ii,jj) = V(1:nVar,ii,jj)
          END DO
        END IF
      END DO
    END DO
  CASE(2)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          CALL MUSCL(&
                    V(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        END IF
      END DO
    END DO
  CASE(3,4)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          CALL WENO_XDIR(&
                    V(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1),&
                    ReconstructionFix)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "ReconstructionFix not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

IF ((Reconstruction .EQ. 1) .OR. (Reconstruction .EQ. 2)) THEN
  RETURN
END IF

SELECT CASE (ReconstructionFix)
  CASE(1)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj) = V(1:nVar,ii,jj)
              WP(1:nVar,iGP,ii,jj) = V(1:nVar,ii,jj)
          END DO
        END IF
      END DO
    END DO
  CASE(2)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          CALL MUSCL(&
                    V(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        END IF
      END DO
    END DO
  CASE(3,4)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          CALL WENO_YDIR(&
                    V(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2),&
                    ReconstructionFix)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "ReconstructionFix not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE MUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: sm, sp, slope
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    sp    = (Q(iVar,+1) - Q(iVar,+0))/dx
    sm    = (Q(iVar,+0) - Q(iVar,-1))/dx
    slope = MINMOD(sm,sp)

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD(x,y)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x, y
REAL            :: MINMOD
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

MINMOD = 0.5*(SIGN(1.0,x) + SIGN(1.0,y))*MIN(ABS(x),ABS(y))

!-------------------------------------------------------------------------------!
END FUNCTION MINMOD
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO_XDIR(V,WM,WP,dx,WhichReconstruction)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: dx
REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
INTEGER            :: iVar, ii, jj, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE(WhichReconstruction)
  CASE(3)
    DO iVar=1,nVar
      DO jj=-nGhosts,nGhosts
        CALL WENO3_FirstSweep(&
          V(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
      END DO
      CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    END DO
  CASE(4)
    DO iVar=1,nVar
      DO jj=-nGhosts,nGhosts
        CALL WENO5_FirstSweep(&
          V(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
      END DO
      CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO_XDIR
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO_YDIR(V,WM,WP,dy,WhichReconstruction)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: dy
REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
INTEGER            :: iVar, ii, jj, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE(WhichReconstruction)
  CASE(3)
    DO iVar=1,nVar
      DO ii=-nGhosts,nGhosts
        CALL WENO3_FirstSweep(&
          V(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
      END DO
      CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    END DO
  CASE(4)
    DO iVar=1,nVar
      DO ii=-nGhosts,nGhosts
        CALL WENO5_FirstSweep(&
          V(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
      END DO
      CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO_YDIR
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO3_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 2.0/3.0
gamma2 = 1.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

W1 = 0.5*(    Q(-1) + Q(+0))
W2 = 0.5*(3.0*Q(+0) - Q(+1))
WM = omega1*W1 + omega2*W2


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/3.0
gamma2 = 2.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1 = 0.5*(-Q(-1) + 3.0*Q(+0))
W2 = 0.5*( Q(+0) +     Q(+1))
WP = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO3_SecondSweep(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(SQRT(3.0)*Q(-1) + 6.0*Q(+0) - SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(SQRT(3.0)*Q(+0) + 6.0*Q(+0) - SQRT(3.0)*Q(+1))
W(1) = omega1*W1 + omega2*W2


!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(-SQRT(3.0)*Q(-1) + 6.0*Q(+0) + SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(-SQRT(3.0)*Q(+0) + 6.0*Q(+0) + SQRT(3.0)*Q(+1))
W(2) = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3_SecondSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 3.0/10.0
gamma2 = 3.0/5.0
gamma3 = 1.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(    -Q(-2) + 5.0*Q(-1) + 2.0*Q(+0))
W2 = (1.0/6.0)*( 2.0*Q(-1) + 5.0*Q(+0) -     Q(+1))
W3 = (1.0/6.0)*(11.0*Q(+0) - 7.0*Q(+1) + 2.0*Q(+2))
WM = omega1*W1 + omega2*W2 + omega3*W3


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/10.0
gamma2 = 3.0/5.0
gamma3 = 3.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(2.0*Q(-2) - 7.0*Q(-1) + 11.0*Q(+0))
W2 = (1.0/6.0)*(   -Q(-1) + 5.0*Q(+0) +  2.0*Q(+1))
W3 = (1.0/6.0)*(2.0*Q(+0) + 5.0*Q(+1) -      Q(+2))
WP = omega1*W1 + omega2*W2 + omega3*W3

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep2nGPs(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = (210.0 + SQRT(3.0))/1080.0
gamma2 = 11.0/18.0
gamma3 = (210.0 - SQRT(3.0))/1080.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = -     SQRT(3.0)*Q(-2) &
         + 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0)
W2     = + 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 1.0*SQRT(3.0)*Q(+1)
W3     = +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0) &
         - 4.0*SQRT(3.0)*Q(+1) &
         +     SQRT(3.0)*Q(+2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0


!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = (210.0 - SQRT(3.0))/1080.0
gamma2 = 11.0/18.0
gamma3 = (210.0 + SQRT(3.0))/1080.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = + 1.0*SQRT(3.0)*Q(-2) &
         - 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0)
W2     = - 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         +     SQRT(3.0)*Q(+1)
W3     = +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0) &
         + 4.0*SQRT(3.0)*Q(+1) &
         -     SQRT(3.0)*Q(+2)
W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep2nGPs
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep3nGPs(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!------------------------------!
! Point: x_{j-1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
gamma1 = (71.0*SQRT(15.0))/5240.0 + 126.0/655.0
gamma2 = 403.0/655.0
gamma3 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = (1.0/30.0 - SQRT(15.)/20.0)*Q(-2) &
         + (SQRT(15.)/5. - 1./15.)*Q(-1) &
         + (31./30. - (3.*SQRT(15.))/20.)*Q(+0)
W2     = + (SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. - SQRT(15.)/20.)*Q(+1)
W3     = + (3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         - (SQRT(15.)/5. + 1./15.)*Q(+1) &
         + (SQRT(15.)/20. + 1./30.)*Q(+2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)


!------------------------------!
! Point: x_{j}                 !
!------------------------------!

! Linear Weights
gamma1 = -9./80.
gamma2 = 49./40.
gamma3 = -9./80.

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = -1./24. *Q(-2) &
         +1./12. *Q(-1) &
         +23./24. *Q(+0) 
W2     = - 1./24. *Q(-1) &
         + 13./12.*Q(+0) &
         - 1./24. *Q(+1)
W3     = + 23./24.*Q(+0) &
         + 1. /12.*Q(+1) &
         - 1./24. *Q(+2)

W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)



!------------------------------!
! Point: x_{j+1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
gamma1 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0
gamma2 = 403.0/655.0
gamma3 = +(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = (1.0/30.0 + SQRT(15.)/20.0)*Q(-2) &
         - (SQRT(15.)/5. + 1./15.)*Q(-1) &
         + (31./30. + (3.*SQRT(15.))/20.)*Q(+0)
W2     = + (-SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. + SQRT(15.)/20.)*Q(+1)
W3     = + (-3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         + (SQRT(15.)/5. - 1./15.)*Q(+1) &
         + (-SQRT(15.)/20. + 1./30.)*Q(+2)
W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep3nGPs
!===============================================================================!
!
!!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep(Q,W) !4nGPS
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.2658420974778319
gamma2 =0.6112504900322486
gamma3 =0.1229074124899195

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = -0.1642562761719537 * Q(-2)+0.7590807081409336 * Q(-1)+0.4051755680310201 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.2663118796250726 * Q(-1)+0.8979443965468811 * Q(0)-0.1642562761719537 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.6968800354220990 * Q(0)-0.9631919150471715 * Q(1)+0.2663118796250726 * Q(2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

gamma1 =0.1281641584355059
gamma2 =0.5219691498503563
gamma3 =0.3498666917141379

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)


W1 = -0.1122135388132497 * Q(-2)+0.3944175994189276 * Q(-1)+0.7177959393943222 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.0577769829791784 * Q(-1)+1.0544365558340714 * Q(0)-0.1122135388132497 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.2277675047716066 * Q(0)-0.2855444877507849 * Q(1)+0.0577769829791784 * Q(2)

W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)


!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.3498666917141379
gamma2 =0.5219691498503563
gamma3 =0.1281641584355059

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial

W1 = +0.0577769829791784 * Q(-2)-0.2855444877507849 * Q(-1)+1.2277675047716066 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1122135388132497 * Q(-1)+1.0544365558340714 * Q(0)+0.0577769829791784 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.7177959393943222 * Q(0)+0.3944175994189276 * Q(1)-0.1122135388132497 * Q(2)

W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.1229074124899195
gamma2 =0.6112504900322486
gamma3 =0.2658420974778319


alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = +0.2663118796250726 * Q(-2)-0.9631919150471715 * Q(-1)+1.6968800354220990 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1642562761719537 * Q(-1)+0.8979443965468811 * Q(0)+0.2663118796250726 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.4051755680310201 * Q(0)+0.7590807081409336 * Q(1)-0.1642562761719537 * Q(2)
W(4)   = (omega1*W1 + omega2*W2 + omega3*W3)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
END MODULE MOD_Reconstruction
!-------------------------------------------------------------------------------!
