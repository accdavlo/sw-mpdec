!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE ExactFunctionWB
  MODULE PROCEDURE ExactFunctionWB
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE Bathymetry   
  MODULE PROCEDURE Bathymetry   
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: ExactFunctionWB
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: Bathymetry
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: xc(2), xm(2), r, r0, hl, hr, r2, r20
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!
!*OUR VARIABLES
REAL               :: Omega, Jamma, u_inf, v_inf, h_inf, DeltaH
INTEGER            :: power




Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  CASE(1) !*UNSTEADY SMOOTH VORTEX
    u_inf = 2.
    v_inf = 3.
    H_inf=1.
    r0 = 1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = (xc(1)**2 + xc(2)**2)
    Omega = sqrt(2.*Gravity*hDerivSmoothAuxiliary(r))
    
    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. 1) THEN
      Prim(1) = hSmoothAuxiliary(r)
      Prim(2) = Prim(2)+Omega*(+xc(2))
      Prim(3) = Prim(3)+Omega*(-xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)

  CASE(2) !*LAKE AT REST
    Prim(1) = 1. - 0.1 * SIN(2.*PI*x(1)) * COS(2.*PI*x(2))
    Prim(2) = 0.
    Prim(3) = 0.

    CALL PrimToCons(Prim,Cons)

  CASE(3) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST 
    Prim(1) = MAX(0.7-Bathymetry(x),MIN_DEPTH) 
    r2 = ((x(1)+2)**2 +(x(2)-0.5)**2)*9.
    IF (r2<1) THEN
      Prim(1) = Prim(1) + 0.05*EXP(1.-1./(1-r2)**2.)
    ENDIF
    Prim(2:3) = 0.0

    CALL PrimToCons(Prim,Cons)

  CASE(4) ! *CIRCULAR DAM BREAK 1
    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm(1)
    xc(2) = x(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    r0    = 7.

    Prim(1) = MIN_DEPTH
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF (r .LE. r0) THEN
      Prim(1) = 2.5
      Prim(2) = 0.0
      Prim(3) = 0.0
    END IF
    CALL PrimToCons(Prim,Cons)

  CASE(5) ! *CIRCULAR DAM BREAK 2
    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm(1)
    xc(2) = x(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    r0    = 15.

    Prim(1) = 0.5 
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF (r .LE. r0) THEN
      Prim(1) = 10.
      Prim(2) = 0.0
      Prim(3) = 0.0
    END IF
    CALL PrimToCons(Prim,Cons)

  CASE(6) ! *WAVE OVER DRY ISLAND 
    r2 = (x(1)+2)**2.
    Prim(1) = MAX(0.7-Bathymetry(x),MIN_DEPTH) 
    IF (r2<1) THEN
      Prim(1) =Prim(1) + 0.5*EXP(1.-1./(1-r2)**2.)
    ENDIF
    IF (Prim(1)>MIN_DEPTH*10) THEN
      Prim(2) = 1.0
    ELSE
      Prim(2) = 0.0
    ENDIF
    Prim(3) = 0.0

    CALL PrimToCons(Prim,Cons)

  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
CONTAINS
   
  REAL FUNCTION hSmoothAuxiliary(x)
     REAL, INTENT(IN) :: x

     hSmoothAuxiliary=1.-0.1*exp(-1./atan(1.-x)**3.)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary(x)
     REAL, INTENT(IN) :: x

     hDerivSmoothAuxiliary=3.*0.1*exp(1./atan(x - 1.)**3.)/(atan(x - 1.)**4*((x - 1.)**2 + 1.))

  END FUNCTION

END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ExactFunctionWB(WhichInitialCondition,x,Cons)
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI 
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage

SELECT CASE (WhichInitialCondition)
  
  CASE(2) !*LAKE AT REST
  
    Prim(1) = 1. - 0.1 * SIN(2.*PI*x(1)) * COS(2.*PI*x(2))
    Prim(2) = 0.
    Prim(3) = 0.

    CALL PrimToCons(Prim,Cons)

  CASE(3) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST 
    Prim(1) = MAX(0.7-Bathymetry(x),MIN_DEPTH) 
    Prim(2:3) = 0.0

    CALL PrimToCons(Prim,Cons)

  CASE DEFAULT
    ErrorMessage = "Exact WB function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

END SUBROUTINE ExactFunctionWB
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep 
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: SW(1:nVar,nGPs,nGPs,nElemsX,nElemsY)
REAL             :: Vtemp(1:nVar,-nGhosts:nElemsX+nGhosts+1,1:nElemsY,1:nGPs)
REAL             :: Vtemp2(1:nVar,1:nGPs,1:nGPs,1:nElemsX,1:nElemsY)
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

S = 0.0 ! S(1:nVar,nElemsX,nElemsY)

!int_{iixjj} S(1:nVar,ii,jj) dxdy

IF (BathymetryFlag .GT. 0) THEN
  SELECT CASE (Reconstruction)
    CASE(1,2)
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          SW(1:nVar,nGPs,nGPs,ii,jj) = SourceFunc( U(1:nVar,ii,jj) , MeshBary(:,ii,jj) )
        END DO
      END DO
    CASE(3)
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO3_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO3_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              SW(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(4)
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO5_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO5_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              SW(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT


  DO jj=1,nElemsY
    DO ii=1,nElemsX
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          S(1:nVar,ii,jj) = S(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * SW(1:nVar,iGP,jGP,ii,jj) 
        END DO 
      END DO 
    END DO
  END DO

END IF
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION SourceFunc(Q,X) RESULT(Source_SW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
IMPLICIT NONE
REAL, DIMENSION(1:nVar), INTENT(IN)  :: Q 
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL, DIMENSION(1:nVar) :: Source_SW 

Source_SW(1) = 0.
Source_SW(2) = -Gravity*Q(1)*Bathymetry_X(X) 
Source_SW(3) = -Gravity*Q(1)*Bathymetry_Y(X) 

!-------------------------------------------------------------------------------!
END FUNCTION SourceFunc
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Bathymetry(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (BathymetryFlag)
  CASE (1) ! used for *LAKE AT REST
    Bathymetry = 0.1 * SIN(2.*PI*X(1)) * COS(2.*PI*X(2)) 

  CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND 
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry = EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry = 0.0
    ENDIF

  CASE DEFAULT
   Bathymetry = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Bathymetry
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Bathymetry_X(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (BathymetryFlag)
  CASE (1) ! used for *LAKE AT REST
    Bathymetry_X = 0.1 * 2. * PI * COS(2.*PI*X(1)) * COS(2.*PI*X(2)) 

  CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry_X = -2.0*X(1)/(1.0-r2)**2.*EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry_X = 0.0
    ENDIF

  CASE DEFAULT
   Bathymetry_X = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Bathymetry_X
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Bathymetry_Y(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (BathymetryFlag)
  CASE (1) ! used for *LAKE AT REST
    Bathymetry_Y = - 0.1 * 2. * PI * SIN(2.*PI*X(1)) * SIN(2.*PI*X(2))

  CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry_Y = -2.0*X(2)/(1.0-r2)**2.*EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry_Y = 0.0
    ENDIF

  CASE DEFAULT
   Bathymetry_Y = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Bathymetry_Y
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Upper Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO

!  CASE DEFAULT
!    ErrorMessage = "Boundary condition not implemented"
!    WRITE(*,*) ErrorMessage
!    STOP
END SELECT


!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO

    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
!  CASE DEFAULT
!    ErrorMessage = "Boundary condition not implemented"
!    WRITE(*,*) ErrorMessage
!    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)           :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowest
REAL,INTENT(OUT),OPTIONAL :: fastest
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                      :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

IF(PRESENT(slowest)) THEN
  slowest = vx - SQRT(Gravity*h)
END IF

IF(PRESENT(fastest)) THEN
  fastest = vx + SQRT(Gravity*h)
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

fastestx = vx + SQRT(Gravity*h)
fastesty = vy + SQRT(Gravity*h)

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, hvx, hvy, ht
!-------------------------------------------------------------------------------!

h    = Cons(1)
hvx  = Cons(2)
hvy  = Cons(3)

!IF (h .LT. MIN_DEPTH) THEN
!  h = MIN_DEPTH
!END IF

ht = h + MIN_DEPTH/h

Prim(1) = h
Prim(2) = hvx/ht
Prim(3) = hvy/ht

!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h   = Prim(1)
vx  = Prim(2)
vy  = Prim(3)

!IF (h .LT. MIN_DEPTH) THEN
!  h = MIN_DEPTH
!END IF

Cons(1) = h
Cons(2) = h*vx
Cons(3) = h*vy

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

!IF (h .LT. MIN_DEPTH) THEN
!  h = MIN_DEPTH
!END IF

Flux(1) = h*vx
Flux(2) = h*vx**2 + 0.5*Gravity*h**2
Flux(3) = h*vx*vy

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims,1:nGPs)
REAL,INTENT(IN)  :: TangVect(1:nDims,1:nGPs)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  ! Rotating the vector quantities       !

  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1,iGP)*PrimL(2,iGP) + NormVect(2,iGP)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1,iGP)*PrimL(2,iGP) + TangVect(2,iGP)*PrimL(3,iGP)

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1,iGP)*PrimR(2,iGP) + NormVect(2,iGP)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1,iGP)*PrimR(2,iGP) + TangVect(2,iGP)*PrimR(3,iGP)

  CALL PrimToCons(PrimLL(1:nVar,iGP),ConsLL(1:nVar,iGP))
  CALL PrimToCons(PrimRR(1:nVar,iGP),ConsRR(1:nVar,iGP))  

  CALL RiemannSolverByRusanov(&
    ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
    PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))

  ! Rotating back the vector quantities  !
  Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
                + TangVect(1:nDims,iGP)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastestL, fastestR
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)
CALL WaveSpeeds1D(PrimL,fastest=fastestL)
CALL WaveSpeeds1D(PrimR,fastest=fastestR)

LambdaMax = MAX(ABS(fastestL),ABS(fastestR))

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
