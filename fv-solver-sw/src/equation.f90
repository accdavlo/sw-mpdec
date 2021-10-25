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

  CASE(100) !Island
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

  CASE(101) !Island at rest
    Prim(1) = MAX(0.7-Bathymetry(x),MIN_DEPTH) 
    Prim(2:3) = 0.0

    CALL PrimToCons(Prim,Cons)

  CASE(102) !Island at rest with perturbation
    Prim(1) = MAX(0.7-Bathymetry(x),MIN_DEPTH) 
    r2 = ((x(1)+2)**2 +(x(2)-0.5)**2)*9.
    IF (r2<1) THEN
      Prim(1) = Prim(1) + 0.05*EXP(1.-1./(1-r2)**2.)
    ENDIF
    Prim(2:3) = 0.0

    CALL PrimToCons(Prim,Cons)

  CASE(200)
    ! Constant State
    Prim(1) = 1.0
    Prim(2) = 0.0
    Prim(3) = 0.0
    CALL PrimToCons(Prim,Cons)
  
  CASE(211)
    ! Dam Break
    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)

    hl    = 2.5
    hr    = MIN_DEPTH!0.5


    IF (x(1) .LE. xm(1)) THEN
      Prim(1) = hl
      Prim(2) = 0.0
      Prim(3) = 0.0
    ELSE
      Prim(1) = hr
      Prim(2) = 0.0
      Prim(3) = 0.0
    END IF

    CALL PrimToCons(Prim,Cons)
  CASE(212)
    ! Circular Dam Break
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
  CASE(213)
    ! Double Dam Break with obstacle

    Prim(1) = MIN_DEPTH
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF ((x(1) .LE. 5.) .OR. (x(1) .GE. 35.)) THEN
      Prim(1) = 5.
    END IF
    CALL PrimToCons(Prim,Cons)
  CASE(214)
    ! Mario's Circular Dam Break
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
  CASE(221)
    ! Double Dam Break with obstacle

    Prim(1) = MIN_DEPTH
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF ((x(1) .LE. 6.) .OR. (x(1) .GE. 34.)) THEN
      Prim(1) = 5.
    END IF
    CALL PrimToCons(Prim,Cons)
  CASE(222)
    ! Double Dam Break with obstacle

    Prim(1) = MIN_DEPTH
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF ((x(2) .LE. 6.) .OR. (x(2) .GE. 34.)) THEN
      Prim(1) = 5.
    END IF
    CALL PrimToCons(Prim,Cons)

  CASE(223)
    ! Double Dam Break with obstacle

    Prim(1) = MIN_DEPTH
    Prim(2) = 0.0
    Prim(3) = 0.0

    IF ((x(2) .LE. 26.) .AND. (x(2) .GE. 14.)) THEN
      Prim(1) = 5.
    END IF
    CALL PrimToCons(Prim,Cons)

  CASE(300) !*STEADY VORTEX
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    Omega = PI / r0
    DeltaH = 0.1

    Jamma=SQRT(Gravity)*PI*2.*SQRT(DeltaH)/SQRT(3.*PI**2 - 16.)/r0

    u_inf=0.
    v_inf=0.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction(Omega*r)-hfunction(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)

  CASE(301) !*UNSTEADY VORTEX
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    Omega = PI / r0

    Jamma=SQRT(10.*Gravity)*PI/(5.*SQRT(3.*PI**2 - 16.))/r0
    u_inf=1.
    v_inf=1.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)


    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction(Omega*r)-hfunction(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)


  CASE(302) !*STEADY VORTEX p=2
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    DeltaH = 0.1
    Omega = PI / r0

    Jamma=(12.*SQRT(DeltaH)*SQRT(Gravity)*PI)/(SQRT(315.*PI**2. - 2048.))/r0
    u_inf=0.
    v_inf=0.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    r0    = PI/Omega

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction2(Omega*r)-hfunction2(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))**2.*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))**2.*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)

  CASE(303) !*STEADY VORTEX p=3
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    Omega = PI / r0
    DeltaH = 0.1

    Jamma=(60.*SQRT(DeltaH)*SQRT(2.*Gravity)*PI)/SQRT(51975.*PI**2. - 367616.)/r0
    u_inf=0.
    v_inf=0.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction3(Omega*r)-hfunction3(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))**3.*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))**3.*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)


  CASE(333) !*UNSTEADY VORTEX p=3
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    Omega = PI / r0
    DeltaH = 0.1

    Jamma=(60.*SQRT(DeltaH)*SQRT(2.*Gravity)*PI)/SQRT(51975.*PI**2. - 367616.)/r0
    u_inf=1.
    v_inf=1.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction3(Omega*r)-hfunction3(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))**3.*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))**3.*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)


  CASE(304) !*STEADY VORTEX p=4 NOt implemented!!!
    !*PARAMETERS OF THE TEST
    r0 = 0.25
    Omega = PI / r0

    Jamma=(504.*SQRT(2*Gravity)*PI)/SQRT(4583103525.*PI**2 - 35168714752.)/r0
    u_inf=0.
    v_inf=0.
    H_inf=1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. r0) THEN
      Prim(1) = Prim(1)+1./Gravity*(Jamma/Omega)**2*(hfunction3(Omega*r)-hfunction3(PI))
      Prim(2) = Prim(2)+Jamma*(1+COS(Omega*r))**4.*(-xc(2))
      Prim(3) = Prim(3)+Jamma*(1+COS(Omega*r))**4.*(+xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)

  CASE(306,307,308,309) !*STEADY VORTEX smooth exp(-1/(1-r2)^p)
    !*PARAMETERS OF THE TEST
    r20 = 0.25**2.
    u_inf=0.
    v_inf=0.
    H_inf=1.
    DeltaH = 0.1
    SELECT CASE(WhichInitialCondition)
      CASE(306)
        power=2
      CASE(307)
        power=3
      CASE(308)
        power=4
      CASE(309)
        power=5
    END SELECT
    Jamma =exp(1.)*DeltaH

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r2     = (xc(1)**2 + xc(2)**2)/r20
    

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r2 .LT. 1) THEN
      Omega = sqrt(2.*Gravity*Jamma*hDerivSmoothAux(r2,power)/r20)
      Prim(1) = Prim(1) - Jamma*hSmoothAux(r2,power)
      Prim(2) = Prim(2)+Omega*(+xc(2))
      Prim(3) = Prim(3)+Omega*(-xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)


  CASE(292,293,294,295,299) !*STEADY VORTEX smooth exp(-1/atan(1-r2)^p)
    !*PARAMETERS OF THE TEST
    r20 = 0.25**2.
    u_inf=0.
    v_inf=0.
    H_inf=1.
    DeltaH = 0.1
    SELECT CASE(WhichInitialCondition)
      CASE(292)
        power=2
      CASE(293)
        power=3
      CASE(294)
        power=4
      CASE(295)
        power=5
      CASE(299)
        power=9
    END SELECT
    Jamma =exp((4./PI)**power)*DeltaH

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r2     = (xc(1)**2 + xc(2)**2)/r20
    

    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r2 .LT. 1) THEN
      Omega = sqrt(2.*Gravity*Jamma*hDerivSmoothTanAux(r2,power)/r20)
      Prim(1) = Prim(1) - Jamma*hSmoothTanAux(r2,power)
      Prim(2) = Prim(2)+Omega*(+xc(2))
      Prim(3) = Prim(3)+Omega*(-xc(1))
    END IF

    CALL PrimToCons(Prim,Cons)


  CASE(280,281,282,283,284) !*STEADY VORTEX smooth exp(-1/atan(1-r2)^p)
    !*PARAMETERS OF THE TEST

    SELECT CASE(WhichInitialCondition)
      CASE(283,284)
        u_inf=1.
        v_inf=1.
      CASE DEFAULT 
        u_inf=0.
        v_inf=0.
    END SELECT

    H_inf=1.
    DeltaH = 0.01
    SELECT CASE(WhichInitialCondition)
      CASE(280)
        r0 = 0.2
        Jamma =(2.*SQRT(15.*DeltaH))/3.
      CASE(281)
        r0 = 0.15
        Jamma = (2.*SQRT(30.*DeltaH))/3
      CASE(282,283)
        r0=0.1
        Jamma = 2.*SQRT(5.*DeltaH)
      CASE(284)
        r0=0.07
        Jamma = 2.*SQRT(Gravity/r0**2.*DeltaH)
    END SELECT

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = SQRT(xc(1)**2 + xc(2)**2)
    
    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    Omega = Jamma*EXP(-(r/r0)**2)
    Prim(1) = Prim(1) - r0**2/4./Gravity*Omega**2
    Prim(2) = Prim(2)+Omega*(+xc(2))
    Prim(3) = Prim(3)+Omega*(-xc(1))

    CALL PrimToCons(Prim,Cons)

  CASE(310) !*STEADY SMOOTH VORTEX
    !*PARAMETERS OF THE TEST
    u_inf=0.
    v_inf=0.
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
  CASE(311,3110) !*UNSTEADY SMOOTH VORTEX
    !*PARAMETERS OF THE TEST
    SELECT CASE (WhichInitialCondition)
    CASE (311)
      u_inf=3.
      v_inf=3.
    CASE (3110)
      u_inf = 2.
      v_inf = 3.
    END SELECT   
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

  CASE(319) !*SMOOTH TEST 
  
    Prim(1) = 10. + EXP(SIN(2.*PI*x(1)))*COS(x(2)*PI*2.)
    Prim(2) = SIN(COS(2.*PI*x(1)))*SIN(2.*PI*x(2))/Prim(1)
    Prim(3) = COS(SIN(2.*PI*x(2)))*COS(2.*PI*x(1))/Prim(1)


    CALL PrimToCons(Prim,Cons)

  CASE(320) !*SMOOTH TEST 
  
    Prim(1) = 1.+0.1*SIN(x(1)*PI*2.)*COS(x(2)*PI*2.)
    Prim(2) = 0.
    Prim(3) = 0.


    CALL PrimToCons(Prim,Cons)
  CASE(321) !*LAKE AT REST
  
    Prim(1) = 1. - 0.1 * SIN(2.*PI*x(1)) * COS(2.*PI*x(2))
    Prim(2) = 0.
    Prim(3) = 0.

    ! IF ((x(1).GT.0.5).AND.(x(1).LT.0.6).AND.(x(2).GT.0.5).AND.(x(2).LT.0.6)) THEN
    !   Prim(1) = Prim(1) + 0.05
    ! ENDIF

    CALL PrimToCons(Prim,Cons)
  CASE(400) !*LAKE AT REST NO BATHYMETRY
  
    Prim(1) = 1.
    Prim(2) = 0.
    Prim(3) = 0.
    IF ((x(1)-0.5)**2.+(x(2)-0.5)**2.<0.01) THEN
      Prim(1)=Prim(1) + 0.1
    ENDIF


    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
CONTAINS
   
  REAL FUNCTION hfunction(x)
     REAL, INTENT(IN) :: x

     hfunction=2.*COS(x)+2.*x*SIN(x)+1./8.*COS(2*x)+0.25*x*SIN(2*x)+0.75*x**2

  END FUNCTION

  REAL FUNCTION hfunction2(r)
     REAL, INTENT(IN) :: r

     hfunction2=(20.*cos(r))/3. + (27.*cos(r)**2.)/16. + (4.*cos(r)**3)/9.+ cos(r)**4/16. + (20.*r*sin(r))/3. &
     + (35.*r**2)/16. + (27.*r*cos(r)*sin(r))/8. + (4.*r*cos(r)**2*sin(r))/3. + (r*cos(r)**3*sin(r))/4.


  END FUNCTION

  REAL FUNCTION hfunction3(r)
     REAL, INTENT(IN) :: r

     hfunction3=(338.*cos(r))/15. + (215.*cos(r)**2)/32. + (124.*cos(r)**3)/45. + (95.*cos(r)**4)/96. + (6.*cos(r)**5)/25. &
      + cos(r)**6/36. + (338.*r*sin(r))/15. + (231.*r**2.)/32. + (215.*r*cos(r)*sin(r))/16. &
      + (124.*r*cos(r)**2*sin(r))/15. + (95.*r*cos(r)**3*sin(r))/24. + (6.*r*cos(r)**4*sin(r))/5. + (r*cos(r)**5*sin(r))/6.

  END FUNCTION



  REAL FUNCTION hfunction4(r)
     REAL, INTENT(IN) :: r

     hfunction4=(552.*cos(r))/7. + (6307.*cos(r)**2)/256. + (248.*cos(r)**3)/21. + (1505.*cos(r)**4)/256. + (88.*cos(r)**5)/35. &
      + (77.*cos(r)**6)/96. + (8.*cos(r)**7)/49. + cos(r)**8/64. - (6435.*pi**2)/256. + (552.*r*sin(r))/7. + (6435.*r**2)/256. &
      + (6307.*r*cos(r)*sin(r))/128. + (248.*r*cos(r)**2*sin(r))/7. + (1505.*r*cos(r)**3*sin(r))/64. + (88.*r*cos(r)**4*sin(r))/7. &
      + (77.*r*cos(r)**5*sin(r))/16. + (8.*r*cos(r)**6*sin(r))/7. + (r*cos(r)**7*sin(r))/8. + 45578./735.

  END FUNCTION

  REAL FUNCTION hSmoothAuxiliary(x)
     REAL, INTENT(IN) :: x

     hSmoothAuxiliary=1.-0.1*exp(-1./atan(1.-x)**3.)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary(x)
     REAL, INTENT(IN) :: x

     hDerivSmoothAuxiliary=3.*0.1*exp(1./atan(x - 1.)**3.)/(atan(x - 1.)**4*((x - 1.)**2 + 1.))

  END FUNCTION

  REAL FUNCTION hSmoothAux(r2,power)
     REAL, INTENT(IN)   :: r2
     INTEGER, INTENT(IN):: power

     hSmoothAux=exp(-1./(1.-r2)**power)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAux(r2,power)
     REAL, INTENT(IN)   :: r2
     INTEGER, INTENT(IN):: power

     hDerivSmoothAux=power*exp(-1./(1.-r2)**power)*(1./(1.-r2)**(power+1))

  END FUNCTION


  REAL FUNCTION hSmoothTanAux(r2,power)
     REAL, INTENT(IN)   :: r2
     INTEGER, INTENT(IN):: power

     hSmoothTanAux=exp(-1./atan(1.-r2)**power)

  END FUNCTION

  REAL FUNCTION hDerivSmoothTanAux(r2,power)
     REAL, INTENT(IN)   :: r2
     INTEGER, INTENT(IN):: power

     hDerivSmoothTanAux=power*exp(-1./atan(1.-r2)**power)/(atan(1.-r2)**(power+1)*(1.+(1.-r2)**2))

  END FUNCTION


  REAL FUNCTION hSmoothAuxiliary1(r2)
     REAL, INTENT(IN) :: r2

     hSmoothAuxiliary1=exp(-1./(1.-r2)**2.)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary1(r2)
     REAL, INTENT(IN) :: r2

     hDerivSmoothAuxiliary1=-(2.*exp(-1./(r2 - 1.)**2.))/(r2 - 1.)**3

  END FUNCTION

  REAL FUNCTION hSmoothAuxiliary2(r2)
     REAL, INTENT(IN) :: r2

     hSmoothAuxiliary2=exp(-1./atan(1.-r2)**2)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary2(r2)
     REAL, INTENT(IN) :: r2

     hDerivSmoothAuxiliary2=-(2.*exp(-1./atan(r2 - 1.)**2))/(atan(r2 - 1.)**3*((r2 - 1.)**2 + 1.))

  END FUNCTION


  REAL FUNCTION hSmoothAuxiliary3(r2)
     REAL, INTENT(IN) :: r2

     hSmoothAuxiliary3=exp(-1./atan(1.-r2)**3)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary3(r2)
     REAL, INTENT(IN) :: r2

     hDerivSmoothAuxiliary3=(3*exp(1./atan(r2 - 1.)**3))/(atan(r2 - 1.)**4*((r2 - 1.)**2 + 1.))

  END FUNCTION


  REAL FUNCTION hSmoothAuxiliary4(r2)
     REAL, INTENT(IN) :: r2

     hSmoothAuxiliary4=exp(-1./atan(1.-r2)**4)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary4(r2)
     REAL, INTENT(IN) :: r2

     hDerivSmoothAuxiliary4=-(4*exp(-1./atan(r2 - 1.)**4))/(atan(r2 - 1.)**5*((r2 - 1.)**2 + 1.))

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
  
  CASE(321) !*LAKE AT REST
  
    Prim(1) = 1. - 0.1 * SIN(2.*PI*x(1)) * COS(2.*PI*x(2))
    Prim(2) = 0.
    Prim(3) = 0.

    CALL PrimToCons(Prim,Cons)

  CASE(101,102) !Island at rest
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
  CASE (1)
   Bathymetry = 0.1 * SIN(2.*PI*X(1)) * COS(2.*PI*X(2)) 
  CASE(2)
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry = EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry = 0.0
    ENDIF
  CASE(3)
      Bathymetry = 2.-SIN(2.*PI*X(1)) - COS(2.*PI*X(2)) 
  CASE(4)
      Bathymetry = 0.
      IF (x(1) .GE. 13. .AND. x(1) .LE. 14) THEN
        Bathymetry = x(1)-13.
      ELSE IF (x(1) .GE. 14. .AND. x(1) .LE. 16) THEN
        Bathymetry = 1.
      ELSE IF (x(1) .GE. 16. .AND. x(1) .LE. 17) THEN
        Bathymetry = 1.-(x(1)-16.)
      ELSE IF (x(1) .GE. 23. .AND. x(1) .LE. 24) THEN
        Bathymetry = x(1)-23.
      ELSE IF (x(1) .GE. 24. .AND. x(1) .LE. 26) THEN
        Bathymetry = 1.
      ELSE IF (x(1) .GE. 26. .AND. x(1) .LE. 27) THEN
        Bathymetry = 1.-(x(1)-26.)
      END IF


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
  CASE (1)
   Bathymetry_X = 0.1 * 2. * PI * COS(2.*PI*X(1)) * COS(2.*PI*X(2)) 
  CASE(2)
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry_X = -2.0*X(1)/(1.0-r2)**2.*EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry_X = 0.0
    ENDIF
  CASE(3)
      Bathymetry_X = -2.*PI*COS(2.*PI*X(1))
  CASE(4)
      Bathymetry_X = 0.
      IF (x(1) .GE. 13. .AND. x(1) .LE. 14) THEN
        Bathymetry_X = 1.
      ELSE IF (x(1) .GE. 14. .AND. x(1) .LE. 16) THEN
        Bathymetry_X = 0.
      ELSE IF (x(1) .GE. 16. .AND. x(1) .LE. 17) THEN
        Bathymetry_X = -1.
      ELSE IF (x(1) .GE. 23. .AND. x(1) .LE. 24) THEN
        Bathymetry_X = 1.
      ELSE IF (x(1) .GE. 24. .AND. x(1) .LE. 26) THEN
        Bathymetry_X = 0.
      ELSE IF (x(1) .GE. 26. .AND. x(1) .LE. 27) THEN
        Bathymetry_X = -1.
      END IF
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
  CASE (1)
   Bathymetry_Y = - 0.1 * 2. * PI * SIN(2.*PI*X(1)) * SIN(2.*PI*X(2))
  CASE(2)
    r2=X(1)**2.+X(2)**2.
    IF (r2<1) THEN
      Bathymetry_Y = -2.0*X(2)/(1.0-r2)**2.*EXP(1.0-1.0/(1.0-r2))
    ELSE
      Bathymetry_Y = 0.0
    ENDIF
  CASE(3)
      Bathymetry_Y = 2.*PI*SIN(2.*PI*X(2)) 
  CASE(4)
      Bathymetry_Y = 0.
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
