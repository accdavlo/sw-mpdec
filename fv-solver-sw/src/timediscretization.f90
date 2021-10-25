!===============================================================================!
MODULE MOD_TimeDiscretization
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE TimeDiscretization
  MODULE PROCEDURE TimeDiscretization
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: TimeDiscretization
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
SUBROUTINE TimeDiscretization()
!-------------------------------------------------------------------------------!
USE MOD_Output,             ONLY: WriteMeshToDisk
USE MOD_Output,             ONLY: WriteSolutionToDisk
USE MOD_Equation,           ONLY: TimeStep
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: tEnd
USE MOD_FiniteVolume2D_vars,ONLY: dt_Analyze
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: tGlobal 
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL :: t, tAnalyze, dt_min
!-------------------------------------------------------------------------------!

t = 0.0
dt_Analyze = tEnd/REAL(nOutputFiles)
tAnalyze   = t + dt_Analyze
IF (WhichOutput .EQ. 1 .OR. WhichOutput .EQ. 3) THEN
  CALL WriteMeshToDisk()
END IF
CALL WriteSolutionToDisk(t)
IF (t .EQ. tEnd) THEN
  RETURN
END IF

DO
  dt_min = TimeStep()
  dt = MIN(MIN(dt_min,tAnalyze-t),MIN(dt_min,tEnd-t))
  
  SELECT CASE(timescheme)
    CASE(1)
      CALL TimeDiscretizationByForwardEuler(t)
    CASE(2)
      CALL TimeDiscretizationBySSPRK4(t)
    CASE(3)
      CALL TimeDiscretizationByRK65(t)
    CASE(4)
      CALL TimeDiscretizationByDeC5(t)
#ifdef PATANKAR
    CASE(5)
      CALL TimeDiscretizationByMPDeC5(t)
    CASE(6)
      CALL TimeDiscretizationByMPEuler(t)
    CASE(7)
      CALL TimeDiscretizationByMPDeC2(t)
#endif
    CASE default
      WRITE(*,*) "Time discretization not implemented"
      STOP
  END SELECT
  
  t = t + dt
  tGlobal = t
  IF (ABS(tAnalyze-t) .LT. 1.0E-10) THEN
    CALL WriteSolutionToDisk(t)
    tAnalyze = tAnalyze + dt_Analyze
    IF (tAnalyze .GT. tEnd) THEN
      tAnalyze = tEnd
    END IF
  END IF
  IF (ABS(t-tEnd) .LT. 1.0E-10) THEN
    EXIT
  END IF
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretization
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! First Stage        !
!--------------------!
tStage = t + 0.0*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = &
    1.00000000000000*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.39175222700392*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K1(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 0.39175222700392*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = &
    0.44437049406734*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.55562950593266*K1(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.36841059262959*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Third Stage        !
!--------------------!
tStage = t + 0.58607968896780*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = &
    0.62010185138540*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.37989814861460*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.25189177424738*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K3(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fourth Stage       !
!--------------------!
tStage = t + 0.474542364687*dt
CALL FVTimeDerivative(tStage)
K4(1:nVar,1:nElemsX,1:nElemsY) = &
    0.17807995410773*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.82192004589227*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.54497475021237*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K4(1:nVar,1:nElemsX,1:nElemsY)
K5(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fifth Stage        !
!--------------------!
tStage  = t + 0.93501063100924*dt
CALL FVTimeDerivative(tStage)
U(1:nVar,1:nElemsX,1:nElemsY)  = &
    0.00683325884039*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.51723167208978*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.12759831133288*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.34833675773694*K4(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.08460416338212*K5(1:nVar,1:nElemsX,1:nElemsY)*dt &
  + 0.22600748319395*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK4
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByForwardEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Forward Euler      !
!--------------------!
tStage = t 
CALL FVTimeDerivative(tStage)

U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByForwardEuler
!===============================================================================!
!
!!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByRK65(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: UN0!U^(0)
USE MOD_FiniteVolume2D_vars,ONLY: K0 !F(U^(0))
USE MOD_FiniteVolume2D_vars,ONLY: K1 !F(U^(1))
USE MOD_FiniteVolume2D_vars,ONLY: K2 !F(U^(2))
USE MOD_FiniteVolume2D_vars,ONLY: K3 !F(U^(3))
USE MOD_FiniteVolume2D_vars,ONLY: K4 !F(U^(4))
USE MOD_FiniteVolume2D_vars,ONLY: K5 !F(U^(5))
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(6,5) :: ARK
REAL, DIMENSION(6)   :: bRK, cRK
REAL                 :: tStage
INTEGER              :: ii, jj, nStages
!-------------------------------------------------------------------------------!
nStages = 6
ARK = reshape((/ 0., 0.25, 0.125, 0., 0.1875, -0.42857142857142855, &
                 0., 0.,   0.125, 0., -0.375, 1.1428571428571428,   &
                 0., 0.,   0.,   0.5,  0.375, 0.8571428571428571, &
                 0., 0.,   0.,   0. , 0.5625, -1.7142857142857142, &
                 0., 0.,   0.,   0. , 0.    , 1.1428571428571428 /), shape(ARK))
bRK = (/ 0.077777777777777778, 0., 0.35555555555555556, 0.13333333333333333, 0.35555555555555556, 0.077777777777777778/)
cRK = (/ 0.0, 0.25, 0.25, 0.5, 0.75, 1.0 /)

!--------------------!
! Zero  Stage        !
!--------------------!

UN0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
tStage = t + cRK(1)*dt
CALL FVTimeDerivative(tStage)
K0(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! First Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +ARK(2,1)*K0(1:nVar,1:nElemsX,1:nElemsY)*dt
tStage = t + cRK(2)*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Second Stage       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(3,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(3,2)*K1(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(3)*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Third Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(4,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(4,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(4,3)*K2(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(4)*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Fourth Stage       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(5,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(5,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(5,3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(5,4)*K3(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(5)*dt
CALL FVTimeDerivative(tStage)
K4(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Fifth Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(6,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(6,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,4)*K3(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,5)*K4(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(6)*dt
CALL FVTimeDerivative(tStage)
K5(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Final Update       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(bRK(1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + bRK(2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(4)*K3(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(5)*K4(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(6)*K5(1:nVar,1:nElemsX,1:nElemsY) )*dt
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByRK65
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC5(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr = 5, MSteps = 4
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC5
!===============================================================================!
!
!!
#ifdef PATANKAR

!===============================================================================!
SUBROUTINE TimeDiscretizationByMPDeC5(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr = 5, MSteps = 4
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPDeC5
!===============================================================================!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByMPDeC2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr = 2, MSteps = 2
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, 0.5000000000000000000, 0.0000000000000000000, 0.5000000000000000000 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPDeC2
!===============================================================================!
!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByMPEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr = 1, MSteps = 2
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, 0.5000000000000000000, 0.0000000000000000000, 0.5000000000000000000 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPEuler
!===============================================================================!
!
!
#ifdef PATANKAR
! JACOBI: mass is split into diagonal (vector) and off diagonal (CRS)
!===============================================================================!
SUBROUTINE MODIFIED_PATANKAR_MATRIX_INVERSION( m_substep, MSteps, thetaDeC, hp, ha)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars
USE MOD_Mesh,               ONLY: GlobalElem
USE MOD_JacobiIteration,    ONLY: jacobi
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: m_substep
INTEGER,INTENT(IN) :: MSteps
REAL, DIMENSION(MSteps,MSteps), INTENT(IN) :: thetaDeC
REAL, DIMENSION(MSteps,nElemsX,nElemsY), INTENT(IN) :: hp
REAL, DIMENSION(nElemsX,nElemsY), INTENT(INOUT) :: ha
REAL, DIMENSION(NNZsparse) :: Mass !OffDiagonal Matrix
REAL, DIMENSION(nElemsX*nElemsY) :: Diagonal
REAL, DIMENSION(nElemsX*nElemsY) :: rhs, haVec
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER              :: ii, rr, jj,K,L,KKSparseVect, KLSparseVect, iRowStart, iiLPer, jjLPer
!-------------------------------------------------------------------------------!

Mass = 0.
Diagonal = 1.
KLSparseVect=0
iRowStart = 0

DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)

    KKSparseVect = SparseIndexMat(K,3)
    
    L = GlobalElem(ii,jj-1)
    iiLPer = MODULO(ii-1,nElemsX)+1
    jjLPer = MODULO(jj-2,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,1)

    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO
    


    L = GlobalElem(ii-1,jj)           
    iiLPer = MODULO(ii-2,nElemsX)+1  
    jjLPer = MODULO(jj-1,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,2)

    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO


    
    L = GlobalElem(ii+1,jj)
    iiLPer = MODULO(ii,  nElemsX)+1
    jjLPer = MODULO(jj-1,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,4)


    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO



    L = GlobalElem(ii,jj+1)
    iiLPer = MODULO(ii-1,nElemsX)+1
    jjLPer = MODULO(jj,  nElemsY)+1

    KLSparseVect = SparseIndexMat(K,5)


    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO


    rhs(GlobalElem(ii,jj))=hp(1,ii,jj)
    
  END DO
END DO

haVec=0.


CALL jacobi(RowStart,ColumnsVector,Mass, Diagonal, rhs, haVec)


DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)
    ha(ii,jj)=haVec(K)
  END DO
END DO


!
!-------------------------------------------------------------------------------!
END SUBROUTINE MODIFIED_PATANKAR_MATRIX_INVERSION
!===============================================================================!
!
#endif
!===============================================================================!
REAL FUNCTION  PATANKAR_DIV( prod , h )
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH
IMPLICIT NONE 
REAL, INTENT(IN) :: prod, h

IF ( h .LT. MIN_DEPTH ) THEN
  PATANKAR_DIV = 0.
ELSE
  PATANKAR_DIV = 2. * h * prod / ( h*h + MAX( h*h , MIN_DEPTH ) )
ENDIF

END  FUNCTION  PATANKAR_DIV 
!===============================================================================!
#endif
!
!===============================================================================!
END MODULE MOD_TimeDiscretization
!-------------------------------------------------------------------------------!
