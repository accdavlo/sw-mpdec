!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
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
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: timescheme 
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
INTEGER            :: iarg, nargs  
CHARACTER(len=32)  :: arg
!-------------------------------------------------------------------------------!

InitialCondition = 1 

nargs = command_argument_COUNT()
IF (nargs > 0) THEN
   CALL get_command_ARGUMENT(1, arg)
   READ(arg, *) iarg
   InitialCondition = iarg
END IF

SELECT CASE(InitialCondition)
  CASE(1) !*UNSTEADY SMOOTH VORTEX
    TEnd    = 0.1
    Gravity = 9.81
    nElemsX = 512
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/3.0,3.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    BathymetryFlag = 0       

  CASE(2) !*LAKE AT REST
    TEnd    = 0.1
    Gravity = 9.81
    nElemsX = 100 
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    BathymetryFlag = 1     

  CASE(3) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST 
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 400
    nElemsY = 120
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2      

  CASE(4) ! *CIRCULAR DAM BREAK 1 
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/40.0,40.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 0      

  CASE(5) ! *CIRCULAR DAM BREAK 2 
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 200
    nElemsY = 200
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/50.0,50.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 0      

  CASE(6) ! *WAVE OVER DRY ISLAND
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 400
    nElemsY = 120
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2      

  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

nargs = command_argument_COUNT()
IF (nargs > 1) THEN
   CALL get_command_ARGUMENT(2, arg)
   READ(arg, *) iarg
   nElemsX = iarg
   nElemsY = nElemsX
END IF

PRINT*, "nElemsX = ", nElemsX, ", nElemsY = ", nElemsY 

CFL      = 1.0

Reconstruction    = 4 ! 1 first order FV, 2 MUSCL, 3 WENO3, 4 WENO5
ReconstructionFix = Reconstruction

timescheme = 5   ! 1 explicit euler,  2 SSPRK64,  3 RK65,  4 DeC5,  5 PatankarDeC5,  6 PatankarEuler,  7 mPDeC2

WhichOutput  = 2 ! 1 Octave, 2 Tecplot, 3 Both
nOutputFiles = 4

VarNameVisu(1) = "Depth"
VarNameVisu(2) = "VelocityX"
VarNameVisu(3) = "VelocityY"
VarNameVisu(4) = "Bathymetry"

PRINT*, "Test = ",InitialCondition
PRINT*, "Time Scheme = ", timescheme
PRINT*, "Reconstruction = ", Reconstruction
PRINT*, "ReconstructionFix = ", ReconstructionFix
#ifdef PATANKAR
  PRINT*, "Patankar scheme"
#else
  PRINT*, "no patankar"
#endif

#ifdef WELLBALANCED
  PRINT*, "Well balanced"
#endif



!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
