!===============================================================================!
MODULE MOD_FiniteVolume2D
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeFiniteVolume
  MODULE PROCEDURE InitializeFiniteVolume
END INTERFACE

INTERFACE FillInitialConditions
  MODULE PROCEDURE FillInitialConditions
END INTERFACE

INTERFACE InitializeWBVariables
  MODULE PROCEDURE InitializeWBVariables
END INTERFACE

INTERFACE FVTimeDerivative
  MODULE PROCEDURE FVTimeDerivative
END INTERFACE

INTERFACE FinalizeFiniteVolume
  MODULE PROCEDURE FinalizeFiniteVolume
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeFiniteVolume
PUBLIC :: FillInitialConditions
PUBLIC :: InitializeWBVariables
PUBLIC :: FVTimeDerivative
PUBLIC :: FinalizeFiniteVolume
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
SUBROUTINE InitializeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: UN0
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
USE MOD_FiniteVolume2D_vars,ONLY: Up
USE MOD_FiniteVolume2D_vars,ONLY: Ua
USE MOD_FiniteVolume2D_vars,ONLY: FUp
USE MOD_FiniteVolume2D_vars,ONLY: Bath
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: timescheme 


#ifdef PATANKAR 
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse
USE MOD_FiniteVolume2D_vars,ONLY: NGlobalRows
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp
USE MOD_FiniteVolume2D_vars,ONLY: DestUp 
USE MOD_FiniteVolume2D_vars,ONLY: SparseIndexMat
USE MOD_Mesh,               ONLY: GlobalElem
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iSparseVect, iRowStart, K, L
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1) ! NONE
    nGhosts = 1
    nGPs    = 1
  CASE(2) ! MUSCL
    nGhosts = 1
    nGPs    = 1
  CASE(3) ! WENO3
    nGhosts = 1
    nGPs    = 2
  CASE(4) ! WENO5
    nGhosts = 2
    nGPs    = 4
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(MeshNodes(1:nDims,       0:nElemsX,0:nElemsY))
ALLOCATE(MeshBary (1:nDims,       1:nElemsX,1:nElemsY))
ALLOCATE(MeshGP   (1:nDims,       1:nElemsX,1:nElemsY,1:nGPs,1:nGPs))
ALLOCATE(WeightsGP(1:nGPs, 1:nGPs))
ALLOCATE(WeightsGPBnd(1:nGPs))
ALLOCATE(NormVectX(1:nDims,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(TangVectX(1:nDims,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(NormVectY(1:nDims,1:nGPs,1:nElemsX,0:nElemsY))
ALLOCATE(TangVectY(1:nDims,1:nGPs,1:nElemsX,0:nElemsY))

ALLOCATE( U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE( V(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(Ut(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Bath(1:nElemsX,1:nElemsY))
ALLOCATE(WM(1:nVar,1:nGPs,0:nElemsX+1,0:nElemsY+1))
ALLOCATE(WP(1:nVar,1:nGPs,0:nElemsX+1,0:nElemsY+1))
ALLOCATE( S(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FX(1:nVar,0:nElemsX,1:nElemsY))
ALLOCATE(FY(1:nVar,1:nElemsX,0:nElemsY))
ALLOCATE(FluxX(1:nVar,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(FluxY(1:nVar,1:nGPs,1:nElemsX,0:nElemsY))
ALLOCATE(Ind(1:2,0:nElemsX+1,0:nElemsY+1))

ALLOCATE(UN0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K1(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K2(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K3(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K4(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K5(1:nVar,1:nElemsX,1:nElemsY))

ALLOCATE(Ua(1:4,1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Up(1:4,1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FUp(1:4,1:nVar,1:nElemsX,1:nElemsY))

#ifdef WELLBALANCED
ALLOCATE( UtWB(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE( SWB(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FXWB(1:nVar,0:nElemsX,1:nElemsY))
ALLOCATE(FYWB(1:nVar,1:nElemsX,0:nElemsY))
#endif

#ifdef PATANKAR
NNZsparse = 5*nElemsX*nElemsY
NGlobalRows = nElemsX*nElemsY
ALLOCATE(ProductionSparse(1:NNZsparse))
ALLOCATE(DestructionSparse(1:NNZsparse))
ALLOCATE(ColumnsVector(1:NNZsparse))
ALLOCATE(RowStart(1:NGlobalRows+1))
ALLOCATE(ProdUp(1:4,1:NNZsparse))
ALLOCATE(DestUp(1:4,1:NNZsparse))
ALLOCATE(SparseIndexMat(1:nElemsX*nElemsY,1:5))
#endif

U  = 0.0
V  = 0.0
S  = 0.0
Ut = 0.0
FX = 0.0
FY = 0.0
WM = 0.0
WP = 0.0
FluxX = 0.0
FluxY = 0.0

Ind = .FALSE.

UN0 = 0.0
K0 = 0.0
K1 = 0.0
K2 = 0.0
K3 = 0.0
K4 = 0.0
K5 = 0.0

Bath = 0.

Ua  = 0.
Up  = 0.
FUp = 0.

#ifdef WELLBALANCED
UtWB = 0.
SWB  = 0.0
FXWB = 0.0
FYWB = 0.0
#endif



#ifdef PATANKAR
iSparseVect = 0
iRowStart=0

DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)
    !K = (ii-1)*NelemsY + jj
    iSparseVect = iSparseVect + 1
    iRowStart = iRowStart +1

    RowStart(iRowStart) = iSparseVect
    !*PRINT*, "iRowStart", iRowStart, "RowStart(iRowStart)", RowStart(iRowStart)


    L = GlobalElem(ii,jj-1)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,1) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii-1,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,2) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,3) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii+1,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,4) = iSparseVect



    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj+1)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,5) = iSparseVect


  END DO
END DO

RowStart(nElemsX*nElemsY+1) = RowStart(nElemsX*nElemsY)+5 
!*PRINT*, "iRowStart", nElemsX*nElemsY+1, "RowStart(iRowStart)", RowStart(nElemsX*nElemsY+1)


#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FillInitialConditions()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: Bathymetry   
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Bath
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp
REAL, DIMENSION(nGPs, nGPs)       :: Btemp
INTEGER :: ii, jj, iGP, jGP
!-------------------------------------------------------------------------------!

U     = 0.
V     = 0.
Utemp = 0.
Btemp = 0.
DO jj=1,nElemsY
  DO ii=1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        ! compute cell average of conservative variables
        CALL ExactFunction(&
          InitialCondition,0.0,MeshGP(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))
        U(1:nVar, ii, jj) = U(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Utemp(1:nVar,iGP,jGP)

        ! compute cell average of bathymetry            
        Btemp(iGP,jGP) = Bathymetry(MeshGP(:,ii,jj,iGP,jGP))
        Bath(ii,jj)    = Bath(ii,jj) + WeightsGP(iGP,jGP) * Btemp(iGP,jGP)
      END DO
    END DO
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FillInitialConditions
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE InitializeWBVariables()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_Equation           ,ONLY: ExactFunctionWB
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp
INTEGER :: ii, jj, iGP, jGP
!-------------------------------------------------------------------------------!

! compute UWB
Utemp = 0.
DO jj=1,nElemsY
  DO ii=1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunctionWB(InitialCondition,MeshGP(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))
        U(1:nVar, ii, jj) = U(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Utemp(1:nVar,iGP,jGP)
      END DO
    END DO
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

CALL FVTimeDerivative(0.)

UtWB = Ut
FXWB = FX
FYWB = FY
SWB  = S


!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeWBVariables
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FVTimeDerivative(t)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Reconstruction, ONLY: ReconstructionX
USE MOD_Reconstruction, ONLY: ReconstructionY
USE MOD_Reconstruction, ONLY: ReconstructionFixX
USE MOD_Reconstruction, ONLY: ReconstructionFixY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

CALL BoundaryConditions(t)

CALL ShocksIndicatorX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL PositivityLimiterX()
CALL NumericalFluxFX()

CALL ShocksIndicatorY()
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL PositivityLimiterY()
CALL NumericalFluxFY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()

!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE UpdateTimeDerivative()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Ut, S
USE MOD_Mesh,               ONLY: GlobalElem
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB

#ifdef PATANKAR
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, K, L, iSparseVect, iRowStart
!-------------------------------------------------------------------------------!

!Production = 0.
!Destruction = 0.

#ifdef PATANKAR
  iSparseVect = 0
  iRowStart=0
#endif  

#ifdef WELLBALANCED
  FX=FX-FXWB
  FY=FY-FYWB
  S =S -SWB
#endif

DO jj=1,nElemsY
  DO ii=1,nElemsX
#ifdef PATANKAR   
    K = GlobalElem(ii,jj)
    
    iSparseVect = iSparseVect + 1
    iRowStart = iRowStart +1

    L = GlobalElem(ii,jj-1)
    ProductionSparse(iSparseVect) = max(FY(1,ii+0,jj-1),0.)/Mesh_DX(2)
    DestructionSparse(iSparseVect) = -min(FY(1,ii+0,jj-1),0.)/Mesh_DX(2)

    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii-1,jj)
    ProductionSparse(iSparseVect) =   max(FX(1,ii-1,jj+0),0.)/Mesh_DX(1)
    DestructionSparse(iSparseVect) = -min(FX(1,ii-1,jj+0),0.)/Mesh_DX(1)


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj)
    ProductionSparse(iSparseVect) = 0.
    DestructionSparse(iSparseVect) = 0.


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii+1,jj)
    ProductionSparse(iSparseVect) =  max(-FX(1,ii+0,jj+0),0.)/Mesh_DX(1)
    DestructionSparse(iSparseVect) = -min(-FX(1,ii+0,jj+0),0.)/Mesh_DX(1)

    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj+1)
    ProductionSparse(iSparseVect) = max(-FY(1,ii+0,jj+0),0.)/Mesh_DX(2)
    DestructionSparse(iSparseVect) = -min(-FY(1,ii+0,jj+0),0.)/Mesh_DX(2)

    Ut(2:nVar,ii,jj) = S(2:nVar,ii,jj) &
                     - (FX(2:nVar,ii+0,jj+0)-FX(2:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(2:nVar,ii+0,jj+0)-FY(2:nVar,ii+0,jj-1))/Mesh_DX(2)

#else
    Ut(1:nVar,ii,jj) = S(1:nVar,ii,jj) &
                     - (FX(1:nVar,ii+0,jj+0)-FX(1:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(1:nVar,ii+0,jj+0)-FY(1:nVar,ii+0,jj-1))/Mesh_DX(2)

#endif
  END DO !ii
END DO !jj


!-------------------------------------------------------------------------------!
END SUBROUTINE UpdateTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFX()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FX    = 0.0
FluxX = 0.0

DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii+0,jj),&
                       WM(1:nVar,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims,1:nGPs,ii,jj),&
                       TangVectX(1:nDims,1:nGPs,ii,jj),&
                       FluxX(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

DO iGP = 1,nGPs
  FX(:,:,:) = FX(:,:,:) + weightsGPBnd(iGP)*FluxX(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFY()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FY    = 0.0
FluxY = 0.0

DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii,jj+0),&
                       WM(1:nVar,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims,1:nGPs,ii,jj),&
                       TangVectY(1:nDims,1:nGPs,ii,jj),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

DO iGP = 1,nGPs
  FY(:,:,:) = FY(:,:,:) + weightsGPBnd(iGP)*FluxY(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PositivityLimiterX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: wLobatto    
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH   

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiX, theta, mmin 
!-------------------------------------------------------------------------------!


DO jj=1,nElemsY
  DO ii=0,nElemsX
     alpha = 0.  
     beta  = 0.
     DO jGP = 1,nGPs
        alpha = alpha + WeightsGPBnd(jGP)*WM(1,jGP,ii+0,jj) 
        beta  = beta  + WeightsGPBnd(jGP)*WP(1,jGP,ii+0,jj)
     END DO
     csiX = ( V(1,ii,jj) - wLobatto*alpha - wLobatto*beta )/( 1. - 2.*wLobatto ) 

     DO jGP = 1,nGPs
        mmin = MIN( csiX , WM(1,jGP,ii+0,jj), WP(1,jGP,ii+0,jj))

        IF ( V(1,ii,jj) .EQ. mmin ) THEN
          theta = 1. 
        ELSE
          theta = MIN( 1. , ABS( (V(1,ii,jj)-MIN_DEPTH)/(V(1,ii,jj)-mmin) ) )
        ENDIF

 
        WP(1,jGP,ii+0,jj) = V(1,ii,jj) + theta * ( WP(1,jGP,ii+0,jj) - V(1,ii,jj) ) 
        WM(1,jGP,ii+0,jj) = V(1,ii,jj) + theta * ( WM(1,jGP,ii+0,jj) - V(1,ii,jj) )
     END DO

  END DO
END DO



END SUBROUTINE PositivityLimiterX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PositivityLimiterY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: wLobatto    
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DEPTH   

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!


DO jj=0,nElemsY
  DO ii=1,nElemsX
     alpha = 0.  
     beta  = 0.
     DO jGP = 1,nGPs
        alpha = alpha + WeightsGPBnd(jGP)*WM(1,jGP,ii,jj+0) 
        beta  = beta  + WeightsGPBnd(jGP)*WP(1,jGP,ii,jj+0)
     END DO
     csiY = ( V(1,ii,jj) - wLobatto*alpha - wLobatto*beta )/( 1. - 2.*wLobatto ) 

     DO jGP = 1,nGPs
        mmin = MIN( csiY , WM(1,jGP,ii,jj+0), WP(1,jGP,ii,jj+0))

        IF ( V(1,ii,jj) .EQ. mmin ) THEN
          theta = 1. 
        ELSE
          theta = MIN( 1. , ABS( (V(1,ii,jj)-MIN_DEPTH)/(V(1,ii,jj)-mmin) ) )
        ENDIF
 
        WP(1,jGP,ii,jj+0) = V(1,ii,jj) + theta * ( WP(1,jGP,ii,jj+0) - V(1,ii,jj) ) 
        WM(1,jGP,ii,jj+0) = V(1,ii,jj) + theta * ( WM(1,jGP,ii,jj+0) - V(1,ii,jj) )
     END DO

  END DO
END DO


END SUBROUTINE PositivityLimiterY
!===============================================================================!
!
!
!
!===============================================================================!
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FinalizeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: UN0
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
USE MOD_FiniteVolume2D_vars,ONLY: Ua
USE MOD_FiniteVolume2D_vars,ONLY: Up
USE MOD_FiniteVolume2D_vars,ONLY: FUp
USE MOD_FiniteVolume2D_vars,ONLY: Bath
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd

#ifdef PATANKAR 
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp
USE MOD_FiniteVolume2D_vars,ONLY: DestUp
USE MOD_FiniteVolume2D_vars,ONLY: SparseIndexMat
#endif   


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

DEALLOCATE(MeshNodes)
DEALLOCATE(MeshBary)
DEALLOCATE(MeshGP)
DEALLOCATE(WeightsGP)
DEALLOCATE(WeightsGPBnd)
DEALLOCATE(NormVectX)
DEALLOCATE(TangVectX)
DEALLOCATE(NormVectY)
DEALLOCATE(TangVectY)
DEALLOCATE(U)
DEALLOCATE(Bath)
DEALLOCATE(V)
DEALLOCATE(Ut)
DEALLOCATE(WM)
DEALLOCATE(WP)
DEALLOCATE(S)
DEALLOCATE(FX)
DEALLOCATE(FY)
DEALLOCATE(FluxX)
DEALLOCATE(FluxY)
DEALLOCATE(UN0)
DEALLOCATE(K0)
DEALLOCATE(K1)
DEALLOCATE(K2)
DEALLOCATE(K3)
DEALLOCATE(K4)
DEALLOCATE(K5)
    
DEALLOCATE(Ua)
DEALLOCATE(Up)
DEALLOCATE(FUp)

#ifdef WELLBALANCED
DEALLOCATE(UtWB)
DEALLOCATE(FXWB)
DEALLOCATE(FYWB)
DEALLOCATE(SWB)
#endif

#ifdef PATANKAR
DEALLOCATE(ProductionSparse)
DEALLOCATE(DestructionSparse)
DEALLOCATE(ColumnsVector)
DEALLOCATE(RowStart)

DEALLOCATE(ProdUp)
DEALLOCATE(DestUp)
DEALLOCATE(SparseIndexMat)
#endif   

!-------------------------------------------------------------------------------!
END SUBROUTINE FinalizeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_FiniteVolume2D
!-------------------------------------------------------------------------------!
