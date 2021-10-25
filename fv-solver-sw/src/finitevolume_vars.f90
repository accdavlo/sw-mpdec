!===============================================================================!
MODULE MOD_FiniteVolume2D_vars
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PUBLIC
!-------------------------------------------------------------------------------!
! >> GLOBAL VARIABLES                                                           !
!-------------------------------------------------------------------------------!
INTEGER,PARAMETER   :: nVar  = 3
INTEGER,PARAMETER   :: nDims = 2
REAL                :: MESH_SX(1:nDims)
REAL                :: MESH_X0(1:nDims)
REAL                :: MESH_X1(1:nDims)
REAL                :: MESH_DX(1:nDims)
INTEGER             :: nGPs
INTEGER             :: nGhosts
INTEGER             :: nElemsX 
INTEGER             :: nElemsY

#ifdef PATANKAR
INTEGER             :: NNZsparse
INTEGER             :: NGlobalRows
INTEGER,ALLOCATABLE :: RowStart(:)
INTEGER,ALLOCATABLE :: ColumnsVector(:)
REAL, ALLOCATABLE   :: ProductionSparse(:)
REAL, ALLOCATABLE   :: DestructionSparse(:)
REAL,ALLOCATABLE    :: ProdUp(:,:)
REAL,ALLOCATABLE    :: DestUp(:,:)
INTEGER, ALLOCATABLE   :: SparseIndexMat(:,:)
#endif

REAL,ALLOCATABLE    :: MeshNodes(:,:,:)
REAL,ALLOCATABLE    :: MeshBary(:,:,:)
REAL,ALLOCATABLE    :: MeshGP(:,:,:,:,:)
REAL,ALLOCATABLE    :: WeightsGP(:,:)
REAL,ALLOCATABLE    :: WeightsGPBnd(:)
REAL,ALLOCATABLE    :: NormVectX(:,:,:,:)
REAL,ALLOCATABLE    :: NormVectY(:,:,:,:)
REAL,ALLOCATABLE    :: TangVectX(:,:,:,:)
REAL,ALLOCATABLE    :: TangVectY(:,:,:,:)

REAL,ALLOCATABLE    :: U(:,:,:)
REAL,ALLOCATABLE    :: V(:,:,:)
REAL,ALLOCATABLE    :: S(:,:,:)
REAL,ALLOCATABLE    :: SWB(:,:,:)
REAL,ALLOCATABLE    :: Ut(:,:,:)
REAL,ALLOCATABLE    :: FX(:,:,:)
REAL,ALLOCATABLE    :: FY(:,:,:)
REAL,ALLOCATABLE    :: FXWB(:,:,:)
REAL,ALLOCATABLE    :: FYWB(:,:,:)
REAL,ALLOCATABLE    :: WM(:,:,:,:)
REAL,ALLOCATABLE    :: WP(:,:,:,:)
REAL,ALLOCATABLE    :: FluxX(:,:,:,:)
REAL,ALLOCATABLE    :: FluxY(:,:,:,:)
LOGICAL,ALLOCATABLE :: Ind(:,:,:)

REAL,ALLOCATABLE    :: UN0(:,:,:)
REAL,ALLOCATABLE    :: K0(:,:,:)
REAL,ALLOCATABLE    :: K1(:,:,:)
REAL,ALLOCATABLE    :: K2(:,:,:)
REAL,ALLOCATABLE    :: K3(:,:,:)
REAL,ALLOCATABLE    :: K4(:,:,:)
REAL,ALLOCATABLE    :: K5(:,:,:)
REAL,ALLOCATABLE    :: FUp(:,:,:,:)
REAL,ALLOCATABLE    :: Up(:,:,:,:)
REAL,ALLOCATABLE    :: Ua(:,:,:,:)

REAL,ALLOCATABLE    :: UtWB(:,:,:)
REAL,ALLOCATABLE    :: Bath(:,:)

INTEGER,PARAMETER   :: UNIT_FILE = 123
INTEGER             :: WhichOutput
INTEGER             :: nOutputFiles
INTEGER             :: InitialCondition
INTEGER             :: BoundaryConditionsType(4)
REAL                :: PrimRefState1(1:nVar)
REAL                :: PrimRefState2(1:nVar)
REAL                :: PrimRefState3(1:nVar)
REAL                :: PrimRefState4(1:nVar)

REAL                :: t
REAL                :: tGlobal
REAL                :: dt
REAL                :: dt_Analyze
REAL                :: CFL
REAL                :: tEnd
REAL                :: Gravity
REAL                :: LambdaMaxX
REAL                :: LambdaMaxY
INTEGER             :: BathymetryFlag

INTEGER             :: Reconstruction
INTEGER             :: ReconstructionFix
INTEGER             :: timescheme 
REAL,PARAMETER      :: wLobatto = 1./12.
REAL,PARAMETER      :: WENOEPS = 1.0E-24
INTEGER,PARAMETER   :: WENOEXP = 2.0


REAL,PARAMETER      :: PI           = ACOS(-1.0)
REAL,PARAMETER      :: EPS          = 1.0E-6
REAL,PARAMETER      :: ACCURACY     = 1.0E-30
REAL,PARAMETER      :: MIN_DEPTH    = 1.0E-6
REAL,PARAMETER      :: MIN_SPEED    = 0.0
REAL,PARAMETER      :: MIN_TIMESTEP = 1.0E-30


CHARACTER(LEN=255)  :: VarNameVisu(1:nVar+1)



!-------------------------------------------------------------------------------!
END MODULE MOD_FiniteVolume2D_vars
!===============================================================================!
