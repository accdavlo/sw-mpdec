!===============================================================================!
PROGRAM FiniteVolume2D
!-------------------------------------------------------------------------------!
USE MOD_Mesh,              ONLY: BuildMesh
USE MOD_Parameters,        ONLY: InitializeParameters
USE MOD_FiniteVolume2D,    ONLY: FillInitialConditions
USE MOD_FiniteVolume2D,    ONLY: InitializeWBVariables
USE MOD_FiniteVolume2D,    ONLY: InitializeFiniteVolume
USE MOD_FiniteVolume2D,    ONLY: FinalizeFiniteVolume
USE MOD_TimeDiscretization,ONLY: TimeDiscretization
USE MOD_PostProcessing    ,ONLY: PostProcessing   

 
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!

CALL InitializeParameters()
CALL InitializeFiniteVolume()
CALL BuildMesh()

#ifdef WELLBALANCED
CALL InitializeWBVariables()
#endif

CALL FillInitialConditions()
CALL TimeDiscretization()
CALL PostProcessing()
CALL FinalizeFiniteVolume()

!-------------------------------------------------------------------------------!
END PROGRAM FiniteVolume2D
!===============================================================================!
