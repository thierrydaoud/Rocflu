! ******************************************************************************
!
! Purpose: Projected feedback force of the actual cells, as last projected by ppiclF.
!
! ******************************************************************************

SUBROUTINE PICL_FillFeedback(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters

!DEC$ NOFREEFORM
#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!DEC$ FREEFORM

  IMPLICIT NONE

  TYPE(t_region), POINTER :: pRegion

  INTEGER(KIND=4) :: icg

  DO icg = 1,pRegion%grid%nCells
    CALL ppiclf_solve_GetProFld(icg,PPICLF_P_JFX,pRegion%mixt%piclFeedback(1,icg))
    CALL ppiclf_solve_GetProFld(icg,PPICLF_P_JFY,pRegion%mixt%piclFeedback(2,icg))
    CALL ppiclf_solve_GetProFld(icg,PPICLF_P_JFZ,pRegion%mixt%piclFeedback(3,icg))
  END DO ! icg

END SUBROUTINE PICL_FillFeedback
