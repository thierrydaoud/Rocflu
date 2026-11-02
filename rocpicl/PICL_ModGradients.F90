! ******************************************************************************
!
! Purpose: Cell gradients of the fields used by the particle models.
!
! Description: Least-squares gradient with the wall constraints of one
!   variable code (or none), then optional WENO-XYZ and positivity limiter.
!
! ******************************************************************************

MODULE PICL_ModGradients

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsWrapper
  USE RFLU_ModLimiters, ONLY: RFLU_LimitGradCellsSimple
  USE RFLU_ModWENO, ONLY: RFLU_WENOGradCellsXYZWrapper

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PICL_ComputeGradCells, &
            PICL_CONSTR_NONE

! No wall constraint (not a V_MIXT_* code)
  INTEGER, PARAMETER :: PICL_CONSTR_NONE = 0

  CONTAINS

! ******************************************************************************
!
! Purpose: Gradients of var(iBegVar:iEndVar,:) into grad(:,1:iEndVar-iBegVar+1,:).
!
! Input:
!   constrCode  V_MIXT_* code whose wall constraints apply, or PICL_CONSTR_NONE
!   smooth      WENO-XYZ on the gradients
!   positive    zero a cell's gradient where a face value would be <= 0
!
! ******************************************************************************

  SUBROUTINE PICL_ComputeGradCells(pRegion,iBegVar,iEndVar,var,grad, &
                                   constrCode,smooth,positive)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: constrCode,iBegVar,iEndVar
    LOGICAL, INTENT(IN) :: positive,smooth
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

    INTEGER :: errorFlag,nGrad
    INTEGER, DIMENSION(:), POINTER :: constrInfo,signInfo
    TYPE(t_global), POINTER :: global

    global => pRegion%global

    CALL RegisterFunction(global,'PICL_ComputeGradCells',__FILE__)

    nGrad = iEndVar - iBegVar + 1

    ALLOCATE(constrInfo(iBegVar:iEndVar),signInfo(iBegVar:iEndVar),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'constrInfo')
    END IF ! global%error

    constrInfo = constrCode
    signInfo   = VAR_INFO_POS

    CALL RFLU_ComputeGradCellsWrapper(pRegion,iBegVar,iEndVar,1,nGrad, &
                                      constrInfo,var,grad)

    IF ( smooth .EQV. .TRUE. ) THEN
      CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,nGrad,grad)
    END IF ! smooth

    IF ( positive .EQV. .TRUE. ) THEN
      CALL RFLU_LimitGradCellsSimple(pRegion,iBegVar,iEndVar,1,nGrad, &
                                     var,signInfo,grad)
    END IF ! positive

    DEALLOCATE(constrInfo,signInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'constrInfo')
    END IF ! global%error

    CALL DeregisterFunction(global)

  END SUBROUTINE PICL_ComputeGradCells

END MODULE PICL_ModGradients
