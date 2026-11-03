! ******************************************************************************
!
! Purpose: Add div(phi_g rho_g R_sg) and div(phi_g Q_sg) to the residual. Their virtual
!   cells must already hold the values of the cells they mirror.
!
! ******************************************************************************

SUBROUTINE PICL_AddSubgridDivergence(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE PICL_ModGradients, ONLY: PICL_ComputeGradCells, PICL_CONSTR_NONE

  IMPLICIT NONE

  TYPE(t_region), POINTER :: pRegion

  INTEGER :: icg
  REAL(RFREAL) :: divQsg,divRsgX,divRsgY,divRsgZ,vol
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: gradQsg,gradRsg
  TYPE(t_global), POINTER :: global

  global => pRegion%global

  CALL RegisterFunction(global,'PICL_AddSubgridDivergence',__FILE__)

  gradRsg => pRegion%mixt%piclGradPhiRsg
  gradQsg => pRegion%mixt%piclGradPhiQsg

  ! WENO would flatten the gradients at fronts
  CALL PICL_ComputeGradCells(pRegion,1,9,pRegion%mixt%piclPhiRsg,gradRsg, &
                             PICL_CONSTR_NONE,.FALSE.,.FALSE.)
  CALL PICL_ComputeGradCells(pRegion,1,3,pRegion%mixt%piclPhiQsg,gradQsg, &
                             PICL_CONSTR_NONE,.FALSE.,.FALSE.)

  DO icg = 1,pRegion%grid%nCells
    vol = pRegion%grid%vol(icg)

    divRsgX = vol*(gradRsg(XCOORD,1,icg) + gradRsg(YCOORD,2,icg) + gradRsg(ZCOORD,3,icg))
    divRsgY = vol*(gradRsg(XCOORD,4,icg) + gradRsg(YCOORD,5,icg) + gradRsg(ZCOORD,6,icg))
    divRsgZ = vol*(gradRsg(XCOORD,7,icg) + gradRsg(YCOORD,8,icg) + gradRsg(ZCOORD,9,icg))
    divQsg  = vol*(gradQsg(XCOORD,1,icg) + gradQsg(YCOORD,2,icg) + gradQsg(ZCOORD,3,icg))

    IF ( divQsg /= divQsg ) THEN
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,'PPICLF:Broken QSG')
    END IF ! divQsg

    pRegion%mixt%DivPhiRsg(XCOORD,icg) = divRsgX
    pRegion%mixt%DivPhiRsg(YCOORD,icg) = divRsgY
    pRegion%mixt%DivPhiRsg(ZCOORD,icg) = divRsgZ

    pRegion%mixt%rhs(CV_MIXT_XMOM,icg) = pRegion%mixt%rhs(CV_MIXT_XMOM,icg) + divRsgX
    pRegion%mixt%rhs(CV_MIXT_YMOM,icg) = pRegion%mixt%rhs(CV_MIXT_YMOM,icg) + divRsgY
    pRegion%mixt%rhs(CV_MIXT_ZMOM,icg) = pRegion%mixt%rhs(CV_MIXT_ZMOM,icg) + divRsgZ
    pRegion%mixt%rhs(CV_MIXT_ENER,icg) = pRegion%mixt%rhs(CV_MIXT_ENER,icg) + divQsg
  END DO ! icg

  CALL DeregisterFunction(global)

END SUBROUTINE PICL_AddSubgridDivergence
