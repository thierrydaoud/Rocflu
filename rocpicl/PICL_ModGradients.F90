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
            PICL_ComputeGradCellsQuad, &
            PICL_GradientNeeds, &
            PICL_CONSTR_NONE

! Gradient weights of the quadratic fit, one set per region on this rank
  TYPE t_quadSet
    INTEGER :: iRegionGlobal
    INTEGER, ALLOCATABLE :: beg(:)
    REAL(RFREAL), ALLOCATABLE :: wts(:,:)
  END TYPE t_quadSet

  INTEGER, SAVE :: nQuadSets = 0
  TYPE(t_quadSet), ALLOCATABLE, SAVE :: quadSets(:)

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

! ******************************************************************************
!
! Purpose: Which gradients the particle models use, from the input flags.
!
! Description: Same conditions as the interpolated slots in
!   libpicl/user_files/ppiclf_user_Layout.f; change both together.
!
! Output:
!   needPres  grad p: pressure gradient, Parmar added mass, baroclinic term
!   needDens  grad(phi rho): added mass, viscous unsteady, vorticity terms
!   needVel   grad u, v, w: Du/Dt (added mass, viscous unsteady), vorticity
!   needRhog  grad rho_g: added mass, viscous unsteady, POST output
!
! ******************************************************************************

  SUBROUTINE PICL_GradientNeeds(global,needPres,needDens,needVel,needRhog)

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) :: needDens,needPres,needRhog,needVel
    TYPE(t_global), POINTER :: global

    needPres = (global%piclPgFlag == 1) .OR. (global%piclAmFlag == 1) .OR. &
               (global%piclCollisionFlag >= 2)
    needDens = (global%piclAmFlag >= 1) .OR. &
               (global%piclViscousUnsteady >= 1) .OR. &
               (global%piclCollisionFlag >= 2)
    needVel  = (global%piclAmFlag >= 1) .OR. &
               (global%piclViscousUnsteady >= 1) .OR. &
               (global%piclCollisionFlag >= 2)
    needRhog = (global%piclAmFlag >= 1) .OR. &
               (global%piclViscousUnsteady >= 1) .OR. &
               (global%postWriteGradRho .EQV. .TRUE.)

  END SUBROUTINE PICL_GradientNeeds

! ******************************************************************************
!
! Purpose: Least-squares cell gradients for the particle models.
!
! Description: Fit of the gradient plus the second derivatives that the stencil
!   resolves, weighted by 1/|d|^2. In a symmetric stencil the gradient equals the
!   linear fit; in one-sided stencils (walls, outflow, incomplete periodic
!   corners) the curvature no longer leaks into the gradient. No wall
!   constraints, WENO or limiter.
!
! Input:
!   iBegVar,iEndVar  variables of var to differentiate
!   iBegGrad         first gradient slot of grad
!
! Output:
!   grad(:,iBegGrad:iBegGrad+iEndVar-iBegVar,1:nCells)
!
! ******************************************************************************

  SUBROUTINE PICL_ComputeGradCellsQuad(pRegion,iBegVar,iEndVar,iBegGrad,var,grad)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iBegGrad,iBegVar,iEndVar
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

    INTEGER :: icg,icg2,iGrad,iSet,isl,iVar,k
    REAL(RFREAL) :: dVar
    TYPE(t_global), POINTER :: global

    global => pRegion%global

    CALL RegisterFunction(global,'PICL_ComputeGradCellsQuad',__FILE__)

    iSet = PICL_QuadWeightSet(pRegion)

    DO icg = 1,pRegion%grid%nCells
      k = quadSets(iSet)%beg(icg)

      DO iVar = iBegVar,iEndVar
        iGrad = iBegGrad + iVar - iBegVar
        grad(XCOORD:ZCOORD,iGrad,icg) = 0.0_RFREAL
      END DO ! iVar

      DO isl = 1,pRegion%grid%c2cs(icg)%nCellMembs
        icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)

        DO iVar = iBegVar,iEndVar
          iGrad = iBegGrad + iVar - iBegVar
          dVar  = var(iVar,icg2) - var(iVar,icg)
          grad(XCOORD:ZCOORD,iGrad,icg) = grad(XCOORD:ZCOORD,iGrad,icg) &
                                        + quadSets(iSet)%wts(XCOORD:ZCOORD,k+isl-1)*dVar
        END DO ! iVar
      END DO ! isl
    END DO ! icg

    CALL DeregisterFunction(global)

  END SUBROUTINE PICL_ComputeGradCellsQuad

! ******************************************************************************
!
! Purpose: Index of this region's weight set, built on first use.
!
! Description: Per cell, the rows of the fit that give the gradient. Columns in
!   order dx, dy, dz, dx^2/2, dy^2/2, dz^2/2, dx dy, dx dz, dy dz. A direction
!   in which the stencil extends less than FLAT_TOL of its size (round-off, as
!   in 2D) is left out; modified Gram-Schmidt drops a column that depends on
!   the earlier ones (a curvature the stencil cannot separate from the gradient).
!
! ******************************************************************************

  INTEGER FUNCTION PICL_QuadWeightSet(pRegion)

    IMPLICIT NONE

    TYPE(t_region), POINTER :: pRegion

    INTEGER, PARAMETER :: NCOLS = 9
    REAL(RFREAL), PARAMETER :: DROP_TOL = 1.0E-6_RFREAL
    REAL(RFREAL), PARAMETER :: FLAT_TOL = 1.0E-8_RFREAL

    INTEGER :: errorFlag,i,icg,icg2,iSet,isl,j,nm,nMembsTot
    LOGICAL :: keep(NCOLS)
    LOGICAL :: flat(3)
    REAL(RFREAL) :: colNorm,d(3),dmax,dot,ext(3),sw
    REAL(RFREAL), ALLOCATABLE :: a(:,:),q(:,:),r(:,:),rinv(:,:)
    TYPE(t_global), POINTER :: global
    TYPE(t_quadSet), ALLOCATABLE :: setsTemp(:)

    global => pRegion%global

    DO iSet = 1,nQuadSets
      IF ( quadSets(iSet)%iRegionGlobal == pRegion%iRegionGlobal ) THEN
        PICL_QuadWeightSet = iSet
        RETURN
      END IF ! quadSets
    END DO ! iSet

    IF ( nQuadSets > 0 ) THEN
      ALLOCATE(setsTemp(nQuadSets))
      setsTemp(1:nQuadSets) = quadSets(1:nQuadSets)
      DEALLOCATE(quadSets)
      ALLOCATE(quadSets(nQuadSets+1))
      quadSets(1:nQuadSets) = setsTemp(1:nQuadSets)
      DEALLOCATE(setsTemp)
    ELSE
      ALLOCATE(quadSets(1))
    END IF ! nQuadSets
    nQuadSets = nQuadSets + 1
    iSet = nQuadSets
    quadSets(iSet)%iRegionGlobal = pRegion%iRegionGlobal

    nMembsTot = 0
    DO icg = 1,pRegion%grid%nCells
      nMembsTot = nMembsTot + pRegion%grid%c2cs(icg)%nCellMembs
    END DO ! icg

    ALLOCATE(quadSets(iSet)%beg(pRegion%grid%nCells), &
             quadSets(iSet)%wts(3,nMembsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'quadSets')
    END IF ! global%error

    nMembsTot = 0
    DO icg = 1,pRegion%grid%nCells
      nm = pRegion%grid%c2cs(icg)%nCellMembs
      quadSets(iSet)%beg(icg) = nMembsTot + 1

      ALLOCATE(a(nm,NCOLS),q(nm,NCOLS),r(NCOLS,NCOLS),rinv(NCOLS,NCOLS))

! Directions in which the stencil has no extent
      ext  = 0.0_RFREAL
      dmax = 0.0_RFREAL
      DO isl = 1,nm
        icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
        d(1:3) = pRegion%grid%cofg(1:3,icg2) - pRegion%grid%cofg(1:3,icg)
        ext(1:3) = MAX(ext(1:3),ABS(d(1:3)))
        dmax = MAX(dmax,SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3)))
      END DO ! isl
      flat(1:3) = ext(1:3) < FLAT_TOL*dmax

! Weighted design matrix, sqrt(w) = 1/|d|
      DO isl = 1,nm
        icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
        d(1:3) = pRegion%grid%cofg(1:3,icg2) - pRegion%grid%cofg(1:3,icg)
        sw = 1.0_RFREAL/SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
        WHERE ( flat(1:3) ) d(1:3) = 0.0_RFREAL
        a(isl,1) = sw*d(1)
        a(isl,2) = sw*d(2)
        a(isl,3) = sw*d(3)
        a(isl,4) = sw*0.5_RFREAL*d(1)*d(1)
        a(isl,5) = sw*0.5_RFREAL*d(2)*d(2)
        a(isl,6) = sw*0.5_RFREAL*d(3)*d(3)
        a(isl,7) = sw*d(1)*d(2)
        a(isl,8) = sw*d(1)*d(3)
        a(isl,9) = sw*d(2)*d(3)
      END DO ! isl

! Modified Gram-Schmidt, dropping dependent columns: a = q r
      r = 0.0_RFREAL
      DO j = 1,NCOLS
        q(:,j)  = a(:,j)
        colNorm = SQRT(SUM(a(:,j)*a(:,j)))
        DO i = 1,j-1
          IF ( keep(i) .EQV. .TRUE. ) THEN
            dot    = SUM(q(:,i)*q(:,j))
            r(i,j) = dot
            q(:,j) = q(:,j) - dot*q(:,i)
          END IF ! keep
        END DO ! i
        r(j,j)  = SQRT(SUM(q(:,j)*q(:,j)))
        keep(j) = (colNorm > 0.0_RFREAL) .AND. (r(j,j) > DROP_TOL*colNorm)
        IF ( keep(j) .EQV. .TRUE. ) THEN
          q(:,j) = q(:,j)/r(j,j)
        ELSE
          q(:,j) = 0.0_RFREAL
          r(:,j) = 0.0_RFREAL
        END IF ! keep
      END DO ! j

! Rows 1-3 of r^-1 (kept columns only), back substitution
      rinv = 0.0_RFREAL
      DO j = NCOLS,1,-1
        IF ( keep(j) .EQV. .FALSE. ) CYCLE
        rinv(j,j) = 1.0_RFREAL/r(j,j)
        DO i = j-1,1,-1
          IF ( keep(i) .EQV. .FALSE. ) CYCLE
          rinv(i,j) = -SUM(r(i,i+1:j)*rinv(i+1:j,j))/r(i,i)
        END DO ! i
      END DO ! j

! Gradient weights: rows 1-3 of r^-1 q^T, times sqrt(w)
      DO isl = 1,nm
        icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
        d(1:3) = pRegion%grid%cofg(1:3,icg2) - pRegion%grid%cofg(1:3,icg)
        sw = 1.0_RFREAL/SQRT(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
        DO i = 1,3
          quadSets(iSet)%wts(i,nMembsTot+isl) = sw*SUM(rinv(i,1:NCOLS)*q(isl,1:NCOLS))
        END DO ! i
      END DO ! isl

      DEALLOCATE(a,q,r,rinv)
      nMembsTot = nMembsTot + nm
    END DO ! icg

    PICL_QuadWeightSet = iSet

  END FUNCTION PICL_QuadWeightSet

END MODULE PICL_ModGradients
