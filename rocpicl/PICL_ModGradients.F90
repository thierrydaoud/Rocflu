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
            PICL_CONSTR_NONE, &
            PICL_FIELD_SCALAR, &
            PICL_FIELD_VECTOR

! Gradient weights of the quadratic fit, one set per region on this rank.
! Cells with slip-wall faces (iWall > 0) also carry the weights of the fit
! with the stencil mirrored across those walls, from begWall(iWall): wtsEven
! per member for scalars; for a vector, tensors wtsVec(dir,a,b) per member,
! wtsVecSelf for the cell and wtsSum (the sum of all row weights)
  TYPE t_quadSet
    INTEGER :: iRegionGlobal
    INTEGER, ALLOCATABLE :: beg(:),begWall(:),iWall(:)
    REAL(RFREAL), ALLOCATABLE :: wts(:,:),wtsEven(:,:),wtsSum(:,:), &
                                 wtsVec(:,:,:,:),wtsVecSelf(:,:,:,:)
  END TYPE t_quadSet

  INTEGER, SAVE :: nQuadSets = 0
  TYPE(t_quadSet), ALLOCATABLE, SAVE :: quadSets(:)

! No wall constraint (not a V_MIXT_* code)
  INTEGER, PARAMETER :: PICL_CONSTR_NONE = 0

! Mirror rule at slip walls: scalar (even) or velocity vector (reflected)
  INTEGER, PARAMETER :: PICL_FIELD_SCALAR = 1, &
                        PICL_FIELD_VECTOR = 2

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
!   linear fit; in one-sided stencils (outflow, incomplete periodic corners) the
!   curvature no longer leaks into the gradient. At slip walls the stencil is
!   mirrored across the wall (as reflective ghost cells): a scalar keeps its
!   value in the image, a velocity is reflected, R u with R = I - 2 n n^T. No
!   WENO or limiter.
!
! Input:
!   iBegVar,iEndVar  variables of var to differentiate
!   iBegGrad         first gradient slot of grad
!   fieldKind        PICL_FIELD_SCALAR, or PICL_FIELD_VECTOR for the three
!                    velocity components iBegVar..iBegVar+2
!
! Output:
!   grad(:,iBegGrad:iBegGrad+iEndVar-iBegVar,1:nCells)
!
! ******************************************************************************

  SUBROUTINE PICL_ComputeGradCellsQuad(pRegion,iBegVar,iEndVar,iBegGrad,var, &
                                       grad,fieldKind)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: fieldKind,iBegGrad,iBegVar,iEndVar
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

    INTEGER :: ia,ib,icg,icg2,iGrad,iSet,isl,iVar,iw,k,kw
    REAL(RFREAL) :: dVar
    TYPE(t_global), POINTER :: global

    global => pRegion%global

    CALL RegisterFunction(global,'PICL_ComputeGradCellsQuad',__FILE__)

    IF ( fieldKind == PICL_FIELD_VECTOR .AND. iEndVar-iBegVar /= 2 ) THEN
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! fieldKind

    iSet = PICL_QuadWeightSet(pRegion)

    DO icg = 1,pRegion%grid%nCells
      k  = quadSets(iSet)%beg(icg)
      iw = quadSets(iSet)%iWall(icg)

      DO iVar = iBegVar,iEndVar
        iGrad = iBegGrad + iVar - iBegVar
        grad(XCOORD:ZCOORD,iGrad,icg) = 0.0_RFREAL
      END DO ! iVar

      IF ( iw == 0 ) THEN
        DO isl = 1,pRegion%grid%c2cs(icg)%nCellMembs
          icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)

          DO iVar = iBegVar,iEndVar
            iGrad = iBegGrad + iVar - iBegVar
            dVar  = var(iVar,icg2) - var(iVar,icg)
            grad(XCOORD:ZCOORD,iGrad,icg) = grad(XCOORD:ZCOORD,iGrad,icg) &
                                          + quadSets(iSet)%wts(XCOORD:ZCOORD,k+isl-1)*dVar
          END DO ! iVar
        END DO ! isl
      ELSE IF ( fieldKind == PICL_FIELD_SCALAR ) THEN
        kw = quadSets(iSet)%begWall(iw)
        DO isl = 1,pRegion%grid%c2cs(icg)%nCellMembs
          icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)

          DO iVar = iBegVar,iEndVar
            iGrad = iBegGrad + iVar - iBegVar
            dVar  = var(iVar,icg2) - var(iVar,icg)
            grad(XCOORD:ZCOORD,iGrad,icg) = grad(XCOORD:ZCOORD,iGrad,icg) &
                                          + quadSets(iSet)%wtsEven(XCOORD:ZCOORD,kw+isl-1)*dVar
          END DO ! iVar
        END DO ! isl
      ELSE
! grad u_a = sum over rows of w (R u_owner - u)_a
        DO ia = 1,3
          iGrad = iBegGrad + ia - 1
          grad(XCOORD:ZCOORD,iGrad,icg) = -quadSets(iSet)%wtsSum(XCOORD:ZCOORD,iw) &
                                          *var(iBegVar+ia-1,icg)
          DO ib = 1,3
            grad(XCOORD:ZCOORD,iGrad,icg) = grad(XCOORD:ZCOORD,iGrad,icg) &
              + quadSets(iSet)%wtsVecSelf(XCOORD:ZCOORD,ia,ib,iw)*var(iBegVar+ib-1,icg)
          END DO ! ib
        END DO ! ia

        DO isl = 1,pRegion%grid%c2cs(icg)%nCellMembs
          icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
          kw   = quadSets(iSet)%begWall(iw) + isl - 1

          DO ia = 1,3
            iGrad = iBegGrad + ia - 1
            DO ib = 1,3
              grad(XCOORD:ZCOORD,iGrad,icg) = grad(XCOORD:ZCOORD,iGrad,icg) &
                + quadSets(iSet)%wtsVec(XCOORD:ZCOORD,ia,ib,kw)*var(iBegVar+ib-1,icg2)
            END DO ! ib
          END DO ! ia
        END DO ! isl
      END IF ! iw
    END DO ! icg

    CALL DeregisterFunction(global)

  END SUBROUTINE PICL_ComputeGradCellsQuad

! ******************************************************************************
!
! Purpose: Gradient rows of the weighted quadratic fit for offsets d(1:3,1:nr).
!
! Description: Columns in order dx, dy, dz, dx^2/2, dy^2/2, dz^2/2, dx dy,
!   dx dz, dy dz, weights 1/|d|^2. A direction in which the offsets extend less
!   than FLAT_TOL of their size (round-off, as in 2D) is left out; modified
!   Gram-Schmidt drops a column that depends on the earlier ones (a curvature
!   the stencil cannot separate from the gradient).
!
! Output:
!   w(1:3,isl): gradient = sum over isl of w(:,isl)*(value at row isl - value
!               at the cell)
!
! ******************************************************************************

  SUBROUTINE PICL_QuadFit(nr,d,w)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nr
    REAL(RFREAL), INTENT(IN) :: d(3,nr)
    REAL(RFREAL), INTENT(OUT) :: w(3,nr)

    INTEGER, PARAMETER :: NCOLS = 9
    REAL(RFREAL), PARAMETER :: DROP_TOL = 1.0E-6_RFREAL
    REAL(RFREAL), PARAMETER :: FLAT_TOL = 1.0E-8_RFREAL

    INTEGER :: i,isl,j
    LOGICAL :: flat(3),keep(NCOLS)
    REAL(RFREAL) :: colNorm,dd(3),dmax,dot,ext(3),sw
    REAL(RFREAL) :: a(nr,NCOLS),q(nr,NCOLS),r(NCOLS,NCOLS),rinv(NCOLS,NCOLS)

! Directions in which the stencil has no extent
    ext  = 0.0_RFREAL
    dmax = 0.0_RFREAL
    DO isl = 1,nr
      dd(1:3) = d(1:3,isl)
      ext(1:3) = MAX(ext(1:3),ABS(dd(1:3)))
      dmax = MAX(dmax,SQRT(dd(1)*dd(1) + dd(2)*dd(2) + dd(3)*dd(3)))
    END DO ! isl
    flat(1:3) = ext(1:3) < FLAT_TOL*dmax

! Weighted design matrix, sqrt(w) = 1/|d|
    DO isl = 1,nr
      dd(1:3) = d(1:3,isl)
      sw = 1.0_RFREAL/SQRT(dd(1)*dd(1) + dd(2)*dd(2) + dd(3)*dd(3))
      WHERE ( flat(1:3) ) dd(1:3) = 0.0_RFREAL
      a(isl,1) = sw*dd(1)
      a(isl,2) = sw*dd(2)
      a(isl,3) = sw*dd(3)
      a(isl,4) = sw*0.5_RFREAL*dd(1)*dd(1)
      a(isl,5) = sw*0.5_RFREAL*dd(2)*dd(2)
      a(isl,6) = sw*0.5_RFREAL*dd(3)*dd(3)
      a(isl,7) = sw*dd(1)*dd(2)
      a(isl,8) = sw*dd(1)*dd(3)
      a(isl,9) = sw*dd(2)*dd(3)
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
    DO isl = 1,nr
      dd(1:3) = d(1:3,isl)
      sw = 1.0_RFREAL/SQRT(dd(1)*dd(1) + dd(2)*dd(2) + dd(3)*dd(3))
      DO i = 1,3
        w(i,isl) = sw*SUM(rinv(i,1:NCOLS)*q(isl,1:NCOLS))
      END DO ! i
    END DO ! isl

  END SUBROUTINE PICL_QuadFit

! ******************************************************************************
!
! Purpose: Index of this region's weight set, built on first use.
!
! Description: Per cell, the weights of the fit (PICL_QuadFit) over its stencil.
!   For a cell with slip-wall faces, also the fit over the stencil and the cell
!   itself imaged across every combination of its wall planes; each image row
!   is folded onto the cell it images (the image of the cell itself drops out
!   for scalars).
!
! ******************************************************************************

  INTEGER FUNCTION PICL_QuadWeightSet(pRegion)

    USE ModBndPatch, ONLY: t_patch

    IMPLICIT NONE

    TYPE(t_region), POINTER :: pRegion

    INTEGER, PARAMETER :: NWALLMAX = 3

    INTEGER :: errorFlag,ia,ib,ic,icg,icg2,ifl,iPatch,iRow,iSet,isl,iw,j, &
               nComb,nm,nMembsTot,nMembsWall,nr,nWallCells
    INTEGER, ALLOCATABLE :: nWall(:),owner(:)
    LOGICAL :: dupl
    REAL(RFREAL) :: hh(3,3),x(3)
    REAL(RFREAL), ALLOCATABLE :: d(:,:),refl(:,:,:),w(:,:),wallN(:,:,:), &
                                 wallX(:,:,:)
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch
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

! Slip-wall planes (point, unit normal) of each real cell
    ALLOCATE(nWall(pRegion%grid%nCells),wallN(3,NWALLMAX,pRegion%grid%nCells), &
             wallX(3,NWALLMAX,pRegion%grid%nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nWall')
    END IF ! global%error

    nWall = 0
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
      IF ( pPatch%bcType /= BC_SLIPWALL ) CYCLE

      DO ifl = 1,pPatch%nBFaces
        icg  = pPatch%bf2c(ifl)
        dupl = .FALSE.
        DO iw = 1,nWall(icg)
          IF ( SUM(wallN(1:3,iw,icg)*pPatch%fn(XCOORD:ZCOORD,ifl)) > 0.999999_RFREAL ) THEN
            dupl = .TRUE.
          END IF ! SUM
        END DO ! iw
        IF ( dupl .OR. nWall(icg) == NWALLMAX ) CYCLE
        nWall(icg) = nWall(icg) + 1
        wallN(1:3,nWall(icg),icg) = pPatch%fn(XCOORD:ZCOORD,ifl)
        wallX(1:3,nWall(icg),icg) = pPatch%fc(XCOORD:ZCOORD,ifl)
      END DO ! ifl
    END DO ! iPatch

    nMembsTot  = 0
    nMembsWall = 0
    nWallCells = 0
    DO icg = 1,pRegion%grid%nCells
      nMembsTot = nMembsTot + pRegion%grid%c2cs(icg)%nCellMembs
      IF ( nWall(icg) > 0 ) THEN
        nWallCells = nWallCells + 1
        nMembsWall = nMembsWall + pRegion%grid%c2cs(icg)%nCellMembs
      END IF ! nWall
    END DO ! icg

    ALLOCATE(quadSets(iSet)%beg(pRegion%grid%nCells), &
             quadSets(iSet)%iWall(pRegion%grid%nCells), &
             quadSets(iSet)%begWall(MAX(nWallCells,1)), &
             quadSets(iSet)%wts(3,nMembsTot), &
             quadSets(iSet)%wtsEven(3,MAX(nMembsWall,1)), &
             quadSets(iSet)%wtsVec(3,3,3,MAX(nMembsWall,1)), &
             quadSets(iSet)%wtsVecSelf(3,3,3,MAX(nWallCells,1)), &
             quadSets(iSet)%wtsSum(3,MAX(nWallCells,1)),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'quadSets')
    END IF ! global%error

    quadSets(iSet)%wtsEven    = 0.0_RFREAL
    quadSets(iSet)%wtsVec     = 0.0_RFREAL
    quadSets(iSet)%wtsVecSelf = 0.0_RFREAL
    quadSets(iSet)%wtsSum     = 0.0_RFREAL

    nMembsTot  = 0
    nMembsWall = 0
    nWallCells = 0
    DO icg = 1,pRegion%grid%nCells
      nm = pRegion%grid%c2cs(icg)%nCellMembs
      quadSets(iSet)%beg(icg)   = nMembsTot + 1
      quadSets(iSet)%iWall(icg) = 0

      ALLOCATE(d(3,nm),w(3,nm))
      DO isl = 1,nm
        icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
        d(1:3,isl) = pRegion%grid%cofg(1:3,icg2) - pRegion%grid%cofg(1:3,icg)
      END DO ! isl
      CALL PICL_QuadFit(nm,d,w)
      quadSets(iSet)%wts(1:3,nMembsTot+1:nMembsTot+nm) = w(1:3,1:nm)
      DEALLOCATE(d,w)

      IF ( nWall(icg) > 0 ) THEN
        nWallCells = nWallCells + 1
        quadSets(iSet)%iWall(icg) = nWallCells
        quadSets(iSet)%begWall(nWallCells) = nMembsWall + 1

! Rows: the members (owner isl, R = I), then the cell (owner 0) and the members
! imaged across each non-empty combination ic of the wall planes
        nComb = 2**nWall(icg) - 1
        nr    = nm + (nm+1)*nComb
        ALLOCATE(d(3,nr),w(3,nr),owner(nr),refl(3,3,nr),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'d')
        END IF ! global%error

        iRow = 0
        DO ic = 0,nComb
          DO isl = 0,nm
            IF ( ic == 0 .AND. isl == 0 ) CYCLE
            IF ( isl == 0 ) THEN
              icg2 = icg
            ELSE
              icg2 = pRegion%grid%c2cs(icg)%cellMembs(isl)
            END IF ! isl

            x(1:3) = pRegion%grid%cofg(1:3,icg2)
            hh = 0.0_RFREAL
            DO ia = 1,3
              hh(ia,ia) = 1.0_RFREAL
            END DO ! ia
            DO j = 1,nWall(icg)
              IF ( BTEST(ic,j-1) .EQV. .FALSE. ) CYCLE
              x(1:3) = x(1:3) - 2.0_RFREAL*wallN(1:3,j,icg) &
                       *SUM((x(1:3) - wallX(1:3,j,icg))*wallN(1:3,j,icg))
              DO ib = 1,3
                hh(1:3,ib) = hh(1:3,ib) - 2.0_RFREAL*wallN(1:3,j,icg) &
                             *SUM(wallN(1:3,j,icg)*hh(1:3,ib))
              END DO ! ib
            END DO ! j

            iRow = iRow + 1
            d(1:3,iRow)      = x(1:3) - pRegion%grid%cofg(1:3,icg)
            owner(iRow)      = isl
            refl(1:3,1:3,iRow) = hh(1:3,1:3)
          END DO ! isl
        END DO ! ic

        CALL PICL_QuadFit(nr,d,w)

        DO iRow = 1,nr
          quadSets(iSet)%wtsSum(1:3,nWallCells) = quadSets(iSet)%wtsSum(1:3,nWallCells) &
                                                + w(1:3,iRow)
          IF ( owner(iRow) > 0 ) THEN
            iw = nMembsWall + owner(iRow)
            quadSets(iSet)%wtsEven(1:3,iw) = quadSets(iSet)%wtsEven(1:3,iw) + w(1:3,iRow)
          END IF ! owner
          DO ia = 1,3
            DO ib = 1,3
              IF ( owner(iRow) > 0 ) THEN
                iw = nMembsWall + owner(iRow)
                quadSets(iSet)%wtsVec(1:3,ia,ib,iw) = quadSets(iSet)%wtsVec(1:3,ia,ib,iw) &
                                                    + w(1:3,iRow)*refl(ia,ib,iRow)
              ELSE
                quadSets(iSet)%wtsVecSelf(1:3,ia,ib,nWallCells) = &
                  quadSets(iSet)%wtsVecSelf(1:3,ia,ib,nWallCells) + w(1:3,iRow)*refl(ia,ib,iRow)
              END IF ! owner
            END DO ! ib
          END DO ! ia
        END DO ! iRow

        DEALLOCATE(d,w,owner,refl,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'d')
        END IF ! global%error

        nMembsWall = nMembsWall + nm
      END IF ! nWall

      nMembsTot = nMembsTot + nm
    END DO ! icg

    DEALLOCATE(nWall,wallN,wallX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nWall')
    END IF ! global%error

    PICL_QuadWeightSet = iSet

  END FUNCTION PICL_QuadWeightSet

END MODULE PICL_ModGradients
