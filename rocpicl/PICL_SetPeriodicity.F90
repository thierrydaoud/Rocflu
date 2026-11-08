! ******************************************************************************
!
! Purpose: Particle periodicity from the periodic patches of the fluid
!   (<case>.bc, transforms built by rflupart).
!
! Description: Each BC_PERIODIC patch is a translation (linear periodicity,
!   along x, y or z) or a rotation about a coordinate axis (rotational
!   periodicity). The directions, axis and sector angle are reduced over all
!   ranks, because a rank holds only the patches of its own region. For a
!   rotation, each of the two planes gets its outward unit normal and the
!   rotation that takes it onto the related patch (the patch transform tm,
!   applied as MATMUL(tm,x)): a particle past a plane re-enters through the
!   other one.
!
! Input:
!   pRegion      Pointer to region
!   domainRange  Fluid domain extent (min, max in x, y, z) over all ranks
!
! Output:
!   planeNormal    Outward unit normal of each rotational periodic plane
!   planeRotation  Rotation taking each plane onto the other one
!   planeRadial    Unit vector in each plane, perpendicular to the axis,
!                  pointing away from it
!   axisVector     Unit vector of the rotation axis
!   global%piclPeriodicX/Y/ZFlag, global%piclRotPeriodicAxis,
!   global%piclRotPeriodicAngle
!
! ******************************************************************************

SUBROUTINE PICL_SetPeriodicity(pRegion,domainRange,planeNormal,planeRotation, &
                               planeRadial,axisVector)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

  TYPE(t_region), POINTER :: pRegion
  REAL(RFREAL), INTENT(IN) :: domainRange(2,3)
  REAL(RFREAL), INTENT(OUT) :: planeNormal(3,2),planeRotation(3,3,2), &
                               planeRadial(3,2),axisVector(3)

  INTEGER :: axisLocal(2),axisAll(2),dirLocal(3),dirAll(3),errorFlag,ifc, &
             iDir,iPatch,iPlane,l
  REAL(RFREAL) :: angle,angleLocal(2),angleAll(2),area,inside,insideAll, &
                  length,normal(3),shift(3),shiftLen,center(3),axisLocal3(3)
  REAL(RFREAL) :: planeLocal(15,2),planeAll(15,2)
  REAL(RFREAL), PARAMETER :: TOL = 1.0E-9_RFREAL
  CHARACTER(80) :: msg
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

  global => pRegion%global

  CALL RegisterFunction(global,'PICL_SetPeriodicity',__FILE__)

  length = MAXVAL(domainRange(2,:) - domainRange(1,:))

! Classify the local periodic patches. axisLocal and angleLocal hold
! (max, -min) so that one MPI_MAX gives both extremes. planeLocal holds,
! for the patch with the lower (1) and higher (2) global number of the
! rotational pair, its normal, rotation (column by column) and radial
! direction.
  dirLocal   = 0
  axisLocal  = (/ 0, -4 /)
  angleLocal = (/ -1.0_RFREAL, -HUGE(1.0_RFREAL) /)
  planeLocal = -HUGE(1.0_RFREAL)

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)
    IF ( pPatch%bcType /= BC_PERIODIC ) CYCLE

    IF ( MAXVAL(ABS(pPatch%tm(XCOORD:ZCOORD,XCOORD:ZCOORD) &
                    - IdentityMatrix())) < TOL ) THEN
      shift    = pPatch%tm(XCOORD:ZCOORD,XYZMAG)
      shiftLen = SQRT(SUM(shift**2))
      iDir     = MAXLOC(ABS(shift),1)
      IF ( SUM(ABS(shift)) - ABS(shift(iDir)) > TOL*shiftLen ) THEN
        WRITE(msg,'(A,I0,A)') 'periodic patch ',pPatch%iPatchGlobal, &
              ': translation not along x, y or z'
        CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(msg))
      END IF ! SUM
      IF ( ABS(shiftLen - (domainRange(2,iDir) - domainRange(1,iDir))) > &
           TOL*shiftLen ) THEN
        WRITE(msg,'(A,I0,A)') 'periodic patch ',pPatch%iPatchGlobal, &
              ': translation differs from the domain length'
        CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(msg))
      END IF ! ABS
      dirLocal(iDir) = 1
    ELSE
      angle = ABS(pPatch%angleRelated)
      axisLocal  = (/ MAX(axisLocal(1),pPatch%axisRelated), &
                      MAX(axisLocal(2),-pPatch%axisRelated) /)
      angleLocal = (/ MAX(angleLocal(1),angle), MAX(angleLocal(2),-angle) /)
      IF ( pPatch%nBFaces == 0 ) CYCLE

! --- Outward normal: area-weighted normal of the faces whose cell lies
!     behind them (a region can also hold faces of the patch rotated onto
!     the related plane, with inward normals). The patch must be flat and
!     its plane must contain the origin (the rotation axis)
      normal = 0.0_RFREAL
      center = 0.0_RFREAL
      area   = 0.0_RFREAL
      DO ifc = 1,pPatch%nBFaces
        IF ( .NOT. Outward(ifc) ) CYCLE
        normal = normal + pPatch%fn(XCOORD:ZCOORD,ifc)*pPatch%fn(XYZMAG,ifc)
        center = center + pPatch%fc(XCOORD:ZCOORD,ifc)*pPatch%fn(XYZMAG,ifc)
        area   = area + pPatch%fn(XYZMAG,ifc)
      END DO ! ifc
      IF ( area == 0.0_RFREAL ) CYCLE
      normal = normal/SQRT(SUM(normal**2))
      DO ifc = 1,pPatch%nBFaces
        IF ( .NOT. Outward(ifc) ) CYCLE
        IF ( MAXVAL(ABS(pPatch%fn(XCOORD:ZCOORD,ifc) - normal)) > 1.0E-6_RFREAL &
             .OR. ABS(DOT_PRODUCT(pPatch%fc(XCOORD:ZCOORD,ifc),normal)) > &
                  1.0E-6_RFREAL*length ) THEN
          WRITE(msg,'(A,I0,A)') 'rotational periodic patch ', &
                pPatch%iPatchGlobal,': not a plane through the origin'
          CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(msg))
        END IF ! MAXVAL
      END DO ! ifc

      iPlane = 1
      IF ( pPatch%iPatchGlobal > pPatch%iPatchRelated ) iPlane = 2
      planeLocal(1:3,iPlane)  = normal
      planeLocal(4:12,iPlane) = RESHAPE(pPatch%tm(XCOORD:ZCOORD,XCOORD:ZCOORD), &
                                        (/ 9 /))
      axisLocal3 = 0.0_RFREAL
      axisLocal3(pPatch%axisRelated) = 1.0_RFREAL
      center = center/area
      center = center - DOT_PRODUCT(center,axisLocal3)*axisLocal3
      planeLocal(13:15,iPlane) = center/SQRT(SUM(center**2))
    END IF ! MAXVAL
  END DO ! iPatch

  CALL MPI_Allreduce(dirLocal,dirAll,3,MPI_INTEGER,MPI_MAX,global%mpiComm, &
                     errorFlag)
  CALL MPI_Allreduce(axisLocal,axisAll,2,MPI_INTEGER,MPI_MAX,global%mpiComm, &
                     errorFlag)
  CALL MPI_Allreduce(angleLocal,angleAll,2,MPI_RFREAL,MPI_MAX, &
                     global%mpiComm,errorFlag)
  CALL MPI_Allreduce(planeLocal,planeAll,30,MPI_RFREAL,MPI_MAX, &
                     global%mpiComm,errorFlag)

  global%piclPeriodicXFlag = dirAll(1)
  global%piclPeriodicYFlag = dirAll(2)
  global%piclPeriodicZFlag = dirAll(3)
  global%piclRotPeriodicAxis  = 0
  global%piclRotPeriodicAngle = 0.0_RFREAL
  planeNormal   = 0.0_RFREAL
  planeRotation = 0.0_RFREAL
  planeRadial   = 0.0_RFREAL
  axisVector    = 0.0_RFREAL

  IF ( axisAll(1) > 0 ) THEN
    IF ( axisAll(1) /= -axisAll(2) .OR. &
         angleAll(1) + angleAll(2) > TOL*angleAll(1) ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
           'rotational periodic patches with different axes or angles')
    END IF ! axisAll
    DO l = 1,3
      IF ( l /= axisAll(1) .AND. dirAll(l) == 1 ) THEN
        CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
             'linear periodicity in the plane of the rotational periodicity')
      END IF ! l
    END DO ! l
    IF ( MINVAL(planeAll) < -0.5_RFREAL*HUGE(1.0_RFREAL) ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
           'rotational periodicity needs exactly one pair of periodic patches')
    END IF ! MINVAL
    planeNormal   = planeAll(1:3,:)
    planeRotation = RESHAPE(planeAll(4:12,:),(/ 3,3,2 /))
    planeRadial   = planeAll(13:15,:)
    axisVector(axisAll(1)) = 1.0_RFREAL

! - Every cell centroid inside both planes: catches an inward normal and
!   sectors wider than 180 degrees
    inside = -HUGE(1.0_RFREAL)
    DO iPlane = 1,2
      inside = MAX(inside,MAXVAL(MATMUL(planeNormal(:,iPlane), &
               pRegion%grid%cofg(XCOORD:ZCOORD,1:pRegion%grid%nCells))))
    END DO ! iPlane
    CALL MPI_Allreduce(inside,insideAll,1,MPI_RFREAL,MPI_MAX, &
                       global%mpiComm,errorFlag)
    IF ( insideAll > 0.0_RFREAL ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
           'cells outside a rotational periodic plane (sector wider than 180 deg?)')
    END IF ! insideAll

    global%piclRotPeriodicAxis  = axisAll(1)
    global%piclRotPeriodicAngle = angleAll(1)
  END IF ! axisAll

  CALL DeregisterFunction(global)

CONTAINS

  LOGICAL FUNCTION Outward(ifc)
    INTEGER, INTENT(IN) :: ifc
    Outward = DOT_PRODUCT(pRegion%grid%cofg(XCOORD:ZCOORD,pPatch%bf2c(ifc)) &
                          - pPatch%fc(XCOORD:ZCOORD,ifc), &
                          pPatch%fn(XCOORD:ZCOORD,ifc)) < 0.0_RFREAL
  END FUNCTION Outward

  FUNCTION IdentityMatrix() RESULT(eye)
    REAL(RFREAL) :: eye(3,3)
    eye = 0.0_RFREAL
    eye(1,1) = 1.0_RFREAL
    eye(2,2) = 1.0_RFREAL
    eye(3,3) = 1.0_RFREAL
  END FUNCTION IdentityMatrix

END SUBROUTINE PICL_SetPeriodicity
