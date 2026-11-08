! ******************************************************************************
!
! Purpose: Particle periodicity from the periodic patches of the fluid
!   (<case>.bc, transforms built by rflupart).
!
! Description: Each BC_PERIODIC patch is a translation (linear periodicity,
!   along x, y or z) or a rotation about a coordinate axis (rotational
!   periodicity). The directions, axis and sector angle are reduced over all
!   ranks, because a rank holds only the patches of its own region.
!
! Input:
!   pRegion      Pointer to region
!   domainRange  Fluid domain extent (min, max in x, y, z) over all ranks
!
! Output: global%piclPeriodicX/Y/ZFlag, global%piclRotPeriodicAxis,
!   global%piclRotPeriodicAngle
!
! ******************************************************************************

SUBROUTINE PICL_SetPeriodicity(pRegion,domainRange)

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

  INTEGER :: axisLocal(2),axisAll(2),dirLocal(3),dirAll(3),errorFlag,iDir, &
             iPatch,l
  REAL(RFREAL) :: angle,angleLocal(2),angleAll(2),shift(3),shiftLen
  REAL(RFREAL), PARAMETER :: TOL = 1.0E-9_RFREAL
  CHARACTER(80) :: msg
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

  global => pRegion%global

  CALL RegisterFunction(global,'PICL_SetPeriodicity',__FILE__)

! Classify the local periodic patches. axisLocal and angleLocal hold
! (max, -min) so that one MPI_MAX gives both extremes.
  dirLocal   = 0
  axisLocal  = (/ 0, -4 /)
  angleLocal = (/ -1.0_RFREAL, -HUGE(1.0_RFREAL) /)

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
    END IF ! MAXVAL
  END DO ! iPatch

  CALL MPI_Allreduce(dirLocal,dirAll,3,MPI_INTEGER,MPI_MAX,global%mpiComm, &
                     errorFlag)
  CALL MPI_Allreduce(axisLocal,axisAll,2,MPI_INTEGER,MPI_MAX,global%mpiComm, &
                     errorFlag)
  CALL MPI_Allreduce(angleLocal,angleAll,2,MPI_RFREAL,MPI_MAX, &
                     global%mpiComm,errorFlag)

  global%piclPeriodicXFlag = dirAll(1)
  global%piclPeriodicYFlag = dirAll(2)
  global%piclPeriodicZFlag = dirAll(3)
  global%piclRotPeriodicAxis  = 0
  global%piclRotPeriodicAngle = 0.0_RFREAL

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
    global%piclRotPeriodicAxis  = axisAll(1)
    global%piclRotPeriodicAngle = angleAll(1)
  END IF ! axisAll

  CALL DeregisterFunction(global)

CONTAINS

  FUNCTION IdentityMatrix() RESULT(eye)
    REAL(RFREAL) :: eye(3,3)
    eye = 0.0_RFREAL
    eye(1,1) = 1.0_RFREAL
    eye(2,2) = 1.0_RFREAL
    eye(3,3) = 1.0_RFREAL
  END FUNCTION IdentityMatrix

END SUBROUTINE PICL_SetPeriodicity
