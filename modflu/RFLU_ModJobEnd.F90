! ******************************************************************************
!
! Purpose: Stop the time stepping in time for the end-of-run dump before the
!   Slurm job ends.
!
! Description: Rank 0 reads the end of the job from SLURM_JOB_END_TIME (Unix
!   time) and turns it into MPI_Wtime units. Each step it compares the time
!   left with a margin measured from the run: twice the slowest step plus the
!   slowest dump so far, plus 2 min; at least 10 min until a dump has been
!   timed. Without SLURM_JOB_END_TIME there is no such stop. The stop line
!   is printed at every verbosity: it is the only sign that the run needs a
!   restart.
!
! Notes: None.
!
! ******************************************************************************

MODULE RFLU_ModJobEnd

  USE ISO_C_BINDING, ONLY: C_LONG,C_PTR,C_NULL_PTR

  USE ModDataTypes
  USE ModParameters
  USE ModMPI
  USE ModGlobal, ONLY: t_global

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_JobEndInit, &
            RFLU_JobEndNear, &
            RFLU_JobEndDumpTime

  REAL(RFREAL), PARAMETER :: SLACK = 120.0_RFREAL, &
                             FLOOR_UNTIMED = 600.0_RFREAL

! Rank 0 only: the job end in MPI_Wtime units and the slowest step and dump
  LOGICAL, SAVE :: haveJobEnd = .FALSE., dumpTimed = .FALSE.
  REAL(RFREAL), SAVE :: jobEnd = 0.0_RFREAL, slowestStep = 0.0_RFREAL, &
                        slowestDump = 0.0_RFREAL

  INTERFACE
    FUNCTION c_time(t) BIND(C,NAME='time')
      IMPORT :: C_LONG,C_PTR
      TYPE(C_PTR), VALUE :: t
      INTEGER(C_LONG) :: c_time
    END FUNCTION c_time
  END INTERFACE

CONTAINS

! ******************************************************************************
!
! Purpose: Read the end of the Slurm job (rank 0).
!
! ******************************************************************************

  SUBROUTINE RFLU_JobEndInit(global)

    TYPE(t_global), POINTER :: global

    CHARACTER(32) :: text
    INTEGER :: length,status,errorFlag
    INTEGER(8) :: endUnix,nowUnix

    haveJobEnd = .FALSE.
    dumpTimed  = .FALSE.
    slowestStep = 0.0_RFREAL
    slowestDump = 0.0_RFREAL

    IF ( global%myProcid /= MASTERPROC ) RETURN

    CALL GET_ENVIRONMENT_VARIABLE('SLURM_JOB_END_TIME',text,length,status)
    IF ( status /= 0 .OR. length == 0 ) RETURN

    READ(text,*,IOSTAT=errorFlag) endUnix
    nowUnix = INT(c_time(C_NULL_PTR),8)
    IF ( errorFlag /= 0 .OR. endUnix <= nowUnix ) RETURN

    haveJobEnd = .TRUE.
    jobEnd = MPI_Wtime() + REAL(endUnix - nowUnix,RFREAL)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,4A)') SOLVER_NAME,'Slurm job ends in ', &
        FixedText(REAL(endUnix - nowUnix,RFREAL)/3600.0_RFREAL,1), &
        ' h: the run stops with a dump before it'
    END IF ! global%verbLevel

  END SUBROUTINE RFLU_JobEndInit

! ******************************************************************************
!
! Purpose: Whether to stop now, so that the end-of-run dump fits before the
!   job ends (collective).
!
! Input:
!   global      Global data
!   stepStart   MPI_Wtime at the start of this step
!
! ******************************************************************************

  LOGICAL FUNCTION RFLU_JobEndNear(global,stepStart)

    TYPE(t_global), POINTER :: global
    REAL(RFREAL), INTENT(IN) :: stepStart

    INTEGER :: errorFlag
    REAL(RFREAL) :: left,margin,now

    RFLU_JobEndNear = .FALSE.

    IF ( global%myProcid == MASTERPROC .AND. haveJobEnd ) THEN
      now = MPI_Wtime()
      slowestStep = MAX(slowestStep,now - stepStart)
      margin = 2.0_RFREAL*(slowestStep + slowestDump) + SLACK
      IF ( dumpTimed .EQV. .FALSE. ) margin = MAX(margin,FLOOR_UNTIMED)
      left = jobEnd - now

      IF ( left <= margin ) THEN
        RFLU_JobEndNear = .TRUE.
        WRITE(STDOUT,'(A,1X,9A)') SOLVER_NAME, &
          'Stopping for the Slurm job end: ',FixedText(left/60.0_RFREAL,1), &
          ' min left, margin ',FixedText(margin/60.0_RFREAL,1), &
          ' min (slowest step ',FixedText(slowestStep,2), &
          ' s, slowest dump ',FixedText(slowestDump,2),' s)'
      END IF ! left
    END IF ! global%myProcid

    CALL MPI_Bcast(RFLU_JobEndNear,1,MPI_LOGICAL,MASTERPROC, &
                   global%mpiComm,errorFlag)

  END FUNCTION RFLU_JobEndNear

! ******************************************************************************
!
! Purpose: Record the time of a dump, the slowest rank's (collective).
!
! ******************************************************************************

  SUBROUTINE RFLU_JobEndDumpTime(global,seconds)

    TYPE(t_global), POINTER :: global
    REAL(RFREAL), INTENT(IN) :: seconds

    INTEGER :: errorFlag
    REAL(RFREAL) :: slowest

    CALL MPI_Reduce(seconds,slowest,1,MPI_RFREAL,MPI_MAX,MASTERPROC, &
                    global%mpiComm,errorFlag)

    IF ( global%myProcid == MASTERPROC ) THEN
      slowestDump = MAX(slowestDump,slowest)
      dumpTimed = .TRUE.
    END IF ! global%myProcid

  END SUBROUTINE RFLU_JobEndDumpTime

! ******************************************************************************
!
! Purpose: A number with the given decimals, leading zero kept, no blanks.
!
! ******************************************************************************

  FUNCTION FixedText(value,decimals) RESULT(text)

    REAL(RFREAL), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: decimals
    CHARACTER(:), ALLOCATABLE :: text

    CHARACTER(32) :: buffer,form

    WRITE(form,'(A,I0,A)') '(F20.',decimals,')'
    WRITE(buffer,form) value
    text = TRIM(ADJUSTL(buffer))

  END FUNCTION FixedText

END MODULE RFLU_ModJobEnd
