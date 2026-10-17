! ******************************************************************************
!
! Purpose: Raw dumps of the state in restart_data/: the fluid, one file per
!   region and step, restart_data/fluid_<step>_<region>.dat (8 and 5 digits),
!   and with particles one file per step for all ranks,
!   restart_data/particles_<step>.dat.
!
! Description: The dump holds the state as it is in memory, real cells only,
!   in double precision, without conversion. It is written with stream access
!   (raw bytes, no Fortran record markers, native byte order):
!
!     CHARACTER(8)  'RFLUDUMP'
!     INTEGER(4)    format version (1)
!     INTEGER(8)    step
!     REAL(8)       time
!     INTEGER(4)    region, nCells, nCv
!     INTEGER(4)    form of the gas state: 1 = pure (rho, rho*u, ..., rho*E),
!                   2 = volume-weighted ((1-phi)*rho, ...; PICL runs)
!     INTEGER(4)    1 if the PICL fields follow, else 0
!     REAL(8)       cv(iCv,1:nCells) for iCv = 1..nCv, one variable after
!                   the other
!     REAL(8)       PICL only: piclVF(1:nCells), then piclgradRhog(i,1,1:nCells)
!                   for i = 1..3
!
!   The particle dump (written by ppiclf_io_WriteRestart) holds, with the
!   same stream layout:
!
!     CHARACTER(8)  'PPICLDMP'
!     INTEGER(4)    format version (1)
!     INTEGER(8)    step
!     REAL(8)       time
!     INTEGER(4)    nRanks
!     INTEGER(8)    nParticles (all ranks)
!     INTEGER(4)    lrs, lrp, lip (slots of y, rprop, iprop), lrp3 (history
!                   values per particle, 0 without VISCOUSUNSTEADY), nHist
!                   (history steps), nTimeBH (history steps filled)
!     REAL(8)       TimeBH(1:nHist)
!     INTEGER(4)    bin counts (3); REAL(8) bin box (xmin,xmax,ymin,...; 6)
!                   and bin size (3)
!     INTEGER(4)    nOut, then per output field: INTEGER(4) source (1 y,
!                   2 rprop, 3 iprop), INTEGER(4) slot, CHARACTER(24) name
!                   (ppiclf_user_OutputFields: the par*.vtu fields)
!     INTEGER(8)    particles on each rank (nRanks)
!     then one block per rank, in rank order, each for that rank's
!     particles: REAL(8) y(lrs,:), REAL(8) rprop(lrp,:), INTEGER(4)
!     iprop(lip,:), REAL(8) rprop3(lrp3,:)
!
!   Each file is written as <name>.tmp and renamed when it is complete, so
!   an interrupted run never leaves a partial dump under a final name.
!
! ******************************************************************************

MODULE RFLU_ModRestartData

  USE ISO_C_BINDING, ONLY: C_CHAR,C_INT,C_NULL_CHAR

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModMPI
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_WriteRestartData, &
            RFLU_CreateRestartDataDir, &
            RFLU_FinishRestartFile

  INTEGER, PARAMETER, PUBLIC :: RESTART_FORM_PURE = 1, &
                                RESTART_FORM_VOLUME_WEIGHTED = 2

  INTEGER, PARAMETER :: RESTART_DATA_VERSION = 1

  INTERFACE
    FUNCTION c_mkdir(path,mode) BIND(C,NAME='mkdir')
      IMPORT :: C_CHAR,C_INT
      CHARACTER(KIND=C_CHAR) :: path(*)
      INTEGER(C_INT), VALUE :: mode
      INTEGER(C_INT) :: c_mkdir
    END FUNCTION c_mkdir

    FUNCTION c_rename(oldPath,newPath) BIND(C,NAME='rename')
      IMPORT :: C_CHAR,C_INT
      CHARACTER(KIND=C_CHAR) :: oldPath(*),newPath(*)
      INTEGER(C_INT) :: c_rename
    END FUNCTION c_rename
  END INTERFACE

CONTAINS

! ******************************************************************************
!
! Purpose: Write the fluid dump of one region.
!
! Input:
!   pRegion     Region data
!
! Notes: rfluinit writes the gas state as initialized (pure, step 0); rflump
!   with particles holds the volume-weighted state and the PICL fields.
!
! Output: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteRestartData(pRegion)

    TYPE(t_region), POINTER :: pRegion

    CHARACTER(CHRLEN) :: dirName,fileName
    INTEGER :: errorFlag,form,i,iCv,iFile,nCells,withPicl
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

    global => pRegion%global
    pGrid  => pRegion%grid

    CALL RegisterFunction(global,'RFLU_WriteRestartData',__FILE__)

    nCells = pGrid%nCells

    form     = RESTART_FORM_PURE
    withPicl = 0
#ifdef PICL
    IF ( (global%piclUsed .EQV. .TRUE.) .AND. &
         (global%moduleType /= MODULE_TYPE_INIT) ) THEN
      form     = RESTART_FORM_VOLUME_WEIGHTED
      withPicl = 1
    END IF ! global%piclUsed
#endif

    CALL RFLU_CreateRestartDataDir(global,dirName)

    WRITE(fileName,'(A,I8.8,A,I5.5,A)') TRIM(dirName)//'/fluid_', &
          global%currentStep,'_',pRegion%iRegionGlobal,'.dat'

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,1X,A)') SOLVER_NAME,'Writing fluid dump', &
                                    TRIM(fileName)
    END IF ! global%verbLevel

    OPEN(NEWUNIT=iFile,FILE=TRIM(fileName)//'.tmp',FORM='UNFORMATTED', &
         ACCESS='STREAM',STATUS='REPLACE',IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(fileName)//'.tmp')
    END IF ! errorFlag

    WRITE(iFile) 'RFLUDUMP',INT(RESTART_DATA_VERSION,4), &
                 INT(global%currentStep,8),REAL(global%currentTime,8), &
                 INT(pRegion%iRegionGlobal,4),INT(nCells,4), &
                 INT(pRegion%mixtInput%nCv,4),INT(form,4),INT(withPicl,4)

    DO iCv = 1,pRegion%mixtInput%nCv
      WRITE(iFile) REAL(pRegion%mixt%cv(iCv,1:nCells),8)
    END DO ! iCv

#ifdef PICL
    IF ( withPicl == 1 ) THEN
      WRITE(iFile) REAL(pRegion%mixt%piclVF(1:nCells),8)
      DO i = 1,3
        WRITE(iFile) REAL(pRegion%mixt%piclgradRhog(i,1,1:nCells),8)
      END DO ! i
    END IF ! withPicl
#endif

    CLOSE(iFile,IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(fileName)//'.tmp')
    END IF ! errorFlag

    CALL RFLU_FinishRestartFile(global,fileName)

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteRestartData

! ******************************************************************************
!
! Purpose: Create restart_data/ (on the calling rank; an existing directory
!   is not an error, so every rank may call it without a barrier).
!
! Input:
!   global      Global data
!
! Output:
!   dirName     The directory name, without a trailing /
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateRestartDataDir(global,dirName)

    TYPE(t_global), POINTER :: global
    CHARACTER(*), INTENT(OUT) :: dirName

    INTEGER(C_INT) :: mkdirStatus

    dirName = TRIM(global%outDir)//'restart_data'
    mkdirStatus = c_mkdir(TRIM(dirName)//C_NULL_CHAR,INT(O'755',C_INT))

  END SUBROUTINE RFLU_CreateRestartDataDir

! ******************************************************************************
!
! Purpose: Give a dump written as <fileName>.tmp its final name.
!
! Input:
!   global      Global data
!   fileName    Final name of the dump
!
! Output: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_FinishRestartFile(global,fileName)

    TYPE(t_global), POINTER :: global
    CHARACTER(*), INTENT(IN) :: fileName

    IF ( c_rename(TRIM(fileName)//'.tmp'//C_NULL_CHAR, &
                  TRIM(fileName)//C_NULL_CHAR) /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(fileName))
    END IF ! c_rename

  END SUBROUTINE RFLU_FinishRestartFile


END MODULE RFLU_ModRestartData
