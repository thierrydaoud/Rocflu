! ******************************************************************************
!
! Purpose: Raw dumps of the fluid state in restart_data/, one file per region
!   and step: restart_data/fluid_<step>_<region>.dat (8 and 5 digits).
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
!   The file is written as <name>.tmp and renamed when it is complete, so an
!   interrupted run never leaves a partial dump under a final name.
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
  PUBLIC :: RFLU_WriteRestartData

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

! Every rank creates the directory (no barrier needed; an existing
! directory is not an error)
    dirName = TRIM(global%outDir)//'restart_data'
    errorFlag = c_mkdir(TRIM(dirName)//C_NULL_CHAR,INT(O'755',C_INT))

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

    IF ( c_rename(TRIM(fileName)//'.tmp'//C_NULL_CHAR, &
                  TRIM(fileName)//C_NULL_CHAR) /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(fileName))
    END IF ! c_rename

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteRestartData

END MODULE RFLU_ModRestartData
