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
            RFLU_ReadRestartData, &
            RFLU_ChooseStartStep, &
            RFLU_CleanRestartData, &
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

    FUNCTION c_unlink(path) BIND(C,NAME='unlink')
      IMPORT :: C_CHAR,C_INT
      CHARACTER(KIND=C_CHAR) :: path(*)
      INTEGER(C_INT) :: c_unlink
    END FUNCTION c_unlink
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
! Purpose: Name of the fluid dump of a region at a step.
!
! ******************************************************************************

  FUNCTION FluidDumpName(global,step,iRegionGlobal) RESULT(fileName)

    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(IN) :: step,iRegionGlobal
    CHARACTER(CHRLEN) :: fileName

    WRITE(fileName,'(A,I8.8,A,I5.5,A)') TRIM(global%outDir)// &
          'restart_data/fluid_',step,'_',iRegionGlobal,'.dat'

  END FUNCTION FluidDumpName

! ******************************************************************************
!
! Purpose: Choose the step the run starts from (START_STEP) and set the start
!   step and time.
!
! Description: START_STEP 0 starts a new run from rfluinit's step 0 dump;
!   restarting from a later dump (START_STEP -1 or N) is not available yet.
!   The time is read from the header of the first region's dump by rank 0 and
!   broadcast.
!
! Input:
!   global      Global data
!
! Output: None (global%currentStep, global%currentTime set).
!
! ******************************************************************************

  SUBROUTINE RFLU_ChooseStartStep(global)

    TYPE(t_global), POINTER :: global

    CHARACTER(8) :: magic
    CHARACTER(CHRLEN) :: fileName,stepString
    INTEGER :: errorFlag,iFile
    INTEGER(4) :: version
    INTEGER(8) :: step
    REAL(8) :: startTime

    CALL RegisterFunction(global,'RFLU_ChooseStartStep',__FILE__)

    IF ( global%startStep /= 0 ) THEN
      WRITE(stepString,'(I0)') global%startStep
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
        'START_STEP '//TRIM(stepString)//': restarting from a dump is not '// &
        'available yet; use START_STEP 0 (a new run)')
    END IF ! global%startStep

    global%currentStep = 0

    IF ( global%myProcid == MASTERPROC ) THEN
! A single-region run has only region 0 (the serial region)
      IF ( global%nRegions == 1 ) THEN
        fileName = FluidDumpName(global,0,0)
      ELSE
        fileName = FluidDumpName(global,0,1)
      END IF ! global%nRegions
      OPEN(NEWUNIT=iFile,FILE=TRIM(fileName),FORM='UNFORMATTED', &
           ACCESS='STREAM',STATUS='OLD',IOSTAT=errorFlag)
      IF ( errorFlag /= 0 ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(fileName)// &
          ' (the step 0 dump of rfluinit: run rfluinit first)')
      END IF ! errorFlag
      READ(iFile,IOSTAT=errorFlag) magic,version,step,startTime
      CLOSE(iFile)
      IF ( errorFlag /= 0 .OR. magic /= 'RFLUDUMP' ) THEN
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(fileName))
      END IF ! errorFlag
    END IF ! global%myProcid

    CALL MPI_Bcast(startTime,1,MPI_DOUBLE_PRECISION,MASTERPROC, &
                   global%mpiComm,errorFlag)
    global%currentTime = REAL(startTime,RFREAL)

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ChooseStartStep

! ******************************************************************************
!
! Purpose: Remove what an earlier run in this directory left behind, so that
!   the run's results replace it.
!
! Description: Rank 0 removes the dumps in restart_data/ newer than the start
!   step and any partial (.tmp) dump. A new run (start step 0) also removes
!   the old-style output of earlier runs: SOL_* (except rfluinit's SOL_0, which
!   still holds the boundary data read at start-up), PARAVIEW_*, par*.vtu,
!   bin*.vtu and the .rin file. The number removed is in
!   global%nFilesRemoved (for the Setup box). Called by every rank.
!
! Input:
!   global      Global data
!
! Output: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CleanRestartData(global)

    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN) :: dirName,entryName,listName
    CHARACTER(512) :: command
    INTEGER :: errorFlag,iFile,nRemoved,step

    CALL RegisterFunction(global,'RFLU_CleanRestartData',__FILE__)

    nRemoved = 0

    IF ( global%myProcid == MASTERPROC ) THEN
      CALL RFLU_CreateRestartDataDir(global,dirName)
      listName = TRIM(dirName)//'/.listing'

! Dumps newer than the start step, and partial dumps
      CALL EXECUTE_COMMAND_LINE('ls -1 '//TRIM(dirName)//' > '// &
                                TRIM(listName),EXITSTAT=errorFlag)
      OPEN(NEWUNIT=iFile,FILE=TRIM(listName),FORM='FORMATTED', &
           STATUS='OLD',IOSTAT=errorFlag)
      IF ( errorFlag == 0 ) THEN
        DO
          READ(iFile,'(A)',IOSTAT=errorFlag) entryName
          IF ( errorFlag /= 0 ) EXIT

          step = -1
          IF ( entryName(1:6) == 'fluid_' ) THEN
            READ(entryName(7:14),'(I8)',IOSTAT=errorFlag) step
          ELSE IF ( entryName(1:10) == 'particles_' ) THEN
            READ(entryName(11:18),'(I8)',IOSTAT=errorFlag) step
          END IF ! entryName
          IF ( errorFlag /= 0 ) step = -1

          IF ( (step > global%currentStep) .OR. &
               ((step >= 0) .AND. (INDEX(entryName,'.tmp') > 0)) ) THEN
            IF ( c_unlink(TRIM(dirName)//'/'//TRIM(entryName)// &
                          C_NULL_CHAR) == 0 ) THEN
              nRemoved = nRemoved + 1
            END IF ! c_unlink
          END IF ! step
        END DO
        CLOSE(iFile,STATUS='DELETE')
      END IF ! errorFlag

! Old-style output of earlier runs
      IF ( global%currentStep == 0 ) THEN
        WRITE(command,'(9A)') 'find '//TRIM(global%outDir)//' -maxdepth 1 ', &
          '\( -name "SOL_*" ! -name "SOL_0.00000E+00" ', &
          '-o -name "PARAVIEW_*" ', &
          '-o -name "par[0-9]*.vtu" -o -name "bin[0-9]*.vtu" ', &
          '-o -name "'//TRIM(global%caseName)//'.rin" \) ', &
          '-print -exec rm -rf {} + > ',TRIM(listName)
        CALL EXECUTE_COMMAND_LINE(TRIM(command),EXITSTAT=errorFlag)
        OPEN(NEWUNIT=iFile,FILE=TRIM(listName),FORM='FORMATTED', &
             STATUS='OLD',IOSTAT=errorFlag)
        IF ( errorFlag == 0 ) THEN
          DO
            READ(iFile,'(A)',IOSTAT=errorFlag) entryName
            IF ( errorFlag /= 0 ) EXIT
            nRemoved = nRemoved + 1
          END DO
          CLOSE(iFile,STATUS='DELETE')
        END IF ! errorFlag
      END IF ! global%currentStep
    END IF ! global%myProcid

    global%nFilesRemoved = nRemoved

! Nothing may be written before rank 0 has finished
    CALL MPI_Barrier(global%mpiComm,errorFlag)

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CleanRestartData

! ******************************************************************************
!
! Purpose: Read the fluid dump of one region at the start step into the gas
!   state (real cells; the virtual cells are filled by the exchange and the
!   boundary conditions that follow).
!
! Input:
!   pRegion     Region data
!
! Output: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadRestartData(pRegion)

    TYPE(t_region), POINTER :: pRegion

    CHARACTER(8) :: magic
    CHARACTER(CHRLEN) :: fileName
    INTEGER :: errorFlag,iCv,iFile
    INTEGER(4) :: dumpRegion,form,nCellsDump,nCvDump,version,withPicl
    INTEGER(8) :: step
    REAL(8) :: dumpTime
    REAL(8), ALLOCATABLE :: values(:)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

    global => pRegion%global
    pGrid  => pRegion%grid

    CALL RegisterFunction(global,'RFLU_ReadRestartData',__FILE__)

    fileName = FluidDumpName(global,global%currentStep,pRegion%iRegionGlobal)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,1X,A)') SOLVER_NAME,'Reading fluid dump', &
                                    TRIM(fileName)
    END IF ! global%verbLevel

    OPEN(NEWUNIT=iFile,FILE=TRIM(fileName),FORM='UNFORMATTED', &
         ACCESS='STREAM',STATUS='OLD',IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(fileName))
    END IF ! errorFlag

    READ(iFile,IOSTAT=errorFlag) magic,version,step,dumpTime,dumpRegion, &
                                 nCellsDump,nCvDump,form,withPicl
    IF ( errorFlag /= 0 .OR. magic /= 'RFLUDUMP' .OR. version /= 1 ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(fileName))
    END IF ! errorFlag
    IF ( nCellsDump /= pGrid%nCells .OR. &
         nCvDump /= pRegion%mixtInput%nCv ) THEN
      CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,TRIM(fileName)// &
        ': cells or variables differ from the grid of this region')
    END IF ! nCellsDump
    IF ( form /= RESTART_FORM_PURE ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(fileName)// &
        ': a dump written by rflump; restarting from it is not available yet')
    END IF ! form

    ALLOCATE(values(pGrid%nCells))
    DO iCv = 1,pRegion%mixtInput%nCv
      READ(iFile,IOSTAT=errorFlag) values
      IF ( errorFlag /= 0 ) THEN
        CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
      END IF ! errorFlag
      pRegion%mixt%cv(iCv,1:pGrid%nCells) = REAL(values,RFREAL)
    END DO ! iCv
    DEALLOCATE(values)

    CLOSE(iFile)

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadRestartData

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
