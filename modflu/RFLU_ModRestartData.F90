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
!     INTEGER(4)    format version (4)
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
!     INTEGER(4)    number of NSCBC patches (non-reflecting boundaries, whose
!                   face state evolves in time), then per patch: INTEGER(4)
!                   global patch number, nBFaces; REAL(8) face state
!                   cv(iCv,1:nBFaces) for iCv = 1..nCv
!
!   The particle dump (written by ppiclf_io_WriteRestart) holds, with the
!   same stream layout:
!
!     CHARACTER(8)  'PPICLDMP'
!     INTEGER(4)    format version (2)
!     INTEGER(8)    step
!     REAL(8)       time
!     INTEGER(4)    nRanks
!     INTEGER(8)    nParticles (all ranks)
!     INTEGER(4)    lrs, lrp, lip (slots of y, rprop, iprop), lrp3 (history
!                   values per particle, 0 without VISCOUSUNSTEADY), lrpPro
!                   (projected fields), nHist (history steps), nTimeBH
!                   (history steps filled)
!     REAL(8)       TimeBH(1:nHist)
!     INTEGER(4)    bin counts (3); REAL(8) bin box (xmin,xmax,ymin,...; 6)
!                   and bin size (3)
!     INTEGER(4)    nOut, then per output field: INTEGER(4) source (1 y,
!                   2 rprop, 3 iprop), INTEGER(4) slot, CHARACTER(24) name
!                   (ppiclf_user_OutputFields: the par*.vtu fields)
!     INTEGER(8)    particles on each rank (nRanks), fluid cells on each rank
!                   (nRanks)
!     then one block per rank, in rank order: for that rank's particles
!     REAL(8) y(lrs,:), REAL(8) rprop(lrp,:), INTEGER(4) iprop(lip,:),
!     REAL(8) rprop3(lrp3,:), REAL(8) feedbk(lrpPro,:) (the feedback to the
!     gas); then REAL(8) the fields projected to its fluid cells, cells x
!     lrpPro. At the start of a stage the gas reads the projected fields and
!     ppiclF then projects feedbk of the stage before: a restart needs both
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
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_WriteRestartData, &
            RFLU_ReadRestartData, &
            RFLU_ChooseStartStep, &
            RFLU_CleanRestartData, &
            RFLU_CreateRestartDataDir, &
            RFLU_FinishRestartFile, &
            RFLU_ParticleDumpName, &
            RFLU_CompleteDumpSteps, &
            RFLU_FluidDumpTime

  INTEGER, PARAMETER, PUBLIC :: RESTART_FORM_PURE = 1, &
                                RESTART_FORM_VOLUME_WEIGHTED = 2

  INTEGER, PARAMETER :: RESTART_DATA_VERSION = 4

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
    INTEGER :: errorFlag,form,i,iCv,iFile,iPatch,nCells,nPatchesNSCBC,withPicl
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

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

! NSCBC face states (HasFaceStates)
    nPatchesNSCBC = 0
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      IF ( HasFaceStates(global,pPatch) ) nPatchesNSCBC = nPatchesNSCBC + 1
    END DO ! iPatch
    WRITE(iFile) INT(nPatchesNSCBC,4)
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      IF ( HasFaceStates(global,pPatch) ) THEN
        WRITE(iFile) INT(pPatch%iPatchGlobal,4),INT(pPatch%nBFaces,4)
        DO iCv = 1,pRegion%mixtInput%nCv
          WRITE(iFile) REAL(pPatch%mixt%cv(iCv,1:pPatch%nBFaces),8)
        END DO ! iCv
      END IF ! bcKind
    END DO ! iPatch

    CLOSE(iFile,IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(fileName)//'.tmp')
    END IF ! errorFlag

    CALL RFLU_FinishRestartFile(global,fileName)

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteRestartData

! ******************************************************************************
!
! Purpose: Whether a patch is an NSCBC patch whose face states exist: always
!   in rflump; in rfluinit only for the from-scratch and hard-coded initial
!   states (the other ones do not create them).
!
! ******************************************************************************

  LOGICAL FUNCTION HasFaceStates(global,pPatch)

    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch

    HasFaceStates = .FALSE.
    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      IF ( global%moduleType /= MODULE_TYPE_INIT ) THEN
        HasFaceStates = .TRUE.
      ELSE
        HasFaceStates = (global%initFlowFlag == INITFLOW_FROMSCRATCH) .OR. &
                        (global%initFlowFlag == INITFLOW_FROMHARDCODE)
      END IF ! global%moduleType
    END IF ! pPatch%bcKind

  END FUNCTION HasFaceStates

! ******************************************************************************
!
! Purpose: Number of NSCBC patches of a region.
!
! ******************************************************************************

  INTEGER FUNCTION CountNSCBC(pRegion)

    TYPE(t_region), POINTER :: pRegion

    INTEGER :: iPatch

    CountNSCBC = 0
    DO iPatch = 1,pRegion%grid%nPatches
      IF ( pRegion%patches(iPatch)%bcKind == BC_KIND_NSCBC ) THEN
        CountNSCBC = CountNSCBC + 1
      END IF ! bcKind
    END DO ! iPatch

  END FUNCTION CountNSCBC

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
! Purpose: Name of the particle dump at a step.
!
! ******************************************************************************

  FUNCTION RFLU_ParticleDumpName(global,step) RESULT(fileName)

    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(IN) :: step
    CHARACTER(CHRLEN) :: fileName

    WRITE(fileName,'(A,I8.8,A)') TRIM(global%outDir)// &
          'restart_data/particles_',step,'.dat'

  END FUNCTION RFLU_ParticleDumpName

! ******************************************************************************
!
! Purpose: Whether every dump a run needs to start from a step exists: the
!   fluid dump of every region and, with particles after step 0, the
!   particle dump (step 0 particles come from points.dat).
!
! ******************************************************************************

  LOGICAL FUNCTION DumpsComplete(global,step)

    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(IN) :: step

    INTEGER :: iReg,iRegLow,iRegUpp
    LOGICAL :: fileExists

! A single-region run has only region 0 (the serial region)
    IF ( global%nRegions == 1 ) THEN
      iRegLow = 0
      iRegUpp = 0
    ELSE
      iRegLow = 1
      iRegUpp = global%nRegions
    END IF ! global%nRegions

    DumpsComplete = .TRUE.
    DO iReg = iRegLow,iRegUpp
      INQUIRE(FILE=TRIM(FluidDumpName(global,step,iReg)),EXIST=fileExists)
      DumpsComplete = DumpsComplete .AND. fileExists
    END DO ! iReg
    IF ( (global%piclUsed .EQV. .TRUE.) .AND. (step > 0) ) THEN
      INQUIRE(FILE=TRIM(RFLU_ParticleDumpName(global,step)),EXIST=fileExists)
      DumpsComplete = DumpsComplete .AND. fileExists
    END IF ! global%piclUsed

  END FUNCTION DumpsComplete

! ******************************************************************************
!
! Purpose: The entries of restart_data/ (on the calling rank).
!
! Output:
!   entries     Names, without the directory
!
! ******************************************************************************

  SUBROUTINE ListRestartData(global,entries)

    TYPE(t_global), POINTER :: global
    CHARACTER(CHRLEN), ALLOCATABLE, INTENT(OUT) :: entries(:)

    CHARACTER(CHRLEN) :: dirName,entryName,listName
    INTEGER :: errorFlag,iFile,nEntries

    CALL RFLU_CreateRestartDataDir(global,dirName)
    listName = TRIM(dirName)//'/.listing'
    CALL EXECUTE_COMMAND_LINE('ls -1 '//TRIM(dirName)//' > '// &
                              TRIM(listName),EXITSTAT=errorFlag)

    nEntries = 0
    OPEN(NEWUNIT=iFile,FILE=TRIM(listName),FORM='FORMATTED', &
         STATUS='OLD',IOSTAT=errorFlag)
    IF ( errorFlag == 0 ) THEN
      DO
        READ(iFile,'(A)',IOSTAT=errorFlag) entryName
        IF ( errorFlag /= 0 ) EXIT
        nEntries = nEntries + 1
      END DO
      ALLOCATE(entries(nEntries))
      REWIND(iFile)
      DO nEntries = 1,SIZE(entries)
        READ(iFile,'(A)') entries(nEntries)
      END DO
      CLOSE(iFile,STATUS='DELETE')
    ELSE
      ALLOCATE(entries(0))
    END IF ! errorFlag

  END SUBROUTINE ListRestartData

! ******************************************************************************
!
! Purpose: Step of a dump from its name (-1 if the name is not a dump's).
!
! ******************************************************************************

  SUBROUTINE RFLU_CompleteDumpSteps(global,steps)

! Purpose: The steps of restart_data/ with every dump a run needs to start
!   from them (DumpsComplete), in increasing order (on the calling rank).

    TYPE(t_global), POINTER :: global
    INTEGER, ALLOCATABLE, INTENT(OUT) :: steps(:)

    CHARACTER(CHRLEN), ALLOCATABLE :: entries(:)
    INTEGER :: iEntry,nSteps,step
    INTEGER, ALLOCATABLE :: found(:)

    CALL ListRestartData(global,entries)
    ALLOCATE(found(SIZE(entries)))
    nSteps = 0
    DO iEntry = 1,SIZE(entries)
      IF ( entries(iEntry)(1:6) /= 'fluid_' ) CYCLE
      IF ( INDEX(entries(iEntry),'.tmp') > 0 ) CYCLE
      step = DumpStep(entries(iEntry))
      IF ( step < 0 ) CYCLE
      IF ( ANY(found(1:nSteps) == step) ) CYCLE
      IF ( DumpsComplete(global,step) .EQV. .FALSE. ) CYCLE
      nSteps = nSteps + 1
      found(nSteps) = step
    END DO ! iEntry

    ALLOCATE(steps(nSteps))
    steps = found(1:nSteps)
    CALL SortIntegers(steps)

  END SUBROUTINE RFLU_CompleteDumpSteps

  SUBROUTINE SortIntegers(a)

    INTEGER, INTENT(INOUT) :: a(:)

    INTEGER :: i,j,v

    DO i = 2,SIZE(a)
      v = a(i)
      j = i - 1
      DO WHILE ( j >= 1 )
        IF ( a(j) <= v ) EXIT
        a(j+1) = a(j)
        j = j - 1
      END DO
      a(j+1) = v
    END DO

  END SUBROUTINE SortIntegers

  REAL(8) FUNCTION RFLU_FluidDumpTime(global,step)

! Purpose: The time in the header of the first region's fluid dump of a step.

    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(IN) :: step

    CHARACTER(8) :: magic
    CHARACTER(CHRLEN) :: fileName
    INTEGER :: errorFlag,iFile
    INTEGER(4) :: version
    INTEGER(8) :: dumpStep8

    IF ( global%nRegions == 1 ) THEN
      fileName = FluidDumpName(global,step,0)
    ELSE
      fileName = FluidDumpName(global,step,1)
    END IF ! global%nRegions
    RFLU_FluidDumpTime = 0.0D0
    OPEN(NEWUNIT=iFile,FILE=TRIM(fileName),FORM='UNFORMATTED', &
         ACCESS='STREAM',STATUS='OLD',IOSTAT=errorFlag)
    IF ( errorFlag == 0 ) THEN
      READ(iFile,IOSTAT=errorFlag) magic,version,dumpStep8,RFLU_FluidDumpTime
      CLOSE(iFile)
    END IF ! errorFlag
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
    END IF ! errorFlag

  END FUNCTION RFLU_FluidDumpTime

  INTEGER FUNCTION DumpStep(entryName)

    CHARACTER(*), INTENT(IN) :: entryName

    INTEGER :: errorFlag

    DumpStep = -1
    errorFlag = 0
    IF ( entryName(1:6) == 'fluid_' ) THEN
      READ(entryName(7:14),'(I8)',IOSTAT=errorFlag) DumpStep
    ELSE IF ( entryName(1:10) == 'particles_' ) THEN
      READ(entryName(11:18),'(I8)',IOSTAT=errorFlag) DumpStep
    END IF ! entryName
    IF ( errorFlag /= 0 ) DumpStep = -1

  END FUNCTION DumpStep

! ******************************************************************************
!
! Purpose: Choose the step the run starts from (START_STEP) and set the start
!   step and time.
!
! Description: Rank 0 decides and broadcasts. START_STEP 0: a new run from
!   rfluinit's step 0 dump. -1: the latest step with every dump the run needs
!   (DumpsComplete). N: step N, which must have them. The time is read from
!   the header of the first region's dump.
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
    CHARACTER(CHRLEN), ALLOCATABLE :: entries(:)
    INTEGER :: chosenStep,errorFlag,iEntry,iFile,step
    INTEGER(4) :: version
    INTEGER(8) :: dumpStep8
    REAL(8) :: startTime

    CALL RegisterFunction(global,'RFLU_ChooseStartStep',__FILE__)

    chosenStep = 0
    startTime  = 0.0D0

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%startStep == 0 ) THEN
        IF ( DumpsComplete(global,0) .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_FILE_EXIST,__LINE__, &
            'restart_data/ has no step 0 fluid dump of every region: run '// &
            'rfluinit first')
        END IF ! DumpsComplete
      ELSE IF ( global%startStep == -1 ) THEN
        chosenStep = -1
        CALL ListRestartData(global,entries)
        DO iEntry = 1,SIZE(entries)
          step = DumpStep(entries(iEntry))
          IF ( (step > chosenStep) .AND. (step > 0) .AND. &
               (INDEX(entries(iEntry),'.tmp') == 0) ) THEN
            IF ( DumpsComplete(global,step) .EQV. .TRUE. ) chosenStep = step
          END IF ! step
        END DO ! iEntry
        IF ( chosenStep < 0 ) THEN
          CALL ErrorStop(global,ERR_FILE_EXIST,__LINE__, &
            'START_STEP -1: restart_data/ has no complete dump after step 0 '// &
            'to restart from')
        END IF ! chosenStep
      ELSE
        chosenStep = global%startStep
        IF ( DumpsComplete(global,chosenStep) .EQV. .FALSE. ) THEN
          WRITE(stepString,'(I0)') chosenStep
          CALL ErrorStop(global,ERR_FILE_EXIST,__LINE__, &
            'START_STEP '//TRIM(stepString)//': restart_data/ does not have '// &
            'the dumps of every region (and the particle dump) of that step')
        END IF ! DumpsComplete
      END IF ! global%startStep

! A single-region run has only region 0 (the serial region)
      IF ( global%nRegions == 1 ) THEN
        fileName = FluidDumpName(global,chosenStep,0)
      ELSE
        fileName = FluidDumpName(global,chosenStep,1)
      END IF ! global%nRegions
      OPEN(NEWUNIT=iFile,FILE=TRIM(fileName),FORM='UNFORMATTED', &
           ACCESS='STREAM',STATUS='OLD',IOSTAT=errorFlag)
      IF ( errorFlag == 0 ) THEN
        READ(iFile,IOSTAT=errorFlag) magic,version,dumpStep8,startTime
        CLOSE(iFile)
      END IF ! errorFlag
      IF ( errorFlag /= 0 .OR. magic /= 'RFLUDUMP' ) THEN
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(fileName))
      END IF ! errorFlag
    END IF ! global%myProcid

    CALL MPI_Bcast(chosenStep,1,MPI_INTEGER,MASTERPROC,global%mpiComm, &
                   errorFlag)
    CALL MPI_Bcast(startTime,1,MPI_DOUBLE_PRECISION,MASTERPROC, &
                   global%mpiComm,errorFlag)
    global%currentStep = chosenStep
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
!   the old-style output of earlier runs: SOL_* (except rfluinit's SOL_0),
!   PARAVIEW_*, par*.vtu, bin*.vtu and the .rin file; a restart leaves it (new
!   par/bin files continue its numbering). The number removed is in
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
    CHARACTER(CHRLEN), ALLOCATABLE :: entries(:)
    CHARACTER(512) :: command
    INTEGER :: errorFlag,iEntry,iFile,nRemoved,step

    CALL RegisterFunction(global,'RFLU_CleanRestartData',__FILE__)

    nRemoved = 0

    IF ( global%myProcid == MASTERPROC ) THEN
      CALL RFLU_CreateRestartDataDir(global,dirName)
      listName = TRIM(dirName)//'/.listing'

! Dumps newer than the start step, and partial dumps
      CALL ListRestartData(global,entries)
      DO iEntry = 1,SIZE(entries)
        step = DumpStep(entries(iEntry))
        IF ( (step > global%currentStep) .OR. &
             ((step >= 0) .AND. (INDEX(entries(iEntry),'.tmp') > 0)) ) THEN
          IF ( c_unlink(TRIM(dirName)//'/'//TRIM(entries(iEntry))// &
                        C_NULL_CHAR) == 0 ) THEN
            nRemoved = nRemoved + 1
          END IF ! c_unlink
        END IF ! step
      END DO ! iEntry

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
! Purpose: Read the fluid dump of one region at the start step: the gas state
!   (real cells; the virtual cells are filled by the exchange and the boundary
!   conditions that follow), the PICL fields and the NSCBC face states.
!
! Description: Step 0 (rfluinit) holds the pure gas state; a later dump holds
!   the state of rflump, volume-weighted with particles, which is used as is.
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
    INTEGER :: errorFlag,expectedForm,i,iCv,iFile,iPatch,iPatchDump,iVar
    INTEGER(4) :: dumpRegion,form,nBFacesDump,nCellsDump,nCvDump, &
                  nPatchesNSCBC,patchGlobalDump,version,withPicl
    INTEGER(8) :: step
    REAL(8) :: dumpTime
    REAL(8), ALLOCATABLE :: values(:)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

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
    IF ( errorFlag /= 0 .OR. magic /= 'RFLUDUMP' .OR. &
         version /= RESTART_DATA_VERSION ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(fileName)// &
        ': not a fluid dump of this version (rerun rfluinit for step 0)')
    END IF ! errorFlag
    IF ( nCellsDump /= pGrid%nCells .OR. &
         nCvDump /= pRegion%mixtInput%nCv ) THEN
      CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,TRIM(fileName)// &
        ': cells or variables differ from the grid of this region')
    END IF ! nCellsDump
    expectedForm = RESTART_FORM_PURE
    IF ( (global%currentStep > 0) .AND. (global%piclUsed .EQV. .TRUE.) ) THEN
      expectedForm = RESTART_FORM_VOLUME_WEIGHTED
    END IF ! global%currentStep
    IF ( form /= expectedForm ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(fileName)// &
        ': written by a run with particles on/off differently from this one')
    END IF ! form

    ALLOCATE(values(pGrid%nCells))
    DO iCv = 1,pRegion%mixtInput%nCv
      READ(iFile,IOSTAT=errorFlag) values
      IF ( errorFlag /= 0 ) THEN
        CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
      END IF ! errorFlag
      pRegion%mixt%cv(iCv,1:pGrid%nCells) = REAL(values,RFREAL)
    END DO ! iCv

! PICL fields: 4 arrays (piclVF, piclgradRhog x, y, z)
    IF ( withPicl == 1 ) THEN
      DO iVar = 1,4
        READ(iFile,IOSTAT=errorFlag) values
        IF ( errorFlag /= 0 ) THEN
          CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
        END IF ! errorFlag
#ifdef PICL
        IF ( iVar == 1 ) THEN
          pRegion%mixt%piclVF(1:pGrid%nCells) = REAL(values,RFREAL)
        ELSE
          i = iVar - 1
          pRegion%mixt%piclgradRhog(i,1,1:pGrid%nCells) = REAL(values,RFREAL)
        END IF ! iVar
#endif
      END DO ! iVar
    END IF ! withPicl

    DEALLOCATE(values)

! NSCBC face states, matched to the local patches by global patch number
    READ(iFile,IOSTAT=errorFlag) nPatchesNSCBC
    IF ( errorFlag /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
    END IF ! errorFlag
    DO iPatchDump = 1,nPatchesNSCBC
      READ(iFile,IOSTAT=errorFlag) patchGlobalDump,nBFacesDump
      pPatch => NULL()
      DO iPatch = 1,pGrid%nPatches
        IF ( pRegion%patches(iPatch)%iPatchGlobal == patchGlobalDump ) THEN
          pPatch => pRegion%patches(iPatch)
        END IF ! iPatchGlobal
      END DO ! iPatch
      IF ( errorFlag /= 0 .OR. (.NOT. ASSOCIATED(pPatch)) ) THEN
        CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName)// &
          ': NSCBC patch not found in this region')
      END IF ! errorFlag
      IF ( pPatch%bcKind /= BC_KIND_NSCBC .OR. &
           pPatch%nBFaces /= nBFacesDump ) THEN
        CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(fileName)// &
          ': the NSCBC patches differ from the .bc file and grid')
      END IF ! pPatch%bcKind
      ALLOCATE(values(nBFacesDump))
      DO iCv = 1,pRegion%mixtInput%nCv
        READ(iFile,IOSTAT=errorFlag) values
        IF ( errorFlag /= 0 ) THEN
          CALL ErrorStop(global,ERR_FILE_READ,__LINE__,TRIM(fileName))
        END IF ! errorFlag
        pPatch%mixt%cv(iCv,1:nBFacesDump) = REAL(values,RFREAL)
      END DO ! iCv
      pPatch%mixt%cvState = CV_MIXT_STATE_CONS
      DEALLOCATE(values)
    END DO ! iPatchDump
    IF ( nPatchesNSCBC /= CountNSCBC(pRegion) ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,TRIM(fileName)// &
        ': no face states for every NSCBC patch of this region')
    END IF ! nPatchesNSCBC

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
