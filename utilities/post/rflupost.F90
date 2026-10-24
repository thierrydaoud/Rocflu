! ******************************************************************************
!
! Purpose: Post program: convert the dumps in restart_data/ into files for
!   ParaView in post/.
!
! Description: For each step converted:
!     post/p<k>/<step>.vtu          fluid of region k+1 (k = 0 with one region)
!     post/particles/<step>.vtu     the particles (with particles)
!     post/bins/<step>.vtu          the ppiclF bins (with particles)
!     post/root/collection_<step>.vtm   the three above as one dataset
!   and post/collection.pvd lists every step in post/ with its time (open it
!   in ParaView), post/fluid.pvd, post/particles.pvd and post/bins.pvd each
!   one kind (open them together to show, hide and filter each on its own). The fluid fields are chosen with the *_WRT keys of # POST;
!   pressure, temperature and sound speed come from the solver's own
!   RFLU_SetDependentVars. The particle fields are those of the dump's
!   output-field table (PARTICLE_WRT, GAS_WRT, ... of # PICL).
!
!   The steps: 'new' (those without post/root/collection_<step>.vtm), 'all',
!   'last', or a comma-separated list. Runs on one rank, region after region.
!
! Input:
!   caseString    Case name
!   verbLevel     Verbosity
!   stepsString   Steps to convert
!
! Output: None.
!
! ******************************************************************************

SUBROUTINE rflupost(caseString,verbLevel,stepsString)

  USE ModError
  USE ModDataTypes
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE RFLU_ModBoundXvUtils, ONLY: RFLU_BXV_CreateVarsCv, &
                                  RFLU_BXV_CreateVarsDv, &
                                  RFLU_BXV_DestroyVarsCv, &
                                  RFLU_BXV_DestroyVarsDv
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModRenumberings
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRegionMapping
  USE RFLU_ModRestartData, ONLY: RFLU_CompleteDumpSteps, &
                                 RFLU_FluidDumpTime, &
                                 RFLU_ParticleDumpName, &
                                 RFLU_ReadRestartData
  USE RFLU_ModPostWrite
  USE RFLU_ModPostParticles

  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, &
                           RFLU_BuildDataStruct, &
                           RFLU_CreateGrid, &
                           RFLU_DeallocMemSolWrapper, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_PrintHeader, &
                           RFLU_SetDependentVars, &
                           RFLU_SetGasVars, &
                           RFLU_SetModuleType, &
                           RFLU_SetVarInfoWrapper, &
                           RFLU_WriteVersionString

  IMPLICIT NONE

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString,stepsString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: binFile,casename,fileName,particleFile,regionDir
  CHARACTER(CHRLEN), ALLOCATABLE :: binFiles(:,:),fluidFiles(:,:), &
                                    particleFiles(:,:),regionFiles(:), &
                                    stepFiles(:)
  INTEGER :: errorFlag,iReg,iRegLow,iRegUpp,iStep,nParticleSteps,nSteps,step
  INTEGER, ALLOCATABLE :: available(:),steps(:)
  LOGICAL :: fileExists
  REAL(8), ALLOCATABLE :: particleTimes(:),times(:)
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start: as rfluinit, every region of the case on this rank
! ******************************************************************************

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag

  casename = caseString(1:LEN(caseString))
  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

  CALL RegisterFunction(global,'rflupost',__FILE__)

  CALL RFLU_SetModuleType(global,MODULE_TYPE_POSTPROC)

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_WriteVersionString(global)
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global%myProcid

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.)

! A single-region case has only region 0 (the serial region)
  IF ( global%nRegions == 1 ) THEN
    iRegLow = 0
    iRegUpp = 0
  ELSE
    iRegLow = 1
    iRegUpp = global%nRegions
  END IF ! global%nRegions

! ******************************************************************************
! Steps to convert
! ******************************************************************************

  CALL RFLU_CompleteDumpSteps(global,available)
  CALL ChooseSteps(TRIM(stepsString),available,steps)
  nSteps = SIZE(steps)

  CALL EXECUTE_COMMAND_LINE('mkdir -p post/root post/particles post/bins', &
                            EXITSTAT=errorFlag)

  WRITE(STDOUT,'(A,I0,A,I0,A)') 'rflupost: converting ',nSteps,' of ', &
                                SIZE(available),' dump steps to post/'

! ******************************************************************************
! Fluid: region after region, each step of the region
! ******************************************************************************

  ALLOCATE(regionFiles(iRegUpp-iRegLow+1))

  IF ( nSteps > 0 ) THEN
    DO iReg = iRegLow,iRegUpp
      pRegion => levels(1)%regions(iReg)

      regionDir = RegionDirName(iReg)
      CALL EXECUTE_COMMAND_LINE('mkdir -p post/'//TRIM(regionDir), &
                                EXITSTAT=errorFlag)

      CALL SetUpRegion(pRegion)
      DO iStep = 1,nSteps
        CALL WriteRegion(pRegion,steps(iStep),regionDir)
      END DO ! iStep
      CALL DestroyRegion(pRegion)
    END DO ! iReg
  END IF ! nSteps

! ******************************************************************************
! Particles, bins and the dataset of each step
! ******************************************************************************

  DO iReg = iRegLow,iRegUpp
    WRITE(regionFiles(iReg-iRegLow+1),'(A,A,A)') '../',TRIM(RegionDirName(iReg)), &
          '/'
  END DO ! iReg

  DO iStep = 1,nSteps
    step = steps(iStep)
    CALL WriteParticlesAndBins(step,fileExists)
    particleFile = ''
    binFile      = ''
    IF ( fileExists ) THEN
      particleFile = StepFileName('../particles/',step)
      binFile      = StepFileName('../bins/',step)
    END IF ! fileExists
    CALL WriteStepVTM(TRIM(StepVTMName(step)),PrefixStep(regionFiles,step), &
                      TRIM(particleFile),TRIM(binFile))
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,I0,A,ES12.5)') 'rflupost:   step ',step,', t = ', &
                                      RFLU_FluidDumpTime(global,step)
    END IF ! global%verbLevel
  END DO ! iStep

! ******************************************************************************
! The series: every step that has its dataset in post/. collection.pvd has
! everything; fluid.pvd, particles.pvd and bins.pvd one kind each, so that
! ParaView shows them as separate objects
! ******************************************************************************

  ALLOCATE(times(SIZE(available)),stepFiles(SIZE(available)), &
           fluidFiles(SIZE(regionFiles),SIZE(available)), &
           particleTimes(SIZE(available)),particleFiles(1,SIZE(available)), &
           binFiles(1,SIZE(available)))
  nSteps = 0
  nParticleSteps = 0
  DO iStep = 1,SIZE(available)
    step = available(iStep)
    INQUIRE(FILE=TRIM(StepVTMName(step)),EXIST=fileExists)
    IF ( .NOT. fileExists ) CYCLE
    nSteps = nSteps + 1
    times(nSteps) = RFLU_FluidDumpTime(global,step)
    stepFiles(nSteps) = StepFileName('root/collection_',step)
    stepFiles(nSteps) = stepFiles(nSteps)(1:LEN_TRIM(stepFiles(nSteps))-4)// &
                        '.vtm'
    DO iReg = iRegLow,iRegUpp
      fluidFiles(iReg-iRegLow+1,nSteps) = &
        StepFileName(TRIM(RegionDirName(iReg))//'/',step)
    END DO ! iReg
    INQUIRE(FILE=TRIM(StepFileName('post/particles/',step)),EXIST=fileExists)
    IF ( fileExists ) THEN
      nParticleSteps = nParticleSteps + 1
      particleTimes(nParticleSteps) = times(nSteps)
      particleFiles(1,nParticleSteps) = StepFileName('particles/',step)
      binFiles(1,nParticleSteps)      = StepFileName('bins/',step)
    END IF ! fileExists
  END DO ! iStep
  CALL WriteSeriesPVD('post/collection.pvd',times(1:nSteps), &
                      RESHAPE(stepFiles(1:nSteps),(/1,nSteps/)))
  CALL WriteSeriesPVD('post/fluid.pvd',times(1:nSteps),fluidFiles(:,1:nSteps))
  IF ( nParticleSteps > 0 ) THEN
    CALL WriteSeriesPVD('post/particles.pvd',particleTimes(1:nParticleSteps), &
                        particleFiles(:,1:nParticleSteps))
    CALL WriteSeriesPVD('post/bins.pvd',particleTimes(1:nParticleSteps), &
                        binFiles(:,1:nParticleSteps))
  END IF ! nParticleSteps
  WRITE(STDOUT,'(A,I0,A)') 'rflupost: post/collection.pvd and fluid.pvd list ', &
                           nSteps,' steps'
  IF ( nParticleSteps > 0 ) THEN
    WRITE(STDOUT,'(A,I0,A)') 'rflupost: post/particles.pvd and bins.pvd list ', &
                             nParticleSteps,' steps'
  END IF ! nParticleSteps

  CALL DeregisterFunction(global)

CONTAINS

! ******************************************************************************
! Names
! ******************************************************************************

! post/ directory of a region: p<k>, k = region - 1 (0 for a single region)
  FUNCTION RegionDirName(iRegion) RESULT(dirName)
    INTEGER, INTENT(IN) :: iRegion
    CHARACTER(CHRLEN) :: dirName
    WRITE(dirName,'(A,I0)') 'p',MAX(iRegion-1,0)
  END FUNCTION RegionDirName

  FUNCTION StepFileName(prefix,stepNumber) RESULT(name)
    CHARACTER(*), INTENT(IN) :: prefix
    INTEGER, INTENT(IN) :: stepNumber
    CHARACTER(CHRLEN) :: name
    WRITE(name,'(A,I8.8,A)') prefix,stepNumber,'.vtu'
  END FUNCTION StepFileName

  FUNCTION StepVTMName(stepNumber) RESULT(name)
    INTEGER, INTENT(IN) :: stepNumber
    CHARACTER(CHRLEN) :: name
    WRITE(name,'(A,I8.8,A)') 'post/root/collection_',stepNumber,'.vtm'
  END FUNCTION StepVTMName

  FUNCTION PrefixStep(dirs,stepNumber) RESULT(names)
    CHARACTER(*), INTENT(IN) :: dirs(:)
    INTEGER, INTENT(IN) :: stepNumber
    CHARACTER(CHRLEN) :: names(SIZE(dirs))
    INTEGER :: i
    DO i = 1,SIZE(dirs)
      names(i) = StepFileName(TRIM(dirs(i)),stepNumber)
    END DO ! i
  END FUNCTION PrefixStep

! ******************************************************************************
! The steps to convert: 'new', 'all', 'last' or 'N,M,...'
! ******************************************************************************

  SUBROUTINE ChooseSteps(selection,available,chosen)

    CHARACTER(*), INTENT(IN) :: selection
    INTEGER, INTENT(IN) :: available(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: chosen(:)

    CHARACTER(LEN(selection)+1) :: rest
    INTEGER :: comma,i,n,wanted
    INTEGER :: picked(SIZE(available))
    LOGICAL :: done

    n = 0
    SELECT CASE ( selection )
      CASE ( 'all' )
        n = SIZE(available)
        picked(1:n) = available
      CASE ( 'last' )
        IF ( SIZE(available) > 0 ) THEN
          n = 1
          picked(1) = available(SIZE(available))
        END IF ! SIZE
      CASE ( 'new' )
        DO i = 1,SIZE(available)
          INQUIRE(FILE=TRIM(StepVTMName(available(i))),EXIST=done)
          IF ( .NOT. done ) THEN
            n = n + 1
            picked(n) = available(i)
          END IF ! done
        END DO ! i
      CASE DEFAULT
        rest = selection//','
        DO WHILE ( LEN_TRIM(rest) > 0 )
          comma = INDEX(rest,',')
          READ(rest(1:comma-1),*,IOSTAT=errorFlag) wanted
          IF ( errorFlag /= 0 .OR. .NOT. ANY(available == wanted) ) THEN
            WRITE(STDOUT,'(A,A,A)') 'rflupost: step ',rest(1:comma-1), &
              ' has no complete dump in restart_data/'
            STOP 1
          END IF ! errorFlag
          n = n + 1
          picked(MIN(n,SIZE(picked))) = wanted
          rest = rest(comma+1:)
        END DO
    END SELECT ! selection

    ALLOCATE(chosen(n))
    chosen = picked(1:n)

  END SUBROUTINE ChooseSteps

! ******************************************************************************
! Grid and arrays of a region
! ******************************************************************************

  SUBROUTINE SetUpRegion(pRegion)

    TYPE(t_region), POINTER :: pRegion

    CALL RFLU_ReadDimensions(pRegion)
    CALL RFLU_CreateGrid(pRegion)
    IF ( pRegion%grid%nPatches > 0 ) THEN
      CALL RFLU_ReadBCInputFileWrapper(pRegion)
    END IF ! pRegion%grid%nPatches
    CALL RFLU_ReadGridWrapper(pRegion)

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)

! Cell volumes, with particles: phi at step 0 is the projected particle
! volume over the cell volume (as in PICL_InitSolver)
    IF ( global%piclUsed .EQV. .TRUE. ) THEN
      CALL RFLU_CreateBVertexLists(pRegion)
      CALL RFLU_BuildBVertexLists(pRegion)
      CALL RFLU_CreateFaceList(pRegion)
      CALL RFLU_BuildFaceList(pRegion)
      CALL RFLU_RenumberBFaceLists(pRegion)
      CALL RFLU_CreateGeometry(pRegion)
      CALL RFLU_BuildGeometry(pRegion)
    END IF ! global%piclUsed

    CALL RFLU_AllocMemSolWrapper(pRegion)
    CALL RFLU_SetVarInfoWrapper(pRegion)

! The dump holds the NSCBC face states
    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_CreateVarsCv(pRegion)
      CALL RFLU_BXV_CreateVarsDv(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC

  END SUBROUTINE SetUpRegion

  SUBROUTINE DestroyRegion(pRegion)

    TYPE(t_region), POINTER :: pRegion

    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_DestroyVarsCv(pRegion)
      CALL RFLU_BXV_DestroyVarsDv(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC
    CALL RFLU_DeallocMemSolWrapper(pRegion)
    IF ( global%piclUsed .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyGeometry(pRegion)
      CALL RFLU_DestroyFaceList(pRegion)
      CALL RFLU_DestroyBVertexLists(pRegion)
    END IF ! global%piclUsed
    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyGrid(pRegion)

  END SUBROUTINE DestroyRegion

! ******************************************************************************
! Fluid file of a region at a step
! ******************************************************************************

  SUBROUTINE WriteRegion(pRegion,stepNumber,dirName)

    TYPE(t_region), POINTER :: pRegion
    INTEGER, INTENT(IN) :: stepNumber
    CHARACTER(*), INTENT(IN) :: dirName

    INTEGER :: c,icg,icl,iCellType,k,nCells,nFields,nVertCell
    INTEGER(1), ALLOCATABLE :: cellTypes(:)
    INTEGER(8) :: nConn
    INTEGER(8), ALLOCATABLE :: connectivity(:),offsets(:)
    REAL(8), ALLOCATABLE :: rhoG(:),phi(:),work(:,:)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_postField), ALLOCATABLE :: fields(:)
    TYPE(t_particleDump) :: particleDump

    pGrid  => pRegion%grid
    nCells = pGrid%nCells

! State of the step: gas and PICL fields, NSCBC faces
    global%currentStep = stepNumber
    CALL RFLU_ReadRestartData(pRegion)
    pRegion%mixt%cvState = CV_MIXT_STATE_CONS
    CALL RFLU_SetGasVars(pRegion,1,nCells)
    CALL RFLU_SetDependentVars(pRegion,1,nCells)

! Particle volume fraction: in the dump after step 0; at step 0 (rfluinit's
! dump has no PICL fields) the one projected at the start of the run, from
! the particle dump (the block of the rank of this region)
    ALLOCATE(phi(nCells),rhoG(nCells))
    phi = 0.0D0
#ifdef PICL
    IF ( global%piclUsed .EQV. .TRUE. ) THEN
      IF ( stepNumber > 0 ) THEN
        phi = pRegion%mixt%piclVF(1:nCells)
      ELSE
        INQUIRE(FILE=TRIM(RFLU_ParticleDumpName(global,0)),EXIST=fileExists)
        IF ( fileExists ) THEN
          CALL ReadParticleDump(TRIM(RFLU_ParticleDumpName(global,0)), &
                                particleDump)
          k = MAX(pRegion%iRegionGlobal,1)
          IF ( k <= particleDump%nRanks ) THEN
            IF ( SIZE(particleDump%ranks(k)%phi) == nCells ) THEN
              phi = particleDump%ranks(k)%phi/pGrid%vol(1:nCells)
            END IF ! SIZE
          END IF ! k
        END IF ! fileExists
      END IF ! stepNumber
    END IF ! global%piclUsed
#endif

! Gas density: the dump after step 0 holds (1-phi)*rho with particles
    rhoG = pRegion%mixt%cv(CV_MIXT_DENS,1:nCells)
    IF ( (global%piclUsed .EQV. .TRUE.) .AND. (stepNumber > 0) ) THEN
      rhoG = rhoG/(1.0D0 - phi)
    END IF ! global%piclUsed

! Fields, as chosen with the *_WRT keys (names as in the old ParaView files)
    nFields = 0
    ALLOCATE(work(3,nCells))
    IF ( global%postWriteDensity ) THEN
      CALL AddPostField(fields,nFields,'Gas Density',RESHAPE(rhoG,(/1,nCells/)))
    END IF ! global%postWriteDensity
    IF ( global%postWriteVelocity ) THEN
      DO c = 1,3
        work(c,:) = pRegion%mixt%cv(CV_MIXT_XMOM+c-1,1:nCells)/ &
                    pRegion%mixt%cv(CV_MIXT_DENS,1:nCells)
      END DO ! c
      CALL AddPostField(fields,nFields,'Velocity',work)
    END IF ! global%postWriteVelocity
    IF ( global%postWritePressure ) THEN
      CALL AddPostField(fields,nFields,'Pressure', &
        RESHAPE(REAL(pRegion%mixt%dv(DV_MIXT_PRES,1:nCells),8),(/1,nCells/)))
    END IF ! global%postWritePressure
    IF ( global%postWriteTemperature ) THEN
      CALL AddPostField(fields,nFields,'Temperature', &
        RESHAPE(REAL(pRegion%mixt%dv(DV_MIXT_TEMP,1:nCells),8),(/1,nCells/)))
    END IF ! global%postWriteTemperature
    IF ( global%postWriteSoundSpeed ) THEN
      CALL AddPostField(fields,nFields,'Speed of Sound', &
        RESHAPE(REAL(pRegion%mixt%dv(DV_MIXT_SOUN,1:nCells),8),(/1,nCells/)))
    END IF ! global%postWriteSoundSpeed
    IF ( global%postWriteEnergy ) THEN
      CALL AddPostField(fields,nFields,'Energy', &
        RESHAPE(pRegion%mixt%cv(CV_MIXT_ENER,1:nCells)/ &
                pRegion%mixt%cv(CV_MIXT_DENS,1:nCells),(/1,nCells/)))
    END IF ! global%postWriteEnergy
#ifdef PICL
    IF ( global%piclUsed .EQV. .TRUE. ) THEN
      IF ( global%postWritePhi ) THEN
        CALL AddPostField(fields,nFields,'Particle Volume Fraction', &
                          RESHAPE(phi,(/1,nCells/)))
      END IF ! global%postWritePhi
      IF ( global%postWriteGradRho ) THEN
        DO c = 1,3
          work(c,:) = pRegion%mixt%piclgradRhog(c,1,1:nCells)
        END DO ! c
        CALL AddPostField(fields,nFields,'Gradient Gas Density',work)
      END IF ! global%postWriteGradRho
    END IF ! global%piclUsed
#endif

! Cells of every type, in Rocflu's cell order (VTK tetra 10, hexahedron 12,
! wedge 13, pyramid 14; vertices in Rocflu's order)
    ALLOCATE(cellTypes(nCells),offsets(nCells))
    nConn = 0
    DO icg = 1,nCells
      SELECT CASE ( pGrid%cellGlob2Loc(1,icg) )
        CASE ( CELL_TYPE_TET )
          nConn = nConn + 4
        CASE ( CELL_TYPE_HEX )
          nConn = nConn + 8
        CASE ( CELL_TYPE_PRI )
          nConn = nConn + 6
        CASE ( CELL_TYPE_PYR )
          nConn = nConn + 5
      END SELECT ! pGrid%cellGlob2Loc
    END DO ! icg
    ALLOCATE(connectivity(nConn))
    nConn = 0
    DO icg = 1,nCells
      iCellType = pGrid%cellGlob2Loc(1,icg)
      icl       = pGrid%cellGlob2Loc(2,icg)
      SELECT CASE ( iCellType )
        CASE ( CELL_TYPE_TET )
          nVertCell = 4
          cellTypes(icg) = 10
          connectivity(nConn+1:nConn+4) = pGrid%tet2v(1:4,icl) - 1
        CASE ( CELL_TYPE_HEX )
          nVertCell = 8
          cellTypes(icg) = 12
          connectivity(nConn+1:nConn+8) = pGrid%hex2v(1:8,icl) - 1
        CASE ( CELL_TYPE_PRI )
          nVertCell = 6
          cellTypes(icg) = 13
          connectivity(nConn+1:nConn+6) = pGrid%pri2v(1:6,icl) - 1
        CASE ( CELL_TYPE_PYR )
          nVertCell = 5
          cellTypes(icg) = 14
          connectivity(nConn+1:nConn+5) = pGrid%pyr2v(1:5,icl) - 1
      END SELECT ! iCellType
      nConn = nConn + nVertCell
      offsets(icg) = nConn
    END DO ! icg

    fileName = StepFileName('post/'//TRIM(dirName)//'/',stepNumber)
    CALL WriteUnstructuredVTU(TRIM(fileName),RFLU_FluidDumpTime(global, &
                              stepNumber),stepNumber, &
                              REAL(pGrid%xyz(1:3,1:pGrid%nVert),8), &
                              connectivity,offsets,cellTypes,fields,nFields, &
                              .FALSE.)

  END SUBROUTINE WriteRegion

! ******************************************************************************
! Particle and bin files of a step (fileWritten: whether there were any)
! ******************************************************************************

  SUBROUTINE WriteParticlesAndBins(stepNumber,fileWritten)

    INTEGER, INTENT(IN) :: stepNumber
    LOGICAL, INTENT(OUT) :: fileWritten

    INTEGER :: i,ib,ix,iy,iz,k,nBinCells,nFields,nx,ny,nz
    INTEGER(1), ALLOCATABLE :: cellTypes(:)
    INTEGER(8) :: n
    INTEGER(8), ALLOCATABLE :: connectivity(:),offsets(:)
    REAL(8), ALLOCATABLE :: points(:,:),values(:,:)
    TYPE(t_particleDump) :: dump
    TYPE(t_postField), ALLOCATABLE :: fields(:)

    fileWritten = .FALSE.
    IF ( global%piclUsed .EQV. .FALSE. ) RETURN
    INQUIRE(FILE=TRIM(RFLU_ParticleDumpName(global,stepNumber)), &
            EXIST=fileWritten)
    IF ( .NOT. fileWritten ) RETURN

    CALL ReadParticleDump(TRIM(RFLU_ParticleDumpName(global,stepNumber)),dump)

! Particles: one vertex each, the fields of the output-field table
    n = dump%nParticles
    ALLOCATE(cellTypes(n),offsets(n),connectivity(n),values(1,n))
    cellTypes = 1
    DO i = 1,INT(n)
      connectivity(i) = i - 1
      offsets(i)      = i
    END DO ! i
    nFields = 0
    DO k = 1,dump%nOut
      SELECT CASE ( dump%outSource(k) )
        CASE ( OUT_SOURCE_Y )
          values(1,:) = dump%y(dump%outSlot(k),:)
        CASE ( OUT_SOURCE_RPROP )
          values(1,:) = dump%rprop(dump%outSlot(k),:)
        CASE ( OUT_SOURCE_TAG )
          values(1,:) = REAL(dump%iprop(dump%outSlot(k),:),8)
      END SELECT ! dump%outSource
      CALL AddPostField(fields,nFields,TRIM(dump%outName(k)),values)
    END DO ! k
    CALL WriteUnstructuredVTU(TRIM(StepFileName('post/particles/',stepNumber)), &
                              dump%time,stepNumber,dump%y(1:3,:), &
                              connectivity,offsets,cellTypes,fields,nFields, &
                              .TRUE.)
    DEALLOCATE(cellTypes,offsets,connectivity,values,fields)

! Bins: the box split into nx x ny x nz hexahedra; bin (ix,iy,iz) belongs to
! rank ix + nx*(iy + ny*iz) (ppiclf_comm_FindParticle)
    nx = dump%nBins(1)
    ny = dump%nBins(2)
    nz = dump%nBins(3)
    nBinCells = nx*ny*nz
    ALLOCATE(points(3,(nx+1)*(ny+1)*(nz+1)),cellTypes(nBinCells), &
             offsets(nBinCells),connectivity(8*nBinCells),values(1,nBinCells))
    DO iz = 0,nz
      DO iy = 0,ny
        DO ix = 0,nx
          i = 1 + ix + (nx+1)*(iy + (ny+1)*iz)
          points(1,i) = dump%binBox(1) + ix*dump%binSize(1)
          points(2,i) = dump%binBox(3) + iy*dump%binSize(2)
          points(3,i) = dump%binBox(5) + iz*dump%binSize(3)
        END DO ! ix
      END DO ! iy
    END DO ! iz
    ib = 0
    DO iz = 0,nz-1
      DO iy = 0,ny-1
        DO ix = 0,nx-1
          ib = ib + 1
          k = ix + (nx+1)*(iy + (ny+1)*iz)
          connectivity(8*ib-7:8*ib) = (/ k, k+1, k+1+(nx+1), k+(nx+1), &
            k+(nx+1)*(ny+1), k+1+(nx+1)*(ny+1), k+1+(nx+1)+(nx+1)*(ny+1), &
            k+(nx+1)+(nx+1)*(ny+1) /)
          offsets(ib)  = 8*ib
          values(1,ib) = REAL(ix + nx*(iy + ny*iz),8)
        END DO ! ix
      END DO ! iy
    END DO ! iz
    cellTypes = 12
    nFields = 0
    CALL AddPostField(fields,nFields,'rank',values)
    CALL WriteUnstructuredVTU(TRIM(StepFileName('post/bins/',stepNumber)), &
                              dump%time,stepNumber,points,connectivity, &
                              offsets,cellTypes,fields,nFields,.FALSE.)

  END SUBROUTINE WriteParticlesAndBins

END SUBROUTINE rflupost
