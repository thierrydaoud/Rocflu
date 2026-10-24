! ******************************************************************************
!
! Purpose: Print the Setup box of the run log: grid, particles, models,
!   ppiclF settings and time stepping, two items per line.
!
! Description: Called once at the end of the initialization, on every rank
!   (the cell counts are summed over the ranks); rank 0 prints. Particle
!   items come from global%picl*, filled by PICL_InitSolver.
!
! Input:
!   regions     Data of the local regions
!
! Output: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PrintSetup(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: MAX_ITEMS = 24

  CHARACTER(15) :: labels(MAX_ITEMS)
  CHARACTER(35) :: values(MAX_ITEMS)
  CHARACTER(105) :: titleLine
  CHARACTER(35) :: stepText
  INTEGER :: errorFlag,iReg,nItems
  INTEGER :: cellCounts(5),cellCountsGlobal(5)
  REAL(RFREAL) :: endTime
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

! Cells of each type, summed over the ranks: all, tetrahedra, hexahedra,
! prisms, pyramids
  cellCounts = 0
  DO iReg = 1,global%nRegionsLocal
    pGrid => regions(iReg)%grid
    cellCounts = cellCounts + (/ pGrid%nCells,pGrid%nTets,pGrid%nHexs, &
                                 pGrid%nPris,pGrid%nPyrs /)
  END DO ! iReg
  CALL MPI_Allreduce(cellCounts,cellCountsGlobal,5,MPI_INTEGER,MPI_SUM, &
                     global%mpiComm,errorFlag)

  IF ( global%myProcid /= MASTERPROC ) RETURN

  WRITE(STDOUT,'(A)') ''
  CALL PrintRule
  titleLine = 'Setup '//TRIM(global%casename)//':'
  WRITE(STDOUT,'(3A)') '| ',titleLine,' |'
  CALL PrintRule

! ==============================================================================
! Case: grid and particles
! ==============================================================================

  nItems = 0
  CALL AddItem('Cells',TRIM(IntText(cellCountsGlobal(1)))//' ('// &
               TRIM(CellTypesText(cellCountsGlobal(2:5)))//')')
  CALL AddItem('Regions',TRIM(IntText(global%nRegions))//' on '// &
               TRIM(IntText(global%nProcs))//' ranks')
  IF ( global%piclUsed ) THEN
    CALL AddItem('Particles',IntText(global%piclNParticles))
    IF ( global%piclDpMin == global%piclDpMax ) THEN
      CALL AddItem('Diameter',TRIM(RealText(global%piclDpMax))//' m')
    ELSE
      CALL AddItem('Diameter',TRIM(RealText(global%piclDpMin))//' to '// &
                   TRIM(RealText(global%piclDpMax))//' m')
    END IF ! global%piclDpMin
    CALL AddItem('Density',TRIM(RealText(global%piclDensity))//' kg/m3')
    CALL AddItem('Cp',TRIM(RealText(global%piclSpecificHeat))//' J/(kg K)')
    CALL AddItem('Temperature',TRIM(RealText(global%piclTemp))//' K')
    CALL AddItem('Periodic',PeriodicText())
  END IF ! global%piclUsed
  CALL PrintGroup

! ==============================================================================
! Particle models
! ==============================================================================

  IF ( global%piclUsed ) THEN
    CALL AddItem('Motion',PickName(global%piclVelocityModel, &
                 (/ 'frozen (VELOCITYMODEL 0)','moving                  ' /)))
    CALL AddItem('QS force',FlagName(global%piclQsFlag, &
                 (/ 'Gidaspow','Parmar  ','Osnes   ' /)))
    CALL AddItem('QS fluct.',FlagName(global%piclQsFluctFlag, &
                 (/ 'Lattanzi','Osnes   ' /)))
    CALL AddItem('Press. grad.',OnOff(global%piclPgFlag))
    CALL AddItem('Added mass',FlagName(global%piclAmFlag, &
                 (/ 'Parmar','Briney' /)))
    IF ( global%piclViscousUnsteady > 0 ) THEN
      CALL AddItem('Visc. unst.',TRIM(FlagName(global%piclViscousUnsteady, &
                   (/ 'trapezoidal, uniform dt ','trapezoidal, variable dt', &
                      'modified trapezoidal    ','Hinsberg                ' /)))// &
                   ', '//TRIM(IntText(global%piclVUSteps))//' steps')
    ELSE
      CALL AddItem('Visc. unst.','off')
    END IF ! global%piclViscousUnsteady
    CALL AddItem('Collisions',FlagName(global%piclCollisionFlag, &
                 (/ 'normal                ','+ tangential, spin    ', &
                    '+ rolling, hydro torq.','+ Saffman, Magnus lift' /)))
    CALL AddItem('Heat transfer',FlagName(global%piclHeatTransferFlag, &
                 (/ 'Stokes       ','Ranz-Marshall','Gunn         ', &
                    'Fox          ','Sun          ' /)))
    CALL AddItem('DU heat tr.',OnOff(global%piclHTUnsteadyFlag))
    CALL AddItem('Feedback',PickName(global%piclFeedbackFlag, &
                 (/ 'one-way','two-way' /)))
    CALL PrintGroup
  END IF ! global%piclUsed

! ==============================================================================
! Numerics: ppiclF and time stepping
! ==============================================================================

  IF ( global%piclUsed ) THEN
    CALL AddItem('Bins',TRIM(IntText(global%piclNBins(1)))//' x '// &
                 TRIM(IntText(global%piclNBins(2)))//' x '// &
                 TRIM(IntText(global%piclNBins(3))))
    CALL AddItem('Filter width',FilterText())
    CALL AddItem('Neighb. width', &
                 TRIM(RealText(global%piclNeighborWidthUsed))//' m')
    CALL AddItem('Neighb. search',PickName(global%piclSBNearFlag, &
                 (/ 'Zwick         ','Avery sub-bins' /)))
  END IF ! global%piclUsed

  endTime = global%restartTime &
          + MIN(global%maxTime - global%restartTime,global%runTime)
  CALL AddItem('CFL',RealText(regions(1)%mixtInput%cfl))
  CALL AddItem('dt max',TRIM(RealText(global%dtImposed))//' s')
  CALL AddItem('Start time',TRIM(RealText(global%restartTime))//' s')
  CALL AddItem('End time',TRIM(RealText(endTime))//' s')
  IF ( global%currentStep == 0 ) THEN
    CALL AddItem('Start step','0 (new run)')
  ELSE
    WRITE(stepText,'(I0,A)') global%currentStep,' (restart)'
    CALL AddItem('Start step',TRIM(stepText))
  END IF ! global%currentStep
  IF ( global%nFilesRemoved > 0 ) THEN
    WRITE(stepText,'(I0,A)') global%nFilesRemoved,' of an earlier run'
    CALL AddItem('Files removed',TRIM(stepText))
  END IF ! global%nFilesRemoved
  CALL AddItem('Print every',TRIM(RealText(global%printTime))//' s')
  CALL AddItem('Write every',TRIM(RealText(global%writeTime))//' s')
  CALL PrintGroup

  WRITE(STDOUT,'(A)') ''

CONTAINS

  SUBROUTINE PrintRule
    WRITE(STDOUT,'(3A)') '+',REPEAT('-',107),'+'
  END SUBROUTINE PrintRule

  SUBROUTINE AddItem(label,value)
    CHARACTER(*), INTENT(IN) :: label,value
    nItems = nItems + 1
    labels(nItems) = label
    values(nItems) = value
  END SUBROUTINE AddItem

! Items two per line ('* label' padded to 15, value padded to 35), then a rule
  SUBROUTINE PrintGroup
    CHARACTER(52) :: leftColumn
    CHARACTER(53) :: rightColumn
    INTEGER :: i

    DO i = 1,nItems,2
      leftColumn = '* '//labels(i)//values(i)
      rightColumn = ''
      IF ( i < nItems ) rightColumn = '* '//labels(i+1)//values(i+1)
      WRITE(STDOUT,'(4A)') '| ',leftColumn,rightColumn,' |'
    END DO ! i
    IF ( nItems > 0 ) CALL PrintRule
    nItems = 0
  END SUBROUTINE PrintGroup

  FUNCTION IntText(n) RESULT(text)
    INTEGER, INTENT(IN) :: n
    CHARACTER(16) :: text
    WRITE(text,'(I0)') n
  END FUNCTION IntText

! Plain decimals for 0.01 <= |x| < 1e5 (trailing zeros dropped), else E format
  FUNCTION RealText(x) RESULT(text)
    REAL(RFREAL), INTENT(IN) :: x
    CHARACTER(16) :: text
    INTEGER :: last

    IF ( x == 0.0_RFREAL ) THEN
      text = '0'
    ELSE IF ( ABS(x) >= 1.0E-2_RFREAL .AND. ABS(x) < 1.0E5_RFREAL ) THEN
      WRITE(text,'(F16.4)') x
      text = ADJUSTL(text)
      last = LEN_TRIM(text)
      DO WHILE ( text(last:last) == '0' )
        last = last - 1
      END DO
      IF ( text(last:last) == '.' ) last = last - 1
      text = text(1:last)
    ELSE
      WRITE(text,'(ES10.3)') x
      text = ADJUSTL(text)
    END IF ! x
  END FUNCTION RealText

  FUNCTION OnOff(flag) RESULT(text)
    INTEGER, INTENT(IN) :: flag
    CHARACTER(3) :: text
    text = 'off'
    IF ( flag > 0 ) text = 'on'
  END FUNCTION OnOff

! names(flag+1) for a 0/1 setting
  FUNCTION PickName(flag,names) RESULT(text)
    INTEGER, INTENT(IN) :: flag
    CHARACTER(*), INTENT(IN) :: names(:)
    CHARACTER(35) :: text
    text = names(MAX(1,MIN(SIZE(names),flag+1)))
  END FUNCTION PickName

! 'off' for flag 0, else 'names(flag) (flag)'
  FUNCTION FlagName(flag,names) RESULT(text)
    INTEGER, INTENT(IN) :: flag
    CHARACTER(*), INTENT(IN) :: names(:)
    CHARACTER(35) :: text
    IF ( flag <= 0 ) THEN
      text = 'off'
    ELSE IF ( flag > SIZE(names) ) THEN
      text = 'flag '//TRIM(IntText(flag))
    ELSE
      text = TRIM(names(flag))//' ('//TRIM(IntText(flag))//')'
    END IF ! flag
  END FUNCTION FlagName

  FUNCTION CellTypesText(counts) RESULT(text)
    INTEGER, INTENT(IN) :: counts(4)
    CHARACTER(CHRLEN) :: text
    CHARACTER(10), PARAMETER :: typeNames(4) = (/ 'tetrahedra','hexahedra ', &
                                                  'prisms    ','pyramids  ' /)
    INTEGER :: i,nTypes

    nTypes = COUNT(counts > 0)
    text = ''
    DO i = 1,4
      IF ( counts(i) == 0 ) CYCLE
      IF ( LEN_TRIM(text) > 0 ) text = TRIM(text)//', '
      IF ( nTypes == 1 ) THEN
        text = TRIM(text)//TRIM(typeNames(i))
      ELSE
        text = TRIM(text)//TRIM(IntText(counts(i)))//' '//TRIM(typeNames(i))
      END IF ! nTypes
    END DO ! i
  END FUNCTION CellTypesText

  FUNCTION PeriodicText() RESULT(text)
    CHARACTER(16) :: text
    text = ''
    IF ( global%piclPeriodicXFlag == 1 ) text = TRIM(text)//', x'
    IF ( global%piclPeriodicYFlag == 1 ) text = TRIM(text)//', y'
    IF ( global%piclPeriodicZFlag == 1 ) text = TRIM(text)//', z'
    IF ( LEN_TRIM(text) == 0 ) THEN
      text = 'none'
    ELSE
      text = text(3:)
    END IF ! LEN_TRIM
  END FUNCTION PeriodicText

! One value if the filter width is the same in x, y and z, else all three
  FUNCTION FilterText() RESULT(text)
    CHARACTER(35) :: text
    REAL(RFREAL) :: widths(3)
    widths = global%piclFilterWidth
    IF ( MAXVAL(widths) - MINVAL(widths) <= 1.0E-6_RFREAL*MAXVAL(widths) ) THEN
      text = TRIM(RealText(widths(1)))//' m'
    ELSE
      text = TRIM(RealText(widths(1)))//' x '//TRIM(RealText(widths(2)))// &
             ' x '//TRIM(RealText(widths(3)))//' m'
    END IF ! MAXVAL
  END FUNCTION FilterText

END SUBROUTINE RFLU_PrintSetup
