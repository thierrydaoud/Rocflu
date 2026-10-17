! ******************************************************************************
!
! Purpose: Write the particle dump of this step, restart_data/particles_<step>.dat
!   (layout in modflu/RFLU_ModRestartData.F90).
!
! Description: Collective over all ranks: ppiclf_io_WriteRestart writes one file
!   for all ranks; rank 0 then gives it its final name. In its own file (as
!   PICL_WriteVTU) so that only programs that link ppiclF pull it in.
!
! Input:
!   pRegion     Region data (for global)
!
! Output: None.
!
! ******************************************************************************

SUBROUTINE PICL_WriteRestartData(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModMPI
  USE RFLU_ModRestartData, ONLY: RFLU_CreateRestartDataDir, &
                                 RFLU_FinishRestartFile

  IMPLICIT NONE

  TYPE(t_region), POINTER :: pRegion

  CHARACTER(CHRLEN) :: dirName,fileName
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

  global => pRegion%global

  CALL RegisterFunction(global,'PICL_WriteRestartData',__FILE__)

  CALL RFLU_CreateRestartDataDir(global,dirName)
  WRITE(fileName,'(A,I8.8,A)') TRIM(dirName)//'/particles_', &
        global%currentStep,'.dat'

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,A)') SOLVER_NAME,'Writing particle dump', &
                                  TRIM(fileName)
  END IF ! global%verbLevel

  CALL ppiclf_io_WriteRestart(TRIM(fileName)//'.tmp', &
                              INT(global%currentStep,8), &
                              REAL(global%currentTime,8))

! Every rank has closed the file; rank 0 gives it its final name
  CALL MPI_Barrier(global%mpiComm,errorFlag)
  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_FinishRestartFile(global,fileName)
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE PICL_WriteRestartData
