! ******************************************************************************
!
! Purpose: Read a particle dump (restart_data/particles_<step>.dat, layout in
!   modflu/RFLU_ModRestartData.F90) for the post program.
!
! Description: The whole file is read on one rank: the particles of all ranks
!   one after the other, the output-field table, the bin layout and the
!   particle volume fraction projected to each rank's fluid cells (the first
!   projected field, P_JPHIP).
!
! ******************************************************************************

MODULE RFLU_ModPostParticles

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_particleDump, ReadParticleDump

  INTEGER, PARAMETER, PUBLIC :: OUT_SOURCE_Y = 1, &
                                OUT_SOURCE_RPROP = 2, &
                                OUT_SOURCE_TAG = 3

! The volume fraction projected to one rank's fluid cells
  TYPE t_rankCells
    REAL(8), ALLOCATABLE :: phi(:)
  END TYPE t_rankCells

  TYPE t_particleDump
    INTEGER(8) :: step,nParticles
    REAL(8) :: time
    INTEGER :: nRanks,nBins(3)
    REAL(8) :: binBox(6),binSize(3)
! Output fields: source (OUT_SOURCE_*), slot and name
    INTEGER :: nOut
    INTEGER, ALLOCATABLE :: outSource(:),outSlot(:)
    CHARACTER(24), ALLOCATABLE :: outName(:)
! Per particle, in rank order
    REAL(8), ALLOCATABLE :: y(:,:),rprop(:,:)
    INTEGER(4), ALLOCATABLE :: iprop(:,:)
    TYPE(t_rankCells), ALLOCATABLE :: ranks(:)
  END TYPE t_particleDump

CONTAINS

  SUBROUTINE ReadParticleDump(fileName,dump)

    CHARACTER(*), INTENT(IN) :: fileName
    TYPE(t_particleDump), INTENT(OUT) :: dump

    CHARACTER(8) :: magic
    INTEGER :: errorFlag,iFile,iRank,k
    INTEGER(4) :: version,nRanks4,lrs,lrp,lip,lrp3,lrpPro,nHist,nTimeBH, &
                  nBins4(3),nOut4,source4,slot4
    INTEGER(8) :: first,n
    INTEGER(8), ALLOCATABLE :: counts(:),cellCounts(:)
    REAL(8), ALLOCATABLE :: skipped(:),projected(:,:)

    OPEN(NEWUNIT=iFile,FILE=fileName,FORM='UNFORMATTED',ACCESS='STREAM', &
         STATUS='OLD',IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      WRITE(*,'(A,A)') 'rflupost: cannot read ',TRIM(fileName)
      STOP 1
    END IF ! errorFlag

    READ(iFile) magic,version,dump%step,dump%time,nRanks4,dump%nParticles, &
                lrs,lrp,lip,lrp3,lrpPro,nHist,nTimeBH
    IF ( magic /= 'PPICLDMP' .OR. version /= 2 ) THEN
      WRITE(*,'(A,A)') 'rflupost: not a version 2 particle dump: ', &
                       TRIM(fileName)
      STOP 1
    END IF ! magic
    dump%nRanks = nRanks4

    IF ( nHist > 0 ) THEN
      ALLOCATE(skipped(nHist))
      READ(iFile) skipped
      DEALLOCATE(skipped)
    END IF ! nHist
    READ(iFile) nBins4,dump%binBox,dump%binSize
    dump%nBins = nBins4
    READ(iFile) nOut4
    dump%nOut = nOut4
    ALLOCATE(dump%outSource(nOut4),dump%outSlot(nOut4),dump%outName(nOut4))
    DO k = 1,nOut4
      READ(iFile) source4,slot4,dump%outName(k)
      dump%outSource(k) = source4
      dump%outSlot(k)   = slot4
    END DO ! k

    ALLOCATE(counts(nRanks4),cellCounts(nRanks4))
    READ(iFile) counts,cellCounts

    ALLOCATE(dump%y(lrs,dump%nParticles),dump%rprop(lrp,dump%nParticles), &
             dump%iprop(lip,dump%nParticles),dump%ranks(nRanks4))
    first = 1
    DO iRank = 1,nRanks4
      n = counts(iRank)
      IF ( n > 0 ) THEN
        READ(iFile) dump%y(:,first:first+n-1)
        READ(iFile) dump%rprop(:,first:first+n-1)
        READ(iFile) dump%iprop(:,first:first+n-1)
        ALLOCATE(skipped((lrp3 + lrpPro)*n))
        IF ( SIZE(skipped) > 0 ) READ(iFile) skipped   ! history, feedback
        DEALLOCATE(skipped)
      END IF ! n
      first = first + n

      ALLOCATE(projected(cellCounts(iRank),lrpPro))
      IF ( SIZE(projected) > 0 ) READ(iFile) projected
      ALLOCATE(dump%ranks(iRank)%phi(cellCounts(iRank)))
      IF ( (cellCounts(iRank) > 0) .AND. (lrpPro > 0) ) THEN
        dump%ranks(iRank)%phi = projected(:,1)
      END IF ! cellCounts
      DEALLOCATE(projected)
    END DO ! iRank

    CLOSE(iFile)

  END SUBROUTINE ReadParticleDump

END MODULE RFLU_ModPostParticles
