! ******************************************************************************
!
! Purpose: Read a particle dump (restart_data/particles_<step>.dat, layout in
!   modflu/RFLU_ModRestartData.F90) for the post program.
!
! Description: The whole file is read on one rank: the particles of all ranks
!   one after the other, the output-field table, the bin layout and the
!   particle volume fraction projected to each rank's fluid cells (the first
!   projected field, P_JPHIP). Versions 3 and 4 are read.
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
                  nBins4(3),nOut4,source4,slot4,nSaved
    INTEGER(4), ALLOCATABLE :: cells(:),saved(:),sharedTags(:),storedF(:), &
                               storedH(:),storedI(:),storedR(:),storedY(:)
    REAL(8), ALLOCATABLE :: props(:,:),sharedF(:),sharedH(:),sharedR(:), &
                            sharedY(:)
    INTEGER(8) :: first,n
    INTEGER(8), ALLOCATABLE :: counts(:),cellCounts(:),nonzeroCounts(:)
    REAL(8), ALLOCATABLE :: skipped(:),projected(:,:)

    OPEN(NEWUNIT=iFile,FILE=fileName,FORM='UNFORMATTED',ACCESS='STREAM', &
         STATUS='OLD',IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      WRITE(*,'(A,A)') 'rflupost: cannot read ',TRIM(fileName)
      STOP 1
    END IF ! errorFlag

    READ(iFile) magic,version,dump%step,dump%time,nRanks4,dump%nParticles, &
                lrs,lrp,lip,lrp3,lrpPro,nHist,nTimeBH
    IF ( magic /= 'PPICLDMP' .OR. (version /= 3 .AND. version /= 4) ) THEN
      WRITE(*,'(A,A)') 'rflupost: not a version 3 or 4 particle dump: ', &
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

! The rprop slots in the file (needed by a restart, or output)
    READ(iFile) nSaved
    ALLOCATE(saved(nSaved))
    IF ( nSaved > 0 ) READ(iFile) saved

! Version 4: which columns are per particle (stored 1) and the value of the
! others; version 3: every column per particle
    ALLOCATE(storedY(lrs),storedR(nSaved),storedI(lip),storedH(lrp3), &
             storedF(lrpPro),sharedY(lrs),sharedR(nSaved),sharedTags(lip), &
             sharedH(lrp3),sharedF(lrpPro))
    storedY = 1
    storedR = 1
    storedI = 1
    storedH = 1
    storedF = 1
    IF ( version == 4 ) THEN
      READ(iFile) storedY,sharedY
      READ(iFile) storedR,sharedR
      READ(iFile) storedI,sharedTags
      READ(iFile) storedH,sharedH
      READ(iFile) storedF,sharedF
    END IF ! version

    ALLOCATE(counts(nRanks4),cellCounts(nRanks4),nonzeroCounts(nRanks4))
    IF ( version == 4 ) THEN
      READ(iFile) counts,cellCounts,nonzeroCounts
    ELSE
      READ(iFile) counts,cellCounts
    END IF ! version

    ALLOCATE(dump%y(lrs,dump%nParticles),dump%rprop(lrp,dump%nParticles), &
             dump%iprop(lip,dump%nParticles),dump%ranks(nRanks4))
    dump%rprop = 0.0D0
    first = 1
    DO iRank = 1,nRanks4
      n = counts(iRank)
      IF ( n > 0 ) THEN
        CALL ReadColumns(iFile,storedY,sharedY,dump%y(:,first:first+n-1))
        ALLOCATE(props(nSaved,n))
        CALL ReadColumns(iFile,storedR,sharedR,props)
        dump%rprop(saved,first:first+n-1) = props
        DEALLOCATE(props)
        CALL ReadTags(iFile,storedI,sharedTags,dump%iprop(:,first:first+n-1))
        ALLOCATE(skipped((COUNT(storedH == 1) + COUNT(storedF == 1))*n))
        IF ( SIZE(skipped) > 0 ) READ(iFile) skipped   ! history, feedback
        DEALLOCATE(skipped)
      END IF ! n
      first = first + n

! Version 3: every cell of the rank; 4: the cells not 0, then their values
      ALLOCATE(dump%ranks(iRank)%phi(cellCounts(iRank)))
      dump%ranks(iRank)%phi = 0.0D0
      IF ( version == 3 ) THEN
        ALLOCATE(projected(cellCounts(iRank),lrpPro))
        IF ( SIZE(projected) > 0 ) READ(iFile) projected
        IF ( (cellCounts(iRank) > 0) .AND. (lrpPro > 0) ) THEN
          dump%ranks(iRank)%phi = projected(:,1)
        END IF ! cellCounts
      ELSE
        ALLOCATE(cells(nonzeroCounts(iRank)), &
                 projected(nonzeroCounts(iRank),lrpPro))
        IF ( SIZE(cells) > 0 ) READ(iFile) cells
        IF ( SIZE(projected) > 0 ) READ(iFile) projected
        IF ( (SIZE(cells) > 0) .AND. (lrpPro > 0) ) THEN
          dump%ranks(iRank)%phi(cells) = projected(:,1)
        END IF ! SIZE
        DEALLOCATE(cells)
      END IF ! version
      DEALLOCATE(projected)
    END DO ! iRank

    CLOSE(iFile)

  END SUBROUTINE ReadParticleDump

! The stored columns of a (columns x particles) from the file, the others set
! to their shared value
  SUBROUTINE ReadColumns(iFile,stored,shared,a)

    INTEGER, INTENT(IN) :: iFile
    INTEGER(4), INTENT(IN) :: stored(:)
    REAL(8), INTENT(IN) :: shared(:)
    REAL(8), INTENT(INOUT) :: a(:,:)

    INTEGER :: j,k
    REAL(8), ALLOCATABLE :: packed(:,:)

    ALLOCATE(packed(COUNT(stored == 1),SIZE(a,2)))
    IF ( SIZE(packed) > 0 ) READ(iFile) packed
    k = 0
    DO j = 1,SIZE(stored)
      IF ( stored(j) == 1 ) THEN
        k = k + 1
        a(j,:) = packed(k,:)
      ELSE
        a(j,:) = shared(j)
      END IF ! stored
    END DO ! j

  END SUBROUTINE ReadColumns

! ReadColumns for the integer columns (iprop)
  SUBROUTINE ReadTags(iFile,stored,shared,a)

    INTEGER, INTENT(IN) :: iFile
    INTEGER(4), INTENT(IN) :: stored(:),shared(:)
    INTEGER(4), INTENT(INOUT) :: a(:,:)

    INTEGER :: j,k
    INTEGER(4), ALLOCATABLE :: packed(:,:)

    ALLOCATE(packed(COUNT(stored == 1),SIZE(a,2)))
    IF ( SIZE(packed) > 0 ) READ(iFile) packed
    k = 0
    DO j = 1,SIZE(stored)
      IF ( stored(j) == 1 ) THEN
        k = k + 1
        a(j,:) = packed(k,:)
      ELSE
        a(j,:) = shared(j)
      END IF ! stored
    END DO ! j

  END SUBROUTINE ReadTags

END MODULE RFLU_ModPostParticles
