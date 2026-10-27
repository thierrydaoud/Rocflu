!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_BinCounts(length,minLen,nRanks,nBins)
!
! Number of bins in x, y, z for a bin box of the given lengths: start
! from one bin and keep splitting the direction with the longest bins,
! among the directions where one more bin keeps the total within nRanks
! and the bins at least minLen long. Every step keeps the layout valid.
! No MPI and no global data (tested by tests/bin_counts).
!
      IMPLICIT NONE
!
! Input:
!
      REAL*8    length(3), minLen(3)
      INTEGER*4 nRanks
!
! Output:
!
      INTEGER*4 nBins(3)
!
! Internal:
!
      INTEGER*4 i, longest, total
      REAL*8    maxLen
!
      nBins(1) = 1
      nBins(2) = 1
      nBins(3) = 1
      DO
        total   = nBins(1)*nBins(2)*nBins(3)
        longest = 0
        maxLen  = 0.0D0
        DO i = 1,3
          IF((total/nBins(i))*(nBins(i)+1) .GT. nRanks) CYCLE
          IF(length(i)/(nBins(i)+1) .LT. minLen(i)) CYCLE
          IF(length(i)/nBins(i) .GT. maxLen) THEN
            maxLen  = length(i)/nBins(i)
            longest = i
          END IF
        END DO
        IF(longest .EQ. 0) EXIT
        nBins(longest) = nBins(longest) + 1
      END DO

      RETURN
      END
