!
! Unit test of ppiclf_comm_BinCounts (libpicl/source/ppiclf_bincounts.f),
! run by run_tests.py as unit/bin-counts.
!   1. Golden table: named bin boxes and the bin counts they must get.
!   2. Scan: every layout over many boxes and rank counts is valid
!      (total <= ranks, bins at least the minimum length, >= 1 bin).
! Prints the failures; exit code 1 if any.
!
      PROGRAM bin_counts_test
      IMPLICIT NONE
      INTEGER*4 NGOLD
      PARAMETER (NGOLD = 10)
      REAL*8    minLen(3), length(3), golden(3,NGOLD), thick(3)
      INTEGER*4 nRanks(NGOLD), expected(3,NGOLD), nBins(3), k, ix, iy,
     >          iz, ip, np, nFail, nCases
      CHARACTER*24 label(NGOLD)

      ! Filter and neighbour widths of the wedge case
      minLen = (/8.703D-4, 3.782D-4, 4.0D-4/)

      label(1)  = 'wedge bed at t = 0'
      golden(:,1)  = (/5.848D-3, 4.791D-3, 1.368D-3/)
      nRanks(1)    = 32
      expected(:,1)  = (/6, 5, 1/)
      label(2)  = 'wedge bed at 11 us'
      golden(:,2)  = (/32.631D-3, 5.673D-3, 1.379D-3/)
      nRanks(2)    = 32
      expected(:,2)  = (/16, 2, 1/)
      label(3)  = 'wedge bed at the crash'
      golden(:,3)  = (/36.7D-3, 5.81D-3, 1.378D-3/)
      nRanks(3)    = 32
      expected(:,3)  = (/16, 2, 1/)
      label(4)  = 'same, 64 ranks'
      golden(:,4)  = (/36.7D-3, 5.81D-3, 1.378D-3/)
      nRanks(4)    = 64
      expected(:,4)  = (/21, 3, 1/)
      label(5)  = 'cube'
      golden(:,5)  = (/20.0D-3, 20.0D-3, 20.0D-3/)
      nRanks(5)    = 64
      expected(:,5)  = (/4, 4, 4/)
      label(6)  = 'slab'
      golden(:,6)  = (/40.0D-3, 40.0D-3, 4.0D-3/)
      nRanks(6)    = 64
      expected(:,6)  = (/8, 8, 1/)
      label(7)  = 'long bed'
      golden(:,7)  = (/100.0D-3, 5.0D-3, 5.0D-3/)
      nRanks(7)    = 128
      expected(:,7)  = (/32, 2, 2/)
      label(8)  = 'one rank'
      golden(:,8)  = (/36.7D-3, 5.81D-3, 1.378D-3/)
      nRanks(8)    = 1
      expected(:,8)  = (/1, 1, 1/)
      label(9)  = 'box near the minimum'
      golden(:,9)  = (/0.9D-3, 0.4D-3, 0.41D-3/)
      nRanks(9)    = 32
      expected(:,9)  = (/1, 1, 1/)
      label(10) = 'prime rank count'
      golden(:,10) = (/20.0D-3, 20.0D-3, 20.0D-3/)
      nRanks(10)   = 7
      expected(:,10) = (/3, 2, 1/)

      nFail = 0
      DO k = 1,NGOLD
        CALL ppiclf_comm_BinCounts(golden(:,k),minLen,nRanks(k),nBins)
        IF (ANY(nBins .NE. expected(:,k))) THEN
          nFail = nFail + 1
          PRINT '(A,A,A,I4,A,3I4,A,3I4)', 'golden ', TRIM(label(k)),
     >      ', ', nRanks(k), ' ranks: got', nBins, ', expected',
     >      expected(:,k)
        END IF
      END DO

      ! Lengths from 2 mm to about 0.8 m in x, 0.8 mm to 5.6 mm in y
      thick = (/1.4D-3, 4.0D-3, 20.0D-3/)
      nCases = 0
      DO ip = 0,8
        np = 2**ip + MOD(ip,3)
        DO ix = 0,300
          DO iy = 0,40
            DO iz = 1,3
              length = (/2.0D-3*1.02D0**ix, 0.8D-3*1.05D0**iy,
     >                   thick(iz)/)
              CALL ppiclf_comm_BinCounts(length,minLen,np,nBins)
              nCases = nCases + 1
              IF (PRODUCT(nBins) .GT. np .OR. ANY(nBins .LT. 1) .OR.
     >            ANY(length/nBins .LT. minLen)) THEN
                nFail = nFail + 1
                IF (nFail .LE. 10) PRINT '(A,3ES11.3,A,I5,A,3I5)',
     >            'invalid layout for', length, ',', np, ' ranks:',
     >            nBins
              END IF
            END DO
          END DO
        END DO
      END DO

      PRINT '(I0,A,I0,A,I0,A)', NGOLD, ' golden boxes, ', nCases,
     >  ' scanned layouts, ', nFail, ' failures'
      IF (nFail .GT. 0) CALL EXIT(1)
      END
