!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
! Allocates the viscous-unsteady history (module ppiclf_unsteady) with
! ppiclf_nUnsteadyData steps, set by Rocflu from the input: 0 unless
! VISCOUSUNSTEADY >= 1. Called again on restart; only zeroes then.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_InitZero
      use ppiclf_data
!
      use ppiclf_unsteady
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 n, ierr
!
! Code:
!
      if (.not. allocated(ppiclf_rprop3)) then
         n = max(ppiclf_nUnsteadyData, 0)
         ! history packed per particle: drudt Mixt/Plag (3 each),
         ! dTdt Mixt/Plag, T Mixt/Plag
         ppiclf_lrp3 = 10*n
         allocate(ppiclf_rprop3(ppiclf_lrp3,ppiclf_maxpart),
     >            ppiclf_TimeBH(n),
     >            ppiclf_drudtPlag(3,n,ppiclf_maxpart),
     >            ppiclf_drudtMixt(3,n,ppiclf_maxpart),
     >            ppiclf_dTdtPlag(n,ppiclf_maxpart),
     >            ppiclf_dTdtMixt(n,ppiclf_maxpart),
     >            ppiclf_TMixt(n,ppiclf_maxpart),
     >            ppiclf_TPlag(n,ppiclf_maxpart), stat=ierr)
         if (ierr .ne. 0) call ppiclf_exittr(
     >      'InitZero: cannot allocate the unsteady history$',0.0d0,n)
      endif

      ppiclf_rprop3    = 0.0d0
      ppiclf_TimeBH    = 0.0d0
      ppiclf_drudtMixt = 0.0d0
      ppiclf_drudtPlag = 0.0d0
      ppiclf_dTdtMixt  = 0.0d0
      ppiclf_dTdtPlag  = 0.0d0
      ppiclf_TMixt     = 0.0d0
      ppiclf_TPlag     = 0.0d0

      return
      end
