!-----------------------------------------------------------------------
!
! History of the viscous-unsteady and diffusive-unsteady heat transfer
! models, sized at run time (ppiclf_data_Allocate) from the input:
! ppiclf_nUnsteadyData past steps when VISCOUSUNSTEADY >= 1, none
! otherwise. ppiclf_rprop3 packs the per-particle history so it moves
! with the particle between ranks (ppiclf_user_plag2prop/prop2plag).
!
! Must come first in ppiclf.f (before ppiclf_data.f): the routines
! below USE it.
!
!-----------------------------------------------------------------------
      module ppiclf_unsteady
      implicit none
      ! number of rprop3 entries per particle (10*ppiclf_nUnsteadyData)
      integer*4 ppiclf_lrp3
      real*8, allocatable :: ppiclf_rprop3(:,:)
      real*8, allocatable :: ppiclf_TimeBH(:)
      real*8, allocatable :: ppiclf_drudtPlag(:,:,:)
      real*8, allocatable :: ppiclf_drudtMixt(:,:,:)
      real*8, allocatable :: ppiclf_dTdtPlag(:,:)
      real*8, allocatable :: ppiclf_dTdtMixt(:,:)
      real*8, allocatable :: ppiclf_TMixt(:,:)
      real*8, allocatable :: ppiclf_TPlag(:,:)
      end module ppiclf_unsteady
