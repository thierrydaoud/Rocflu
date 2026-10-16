#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
! The particle arrays (y, y1, ydot, ydotc, rprop*, feedbk, iprop) are
! allocatable, in module ppiclf_data (ppiclf_data.f)

      COMMON /PPICLF_SLN_CURRENT_N/ PPICLF_NPART
      INTEGER*4 PPICLF_NPART

! The viscous-unsteady history and ppiclf_rprop3 are allocatable, in
! module ppiclf_unsteady (ppiclf_unsteady.f)
