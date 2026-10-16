! Slot numbers of the per-particle and per-cell field arrays:
! PPICLF_J* (y), PPICLF_R_* (rprop) and PPICLF_P_* (projected fields).
! They are listed in ppiclf_slots.txt; CMake numbers them and writes
! PPICLF_LAYOUT.h into the build directory. To add a slot, add its name
! to ppiclf_slots.txt. The array sizes (ppiclf_lrs, ppiclf_lrp,
! ppiclf_lrp_int, ppiclf_lrp_pro) and the lists of interpolated and
! ghost slots are set at run time by ppiclf_user_Layout
! (ppiclf_user_Layout.f).

! The viscous-unsteady history (PPICLF_VU, PPICLF_LRP3) is sized at run
! time from the input (VISCOUSUNSTEADY, VUSTEPS); see ppiclf_unsteady.f

! Particle, ghost, cell and wall-triangle counts (formerly PPICLF_LPART,
! PPICLF_LPART_GP, PPICLF_LEE, PPICLF_LWALL) are not set here: ppiclF
! sizes those arrays at run time and grows them (ppiclf_data.f)

#include "PPICLF_LAYOUT.h"
