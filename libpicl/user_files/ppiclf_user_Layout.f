!-----------------------------------------------------------------------
!
! Field layout of a run: how many entries each particle carries in y
! and rprop, how many gas fields are interpolated to the particles and
! how many particle fields are projected to the gas, and which rprop
! slots are interpolated or sent with ghost particles.
!
! Called once by ppiclf_solve_InitLayout, after the PICL flags of the
! input file are set (/RFLU_ppiclF/, PPICLF_USER_COMMON.h) and before
! any particle array is allocated. The slot numbers are the PPICLF_J*,
! PPICLF_R_* and PPICLF_P_* names of PPICLF_USER.h.
!
!   ppiclf_lrs     : y, y1, ydot, ydotc
!   ppiclf_lrp     : rprop
!   ppiclf_lrp_pro : projected fields (PPICLF_P_*)
!   ppiclf_solve_AddInterpSlot(jp): rprop slot jp is interpolated from
!       the gas every stage; rocpicl/PICL_TEMP_Runge*.F90 hands it over
!       with ppiclf_solve_InterpFieldUser(jp, field)
!   ppiclf_solve_AddGhostSlot(jp): rprop slot jp is sent with ghost
!       particles, because ppiclf_user_EvalNearestNeighbor reads it
!       from the neighbour (rpropj)
!
!-----------------------------------------------------------------------
      subroutine ppiclf_user_Layout
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
      integer*4 jp
!
      ppiclf_lrs     = PPICLF_JOZ
      ppiclf_lrp     = PPICLF_R_RTerm2Perp
      ppiclf_lrp_pro = PPICLF_P_JAlphaPT33
!
! Gas at the particle: every slot from PPICLF_R_JRHOF to PPICLF_R_JSDOZ
!
      do jp = PPICLF_R_JRHOF, PPICLF_R_JSDOZ
         call ppiclf_solve_AddInterpSlot(jp)
      end do
!
! Read from the neighbour by the collision and binary added-mass
! models
!
      call ppiclf_solve_AddGhostSlot(PPICLF_R_JDP)
      call ppiclf_solve_AddGhostSlot(PPICLF_R_JRHOP)
      call ppiclf_solve_AddGhostSlot(PPICLF_R_JVOLP)
      call ppiclf_solve_AddGhostSlot(PPICLF_R_WDOTX)
      call ppiclf_solve_AddGhostSlot(PPICLF_R_WDOTY)
      call ppiclf_solve_AddGhostSlot(PPICLF_R_WDOTZ)

      return
      end
