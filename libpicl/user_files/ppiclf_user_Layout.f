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
!       with ppiclf_solve_InterpFieldUser(jp, field). A slot that is not
!       listed is not interpolated and stays zero: list every slot a
!       model that is on reads.
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
      ppiclf_lrs = PPICLF_JOZ
      ppiclf_lrp = PPICLF_R_RTerm2Perp
!
! Projected fields: volume fraction, force and energy feedback, and the
! fields the probes read (PPICLF_P_JPHIPD to PPICLF_P_JPHIPT); the
! pseudo-turbulence stresses and heat flux only with PSEUDOTURB
!
      if (pseudoTurb_flag .gt. 0) then
         ppiclf_lrp_pro = PPICLF_P_JAlphaPT33
      else
         ppiclf_lrp_pro = PPICLF_P_JPHIPT
      endif
!
! Gas at the particle, interpolated every stage. Only the fields a model
! that is on reads; Rocflu hands over all of them and ppiclF ignores the
! others (ppiclf_solve_InterpFieldUser).
!
!     read by the drag, heat transfer and every other model
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JRHOF)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JPHIP)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JUX)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JUY)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JUZ)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JCS)
      call ppiclf_solve_AddInterpSlot(PPICLF_R_JT)
!     pressure gradient: pressure-gradient force, Parmar added mass
      if (pg_flag .eq. 1 .or. am_flag .eq. 1) then
         do jp = PPICLF_R_JDPDX, PPICLF_R_JDPDZ
            call ppiclf_solve_AddInterpSlot(jp)
         end do
      endif
!     Du/Dt, d(rho phi_g)/dt, grad(rho phi_g), grad(rho): added mass
!     and viscous unsteady force
      if (am_flag .ge. 1 .or. ViscousUnsteady_flag .ge. 1) then
         do jp = PPICLF_R_JSDRX, PPICLF_R_JPGCZ
            call ppiclf_solve_AddInterpSlot(jp)
         end do
         do jp = PPICLF_R_JRHOGX, PPICLF_R_JRHOGZ
            call ppiclf_solve_AddInterpSlot(jp)
         end do
      endif
!     vorticity and its substantial derivative: torque and lift
!     (COLLISIONS 2 and above)
      if (collisional_flag .ge. 2) then
         do jp = PPICLF_R_JXVOR, PPICLF_R_JZVOR
            call ppiclf_solve_AddInterpSlot(jp)
         end do
         do jp = PPICLF_R_JSDOX, PPICLF_R_JSDOZ
            call ppiclf_solve_AddInterpSlot(jp)
         end do
      endif
!     divergence of the viscous stress: pressure-gradient force with
!     the Navier-Stokes flow model
      if (pg_flag .eq. 1 .and. flow_model .eq. 1) then
         do jp = PPICLF_R_JDPVDX, PPICLF_R_JDPVDZ
            call ppiclf_solve_AddInterpSlot(jp)
         end do
      endif
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
