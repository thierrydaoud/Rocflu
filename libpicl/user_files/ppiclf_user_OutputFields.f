!-----------------------------------------------------------------------
!
! The fields written to the particle files (par*.vtu), and their names.
! Positions are always written (as the points of the file). Each call
! adds one field: a slot of ppiclf_y (PPICLF_OUT_Y), ppiclf_rprop
! (PPICLF_OUT_RPROP) or ppiclf_iprop (PPICLF_OUT_TAG).
!
! A slot that is not listed here is still computed and used; it is
! only not written. Groups beyond the default are switched on in the
! PICL section of the input file:
!   OUTPARTICLE 1 : temperature, spin, density and volume of the particle
!   OUTGAS      1 : gas density, velocity, sound speed and temperature
!                   interpolated to the particle
!   OUTGRADGAS  1 : gas gradients and derivatives interpolated to the
!                   particle (for checking the interpolation)
!   OUTMODEL    1 : states and diagnostics of the models that are on
!
! Called once, before the first particle file is written.
!
!-----------------------------------------------------------------------
      subroutine ppiclf_user_OutputFields
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Default: particle velocity, volume fraction, diameter and ID
!
      call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JVX, 'vp_x')
      call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JVY, 'vp_y')
      call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JVZ, 'vp_z')
      ! particle volume fraction: the projected particle volume,
      ! interpolated back to the particle
      call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP, PPICLF_R_JPHIP,
     >                              'phi_p')
      call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP, PPICLF_R_JDP,
     >                              'dp')
      ! ID: index on the rank that created the particle, that rank,
      ! and the cycle it was created on
      call ppiclf_io_AddOutputField(PPICLF_OUT_TAG, 1, 'id')
      call ppiclf_io_AddOutputField(PPICLF_OUT_TAG, 2, 'id_rank')
      call ppiclf_io_AddOutputField(PPICLF_OUT_TAG, 3, 'id_cycle')
!
! OUTPARTICLE
!
      if (out_particle_flag .ge. 1) then
         call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JT, 'Tp')
         call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JOX,
     >                                 'omp_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JOY,
     >                                 'omp_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_Y, PPICLF_JOZ,
     >                                 'omp_z')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHOP, 'rhop')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JVOLP, 'volp')
      endif
!
! OUTGAS
!
      if (out_gas_flag .ge. 1) then
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHOF, 'rho_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JUX, 'ug_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JUY, 'ug_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JUZ, 'ug_z')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JCS, 'c_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JT, 'T_g')
      endif
!
! OUTGRADGAS
!
      if (out_gradgas_flag .ge. 1) then
         ! pressure gradient
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPDX, 'dpdx_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPDY, 'dpdy_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPDZ, 'dpdz_g')
         ! substantial derivative of the gas velocity, Du/Dt
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDRX, 'DuDt_g_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDRY, 'DuDt_g_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDRZ, 'DuDt_g_z')
         ! d(rho phi_g)/dt and grad(rho phi_g)
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHSR, 'drhophidt_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JPGCX, 'drhophidx_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JPGCY, 'drhophidy_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JPGCZ, 'drhophidz_g')
         ! vorticity
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JXVOR, 'vort_g_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JYVOR, 'vort_g_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JZVOR, 'vort_g_z')
         ! gas density gradient
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHOGX, 'drhodx_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHOGY, 'drhody_g')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JRHOGZ, 'drhodz_g')
         ! divergence of the viscous stress
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPVDX, 'divtau_g_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPVDY, 'divtau_g_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JDPVDZ, 'divtau_g_z')
         ! substantial derivative of the vorticity
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDOX, 'DvortDt_g_x')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDOY, 'DvortDt_g_y')
         call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_JSDOZ, 'DvortDt_g_z')
      endif
!
! OUTMODEL: only the models that are on
!
      if (out_model_flag .ge. 1) then
         if (qs_fluct_flag .ge. 1) then
            ! fluctuating quasi-steady force
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_FLUCTFX, 'Ffluct_x')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_FLUCTFY, 'Ffluct_y')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_FLUCTFZ, 'Ffluct_z')
         endif
         if (am_flag .eq. 2) then
            ! density-weighted relative acceleration (binary added mass)
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_WDOTX, 'Wdot_x')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_WDOTY, 'Wdot_y')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_WDOTZ, 'Wdot_z')
         endif
         if (pseudoTurb_flag .ge. 1) then
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_XIPAR, 'xi_par')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_XIPERP, 'xi_perp')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_XIT, 'xi_T')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_RTildePar, 'Rtilde_par')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                               PPICLF_R_RTildePerp, 'Rtilde_perp')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_R11, 'R11')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_R22, 'R22')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_R33, 'R33')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_KMehra, 'k_Mehra')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                    PPICLF_R_KOsnes, 'k_Osnes')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_RparMehra, 'Rpar_Mehra')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                               PPICLF_R_RperpMehra, 'Rperp_Mehra')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                                 PPICLF_R_RTerm2Par, 'Rterm2_par')
            call ppiclf_io_AddOutputField(PPICLF_OUT_RPROP,
     >                               PPICLF_R_RTerm2Perp, 'Rterm2_perp')
         endif
      endif

      return
      end
