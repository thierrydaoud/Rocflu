!-----------------------------------------------------------------------
!
! Sam - Code for solving ydot = F(y,t)
!
! Set internal flags
!   These flags allows the user to turn on or
!      off various force terms, or to select
!      between different versions of each force model
!   To keep particles where they are, set velocity_model = 0
!   To turn off a force set the flag to zero
!   For two-way coupling set collisional_flag = 0
!   For four-way coupling set collisional_flag = 1
!   To turn off feedback force set feedback_flag = 0
!   To use fluctuations turn corresponding flag = 1
!
!
!   velocity_model = 0: particles keep their position and velocity
!         (and temperature); forces are still computed and the
!         feedback force can be on
!                  = 1: particles move under the forces
!
!   qs_flag = 3  ! none = 0; Gidaspow = 1; Parmar = 2; Osnes = 3
!   am_flag = 1  ! Parmar = 1; Briney = 2
!   pg_flag = 1
!   collisional_flag = 1 ! two way coupled = 0 ; four-way = 1
!
!   heattransfer_flag = 1
!
!   ViscousUnsteady_flag = 0 no viscous unsteady drag
!                        = 1 history kernal for visc. unsteady drag
!
!   feedback_flag = 1
!
!   qs_fluct_flag = 1 ! None = 0 ; Lattanzi = 1 ; Osnes = 2 
!
!
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_SetYdot
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 :: ppiclf_rcp_part, ppiclf_p0
      integer :: ppiclf_moveparticle, ierr
      common /RFLU_ppiclf_misc01/ ppiclf_rcp_part
      common /RFLU_ppiclf_misc03/ ppiclf_p0, ppiclf_moveparticle
      real*8 fqsx, fqsy, fqsz
      real*8 fqsforce
      real*8 fqs_fluct(3)
      real*8 xi_par, xi_perp, xi_T, Rsg(3,3), Tsg(3), alpha_PT(3,3)
      real*8 famx, famy, famz 
      real*8 fdpdx, fdpdy, fdpdz
      real*8 fdpvdx, fdpvdy, fdpvdz
      real*8 fcx, fcy, fcz
      real*8 fbx, fby, fbz 
      real*8 fvux, fvuy, fvuz

      real*8 ug, vg, wg

      real*8 beta,cd

      real*8 factor, rmass_add

      real*8 gkern

! Needed for Added mass calculation
      integer*4 j, l
      real*8 SDrho, maxFilter

      real*8 vgradrhog
      integer*4 i, n, ic, k
      integer*4 store_forces

! Needed for heat transfer
      real*8 qq, rmass_therm, temp, qq_du, Nuss

      real*8    Pres

! Needed for angular velocity
      real*8 taux, tauy, tauz, rmass_omega,
     >       taux_hydro, tauy_hydro, tauz_hydro 
      real*8 liftx, lifty, liftz

! Finite Diff Material derivative Variables
      integer*4 nstage, iStage
      integer*4 icallb
      save      icallb
      data      icallb /0/

! Print Data to file
      LOGICAL I_EXIST 
      Character(LEN=25) str 
      integer*4 f_dump
      save      f_dump  
      data      f_dump /1/

      logical exist_file
!
!-----------------------------------------------------------------------
!   

      ! Avery added 10/10/2024 for subbin nearest neighbor search
      
      ! Subbin s holds SBin_list(SBin_start(s):SBin_start(s+1)-1),
      ! ghosts as negative indices; filled by ppiclf_user_subbinMap
      INTEGER*4 SBin_start( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) ))
      INTEGER*4 SBin_list(ppiclf_npart+ppiclf_npart_gp)
      INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin


!
! Code:
!

      icallb = icallb + 1
      nstage = 3
      iStage = mod(icallb,nstage)
      if (iStage .eq. 0) iStage = 3  

      ! Print dt and time every time step (from -v 2 on)
      if (ppiclf_nid==0 .and. ppiclf_verbose>0) then
        if (iStage .eq. 1) then
          write(6,'(a,2x,2(1pe14.6),2x,i3)') '*** PPICLF dt, time = ',
     >      ppiclf_dt,ppiclf_time
        endif
      endif

      rpi         = acos(-1.0d0) ! 3.1415...
      rcp_part    = ppiclf_rcp_part 
      rpr         = 0.70d0 ! Prandtl Number
      rcp_fluid   = 1004.64d0 ! Fluid Specific Heat at const pressure

      fac = ppiclf_rk3dtrk(iStage)*ppiclf_dt

      OneThird = 1.0d0/3.0d0

!
!-----------------------------------------------------------------------
!
! Reset arrays for Viscous Unsteady Force
!
      if (ViscousUnsteady_flag>=1) then
         call ppiclf_user_prop2plag
      endif
!
!-----------------------------------------------------------------------
!
! Avery added 10/10/2024 - Map particles to subbins if collisional force, 
! Briney Added Mass force, QS fluctation force, or pseudo turbulence is flagged
!
      if((am_flag==2).or.(collisional_flag>=1)
     >                  .or.(qs_fluct_flag>=1)
     >                .or.(pseudoTurb_flag .gt. 0)) then
        call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin 
     >                               ,SBin_start ,SBin_list)
        if (collisional_flag>=1 .and. ppiclf_nwall>0
     >      .and. ppiclf_npart>0)
     >     call ppiclf_solve_WallSubbins(i_Bin, n_SBin, tot_SBin)
      endif ! Collisions, QS Fluct, Briney AM, or pseudoTurb flags on

!
!-----------------------------------------------------------------------
!
!
! Evaluate ydot, the rhs of the equations of motion
! for the particles
!

      do i=1,ppiclf_npart

         ! Choose viscosity law
         if (rmu_flag==rmu_fixed_param) then
            ! Constant viscosity law
            rmu = rmu_ref
         elseif (rmu_flag==rmu_suth_param) then
            ! Sutherland law
            temp    = ppiclf_rprop(PPICLF_R_JT,i)
            rmu     = rmu_ref*sqrt(temp/tref)
     >                   *(1.0d0+suth/tref)/(1.0d0+suth/temp)
         else
             call ppiclf_exittr('Unknown viscosity law$', 0.0d0, 0)
         endif
         rkappa = rcp_fluid*rmu/rpr


         ! Useful values
         rmass  = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >              *ppiclf_rprop(PPICLF_R_JRHOP,i)
         vx     = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
         vy     = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
         vz     = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
         vmag   = sqrt(vx*vx + vy*vy + vz*vz)
         rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)
         dp     = ppiclf_rprop(PPICLF_R_JDP,i)
         rep    = vmag*dp*rhof/rmu
         rphip  = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphif  = 1.0d0-ppiclf_rprop(PPICLF_R_JPHIP,i)
         rem    = (1.0d0-rphip)*rep   ! mean slip Reynolds Number
         asndf  = ppiclf_rprop(PPICLF_R_JCS,i)
         rmachp = vmag/asndf
         rhop   = ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Every force is evaluated at any slip, as in MFC. Only the
         ! quasi-steady drag (Step 2) is skipped at zero slip, where its
         ! correlations are singular and the force is zero. Previously
         ! particles with vmag <= 1e-3 skipped all forces, ydot, the
         ! feedback and the neighbor search, so resting particles kept
         ! stale values, never collided and never exchanged heat.

         ! No floor on rep: the models that use it have finite limits
         ! as Re -> 0 (the QS models compute their own Re)

         ! Redefine volume fractions
         ! Need to make sure phi_p + phi_f = 1
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphip = min(rphip,0.62d0)
         rphif = 1.0d0-rphip
         !*** AVERY - we should rethink this limit of 62% pVF ***

         ! TLJ: Needed for viscous unsteady force
         rnu = rmu/rhof

         ! Zero out for each particle i
         famx = 0.0d0; famy = 0.0d0; famz = 0.0d0; rmass_add = 0.0d0;
         Fam(1) = 0.0d0; Fam(2) = 0.0d0; Fam(3) = 0.0d0
         FamUnary(1)=0.0d0;FamUnary(2)=0.0d0;FamUnary(3)=0.0d0;
         FamBinary(1)=0.0d0;FamBinary(2)=0.0d0;FamBinary(3)=0.0d0;
         Wdot_neighbor_mean(1) = 0.0d0; Wdot_neighbor_mean(2) = 0.0d0;
         Wdot_neighbor_mean(3) = 0.0d0; nneighbors = 0.0d0
         fqsx = 0.0d0; fqsy = 0.0d0; fqsz = 0.0d0; beta = 0.0d0;
         fqs_fluct(1)=0.0d0;fqs_fluct(2)=0.0d0;fqs_fluct(3)=0.0d0;
         fdpdx = 0.0d0; fdpdy = 0.0d0; fdpdz = 0.0d0;
         fcx = 0.0d0; fcy = 0.0d0; fcz = 0.0d0;
         taux = 0.0d0; tauy = 0.0d0; tauz = 0.0d0;
         taux_hydro = 0.0d0; tauy_hydro = 0.0d0; tauz_hydro = 0.0d0;
         liftx = 0.0d0; lifty = 0.0d0; liftz = 0.0d0;
         fvux = 0.0d0; fvuy = 0.0d0; fvuz = 0.0d0;
         qq=0.0d0; qq_du=0.0d0; Nuss = 0.0d0
         upmean = 0.0; vpmean = 0.0; wpmean = 0.0;
         u2pmean = 0.0; v2pmean = 0.0; w2pmean = 0.0;
         fdpvdx = 0.0d0; fdpvdy = 0.0d0; fdpvdz = 0.0d0;
         !--- Added for PseudoTurbulence
         Rsg = 0.0d0; Tsg = 0.0d0; alpha_PT = 0.0d0;
         xi_par = 0.0d0; xi_perp = 0.0d0; xi_T=0.0d0

!
! Step 1a: New Added-Mass model of Briney
!
         ! 07/15/2024 - If am_flag = 2, then we need to call
         !   the Unary term before we make any calls to nearest
         !   neighbor
         ! 06/05/2024 - Thierry - For each particle i, initialize
         ! variables to be used in nearest neighbors to zero
         ! before looping over particle j (j neq i)
         ! Briney Added Mass flag
         if (am_flag == 2) then 
            ! 07/14/24 - Thierry - If Briney Algorithm flag and fluct_flag
            !   are ON -> evaluate added-mass unary term before evaluating
            !   neighbor-induced acceleration in EvalNearestNeighbor
            call ppiclf_user_AM_Briney_Unary(i,iStage,
     >           famx,famy,famz,rmass_add)
         endif ! end am_flag = 2

!
! Step 1b: Call NearestNeighbor if particles i and j interact
!
         if ((am_flag==2).or.(collisional_flag>=1)
     >          .or.(qs_fluct_flag>=1)) then

         !AVERY - we should fix vmag ~0 bug and remove conditional check
         if ((qs_fluct_flag>=1) .and. (vmag .gt. 1.d-8)) then
            ! Compute mean for particle i
            !    add neighbor particle j afterward
            ! Box filter is used if qs_fluct_filter_flag=0
            !   The box filter used here is a simple cube centered
            !     at particle i with half-width dist2 (see
            !     ppiclf_user_EvalNearestNeighbor.f for definition)
            !   We use a simple arithmetic mean
            !   phipmean is not used
            ! Gaussian filter is used if qs_fluct_filter_flag=1
            !   We use the value of the Gaussian times the volume
            !     of the particle to get the filtered particle volume

            if (qs_fluct_filter_flag==0) then
               ! box filter
               phipmean = ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = ppiclf_y(PPICLF_JVX,i)
               vpmean   = ppiclf_y(PPICLF_JVY,i)
               wpmean   = ppiclf_y(PPICLF_JVZ,i)
               u2pmean  = upmean**2
               v2pmean  = vpmean**2
               w2pmean  = wpmean**2
               icpmean  = 1
            else if (qs_fluct_filter_flag==1) then
               ! gaussian kernel
               ! r = 0
               !***AVERY - verify max filter dimension should be used***
               maxFilter = MAX(ppiclf_filter(1),ppiclf_filter(2),
     >                         ppiclf_filter(3))
               gkern = sqrt(rpi*maxFilter**2/
     >                (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
               !***AVERY - I think phipmean is using wrong index ***
               ! Also, none of the below are means...
               phipmean = gkern*ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = gkern*ppiclf_y(PPICLF_JVX,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               vpmean   = gkern*ppiclf_y(PPICLF_JVY,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               wpmean   = gkern*ppiclf_y(PPICLF_JVZ,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               u2pmean  = gkern*(ppiclf_y(PPICLF_JVX,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               v2pmean  = gkern*(ppiclf_y(PPICLF_JVY,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               w2pmean  = gkern*(ppiclf_y(PPICLF_JVZ,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               icpmean = 1
            end if
         end if

         ! add neighbors
         ! AVERY - always use subbins
         IF ( 1 .EQ. 1 ) THEN !sbNearest_flag .EQ. 1) THEN
            CALL ppiclf_solve_NearestNeighborSB(
     >           i,tot_SBin,SBin_start,SBin_list,n_SBin)
            IF (collisional_flag>=1 .AND. ppiclf_nwall>0)
     >         CALL ppiclf_solve_WallContacts(i)
         END IF

         end if ! end Step 1b; nearestneighbor


!
! Step 2: Force component quasi-steady
!
         ! Same numbering as MFC (particle_params%qs_force). Before
         ! 2026-10, 1 = Parmar, 2 = Osnes, 3 = Modified Parmar
         ! (removed) and 4 = Gidaspow; 4 now stops with an error.
         ! Zero at zero slip, where the correlations are singular
         if (qs_flag==0 .or. vmag <= 0.0d0) then
           beta = 0.0d0
         else if (qs_flag==1) then
           call ppiclf_user_QS_Gidaspow(i,beta)
         else if (qs_flag==2) then
           call ppiclf_user_QS_Parmar(i,beta)
         else if (qs_flag==3) then
           call ppiclf_user_QS_Osnes (i,beta)
         else
           print*, "***PPICLF: Error in QS Model Selection!",
     >             " QSFORCE must be 0 (none), 1 (Gidaspow),",
     >             " 2 (Parmar) or 3 (Osnes)"
           call ppiclf_exittr('Wrong QS Model Choice$', 0.0, 0)
         endif

         fqsx = beta*vx
         fqsy = beta*vy
         fqsz = beta*vz

!
! Step 3: Force fluctuation for quasi-steady force
!
         ! Note: QS fluctuations needs nearest neighbors,
         !   and is called above in Step 1b
         if (qs_fluct_flag==1) then
            call ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
         elseif (qs_fluct_flag==2) then
            call ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct)
         endif

         ! Add fluctuation part to quasi-steady force
         fqsx = fqsx + fqs_fluct(1)
         fqsy = fqsy + fqs_fluct(2)
         fqsz = fqsz + fqs_fluct(3)

         ! Store quasi-steady fluctuating force
#if PPICLF_R_FLUCTFX > 0
         ppiclf_rprop(PPICLF_R_FLUCTFX,i) = fqs_fluct(1)
         ppiclf_rprop(PPICLF_R_FLUCTFY,i) = fqs_fluct(2)
         ppiclf_rprop(PPICLF_R_FLUCTFZ,i) = fqs_fluct(3)
#endif
         
!
! Step 4: Force component added mass
!
         if (am_flag == 1) then 
            call ppiclf_user_AM_Parmar(i,iStage,
     >           famx,famy,famz,rmass_add)

!-----------------------------------------------------------------------
!Thierry - Added Mass code continues here
         
         elseif (am_flag == 2) then 

            ! Thierry - binary_model.f90 evaluates the terms
            ! in the folllowing order:
            !   (1) Unary Term
            !   (2) Evaluates Neighbor Acceleration
            !   (3) Binary Term
            ! We replicate that here by calling them in the same order
            ! Unary and Binary calculations are now under 
            !   two separate subroutines
            ! Thierry - need to make sure NearestNeighbor is called
            !    if fluct_flag = 0 (ie, no QS fluctuations)

            ! Binary subroutine only valid when number of neighbors .gt. 0
            if (nneighbors .gt. 0) then
               call ppiclf_user_AM_Briney_Binary(i,iStage,
     >              famx,famy,famz,rmass_add)
               FamBinary(1) = famx - FamUnary(1)
               FamBinary(2) = famy - FamUnary(2)
               FamBinary(3) = famz - FamUnary(3)
            else
            ! if particle has no neighbors, need to multiply added mass forces
            ! by volume, as this is taken care of in Binary subroutine
               famx = famx*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famy = famy*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famz = famz*ppiclf_rprop(PPICLF_R_JVOLP,i)
            endif
         endif


!-----------------------------------------------------------------------

!
! Step 5: Force component pressure gradient
!
         if (pg_flag == 1) then
            fdpdx = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDX,i)
            fdpdy = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDY,i)
            fdpdz = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDZ,i)

            if (flow_model == 1) then ! Navier-Stokes Flow Model
               fdpvdx = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDX,i)
               fdpvdy = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDY,i)
               fdpvdz = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDZ,i)
            endif ! flow_model

            fdpdx = fdpdx + fdpvdx
            fdpdy = fdpdy + fdpvdy
            fdpdz = fdpdz + fdpvdz
         endif ! end pg_flag = 1

!
! Step 6: Force component collisional force, ie, particle-particle
!
         if (collisional_flag >= 1) then
            ! Collision force:
            !  A discrete numerical model for granular assemblies
            !  - Cundall and Strack (1979)
            !  - Geotechnique

            ! Sam - STILL NEED TO VALIDATE COLLISION FORCE
            ! Sam - Step 1b already calls nearest neighbor
            
            fcx  = ppiclf_ydotc(PPICLF_JVX,i)
            fcy  = ppiclf_ydotc(PPICLF_JVY,i)
            fcz  = ppiclf_ydotc(PPICLF_JVZ,i) 

         endif ! collisional_flag >= 1


!
! Step 7: Viscous unsteady force with history kernel
!
! Step 7: Viscous-unsteady force with history kernel
!         Diffusive-unsteady heat transfer is also calculated in this call
         ! 1 == original Trapezoidal code for uniform dt
         ! 2 == original Trapezoidal code for variable dt
         ! 3 == modified Trapezoidal method for variable dt
         ! 4 == Hinsberg method for variable dt
         if (ViscousUnsteady_flag==1) then
            call ppiclf_user_history_uniformtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==2) then
            call ppiclf_user_history_variabletrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==3) then
            call ppiclf_user_history_modifiedtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==4) then
            call ppiclf_user_history_hinsberg
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         endif

!
! Step 8a: Heat transfer model
!
         rmass_therm = rmass*rcp_part
         qq = 0.0d0

         if (heattransfer_flag >= 1) then
            call ppiclf_user_HT_driver(i,Nuss,qq)
         endif ! heattransfer_flag >= 1

!
! Step 8b : Diffusive Unsteady Heat Transfer model
!           Calculated in same subroutine as Viscous Unsteady
         if (HTUnsteady_flag==0) qq_du = 0.0d0
!
! Step 9a: Angular velocity model
!
         rmass_omega = rmass*dp*dp/10.0d0

         if (collisional_flag >= 2) then
            taux  = ppiclf_ydotc(PPICLF_JOX,i)
            tauy  = ppiclf_ydotc(PPICLF_JOY,i)
            tauz  = ppiclf_ydotc(PPICLF_JOZ,i) 
            call ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz,
     >                            taux_hydro,tauy_hydro,tauz_hydro)
         endif ! collisional_flag >= 2

!
! Step 9b: Saffman and Magnus Lift models
!          Lift models requires gas-phase vorticity and
!          particle angular velocity
!
         if (collisional_flag == 4) then
            call ppiclf_user_Lift_driver(i,iStage,liftx,lifty,liftz)
         endif ! collisional_flag == 4

!
! Step 10: Pseudo-Turbulence (PT) Models
!          Calculate Reynolds Stress Tensor (Rsg), 
!          and PT-Heat Flux (Qsg)

         ! 03/16/2026 - Thierry - I don't think we need all 3 checks for
         ! the future, but I need that for PseudoTurbulence paper
         ! Ideally, the user should either include all therms, or none. 
         if(pseudoTurb_flag .gt. 0) then
!           ! Rsg = 0, Qsg = 0 
!           call ppiclf_user_PseudoTurb(i,iStage,Nuss,fqsx,fqsy,fqsz,
!     >                                 xi_par,xi_perp,xi_T,
!     >                                 Rsg, T_par, alpha_PT)
!         elseif(pseudoTurb_flag==1) then
!           ! Rsg !=0, Qsg = 0
!           call ppiclf_user_PseudoTurb(i,iStage,Nuss,fqsx,fqsy,fqsz,
!     >                                 xi_par,xi_perp,xi_T,
!     >                                 Rsg, T_par, alpha_PT)

!         elseif(pseudoTurb_flag==2) then
           ! All included
           ! The model needs C_D'/<C_D>, the fluctuation over the mean
           ! drag, so it gets the mean quasi-steady drag beta*v: fqsx
           ! already holds the fluctuation (Step 3), and a total drag
           ! near zero made the ratio blow up (negative R and ksg)
           call ppiclf_user_PseudoTurb(i,iStage,Nuss,
     >                                 beta*vx,beta*vy,beta*vz,
     >                                 xi_par,xi_perp,xi_T,
     >                                 Rsg, Tsg, alpha_PT)
         endif

         ! Store normally distributed random variables xi for PseudoTurbulence
#if PPICLF_R_XIPAR > 0
         ppiclf_rprop(PPICLF_R_XIPAR,i)  = xi_par
         ppiclf_rprop(PPICLF_R_XIPERP,i) = xi_perp
         ppiclf_rprop(PPICLF_R_XIT,i)    = xi_T
#endif
!
! Step 10: Set ydot for all PPICLF_SLN number of equations
!
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+famx+fdpdx+fvux+liftx+fcx)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+famy+fdpdy+fvuy+lifty+fcy)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+famz+fdpdz+fvuz+liftz+fcz)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JT,i)  = (qq+qq_du)/rmass_therm
#if PPICLF_JOX > 0
         ppiclf_ydot(PPICLF_JOX,i) = taux/rmass_omega
         ppiclf_ydot(PPICLF_JOY,i) = tauy/rmass_omega
         ppiclf_ydot(PPICLF_JOZ,i) = tauz/rmass_omega
#endif

!
! Update data for viscous unsteady case at every RK3 stage
!
         if (ViscousUnsteady_flag>=1) then
            call ppiclf_user_UpdatePlag(i,iStage)
         endif

!
! Step 11: Feed Back force to the gas phase
!
         ! Comment: ydotc represented the collisional force in the
         !    particle eqautions above. Here, we over-write the
         !    ydotc vectors for the feedback force used in Rocflu.
         !    Note that Rocflu uses a negative of the RHS, and
         !    so ppiclf must respect this odd convention.
         !
         ! Project work done by hydrodynamic forces:
         !   Inter-phase heat transfer and energy coupling in turbulent 
         !   dispersed multiphase flows
         !   - Ling et al. (2016)
         !   - Physics of Fluids
         ! See also for more details
         !   Explosive dispersal of particles in high speed environments
         !   - Durant et al. (2022)
         !   - Journal of Applied Physics


        IF(feedback_flag==0) THEN
          ppiclf_feedbk(PPICLF_P_JFX,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JFY,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JFZ,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JE,i)  = 0.0d0
        END IF

        IF(feedback_flag==1) THEN
          ! Momentum equations feedback terms
          ppiclf_feedbk(PPICLF_P_JFX,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVX,i)*rmass - fcx)
          ppiclf_feedbk(PPICLF_P_JFY,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVY,i)*rmass - fcy)
          ppiclf_feedbk(PPICLF_P_JFZ,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVZ,i)*rmass - fcz)

          IF(pseudoTurb_flag==0) THEN
          ! Energy equation feedback term
            ppiclf_feedbk(PPICLF_P_JE,i) = ppiclf_rprop(PPICLF_R_JSPL,i)
     >       * ( (fqsx+fvux+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >           (fqsy+fvuy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >           (fqsz+fvuz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
     >                  famx*ppiclf_rprop(PPICLF_R_JUX,i) +
     >                  famy*ppiclf_rprop(PPICLF_R_JUY,i) +
     >                  famz*ppiclf_rprop(PPICLF_R_JUZ,i) +
#if PPICLF_JOX > 0
     >                  taux_hydro*ppiclf_y(PPICLF_JOX,i) +
     >                  tauy_hydro*ppiclf_y(PPICLF_JOY,i) +
     >                  tauz_hydro*ppiclf_y(PPICLF_JOZ,i) +
#endif
     >           qq )
          ELSEIF(pseudoTurb_flag .gt. 0) THEN
            ppiclf_feedbk(PPICLF_P_JE,i) = ppiclf_rprop(PPICLF_R_JSPL,i)
     >       * ( (fqsx+fvux+famx+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >           (fqsy+fvuy+famy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >           (fqsz+fvuz+famz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
#if PPICLF_JOX > 0
     >           taux_hydro*ppiclf_y(PPICLF_JOX,i) +
     >           tauy_hydro*ppiclf_y(PPICLF_JOY,i) +
     >           tauz_hydro*ppiclf_y(PPICLF_JOZ,i) +
#endif
     >           qq )
          END IF ! pseudoTurb_flag

          ! Pseudo-turbulence projected fields: allocated only with
          ! PSEUDOTURB (ppiclf_user_Layout)
          IF(pseudoTurb_flag .GT. 0) THEN
            ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Feedback
            ppiclf_feedbk(PPICLF_P_JRSG11,i) = Rsg(1,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG12,i) = Rsg(1,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG13,i) = Rsg(1,3) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG21,i) = Rsg(2,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG22,i) = Rsg(2,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG23,i) = Rsg(2,3) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG31,i) = Rsg(3,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG32,i) = Rsg(3,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JRSG33,i) = Rsg(3,3)  
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)

            ppiclf_feedbk(PPICLF_P_JTSG1,i) = Tsg(1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JTSG2,i) = Tsg(2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JTSG3,i) = Tsg(3) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)

            ! 03/17/2026 - Thierry - Pseudo Turbulent Heat Flux Projection
            ppiclf_feedbk(PPICLF_P_JAlphaPT11,i) = alpha_PT(1,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT12,i) = alpha_PT(1,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT13,i) = alpha_PT(1,3) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT21,i) = alpha_PT(2,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT22,i) = alpha_PT(2,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT23,i) = alpha_PT(2,3) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT31,i) = alpha_PT(3,1) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT32,i) = alpha_PT(3,2) 
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
            ppiclf_feedbk(PPICLF_P_JAlphaPT33,i) = alpha_PT(3,3)  
     >      *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JSPL,i)
          END IF ! pseudoTurb_flag

        END IF ! Feedback flag

        ! Update volume fraction feedback quantities with feedback on or off
        ppiclf_feedbk(PPICLF_P_JPHIP,i) = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >   *ppiclf_rprop(PPICLF_R_JSPL,i)
        ppiclf_feedbk(PPICLF_P_JPHIPD,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JRHOP,i)
        ppiclf_feedbk(PPICLF_P_JPHIPU,i) = 
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVX,i)
        ppiclf_feedbk(PPICLF_P_JPHIPV,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVY,i)
        ppiclf_feedbk(PPICLF_P_JPHIPW,i) = 
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVZ,i)
        ppiclf_feedbk(PPICLF_P_JPHIPT,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JT,i)


! Step 12: Velocity model 0: particles keep their position, velocity,
! temperature and spin. Feedback can still be on though.
!
         if (velocity_model == 0) then
            ppiclf_ydot(PPICLF_JX ,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JY ,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JZ, i)  = 0.0d0
            ppiclf_ydot(PPICLF_JVX,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JVY,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JVZ,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JT,i)   = 0.0d0
#if PPICLF_JOX > 0
            ppiclf_ydot(PPICLF_JOX,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JOY,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JOZ,i)  = 0.0d0
#endif
         endif




      enddo ! do i=1,ppiclf_npart

!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!06/05/2024 - Thierry - Store density-weighted acceleration
!
      ! Briney Added Mass flag
      if (am_flag==2) then 
         do i=1,ppiclf_npart
     
            ! Substantial derivative of density - how rocflu does it 
            SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)
            
            ! material derivative is phi weighted in Rocflu
            ! drho/dt
            SDrho = SDrho / (rphif)  
            vgradrhog = vx * ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >                  vy * ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >                  vz * ppiclf_rprop(PPICLF_R_JRHOGZ,i)
      
            ! Fluid density
            rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)

            vx = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
            vy = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
            vz = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
            ug = ppiclf_rprop(PPICLF_R_JUX,i)
            vg = ppiclf_rprop(PPICLF_R_JUY,i)
            wg = ppiclf_rprop(PPICLF_R_JUZ,i)
            ! Unary added mass solves rho^g d(u^p)/dt implicitly
            ! Binary added mass solves it explicitly and not implicitly
            ! WDOTX = D(rho^g u^g)/Dt - d(rho^g u^p)/dt)
            ! X-acceleration
            ppiclf_rprop(PPICLF_R_WDOTX,i) =
     >                vx*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i)
     >              + ug*vgradrhog
     >              - rhof*ppiclf_ydot(PPICLF_JVX,i)
          
            ! Y-acceleration
            ppiclf_rprop(PPICLF_R_WDOTY,i) =
     >                vy*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i)
     >              + vg*vgradrhog 
     >              - rhof*ppiclf_ydot(PPICLF_JVY,i)

            ! Z-acceleration
            ppiclf_rprop(PPICLF_R_WDOTZ,i) =
     >                vz*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i)
     >              + wg*vgradrhog
     >              - rhof*ppiclf_ydot(PPICLF_JVZ,i)

         enddo
      endif


!
! ----------------------------------------------------------------------
!
!

      !
      ! Reset arrays for Viscous Unsteady Force
      !
      if (ViscousUnsteady_flag>=1) then
         if (iStage==3) call ppiclf_user_ShiftUnsteadyData
         call ppiclf_user_plag2prop
      endif


! ----------------------------------------------------------------------

      return
      end
