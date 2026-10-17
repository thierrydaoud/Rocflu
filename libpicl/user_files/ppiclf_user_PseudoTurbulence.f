!******************************************************************************
! Date Created: March 16, 2026
!      
! Purpose: Calculate subgrid terms for gas momentum and energy
!          equations.
!
! Description: This subroutine calculates:
!              (1) Sub-grid Reynolds Stress Tensor
!              (2) Triple Velocity Correlation
!              (3) Pseudo-Turbulent Heat Flux
!              These terms are then projected back to Rocflu's Momentum
!              and Energu Equations.
!
! Flags: pseudoTurb_flag==0
!        pseudoTurb_flag==1
!        pseudoTurb_flag==2
!
! Input:
!      
! Output: Rsg(3x3), T_par(1x3), alpha_PT(3x3)
!
! Notes: none.
! 
! References:
!     
!      Osnes et al., Pseudo-turbulence models for compressible ﬂow through random
!      particle suspensions, (2025)
!      
!      Zhou et al, Modeling pseudo-turbulent heat flux in gas-solid heat
!      transfer, (2024)
!      
!      Peng et al., Implementation of pseudo-turbulence closures in an
!      Eulerian–Eulerian two-fluid model for non-isothermal gas–solid
!      flow, (2019)
!
!      Sun et al., Pseudo-turbulent heat flux and average gas–phase conduction
!      during gas–solid heat transfer: flow past random fixed particle
!      assemblies, (2016)
!      
!      Mehrabadi et al., Pseudo-turbulent gas-phase velocity
!      fluctuations in homogeneous gas–solid flow: fixed particle
!      assemblies and freely evolving suspensions, (2015)
!      
!******************************************************************************

      subroutine ppiclf_user_PseudoTurb(i,iStage,Nu,fqsx,fqsy,fqsz,
     >                         xi_par,xi_perp,xi_T,
     >                         Rsg, Tsg, alpha_PT)
      use ppiclf_data
!                                                    
      implicit none
!
      include "PPICLF"
!
! Input:
      integer*4 i, iStage
      real*8 Nu, fqsx, fqsy, fqsz
!
! Output:
      real*8 xi_par, xi_perp, xi_T, Rsg(3,3), alpha_PT(3,3), Tsg(3) 
!
! Internal:
      integer*4 m, j
      real*8 aSDE,bq,chi,denum,dW1,dW2,dW3,tF_inv,theta,Z1,Z2,Z3
      real*8 TwoPi
      real*8 bSDE_CD, bSDE_CL, bSDE_CT, CD_prime, CD_average
      real*8 avec(3), bvec(3), cvec(3)
      real*8 eunit(3)
      real*8 s_par, s_perp, s_T, Rmean_par, Rmean_perp, R_par, R_perp
      real*8 R(3,3), Q(3,3), Qt(3,3), Tmean_par(3), T_par(3)
      real*8 k_tilde, k_Mach, b_tilde, b_Mach, b_par, b_perp,
     >       k_Osnes, b_Osnes, KE_mean
      real*8 alpha_fluid, alpha_par, alpha_num, alpha_denum, 
     >       alpha_perp, alpha(3,3)
      real*8 Rpar_Mehra, Rperp_Mehra

! Modeling Constants
      real*8 C1, C2, C3, C4, C5, C6, C7,
     >       D1, D2, D3, D4, D5, D6, D7, D8,
     >       E1, E2, E3, E4
      real*8 F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12,
     >       G1, G2, G3, G4, G5, G6, G7, G8,
     >       H1, H2, H3, H4, H5, H6, H7, H8,
     >       A1, A2, A3, A4
      real*8 D9, D10, D11
      real*8 fit_func
!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      ! Constants taken from Osnes paper, Table 1
      C1 = -1.2152; D1 = -0.0462; E1 = -0.2906;
      C2 = -7.6314; D2 = -0.1068; E2 =  1.1899; 
      C3 =  0.2889; D3 =  0.6793; E3 =  0.5218; 
      C4 =  0.6143; D4 =  1.1461; E4 =  0.0699;
      C5 =  0.3082; D5 = -2.6886; 
                    D6 = -2.1376;
                    D7 =  0.4873;
                    D8 =  0.2395;
                  ! Dr. Bala's terms (and added terms)
       C6 =  2.238; D9  =  0.716
       C7 = -6.918; D10 = -2.140
                    D11 =  1.600

      ! Constants taken from Osnes paper, Table 2
      ! (F11, G4, and G5 are modified)
      F1  = -0.0022; G1 = -0.2867; H1 =  0.4992
      F2  = -0.0219; G2 =  0.2176; H2 = -1.3528
      F3  =  0.0932; G3 =  0.2826; H3 = -0.1358
      F4  = -0.0135; G4 = -0.494 ; H4 = -0.1463
      F5  =  0.0361; G5 =  0.233 ; H5 =  0.2583
      F6  =  0.0403;             ; H6 = -0.3339
      F7  = -0.0761;             ; H7 = -0.0407
      F8  =  0.0599;             ; H8 = -0.0806
      F9  =  0.0164;
      F10 =  0.0453;
      
      F11 = -0.414;
      F12 =  0.230
     
      ! zero out variables  at first
      Rmean_par = 0.0d0; Rmean_perp = 0.0d0
      R = 0.0d0; Rsg = 0.0d0
      Tmean_par = 0.0d0; T_par = 0.0d0
      alpha = 0.0d0; alpha_PT = 0.0d0; alpha_par = 0.0d0;

      !mp = max(0.0d0, min(0.87d0, rmachp))
      !re = max(0.0d0, min(266.0d0, rep))
      mp = rmachp
      re = rep
      rem = (1.0d0-rphip)*re

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!

      ! Granular temperature
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      chi = (1.0 + 2.50*rphip + 4.51*(rphip**2) + 4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      tF_inv = (24.0*rphip*chi/dp) * sqrt(theta/rpi)

! avec : unit vector in main direction (undefined at zero slip, where
! the slip-driven pseudo-turbulence vanishes)
      CD_average = 0.0d0
      if (vmag > 0.0d0) then
        avec = [vx,vy,vz]/vmag

        ! CD_prime has unit of Force
        CD_prime = ppiclf_rprop(PPICLF_R_FLUCTFX,i)*avec(1) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFY,i)*avec(2) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFZ,i)*avec(3)

        ! CD_average has unit of Force, is zero at early time steps
        CD_average = fqsx*avec(1) +
     >               fqsy*avec(2) +
     >               fqsz*avec(3)
      endif

      ! avoiding singularity (also zero slip)
      if(CD_average .lt. 1.d-8) then
        xi_par  = 0.0d0
        xi_perp = 0.0d0
        xi_T    = 0.0d0
        Rsg     = 0.0d0
        Tsg     = 0.0d0
        alpha_PT = 0.0d0
        return
      endif

! Calculate the three orthogonal unit vectors
! The first vector (avec) is vx/vmag, vy/vmag, and vz/vmag
! The second (bvec) is constructued by taking cross-product with eunit
!   Note: if avec is in dir. of e_x=(1,0,0), use e_y=(0,1,0) to get e_z
!       : if avec is in dir. of e_y=(0,1,0), use e_z=(0,0,1) to get e_x
!       : if avec is in dir. of e_z=(0,0,1), use e_x=(1,0,0) to get e_y
! The third (cvec) is cross product of the first two
! written 3/6/24
! bvec, cvec: two orthogonal vectors to avec
      eunit = [1,0,0]
      if (abs(avec(2))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,1,0]
      elseif (abs(avec(1))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,0,1]
      endif

      bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
      bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
      bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
      denum   = max(1.d-8,sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
      bvec    = bvec / denum

      cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
      cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
      cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
      denum   = max(1.d-8,sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
      cvec    = cvec / denum

      ! Generate  Gaussian Random Values
      call ppiclf_op_UniformRandom(ppiclf_rprop(PPICLF_R_SEED,i),
     >                             UnifRnd,6)

      ! Box-Muller transform for generating two independent standard normal
      ! (Gaussian) random variables
      ! Z1 & Z2 are standard normal random variables
      Z1 = sqrt(-2.0d0*log(UnifRnd(1))) * cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3))) * sin(TwoPi*UnifRnd(4))
      Z3 = sqrt(-2.0d0*log(UnifRnd(5))) * cos(TwoPi*UnifRnd(6))

      ! dW1 & dW2 are scaled stochastic amplitudes       
      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2
      dW3 = sqrt(fac)*Z3

c---  Q = [avec | bvec | cvec], 3x3 matrix
c---  avec : unit vector in main direction
c---  bvec, cvec: two orthogonal vectors to avec

        ! 08/15/2025 - Thierry - still need to finalize how to do the
        ! rotation of the R tensor 
      do m=1,3
        Q(m,1) = avec(m)
        Q(m,2) = bvec(m)
        Q(m,3) = cvec(m)
      enddo
  
      Qt = transpose(Q)

      ! Reynolds Subgrid Stress Tensor - Eulerian Mean Model 
                                                                     
      ! Reynolds number and vol fraction dependent k^tilde and b_par 
      ! Mehrabadi's terms
      k_tilde = 2.0*rphip + 2.5*rphip*((1.0-rphip)**3) * 
     >       exp(-rphip*(rem**0.5))                                 
                                                                    
!        b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
!     >          exp(-3.511*phi/(1.0+1.801*exp(-0.005*rem)))

      ! 08/25/2025 - Thierry - Fitted rphip function to better match
      ! formulation with Osnes's low Mach number data
      fit_func = -10.18530152*rphip**3 + 10.94163073*rphip**2
     >            -7.07374862*rphip +  0.38424203
      b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
     >        exp(fit_func/(1.0+1.801*exp(-0.005*rem)))
                                                                       
        ! Mach number correction provided by Osnes                            
!        k_Mach = phi*(C1 + C2*phi + re**C3) * 
!     >        (tanh(C4/C5) + tanh((mp - C4)/C5))

      ! Mach number correction at Re=100, coeff taken from Osnes
      ! cap vol fraction here at 0.3
      k_Mach = min(rphip,0.3)*(C7*min(rphip,0.3) + C6) *
     >      (tanh(C4/C5) + tanh((mp - C4)/C5))

!        b_Mach = (D1 + (re/300.0)*(D2 + D3*re/300.0) +
!     >         phi*(D4 + D5*(re**2/300.0**2) + D6*phi)) *
!     >         (tanh(-D7/D8) - tanh((mp-D7)/D8))

      ! First term was corrected by Dr. Bala (D9, D10, D11)
      b_Mach = rphip*(D9 + D10*rphip + D11*rphip**2) *
     >       (tanh(-D7/D8) - tanh((mp-D7)/D8))

      ! Corrected k^tilde and b_par components                       
      k_Osnes =  k_tilde*(1.0d0 + k_Mach)
      b_Osnes =  b_par  *(1.0d0 + b_Mach)
      b_perp  = -b_Osnes/2.0d0

c--  Multiply by the mean relative flow kinetic energy to dimentionalize      
      KE_mean = 0.5d0 * vmag**2
      k_Osnes = k_Osnes * KE_mean

      Rpar_Mehra  = 2.0d0*rhof*k_tilde*KE_mean
     >                 *(b_par  + 1.0d0/3.0d0)
      Rperp_Mehra  = 2.0d0*rhof*k_tilde*KE_mean
     >                 *(-b_par/2.0  + 1.0d0/3.0d0)

      ! Mean Eulerian Reynolds Subgrid Stress - Parallel Component   
      Rmean_par  = 2.0d0*k_Osnes*(b_Osnes  + 1.0d0/3.0d0)
                                                                     
      ! Mean Eulerian Reynolds Subgrid Stress - Perpendicular Component
      Rmean_perp = 2.0d0*k_Osnes*(b_perp + 1.0d0/3.0d0)

      ! Mean Triple Velocity Correlation
      Tmean_par = E1 + E2*rphip/(E3 + re/300.0) + E4*mp
        
c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize 
      Tmean_par(1)  = Tmean_par(1) * vx * k_Osnes
      Tmean_par(2)  = Tmean_par(2) * vy * k_Osnes
      Tmean_par(3)  = Tmean_par(3) * vz * k_Osnes

c------ Lagrangian Model
  
        ! 08/15/2025 - Ditch A1 per Dr. Bala, and set A2 to this constant per Osnes
      A1 = 0.0d0
      A2 = 0.064

      ! 05/18/2026 - Thierry
      ! This A3 is the old formulation of Osnes, leads to negative Rsg
      ! A3 = G1 + G2/(min(rphip,0.3) + G3) + G4 * mp

      ! New A3 formulation suggested by Osnes
      A3 = (G1 + G2/(min(rphip,0.3) + G3)) * (1.0 + 0.169
     >   *(tanh(G4/G5) - tanh((mp + G4)/G5)))
      
      ! 09/02/2025 - Cap according to Osnes model range
      A4 = H1 + H2*max(0.0d0, min(0.3d0, rphip))
     >        + H3*max(0.0d0, min(0.87d0, mp)) 
     >        + H4*max(30.0d0, min(266.0d0, re))/300.0d0
  
      !s_par = F8 + F9/(rphip + F10) + F11 * mp

      s_par = F8 + F9/(rphip + F10) * (1.0 + 0.085*
     > (tanh(F11/F12) - tanh((mp + F11)/F12)))

      !s_perp = G5/(phi + G6) + (G7*re)/(300.0*phi) + G8
c---   We ditch Osnes expression for s_perp and assume it as big as s_par
      s_perp = s_par
        
      ! 09/02/2025 - Cap according to Osnes model range
      s_T = H5 + H6*max(0.0d0, min(0.3d0, rphip))
     >         + H7*max(0.0d0, min(0.87d0, mp)) 
     >         + H8*max(30.0d0, min(266.0d0, re))/300.0d0

      tF_inv  = (24.0*rphip*chi/dp) * sqrt(theta/rpi)
      aSDE    = tF_inv
      bSDE_CD = s_par *sqrt(2.0*tF_inv)
      bSDE_CL = s_perp*sqrt(2.0*tF_inv)
      bSDE_CT = s_T *sqrt(2.0*tF_inv)

      ! Langevin Model implemented for xi_par, xi_perp, xi_T
      xi_par = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPAR,i)
     >          + bSDE_CD*dW1
      xi_perp = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPERP,i)
     >          + bSDE_CL*dW2
      xi_T = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIT,i)
     >          + bSDE_CT*dW3


      ! Lagrangian Reynolds Subgrid Stress - Parallel Component
      R_par = 1.0 + A1 + A2*CD_prime/max(CD_average, 1.0d-8) + xi_par
  
      ! Lagrangian Reynolds Subgrid Stress - Perpendicular Component
      R_perp = 1.0 + A3*CD_prime/max(CD_average, 1.0d-8) + xi_perp

c--  Multiply Lagrangian Model by Eulerian Model and gas density
      R_par  = R_par  * Rmean_par * rhof
      R_perp = R_perp * Rmean_perp * rhof

c--- R = |R_par,   0   ,   0   |
c---     | 0   , R_perp,   0   |
c---     | 0       0   , R_perp|
  
c---  R matrix only has diagonal components
      R(1,1) = R_par
      R(2,2) = R_perp
      R(3,3) = R_perp
  
c--- Now Rotate the matrix, Rsg = Q . R . Q^T
  
      Rsg = matmul(Q, matmul(R,Qt))

! Diagnostics: built only with MODEL_WRT (ppiclf_slots.txt, ptdiag)
#if PPICLF_R_RTildePar > 0
      ppiclf_rprop(PPICLF_R_RTildePar, i)  = Rmean_par*rhof
      ppiclf_rprop(PPICLF_R_RTildePerp, i) = Rmean_perp*rhof
      ppiclf_rprop(PPICLF_R_R11,i)         = Rsg(1,1)
      ppiclf_rprop(PPICLF_R_R22,i)         = Rsg(2,2)
      ppiclf_rprop(PPICLF_R_R33,i)         = Rsg(3,3)
      ppiclf_rprop(PPICLF_R_KMehra,i)      = 2.0*rhof*k_tilde*KE_mean
      ppiclf_rprop(PPICLF_R_KOsnes,i)      = 2.0*rhof*k_Osnes
      ppiclf_rprop(PPICLF_R_RparMehra,i)   = Rpar_Mehra
      ppiclf_rprop(PPICLF_R_RperpMehra,i)   = Rperp_Mehra
      ppiclf_rprop(PPICLF_R_RTerm2Par,i) = 
     >                    A2*CD_prime/max(CD_average,1.0d-8) 
      ppiclf_rprop(PPICLF_R_RTerm2Perp,i) = 
     >                    A3*CD_prime/max(CD_average,1.0d-8) 
#endif

c--- Osnes Formulation for Triple Velocity Correlation

      T_par = A4 * CD_prime/max(CD_average, 1.0d-8) + xi_T

c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize      
c--  then add mean
      ! check if I mu;itply by avec(1)*vmag instead of vx
      T_par(1) = T_par(1) * vx * k_Osnes
     >           + Tmean_par(1)

      T_par(2) = T_par(2) * vy * k_Osnes
     >           + Tmean_par(2)

      T_par(3) = T_par(3) * vz * k_Osnes
     >           + Tmean_par(3)

      Tsg = matmul(Q, T_par)

      ! Zhou et al.,  Eq. (31)
      ! Parallel component of Pseudo-Turbulent Diffusivity Tensor
      alpha_num = 2.0d0*rem*(rem+1.4d0)*(rpr**2)*exp(-0.002089d0*rem)*
     > (rphif*(-5.11d0*rphip+10.1d0*rphip**2-10.85d0*rphip**3)
     > +1.0d0-exp(-10.96d0*rphip))

      alpha_denum =  3.0d0*rpi*Nu*(1.17d0*rphip-0.2021d0
     > *rphip**(1.0/2.0) +0.08568*rphip**(1.0/4.0))
     > *(rphif**2)*(1.0d0-1.6d0*rphip*rphif-3.0d0*rphip*(rphif**4)
     > *exp((-rem**0.4)*rphip))

      alpha_par = alpha_num/max(alpha_denum, 1.0d-12)

      ! Thermal Diffusivity
      alpha_fluid = rkappa/(rhof * rcp_fluid) 

      ! Multiply by thermal diffusivity
      alpha_par = alpha_par * alpha_fluid
      
      ! Perpendicular Component of Pseudo-Turbulent Diffusivity Tensor
      alpha_perp = (3.0d0*b_perp + 1.0d0)/(3.0d0*b_Osnes + 1.0d0)
     >             * alpha_par

      alpha(1,1) = alpha_par
      alpha(2,2) = alpha_perp
      alpha(3,3) = alpha_perp

      alpha_PT = matmul(Q, matmul(alpha,Qt))

      return
      end
