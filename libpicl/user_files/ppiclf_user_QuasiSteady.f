!-----------------------------------------------------------------------
!
! Quasi-steady drag models. Each returns beta, with the quasi-steady
! force F_qs = beta*(u_f - u_p).
!
! The models are called only for vmag > 0 and compute their own slip
! Reynolds and Mach numbers from vmag. Unlike the shared rep, Re is not
! floored in the drag terms: each term tends to its Stokes limit as
! Re -> 0. Re is only bounded where it enters log(Re) in Parmar's Mach
! interpolation, which is fitted for ln(Re) = 6.48 to 12.21. These
! match the fixed models in MFC (m_particles_EL_kernels.fpp).
!
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for quasi-steady force
!
! Quasi-steady force (Re_p and Ma_p corrections):
!   Improved Drag Correlation for Spheres and Application
!   to Shock-Tube Experiments
!   - Parmar et al. (2010)
!   - AIAA Journal
!
! Quasi-steady force (phi corrections):
!   The Added Mass, Basset, and Viscous Drag Coefficients
!   in Nondilute Bubbly Liquids Undergoing Small-Amplitude
!   Oscillatory Motion
!   - Sangani et al. (1991)
!   - Phys. Fluids A
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Parmar(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 beta, rcd1, phi_corr
      real*8 reqs, lre, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1,
     >   rcd_M2, C(3), fM(3), factor
      real*8 ppiclf_user_lagrange3
      external ppiclf_user_lagrange3
!
! Code:
!
      reqs = vmag*dp*rhof/rmu
      lre  = log(dmax1(reqs,0.1d0)) ! bounds ln(Re) Mach interpolation

      rmacr = 0.6d0 ! Critical Mach number
      rcd_mcr = (1.0d0+0.15d0*reqs**(0.684d0)) +
     >   (reqs/24.0d0)*(0.513d0/(1.0d0+483.0d0/reqs**(0.669d0)))
      if (rmachp .le. rmacr) then
         rcd_std = (1.0d0+0.15d0*reqs**(0.687d0)) +
     >      (reqs/24.0d0)*(0.42d0/(1.0d0+42500.0d0/reqs**(1.16d0)))
         rmach_rat = rmachp/rmacr
         rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
      else if (rmachp .le. 1.0d0) then
         rcd_M1 = (1.0d0+0.118d0*reqs**0.813d0) +
     >      (reqs/24.0d0)*0.69d0/(1.0d0+3550.0d0/reqs**0.793d0)
         C(1) =  6.48d0
         C(2) =  9.28d0
         C(3) = 12.21d0
         fM(1) = -1.884d0 +8.422d0*rmachp -13.70d0*rmachp**2
     >           +8.162d0*rmachp**3
         fM(2) = -2.228d0 +10.35d0*rmachp -16.96d0*rmachp**2
     >           +9.840d0*rmachp**3
         fM(3) =  4.362d0 -16.91d0*rmachp +19.84d0*rmachp**2
     >           -6.296d0*rmachp**3
         factor = ppiclf_user_lagrange3(lre,C,fM)
         rcd1 = rcd_mcr + (rcd_M1-rcd_mcr)*factor
      else if (rmachp .lt. 1.75d0) then
         rcd_M1 = (1.0d0+0.118d0*reqs**0.813d0) +
     >      (reqs/24.0d0)*0.69d0/(1.0d0+3550.0d0/reqs**0.793d0)
         rcd_M2 = (1.0d0+0.107d0*reqs**0.867d0) +
     >      (reqs/24.0d0)*0.646d0/(1.0d0+861.0d0/reqs**0.634d0)
         C(1) =  6.48d0
         C(2) =  8.93d0
         C(3) = 12.21d0
         fM(1) = -2.963d0 +4.392d0*rmachp -1.169d0*rmachp**2
     >           -0.027d0*rmachp**3 -0.233d0*exp((1.0d0-rmachp)/0.011d0)
         fM(2) = -6.617d0 +12.11d0*rmachp -6.501d0*rmachp**2
     >           +1.182d0*rmachp**3 -0.174d0*exp((1.0d0-rmachp)/0.010d0)
         fM(3) = -5.866d0 +11.57d0*rmachp -6.665d0*rmachp**2
     >           +1.312d0*rmachp**3 -0.350d0*exp((1.0d0-rmachp)/0.012d0)
         factor = ppiclf_user_lagrange3(lre,C,fM)
         rcd1 = rcd_M1 + (rcd_M2-rcd_M1)*factor
      else
         rcd1 = (1.0d0+0.107d0*reqs**0.867d0) +
     >      (reqs/24.0d0)*0.646d0/(1.0d0+861.0d0/reqs**0.634d0)
      end if ! rmachp


      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5
      phi_corr = (1.0d0+5.94d0*min(rphip,0.5d0))

      beta = rcd1*3.0d0*rpi*rmu*dp*phi_corr

      return
      end
!-----------------------------------------------------------------------
!
! Quadratic Lagrange interpolant through (C(i), fM(i)), evaluated at x.
!
!-----------------------------------------------------------------------
!
      real*8 function ppiclf_user_lagrange3(x,C,fM)
!
      implicit none
!
      real*8 x, C(3), fM(3)
!
      ppiclf_user_lagrange3 =
     >     fM(1)*(x-C(2))*(x-C(3))/((C(1)-C(2))*(C(1)-C(3)))
     >   + fM(2)*(x-C(1))*(x-C(3))/((C(2)-C(1))*(C(2)-C(3)))
     >   + fM(3)*(x-C(1))*(x-C(2))/((C(3)-C(1))*(C(3)-C(2)))

      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
! Modified March 5, 2024
!
! Subroutine for quasi-steady force
!
! QS Force calculated as a function of Re, Ma and phi
!
! Use Osnes etal (2023) correlations
! A.N. Osnes, M. Vartdal, M. Khalloufi,
!    J. Capecelatro, and S. Balachandar.
! Comprehensive quasi-steady force correlations for compressible flow
!    through random particle suspensions.
! International Journal of Multiphase Flow, Vol. 165, 104485, (2023).
! doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
!
! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura.
! Supersonic and hypersonic drag coefficients for a sphere.
! AIAA Journal, Vol. 59(8), pp. 3261-3274, (2021).
! doi: https://doi.org/10.2514/1.J060153.
!
! NOTE: Re<45 Rarified formula of Loth et al has been redefined by Balachandar
! to avoid singularity as Ma -> 0.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Osnes(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 gamma,Knp,fKn,CD1,s,CD2,
     >   cd_loth,CM,GM,HM,b1,b2,b3,cd,beta,phi
      real*8 sgby2, JMt, reqs, mpqs, rcd_SN
!
! Code:
!
      gamma = 1.4d0
      ! Mach floor only for Loth's O(M) rarefied terms (CD2, JMt),
      ! which are 0/0 at M = 0
      mpqs = dmax1(rmachp,1.0d-6)
      phi  = rphip
      reqs = vmag*dp*rhof/rmu
      rcd_SN = 1.0d0 + 0.15d0*reqs**(0.687d0) ! Schiller-Naumann

      ! Loth's correlation
      if (reqs .le. 45.0d0) then
         ! Rarefied-dominated regime
         Knp = sqrt(0.5d0 * rpi * gamma) * rmachp / reqs
         fKn = 1.0d0 / (1.0d0
     >         + Knp*(2.514d0 + 0.8d0*exp(-0.55d0/dmax1(Knp,0.01d0))))
         CD1 = (24.0d0/reqs)*rcd_SN*fKn
         sgby2 = sqrt(0.5d0 * gamma)
         s = mpqs * sgby2
         ! J_M continuous at M = 1
         if (mpqs <= 1.0d0) then
            JMt = 2.26d0*(mpqs**4) + 0.14d0*mpqs
         else
            JMt = 1.6d0*(mpqs**4) + 0.25d0*(mpqs**3)
     >            + 0.11d0*(mpqs**2) + 0.44d0*mpqs
         end if
!
! Reformulated version of Loth et al. to avoid singularity at mp = 0
!
         CD2 = (1.0d0 + 2.0d0*(s**2)) * exp(-s**2) * mpqs
     >          / ((sgby2**3)*sqrt(rpi))
     >          + (4.0d0*(s**4) + 4.0d0*(s**2) - 1.0d0)
     >          * erf(s) / (2.0d0*(sgby2**4))
     >          + (2.0d0*(mpqs**3) / (3.0d0 * sgby2)) * sqrt(rpi)

         CD2 = CD2 / (1.0d0 + (((CD2/JMt) - 1.0d0)
     >              * sqrt(reqs/45.0d0)))
         cd_loth = (CD1 + CD2) / (1.0d0 + (mpqs**4))
      else
         ! Compression-dominated regime
         ! TLJ: coefficients tweaked to get continuous values
         !      on the two branches at the critical points
         if (mpqs < 1.5d0) then
            CM = 1.65d0 + 0.65d0*tanh(4d0*mpqs - 3.4d0)
         else
            CM = 2.18d0
     >         - 0.12913149918318745d0*tanh(0.9d0*mpqs - 2.7d0)
         end if
         if (mpqs < 0.8d0) then
            GM = 166.0d0*(mpqs**3) + 3.29d0*(mpqs**2)
     >           - 10.9d0*mpqs + 20.d0
         else
            GM = 5.0d0 + 47.809331200000017d0*(mpqs**(-3))
         end if
         if (mpqs < 1.0d0) then
            HM = 0.0239d0*(mpqs**3) + 0.212d0*(mpqs**2)
     >           - 0.074d0*mpqs + 1.d0
         else
            HM = 0.93967777777777772d0 + 1.0d0 / (3.5d0 + (mpqs**5))
         end if

         cd_loth = (24.0d0/reqs)*rcd_SN*HM +
     >      0.42d0*CM/(1.0d0+42500.0d0/reqs**(1.16d0*CM)
     >                 + GM/sqrt(reqs))

      end if

      b1 = 5.81d0*phi/((1.0d0-phi)**2) +
     >     0.48d0*(phi**(1.d0/3.d0))/((1.0d0-phi)**3)

      b2 = ((1.0d0-phi)**2)*(phi**3)*
     >     reqs*(0.95d0+0.61d0*(phi**3)/((1.0d0-phi)**2))

      b3 = dmin1(sqrt(20.0d0*rmachp),1.0d0)*
     >     (5.65d0*phi-22.0d0*(phi**2)+23.4d0*(phi**3))*
     >     (1.0d0+tanh((rmachp-(0.65d0-0.24d0*phi))/0.35d0))

      cd = cd_loth/(1.0d0-phi) + b3
     >   + (24.0d0/reqs)*(1.0d0-phi)*(b1+b2)

      ! Every 24/Re term in cd is multiplied by Re/24 here (Stokes limit)
      beta = 3.0d0*rpi*rmu*dp*(reqs/24.0d0)*cd

      return
      end
!-----------------------------------------------------------------------
!
! Created Sep. 29, 2025
!
! Subroutine for Quasi-Steady Drag Model of Gidaspow
!
! D. Gidaspow, Multiphase Flow and Fluidization (Academic Press, 1994)
!
! Note: Model is provided per cell volume. We convert that to per
! particle using the particle volume fraction and volume. Every term of
! the per-volume coefficient is proportional to phi, so the division by
! phi is done in closed form and a dilute particle (phi -> 0) keeps its
! single-particle drag.
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Gidaspow(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal variables
!
      integer*4 i
      real*8 cd, beta, phifRep, phif, phi, reqs
!
! Code:
!
      phi  = rphip
      phif = dmax1(1.0d0-rphip,0.0001d0)
      reqs = vmag*dp*rhof/rmu

      ! No Re floor: cd*vmag keeps its Stokes limit
      phifRep = phif*reqs

      if(phifRep .lt. 1000.0d0) then
        cd = 24.0d0/phifRep *(1.0d0+0.15d0*(phifRep)**0.687d0)
      else
        cd = 0.44d0
      endif

      ! Exchange coefficient per unit volume, divided by phi
      if(phif .lt. 0.8d0) then
        beta = 150.0d0*phi*rmu/(phif * dp**2)
     >          + 1.75d0*rhof*vmag/dp
      else
        beta = 0.75d0*cd*rhof*vmag/(dp*phif**1.65d0)
      endif

      beta = beta*(rpi*dp**3)/6.0d0

      return
      end
!-----------------------------------------------------------------------
