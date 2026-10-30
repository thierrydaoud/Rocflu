!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer models for non-burning
!    particles
! May update if we include particle combustion
!
! if heattransfer_flag = 0  ignore heat transfer
!                      = 1  Stokes
!                      = 2  Ranz-Marshall (1952)
!                      = 3  Gunn (1977)
!                      = 4  Fox (1978)
!                      = 5  Sun et al. (2015)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_HT_driver(i,Nuss,qq,Nu_bulk,theta)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 qq, Nuss, Q_conv, Nu_bulk, theta

!
! Code:
!
      Q_conv = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      Nuss = 0.0d0; Nu_bulk = 0.0d0; theta = 0.0d0
      if (heattransfer_flag == 1) then
         call HTModel_Stokes(i,Nuss)
      elseif (heattransfer_flag == 2) then
         call HTModel_RM(i,Nuss)
      elseif (heattransfer_flag == 3) then
         call HTModel_Gunn(i,Nuss)
      elseif (heattransfer_flag == 4) then
         call HTModel_Fox(i,Nuss)
      elseif (heattransfer_flag == 5) then
         call HTModel_Sun(i,Nuss,Nu_bulk,theta)
      else
         call ppiclf_exittr('Unknown heat transfer model$', 0.0d0, 0)
      endif

      qq = qq + Q_conv*Nuss


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Stokes limit
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Stokes(i,Nuss)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss
!
! Code:
!
      Nuss = 2.0d0

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Ranz-Marshall
!    define Nusselt number Nu = Nu(Pr,Re)
!
! (1) Ling et al (2012) 
!     "Interaction of a planar shock wave with a dense particle
!     curtain: Modeling and experiments"
!     Physics of Fluids, Vol. 24, 113301
! (2) Durant et al (2022) 
!     "Explosive dispersal of particles in high speed environments"
!     Journal of Applied Physics, Vol. 132, 184902
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_RM(i,Nuss)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss
!
! Code:
!
      Nuss = 2.0d0+0.6d0*(rep**0.5d0)*(rpr**OneThird)

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Gunn
!    define Nusselt number Nu = Nu(Pr,Re,phi)
!
! (1)  Gunn (1978)
!        "Transfer of heat or mass to particles in
!        fixed and fluidised beds"
!        Int. J. Heat Mass Transfer, Vol. 21, pp. 467-476
! (2)  Houim and Oran (2016)
!        "A multiphase model for compressible granular-gaseous
!        flows: formulation and initial tests"
!        J. Fluid Mech., Vol. 789, pp. 166-220
! (3) Boniou and Fox (2023)
!        "Shock-particle-curtain-interaction study with
!        a hyperbolic two-fluid model: Effect of particle
!        force models"
!        Int. Journal of Mutiphase Flow, Vol. 169, 104591
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Gunn(i,Nuss)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 vg
      real*8 Nuss
!
! Code:
!
      ! define bed voidage = ratio of free volume avaliable
      ! for flow to the total volume of bed; aka, volume
      ! fraction of the gas phase
      ! vg = 1.0 - rphip == rphif
      vg = rphif

      ! define Nusselt number Nu = Nu(Pr,Re,phi)
      Nuss = (7.0d0-10.0d0*vg+5.0d0*vg*vg)
     >           *(1.0d0+0.7d0*(rem**0.2d0)*(rpr**OneThird))
     >     + (1.33d0-2.4d0*vg+1.2d0*vg*vg)*(rem**0.7d0)*(rpr**OneThird)

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Fox
!    define Nusselt number Nu = Nu(Pr,Re,M)
!
! (1)  Fox, Rackett, and Nicholls (1978)
!        "Shock wave ignition of magnesium powders"
!        in Proc. 11th Int. Symp. Shock Tubes and Waves.
!        University of Washington Press, Seattle, WA, pp. 262-268
! (2)  Ling, Haselbacher, and Balachandar (2011)
!        "Importance of unsteady contributions to force and heating
!        for particles in compressible flows. Part I: Modeling and
!        anaylsis for shock-particle interaction"
!        Int. Journal of Multiphase Flow, Vol. 37, pp. 1026-1044
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Fox(i,Nuss)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss, Ma_Re
!
! Code:
!
      ! define Nusselt number Nu = Nu(Pr,Re,M)
      ! M/Re = mu/(rho*c*dp) does not depend on the slip, so use it
      ! directly: rmachp/rep is 0/0 at zero slip
      Ma_Re = rmu/(rhof*asndf*dp)
      Nuss = 2.0d0*exp(-rmachp)/(1.0d0+17.0d0*Ma_Re)
     >     + 0.495d0*(rpr**OneThird)*(rep**0.55d0)*
     >       ((1.0d0+0.5d0*exp(-17.0d0*Ma_Re))/1.5d0)

      return
      end
!-----------------------------------------------------------------------
!
! Created October 16, 2025
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: define Nusselt number Nu = Nu(Pr,Re,M)
!
!    B. Sun, S. Tenneti, S. Subramaniam (2015)
!        "Modeling average gas-solid heat transfer using
!        particle-resolved direct numerical simulation"
!     International Journal of Heat and Mass Transfer
!        
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Sun(i,Nuss,Nu_bulk,theta)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss, Nu_bulk, theta, Nu_held, theta_held, Nu_rm, s, w
!
! Code:
!
      ! Fitted for 0.1 <= eps_s <= 0.5 and 1 <= Re_m <= 100
      call HTModel_Sun_NuTheta(rphip, rem, Nu_bulk, theta)

      ! Per unit sphere surface (pi/4) and average fluid T (1/theta)
      if (rphip >= 0.1d0) then
         Nuss = rpi/4.0d0 * Nu_bulk/theta
      else
         ! Below eps_s = 0.1: Sun at 0.1, blended to Ranz-Marshall
         call HTModel_Sun_NuTheta(0.1d0, rem, Nu_held, theta_held)
         Nu_held = rpi/4.0d0 * Nu_held/theta_held
         call HTModel_RM(i,Nu_rm)
         s = rphip/0.1d0
         w = 1.0d0 - s*s*(3.0d0 - 2.0d0*s)
         Nuss = w*Nu_rm + (1.0d0 - w)*Nu_held
      endif

      return
      end
!-----------------------------------------------------------------------
!
! Sun et al. (2015) eq. 27 (Nu per unit perimeter of plane sections,
! bulk fluid T) and eq. 30 (mean scaled fluid temperature)
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Sun_NuTheta(es, Re_mean, Nu_bulk, theta)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 es, Re_mean, Nu_bulk, theta, ef
!
! Code:
!
      ef = 1.0d0 - es
      Nu_bulk = (-0.46d0 + 1.77d0*ef + 0.69d0*ef**2)*(ef**(-3))
     >        + (1.37d0 - 2.4d0*ef + 1.2d0*ef**2)
     >         *(Re_mean**0.7d0)*(rpr**OneThird)
      theta = 1.0d0 - 1.6d0*es*ef
     >      - 3.0d0*es*ef**4*exp(-(Re_mean**0.4d0)*es)

      return
      end
!-----------------------------------------------------------------------
