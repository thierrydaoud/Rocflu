#include "PPICLF_USER.h"
#include "PPICLF_STD.h"


!
! General useage
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, flow_model, pseudoTurb_flag,
     >   HTUnsteady_flag, out_particle_flag, out_gas_flag,
     >   out_gradgas_flag, out_model_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ rmu_ref, tref, suth, ksp, erest,
     >   stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, flow_model, pseudoTurb_flag,
     >   HTUnsteady_flag, out_particle_flag, out_gas_flag,
     >   out_gradgas_flag, out_model_flag

      real*8 rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,
     >   rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,
     >   vx,vy,vz,
     >   rcp_part,rpr,
     >   mp, re, rem, rcp_fluid
      common /RFLU_user/ rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,
     >   rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,
     >   vx,vy,vz,
     >   rcp_part,rpr,
     >   mp,re, rem, rcp_fluid

!
! For misc values
!
      real*8 OneThird
      common /ppiclf_misc01/ OneThird
 

!
! For ppiclf_user_Fluctuations.f
!
      integer*4 icpmean
      real*8 upmean, vpmean, wpmean, phipmean
      real*8 u2pmean, v2pmean, w2pmean
      real*8 UnifRnd(6)
      common /user_fluct01/ icpmean
      common /user_fluct02/ upmean, vpmean, wpmean, phipmean
      common /user_fluct03/ u2pmean, v2pmean, w2pmean
      common /user_fluct04/ UnifRnd



!
! For ppiclf_user_AddedMass.f
!
      integer*4 nneighbors
      real*8 Fam(3), FamUnary(3), FamBinary(3),
     >       Wdot_neighbor_mean(3), R_pair(6,6)

      common /user_AddedMass01/ nneighbors
      common /user_AddedMass02/ Fam, FamUnary, FamBinary,
     >                          Wdot_neighbor_mean, R_pair

!
! For ppiclf_solve_InitAngularPeriodic
!
      integer*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag,
     >          ang_case 
      real*8 ang_per_angle, ang_per_xangle,
     >       ang_per_rin, ang_per_rout,
     >       xrot(3) , vrot(3)
      real*8 x_per_min, x_per_max,
     >  y_per_min, y_per_max, z_per_min, z_per_max

      common /solve_InitAngularPeriodic01/ x_per_flag, y_per_flag, 
     >                                z_per_flag, ang_per_flag,
     >                                ang_case
      common /solve_InitAngularPeriodic02/ ang_per_angle, 
     >                                ang_per_xangle,
     >                                ang_per_rin, 
     >                                ang_per_rout,
     >                                xrot, vrot
      common /solve_InitAngularPeriodic03/ x_per_min, x_per_max,
     >  y_per_min, y_per_max, z_per_min, z_per_max
     
