      SUBROUTINE ppiclf_solve_Initialize(xi1,xpmin,xpmax,
     >           yi1,ypmin,ypmax,zi1,zpmin,zpmax,
     >           ai1,apa,apxa,aprin,aprout)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 xi1, yi1, zi1, ai1
      REAL*8    xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,
     >          apa,apxa,aprin,aprout, pi, angled

      ! Called by rocpicl/PICL_TEMP_InitSolver.F90
      ! xdrange adjusts the bin boundaries to ensure they aren't 
      ! larger than the cartesian fluid domain extremes.

      ! Linear X-Periodicity
      ppiclf_xdrange(1,1) = xpmin
      ppiclf_xdrange(2,1) = xpmax
      x_per_flag = xi1
      IF(x_per_flag.EQ.1) THEN
        IF(xpmin .ge. xpmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'xpMin > xpMax'
          CALL ppiclf_exittr('PeriodicX must have xmin < xmax$',xpmin,0)
        END IF
        ppiclf_linperiodic(1) = .TRUE.
      END IF

      ! Linear Y-Periodicity
      y_per_flag = yi1
      ppiclf_xdrange(1,2) = ypmin
      ppiclf_xdrange(2,2) = ypmax
      IF(y_per_flag.EQ.1) THEN
        IF(ypmin .ge. ypmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'ypMin > ypMax'
          CALL ppiclf_exittr('PeriodicX must have ymin < ymax$',ypmin,0)
        END IF
        ppiclf_linperiodic(2) = .TRUE.
      END IF

      ! Linear Z-Periodicity
      z_per_flag = zi1
      ppiclf_xdrange(1,3) = zpmin
      ppiclf_xdrange(2,3) = zpmax
      IF(z_per_flag.EQ.1) THEN
        IF(zpmin .ge. zpmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'ypMin > ypMax'
          CALL ppiclf_exittr('PeriodicZ must have zmin < zmax$',zpmin,0)
        END IF
        ppiclf_linperiodic(3) = .TRUE.
      END IF

!*** THIS CODE WILL CHANGE
!      ! Angular Periodicity
!      ang_per_flag = ai1
!      IF(ang_per_flag.EQ.1) THEN
!        ppiclf_linperiodic(1) = .TRUE. ! X-Periodicity
!        ppiclf_linperiodic(2) = .TRUE. ! Y-Periodicity
!        ang_per_angle  = apa
!        ang_per_xangle = apxa
!        ang_per_rin    = aprin
!        ang_per_rout   = aprout
!      END IF
!
!      ! User cannot initialize X/Y-Periodicity with Angular Periodicity
!      IF(((x_per_flag.EQ.1).OR.(y_per_flag.EQ.1))
!     >                     .AND.(ang_per_flag.EQ.1))
!     >   CALL ppiclf_exittr('PPICLF: Invalid Periodicity choice$',0,0)
!
!      ! Thierry - compute ang_case
!
!      pi = ACOS(-1.0)
!      angled = ang_per_angle * 180.0d0 / pi ! store angle value in degrees
!
!      IF(ang_per_flag.EQ.0) THEN
!         ang_case = 0 ! standard geometry
!      ELSE
!         IF(angled .lt. 90.0)        ang_case = 1 ! general wedge
!         IF(NINT(angled) .EQ. 90.0)  ang_case = 2 ! quarter cylinder
!         IF(NINT(angled) .EQ. 180.0) ang_case = 3 ! half cylinder
!      END IF
!
!      IF(ppiclf_nid.EQ.0 .AND. ang_case.NE.0) THEN
!         PRINT*, " "
!         PRINT*, " ======================================="
!         PRINT*, " "
!         PRINT*, "!!! PPICLF Angular Periodicity Initialized !!!!"
!         PRINT*, "  Angular periodicity flag =", ang_per_flag
!         PRINT*, "  Init Angular- angle =", ang_per_angle
!         PRINT*, "  Init Angular- angled =", angled
!         PRINT*, "  Init Angular- nint(angled) =", NINT(angled)
!         PRINT*, "  Init Angular- ang_case =", ang_case
!         PRINT*, " "
!         PRINT*, " ======================================="
!         PRINT*, " "
!      END IF
!! *** END CHANGE ***
!
      RETURN
      END
!
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitLayout(nstate, nprop)
      use ppiclf_data
!
! Called from rocpicl/PICL_TEMP_InitSolver.F90 once, after the PICL
! flags are set (/RFLU_ppiclF/) and before any particle is created.
! ppiclf_user_Layout sets the field sizes (module ppiclf_data) and
! lists the interpolated and ghost slots; this checks them and returns
! the sizes of y (nstate) and rprop (nprop), which the caller needs to
! hand particles to ppiclf_solve_InitParticle.
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Output:
!
      INTEGER*4 nstate, nprop
!
! Internal:
!
      INTEGER*4 k
!
      IF(ppiclf_lrs .GT. 0)
     >   CALL ppiclf_exittr('InitLayout called twice$',0.0d0,0)

      CALL ppiclf_user_Layout

      IF(ppiclf_lrs .LT. PPICLF_JVZ)
     >   CALL ppiclf_exittr('y must hold position and velocity$'
     >                     ,0.0d0,ppiclf_lrs)
      IF(ppiclf_lrp .LE. 0)
     >   CALL ppiclf_exittr('rprop has no slots$',0.0d0,ppiclf_lrp)
      IF(ppiclf_lrp_pro .LT. 0)
     >   CALL ppiclf_exittr('Invalid number of projected fields$'
     >                     ,0.0d0,ppiclf_lrp_pro)
      ! a layout without interpolated or ghost slots still gets the
      ! (empty) lists
      IF(.NOT. ALLOCATED(ppiclf_int_slots))
     >   ALLOCATE(ppiclf_int_slots(ppiclf_lrp))
      IF(.NOT. ALLOCATED(ppiclf_ghost_slots))
     >   ALLOCATE(ppiclf_ghost_slots(ppiclf_lrp))
      ALLOCATE(ppiclf_int_given(max(ppiclf_lrp_int,1)))
      ppiclf_int_given  = .FALSE.
      ppiclf_int_ngiven = 0
      ppiclf_lrp_gp = ppiclf_lrs + ppiclf_lrp_ghost

      nstate = ppiclf_lrs
      nprop  = ppiclf_lrp

      IF(ppiclf_nid .EQ. 0) THEN
        WRITE(6,'(A)') ' ppiclF field sizes (ppiclf_user_Layout):'
        WRITE(6,'(A,I4)') '   y     (state)        : ', ppiclf_lrs
        WRITE(6,'(A,I4)') '   rprop (properties)   : ', ppiclf_lrp
        WRITE(6,'(A,I4)') '   interpolated fields  : ', ppiclf_lrp_int
        WRITE(6,'(A,I4)') '   projected fields     : ', ppiclf_lrp_pro
        WRITE(6,'(A,I4)') '   rprop slots of ghosts: ',
     >                    ppiclf_lrp_ghost
        WRITE(6,'(A,64I3)') '   interpolated slots   :',
     >                    (ppiclf_int_slots(k), k=1,ppiclf_lrp_int)
        WRITE(6,'(A,64I3)') '   ghost slots          :',
     >                    (ppiclf_ghost_slots(k), k=1,ppiclf_lrp_ghost)
      END IF

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_AddInterpSlot(jp)
      use ppiclf_data
!
! Called from ppiclf_user_Layout, after ppiclf_lrp is set: rprop slot
! jp is interpolated from the gas every stage. The fluid solver hands
! the field over with ppiclf_solve_InterpFieldUser(jp, field).
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
      INTEGER*4 jp, k
!
      IF(jp .LE. 0 .OR. jp .GT. ppiclf_lrp)
     >   CALL ppiclf_exittr('Interpolated slot outside rprop$',0.0d0,jp)
      IF(.NOT. ALLOCATED(ppiclf_int_slots))
     >   ALLOCATE(ppiclf_int_slots(ppiclf_lrp))
      DO k=1,ppiclf_lrp_int
        IF(ppiclf_int_slots(k) .EQ. jp)
     >    CALL ppiclf_exittr('Interpolated slot listed twice$',0.0d0,jp)
      END DO
      ppiclf_lrp_int = ppiclf_lrp_int + 1
      ppiclf_int_slots(ppiclf_lrp_int) = jp

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_AddGhostSlot(jp)
      use ppiclf_data
!
! Called from ppiclf_user_Layout, after ppiclf_lrp is set: rprop slot
! jp is sent with ghost particles, because the neighbour routine
! (ppiclf_user_EvalNearestNeighbor) reads it from the neighbour
! (rpropj). A ghost's other rprop slots are zero in that routine.
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
      INTEGER*4 jp, k
!
      IF(jp .LE. 0 .OR. jp .GT. ppiclf_lrp)
     >   CALL ppiclf_exittr('Ghost slot outside rprop$',0.0d0,jp)
      IF(.NOT. ALLOCATED(ppiclf_ghost_slots))
     >   ALLOCATE(ppiclf_ghost_slots(ppiclf_lrp))
      DO k=1,ppiclf_lrp_ghost
        IF(ppiclf_ghost_slots(k) .EQ. jp)
     >    CALL ppiclf_exittr('Ghost slot listed twice$',0.0d0,jp)
      END DO
      ppiclf_lrp_ghost = ppiclf_lrp_ghost + 1
      ppiclf_ghost_slots(ppiclf_lrp_ghost) = jp

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
     > bind(C, name="ppiclc_solve_InitParticle")
#else
      SUBROUTINE ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
#endif
      use ppiclf_data
!
! Called from rocpicl/PICL_TEMP_InitSolver.F90

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
! Input: 
!
      INTEGER*4  imethod ! From rocpicl: 2 (Same RK3 as Rocflu)
      INTEGER*4  ndim    ! From rocpicl: 3
      INTEGER*4  iendian ! From rocpicl: 0
      INTEGER*4  npart
      INTEGER*4  l
      REAL*8     y(*)
      REAL*8     rprop(*)
      REAL*8     filt2(3)
      REAL*8     filt3

      IF(.NOT. PPICLF_LCOMM)
     >CALL ppiclf_exittr('InitMPI must be before InitParticle$',0.0d0
     >   ,ppiclf_nid)
      IF(PPICLF_OVERLAP)
     >CALL ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)

      CALL ppiclf_prints('*Begin InitParticle$')
      CALL ppiclf_prints('   *Begin InitParam$')

      CALL ppiclf_solve_InitParam(imethod,ndim,iendian)
      DO l = 1,3
        ppiclf_filter(l) = filt2(l)
      END DO
      ppiclf_nndist = filt3
      
      CALL ppiclf_prints('    End InitParam$')

      IF(.NOT. PPICLF_RESTART) THEN
        CALL ppiclf_prints('   *Begin InitZero$')
        CALL ppiclf_solve_InitZero
        CALL ppiclf_prints('   *End InitZero$')
        CALL ppiclf_prints('   *Begin AddParticles$')
        CALL ppiclf_solve_AddParticles(npart,y,rprop)
        CALL ppiclf_prints('   *End AddParticles$')
      END IF

      CALL ppiclf_prints(' End InitParticle$')
!
      PPICLF_LINIT = .TRUE.

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitParam(imethod,ndim,iendian)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4  imethod
      INTEGER*4  ndim
      INTEGER*4  iendian
!
      IF(imethod .EQ. 0 .OR. imethod .GE. 3 .OR. imethod .LE. -2)
     >   CALL ppiclf_exittr('Invalid integration method$',0.0d0,imethod)
      IF(ndim .NE. 3) THEN
        PRINT*, 'ERROR: Rocpicl is written for 3D grids.'
        CALL ppiclf_exittr('Invalid problem dimension$',0.0d0,ndim)
      END IF
      IF(iendian .LT. 0 .OR. iendian .GT. 1)
     >   CALL ppiclf_exittr('Invalid Endian$',0.0d0,iendian)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian

      ppiclf_linperiodic(1) = .FALSE.    
      ppiclf_linperiodic(2) = .FALSE.   
      ppiclf_linperiodic(3) = .FALSE. 

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0d0
      ppiclf_time   = 0.0d0

!      ppiclf_readytosolve = .FALSE.
      ppiclf_overlap      = .FALSE.
      ppiclf_linit        = .FALSE.
      ppiclf_lintp        = .FALSE.
      ppiclf_lproj        = .FALSE.
      ppiclf_binchanged   = .TRUE.
      ppiclf_printbinvtu  = .TRUE.
      IF(ppiclf_lrp_int .GT. 0) ppiclf_lintp = .TRUE.
      IF(ppiclf_lrp_pro .GT. 0) ppiclf_lproj = .TRUE.

      ppiclf_xdrange(1,1) = -1E20
      ppiclf_xdrange(2,1) =  1E20
      ppiclf_xdrange(1,2) = -1E20
      ppiclf_xdrange(2,2) =  1E20
      ppiclf_xdrange(1,3) = -1E20
      ppiclf_xdrange(2,3) =  1E20

      ppiclf_filter(1)        = 0.0D0
      ppiclf_filter(2)        = 0.0D0
      ppiclf_filter(3)        = 0.0D0
      ppiclf_nndist           = 0.0D0
      ppiclf_interp_dchk(1) = 0.0D0
      ppiclf_interp_dchk(2) = 0.0D0
      ppiclf_interp_dchk(3) = 0.0D0
 
      ppiclf_n_bins(1) = 1
      ppiclf_n_bins(2) = 1
      ppiclf_n_bins(3) = 1

      ppiclf_binb(1) = 0.0
      ppiclf_binb(2) = 0.0
      ppiclf_binb(3) = 0.0
      ppiclf_binb(4) = 0.0
      ppiclf_binb(5) = 0.0
      ppiclf_binb(6) = 0.0

      ppiclf_previousbinb(1) =  1.0E9
      ppiclf_previousbinb(2) = -1.0E9
      ppiclf_previousbinb(3) =  1.0E9
      ppiclf_previousbinb(4) = -1.0E9
      ppiclf_previousbinb(5) =  1.0E9
      ppiclf_previousbinb(6) = -1.0E9

      ppiclf_bins_dx(1) = 1.0
      ppiclf_bins_dx(2) = 1.0
      ppiclf_bins_dx(3) = 1.0     

      ppiclf_nwall    = 0
      ppiclf_iwallm   = 0

      ppiclf_int_ngiven = 0
      IF(ALLOCATED(ppiclf_int_given)) ppiclf_int_given = .FALSE.

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitZero
      use ppiclf_data
      use ppiclf_unsteady
!
! Allocates ppiclF's arrays on the first call and zeroes them (called
! again on restart, after the particle file is read)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
      ! Allocate the arrays of module ppiclf_data on the first call.
      ! They start with room for one entry and grow as the run needs
      ! (ppiclf_data_Grow*).
      IF (ppiclf_maxpart .EQ. 0) THEN
        ppiclf_maxpart  = 1
        ppiclf_maxghost = 1
        ppiclf_maxcell  = 1
        ppiclf_maxfv    = 1
        ppiclf_maxwall  = 1
        CALL ppiclf_data_Allocate
      END IF

      ! zero'ing real and ghost particle properties
      ppiclf_y     = 0.0D0
      ppiclf_ydot  = 0.0D0
      ppiclf_ydotc = 0.0D0
      ppiclf_y1    = 0.0D0
      ppiclf_rprop  = 0.0D0
      ppiclf_iprop  = 0
      ppiclf_feedbk = 0.0D0
      ppiclf_iprop_gp = 0
      ppiclf_rprop_gp = 0.0D0

      ppiclf_npart = 0
      ! zero'ing grid properties for interpolation
      ppiclf_int_fld = 0.0D0

      ! viscous-unsteady history
      ppiclf_rprop3    = 0.0D0
      ppiclf_TimeBH    = 0.0D0
      ppiclf_drudtMixt = 0.0D0
      ppiclf_drudtPlag = 0.0D0
      ppiclf_dTdtMixt  = 0.0D0
      ppiclf_dTdtPlag  = 0.0D0
      ppiclf_TMixt     = 0.0D0
      ppiclf_TPlag     = 0.0D0

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
#else
      SUBROUTINE ppiclf_solve_AddParticles(npart,y,rprop)
#endif
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4  npart
      REAL*8     y(*)
      REAL*8     rprop(*)
!
! Internal:
!
      INTEGER*4 ppiclf_iglsum,ntotal,i,j
      external ppiclf_iglsum
!

      CALL ppiclf_prints('   *Begin AddParticles$')

      IF(npart .lt. 0)
     >   CALL ppiclf_exittr('Invalid number of particles$',
     >                      0.0D0,npart)
      CALL ppiclf_data_GrowPart(ppiclf_npart + npart)

      CALL ppiclf_printsi('      -Begin copy particles$',npart)

      CALL ppiclf_copy(ppiclf_y(1,ppiclf_npart+1),
     >                 y,
     >                 npart*ppiclf_lrs)
      CALL ppiclf_copy(ppiclf_rprop(1,ppiclf_npart+1),
     >                 rprop, npart*ppiclf_lrp)

      ppiclf_npart = ppiclf_npart + npart

      CALL ppiclf_printsi('      -Begin copy particles$',ppiclf_npart)

      IF(.NOT. PPICLF_RESTART) THEN
         CALL ppiclf_prints('      -Begin ParticleTag$')
            CALL ppiclf_solve_SetParticleTag(npart)
         CALL ppiclf_prints('       End ParticleTag$')
      END IF

      CALL ppiclf_prints('    End AddParticles$')

      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_SetParticleTag(npart)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      INTEGER*4 npart
! 
! Internal: 
! 
      INTEGER*4 i
!
      DO i=ppiclf_npart-npart+1,ppiclf_npart
         ppiclf_iprop(1,i) = i
         ppiclf_iprop(2,i) = ppiclf_nid
         ppiclf_iprop(3,i) = ppiclf_cycle
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_NearestNeighborSB(i,SBt,SBs,SBl,SBn)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input:
!
      INTEGER*4 i, SBt, SBn(3)
      ! Particles of subbin s: SBl(SBs(s):SBs(s+1)-1); ghosts negative
      INTEGER*4 SBs(0:SBt), SBl(*)
! 
! Internal: 
! 
      REAL*8 xp(3), bin_xMin(3), distSQ, xdistSQ, ydistSQ, zdistSQ,
     >       dist_total
      ! a ghost's rprop, unpacked to full length (slots a ghost does
      ! not carry stay zero)
      REAL*8 rpropj(ppiclf_lrp)
      INTEGER*4 k, j, jp, l, iSB, jSB, kSB, loopSB, tempSB, iSBin(3)
! 
      distSQ = ppiclf_nndist**2
      rpropj = 0.0D0
      
      ! find ith particle subbin (tempSB)
      DO l = 1,3
        xp(l) = ppiclf_y(l,i)
        bin_xMin(l) = ppiclf_bin_pos(1,l)
      END DO
      DO l = 1,3
        iSBin(l) = FLOOR((xp(l) - (bin_xMin(l)
     >         - ppiclf_nndist))/ppiclf_nndist)
      END DO
      tempSB = iSBin(1) + iSBin(2)*SBn(1) + iSBin(3)*SBn(1)*SBn(2)
      ! Loop through real particles
      DO iSB = 1,3     !to look at -1,current,+1 x-dir subbins
        DO jSB = 1,3   !to look at -1,current,+1 x-dir subbins
          DO kSB = 1,3 !to look at -1,current,+1 x-dir subbins
          ! Loops through 27 adjacent subbins
          loopSB = tempSB + (-2+iSB) + (-2+jSB)*SBn(1) 
     >             + (-2+kSB)*SBn(1)*SBn(2)
          IF (loopSB .GT. -1 .AND. loopSB .LT. SBt) THEN
            DO k = SBs(loopSB),SBs(loopSB+1)-1
              j = SBl(k)
              IF (j .GT. 0) THEN ! Real particle
                IF (j .EQ. i) CYCLE ! Same particle
                xdistSQ = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
                IF (xdistSQ .GT. distSQ) CYCLE
                ydistSQ = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
                IF (ydistSQ .GT. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                zdistSQ = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
                IF (zdistSQ .GT. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                IF (dist_total .GT. distSQ) CYCLE
                CALL ppiclf_user_EvalNearestNeighbor(i,j
     >                                   ,ppiclf_cp_map(1,i)
     >                                   ,ppiclf_cp_map(1+ppiclf_lrs,i)
     >                                   ,ppiclf_cp_map(1,j)
     >                                   ,ppiclf_cp_map(1+ppiclf_lrs,j))
              ELSE IF (j .LT. 0) THEN ! Ghost Particle
                ! Negative was just use for ghost particle indicator
                ! in subbin mapping array. Need to flip sign
                j = - j                 
                xdistSQ = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
                IF (xdistSQ .GT. distSQ) CYCLE
                ydistSQ = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
                IF (ydistSQ .GT. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                IF (ppiclf_ndim .EQ. 3) THEN
                zdistSQ = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
                IF (zdistSQ .GT. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                END IF
                IF (dist_total .GT. distSQ) CYCLE
                jp = -1*j
                DO l = 1,ppiclf_lrp_ghost
                  rpropj(ppiclf_ghost_slots(l)) =
     >              ppiclf_rprop_gp(ppiclf_lrs+l,j)
                END DO
                CALL ppiclf_user_EvalNearestNeighbor(i,jp
     >                                 ,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+ppiclf_lrs,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,rpropj)
              END IF
            END DO !k
          END IF ! if loopSB is valid
        END DO !kSB
      END DO !jSB
      END DO !iSB

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_WriteVTU(time)
     > bind(C, name="ppiclc_solve_WriteVTU")
#else
      SUBROUTINE ppiclf_solve_WriteVTU(time)
#endif
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      REAL*8    time
! 
! Internal:
!
      ppiclf_time   = time

      CALL ppiclf_io_WriteParticleVTU('')
      CALL ppiclf_io_WriteBinVTU('')
      ! Output diagnostics
      CALL ppiclf_io_OutputDiagAll

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
#else
      SUBROUTINE ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
#endif
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      INTEGER*4 istep
      INTEGER*4 iostep
      REAL*8    dt
      REAL*8    time
! 
! Internal:
!
      LOGICAL iout
!
      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      IF(ppiclf_imethod .EQ. 2) THEN
        CALL ppiclf_solve_IntegrateRK3s_Rocflu(iout)
      ELSE
        PRINT*, 'ERROR: Wrong RK selected for rocpicl. Use RK3!'
        CALL ppiclf_exittr('Wrong RK for rocpicl',0.D0,0)
      END IF

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_IntegrateRK3s_Rocflu(iout)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 i, ndum, nstage, istage
      INTEGER*4 icalld
      INTEGER*4 j
      save      icalld
      data      icalld /0/

      EXTERNAL  ppiclf_glsum
      REAL*8    ppiclf_glsum
      REAL*8    fxsum, fysum, fzsum, fxabssum, fyabssum, fzabssum

!
! Output:
!
      LOGICAL iout
!
      icalld = icalld + 1

      ! get rk3 coeffs
      CALL ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = MOD(icalld,nstage)
      if (istage .EQ. 0) istage = 3
      iout = .FALSE.
      if (istage .EQ. nstage) iout = .TRUE.

      ! evaluate ydot
      CALL ppiclf_solve_SetYdot

      !Zero out for first stage
      if (istage .EQ. 1) then
        ppiclf_y1 = 0.0D0
      END IF

      ! The Rocflu RK3 can be found in equation (7) of:
      ! S. Yu. "Runge-Kutta Methods Combined with Compact
      !   Difference Schemes for the Unsteady Euler Equations".
      !   Center for Modeling of Turbulence and Transition.
      !   Research Briefs, 1991.

      DO i = 1,PPICLF_NPART
        DO j = 1,ppiclf_lrs
          ppiclf_y(j,i) =  - ppiclf_rk3coef(1,istage)*ppiclf_y1   (j,i)
     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (j,i)
     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (j,i)
        END DO
      END DO
      
      IF(ppiclf_linperiodic(1) .OR. ppiclf_linperiodic(2) .OR.
     >                             ppiclf_linperiodic(3)) THEN
        CALL ppiclf_solve_PeriodicParticleShift
      END IF

      !Store Current stage RHS for next stage's use
      DO i= 1,PPICLF_NPART
        DO j= 1,ppiclf_lrs
          ppiclf_y1(j,i) =  ppiclf_ydot(j,i)
        END DO
      END DO
      

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_PeriodicParticleShift
      use ppiclf_data
!
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 i  
      REAL*8 per_alpha
!
      DO i=1,ppiclf_npart
!!!!!!!!!!!!!!!!        Rotational Periodicity Starts Here     !!!!!!!!!!!!!!!!!!!!
!            ! currently only supports angular rotation around z-axis
!            ! and only check theta component and not radial
!            ! applying the radial periodicity is very straightforward, but not needed for now
!
!         IF(ang_per_flag .EQ. 1) THEN  ! Angular periodicity
!           
!           ! particle angle w/ x-axis
!           ! per_alpha here is obtained in radians
!           ! ang_per_angle & ang_per_xangle are transformed 
!           !   to radians in PICL_TEMP_InitFlowSolver.F90
!           per_alpha = 
!     >         atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))
!
!           ! check if particle leaving through lower face or upper face of wedge
!           IF((per_alpha .LT. ang_per_xangle) .OR. 
!     >         (per_alpha .GT. (ang_per_xangle + ang_per_angle))) THEN
!             CALL ppiclf_solve_InvokeAngularPeriodic(i, ang_per_flag,
!     >                                                per_alpha, 
!     >                                                ang_per_angle,
!     >                                                ang_per_xangle,
!     >                                                1)
!           END IF ! per_alpha
!         END IF ! ang_per_flag
!
         ! Linear Periodicity Invoked
         IF(ppiclf_linperiodic(1) .OR. ppiclf_linperiodic(2) 
     >                            .OR. ppiclf_linperiodic(3)) THEN
           CALL ppiclf_solve_InvokeLinearPeriodic(i)
         END IF 
      END DO ! i=1,ppiclf_part

      RETURN
      END

!----------------------------------------------------------------------- 

      SUBROUTINE ppiclf_solve_InvokeLinearPeriodic(i)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 i, j
!

! Case 1 - Linear Periodicity in any of 3 directions ; NO Anuglar Periodicity
      IF(((x_per_flag.EQ.1).OR.(y_per_flag.EQ.1).OR.(z_per_flag.EQ.1))
     >   .AND.(ang_per_flag.EQ.0)) THEN

        DO j= 1,3
          ! particle leaving min. periodic face -> move it relative to 
          !                                         max periodic face
          IF(ppiclf_y(j,i) .LT. ppiclf_xdrange(1,j)) THEN
             ppiclf_y(j,i) = ppiclf_xdrange(2,j) - 
     >                    ABS(ppiclf_xdrange(1,j) - ppiclf_y(j,i))
            CYCLE
          END IF

          ! particle leaving max. periodic face -> move it relative to 
          !                                         min periodic face
          IF(ppiclf_y(j,i).GT.ppiclf_xdrange(2,j)) THEN
             ppiclf_y(j,i) = ppiclf_xdrange(1,j) + 
     >                    ABS(ppiclf_y(j,i) - ppiclf_xdrange(2,j))
          END IF
        END DO ! j
      END IF 
 
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_SetRK3Coeff(dt)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      REAL*8 dt
!
      IF(ppiclf_imethod .EQ. 2) THEN
        !BD:Rocflu's rk3 scheme
        !Folowing form
        !rk3(1,:) = Temp storage i.e. Previous Stage RHS
        !rk3(2,:) = Temp storage i.e. Current Stage iteration
        !rk3(3,:) = Temp storage i.e. Current Stage RHS

        ! 11/20/2025 - These are the time fraction steps between
        ! the different RK stages. Multiply by ppiclf_dt to get
        ! the dt at the current stage
        ppiclf_rk3dtrk(1) = 8.0d0/15.0d0
        ppiclf_rk3dtrk(2) = 2.0d0/15.0d0
        ppiclf_rk3dtrk(3) = 5.0d0/15.0d0

        ppiclf_rk3ark(1) = 8.0d0/15.0d0
        ppiclf_rk3ark(2) = 5.0d0/12.0d0
        ppiclf_rk3ark(3) = 0.75d0

        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0
        ppiclf_rk3coef(3,1) = dt*8.0d0/15.0d0
        ppiclf_rk3coef(1,2) = dt*17.0d0/60.0d0
        ppiclf_rk3coef(2,2) = 1.0d0
        ppiclf_rk3coef(3,2) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(1,3) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(2,3) = 1.0d0
        ppiclf_rk3coef(3,3) = dt*3.0d0/4.0d0
!      ELSE
!        !BD:Original Code This follows CMT-nek's rk 3 scheme
!        ppiclf_rk3coef(1,1) = 0.d00
!        ppiclf_rk3coef(2,1) = 1.0d0 
!        ppiclf_rk3coef(3,1) = dt
!        ppiclf_rk3coef(1,2) = 3.0d0/4.0d0
!        ppiclf_rk3coef(2,2) = 1.0d0/4.0d0 
!        ppiclf_rk3coef(3,2) = dt/4.0d0
!        ppiclf_rk3coef(1,3) = 1.0d0/3.0d0
!        ppiclf_rk3coef(2,3) = 2.0d0/3.0d0
!        ppiclf_rk3coef(3,3) = dt*2.0d0/3.0d0
!        !BD: Original Code END
      END IF

      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_SetYdot
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
!      IF(.NOT. PPICLF_READYTOSOLVE)
!     >  CALL ppiclf_solve_InitSolve
      CALL ppiclf_solve_InitSolve
      CALL ppiclf_user_SetYdot

      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitSolve
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! 
! Internal: 
! 
      INTEGER*4 :: i, j
      ! ppiclf_binchanged set in CreateBin
      ! ppiclf_binchanged .TRUE. means
      ! bin coordinates changed
      CALL ppiclf_comm_CreateBin

      ! ppiclf_particleMoved set in FindParticle
      ! ppiclf_particleMoved .EQ. 0 means all particles
      ! stayed in same bin as previous RK Stage.
      CALL ppiclf_comm_FindParticle
      IF(ppiclf_particleMoved .NE. 0 .OR.
     >              ppiclf_binchanged) THEN
        CALL ppiclf_comm_MoveParticle
      END IF

      IF(ppiclf_overlap .AND. ppiclf_binchanged) THEN
        CALL ppiclf_comm_MapOverlapGrid
      END IF

      IF(ppiclf_overlap) THEN
        ! Interpolate fluid solver grid to particle
        CALL ppiclf_solve_InterpParticleGrid
        ! Project particle feedback to fluid solver grid
        CALL ppiclf_solve_ProjectParticleGrid
      END IF

!      IF(ppiclf_gprequired) THEN
      ! Ghost particles are needed 
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
!      END IF

      ! Zero collisions 
      ppiclf_ydotc = 0.0D0

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpParticleGrid
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 j
!
      ! Every interpolated field must have been handed over this stage
      ! (ppiclf_solve_InterpFieldUser), or none: the call from
      ! ppiclf_solve_InitSolve comes before the fluid solver has given
      ! any, and then nothing is interpolated
      IF(ppiclf_int_ngiven .NE. 0 .AND.
     >   ppiclf_int_ngiven .NE. ppiclf_lrp_int) THEN
        DO j=1,ppiclf_lrp_int
          IF(.NOT. ppiclf_int_given(j))
     >      CALL ppiclf_exittr('Interpolated field not given, slot$'
     >                        ,0.0d0,ppiclf_int_slots(j))
        END DO
      END IF

      ! Copies Grid Cell ID for all Rocflu elements that map
      ! to ppiclf domain for GSLIB Transfer.  This copy is from
      ! MapOverlapGrid.
      CALL ppiclf_solve_InitInterp

      ! Makes array (ppiclf_int_fld_input) of all rprop data
      ! for grid cellss that map to ppiclf domain.
      IF(ppiclf_int_ngiven .GT. 0) THEN
        DO j=1,ppiclf_lrp_int
           CALL ppiclf_solve_InterpField(j)
        END DO
      END IF
      
      ! Transfers ppiclf_er_mapc & ppiclf_int_fld for all Rocflu Grid
      ! cells that map to ppiclf domain.
      CALL ppiclf_solve_InterpTupleTransfer

      ! Maps up to 27 closest cell centers to particle (projection
      ! stencil) and interpolates rprop data to the particles
      ! (ppiclf_solve_InterpLS)
      CALL ppiclf_solve_SBParticleToCellMap

      ppiclf_int_ngiven = 0
      ppiclf_int_given = .FALSE.


      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpFieldUser(jp,infld)
      use ppiclf_data
!
! This is called by rocpicl/PICL_TEMP_Runge.F90, once per stage for
! each gas field: infld (one value per fluid cell) is interpolated to
! rprop slot jp if ppiclf_user_Layout listed jp (ppiclf_int_slots), and
! ignored otherwise, because no model that is on reads it. Each listed
! slot must be given once per stage.
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 jp !rprop index
      REAL*8 infld(*) !value to set rprop to
!
! Internal:
!
      INTEGER*4 n, k, kfld
!
      IF(jp .LE. 0 .OR. jp .GT. ppiclf_lrp)
     >   CALL ppiclf_exittr('Interpolated slot outside rprop$',0.0d0,jp)
      kfld = 0
      DO k=1,ppiclf_lrp_int
        IF(ppiclf_int_slots(k) .EQ. jp) kfld = k
      END DO
      ! not needed by this run
      IF(kfld .EQ. 0) RETURN
      IF(ppiclf_int_given(kfld))
     >   CALL ppiclf_exittr('Interpolated field given twice, slot$'
     >                     ,0.0d0,jp)
      ppiclf_int_given(kfld) = .TRUE.
      ppiclf_int_ngiven = ppiclf_int_ngiven + 1

      ! copy to infld internal storage
      n = ppiclf_nFVCells
      CALL ppiclf_copy(ppiclf_int_fld_input(1,kfld)
     >                ,infld(1),n)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitInterp
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 ie
!
      IF(.NOT.ppiclf_overlap)
     >CALL ppiclf_exittr('Cannot interpolate unless overlap grid$',0.0d0
     >                   ,0)
      IF(.NOT.ppiclf_lintp) 
     >CALL ppiclf_exittr('No interpolated fields (ppiclf_user_Layout)$'
     >                   ,0.0d0,0)
      ppiclf_nCells_Interp = ppiclf_nCells_FV2PICL_Orig
      DO ie=1,ppiclf_nCells_Interp
        CALL ppiclf_icopy(ppiclf_cell_map_interp(1,ie),
     >        ppiclf_cell_map_Orig(1,ie), PPICLF_LRMAX)
      END DO
      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpField(j)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4 jp
!
! Internal:
!
      INTEGER*4 n, ie, iee, j
!
      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      DO ie=1,ppiclf_nCells_Interp
         ! iee is the Rocflu element ID from previous MapOverlapGrid
         ! subroutine
         ! j is the rprop index
         iee = ppiclf_cell_map_interp(1,ie) 
         CALL ppiclf_copy(ppiclf_int_fld (j,ie)
     >                   ,ppiclf_int_fld_input(iee,j),1)
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpTupleTransfer
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal: 
!
      REAL*8 Max_CellLen(3)
      INTEGER*4 nkey(2), nl, nii, njj, nrr, ie, l, nin
      LOGICAL partl
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 3
      nrr  = ppiclf_lrp_int
      nkey(1) = 2
      nkey(2) = 1

      ! Grow the cell arrays to what this rank will receive
      CALL ppiclf_comm_NumIncoming(ppiclf_nCells_Interp, nii,
     >                             ppiclf_cell_map_interp, njj, nin)
      CALL ppiclf_data_GrowCell(nin)

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Interp, ppiclf_maxcell ! Columns allocated
     >      ,ppiclf_cell_map_interp, nii      ! Integer communication
     >      ,partl, nl                        ! Logical communication
     >      ,ppiclf_int_fld, nrr              ! Real communication
     >      ,njj)                             ! Proc index to send to
      CALL pfgslib_crystal_tuple_sort(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Interp             ! Amount of columns to sort
     >      ,ppiclf_cell_map_interp,nii       ! Integer data
     >      ,partl,nl                         ! Logical data
     >      ,ppiclf_int_fld,nrr               ! Real data
     >      ,nkey,2)                          ! Sorting order

      ! Find distance check for interpolation (subbin size).
      ! This is 1.8*MaxCellLength, the support radius of the
      ! least-squares interpolation weights. The 27-cell projection
      ! stencil keeps its 1.5*MaxCellLength box (dchk/1.2).
      Max_CellLen(1) = 0.0D0
      Max_CellLen(2) = 0.0D0
      Max_CellLen(3) = 0.0D0
      DO ie = 1,ppiclf_nCells_Interp ! Loop through cells mapped to bin
        DO l = 1,3
          ! Find max cell lengths in all dimensions
          IF(ppiclf_picl_grid(3+l,ie) .GT. Max_CellLen(l))
     >      Max_CellLen(l) = ppiclf_picl_grid(3+l,ie)
        END DO !l
      END DO !ie

      DO l = 1,3
        ppiclf_interp_dchk(l) = Max_CellLen(l)*1.8D0
      END DO

      RETURN
      END

!
!______________________________________________________________________
!
      SUBROUTINE ppiclf_solve_SBParticleToCellMap
      use ppiclf_data

      IMPLICIT NONE

      INCLUDE "PPICLF"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          CellID_nearest(28), partCount
      REAL*8    dSQl, dSQi, dSQ(28), xp(3), dl, dNSQi, dNSQ(28),
     >          CellCenter(3,28), w(27), binblength(3),  
     >          Max_CellLen(3), Max_CellLenSQ(3), dSQchk(3)
      LOGICAL   added, farAway, alreadyMapped
 
      ! Cells of subbin s: SBin_list(SBin_start(s):SBin_start(s+1)-1),
      ! in the order they were mapped. A cell sits in up to 8 subbins
      ! (2 per direction with periodicity). Entries are collected first
      ! (entSB, entID), then sorted by subbin keeping their order.
      INTEGER*4  SBin_start( 0 : (
     > (FLOOR((ppiclf_bins_dx(1) + 2*ppiclf_interp_dchk(1))
     >                       /ppiclf_interp_dchk(1)) + 1) *
     > (FLOOR((ppiclf_bins_dx(2) + 2*ppiclf_interp_dchk(2))
     >                       /ppiclf_interp_dchk(2)) + 1) *
     > (FLOOR((ppiclf_bins_dx(3) + 2*ppiclf_interp_dchk(3))
     >                       /ppiclf_interp_dchk(3)) + 1) ))
      INTEGER*4  SBin_list(8*ppiclf_nCells_Interp),
     >           entSB(8*ppiclf_nCells_Interp),
     >           entID(8*ppiclf_nCells_Interp), nent
       INTEGER*4 n_SBin(3), tot_SBin, i_SBin(3), iTemp_SBin(3),
     >           temp_SBin, iSB, jSB, kSB, loopSB, i_count,
     >           firstSB(3), lastSB(3)  
      REAL*8    bin_Min(3), x_range(3), size_SBin(3)
      LOGICAL   remove
      ! Least-squares interpolation stencil: cells whose support
      ! (1.8 cell lengths per direction) contains the particle.
      ! At most 4 per direction on a uniform hex mesh.
      INTEGER*4 LSMAX
      PARAMETER (LSMAX = 125)
      INTEGER*4 nLS, lsID(LSMAX), m, mpos
      REAL*8    lsDX(3,LSMAX), lsR(3,LSMAX), lsD2(LSMAX), dx(3)
      LOGICAL   inSupp, dupl, skip
      ! lsMark(ie) = ip once cell ie is in particle ip's stencil
      INTEGER*4, ALLOCATABLE :: lsMark(:)
      ! Largest support radius in the stencil and of the current cell,
      ! and the distance window for the same-offset check
      REAL*8    lsRmax(3), rc(3), d2tol
      ! Centre and lengths of the cell of each subbin entry, in subbin
      ! order (read contiguously by the particle loop)
      REAL*8, ALLOCATABLE :: sbCell(:,:)
      LOGICAL   perl(3)
      !***************************************************************

      IF(ppiclf_npart .LT. 1) RETURN
      IF(ppiclf_nCells_Interp .EQ. 0 . AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'ERROR: ',ppiclf_npart, 'Particles mapped to bin:'
     >         ,ppiclf_nid
        PRINT*,'No cells mapped to bin for Interpolation/Projection.'
        CALL ppiclf_exittr('Failure in particle to cell mapping',0.D0,0)
      END IF
 
      DO l = 1,3
        binblength(l) = ppiclf_binb(2*l) - ppiclf_binb((2*l)-1)
        bin_Min(l) = ppiclf_bin_pos(1,l) - ppiclf_interp_dchk(l)
        n_SBin(l) = FLOOR((ppiclf_bins_dx(l) + 2*ppiclf_interp_dchk(l))
     >                                     / ppiclf_interp_dchk(l)) + 1
        ! 27-cell projection stencil box: 1.5*MaxCellLength
        dSQchk(l) = (ppiclf_interp_dchk(l)/1.2D0)**2
        ! SB at bin min boundary
        firstSB(l) = FLOOR((ppiclf_bin_pos(1,l) - bin_Min(l))
     >               / ppiclf_interp_dchk(l))
        ! SB at bin miax boundary
        lastSB(l)  =  FLOOR((ppiclf_bin_pos(2,l) - bin_Min(l))
     >               / ppiclf_interp_dchk(l))
      END DO
 
      nent = 0
      tot_SBin = n_SBin(1)*n_SBin(2)*n_SBin(3)
     
      ! Loop through all elements to map to subbins.
      ! Particles don't need to be mapped, since the particle
      ! subbin is determined in following loop.
      DO ie = 1,ppiclf_nCells_Interp  
        DO l = 1,3
          i_SBin(l) = FLOOR((ppiclf_picl_grid(l,ie) - 
     >                bin_Min(l)) / ppiclf_interp_dchk(l))
        END DO
        ! In the i,j,k loops below, 0 takes care of non-periodic mapping
        ! and 1 takes care of periodic mapping.  If a cell is in corner,
        ! it'll be mapped to 2*2*2=8 subbins.
!*** need to adjust for angular periodicity ***
        DO i = 0,1
          IF(i .EQ. 0) THEN
            iTemp_SBin(1) = i_SBin(1)
            IF(iTemp_SBin(1) .LT. 0) 
     >        iTemp_SBin(1) = 0
            IF(iTemp_SBin(1) .GT. n_SBin(1) - 1) 
     >        iTemp_SBin(1) = n_SBin(1) - 1
          ELSE ! i .EQ. 1
            IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)) THEN
              IF(i_SBin(1) .LE. firstSB(1) + 1) THEN
                iTemp_SBin(1) = lastSB(1) - 1
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE 
              ELSE IF(i_SBin(1) .GE. lastSB(1) - 1) THEN
                iTemp_SBin(1) = firstSB(1) + 1
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE
              ELSE
                CYCLE
              END IF
            ELSE 
              CYCLE
            END IF
          END IF
          DO j = 0,1
            IF(j .EQ. 0) THEN
              iTemp_SBin(2) = i_SBin(2)
              IF(iTemp_SBin(2) .LT. 0) 
     >          iTemp_SBin(2) = 0
               IF(iTemp_SBin(2) .GT. n_SBin(2) - 1) 
     >          iTemp_SBin(2) = n_SBin(2) - 1
            ELSE ! j .EQ. 1
              ! This takes care of periodicity for single processor
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)) THEN
                IF(i_SBin(2) .LE. firstSB(2) + 1) THEN
                  iTemp_SBin(2) = lastSB(2) - 1
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE IF(i_SBin(2) .GE. lastSB(2) - 1) THEN
                  iTemp_SBin(2) = firstSB(2) + 1
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE
                  CYCLE
                END IF
              ELSE 
                CYCLE
              END IF
            END IF
            DO k = 0,1
              IF(k .EQ. 0) THEN
                iTemp_SBin(3) = i_SBin(3)
                IF(iTemp_SBin(3) .LT. 0) 
     >            iTemp_SBin(3) = 0
                IF(iTemp_SBin(3) .GT. n_SBin(3) - 1) 
     >            iTemp_SBin(3) = n_SBin(3) - 1
              ELSE ! k .EQ. 1
                ! This takes care of periodicity for single processor
                IF(ppiclf_linperiodic(3) .AND. 
     >                                      ppiclf_EqualDomain(3)) THEN 
                  IF(i_SBin(3) .LE. firstSB(3) + 1) THEN
                    iTemp_SBin(3) = lastSB(3) - 1
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE IF(i_SBin(3) .GE. lastSB(3) - 1) THEN
                    iTemp_SBin(3) = firstSB(3) + 1
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE
                    CYCLE
                  END IF
                ELSE 
                  CYCLE
                END IF
              END IF
              ! Finally, add the cell to a subbin 
              temp_SBin = iTemp_SBin(1) + iTemp_SBin(2)*n_SBin(1) +
     >                    iTemp_SBin(3)*n_SBin(1)*n_SBin(2)
              nent = nent + 1
              entSB(nent) = temp_SBin
              entID(nent) = ie
            END DO !k
          END DO !j 
        END DO !i
      END DO !ie
      CALL ppiclf_solve_SortSubbins(nent,entSB,entID,tot_SBin,
     >                              SBin_start,SBin_list)
      ALLOCATE(lsMark(ppiclf_nCells_Interp), sbCell(6,MAX(nent,1)))
      lsMark = 0
      DO i_count = 1,nent
        ie = SBin_list(i_count)
        DO l = 1,6
          sbCell(l,i_count) = ppiclf_picl_grid(l,ie)
        END DO
      END DO
      DO l = 1,3
        perl(l) = ppiclf_linperiodic(l) .AND. ppiclf_EqualDomain(l)
      END DO
      partCount = 0
      DO ip=1,ppiclf_npart !Loop all particles in this bin
        nnearest = 0 ! number of nearest elements
        nLS = 0
        DO l=1,3
          lsRmax(l) = 0.0D0
        END DO
        DO ie = 1,28
          CellID_nearest(ie) = -1 ! index of nearest elements
          dSQ(ie) = 1E20 ! distance to center of nearest element
          dNSQ(ie) = 1E20
        END DO !ie
        ! particle centers in all directions
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        DO l = 1,3
          i_SBin(l) = FLOOR((xp(l) - bin_Min(l)) 
     >                /ppiclf_interp_dchk(l))
        END DO
        temp_SBin = i_SBin(1) + i_SBin(2)*n_SBin(1) +
     >              i_SBin(3)*n_SBin(1)*n_SBin(2)
        DO iSB = 1,3      ! -1,+0,+1 subbin in x-dir
          DO jSB = 1,3    ! -1,+0,+1 subbin in y-dir
            DO kSB = 1,3  ! -1,+0,+1 subbin in z-dir
              loopSB = temp_SBin + (-2+iSB) + (-2+jSB)*n_SBin(1)
     >                 + (-2 + kSB)*n_SBin(1)*n_SBin(2)  
              IF(loopSB .GT. (-1) .AND. loopSB .LT. tot_SBin) THEN
                DO i_count = SBin_start(loopSB),SBin_start(loopSB+1)-1
                  ie = SBin_list(i_count) 
                  ! get offset from particle to center (nearest
                  ! periodic image)
                  farAway = .FALSE.
                  inSupp  = .TRUE.
                  skip    = .FALSE.
                  DO l=1,3
                    dl = sbCell(l,i_count) - xp(l)
                    IF (perl(l)) THEN
                      IF(ABS(dl) .GT. binblength(l)-ABS(dl))
     >                  dl = dl - SIGN(binblength(l),dl)
                    END IF
                    dx(l) = dl
                    IF (dl**2 .GT. dSQchk(l)) farAway = .TRUE.
                    rc(l) = 1.8D0*sbCell(3+l,i_count)
                    IF (ABS(dl) .GE. rc(l)) inSupp = .FALSE.
                    ! Outside both the least-squares support and the
                    ! projection box: the cell is not used at all
                    IF (farAway .AND. .NOT. inSupp) THEN
                      skip = .TRUE.
                      EXIT
                    END IF
                  END DO !l
                  IF (skip) CYCLE !i_count
                  dSQi    = 0.0
                  dNSQi   = 0.0D0
                  DO l=1,3
                    dSQl = dx(l)**2
                    dSQi = dSQi + dSQl
                    ! Offset in units of the cell length (projection)
                    dNSQi = dNSQi + dSQl/sbCell(3+l,i_count)**2
                  END DO !l
                  ! Add to the least-squares stencil, once per cell
                  ! (a cell can sit in two subbins with periodicity,
                  ! or be received twice for a two-bin periodic box)
                  IF (inSupp) THEN
                    ! Once per cell: skip a cell already in the stencil
                    ! (a cell can sit in two subbins with periodicity),
                    ! and a copy of one at the same offset (within
                    ! 1e-9 of the support; a cell received twice for a
                    ! two-bin periodic box).
                    dupl = (lsMark(ie) .EQ. ip)
                    ! Find the place that keeps the stencil sorted by
                    ! distance (then offset), so the sums do not depend
                    ! on the order cells arrive in (domain
                    ! decomposition). A same-offset copy has almost the
                    ! same distance: only entries within d2tol of dSQi
                    ! are compared (d2tol bounds the change in distance
                    ! that offsets within 1e-9 of the support allow).
                    IF (.NOT. dupl) THEN
                      d2tol = 0.0D0
                      DO l=1,3
                        d2tol = d2tol + MAX(lsRmax(l),rc(l))**2
                      END DO
                      d2tol = 4.0D-9*d2tol
                      mpos = nLS + 1
                      DO m = nLS,1,-1
                        IF (lsD2(m) .LT. dSQi - d2tol) EXIT
                        IF (ABS(lsD2(m)-dSQi) .LE. d2tol) THEN
                          IF (ABS(lsDX(1,m)-dx(1)) .LT. 1.0D-9*lsR(1,m)
     >                  .AND. ABS(lsDX(2,m)-dx(2)) .LT. 1.0D-9*lsR(2,m)
     >                  .AND. ABS(lsDX(3,m)-dx(3)) .LT. 1.0D-9*lsR(3,m))
     >                      dupl = .TRUE.
                        END IF
                        ! entry m goes after the new cell?
                        IF (lsD2(m) .GT. dSQi) THEN
                          mpos = m
                        ELSE IF (lsD2(m) .EQ. dSQi) THEN
                          IF (lsDX(1,m) .GT. dx(1) .OR.
     >                       (lsDX(1,m) .EQ. dx(1) .AND.
     >                        lsDX(2,m) .GT. dx(2)) .OR.
     >                       (lsDX(1,m) .EQ. dx(1) .AND.
     >                        lsDX(2,m) .EQ. dx(2) .AND.
     >                        lsDX(3,m) .GE. dx(3))) mpos = m
                        END IF
                      END DO
                    END IF
                    IF (.NOT. dupl) THEN
                      nLS = nLS + 1
                      IF (nLS .GT. LSMAX) CALL ppiclf_exittr(
     >                  'Too many cells in interpolation stencil$',
     >                  0.0D0,nLS)
                      DO m = nLS,mpos+1,-1
                        lsID(m) = lsID(m-1)
                        lsD2(m) = lsD2(m-1)
                        DO l=1,3
                          lsDX(l,m) = lsDX(l,m-1)
                          lsR(l,m)  = lsR(l,m-1)
                        END DO
                      END DO
                      m = mpos
                      lsID(m) = ie
                      lsD2(m) = dSQi
                      DO l=1,3
                        lsDX(l,m) = dx(l)
                        lsR(l,m)  = rc(l)
                        lsRmax(l) = MAX(lsRmax(l),rc(l))
                      END DO
                      lsMark(ie) = ip
                    END IF
                  END IF
                  ! skip to next fluid cell if greater than 1.5*max cell
                  ! distance in respective x,y,z direction.
                  IF (farAWAY) CYCLE !i_count
                  ! Sort closest fluid cell centers
                  added = .FALSE.
                  alreadyMapped = .FALSE.
                  DO i = 1,27
                    IF(ie .EQ. CellID_nearest(i)) THEN
                      alreadyMapped = .TRUE.
                      EXIT
                    END IF
                  END DO
                  DO i = 1,27
                    IF(alreadyMapped) EXIT ! go to next cell in SB
                    j = 27 - i + 1
                    IF (dSQi .LT. dSQ(j)) THEN
                      dSQ(j+1) = dSQ(j)
                      dNSQ(j+1) = dNSQ(j)
                      CellID_nearest(j+1) = CellID_nearest(j)
                      DO l=1,3
                        CellCenter(l, j+1) = CellCenter(l, j)
                      END DO
                      dSQ(j) = dSQi
                      dNSQ(j) = dNSQi
                      CellID_nearest(j) = ie
                      DO l=1,3
                        CellCenter(l,j) = ppiclf_picl_grid(l,ie)
                      END DO
                      added = .TRUE.
                    ELSE ! If not within closest cell list
                      EXIT !i
                    END IF
                  END DO !i
                  IF (added) nnearest = nnearest + 1  
                END DO ! i_count
              END IF !SB out of domain
            END DO !kSB
          END DO !jSB
        END DO !iSB
        nnearest = MIN(nnearest,27)
        remove = .FALSE.
        IF(nnearest .LT. 1) remove = .TRUE.

        ! Thierry - 03/11/2026 - This section is not correct 
        ! as it was deleting erroneously particles. Commenting out for
        ! now.

!        ! Remove particle outside of any cell
!        IF(nnearest .LT. 27 .AND. nnearest .GE. 1) THEN
!         ie = CellID_nearest(1)
!         DO l = 1,3
!          ! similar distance check as above, but only for the cell that is
!          ! closest to the particle (1st index in nnearest)
!          ! used ABS since distance isn't squared.
!          ! ppiclf_interp_dchk is set to be 1.5xmax cell length per
!          ! dimension (in SUBROUTINE ppiclf_solve_InterpTupleTransfer)
!          IF(ppiclf_linperiodic(l) .AND.
!     >                            ppiclf_EqualDomain(l)) THEN
!            dl = ABS(MIN((ppiclf_picl_grid(l,ie) - xp(l)), 
!     >             (binblength(l)-ABS(ppiclf_picl_grid(l,ie)
!     >              - xp(l)))))
!          ELSE
!            dl = ABS(ppiclf_picl_grid(l,ie) - xp(l))
!          END IF
!          ! Ensure particle is within 1/2 cell distance of one cell.
!          IF(dl .GT. ppiclf_interp_dchk(l)/1.5D0*0.5D0) remove = .TRUE.
!         END DO
!        END IF

        IF (remove) THEN
          ! Particle is outside of fluid domain.
          ! iprop(8,ip) set to -1 means it will be removed
          ! from ppiclf_y, ppiclf_rprop and rprop3
          ppiclf_iprop(9,ip) = -1
          ppiclf_remove_particle = .TRUE.
          PRINT*, 'part # on proc # removed.', ppiclf_iprop(1,ip),
     >            ppiclf_nid
          remove = .FALSE.
        ELSE
          partCount = partCount + 1
          ! use partCount since ip includes possible removed particles
          ppiclf_nPart2Cell(partCount) = nnearest
          DO i = 1,nnearest
            ppiclf_Part2Cell_map(partCount,i) = CellID_nearest(i) ! Cell ID
            ! Squared particle to cell center offset, scaled by the
            ! cell length in each direction
            ppiclf_Part2Cell_dist(partCount,i) = dNSQ(i)
          END DO !i
          ! A particle outside the mesh may have no cell support
          ! containing it: use the nearest cell value.
          IF (nLS .EQ. 0) THEN
            nLS = 1
            lsID(1) = CellID_nearest(1)
            DO l=1,3
              lsDX(l,1) = 0.0D0
              lsR(l,1)  = 1.0D0
            END DO
          END IF
          ! Interpolate fluid fields to the particle (index ip, before
          ! removed particles are compacted)
          CALL ppiclf_solve_InterpLS(ip,nLS,lsID,lsDX,lsR)
        END IF !nnearest
      END DO !ip

      DEALLOCATE(lsMark, sbCell)

      IF(ppiclf_remove_particle) THEN
        ! Delete particles that are outside of fluid grid
        CALL ppiclf_solve_RemoveParticle
        ppiclf_remove_particle = .FALSE.
      END IF

      RETURN
      END
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_solve_InterpLS(ip,n,id,dx,rs)
      use ppiclf_data
!
! Interpolates the fluid fields to particle ip by a weighted linear
! least-squares fit f ~ a + b.(x - xp) over the n stencil cells.
! Cell c has weight prod_l (1 - (dx_l/rs_l)^2)^6, which vanishes at
! its support radius rs_l = 1.8 cell lengths, so cells enter and leave
! the stencil with zero weight and the interpolated value is
! continuous in xp. The fit is exact for linear fields. Directions
! where the stencil has no spread (one-cell-thick quasi-2D meshes) are
! dropped from the fit. Each value is clamped to the stencil min/max.
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 ip, n, id(n)
      REAL*8    dx(3,n), rs(3,n)
!
! Internal:
!
      INTEGER*4 i, j, k, l, c, r, nc, ipiv, col(4), nfld
      REAL*8    w(n), s(3,n), phi(n), A(4,5), hs(3), wsum, sm(3),
     >          var(3), ak(4), q, val, fmin, fmax, f, tmp, amax
!
      ! Scale offsets by the largest support radius (conditioning)
      DO l = 1,3
        hs(l) = 0.0D0
        DO k = 1,n
          hs(l) = MAX(hs(l),rs(l,k))
        END DO
      END DO

      wsum = 0.0D0
      DO l = 1,3
        sm(l) = 0.0D0
      END DO
      DO k = 1,n
        w(k) = 1.0D0
        DO l = 1,3
          q = dx(l,k)/rs(l,k)
          w(k) = w(k)*(1.0D0 - q*q)**6
          s(l,k) = dx(l,k)/hs(l)
          sm(l) = sm(l) + w(k)*s(l,k)
        END DO
        wsum = wsum + w(k)
      END DO

      ! Keep the constant and every direction with weighted spread
      nc = 1
      col(1) = 0
      DO l = 1,3
        sm(l) = sm(l)/wsum
        var(l) = 0.0D0
        DO k = 1,n
          var(l) = var(l) + w(k)*(s(l,k) - sm(l))**2
        END DO
        var(l) = var(l)/wsum
        IF (var(l) .GT. 1.0D-8) THEN
          nc = nc + 1
          col(nc) = l
        END IF
      END DO

      ! Normal equations M g = e1, M = sum_k w_k a_k a_k^T
      DO r = 1,nc
        DO c = 1,nc+1
          A(r,c) = 0.0D0
        END DO
      END DO
      A(1,nc+1) = 1.0D0
      DO k = 1,n
        ak(1) = 1.0D0
        DO c = 2,nc
          ak(c) = s(col(c),k)
        END DO
        DO r = 1,nc
          DO c = 1,nc
            A(r,c) = A(r,c) + w(k)*ak(r)*ak(c)
          END DO
        END DO
      END DO

      ! Gaussian elimination with partial pivoting
      DO c = 1,nc
        ipiv = c
        amax = ABS(A(c,c))
        DO r = c+1,nc
          IF (ABS(A(r,c)) .GT. amax) THEN
            ipiv = r
            amax = ABS(A(r,c))
          END IF
        END DO
        IF (amax .LE. 1.0D-14*wsum) CALL ppiclf_exittr(
     >    'Singular least-squares interpolation matrix$',0.0D0,ip)
        IF (ipiv .NE. c) THEN
          DO j = c,nc+1
            tmp = A(c,j)
            A(c,j) = A(ipiv,j)
            A(ipiv,j) = tmp
          END DO
        END IF
        DO r = c+1,nc
          tmp = A(r,c)/A(c,c)
          DO j = c,nc+1
            A(r,j) = A(r,j) - tmp*A(c,j)
          END DO
        END DO
      END DO
      DO r = nc,1,-1
        tmp = A(r,nc+1)
        DO c = r+1,nc
          tmp = tmp - A(r,c)*A(c,nc+1)
        END DO
        A(r,nc+1) = tmp/A(r,r)
      END DO

      ! Interpolation weights phi_k = w_k a_k^T g (sum to one)
      DO k = 1,n
        tmp = A(1,nc+1)
        DO c = 2,nc
          tmp = tmp + A(c,nc+1)*s(col(c),k)
        END DO
        phi(k) = w(k)*tmp
      END DO

      ! no field given yet (start-up): nothing to interpolate
      nfld = 0
      IF(ppiclf_int_ngiven .GT. 0) nfld = ppiclf_lrp_int
      DO i = 1,nfld
        j = ppiclf_int_slots(i)
        val  = 0.0D0
        fmin = ppiclf_int_fld(i,id(1))
        fmax = fmin
        DO k = 1,n
          f = ppiclf_int_fld(i,id(k))
          val = val + phi(k)*f
          fmin = MIN(fmin,f)
          fmax = MAX(fmax,f)
        END DO ! k
        IF (isnan(val)) THEN
          PRINT *, 'INTERP NAN: Particle, processor id, ncells', ip,
     >                                   ppiclf_nid,n
          PRINT*, 'Index:',j, 'Value:',val
          CALL ppiclf_exittr('rprop NaN in Interpolate$',0.D0,0)
        END IF
        ppiclf_rprop(j, ip) = MIN(fmax,MAX(fmin,val))
      END DO ! i

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_RemoveParticle
      use ppiclf_data
!
      USE ppiclf_unsteady
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

! Internal:
!
      INTEGER*4 i, icount
      icount = 0
      DO i=1,ppiclf_npart
         IF(ppiclf_iprop(9,i) .NE. -1) THEN
            ! Keep particle - copy the column with particle information
            icount = icount + 1 
            IF(i .NE. icount) THEN
               CALL ppiclf_copy
     >          (ppiclf_y     (1,icount),     ppiclf_y(1,i), ppiclf_lrs)
               CALL ppiclf_copy
     >          (ppiclf_y1    (1,icount),    ppiclf_y1(1,i), ppiclf_lrs)
               CALL ppiclf_copy
     >          (ppiclf_ydot  (1,icount),  ppiclf_ydot(1,i), ppiclf_lrs)
               CALL ppiclf_copy
     >          (ppiclf_ydotc (1,icount), ppiclf_ydotc(1,i), ppiclf_lrs)
               CALL ppiclf_copy
     >          (ppiclf_rprop (1,icount), ppiclf_rprop(1,i), ppiclf_lrp)
               IF(ppiclf_lrp3 .GT. 1) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop3(1,icount),ppiclf_rprop3(1,i),ppiclf_lrp3)
               END IF
               CALL ppiclf_copy(ppiclf_feedbk(1,icount), 
     >                          ppiclf_feedbk(1,i), ppiclf_lrp_pro)
               CALL ppiclf_icopy
     >          (ppiclf_iprop(1,icount) , ppiclf_iprop(1,i), PPICLF_LIP)
            END IF
         ! Else - don't copy particle column if marked for removal
         ! Particles marked for removal if outside fluid domain, which
         ! is found in the particle to cell mapping during interpolation
         END IF
      END DO

      ppiclf_npart = icount

      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_ProjectParticleGrid
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

      ! Internal:
      INTEGER*4 i, j, ip, ie, nCellProj, CellID, nl, nii, njj,
     >          nrr, nkey(2), iee, nin
      REAL*8    CellVol, GaussianConst, dist, w(27), wsum
      LOGICAL   partl 
 
      GaussianConst = 2.305D0 ! Distribution over 2 cell widths
      ppiclf_pro_fld_picl = 0.0d0
      DO ip=1,ppiclf_npart
        ! Update volume fraction for feedback - important for 1st RK
        ! step at time = 0.0
        ppiclf_feedbk(PPICLF_P_JPHIP,ip) =
     >  ppiclf_rprop(PPICLF_R_JVOLP,ip) * ppiclf_rprop(PPICLF_R_JSPL,ip)
        nCellProj = ppiclf_nPart2Cell(ip)
        wsum = 0.0D0
        ! Loop to find individual cell weightings. As in MFC, the
        ! offset in each direction is scaled by the receiving cell's
        ! length in that direction, so elongated cells spread the
        ! particle along their long side too.
        DO i = 1,nCellProj
          CellID = ppiclf_Part2Cell_map(ip,i)
          dist = ppiclf_Part2Cell_dist(ip,i)
          CellVol = ppiclf_picl_grid(7,CellID)
          w(i) = ABS(CellVol*EXP(-GaussianConst*dist))
          wsum = wsum + w(i)
        END DO !i
        DO j=1,ppiclf_lrp_pro
          ! Loop through cells to apply feedback     
          DO i = 1,nCellProj
            CellID = ppiclf_Part2Cell_map(ip,i)
            ppiclf_pro_fld_picl(j,CellID) = 
     >         ppiclf_pro_fld_picl(j,CellID) 
     >         + ppiclf_feedbk(j,ip)*w(i)/wsum
          END DO !i
        END DO !j
      END DO !ip

      ! Now send feedback information to processor that contains 
      ! the cell for the fluid solver

      ppiclf_nCells_Proj = ppiclf_nCells_Interp
      DO i = 1,ppiclf_nCells_Proj
        CALL ppiclf_icopy(ppiclf_cell_map_proj(1,i),
     >         ppiclf_cell_map_interp(1,i),PPICLF_LRMAX)
      END DO

      nl = 0
      nii = PPICLF_LRMAX
      njj = 2 ! original processor with cell for fluid grid
      nrr = ppiclf_lrp_pro
      nkey(1) = 2
      nkey(2) = 1
      ! Grow the cell arrays to what this rank will receive
      CALL ppiclf_comm_NumIncoming(ppiclf_nCells_Proj, nii,
     >                             ppiclf_cell_map_proj, njj, nin)
      CALL ppiclf_data_GrowCell(nin)

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Proj, ppiclf_maxcell ! Columns allocated
     >      ,ppiclf_cell_map_proj, nii      ! Integer communication
     >      ,partl, nl                      ! Logical communication
     >      ,ppiclf_pro_fld_picl, nrr       ! Real communication
     >      ,njj)                           ! Proc index to send to
      CALL pfgslib_crystal_tuple_sort(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Proj             ! Amount of columns to sort
     >      ,ppiclf_cell_map_proj,nii       ! Integer data
     >      ,partl,nl                       ! Logical data
     >      ,ppiclf_pro_fld_picl,nrr        ! Real data
     >      ,nkey,2)                        ! Sorting order

      ppiclf_pro_fld = 0.0d0
      DO ie=1,ppiclf_nCells_Proj
         iee = ppiclf_cell_map_Proj(1,ie)
         DO j=1,ppiclf_lrp_pro
           ! Mapped to the fluid solver domain
           ppiclf_pro_fld(iee,j) = ppiclf_pro_fld(iee,j) +
     >                                    ppiclf_pro_fld_picl(j,ie)
         END DO
      END DO

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_GetProFld(e,m,fld)
      use ppiclf_data
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 e,m
! e - fluid_grid element number
! m - projection property index

!
! Output:
!
      REAL*8 fld
!
      fld = ppiclf_pro_fld(e,m)

      RETURN
      END

!-----------------------------------------------------------------------
!       SUBROUTINE ppiclf_solve_InvokeAngularPeriodic(i,flag,
!     >                                              per_alpha,
!     >                                              angle, xangle,
!     >                                              register)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! :
!! Input: 
!! 
!      ! Thierry -  07/24/24 - modified code begins here
!      ! global variables
!      INTEGER*4 i, flag
!      REAL*8 rin, rout, per_alpha, angle, xangle
!      ! local variables
!      REAL*8 ct, st, ex, ey, ez, local_angle
!      REAL*8 rotmat(3,3) , v(3), x(3)
!      INTEGER*4 register
!
!        ! Thierry - 07/24/24 - modified code begins here
!        ! Implementation of Rotational Periodicity
!        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
!        ! this is invoked when particle is leaving x-axis or y-axis
!       
!!        print*, "!!! Rotational Periodicity Invoked !!!!" 
!          
!        ! use local angle so the value of angle does not get affected globally
!        local_angle = angle
!        ! Thierry 
!        !    (1) sign convention for theta is +ve when measured CCW
!        !           switch angle sign when particle is leaving from 
!        !           upper periodic face
!        !    (2) 0.5 instead of 1.0 to switch angle for ghost algorithm
!        !           since the ghost is being created before the 
!        !           particle is leaving domain
!        if(per_alpha.gt. 0.5*(xangle+angle)) 
!     >    local_angle=-1.0*local_angle
!        
!        ! Half-cylinder case - particle leaving +ve x-axis 
!        !                    - adjust rotation matrix angle accordingly
!        if(ang_case .EQ. 3) then
!          if(per_alpha .lt. xangle) local_angle = 0.0 
!        END if
!
!        ! convert from degrees to radians
!        ct = cos(local_angle)
!        st = sin(local_angle)
!        
!        SELECT CASE(flag)
!          !CASE(1)
!          !  ex = 1.0d0
!          !  ey = 0.0d0
!          !  ez = 0.0d0
!          !  print*, "X-Rotational Axis"
!
!          !CASE(2)
!          !  ex = 0.0d0
!          !  ey = 1.0d0
!          !  ez = 0.0d0
!          !  print*, "Y-Rotational Axis"
!
!          CASE(1)
!            ex = 0.0d0
!            ey = 0.0d0
!            ez = 1.0d0
!!            print*, "Z-Rotational Axis"
!          CASE DEFAULT
!            CALL ppiclf_exittr('Invalid Axis of Rotation!$',0.0d0
!     >         ,ppiclf_nid)
!
!          END SELECT 
!          
!          ! Rotation Matrix calculation
!          rotmat(1,1) = ct + (1.0d0-ct)*ex*ex
!          rotmat(1,2) =      (1.0d0-ct)*ex*ey - st*ez
!          rotmat(1,3) =      (1.0d0-ct)*ex*ez + st*ey
!          
!          rotmat(2,1) =      (1.0d0-ct)*ey*ex + st*ez
!          rotmat(2,2) = ct + (1.0d0-ct)*ey*ey
!          rotmat(2,3) =      (1.0d0-ct)*ey*ez - st*ex
!          
!          rotmat(3,1) =      (1.0d0-ct)*ez*ex - st*ey
!          rotmat(3,2) =      (1.0d0-ct)*ez*ey + st*ex
!          rotmat(3,3) = ct + (1.0d0-ct)*ez*ez
!
!          ! Corrdinates modification
!          x(1) = ppiclf_y(PPICLF_JX,i)
!          x(2) = ppiclf_y(PPICLF_JY,i)
!          x(3) = ppiclf_y(PPICLF_JZ,i)
!          
!          xrot = MATMUL(rotmat, x)
!          
!          ! Velocity vector modification
!
!          v(1) = ppiclf_y(PPICLF_JVX,i)
!          v(2) = ppiclf_y(PPICLF_JVY,i)
!          v(3) = ppiclf_y(PPICLF_JVZ,i)
!
!          vrot = MATMUL(rotmat, v)
!          
!          ! 08/27/24 - Thierry - we add a register variable to 
!          !   choose if we want to register
!          !   the angularly modified variables 
!          ! register = 1 when called from RemoveParticle -> we want to modify the values
!          ! register = 0 when called from AngularCreateGhost -> we don't want to modify values
!          
!          ! register modified values
!          if (register==1) then
!            !print*, "Registering values!" 
!            ppiclf_y(PPICLF_JX,i) = xrot(1)
!            ppiclf_y(PPICLF_JY,i) = xrot(2)
!            ppiclf_y(PPICLF_JZ,i) = xrot(3)
!            
!            ppiclf_y(PPICLF_JVX,i) = vrot(1)
!            ppiclf_y(PPICLF_JVY,i) = vrot(2)
!            ppiclf_y(PPICLF_JVZ,i) = vrot(3)
!            
!          END if 
!       
!      RETURN
!      END
!-----------------------------------------------------------------------
!      SUBROUTINE ppiclf_solve_InitAngularPlane(i,rin, rout,
!     >                                         angle, xangle,
!     >                                         dist1, dist2)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! Inputs 
!!
!      INTEGER*4 i
!      REAL*8 rin, rout, angle, xangle
!! Local Variables:
!!     
!      REAL*8 p1(3), p2(3), p3(3), p4(3), p5(3), p6(3),
!     >       v1(3), v2(3), v3(3), v4(3), n1(3), n2(3),
!     >       A, B, C, D, E, F, G, H, zt, xp, yp, zp
!!
!! Outputs
!      REAL*8 dist1, dist2
!!
!      xp = ppiclf_y(PPICLF_JX,i)
!      yp = ppiclf_y(PPICLF_JY,i)
!      zp = ppiclf_y(PPICLF_JZ,i)
!      
!      zt = ppiclf_binb(6) - ppiclf_binb(5) ! bin thickness in z-direction
!
!      !!! upper plane calculation !!! 
!      
!      ! plane equation
!      ! Ax + By + Cz + D = 0
!
!      ! p1, p2, p3 are 3 points in the upper plane
!      p1 = (/rin, tan(angle - abs(xangle))*rin, 0.0d0/) 
!      p2 = (/rout, tan(angle - abs(xangle))*rout, 0.0d0/) 
!      p3 = (/rin, tan(angle - abs(xangle))*rin, zt/) 
!
!      v1 = p2 - p1 ! vector P1P2
!      v2 = p3 - p1 ! vector P1P3
!      
!      ! upper plane normal vector - n1(A,B,C) = v1 x v2
!      ! cross product calculation
!      A =  v1(2)*v2(3) - v1(3)*v2(2)
!      B = -v1(1)*v2(3) + v1(3)*v2(1)
!      C =  v1(1)*v2(2) - v1(2)*v2(1)
!      n1(1)=A ; n1(2)=B; n1(3)=C
!      
!      ! values of either p1, p2, or p3 can be used to calculate D
!      D = -A*p1(1) - B*p1(2) - C*p1(3)
!      
!      ! P(xp, yp, zp) arbitrary point
!      ! dist = distance between P and upper plane 
!      dist1 = abs(A*xp + B*yp + C*zp + D)
!      dist1 = dist1/sqrt(A**2 + B**2 + C**2)
!      
!      !!! lower plane calculation !!! 
!      ! plane equation
!      ! Ex + Fy + Gz + H = 0
!
!      ! p4, p5, p6 are 3 points in the lower plane
!      p4 = (/rin, -tan(angle - abs(xangle))*rin, 0.0d0/)
!      p5 = (/rout, -tan(angle - abs(xangle))*rout, 0.0d0/)
!      p6 = (/rin, -tan(angle - abs(xangle))*rin, zt/)
!      
!      v3 = p5 - p4 ! vector P4P5
!      v4 = p6 - p4 ! vector P4P6
!      
!      ! lower plane normal vector - n2(E,F,G) = v3 x v4
!      ! cross product calculation
!      E =  v3(2)*v4(3) - v3(3)*v4(2)
!      F = -v3(1)*v4(3) + v3(3)*v4(1)
!      G =  v3(1)*v4(2) - v3(2)*v4(1)
!      n2(1)=E ; n2(2)=F; n2(3)=G
!      
!      ! values of either p4, p5, or p6 can be used to calculate H
!      H = -E*p4(1) - F*p4(2) - G*p4(3)
!
!      ! P(xp, yp, zp) arbitrary point
!      ! dist = distance between P and lower plane 
!      dist2 = abs(E*xp + F*yp + G*zp + H)
!      dist2 = dist2/sqrt(E**2 + F**2 + G**2)
!
!      RETURN
!      END

!-----------------------------------------------------------------------
!
!
!!
!!
!! Maybe this one can be deleted??? - AVERY ***
!!
!!
!      SUBROUTINE ppiclf_solve_InitAngularPeriodic(flag,
!     >              rin, rout, angle, xangle)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Input: 
!! 
!      ! Thierry -  07/24/24 - modified code begings here
!      ! global variables - user input file
!      INTEGER*4 flag
!      REAL*8 rin, rout, angle, xangle
!      ! local variables
!      REAL*8 pi, angled
!
!        ! Thierry - 07/24/24 - modified code begins here
!        ! Implementation of Angular Periodicity
!        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
!        ! this is invoked when particle is leaving x-axis or y-axis
!       
!        ! sign convention for theta is +ve when measured CCW
!        ! switch angle sign when particle is leaving from upper face
!        if (rin .ge. rout)
!     >   CALL ppiclf_exittr('Angular Per must have rin < rout$',rout,0)
!
!            ppiclf_iperiodic(1) = 0 ! X-periodic
!            ppiclf_iperiodic(2) = 0 ! Y-periodic
!
!            SELECT CASE (ang_case)
!              CASE (1) ! general wedge ; 0 <= angle < 90
!                if (ppiclf_nid.EQ.0) print*,"General Wedge Initialized!"
!                ppiclf_xdrange(1,1) = rin  ! Min X-periodic face
!                ppiclf_xdrange(2,1) = rout ! Max X-periodic face
!                ppiclf_xdrange(1,2) = tan(xangle)*rout ! Min Y-periodic face
!                ppiclf_xdrange(2,2) = tan(angle - abs(xangle))*rout ! Max Y-periodic face
!
!              CASE (2) ! quarter cylinder ; angle = 90
!                if (ppiclf_nid.EQ.0)
!     >             print*,"Quarter Cylinder Initialized!"
!                ppiclf_xdrange(1,1) = rin  
!                ppiclf_xdrange(2,1) = rout 
!                ppiclf_xdrange(1,2) = tan(xangle)*rout
!                ppiclf_xdrange(2,2) = rout 
!              
!              CASE (3) ! half cylinder ; angle = 180
!                if (ppiclf_nid.EQ.0)
!     >             print*,"Half Cylinder Initialized!"
!                ppiclf_xdrange(1,1) = -1.0*rout
!                ppiclf_xdrange(2,1) = rout 
!                ppiclf_xdrange(1,2) = tan(xangle)*rout
!                ppiclf_xdrange(2,2) = rout 
!              
!              CASE DEFAULT
!                CALL ppiclf_exittr('Invalid Rotational Case!$',0.0d0
!     >             ,ppiclf_nid)
!              END SELECT
!
!      
!      RETURN
!      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_SortSubbins(nent,entSB,entID,nSB,
     >                                    SBstart,SBlist)
!
! Stable counting sort of nent (subbin, id) entries. On return the ids
! of subbin s are SBlist(SBstart(s):SBstart(s+1)-1), in entry order.
!
      IMPLICIT NONE
!
! Input:
!
      INTEGER*4 nent, nSB, entSB(nent), entID(nent)
!
! Output:
!
      INTEGER*4 SBstart(0:nSB), SBlist(*)
!
! Internal:
!
      INTEGER*4 e, sb, next(0:nSB-1)
!
      DO sb = 0,nSB
        SBstart(sb) = 0
      END DO
      DO e = 1,nent
        SBstart(entSB(e)+1) = SBstart(entSB(e)+1) + 1
      END DO
      SBstart(0) = 1
      DO sb = 1,nSB
        SBstart(sb) = SBstart(sb) + SBstart(sb-1)
      END DO
      DO sb = 0,nSB-1
        next(sb) = SBstart(sb)
      END DO
      DO e = 1,nent
        sb = entSB(e)
        SBlist(next(sb)) = entID(e)
        next(sb) = next(sb) + 1
      END DO

      RETURN
      END
