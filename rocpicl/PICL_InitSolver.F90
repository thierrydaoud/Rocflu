!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: 
!
! Description: none.
!
! Input: 
!
! Output:
!
! Notes: 
!
!******************************************************************************
!
! $Id: PICL_F90,v 1.0 2022/05/08 bdurant Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PICL_InitSolver( pRegion)

!  USE 

!USE ModInterfaces, ONLY: 
 USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModMPI      
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt
  USE ModBndPatch, ONLY: t_patch

!PICL
  USE ModRandom, ONLY: Rand1Uniform,Rand1Normal
  USE RFLU_ModInCellTest

  IMPLICIT NONE
#ifdef PICL
!DEC$ NOFREEFORM
#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!DEC$ FREEFORM
#endif



! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!TYPE(t_region), DIMENSION(:), POINTER :: regions
TYPE(t_global), POINTER :: global
!TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
INTEGER :: errorFlag,icg

#ifdef PICL
   CHARACTER(CHRLEN) :: endString,iFileName,comment
   CHARACTER(12) :: vtuFile,vtuFile1
   LOGICAL :: notfoundFlag, pf_restart,pf_rpInit,pf_settle,&
              wall_exists, fexists
   INTEGER :: i,npart,nCells,vi,vii,ii,jj,kk,loopCounter,ipart,icl
   INTEGER :: PPC,numPclCells,npart_local,i_global,i_global_min,i_global_max,&
              iFile, k, j, l, m
   REAL(RFREAL) :: dp_min,dp_max,rhop,tester,ratio,total_vol,xMinCurt,&
                   xMaxCurt,yMinCurt,yMaxCurt,xMinCell,xMaxCell,yMinCell,&
                   yMaxCell,zMinCell,zMaxCell,x,vFrac,volpclsum,xLoc,yLoc,zLoc,yL, &
                   zpf_factor,xpf_factor,dp,neighborWidth,xp_min,xp_max, &
                   yp_min, yp_max, zp_min, zp_max, MinFluidCells, maxVF
   ! this rank's share of points.dat, handed to ppiclF
   REAL(RFREAL), ALLOCATABLE :: y(:,:), rprop(:,:)
   REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: rocGrid 
   REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: xGrid, yGrid, zGrid
   REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: xData,yData,zData,rData,dumData    
   REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: volp,SPL 
   REAL(RFREAL), DIMENSION(3) :: tpw1,tpw2,tpw3         
   REAL(RFREAL) :: xin, wout
   ! Particle walls from PARTWALL patches
   TYPE(t_patch), POINTER :: pPatch
   INTEGER :: iPatch, ifc, iv, nv, iTri, nTriLoc, nTriTot, nProcs
   INTEGER, DIMENSION(:), ALLOCATABLE :: triCounts, triDispls, partWallPatch
   REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: triLoc, triAll
   REAL(RFREAL) :: xv(3,4)
   ! points.dat line and optional initial velocity
   CHARACTER(512) :: pointsLine
   REAL(RFREAL) :: velp(3), xyzp(3)
   INTEGER :: velErr

   INTEGER :: seed(33), isize, CellVertices

   INTEGER :: velocity_model, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, rmu_flag, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, flow_model, pseudoTurb_flag, &
        HTUnsteady_flag, out_particle_flag, out_gas_flag, &
        out_gradgas_flag, out_model_flag
   REAL(RFREAL) :: rmu_ref, tref, suth, ksp, erest
   COMMON /RFLU_ppiclF/ rmu_ref, tref, suth, ksp, erest, &
        velocity_model, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, rmu_flag, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, flow_model, pseudoTurb_flag, &
        HTUnsteady_flag, out_particle_flag, out_gas_flag, &
        out_gradgas_flag, out_model_flag
   ! entries per particle in y and rprop (ppiclf_user_Layout)
   INTEGER :: nStateSlots, nPropSlots
   REAL(RFREAL) :: ppiclf_rcp_part
   COMMON /RFLU_ppiclf_misc01/ ppiclf_rcp_part

   ! Fluid domain extent (min, max in x, y, z) and periodic directions
   REAL(RFREAL) :: domainRange(2,3)
   LOGICAL :: periodic(3)

   REAL(RFREAL), DIMENSION(3) :: MaxPoint, MinPoint, CellLen, Max_CellLen, &
                                 filter, filter_local
   
   ! 04/04/2025 - TLJ - added min/max grid for periodicity
   INTEGER :: errorFrag
   REAL(RFREAL) :: gridmin,gridmax
#endif

   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_InitSolver.F90,v $ $Revision: 1.1.1.1 $'
  
  global => pRegion%global !pRegion%global

  CALL RegisterFunction( global,'PICL_InitSolver',__FILE__ )


! Set pointers ----------------------------------------------------------------

    !levels(0)=>pLevel    
  !  pRegion => regions(1)         !pLevel%regions(iReg)
    pGrid   => pRegion%grid!pRegion%grid

!MOVING PICL HERE TO AVOID SEG FAULT OF ROCFLU STORED VF
#ifdef PICL

! Setting flags from input file (checked before the run by
! utilities/scripts/case_validator.py: RK3, viscosity, flag values)
! ************************************************************
velocity_model = global%piclVelocityModel
qs_flag = global%piclQsFlag
am_flag = global%piclAmFlag
pg_flag = global%piclPgFlag
collisional_flag = global%piclCollisionFlag
ViscousUnsteady_flag = global%piclViscousUnsteady
heattransfer_flag = global%piclHeatTransferFlag
HTUnsteady_flag = global%piclHTUnsteadyFlag
feedback_flag = global%piclFeedbackFlag
qs_fluct_flag = global%piclQsFluctFlag
sbNearest_flag = global%piclSBNearFlag
pseudoTurb_flag = global%piclPseudoTurbFlag
! optional groups of fields in the particle files
out_particle_flag = global%piclOutParticle
out_gas_flag      = global%piclOutGas
out_gradgas_flag  = global%piclOutGradGas
out_model_flag    = global%piclOutModel

rmu_flag = pRegion%mixtInput%viscModel
rmu_ref = pRegion%mixtInput%refVisc
tref = pRegion%mixtInput%refViscTemp
suth = pRegion%mixtInput%suthCoef
rmu_suth_param = VISC_SUTHR
rmu_fixed_param = VISC_FIXED
flow_model = int(pRegion%mixtInput%flowModel)

! Stop if the input needs particle slots this build does not have (the
! build's layout comes from libpicl/user_files/ppiclf_slots.txt)
#include "PPICLF_LAYOUT_CHECK.inc"

ksp = global%piclKsp
erest = global%piclERest

qs_fluct_filter_flag = global%piclQsFluctFilterFlag
qs_fluct_filter_adapt_flag = global%piclQsFluctFilterAdaptFlag

periodic(1) = global%piclPeriodicXFlag == 1
periodic(2) = global%piclPeriodicYFlag == 1
periodic(3) = global%piclPeriodicZFlag == 1
! Min/max grid coordinates across MPI ranks (l = XCOORD, YCOORD, ZCOORD)
DO l = 1,3
  gridmin = MINVAL(pGrid%xyz(l,1:pGrid%nVert))
  gridmax = MAXVAL(pGrid%xyz(l,1:pGrid%nVert))
  CALL MPI_AllReduce(gridmin,domainRange(1,l),1,MPI_RFREAL,MPI_MIN, &
          global%mpiComm,errorFlag)
  CALL MPI_AllReduce(gridmax,domainRange(2,l),1,MPI_RFREAL,MPI_MAX, &
          global%mpiComm,errorFlag)
END DO


 ! Initialization for viscous unsteady term: ppiclF allocates this many
 ! steps of history (none unless VISCOUSUNSTEADY >= 1)
 ppiclf_nTimeBH = 1
 ppiclf_nUnsteadyData = 0
 IF (global%piclViscousUnsteady >= 1) ppiclf_nUnsteadyData = global%piclVUSteps

! Field sizes and slot lists of this run (ppiclf_user_Layout, from the
! flags above)
CALL ppiclf_solve_InitLayout(nStateSlots, nPropSlots)

! Needed for fluctuations
seed = 1
CALL RANDOM_SEED(put=seed)
CALL RANDOM_SEED(size=isize)


! ************************************************************

! Josh Gillis - Fixed restart probelm
! TLJ - we should probably use the .rin file instead
global%restartFromScratch = .true.
vtuFile1 = 'par00002.vtu'
INQUIRE(FILE=trim(vtuFile1), EXIST=fexists)
IF( global%myProcid == MASTERPROC) then
   write(*,*) fexists, "fexists"
ENDIF
IF(fexists) THEN
   global%restartFromScratch = .false.
   IF ( global%myProcid == MASTERPROC) THEN
      PRINT*, " "
      PRINT*, " ======================================="
      PRINT*, " "
      WRITE(*,*) 'Starting PPICLF Restart'
   END IF
END IF

IF(global%restartFromScratch) THEN
   IF ( global%myProcid == MASTERPROC) then
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading points.dat file...'
   END IF
  
   iFileName = 'points.dat'
   iFile = 0

   ! open data file
   OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
   global%error = errorFlag   
   IF( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
   END IF

   ! check for comments at beginning of file
   ! No comments allowed after beginning of file
   comment = '!'
   DO WHILE (comment(1:1) == '!')
      READ(iFile, '(A)') comment
      comment = ADJUSTL(TRIM(comment))
   END DO

   BACKSPACE(iFile, IOSTAT=ErrorFlag)
  
   READ(iFile,*) npart ! global number of particles
  
   npart_local = CEILING(REAL(npart)/(REAL(global%nProcs)))
   ALLOCATE(y(nStateSlots,MAX(npart_local,1)), &
            rprop(nPropSlots,MAX(npart_local,1)),STAT=errorFlag)
   IF ( errorFlag /= 0 ) THEN
     CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:y,rprop')
   END IF
   y = 0.0_RFREAL
   i = 1
   i_global = 1
   i_global_min = npart_local*global%myProcid
   i_global_max = npart_local*(global%myProcid+1)
   IF(i_global_max > npart) i_global_max = npart

   rprop = 0.0_RFREAL
   dp_max = 0.0D0
   xp_min =  17400000.0
   yp_min =  17400000.0
   zp_min =  17400000.0
   xp_max = -17400000.0
   yp_max = -17400000.0
   zp_max = -17400000.0
   ! Every particle has the density and specific heat of the # PICL input
   rhop = global%piclDensity
   ppiclf_rcp_part = global%piclSpecificHeat

   ! Each line: x y z vx vy vz dp. Every line is read (for the global
   ! bounds), but only this rank's lines go into y: y has exactly
   ! npart_local columns.
   DO i_global=1,npart
      READ(iFile,'(A)') pointsLine
      READ(pointsLine,*,IOSTAT=velErr) xyzp,velp,dp
      IF (velErr /= 0) THEN
         CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
           'points.dat: a particle line must be x y z vx vy vz dp '// &
           '(./rocflu.sh validate shows which line)')
      END IF

      ! These are global max/min's since in 1:npart loop
      dp_max = max(dp_max,dp)
      xp_min = min(xp_min,xyzp(1)-dp/2.0)
      xp_max = max(xp_max,xyzp(1)+dp/2.0)
      yp_min = min(yp_min,xyzp(2)-dp/2.0)
      yp_max = max(yp_max,xyzp(2)+dp/2.0)
      zp_min = min(zp_min,xyzp(3)-dp/2.0)
      zp_max = max(zp_max,xyzp(3)+dp/2.0)
  
      ! if in range for this processor set all the other properties and increment i
      IF((i_global .GT. i_global_min) .AND. (i_global .LE. i_global_max)) THEN
         y(PPICLF_JX,i) = xyzp(1)
         y(PPICLF_JY,i) = xyzp(2)
         y(PPICLF_JZ,i) = xyzp(3)
         y(PPICLF_JVX,i) = velp(1)
         y(PPICLF_JVY,i) = velp(2)
         y(PPICLF_JVZ,i) = velp(3)
         y(PPICLF_JT, i) = global%piclTemp
#if PPICLF_JOX > 0
         y(PPICLF_JOX,i) = 0.0D0
         y(PPICLF_JOY,i) = 0.0D0
         y(PPICLF_JOZ,i) = 0.0D0
#endif
  
         IF ( global%myProcid == MASTERPROC) then
            IF (i==1) THEN
               print*
               print*,'   Density: ',rhop
               print*,'   C_p:     ',ppiclf_rcp_part
               print*,'   T_p:     ',global%piclTemp
               print*
            END IF
         END IF
    
         ! now set properties that are not interpolated from Rocflu onto the particles
         rprop(PPICLF_R_JRHOP,i) = rhop   ! particle density
         rprop(PPICLF_R_JDP,i)   = dp ! particle diameter
         rprop(PPICLF_R_JVOLP,i) = (4.0_RFREAL/3.0_RFREAL)*global%pi*&
                                   (0.5_RFREAL*dp)**3 ! particle volume
         ! Super Particle Loading (Real Number of particles = SPL * number of compuational particles)
         rprop(PPICLF_R_JSPL,i) = 1.0_RFREAL

         i = i + 1
      END IF
   END DO
   npart_local = i - 1

   CALL MPI_ALLREDUCE(xp_min,xp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(xp_max,xp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(yp_min,yp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(yp_max,yp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(zp_min,zp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(zp_max,zp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   IF(global%myProcid == MASTERPROC) THEN
      print*
      print*,'Particle domain boundaries at t=0'
      print*,'x - min, max, dx', xp_min, xp_max, xp_max - xp_min
      print*,'y - min, max, dx', yp_min, yp_max, yp_max - yp_min
      print*,'z - min, max, dx', zp_min, zp_max, zp_max - zp_min
      print*
   END IF


   IF(xp_min < domainRange(1,1) .OR. xp_max > domainRange(2,1) .OR. &
      yp_min < domainRange(1,2) .OR. yp_max > domainRange(2,2) .OR. &
      zp_min < domainRange(1,3) .OR. zp_max > domainRange(2,3)) THEN
     IF(global%myProcid == MASTERPROC) THEN
       WRITE(*,*) 'WARNING - Particles initalized outside of fluid domain'
       WRITE(*,*) 'Particle domain boundaries at t=0'
       WRITE(*,*) 'x - min, max, dx', xp_min, xp_max, xp_max - xp_min
       WRITE(*,*) 'y - min, max, dx', yp_min, yp_max, yp_max - yp_min
       WRITE(*,*) 'z - min, max, dx', zp_min, zp_max, zp_max - zp_min
       WRITE(*,*) 'x fluid min/max', domainRange(:,1)
       WRITE(*,*) 'y fluid min/max', domainRange(:,2)
       WRITE(*,*) 'z fluid min/max', domainRange(:,3)
     END IF
     CALL MPI_Barrier(global%mpiComm,errorFlag)
     !CALL ErrorStop(global,0,459,"rocpicl Init: Particles outside fluid domain")
   END IF

   ! Close points.dat file
   CLOSE(iFile, IOSTAT=errorFlag)
   global%error = errorFlag   
   IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
   END IF ! global%error  
ELSE
!  This is for a restart
   fexists = .true.
   ii=0 !last exists par.vtu file
   DO WHILE (fexists)
      ii = ii + 1
      vtuFile = ''
      WRITE(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
      INQUIRE(FILE=trim(vtuFile), EXIST=fexists)
      IF( global%myProcid == MASTERPROC) THEN
         WRITE(*,*) 'PPICLF par file: ',TRIM(vtuFile),ii,'  ',fexists
      END IF
   END DO
   ii = ii - 1
   vtuFile = ''
   WRITE(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
   IF( global%myProcid == MASTERPROC) THEN
      WRITE(*,*) 'Reading ', vtuFile, len(vtuFile)
   ENDIF

   IF(ii .lt. 0) THEN
      CALL ErrorStop(global,ERR_FILE_EXIST,__LINE__,vtuFile)
   END IF

   ! The particle files hold only the fields chosen for output
   ! (ppiclf_user_OutputFields), not the whole state, so they cannot be
   ! read back. A separate particle restart file is planned.
   CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
     'Particle restart is not available yet: par*.vtu files hold only '// &
     'the output fields. Remove the .rin file to start a new run.')

   npart = -1
   dp_max = -1.0
   CALL ppiclf_io_ReadParticleVTU(trim(vtuFile), ii, npart, dp_max)
   CALL MPI_Allreduce(dp_max,dp_max,1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
   print*,global%myProcid,npart,dp_max,dp_max
   IF( global%myProcid == MASTERPROC) THEN
      print*, " "
      WRITE(*,*) 'Finished PPICLF Restart'
      print*, " "
      print*, " ======================================="
      print*, " "
   END IF
END IF ! global%restartFromScratch

!CALL MPI_Barrier(global%mpiComm,errorFlag)


! User sets up overlap grid:
nCells = pRegion%grid%nCells
 
ALLOCATE(xGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
END IF ! global%error

ALLOCATE(yGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:yGrid')
END IF ! global%error

ALLOCATE(zGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:zGrid')
END IF ! global%error

ALLOCATE(rocGrid(7,nCells),STAT=errorFlag)
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:rocGrid')
END IF ! global%error



DO i = 1, nCells 
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Tetrahedral Cell
    CellVertices = 4
    DO k = 1, CellVertices
      vi = pRegion%grid%tet2v(k,i) 
             xGrid(k,i) = pRegion%grid%xyz(XCOORD,vi) 
             yGrid(k,i) = pRegion%grid%xyz(YCOORD,vi) 
             zGrid(k,i) = pRegion%grid%xyz(ZCOORD,vi) 
    END DO
  ELSE IF(pGrid%cellGlob2Loc(1,i) == 2) THEN
    ! Hexahedral Cell
    CellVertices = 8
    DO k = 1, CellVertices
      vi = pRegion%grid%hex2v(k,i) 
             xGrid(k,i) = pRegion%grid%xyz(XCOORD,vi) 
             yGrid(k,i) = pRegion%grid%xyz(YCOORD,vi) 
             zGrid(k,i) = pRegion%grid%xyz(ZCOORD,vi) 
    END DO !nCells
  ELSE
    CellVertices = 0
    WRITE(*,*) 'ERROR: Rocflupicl only support tetrahedral and hexahedral cell types.'
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:CellLen')
  END IF
END DO
! Find cell lengths
Max_CellLen  = 0.0D0
DO i = 1,nCells
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Tetrahedral Cell
    CellVertices = 4
  ELSE IF(pGrid%cellGlob2Loc(1,i) == 2) THEN
    ! Hexahedral Cell
    CellVertices = 8
  ELSE
    CellVertices = 0
    WRITE(*,*) 'ERROR: Rocflupicl only support tetrahedral and hexahedral cell types.'
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:CellLen')
  END IF
  ! Initialize as zero for each cell
  DO l = 1,3
    MaxPoint(l) = -1.0D10 
    MinPoint(l) =  1.0D10 
    CellLen(l)   =  0.0D0   
  END DO !l
  ! Add all x,y,z cell corners for centroid and find extremes
  DO k = 1,CellVertices
    IF (xGrid(k,i) > MaxPoint(1)) &
      MaxPoint(1) = xGrid(k,i)
    IF (xGrid(k,i) < MinPoint(1)) &
      MinPoint(1) = xGrid(k,i)
    IF (yGrid(k,i) > MaxPoint(2)) &
      MaxPoint(2) = yGrid(k,i)  
    IF (yGrid(k,i) < MinPoint(2)) &
      MinPoint(2) = yGrid(k,i)
    IF (zGrid(k,i) > MaxPoint(3)) &
      MaxPoint(3) = zGrid(k,i)  
    IF (zGrid(k,i) < MinPoint(3)) &
      MinPoint(3) = zGrid(k,i)
  END DO !k
  DO l = 1,3
    ! Find element length in all dimensions
    CellLen(l) = ABS(MaxPoint(l)-MinPoint(l))
    ! Find max lengths for grid cells in particle domain on this processor.
    ! Use cell bounding-box overlap with the particle box, not centroid
    ! inclusion: a narrow cloud (e.g. one particle) may contain no centroid,
    ! which gave Max_CellLen = 0 and a zero filter width.
    IF( (MaxPoint(1) >= xp_min) .AND. (MinPoint(1) <= xp_max) .AND. &
        (MaxPoint(2) >= yp_min) .AND. (MinPoint(2) <= yp_max) .AND. &
        (MaxPoint(3) >= zp_min) .AND. (MinPoint(3) <= zp_max)) THEN
      IF(CellLen(l) .GT. Max_CellLen(l)) Max_CellLen(l) = CellLen(l)
      ! For tets, make a rectangular box around it with maxCellLength equal in all dimensions
      IF(pGrid%cellGlob2Loc(1,i) == 1) Max_CellLen(l) = MAXVAL(Max_CellLen(:))
    END IF
  END DO !l
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Make cell length same in all dimensions for tetrahedral
    ! This is due to their shape not being rectangular
    DO l = 1,3
      CellLen(l) = MAX(CellLen(1),CellLen(2),CellLen(3))
    END DO !l
  END IF

  rocGrid(1,i) = pRegion%grid%cofg(1,i) !Cell Centroid x position
  rocGrid(2,i) = pRegion%grid%cofg(2,i) !Cell Centroid y position
  rocGrid(3,i) = pRegion%grid%cofg(3,i) !Cell Centroid z position
  rocGrid(4,i) = CellLen(1)             !Cell largest dx
  rocGrid(5,i) = CellLen(2)             !Cell largest dy
  rocGrid(6,i) = CellLen(3)             !Cell largest dz
  rocGrid(7,i) = pRegion%grid%vol(i)    !Cell volume
END DO !i

!     CALL MPI_Barrier(global%mpiComm,errorFlag)

MinFluidCells = 2.0D0 !Number of fluid cells for Minimum Bin Size and ppiclf_filter(1:3)
DO l = 1,3
  filter_local(l) = MinFluidCells*Max_CellLen(l) 
  ! Find max x,y,z cell lengths across all MPI ranks (entire rocflu domain)
  CALL MPI_Allreduce(filter_local(l),filter(l),1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
END DO

neighborWidth = 4.0_RFREAL*dp_max ! Minimum value. User can entire larger if desired
IF((neighborWidth .GT. global%piclNeighborWidth) &
    .AND. (global%myProcid == MASTERPROC)) THEN
    WRITE(STDOUT, '(A)') &
        '*** WARNING *** PICL NEIGHBORWIDTH too small, defaulting to 4*dp_max'
END IF
neighborWidth = MAX(neighborWidth, global%piclNeighborWidth)

!CALL MPI_Barrier(global%mpiComm,errorFlag)

IF(global%myProcid == MASTERPROC) THEN
   PRINT*,' '
   PRINT*,'PPICLF: '
   PRINT*,'  Inputs NEIGHBORWIDTH  : ',global%piclNeighborWidth
   PRINT*,'  dp_max (points.dat)   : ',dp_max
   PRINT*,'  ppiclf_filter(1)      : ',filter(1)
   PRINT*,'  ppiclf_filter(2)      : ',filter(2)
   PRINT*,'  ppiclf_filter(3)      : ',filter(3)
   PRINT*,'  ppiclf_nndist         : ',neighborWidth
   PRINT*,' '
END IF

! Sets ppiclf_nndist, ppiclf_filter(1:3), ppiclf_y(), ppiclf_rprop()
!(imethod (RK pick), nDimensions, iendian (IO format), .....)
! On a restart the particles come from the par file and y, rprop are not
! read; ppiclF still takes the arguments
IF (.NOT. ALLOCATED(y)) THEN
  ALLOCATE(y(nStateSlots,1),rprop(nPropSlots,1))
  y = 0.0_RFREAL
  rprop = 0.0_RFREAL
END IF
CALL ppiclf_solve_InitParticle(2,3,0,npart_local,y,rprop,filter,neighborWidth) 
IF(global%myProcid == MASTERPROC) THEN
  PRINT*, 'x fluid min/max', domainRange(:,1)
  PRINT*, 'y fluid min/max', domainRange(:,2)
  PRINT*, 'z fluid min/max', domainRange(:,3)
END IF

! Creates OverlapGrid (fluid cells, domain extent, periodicity) and Calls
! Init Solve
CALL ppiclf_comm_InitOverlapGrid(nCells,rocGrid,domainRange,periodic)
!CALL MPI_Barrier(global%mpiComm,errorFlag)
CALL ppiclf_solve_InitSolve

! Particle walls used to be read from filein.vtk. They now come only from
! the fluid wall patches, so stop rather than silently run without walls.
INQUIRE(FILE='filein.vtk', EXIST=wall_exists)
IF ( wall_exists ) THEN
  CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
    'filein.vtk is no longer read: set PARTWALL 1 on the wall patches '// &
    'in the .bc file and remove filein.vtk')
END IF

! Particle walls from the Rocflu wall patches marked PARTWALL 1 in the .bc
! file. Every rank gets every face (ppiclF bins do not follow the Rocflu
! regions); quadrilaterals are split in two. The normal -fn points into the
! fluid, so these walls push particles back even after their centre has
! crossed the wall (see libpicl/source/ppiclf_wall.f).
nTriLoc = 0
DO iPatch = 1,pGrid%nPatches
  pPatch => pRegion%patches(iPatch)
  IF ( .NOT. pPatch%partWallFlag ) CYCLE
  DO ifc = 1,pPatch%nBFaces
    IF ( pPatch%bf2v(4,ifc) == VERT_NONE ) THEN
      nTriLoc = nTriLoc + 1
    ELSE
      nTriLoc = nTriLoc + 2
    END IF
  END DO
END DO

ALLOCATE(triLoc(12,MAX(nTriLoc,1)))
iTri = 0
DO iPatch = 1,pGrid%nPatches
  pPatch => pRegion%patches(iPatch)
  IF ( .NOT. pPatch%partWallFlag ) CYCLE
  DO ifc = 1,pPatch%nBFaces
    nv = 4
    IF ( pPatch%bf2v(4,ifc) == VERT_NONE ) nv = 3
    DO iv = 1,nv
      xv(1:3,iv) = pGrid%xyz(XCOORD:ZCOORD,pPatch%bv(pPatch%bf2v(iv,ifc)))
    END DO
    iTri = iTri + 1
    triLoc(1:3,iTri)   = xv(1:3,1)
    triLoc(4:6,iTri)   = xv(1:3,2)
    triLoc(7:9,iTri)   = xv(1:3,3)
    triLoc(10:12,iTri) = -pPatch%fn(XCOORD:ZCOORD,ifc)
    IF ( nv == 4 ) THEN
      iTri = iTri + 1
      triLoc(1:3,iTri)   = xv(1:3,1)
      triLoc(4:6,iTri)   = xv(1:3,3)
      triLoc(7:9,iTri)   = xv(1:3,4)
      triLoc(10:12,iTri) = -pPatch%fn(XCOORD:ZCOORD,ifc)
    END IF
  END DO
END DO

CALL MPI_Comm_size(global%mpiComm,nProcs,errorFlag)
ALLOCATE(triCounts(nProcs),triDispls(nProcs))
CALL MPI_Allgather(12*nTriLoc,1,MPI_INTEGER,triCounts,1,MPI_INTEGER, &
                   global%mpiComm,errorFlag)
triDispls(1) = 0
DO k = 2,nProcs
  triDispls(k) = triDispls(k-1) + triCounts(k-1)
END DO
nTriTot = (triDispls(nProcs) + triCounts(nProcs))/12
ALLOCATE(triAll(12,MAX(nTriTot,1)))
CALL MPI_Allgatherv(triLoc,12*nTriLoc,MPI_RFREAL,triAll,triCounts, &
                    triDispls,MPI_RFREAL,global%mpiComm,errorFlag)
DO k = 1,nTriTot
  CALL ppiclf_solve_AddWall(triAll(1:3,k),triAll(4:6,k),triAll(7:9,k), &
                            triAll(10:12,k))
END DO

! Report which global patches are particle walls
ALLOCATE(partWallPatch(global%nPatches))
partWallPatch = 0
DO iPatch = 1,pGrid%nPatches
  pPatch => pRegion%patches(iPatch)
  IF ( pPatch%partWallFlag ) partWallPatch(pPatch%iPatchGlobal) = 1
END DO
CALL MPI_Allreduce(MPI_IN_PLACE,partWallPatch,global%nPatches,MPI_INTEGER, &
                   MPI_MAX,global%mpiComm,errorFlag)
IF ( global%myProcid == MASTERPROC ) THEN
  IF ( nTriTot > 0 ) THEN
    WRITE(*,'(A,I8,A)') ' Particle walls: ',nTriTot, &
      ' triangles from the patches with PARTWALL 1'
    WRITE(*,'(A)',ADVANCE='NO') '   patches:'
    DO k = 1,global%nPatches
      IF ( partWallPatch(k) == 1 ) WRITE(*,'(I5)',ADVANCE='NO') k
    END DO
    WRITE(*,*)
    IF ( global%piclCollisionFlag == 0 ) THEN
      WRITE(*,*) '*** WARNING *** particle walls have no effect with COLLISIONS 0'
    END IF
  ELSE
    WRITE(*,*) 'Particle walls: none (no wall patch has PARTWALL 1)'
  END IF
END IF
DEALLOCATE(triLoc,triAll,triCounts,triDispls,partWallPatch)

! 03/24/2025 - Thierry - store the RocfluMP Flow Model chosen (Euler or NS)
!                        this is used in ppiclF for calculating the pressure gradient
!                        whether with or without the viscous part
!flow_model = 0
!flow_model = int(pRegion%mixtInput%flowModel)

DEALLOCATE(rocGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:rocGrid')
END IF ! global%error


DEALLOCATE(xGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
END IF ! global%error

DEALLOCATE(yGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:yGrid')
END IF ! global%error

DEALLOCATE(zGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:zGrid')
END IF ! global%error

ALLOCATE(volp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:volp')
    END IF ! global%error

DO i=1,pGrid%nCellsTot
        pRegion%mixt%piclVF(i) = 0.0_RFREAL
END DO

IF ( global%myProcid == MASTERPROC) WRITE(*,*) "PFINIT: Calc Init VolP"
DO i = 1, nCells
       CALL ppiclf_solve_GetProFld(i, PPICLF_P_JPHIP, volp(i))
       volp(i) = volp(i)/pRegion%grid%vol(i)
       pRegion%mixt%piclVF(i) = volp(i) 
END DO

maxVF = MAXVAL(volp)
CALL MPI_Allreduce(maxVF,maxVF,1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
! Particle Volume Fraction Sanity Check
IF(maxVF .GT. 0.65 .AND. global%myProcid == MASTERPROC) THEN
  PRINT*, '*** WARNING*** Max particle vf: ', maxVF, & 
         'Adjust grid for a larger minimum cell size or decrease particle diameter'
END IF

! TLJ:
! This section takes as input from utilities/init/RFLU_InitFlowHardCode.F90
!    (r,r*u,r*v,r*w,r*E) and changes to RocfluMP conserved variables
!    (phig*r,phig*r*u,phig*r*v,phig*r*w,phig*r*E), where phig is the gas
!    phase volume fraction and can only be computed after the particles
!    are read in.
! Only actual cells: piclVF of virtual cells is not known yet (it is
! exchanged afterwards in RFLU_InitFlowSolver, together with cv). Converting
! virtual cells with piclVF = 0 left them at the gas-phase state, so the
! two sides of a region border disagreed during the first RK stage.

DO icg = 1,pGrid%nCells
    vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(icg)
    pRegion%mixt%cv(CV_MIXT_DENS,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_DENS,icg)
    pRegion%mixt%cv(CV_MIXT_XMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
    pRegion%mixt%cv(CV_MIXT_YMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ENER,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ENER,icg)
    IF(pRegion%mixt%cv(CV_MIXT_DENS,icg) .le. 0.0) THEN
         WRITE(*,*) "Error: negative density: ",pRegion%mixt%cv(CV_MIXT_DENS,icg)     
         PRINT*, 'From rocpicl/PICL_InitSolver.F90' 
         CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,'PPICLF:init')
    END IF    
END DO ! icg

DEALLOCATE(volp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:zGrid')
    END IF ! global%error

#endif

IF ( global%myProcid == MASTERPROC) then
   print*, ' '
   print*, '***********************************************'
   print*, 'Starting PICL_InitSolver.F90'
   print*, ' '
   print*, 'Velocity Model              = ',global%piclVelocityModel
   print*, 'QS Flag                     = ',global%piclQsFlag
   print*, 'QS Fluct Flag               = ',global%piclQsFluctFlag
   print*, 'Pseudo Turbulence Flag      = ',pseudoTurb_flag
   print*, 'Added Mass Flag             = ',global%piclAmFlag
   print*, 'Pressure Gradient Flag      = ',global%piclPgFlag
   print*, 'Collisions Flag             = ',global%piclCollisionFlag
   print*, 'Viscous-Unsteady Flag       = ',global%piclViscousUnsteady
   print*, 'Average Heat Transfer Flag  = ',global%piclHeatTransferFlag
   print*, 'Unsteady Heat Transfer Flag = ',global%piclHTUnsteadyFlag
   print*, 'Feedback Flag               = ',global%piclFeedbackFlag
   print*, 'Unsteady history steps      = ',ppiclf_nUnsteadyData
   print*, 'SubBin Nearest Flag         = ',global%piclSBNearFlag

   IF (global%piclViscousUnsteady >=1) THEN
      print*,'  Using Viscous unsteady history term'
      print*,'    ppiclf_nTimeBH       = ',ppiclf_nTimeBH
      print*,'    ppiclf_nUnsteadyData = ',ppiclf_nUnsteadyData
   ENDIF

   print*, ' '
   print*, "XLOC = ", y(PPICLF_JX ,1) 
   print*, "YLOC = ", y(PPICLF_JY ,1) 
   print*, "ZLOC = ", y(PPICLF_JZ ,1)   

   print*, ' '
   print*, 'Reading points.dat file...'
   print*, "npart          = ", npart
   print*, 'dp_max         = ', dp_max
   print*, 'rho_dens       = ', rhop
   print*, 'cv_particle    = ', ppiclf_rcp_part
   
   if(periodic(1)) print*, "ppiclF X-periodic (min, max) ", domainRange(:,1)
   if(periodic(2)) print*, "ppiclF Y-periodic (min, max) ", domainRange(:,2)
   if(periodic(3)) print*, "ppiclF Z-periodic (min, max) ", domainRange(:,3)

   print*, ' '
   print*, 'Ending PICL_InitSolver.F90'
   print*, '***********************************************'
ENDIF

! finalize --------------------------------------------------------------------

  ! y(:,1) is printed in the summary above; ppiclF holds its own copy
  IF (ALLOCATED(y)) DEALLOCATE(y,rprop)

  CALL DeregisterFunction( global )

END SUBROUTINE PICL_InitSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************


