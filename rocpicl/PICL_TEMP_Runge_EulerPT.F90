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

SUBROUTINE PICL_TEMP_Runge(pRegion)

!  USE 

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt

  USE RFLU_ModDifferentiationCells
  USE RFLU_ModLimiters, ONLY: RFLU_CreateLimiter, &
                              RFLU_ComputeLimiterBarthJesp, &
                              RFLU_ComputeLimiterVenkat, &
                              RFLU_LimitGradCells, &
                              RFLU_LimitGradCellsSimple, &
                              RFLU_DestroyLimiter
  USE RFLU_ModWENO, ONLY: RFLU_WENOGradCellsWrapper, &
                          RFLU_WENOGradCellsXYZWrapper
#ifdef PICL
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                             RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 
#endif



#ifdef PICL
!DEC$ NOFREEFORM
#include "../libpicl/ppiclF/source/PPICLF_USER.h"
#include "../libpicl/ppiclF/source/PPICLF_STD.h"
!DEC$ FREEFORM
#endif


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

#ifdef PICL
  LOGICAL :: doWrite      
  INTEGER(KIND=4) :: i,piclIO,nCells
  INTEGER :: errorFlag,icg      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime, &
                  temp_dudtMixt, temp_dvdtMixt, temp_dwdtMixt, energydotg, &
                  dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, vFrac, Cp

  REAL(KIND=8), DIMENSION(3) :: ug      
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhoF, uxF, uyF, uzF, csF, tpF, vfP, &
                                             dpxF, dpyF, dpzF, SDRX, SDRY, SDRZ, rhsR, &
                                             pGcX, pGcY, pGcZ, PhiP, YTEMP, &
                                             JFXCell, JFYCell, JFZCell, JFECell, &
                                             domgdx, domgdy, domgdz, drhodx, drhody, drhodz, &
                                             dpvxF, dpvyF, dpvzF, SDOX, SDOY, SDOZ

  REAL(KIND=8), DIMENSION(:,:,:), POINTER :: pGc 
!---------------------------------------------------------------  
! Below added for pseudo turbulence
  INTEGER(KIND=4) :: j
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhog, DivPhiQsg, ParticleDiameter, & 
                                             RelVelMagCell, ReynoldsNoCell, MachNoCell, &
                                             viscosityGas
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JRSGCell, JTSGCell, JAlphaPTCell, JPTHFCell, &
                                               DivPhiRSG, ugas, Qsg, JVpUpCell, urelCell, upCell, &
                                               reynoldsStress
  
  INTEGER(KIND=4) :: dir
  REAL(KIND=8) :: bPar, bPerp, kTilde, MeanKE, rhoK
  REAL(KIND=8), DIMENSION(3,3) :: Q, temp, bij
  REAL(KIND=8), DIMENSION(3) :: u1, u2, u3, u1dot, u2dot, u3dot                                             
  REAL(KIND=8), PARAMETER :: oneThird = 1.0_RFREAL / 3.0_RFREAL                                            
  REAL(KIND=8), PARAMETER :: a = 0.523_RFREAL
  REAL(KIND=8), PARAMETER :: b = 0.305_RFREAL
  REAL(KIND=8), PARAMETER :: c = 0.114_RFREAL
  REAL(KIND=8), PARAMETER :: d = 3.511_RFREAL
  REAL(KIND=8), PARAMETER :: e = 1.801_RFREAL
  REAL(KIND=8), PARAMETER :: f = 0.005_RFREAL
  !---------------------------------------------------------------  
  
  ! TLJ - added for Feedback term - 04/01/2025
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfoPicl
  INTEGER, DIMENSION(:), POINTER :: piclcvInfo
    REAL(KIND=8) :: dodx, dody, dodz,     &
                  omgx, omgy, omgz,     &
                  divu,                 &
                  dprdx, dprdy, dprdz,  &
                  dpdx, dpdy, dpdz,     &
                  phirho, ir, ir2 ,     &
                  dfxdx, dfxdy, dfxdz,  &
                  dfydx, dfydy, dfydz,  &
                  dfzdx, dfzdy, dfzdz   

#endif
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_Runge.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction(global, 'PICL_TEMP_Runge',__FILE__ )

! Set pointers ----------------------------------------------------------------

    !pRegion => regions!pLevel%regions(iReg)
    pGrid   => pRegion%grid

!PPICLF Integration
#ifdef PICL

     piclIO = 100000000
     piclDtMin = REAL(global%dtMin,8)
     piclCurrentTime = REAL(global%currentTime,8)

     ! TLJ - 11/23/2024
     !     - This has now been removed
     doWrite = RFLU_DecideWrite(global)
     !Figure out piclIO call, might need to look into timestepping
     IF((doWrite .EQV. .TRUE.)) piclIO = 1


!PARTICLE stuff possbile needed
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)


!allocate arrays to send to picl
    nCells = pRegion%grid%nCells
    ALLOCATE(rhoF(nCells), uxF(nCells), uyF(nCells), uzF(nCells), csF(nCells), &
             tpF(nCells), vfP(nCells), dpxF(nCells), dpyF(nCells), dpzF(nCells), &
             SDRX(nCells), SDRY(nCells), SDRZ(nCells), rhsR(nCells), pGcX(nCells), pGcY(nCells), &
             pGcZ(nCells), JFXCell(nCells), JFYCell(nCells), JFZCell(nCells), JFECell(nCells), &
             PhiP(nCells), domgdx(nCells), domgdy(nCells), domgdz(nCells), &
             drhodx(nCells), drhody(nCells), drhodz(nCells), dpvxF(nCells), dpvyF(nCells), &
             dpvzF(nCells), SDOX(nCells), SDOY(nCells), SDOZ(nCells), STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    IF(pRegion%mixtInput%axiFlag) THEN
      ALLOCATE(YTEMP(nCells),STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

! Added for pseudo turbulence
!---------------------------------------------------------------  
! Maybe add if statement to only allocate if PseudoTurbulence is used
    ALLOCATE(JRSGCell(9,nCells), JTSGCell(3,nCells), JAlphaPTCell(9,nCells), &
             JPTHFCell(3,nCells), DivPhiRsg(3,nCells), rhog(nCells), ugas(3,nCells), &
             Qsg(3,nCells), DivPhiQsg(nCells), JVpUpCell(3,nCells), urelCell(3,nCells), &
             upCell(3,nCells), ReynoldsNoCell(nCells), ParticleDiameter(nCells), & 
             RelVelMagCell(nCells), MachNoCell(nCells), viscosityGas(nCells), & 
             reynoldsStress(9,nCells) &
            ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
!---------------------------------------------------------------  

pGc => pRegion%mixt%gradCell
    ! 04/01/2025 - TLJ - we need feedback terms and their gradients to
    !       calculate the undisturbed torque component
    ! Internal definitions; some redundancy but just ignore
    ! We do not need energy, but might in the future
    DO i = 1,pRegion%grid%nCells
       JFXCell(i) = 0.0_RFREAL
       JFYCell(i) = 0.0_RFREAL
       JFZCell(i) = 0.0_RFREAL
       JFECell(i) = 0.0_RFREAL
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFX,JFXCell(i))  
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFY,JFYCell(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFZ,JFZCell(i))
       pregion%mixt%piclFeedback(1,i) = JFXCell(i)
       pregion%mixt%piclFeedback(2,i) = JFYCell(i)
       pregion%mixt%piclFeedback(3,i) = JFZCell(i)
    END DO
    ! Now calculate the gradient of the feedback force
    ALLOCATE(varInfoPicl(3), piclcvInfo(3), STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    varInfoPicl(1) = 1
    varInfoPicl(2) = 2
    varInfoPicl(3) = 3
    piclcvInfo = varInfoPicl
    CALL RFLU_ComputeGradCellsWrapper(pRegion,1,3,1,3,varInfoPicl, &
                                      pRegion%mixt%piclFeedback,&
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,3, &
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_LimitGradCellsSimple(pRegion,1,3,1,3, &
                                   pRegion%mixt%piclFeedback,&
                                   piclcvInfo,&
                                   pRegion%mixt%piclgradFeedback)
    DEALLOCATE(varInfoPicl, piclcvInfo, STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

  ! END - TLJ calculating gradient of feedback force

!Fill arrays for interp field
    DO i = 1,pRegion%grid%nCells
!Zero out phip
       PhiP(i) = 0.0_RFREAL

       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
!------------------------------------------------------------------------------------- 
! Brad's Old formulation. Keep it here for reference
!         temp_drudtMixt = -pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,2,i)) & 
!                   +ug(XCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!                                                                                    
!         temp_drvdtMixt = -pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,3,i))&  
!                   +ug(YCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!                                                                                    
!         temp_drwdtMixt = -pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,4,i))&  
!                   +ug(ZCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!------------------------------------------------------------------------------------- 
! IMPORTANT NOTE:
! 05/25/2025 - Thierry - rhs & diss arrays have to be divided by pRegion%grid%vol(i)
!                        when being used for interpolating values from Rocflu to ppiclF
!------------------------------------------------------------------------------------- 
      
       ! 03/11/2025 - Thierry - Du/Dt, Dv/Dt, Dw/Dt (not weighted by phi^g or rho^g)

       temp_dudtMixt  = (-pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)& 
                         +ug(XCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,2,i))

       temp_dvdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)& 
                         +ug(YCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,3,i))
                         
       temp_dwdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)& 
                         +ug(ZCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,4,i))

       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       vfP(i) = vfP(i)/pRegion%grid%vol(i)
       PhiP(i) = vfP(i)
       !VOL Frac cap
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       vfp(i) = PhiP(i)      

       ! TLJ - 02/07/2025 scaled conserved density by gas-phase volume fraction
       vFrac = 1.0_RFREAL - PhiP(i)!pRegion%mixt%piclVF(i)

       rhoF(i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / vFrac
       uxF(i) =  pRegion%mixt%cv(CV_MIXT_XMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uyF(i) =  pRegion%mixt%cv(CV_MIXT_YMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uzF(i) =  pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)

       csF(i) = pRegion%mixt%dv(DV_MIXT_SOUN,i)
       tpF(i) = pRegion%mixt%dv(DV_MIXT_TEMP,i) 
       ! Davin - added pressure to interpolation values 02/22/2025

       dpxF(i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpyF(i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpzF(i) = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz

       dudx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_XVEL,i)
       dudy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_XVEL,i)
       dudz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_XVEL,i)

       dvdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_YVEL,i)
       dvdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_YVEL,i)
       dvdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_YVEL,i)

       dwdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_ZVEL,i)
       dwdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_ZVEL,i)
       dwdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_ZVEL,i)

       domgdx(i) = dwdy - dvdz
       domgdy(i) = dudz - dwdx
       domgdz(i) = dvdx - dudy

       ! 04/01/2025 - TLJ - Calculate the substantial derivative of vorticity
       ! Internal definitions; some redundancy but just ignore
       dodx   = 0.0_RFREAL ! D(Omega_x)/DT
       dody   = 0.0_RFREAL ! D(Omega_y)/DT
       dodz   = 0.0_RFREAL ! D(Omega_z)/DT
       omgx   = dwdy - dvdz ! Omega_x
       omgy   = dudz - dwdx ! Omega_y
       omgz   = dvdx - dudy ! Omega_z
       divu   = dudx + dvdy + dwdz ! u_x+v_y+w_z; divergence of velocity
       dprdx  = pGc(XCOORD,1,i) ! d(rho phi)/dx
       dprdy  = pGc(YCOORD,1,i) ! d(rho phi)/dy
       dprdz  = pGc(ZCOORD,1,i) ! d(rho phi)/dz
       dpdx   = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpdy   = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpdz   = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz
       dfxdx  = pRegion%mixt%piclgradFeedback(XCOORD,1,i) ! dFx/dx
       dfxdy  = pRegion%mixt%piclgradFeedback(YCOORD,1,i) ! dFx/dy
       dfxdz  = pRegion%mixt%piclgradFeedback(ZCOORD,1,i) ! dFx/dz
       dfydx  = pRegion%mixt%piclgradFeedback(XCOORD,2,i) ! dFy/dx
       dfydy  = pRegion%mixt%piclgradFeedback(YCOORD,2,i) ! dFy/dy
       dfydz  = pRegion%mixt%piclgradFeedback(ZCOORD,2,i) ! dFy/dz
       dfzdx  = pRegion%mixt%piclgradFeedback(XCOORD,3,i) ! dFz/dx
       dfzdy  = pRegion%mixt%piclgradFeedback(YCOORD,3,i) ! dFz/dy
       dfzdz  = pRegion%mixt%piclgradFeedback(ZCOORD,3,i) ! dFz/dz
       phirho = pRegion%mixt%cv(CV_MIXT_DENS,i) ! phi_g*rho_g
       ir     = 1.0_RFREAL / phirho
       ir2    = ir*ir
       ! 1. Vortex stretching
       dodx = omgx*dudx + omgy*dudy + omgz*dudz
       dody = omgx*dvdx + omgy*dvdy + omgz*dvdz
       dodz = omgx*dwdx + omgy*dwdy + omgz*dwdz
       ! 2. Vortex dilatation
       dodx = dodx - omgx*divu
       dody = dody - omgy*divu
       dodz = dodz - omgz*divu
       ! 3. Baroclinic
       dodx = dodx + (dprdy*dpdz - dprdz*dpdy)*ir2
       dody = dody + (dprdz*dpdx - dprdx*dpdz)*ir2
       dodz = dodz + (dprdx*dpdy - dprdy*dpdx)*ir2
       ! 4. Torque due to feedback force
       dodx = dodx + (dfzdy - dfydz)*ir
       dody = dody + (dfxdz - dfzdx)*ir
       dodz = dodz + (dfydx - dfxdy)*ir
       ! 5. Misalignment of phi*rho and feedback force
       dodx = dodx + (dprdy*JFZCell(i) - dprdz*JFYCell(i))*ir2
       dody = dody + (dprdz*JFXCell(i) - dprdx*JFZCell(i))*ir2
       dodz = dodz + (dprdx*JFYCell(i) - dprdy*JFXCell(i))*ir2
       ! 6. Add terms and store
       SDOX(i) = dodx
       SDOY(i) = dody
       SDOZ(i) = dodz
       ! End - TLJ - Calculate the substantial derivative of vorticity

       ! Substantial derivative of gas-phase velocity
       SDRX(i) = temp_dudtMixt ! Du/Dt
       SDRY(i) = temp_dvdtMixt ! Dv/Dt
       SDRZ(i) = temp_dwdtMixt ! Dw/Dt

       rhsR(i) = -pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i) ! \p(rho*phi)/\p(t)

       pGcX(i) = pGc(XCOORD,1,i) ! d(rho phi)/dx
       pGcY(i) = pGc(YCOORD,1,i) ! d(rho phi)/dy
       pGcz(i) = pGc(ZCOORD,1,i) ! d(rho phi)/dz

       ! Gradient of rho^g of mixture (not weighted by phi^g!)
       ! Using grad(rhog) directly
       drhodx(i) = pRegion%mixt%piclgradRhog(1,1,i) ! d(rho)/dx
       drhody(i) = pRegion%mixt%piclgradRhog(2,1,i) ! d(rho)/dy
       drhodz(i) = pRegion%mixt%piclgradRhog(3,1,i) ! d(rho)/dz

       ! Viscous term of pressure gradient (divergence of tau)
       dpvxF(i) = pRegion%mixt%diss(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)
       dpvyF(i) = pRegion%mixt%diss(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)
       dpvzF(i) = pRegion%mixt%diss(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)
     END DO

! Interp field calls
! TLJ - interpolates various fluid quantities onto the 
!       the ppiclf particle locations
! Each slot below must be one of the interpolated slots that
! ppiclf_user_Layout lists (AddInterpSlot), given once per stage;
! ppiclF stops if a listed slot is missing or an unlisted one is given.
 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOF,rhoF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,uxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,uyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,uzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDX,dpxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDY,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDZ,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JCS,csF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRX,SDRX)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRY,SDRY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRZ,SDRZ)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHSR,rhsR)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCX,pGcX) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCY,pGcY) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCZ,pGcZ) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JXVOR,domgdx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JYVOR,domgdy)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JZVOR,domgdz)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGX,drhodx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGY,drhody)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGZ,drhodz)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDX,dpvxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDY,dpvyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDZ,dpvzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOX,SDOX)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOY,SDOY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOZ,SDOZ)  

 ! Solve RK stage of time stepping particle solution
     CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)

!FEED BACK TERMS
     JFXCell = 0.0_RFREAL
     JFYCell = 0.0_RFREAL
     JFZCell = 0.0_RFREAL
     JFECell = 0.0_RFREAL
     JRSGCell = 0.0_RFREAL
     JAlphaPTCell = 0.0_RFREAL
     JTSGCell = 0.0_RFREAL
     Qsg = 0.0_RFREAL
     DivPhiRSG = 0.0_RFREAL
     DivPhiQsg = 0.0_RFREAL
     JVpUpCell = 0.0_RFREAL
     upCell = 0.0_RFREAL


!Fill arrays for interp field
IF(global%piclFeedbackFlag == 1) THEN
    DO i = 1,pRegion%grid%nCells
       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

       
      ! Just testing to make sure i can get rid of the duplicate variables
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFX,JFXCell(i))  
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFY,JFYCell(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFZ,JFZCell(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JE,JFECell(i)) 
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       PhiP(i) = vfP(i)/pRegion%grid%vol(i)

       !VOL Frac cap
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       vfp(i) = PhiP(i)      
   !---------------------------------------------------------------------------------------
       ! 07/21/2025 - Thierry - begins here - added for PseudoTurbulence
       if(global%piclPseudoTurbFlag .gt. 0) then
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG11,JRSGCell(1,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG12,JRSGCell(2,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG13,JRSGCell(3,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG21,JRSGCell(4,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG22,JRSGCell(5,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG23,JRSGCell(6,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG31,JRSGCell(7,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG32,JRSGCell(8,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG33,JRSGCell(9,i))

         !JRSGCell(1,i) = JRSGCell(1,i) / (vfP(i)*pRegion%grid%vol(i) + epsilon(1.0_RFREAL))

         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG1,JTSGCell(1,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG2,JTSGCell(2,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG3,JTSGCell(3,i))

         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT11,JAlphaPTCell(1,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT12,JAlphaPTCell(2,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT13,JAlphaPTCell(3,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT21,JAlphaPTCell(4,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT22,JAlphaPTCell(5,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT23,JAlphaPTCell(6,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT31,JAlphaPTCell(7,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT32,JAlphaPTCell(8,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JAlphaPT33,JAlphaPTCell(9,i))

         ! Collect back Particle Volume*Velocity: Vp*Up
         call ppiclf_solve_GetProFld(i,PPICLF_P_JVPVX,JVpUpCell(1,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JVPVY,JVpUpCell(2,i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JVPVZ,JVpUpCell(3,i))

         ! Vp*Up -> Up (velocity of particles per cell)

         upCell(1,i) = JVpUpCell(1,i) / (vfP(i)*pRegion%grid%vol(i) + epsilon(1.0_RFREAL))
         upCell(2,i) = JVpUpCell(2,i) / (vfP(i)*pRegion%grid%vol(i) + epsilon(1.0_RFREAL))
         upCell(3,i) = JVpUpCell(3,i) / (vfP(i)*pRegion%grid%vol(i) + epsilon(1.0_RFREAL))


       endif ! piclPseudoTurbFlag
!---------------------------------------------------------------------------------------
       ! I need to look more into this, don't really understand it
       energydotg = JFECell(i) ! includs KE feedback already

       IF(IsNan(JFXCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PX",i,JFXCell(i),ug(1),ug(2),ug(3)
         write(*,*) "JFY",i,JFYCell(i)
         write(*,*) "JFZ",i,JFZCell(i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PX')
       END IF
       IF(IsNan(JFYCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PY",i,JFYCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
       END IF
       IF(IsNan(JFZCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PZ",i,JFZCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
       END IF
       IF(IsNan(energydotg) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PE",energydotg,i,JFXCell(i),ug(1),JFYCell(i),ug(2),pregion%grid%vol(i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PE')
       END IF
       IF (ANY(IsNan(JRSGCell(:,i))) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-RSG",i,JRSGCell(:,i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken Reynolds SG')
       endif
       IF(ANY(IsNan(JTSGCell(:,i))) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-TSG",i,JTSGCell(:,i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken TSG')
        endif

       IF(ANY(IsNan(JAlphaPTCell(:,i))) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-Alpha_PT",i,JAlphaPTCell(:,i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken Alpha_PT')
        endif

        pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         + JFXCell(i)
        
        pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         + JFYCell(i)

        pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         + JFZCell(i)

        pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         + energydotg
!---------------------------------------------------------------------------------------
! Compute the interpolated Reynolds Stress components before projecting them back.
!---------------------------------------------------------------------------------------

! I think I will move this to the below if statement section very soon

         if(global%piclPseudoTurbFlag .gt. 0) then
         ! Conservative to primitive variables
           rhog(i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / (1.0_RFREAL - PhiP(i))
           ugas(XCOORD,i) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                           /pRegion%mixt%cv(CV_MIXT_DENS,i)

           ugas(YCOORD,i) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                            /pRegion%mixt%cv(CV_MIXT_DENS,i)

           ugas(ZCOORD,i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

!-----------------------------------------------------------------------
! Trying to move Eulerian field Pseudo-Turb from ppiclF to Rocflu
!----------------------------------------------------------------------

            urelCell(XCOORD,i) = ugas(XCOORD,i) - upCell(XCOORD,i)
            urelCell(YCOORD,i) = ugas(YCOORD,i) - upCell(YCOORD,i)
            urelCell(ZCOORD,i) = ugas(ZCOORD,i) - upCell(ZCOORD,i)

            RelVelMagCell(i) = SQRT(urelCell(XCOORD,i)**2 + urelCell(YCOORD,i)**2 + &
                                    urelCell(ZCOORD,i)**2)

            ! Variables hard-coded but need to fix them later
            ParticleDiameter(i) = 115.0E-6
            viscosityGas(i) =  1.716E-5
            ! End variables hard-coded

            ! Getting an error to access the viscosity variable
            ReynoldsNoCell(i) = RelVelMagCell(i)*ParticleDiameter(i)*rhog(i)/&
                                 viscosityGas(i)
                               ! variable below leads to segfault
!                                    pRegion%mixt%tv(TV_MIXT_MUEL,i)
!
!
            MachNoCell(i) = RelVelMagCell(i)/csF(i)

            MeanKE = 0.5_RFREAL * RelVelMagCell(i)**2

            ! Mehrabadi's terms 
            ! still need to declare and allocate the terms below
            kTilde = 2.0*PhiP(i) + 2.5*PhiP(i)*((1.0-PhiP(i))**3) &
                         * exp(-PhiP(i)*SQRT(ReynoldsNoCell(i)))

            rhoK= kTilde*MeanKE*rhog(i)

            bPar = a/(1.0 + b*exp(-c*ReynoldsNoCell(i))) &
                    *exp((-d*PhiP(i))/(1.0+e*exp(-f*ReynoldsNoCell(i))))

            bPerp = -0.5_RFREAL * bPar

            ! Add trace
            bPar  = bPar  + oneThird
            bPerp = bPerp + oneThird

            !u1 = urelCell(:,i)

            ! I just want to try to check if flipping the sign would affect
            u1(1) = upCell(XCOORD,i) - ugas(XCOORD,i)
            u1(2) = upCell(YCOORD,i) - ugas(YCOORD,i)
            u1(3) = upCell(ZCOORD,i) - ugas(ZCOORD,i)
            u1dot = dot_product(u1, u1) + epsilon(1.0_RFREAL)


            ! Generate an orthonormal set based on max slip component
            dir = maxloc(abs(u1), dim = 1)
            select case (dir)
            case (1)
               ! Max slip in direction `1`
               u2    = -(u1(2) / u1dot) * u1
               u2(2) = 1.0_RFREAL + u2(2)
               u2dot = dot_product(u2, u2) + epsilon(1.0_RFREAL)
            case(2)
               ! Max slip in direction `2`
               u2    = -(u1(1) / u1dot) * u1
               u2(1) = 1.0_RFREAL + u2(1)
               u2dot = dot_product(u2, u2) + epsilon(1.0_RFREAL)
            case(3)
               ! Max slip in direction `23`
               u2    = -(u1(1) / u1dot) * u1
               u2(1) = 1.0_RFREAL + u2(1)
               u2dot = dot_product(u2, u2) + epsilon(1.0_RFREAL)
            end select

            ! Right-hand coordinate system for U3 (cross-product)
            u3(1)  = u1(2) * u2(3) - u1(3) * u2(2)
            u3(2)  = u1(3) * u2(1) - u1(1) * u2(3)
            u3(3)  = u1(1) * u2(2) - u1(2) * u2(1)
            u3dot = dot_product(u3, u3) + epsilon(1.0_RFREAL)

            ! Normalize basis vectors
            u1 = u1 / sqrt(u1dot)
            u2 = u2 / sqrt(u2dot)
            u3 = u3 / sqrt(u3dot)

            ! Construct rotation matrices
            Q(:,1) = u1
            Q(:,2) = u2
            Q(:,3) = u3

            ! Direction `1` parallel to slip velocity (bPar)
            ! Multiply diagonal b^dag by Q^T
            temp(1,1) = bPar *Q (1,1)
            temp(1,2) = bPar *Q (2,1)
            temp(1,3) = bPar *Q (3,1)

            temp(2,1) = bPerp * Q(1,2)
            temp(2,2) = bPerp * Q(2,2)
            temp(2,3) = bPerp * Q(3,2)

            temp(3,1) = bPerp * Q(1,3)
            temp(3,2) = bPerp * Q(2,3)
            temp(3,3) = bPerp * Q(3,3)

            ! Multiply Q by b^dag*Q^T (temp) to get bij tensor
            bij(1,1) = Q(1,1) * temp(1,1) + Q(1,2) * temp(2,1) + Q(1,3) * temp(3,1)
            bij(1,2) = Q(1,1) * temp(1,2) + Q(1,2) * temp(2,2) + Q(1,3) * temp(3,2)
            bij(1,3) = Q(1,1) * temp(1,3) + Q(1,2) * temp(2,3) + Q(1,3) * temp(3,3)

            bij(2,1) = Q(2,1) * temp(1,1) + Q(2,2) * temp(2,1) + Q(2,3) * temp(3,1)
            bij(2,2) = Q(2,1) * temp(1,2) + Q(2,2) * temp(2,2) + Q(2,3) * temp(3,2)
            bij(2,3) = Q(2,1) * temp(1,3) + Q(2,2) * temp(2,3) + Q(2,3) * temp(3,3)

            bij(3,1) = Q(3,1) * temp(1,1) + Q(3,2) * temp(2,1) + Q(3,3) * temp(3,1)
            bij(3,2) = Q(3,1) * temp(1,2) + Q(3,2) * temp(2,2) + Q(3,3) * temp(3,2)
            bij(3,3) = Q(3,1) * temp(1,3) + Q(3,2) * temp(2,3) + Q(3,3) * temp(3,3)

            ! Still need to calculate rhoK

            ! Calculate Reynolds Stress
            reynoldsStress(1,i) = 2.0_RFREAL * rhoK * bij(1,1)
            reynoldsStress(2,i) = 2.0_RFREAL * rhoK * bij(1,2)
            reynoldsStress(3,i) = 2.0_RFREAL * rhoK * bij(1,3)
            reynoldsStress(4,i) = reynoldsStress(2,i)
            reynoldsStress(5,i) = 2.0_RFREAL * rhoK * bij(2,2)
            reynoldsStress(6,i) = 2.0_RFREAL * rhoK * bij(2,3)
            reynoldsStress(7,i) = reynoldsStress(3,i)
            reynoldsStress(8,i) = reynoldsStress(6,i)
            reynoldsStress(9,i) = 2.0_RFREAL * rhoK * bij(3,3)
                          
!
!            ! storing for paraview plotting
            pRegion%mixt%MachNoCell(i) = MachNoCell(i)
            pRegion%mixt%ReynoldsNoCell(i) = ReynoldsNoCell(i)
            pRegion%mixt%RsgVol(i) = reynoldsStress(1,i)
            pRegion%mixt%RsgPar(i) = JRSGCell(1,i)

            pRegion%mixt%FCell(XCOORD,i) = JFXCell(i)
            pRegion%mixt%FCell(YCOORD,i) = JFYCell(i)
            pRegion%mixt%FCell(ZCOORD,i) = JFZCell(i)


  
           ! Subgrid Kinetic Energy
           ! K_sg = 1/(2*rhof) * tr(Rsg), dimension of (nCellsTot)
           ! K_sg is added to the Total Gas Energy term
!           pRegion%mixt%piclKsg(i) = 1.0_RFREAL/(2.0_RFREAL*rhog(i)) &
!                                  * (JRSGCell(1,i) + JRSGCell(5,i) + JRSGCell(9,i))
!
           do j=1,9
             ! \phi_g \rho_g R_sg
             pRegion%mixt%piclPhiRSG(j,i) = reynoldsStress(j,i) * (1.0_RFREAL - PhiP(i))
           end do
         
         endif ! piclPseudoTurbFlag

    END DO !nCells
 
! Commenting below out just for now as I'm cleaning things up


       if(global%piclPseudoTurbFlag .gt. 0) then

         ALLOCATE(varInfoPicl(9), piclcvInfo(9), STAT=errorFlag)
         varInfoPicl = [(j, j=1, 9)]
         piclcvInfo = varInfoPicl
!
         CALL RFLU_ComputeGradCellsWrapper(pRegion,1,9,1,9,varInfoPicl, &   
                                           pRegion%mixt%piclPhiRsg,&                
                                           pRegion%mixt%piclGradPhiRsg)       
                                                                           
         CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,9, &                   
                                           pRegion%mixt%piclGradPhiRsg)       
                                                                           
         ! Limiters are usually designed for scalars, not tensors
         ! and can break symmetry and distort anisotropy
         ! I need to do a comparative run between turning on and off
         CALL RFLU_LimitGradCellsSimple(pRegion,1,9,1,9, &                  
                                        pRegion%mixt%piclPhiRsg,&                   
                                        piclcvInfo,&                        
                                        pRegion%mixt%piclGradPhiRsg) 

         DEALLOCATE(varInfoPicl, piclcvInfo, STAT=errorFlag)
!
!         ALLOCATE(varInfoPicl(3), piclcvInfo(3), STAT=errorFlag)
!         varInfoPicl = [(j, j=1, 3)]
!         piclcvInfo = varInfoPicl
!
!         ! Calculate Gradient of T_g (Gas Temperature)
!         CALL RFLU_ComputeGradCellsWrapper(pRegion,DV_MIXT_TEMP,DV_MIXT_TEMP, &
!                                           1,1,varInfoPicl, pRegion%mixt%dv,&                
!                                           pRegion%mixt%piclgradTg)       
!                                                                           
!         CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,1,pRegion%mixt%piclgradTg)       
!                                                                           
!         CALL RFLU_LimitGradCellsSimple(pRegion,DV_MIXT_TEMP,DV_MIXT_TEMP,1,1, &                  
!                                        pRegion%mixt%dv, piclcvInfo, &                        
!                                        pRegion%mixt%piclgradTg)          
!         DEALLOCATE(varInfoPicl, piclcvInfo, STAT=errorFlag)
!
!
!         ! PTHF = alpha_PT \cdot grad(T_g)
!          JPTHFCell(XCOORD,:) =  JAlphaPTCell(1,:)*pRegion%mixt%piclgradTg(XCOORD,1,:nCells) & 
!                                +JAlphaPTCell(2,:)*pRegion%mixt%piclgradTg(YCOORD,1,:nCells) &
!                                +JAlphaPTCell(3,:)*pRegion%mixt%piclgradTg(ZCOORD,1,:nCells)
!
!          JPTHFCell(YCOORD,:) =  JAlphaPTCell(4,:)*pRegion%mixt%piclgradTg(XCOORD,1,:nCells) & 
!                                +JAlphaPTCell(5,:)*pRegion%mixt%piclgradTg(YCOORD,1,:nCells) &
!                                +JAlphaPTCell(6,:)*pRegion%mixt%piclgradTg(ZCOORD,1,:nCells)
!
!          JPTHFCell(ZCOORD,:) =  JAlphaPTCell(7,:)*pRegion%mixt%piclgradTg(XCOORD,1,:nCells) & 
!                                +JAlphaPTCell(8,:)*pRegion%mixt%piclgradTg(YCOORD,1,:nCells) &
!                                +JAlphaPTCell(9,:)*pRegion%mixt%piclgradTg(ZCOORD,1,:nCells)
!
!           ! Qsg : Subgrid Energy Flux dimension (3, nCells)
!           ! Qsg = rhog*Cp*PTHF + rhog/2*Tsg + ug \cdot Rsg
!           ! PTHF = alpha(9,ncells) \cdot grad(T_g)
!           ! Note: Coefficient at constant pressure (Cp) here is set to the usual constant value
!           !       and this is only valid for GASMODEL 1 (TCPERF).
!           !       If GASMODEL 3 (MIXT_TCPERF) or GASMODEL 7 (SPECIES) are activated, then Cp needs to be
!           !       evaluated for each cell such as below
!           ! Cp = pRegion%mixt%gv(GV_MIXT_CP,pRegion%mixtInput%indCp*icell)
!           
!          Cp = 1004.64_RFREAL
!           
!         Qsg(XCOORD,:) = rhog(:)*Cp*JPTHFCell(XCOORD,:) & 
!                       + rhog(:)/2.0_RFREAL * JTSGCell(XCOORD,:)   &
!            + (JRSGCell(1,:)*ugas(XCOORD,:) + JRSGCell(2,:)*ugas(YCOORD,:) + JRSGCell(3,:)*ugas(ZCOORD,:))
!
!         Qsg(YCOORD,:) = rhog(:)*Cp*JPTHFCell(YCOORD,:) &
!                       + rhog(:)/2.0_RFREAL * JTSGCell(YCOORD,:)   & 
!            + (JRSGCell(4,:)*ugas(XCOORD,:) + JRSGCell(5,:)*ugas(YCOORD,:) + JRSGCell(6,:)*ugas(ZCOORD,:))
!
!         Qsg(ZCOORD,:) = rhog(:)*Cp*JPTHFCell(ZCOORD,:) & 
!                       + rhog(:)/2.0_RFREAL * JTSGCell(ZCOORD,:)   &
!            + (JRSGCell(7,:)*ugas(XCOORD,:) + JRSGCell(8,:)*ugas(YCOORD,:) + JRSGCell(9,:)*ugas(ZCOORD,:))
!
!
!         do j=1,3
!           ! Q_sg -> \phi_g Q_sg
!           pRegion%mixt%piclPhiQsg(j,1:nCells) = Qsg(j,:) * (1.0_RFREAL - PhiP(:))
!         enddo
!
!         IF (ANY(IsNan(pRegion%mixt%piclPhiQsg))) THEN
!             write(*,*) "BROKEN piclPhiQsg"
!             CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,'PPICLF:Broken PhiQSG')
!           END IF
!
!         ALLOCATE(varInfoPicl(3), piclcvInfo(3), STAT=errorFlag)
!         varInfoPicl = [(j, j=1, 3)]
!         piclcvInfo = varInfoPicl
!         ! Calculat Gradient of Q_sg (Subgrid Energy Flux)
!
!         CALL RFLU_ComputeGradCellsWrapper(pRegion,1,3,1,3,varInfoPicl, &   
!                                           pRegion%mixt%piclPhiQsg,&                
!                                           pRegion%mixt%piclGradPhiQsg)       
!                                                                           
!         CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,3, &                   
!                                           pRegion%mixt%piclGradPhiQsg)       
!                                                                           
!         CALL RFLU_LimitGradCellsSimple(pRegion,1,3,1,3, &                  
!                                        pRegion%mixt%piclPhiQsg,&                   
!                                        piclcvInfo,&                        
!                                        pRegion%mixt%piclGradPhiQsg)          
!
!         DEALLOCATE(varInfoPicl, piclcvInfo, STAT=errorFlag)
!
!
!!       ! Now compute Div(\phi_g R_sg)
!        ! 07/23/2025 - Thierry Daoud 
!        ! pRegion%mixt%piclgradPhiRSG dimension is (3,9,nCellsTot)
!        ! nCellsTot = no. actual cells (nCells) + no. dummy (ghost?) cells
!        ! When calculating the gradient, you need to use nCellsTot
!        
!! Div (\phi_g R_sg) - comma denotes partial derivative (,3 -> partial / partial x_3)
!! x-direction: Div(\phi_g R_sg),x = (\phi R_11),1 + (\phi R_12),2 + (\phi R_13),3
!! y-direction: Div(\phi_g R_sg),y = (\phi R_21),1 + (\phi R_22),2 + (\phi R_23),3
!! z-direction: Div(\phi_g R_sg),z = (\phi R_31),1 + (\phi R_32),2 + (\phi R_33),3
!
    DO i = 1,pRegion%grid%nCells
         DivPhiRSG(XCOORD,i) = pregion%grid%vol(i) * &
                              (pRegion%mixt%piclGradPhiRsg(XCOORD,1,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,2,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,3,i))

         DivPhiRSG(YCOORD,i) = pregion%grid%vol(i)* &
                              (pRegion%mixt%piclGradPhiRsg(XCOORD,4,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,5,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,6,i))

         DivPhiRSG(ZCOORD,i) = pregion%grid%vol(i) * &
                              (pRegion%mixt%piclGradPhiRsg(XCOORD,7,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,8,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,9,i))

!
!! Div (\phi_g Q_sg) - comma denotes partial derivative (,3 -> partial / partial x_3)
!! Scalar: Div(\phi_g Q_sg) = (\phi Q_1),1 + (\phi Q_2),2 + (\phi Q_3),3
!         DivPhiQsg(i) =  pRegion%mixt%piclGradPhiQsg(XCOORD,1,i) &
!                       + pRegion%mixt%piclGradPhiQsg(YCOORD,2,i) &
!                       + pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i)
!
! IF (IsNan(DivPhiQsg(i)) .EQV. .TRUE.) THEN
!        write(*,*) "BROKEN- DivPhiQSG, i ", DivPhiQsg(i), i
!        write(*,*) "pregion%grid%vol(i)", pregion%grid%vol(i)
!        write(*,*) "pRegion%mixt%piclGradPhiQsg(XCOORD,1,i)", pRegion%mixt%piclGradPhiQsg(XCOORD,1,i)
!        write(*,*) "pRegion%mixt%piclGradPhiQsg(YCOORD,2,i)", pRegion%mixt%piclGradPhiQsg(YCOORD,2,i)
!        write(*,*) "pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i)", pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i)
!        write(*,*) "rhog(i), ugas(1:3,i)", rhog(i), ugas(1:3,i)
!        write(*,*) "JTSGCell(1:3,i) ", JTSGCell(1:3,i)
!        write(*,*) "JRSGCell(1:9,i)", JRSGCell(:,i)
!        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken QSG')
!endif
!
!
         ! Feedback Div(phi Rsg) to the Fluid Momentum Equations
         pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                          + DivPhiRsg(XCOORD,i) 
         
         pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                          + DivPhiRsg(YCOORD,i) 

         pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                          + DivPhiRsg(ZCOORD,i)
!!                       
!!         ! Feedback Div(phi Qsg) to the Fluid Energy Equation
!!         pRegion%mixt%rhs(CV_MIXT_ENER,i) &
!!                          = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
!!                          + DivPhiQsg(i)
!
        ENDDO

       endif ! piclPseudoTurbFlag
END IF ! global%piclFeedbackFlag

!Due to moving particle integration stuff stoping this for now
DO i = 1,pRegion%grid%nCells
!zero out PhiP
       PhiP(i) = 0.0D0
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       vfP(i) = vfP(i)/pRegion%grid%vol(i)
       PhiP(i) = vfP(i)
!VOL Frac Cap
! Should we keep this???
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       pRegion%mixt%piclVF(i) = PhiP(i) 
END DO
!Deallocate arrays


    IF(pRegion%mixtInput%axiFlag) THEN
      DEALLOCATE(YTEMP,STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    END IF

    DEALLOCATE(rhoF, uxF, uyF, uzF, csF, tpF, &
               SDRX, SDRY, SDRZ, rhsR, pGcX, pGcY, pGcZ, &
               JFXCell, JFYCell, JFZCell, JFECell, &
               domgdx, domgdy, domgdz, drhodx, drhody, drhodz, &
               dpvxF, dpvyF, dpvzF, SDOX, SDOY, SDOZ, PhiP, vfP, &
               dpxF, dpyF, dpzF, &
               STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
!---------------------------------------------------------------  
    DEALLOCATE(JRSGCell, JTSGCell, JAlphaPTCell, JPTHFCell, &
               DivPhiRSG, rhog, ugas, Qsg, DivPhiQsg, &
               JVpUpCell, urelCell, &
               upCell, ReynoldsNoCell, ParticleDiameter, &
               RelVelMagCell, MachNoCell, viscosityGas, reynoldsStress  &
               ,STAT=errorFlag)

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

!---------------------------------------------------------------  
#endif
!PPICLF Integration END

! finalize --------------------------------------------------------------------
  CALL DeregisterFunction(global )
END SUBROUTINE PICL_TEMP_Runge

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************
