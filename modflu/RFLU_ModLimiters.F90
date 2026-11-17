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
! ******************************************************************************
!
! Purpose: Suite of routines related to limiter functions.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModLimiters.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModLimiters

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ComputeLimiterBarthJesp, &
            RFLU_LimitGradVelocityBarthJesp, &
            RFLU_ComputeLimiterVenkat, &
            RFLU_CreateLimiter, &
            RFLU_DestroyLimiter, & 
            RFLU_LimitGradCells, &
            RFLU_LimitGradCellsSimple

        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModLimiters.F90,v $ $Revision: 1.1.1.1 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! ******************************************************************************
!
! Purpose: Compute limiters for gradients at cell centers using Barth-
!   Jespersen limiter function.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar	Beginning index in var array
!   iEndVar	Beginning index in var array
!   iBegGrad	Beginning index in grad array
!   iEndGrad	Beginning index in grad array
!   var		Variables to be limited
!   grad	Gradients of variables
!
! Output:
!   lim        	Limiter
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeLimiterBarthJesp(pRegion,iBegVar,iEndVar,iBegGrad, &
                                        iEndGrad,var,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: iBegGrad,iBegVar,iEndGrad,iEndVar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim,var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,errorFlag,icg,ifg,iGrad,iVar
  REAL(RFREAL) :: dx1,dx2,dy1,dy2,dz1,dz2,d1max1,d1max2,d1min1,d1min2,d21, &
                  d22,term,var1,var2,xc,yc,zc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: varMax,varMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLimiterBarthJesp',__FILE__ )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! iEndVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************    
! Allocate temporary memory
! ******************************************************************************    

  ALLOCATE(varMin(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  ALLOCATE(varMax(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************    
! Initialize
! ******************************************************************************    

  DO icg = 1,pGrid%nCellsTot                       
    DO iVar = iBegVar,iEndVar
      varMax(iVar,icg) = var(iVar,icg) 
      varMin(iVar,icg) = var(iVar,icg)
    END DO ! iVar
    
    DO iGrad = iBegGrad,iEndGrad
      lim(iGrad,icg) = 1.0_RFREAL
    END DO ! iGrad
  END DO ! icg 

! ******************************************************************************    
! Compute varMax and varMin 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)                       

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      varMin(iVar,c1) = MIN(varMin(iVar,c1),var2)
      varMin(iVar,c2) = MIN(varMin(iVar,c2),var1)
      varMax(iVar,c1) = MAX(varMax(iVar,c1),var2)
      varMax(iVar,c2) = MAX(varMax(iVar,c2),var1)
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Loop over faces and compute limiter 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)
  
    dx1 = pGrid%cofg(XCOORD,c1) - xc
    dy1 = pGrid%cofg(YCOORD,c1) - yc
    dz1 = pGrid%cofg(ZCOORD,c1) - zc

    dx2 = pGrid%cofg(XCOORD,c2) - xc
    dy2 = pGrid%cofg(YCOORD,c2) - yc
    dz2 = pGrid%cofg(ZCOORD,c2) - zc

    iGrad = iBegGrad

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      d1min1 = varMin(iVar,c1) - var1
      d1max1 = varMax(iVar,c1) - var1
      d1min2 = varMin(iVar,c2) - var2
      d1max2 = varMax(iVar,c2) - var2

      d21 = grad(XCOORD,iGrad,c1)*dx1 &
          + grad(YCOORD,iGrad,c1)*dy1 &
          + grad(ZCOORD,iGrad,c1)*dz1
      d22 = grad(XCOORD,iGrad,c2)*dx2 &
          + grad(YCOORD,iGrad,c2)*dy2 &
          + grad(ZCOORD,iGrad,c2)*dz2
            
      IF ( d21 > 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1max1/d21)
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ELSEIF ( d21 < 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1min1/d21)
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      END IF ! d21       
                     
      IF ( d22 > 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1max2/d22)
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ELSEIF ( d22 < 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1min2/d22)
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      END IF ! d22       

      iGrad = iGrad + 1
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Deallocate temporary memory
! ******************************************************************************    

  DEALLOCATE(varMin,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  DEALLOCATE(varMax,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeLimiterBarthJesp







! ******************************************************************************
!
! Purpose: Compute limiters for gradients at cell centers using Venkatakrishnan
!   limiter function.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar	Beginning index in var array
!   iEndVar	Beginning index in var array
!   iBegGrad	Beginning index in grad array
!   iEndGrad	Beginning index in grad array
!   var		Variables to be limited
!   grad	Gradients of variables
!
! Output:
!   lim        	Limiter
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeLimiterVenkat(pRegion,iBegVar,iEndVar,iBegGrad, &
                                     iEndGrad,var,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: iBegGrad,iBegVar,iEndGrad,iEndVar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim,var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,errorFlag,icg,ifg,iGrad,iVar
  REAL(RFREAL), PARAMETER :: THRD = 1.0_RFREAL/3.0_RFREAL
  REAL(RFREAL), PARAMETER :: TINY = 1.0E-12_RFREAL 
  REAL(RFREAL) :: denom,ds1,ds2,dx1,dx2,dy1,dy2,dz1,dz2,d1max1,d1max2,d1min1, &
                  d1min2,d21,d22,epsq1,epsq2,numer,term,var1,var2,venkatLimK, &
                  xc,yc,zc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: varMax,varMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLimiterVenkat',__FILE__ )

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! TEMPORARY
  venkatLimK = 1.0_RFREAL
! END TEMPORARY

! ******************************************************************************    
! Allocate temporary memory
! ******************************************************************************    

  ALLOCATE(varMin(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  ALLOCATE(varMax(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************    
! Initialize
! ******************************************************************************    

  DO icg = 1,pGrid%nCellsTot                       
    DO iVar = iBegVar,iEndVar
      varMax(iVar,icg) = var(iVar,icg) 
      varMin(iVar,icg) = var(iVar,icg)
    END DO ! iVar
    
    DO iGrad = iBegGrad,iEndGrad
      lim(iGrad,icg) = 1.0_RFREAL
    END DO ! iGrad
  END DO ! icg 

! ******************************************************************************    
! Compute varMax and varMin 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)                       

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      varMin(iVar,c1) = MIN(varMin(iVar,c1),var2)
      varMin(iVar,c2) = MIN(varMin(iVar,c2),var1)
      varMax(iVar,c1) = MAX(varMax(iVar,c1),var2)
      varMax(iVar,c2) = MAX(varMax(iVar,c2),var1)
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Loop over faces and compute limiter 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)
  
    dx1 = pGrid%cofg(XCOORD,c1) - xc
    dy1 = pGrid%cofg(YCOORD,c1) - yc
    dz1 = pGrid%cofg(ZCOORD,c1) - zc

    dx2 = pGrid%cofg(XCOORD,c2) - xc
    dy2 = pGrid%cofg(YCOORD,c2) - yc
    dz2 = pGrid%cofg(ZCOORD,c2) - zc

    ds1   = pGrid%vol(c1)**THRD 
    ds2   = pGrid%vol(c2)**THRD 
    epsq1 = (venkatLimK*ds1)*(venkatLimK*ds1)*(venkatLimK*ds1)
    epsq2 = (venkatLimK*ds2)*(venkatLimK*ds2)*(venkatLimK*ds2)

    iGrad = iBegGrad

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      d1min1 = varMin(iVar,c1) - var1
      d1max1 = varMax(iVar,c1) - var1
      d1min2 = varMin(iVar,c2) - var2
      d1max2 = varMax(iVar,c2) - var2

      d21 = grad(XCOORD,iGrad,c1)*dx1 &
          + grad(YCOORD,iGrad,c1)*dy1 &
          + grad(ZCOORD,iGrad,c1)*dz1
      d22 = grad(XCOORD,iGrad,c2)*dx2 &
          + grad(YCOORD,iGrad,c2)*dy2 &
          + grad(ZCOORD,iGrad,c2)*dz2
            
      IF ( d21 > 0.0_RFREAL ) THEN
        numer = (d1max1*d1max1+epsq1)*d21 + 2*d21*d21*d1max1
        denom = d21*(d1max1*d1max1 + 2*d21*d21 + d1max1*d21 + epsq1)
        term  = numer/(denom+TINY)         
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ELSEIF ( d21 < 0.0_RFREAL ) THEN
        numer = (d1min1*d1min1+epsq1)*d21 + 2*d21*d21*d1min1
        denom = d21*(d1min1*d1min1 + 2*d21*d21 + d1min1*d21 + epsq1)
        term  = numer/(denom+TINY) 
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ENDIF ! d21       
                     
      IF ( d22 > 0.0_RFREAL ) THEN
        numer = (d1max2*d1max2+epsq2)*d22 + 2*d22*d22*d1max2
        denom = d22*(d1max2*d1max2 + 2*d22*d22 + d1max2*d22 + epsq2)
        term  = numer/(denom+TINY) 
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ELSEIF ( d22 < 0.0_RFREAL ) THEN
        numer = (d1min2*d1min2+epsq2)*d22 + 2*d22*d22*d1min2
        denom = d22*(d1min2*d1min2 + 2*d22*d22 + d1min2*d22 + epsq2)
        term  = numer/(denom+TINY) 
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ENDIF ! d22       

      iGrad = iGrad + 1
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Deallocate temporary memory
! ******************************************************************************    

  DEALLOCATE(varMin,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  DEALLOCATE(varMax,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeLimiterVenkat








! ******************************************************************************
!
! Purpose: Limit the cell velocity gradients with Barth-Jespersen along the
!   principal axes of the strain rate, so that the result does not depend on
!   the coordinate frame.
!
! Description: In each cell, q_1, q_2, q_3 are the eigenvectors of (G+G^T)/2,
!   with G the unlimited velocity gradient. The Barth-Jespersen factor alpha_j
!   of the projection q_j.u is found as for a scalar; axes with equal
!   eigenvalues are not unique and share the smaller factor. The gradient
!   becomes Q diag(alpha) Q^T G.
!
! Input:
!   pRegion     Pointer to region
!   iVelVar     Index of x-velocity in var (y- and z-velocity follow)
!   iVelGrad    Index of x-velocity gradient in grad (y and z follow)
!   var         Primitive variables
!   grad        Unlimited gradients
!
! Output:
!   grad        Velocity gradients replaced by the limited ones
!   lim         Velocity entries set to 1, so RFLU_LimitGradCells keeps them
!
! Notes:
!   1. Directional limiting as in Hajduk, Kuzmin and Aizinger, Frame-invariant
!      directional vector limiters for discontinuous Galerkin methods,
!      Ergebnisberichte Angew. Math. 588, TU Dortmund (2018), with the strain
!      rate axes of Maire (2011) and Zeng and Scovazzi (2014).
!   2. Face offsets and bounds as in RFLU_ComputeLimiterBarthJesp; with Q = I
!      the result equals componentwise Barth-Jespersen.
!
! ******************************************************************************

SUBROUTINE RFLU_LimitGradVelocityBarthJesp(pRegion,iVelVar,iVelGrad,var, &
                                           grad,lim)

  IMPLICIT NONE

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iVelGrad,iVelVar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim,var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

! Relative eigenvalue difference below which strain-rate axes count as equal
  REAL(RFREAL), PARAMETER :: EQUAL_EIG_TOL = 1.0E-3_RFREAL

  LOGICAL :: equal12,equal23
  INTEGER :: icg,iComp,ifg,iSide,icl,errorFlag,j
  REAL(RFREAL) :: eigScale,minAlpha,proj
  REAL(RFREAL) :: du(3),offset(3),shift(3)
  REAL(RFREAL) :: gradVel(3,3),rotated(3,3)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: alpha,diffMax,diffMin,eigVals
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: axes
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_LimitGradVelocityBarthJesp',__FILE__)

  pGrid => pRegion%grid

  ALLOCATE(axes(3,3,pGrid%nCellsTot),eigVals(3,pGrid%nCellsTot), &
           alpha(3,pGrid%nCellsTot),diffMin(3,pGrid%nCellsTot), &
           diffMax(3,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'axes')
  END IF ! global%error

! ******************************************************************************
! Principal axes of the strain rate of the unlimited gradient
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    DO iComp = 1,3
      gradVel(iComp,1:3) = grad(XCOORD:ZCOORD,iVelGrad+iComp-1,icg)
    END DO ! iComp

    CALL RFLU_SymmetricEigen3(0.5_RFREAL*(gradVel + TRANSPOSE(gradVel)), &
                              axes(:,:,icg),eigVals(:,icg))

    alpha(:,icg)   = 1.0_RFREAL
    diffMin(:,icg) = 0.0_RFREAL
    diffMax(:,icg) = 0.0_RFREAL
  END DO ! icg

! ******************************************************************************
! Range of the neighbours' velocity along each axis of the cell
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    du(1:3) = var(iVelVar:iVelVar+2,pGrid%f2c(2,ifg)) &
            - var(iVelVar:iVelVar+2,pGrid%f2c(1,ifg))

    DO iSide = 1,2
      icl = pGrid%f2c(iSide,ifg)

      DO j = 1,3
        proj = DOT_PRODUCT(axes(:,j,icl),du)
        IF ( iSide == 2 ) THEN
          proj = -proj
        END IF ! iSide

        diffMin(j,icl) = MIN(diffMin(j,icl),proj)
        diffMax(j,icl) = MAX(diffMax(j,icl),proj)
      END DO ! j
    END DO ! iSide
  END DO ! ifg

! ******************************************************************************
! Barth-Jespersen factor along each axis
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    DO iSide = 1,2
      icl = pGrid%f2c(iSide,ifg)

      offset(1:3) = pGrid%cofg(XCOORD:ZCOORD,icl) - pGrid%fc(XCOORD:ZCOORD,ifg)

      DO iComp = 1,3
        gradVel(iComp,1:3) = grad(XCOORD:ZCOORD,iVelGrad+iComp-1,icl)
      END DO ! iComp

      shift = MATMUL(gradVel,offset)

      DO j = 1,3
        proj = DOT_PRODUCT(axes(:,j,icl),shift)

        IF ( proj > 0.0_RFREAL ) THEN
          alpha(j,icl) = MIN(alpha(j,icl),diffMax(j,icl)/proj)
        ELSE IF ( proj < 0.0_RFREAL ) THEN
          alpha(j,icl) = MIN(alpha(j,icl),diffMin(j,icl)/proj)
        END IF ! proj
      END DO ! j
    END DO ! iSide
  END DO ! ifg

! ******************************************************************************
! Axes with equal eigenvalues share the smaller factor
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    eigScale = EQUAL_EIG_TOL*MAXVAL(ABS(eigVals(:,icg)))

    equal12 = ( eigVals(2,icg) - eigVals(1,icg) <= eigScale )
    equal23 = ( eigVals(3,icg) - eigVals(2,icg) <= eigScale )

    IF ( equal12 .AND. equal23 ) THEN
      alpha(:,icg) = MINVAL(alpha(:,icg))
    ELSE IF ( equal12 ) THEN
      minAlpha = MIN(alpha(1,icg),alpha(2,icg))
      alpha(1:2,icg) = minAlpha
    ELSE IF ( equal23 ) THEN
      minAlpha = MIN(alpha(2,icg),alpha(3,icg))
      alpha(2:3,icg) = minAlpha
    END IF ! equal12
  END DO ! icg

! ******************************************************************************
! Limited gradient Q diag(alpha) Q^T G
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    DO iComp = 1,3
      gradVel(iComp,1:3) = grad(XCOORD:ZCOORD,iVelGrad+iComp-1,icg)
    END DO ! iComp

    rotated = MATMUL(TRANSPOSE(axes(:,:,icg)),gradVel)
    DO j = 1,3
      rotated(j,1:3) = alpha(j,icg)*rotated(j,1:3)
    END DO ! j
    gradVel = MATMUL(axes(:,:,icg),rotated)

    DO iComp = 1,3
      grad(XCOORD:ZCOORD,iVelGrad+iComp-1,icg) = gradVel(iComp,1:3)
    END DO ! iComp

    lim(iVelGrad:iVelGrad+2,icg) = 1.0_RFREAL
  END DO ! icg

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(axes,eigVals,alpha,diffMin,diffMax,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'axes')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_LimitGradVelocityBarthJesp








! ******************************************************************************
!
! Purpose: Eigenvalues and eigenvectors of a symmetric 3x3 matrix.
!
! Description: Cyclic Jacobi rotations until the off-diagonal part vanishes
!   relative to the matrix.
!
! Input:
!   matrix      Symmetric matrix
!
! Output:
!   vectors     Orthonormal eigenvectors (columns)
!   values      Eigenvalues in ascending order, matching the columns
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SymmetricEigen3(matrix,vectors,values)

  IMPLICIT NONE

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: matrix(3,3)
  REAL(RFREAL), INTENT(OUT) :: values(3),vectors(3,3)

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: MAX_SWEEPS = 50

  INTEGER :: i,iMin,ip,iq,iSweep
  REAL(RFREAL) :: cosine,held,offDiag,sine,tangent,theta,total
  REAL(RFREAL) :: column(3),reduced(3,3),rotation(3,3)

! ******************************************************************************
! Start
! ******************************************************************************

  reduced = matrix
  vectors = 0.0_RFREAL
  DO i = 1,3
    vectors(i,i) = 1.0_RFREAL
  END DO ! i

  total = SUM(reduced**2)

  DO iSweep = 1,MAX_SWEEPS
    offDiag = reduced(1,2)**2 + reduced(1,3)**2 + reduced(2,3)**2
    IF ( offDiag <= 1.0E-30_RFREAL*total ) THEN
      EXIT
    END IF ! offDiag

    DO ip = 1,2
      DO iq = ip+1,3
        IF ( reduced(ip,iq) /= 0.0_RFREAL ) THEN
          theta = (reduced(iq,iq) - reduced(ip,ip))/(2.0_RFREAL*reduced(ip,iq))
          IF ( ABS(theta) > 1.0E+150_RFREAL ) THEN
            tangent = 0.5_RFREAL/theta
          ELSE
            tangent = SIGN(1.0_RFREAL,theta)/(ABS(theta) + SQRT(theta**2 + 1.0_RFREAL))
          END IF ! theta
          cosine = 1.0_RFREAL/SQRT(tangent**2 + 1.0_RFREAL)
          sine   = tangent*cosine

          rotation = 0.0_RFREAL
          DO i = 1,3
            rotation(i,i) = 1.0_RFREAL
          END DO ! i
          rotation(ip,ip) =  cosine
          rotation(iq,iq) =  cosine
          rotation(ip,iq) =  sine
          rotation(iq,ip) = -sine

          reduced = MATMUL(TRANSPOSE(rotation),MATMUL(reduced,rotation))
          reduced(ip,iq) = 0.0_RFREAL
          reduced(iq,ip) = 0.0_RFREAL
          vectors = MATMUL(vectors,rotation)
        END IF ! reduced
      END DO ! iq
    END DO ! ip
  END DO ! iSweep

  DO i = 1,3
    values(i) = reduced(i,i)
  END DO ! i

! Sort ascending
  DO i = 1,2
    iMin = i - 1 + MINLOC(values(i:3),DIM=1)
    IF ( iMin /= i ) THEN
      held         = values(i)
      values(i)    = values(iMin)
      values(iMin) = held
      column          = vectors(:,i)
      vectors(:,i)    = vectors(:,iMin)
      vectors(:,iMin) = column
    END IF ! iMin
  END DO ! i

  END SUBROUTINE RFLU_SymmetricEigen3








! ******************************************************************************
!
! Purpose: Create limiter.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iBegGrad	Beginning index of gradient array
!   iEndGrad	Ending index of gradient array
!
! Output: 
!   lim		Limiter
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_CreateLimiter(pRegion,iBegGrad,iEndGrad,lim)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  CALL RegisterFunction(global,'RFLU_CreateLimiter', &
                        __FILE__)

! ******************************************************************************
! Create limiter
! ******************************************************************************

  ALLOCATE(lim(iBegGrad:iEndGrad,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'lim')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreateLimiter







! ******************************************************************************
!
! Purpose: Destroy limiter.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_DestroyLimiter(pRegion,lim)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyLimiter', &
                        __FILE__)

! ******************************************************************************
! Destroy limiter
! ******************************************************************************

  DEALLOCATE(lim,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'lim')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyLimiter







! ******************************************************************************
!
! Purpose: Apply limiter to gradient.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iBegGrad	Beginning index of grad array
!   iEndGrad	Ending index of grad array
!   grad	Gradient array
!   lim		Limiter
!
! Output: 
!   grad	Limited gradient array
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_LimitGradCells(pRegion,iBegGrad,iEndGrad,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,iGrad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  CALL RegisterFunction(global,'RFLU_LimitGradCells', &
                        __FILE__)

! ******************************************************************************
! Apply limiter
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot                       
    DO iGrad = iBegGrad,iEndGrad
      grad(XCOORD,iGrad,icg) = lim(iGrad,icg)*grad(XCOORD,iGrad,icg) 
      grad(YCOORD,iGrad,icg) = lim(iGrad,icg)*grad(YCOORD,iGrad,icg)
      grad(ZCOORD,iGrad,icg) = lim(iGrad,icg)*grad(ZCOORD,iGrad,icg) 
    END DO ! iVar
  END DO ! icg 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_LimitGradCells








! ******************************************************************************
!
! Purpose: Limit gradients of any vector or scalar at cell centers 
!   to ensure positive face quantities.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables 
!   varInfo     Information on variables
!   grad        Gradients of variables at cell centers
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_LimitGradCellsSimple(pRegion,iBegVar,iEndVar,iBegGrad, &
                                     iEndGrad,var,varInfo,grad)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  INTEGER, DIMENSION(:), POINTER :: varInfo    
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,icl,ict,ifg,ifl,iGrad,iPatch,iVar,nFacesPerCell
  INTEGER, DIMENSION(:,:,:), POINTER :: pC2f
  REAL(RFREAL) :: dx,dy,dz,varc,varf,xc,yc,zc,xfc,yfc,zfc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_LimitGradCellsSimple',__FILE__ )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! iEndVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Loop over cells
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    ict = pGrid%cellGlob2Loc(1,icg) ! cell type
    icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

    SELECT CASE ( ict ) 
      CASE ( CELL_TYPE_TET ) 
        nFacesPerCell = 4
        pC2f => pGrid%tet2f              
      CASE ( CELL_TYPE_HEX ) 
        nFacesPerCell = 6
        pC2f => pGrid%hex2f              
      CASE ( CELL_TYPE_PRI ) 
        nFacesPerCell = 5    
        pC2f => pGrid%pri2f                         
      CASE ( CELL_TYPE_PYR ) 
        nFacesPerCell = 5
        pC2f => pGrid%pyr2f                
      CASE DEFAULT  
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict

    xc = pGrid%cofg(XCOORD,icg)
    yc = pGrid%cofg(YCOORD,icg)
    zc = pGrid%cofg(ZCOORD,icg)                    

! ==============================================================================
!   Loop over positive-definite variables
! ==============================================================================     

    iGrad = iBegGrad    

    DO iVar = iBegVar,iEndVar              
      IF ( varInfo(iVar) == VAR_INFO_POS ) THEN 
        varc = var(iVar,icg)

        faceLoop: DO ifl = 1,nFacesPerCell
          iPatch = pC2f(1,ifl,icl)
          ifg    = pC2f(2,ifl,icl)

          IF ( iPatch == 0 ) THEN ! Interior face
            dx = pGrid%fc(XCOORD,ifg) - xc
            dy = pGrid%fc(YCOORD,ifg) - yc
            dz = pGrid%fc(ZCOORD,ifg) - zc                         
          ELSE ! Boundary face
            pPatch => pRegion%patches(iPatch)

            dx = pPatch%fc(XCOORD,ifg) - xc
            dy = pPatch%fc(YCOORD,ifg) - yc
            dz = pPatch%fc(ZCOORD,ifg) - zc                             
          END IF ! pC2f 

          varf = varc + grad(XCOORD,iGrad,icg)*dx &
                      + grad(YCOORD,iGrad,icg)*dy &
                      + grad(ZCOORD,iGrad,icg)*dz

          IF ( varf <= 0.0_RFREAL ) THEN 
            grad(XCOORD,iGrad,icg) = 0.0_RFREAL
            grad(YCOORD,iGrad,icg) = 0.0_RFREAL
            grad(ZCOORD,iGrad,icg) = 0.0_RFREAL

            EXIT faceLoop
          END IF ! varf            

        END DO faceLoop                            
      END IF ! iVar

      iGrad = iGrad + 1        
    END DO ! iVar
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_LimitGradCellsSimple





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModLimiters

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModLimiters.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/04/27 15:10:55  haselbac
! Many bug fixes and clean-up
!
! Revision 1.2  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.1  2005/07/11 19:33:47  mparmar
! Initial revision
!
! ******************************************************************************

