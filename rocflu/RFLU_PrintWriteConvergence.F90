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
! Purpose: Sum the global forces, mass flow and stop flag over all ranks and
!   print the progress line (unsteady) or convergence line (steady).
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintWriteConvergence.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintWriteConvergence(global)

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI

  USE ModTools, ONLY: MakeNonZero

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: timing
  INTEGER :: errorFlag,nFix,percent,stepsSinceLast,worstProc
  REAL(RFREAL) :: endTime,eta,resNorm,secAverage,secElapsed,secPerStep, &
                  timeDone,wallNow
  REAL(RFREAL), DIMENSION(2) :: worstLocal,worstGlobal
  REAL(RFREAL), DIMENSION(6) :: box,worst
  REAL(RFREAL), DIMENSION(7) :: localVals,globalVals

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_PrintWriteConvergence',__FILE__)

! ******************************************************************************
! Reduce data
! ******************************************************************************

  localVals(1) = global%forceX
  localVals(2) = global%forceY
  localVals(3) = global%forceZ
  localVals(4) = global%massIn
  localVals(5) = global%massOut
  localVals(6) = global%stopRun
  localVals(7) = REAL(global%posFixCount,RFREAL)

  CALL MPI_Allreduce(localVals,globalVals,7,MPI_RFREAL,MPI_SUM,global%mpiComm, &
                     errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag

  global%forceX  = globalVals(1)
  global%forceY  = globalVals(2)
  global%forceZ  = globalVals(3)
  global%massIn  = globalVals(4)
  global%massOut = globalVals(5)
  global%stopRun = globalVals(6)

! ******************************************************************************
! Print and write global data
! ******************************************************************************

! ==============================================================================
! Steady flow 
! ==============================================================================

  IF ( (global%flowType == FLOW_STEADY) .AND. & 
       (global%myProcid == MASTERPROC) ) THEN
    IF ( global%currentIter == 1 ) THEN
      resNorm = 1.0_RFREAL
    ELSE
      resNorm = global%residual/MakeNonZero(global%resInit)
      resNorm = MakeNonZero(resNorm)
    END IF ! global%currentIter

    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,I6,1PE13.4,5E13.4)') SOLVER_NAME, & 
        global%currentIter,LOG10(resNorm), &
        global%forceX,global%forceY,global%forceZ, &
        global%massIn,global%massOut
    END IF ! global%verbLevel

! ==============================================================================
! Unsteady flow 
! ==============================================================================

! Progress line, e.g.
!  [ 41%] t =  2.0306E-04 dt =  8.4559E-08 @ step     3179
!         t/step  1.19E-01 s (avg  1.15E-01 s)  ETA 0:06:43 (total 0:11:20)
! (one line). Progress is the simulated time done in this run over the time
! it has to do; the ETA assumes the rest runs at the average rate so far.
! The first line has no timing yet.

  ELSE IF ( (global%flowType == FLOW_UNSTEADY) .AND. & 
            (global%myProcid==MASTERPROC) ) THEN
    endTime  = global%restartTime &
             + MIN(global%maxTime - global%restartTime,global%runTime)
    timeDone = global%currentTime - global%restartTime
    percent  = 100
    IF ( endTime > global%restartTime ) THEN
      percent = INT(100.0_RFREAL*timeDone/(endTime - global%restartTime))
      percent = MAX(0,MIN(100,percent))
    END IF ! endTime

    wallNow        = MPI_Wtime()
    stepsSinceLast = global%iterSinceRestart - global%progressStepLast
    timing = ''
    IF ( global%progressStepLast > 0 .AND. stepsSinceLast > 0 ) THEN
      secElapsed = wallNow - global%progressWallStart
      secPerStep = (wallNow - global%progressWallLast)/stepsSinceLast
      secAverage = secElapsed/global%iterSinceRestart
      eta        = 0.0_RFREAL
      IF ( timeDone > 0.0_RFREAL ) THEN
        eta = MAX(0.0_RFREAL,secElapsed*(endTime - global%currentTime)/timeDone)
      END IF ! timeDone
      WRITE(timing,'(A,ES9.2,A,ES9.2,A)') '  t/step',secPerStep,' s (avg', &
                                          secAverage,' s)'
      timing = TRIM(timing)//'  ETA '//TRIM(FormatWallTime(eta))// &
               ' (total '//TRIM(FormatWallTime(secElapsed + eta))//')'
    END IF ! global%progressStepLast

    WRITE(STDOUT,'(A,I3,A,ES11.4,A,ES11.4,A,I8,A)') ' [',percent,'%] t =', &
      global%currentTime,' dt =',global%dtMin,' @ step',global%currentStep, &
      TRIM(timing)

    global%progressWallLast = wallNow
    global%progressStepLast = global%iterSinceRestart
  END IF ! global%flowType

! ******************************************************************************
! Negative-pressure fixes since the last call, e.g.
!         negative pressure fixed 6820 times, cells in x [ 1.002E-01, 1.002E-01]
!           y [-8.381E-04, 8.381E-04] z [ 1.000E-04, 9.000E-04]  (one line)
!         lowest p -2.523E+03 Pa, set to  7.507E+01 Pa (rho  1.877E-02) at
!           ( 1.002E-01,-3.644E-05, 5.000E-04)  (one line)
! ******************************************************************************

  nFix = NINT(globalVals(7))
  IF ( nFix > 0 ) THEN
    CALL MPI_Allreduce(global%posFixBox,box,6,MPI_RFREAL,MPI_MIN, &
                       global%mpiComm,errorFlag)
    worstLocal(1) = global%posFixWorst(1)
    worstLocal(2) = REAL(global%myProcid,RFREAL)
    CALL MPI_Allreduce(worstLocal,worstGlobal,1,MPI_2DOUBLE_PRECISION, &
                       MPI_MINLOC,global%mpiComm,errorFlag)
    worstProc = NINT(worstGlobal(2))
    worst = global%posFixWorst
    CALL MPI_Bcast(worst,6,MPI_RFREAL,worstProc,global%mpiComm,errorFlag)

    IF ( global%myProcid == MASTERPROC ) THEN
      WRITE(STDOUT,'(8X,A,I0,A,3(1X,A,ES10.3,A,ES10.3,A)/8X,A,ES10.3,A, &
                    &ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A)') &
        'negative pressure fixed ',nFix,' times, cells in', &
        'x [',box(1),',',-box(4),']','y [',box(2),',',-box(5),']', &
        'z [',box(3),',',-box(6),']', &
        'lowest p ',worst(1),' Pa, set to ',worst(2),' Pa (rho ',worst(3), &
        ') at (',worst(4),',',worst(5),',',worst(6),')'
    END IF ! global%myProcid

    global%posFixCount = 0
    global%posFixBox   = HUGE(1.0_RFREAL)
    global%posFixWorst = HUGE(1.0_RFREAL)
  END IF ! nFix

! ******************************************************************************
! Start
! ******************************************************************************

  CALL DeregisterFunction(global)

CONTAINS

! Wall time in seconds as H:MM:SS
  FUNCTION FormatWallTime(seconds) RESULT(hms)
    REAL(RFREAL), INTENT(IN) :: seconds
    CHARACTER(16) :: hms
    INTEGER :: total

    total = NINT(seconds)
    WRITE(hms,'(I0,A,I2.2,A,I2.2)') total/3600,':',MOD(total/60,60),':', &
                                    MOD(total,60)
  END FUNCTION FormatWallTime

END SUBROUTINE RFLU_PrintWriteConvergence

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintWriteConvergence.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.1  2005/04/15 15:07:10  haselbac
! Initial revision
!
! ******************************************************************************

