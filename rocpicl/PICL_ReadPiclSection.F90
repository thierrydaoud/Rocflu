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
! Purpose: read in user input related to numerical procedure.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = numerical switches and parameters.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadNumericsSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PICL_ReadPiclSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: ReadSection  
  
  IMPLICIT NONE

! ... local variables
  INTEGER :: nVals

  INTEGER :: iReg
  INTEGER, PARAMETER :: NVALS_MAX = 27

  CHARACTER(20) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  CALL RegisterFunction( global,'PICL_ReadPiclSection',__FILE__ )

! specify keywords and search for them

  nVals = NVALS_MAX

  keys( 1)  = 'QSFORCE'
  keys( 2)  = 'AMFORCE'
  keys( 3)  = 'PGFORCE'
  keys( 4)  = 'COLLISIONS'
  keys( 5)  = 'HEATTRANSFER'
  keys( 6)  = 'FEEDBACK'
  keys( 7)  = 'VELOCITYMODEL'
  keys( 8)  = 'USED'
  keys( 9)  = 'QSFLUCT'
  keys(10)  = 'TEMP'
  keys(11)  = 'DU HEATTRANSFER'
  keys(12)  = 'NEIGHBORWIDTH'
  keys(13)  = 'KSP'
  keys(14)  = 'EREST'
  keys(15)  = 'QSFLUCTFILT'
  keys(16)  = 'QSFLUCTFILTADAPT'
  keys(17)  = 'VISCOUSUNSTEADY'
  keys(18)  = 'PERIODICX'
  keys(19)  = 'PERIODICY'
  keys(20)  = 'PERIODICZ'
  keys(21)  = 'SUBBIN'
  keys(22)  = 'PSEUDOTURB'
  keys(23)  = 'VUSTEPS'
! Optional groups of fields in the particle files (ppiclf_user_OutputFields).
! No key may be the start of another: ReadSection matches the start of a line.
  keys(24)  = 'OUTPARTICLE'
  keys(25)  = 'OUTGAS'
  keys(26)  = 'OUTGRADGAS'
  keys(27)  = 'OUTMODEL'
 
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) ) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    global%piclQsFlag = NINT(vals(1))
  END IF ! defined
  
  IF ( defined(2) .EQV. .TRUE. ) THEN 
    global%piclAmFlag = NINT(vals(2))
  END IF ! defined 
  
  IF ( defined(3) .EQV. .TRUE. ) THEN 
    global%piclPgFlag = NINT(vals(3))
  END IF ! defined 
  
  IF ( defined(4) .EQV. .TRUE. ) THEN 
    global%piclCollisionFlag = NINT(vals(4))
  END IF ! defined 
  
  IF ( defined(5) .EQV. .TRUE. ) THEN 
    global%piclHeatTransferFlag = NINT(vals(5))
  END IF ! defined

  IF ( defined(6) .EQV. .TRUE. ) THEN 
    global%piclFeedbackFlag = NINT(vals(6))
  END IF ! defined
  
! 0: particles keep their position and velocity; 1: they move
  IF ( defined(7) .EQV. .TRUE. ) THEN 
    global%piclVelocityModel = NINT(vals(7))
    IF ( global%piclVelocityModel < 0 .OR. &
         global%piclVelocityModel > 1 ) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
        'VELOCITYMODEL must be 0 (particles fixed) or 1 (particles move)')
    END IF
  END IF ! defined 

  IF ( defined(8) .EQV. .TRUE. ) THEN 
    IF (NINT(vals(8)) .EQ. 1) THEN
      global%piclUsed = .TRUE.
    END IF
  END IF ! defined 

  IF ( defined(9) .EQV. .TRUE. ) THEN 
    global%piclQSFluctFlag = NINT(vals(9))
  END IF ! defined 

  IF (defined(10) .EQV. .TRUE. ) THEN
    global%piclTemp = vals(10)
  END IF

  IF (defined(11) .EQV. .TRUE. ) THEN
     global%piclHTUnsteadyFlag = NINT(vals(11))
  END IF

  IF (defined(12) .EQV. .TRUE. ) THEN
    global%piclNeighborWidth = vals(12)
  END IF

  IF (defined(13) .EQV. .TRUE. ) THEN
    global%piclKsp = vals(13)
  END IF

  IF (defined(14) .EQV. .TRUE. ) THEN
    global%piclERest = vals(14)
  END IF

  IF (defined(15) .EQV. .TRUE. ) THEN
    global%piclQsFluctFilterFlag = vals(15)
  END IF

  IF (defined(16) .EQV. .TRUE. ) THEN
    global%piclQsFluctFilterAdaptFlag = vals(16)
  END IF

  IF (defined(17) .EQV. .TRUE. ) THEN
    global%piclViscousUnsteady = vals(17)
  END IF

  IF (defined(18) .EQV. .TRUE. ) THEN
    global%piclPeriodicXFlag = NINT(vals(18))
  END IF
  
  IF (defined(19) .EQV. .TRUE. ) THEN
    global%piclPeriodicYFlag = NINT(vals(19))
  END IF
  
  IF (defined(20) .EQV. .TRUE. ) THEN
    global%piclPeriodicZFlag = NINT(vals(20))
  END IF

  IF (defined(21) .EQV. .TRUE. ) THEN
    global%piclSBNearFlag = vals(21)
  END IF

  IF (defined(22) .EQV. .TRUE. ) THEN
    global%piclPseudoTurbFlag = NINT(vals(22))
  END IF

  IF (defined(23) .EQV. .TRUE. ) THEN
    global%piclVUSteps = NINT(vals(23))
  END IF

  IF (defined(24) .EQV. .TRUE. ) global%piclOutParticle = NINT(vals(24))
  IF (defined(25) .EQV. .TRUE. ) global%piclOutGas      = NINT(vals(25))
  IF (defined(26) .EQV. .TRUE. ) global%piclOutGradGas  = NINT(vals(26))
  IF (defined(27) .EQV. .TRUE. ) global%piclOutModel    = NINT(vals(27))

! The unsteady history is allocated only when VISCOUSUNSTEADY >= 1, and the
! diffusive-unsteady heat transfer is computed in the same routines
  IF ( global%piclViscousUnsteady >= 1 .AND. global%piclVUSteps < 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
      'VUSTEPS (steps of viscous-unsteady history) must be at least 1')
  END IF
  IF ( global%piclHTUnsteadyFlag >= 1 .AND. &
       global%piclViscousUnsteady < 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
      'DU HEATTRANSFER needs VISCOUSUNSTEADY >= 1 (shared history)')
  END IF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PICL_ReadPiclSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadNumericsSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/03/27 12:06:46  haselbac
! Added axiFlag
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.10  2006/10/20 21:20:48  mparmar
! Renamed keyword SORDERBF to ORDERBF
!
! Revision 1.9  2006/08/19 15:38:29  mparmar
! Added reading of mixtInput%spaceOrderBFaces
!
! Revision 1.8  2006/04/07 14:37:44  haselbac
! Added new stencilDimens parameters
!
! Revision 1.7  2006/01/06 22:04:19  haselbac
! Added SDIMENS
!
! Revision 1.6  2005/12/25 15:20:21  haselbac
! Added constrained reconstuction input variables, minor bug fixes (NINT instead of ABS)
!
! Revision 1.5  2005/12/24 21:23:35  haselbac
! Added reading of ICT tolerance
!
! Revision 1.4  2005/10/27 18:53:24  haselbac
! Added input parameter for constrained reconstruction
!
! Revision 1.3  2005/07/11 19:22:20  mparmar
! Added reading of RECONST keyword
!
! Revision 1.2  2005/03/09 14:51:41  haselbac
! Added DIMENS variable
!
! Revision 1.1  2004/12/01 16:50:37  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/09/02 02:34:11  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.18  2004/07/08 02:18:37  haselbac
! Added dissFact, cosmetics
!
! Revision 1.17  2003/12/05 16:54:34  haselbac
! Fixed bug in DO loop ranges
!
! Revision 1.16  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.12  2003/04/09 20:42:54  wasistho
! enable k4 = 0.0
!
! Revision 1.11  2002/09/09 14:03:46  haselbac
! mixtInput now under regions
!
! Revision 1.10  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/07/25 00:38:00  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.8  2002/05/04 16:21:58  haselbac
! Bug fix: spaceOrder was set twice
!
! Revision 1.7  2002/03/26 19:06:47  haselbac
! Added ROCFLU functionality
!
! Revision 1.6  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************

