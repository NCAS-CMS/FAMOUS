C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C Use, duplication or disclosure of this code is subject to the
C restrictions as set forth in the contract.
C
C                Meteorological Office
C                London Road
C                BRACKNELL
C                Berkshire UK
C                RG12 2SZ
C 
C If no contract has been raised with this copy of the code, the use,
C duplication or disclosure of it is strictly prohibited.  Permission
C to do so must first be obtained in writing from the Head of Numerical
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C
      Subroutine READWRITD(UNIT,ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Purpose: Reads information for time-step control of WRITD1 from
!          NAMLST file
!
! Current code owner: S.J.Swarbrick
!
!  Model               Modification history:
! version  Date
! -------  ----        --------------------
!   3.4    28/07/94    Original code - S.J.Swarbrick
!  3.5  08/06/95  Add UNIT - moved from INITIAL to READCNTL. RTHBarnes
!  4.1  31/05/96  Add L_WRIT_WAVSTEP for wave sub-model. RTHBarnes. 
!
! Code description:
!   FORTRAN 77 + common FORTRAN 90 extensions. Written to UM
!   programming standards vn. 7.
!
!
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control of WRITD1
!
! Code owner: S.J.Swarbrick
!
! History:
! Version  Date       Comment
! -------  ----       -------
! 3.4      28/07/94   Original code - S.J.Swarbrick
!  4.1  28/07/94  Introduce Wave sub-model.  RTHBarnes.    
!
!                              Switches which
!                              activate WRITD1:
      LOGICAL L_WRIT_ATMSTEP     ! in ATMSTEP
      LOGICAL L_WRIT_DYN         ! in ATMDYN
      LOGICAL L_WRIT_PHY         ! in ATMPHY
      LOGICAL L_WRIT_OCNSTEP     ! in OCNSTEP
      LOGICAL L_WRIT_WAVSTEP     ! in WAV_STEP  
      LOGICAL L_WRIT_INIT        ! in INITDUM
!
      INTEGER WRITD1_TEST_PREV   ! Used within time-step
      INTEGER WRITD1_TEST        !     control mechanism.
      INTEGER T_WRITD1_START     ! Initial time step for writing dump
      INTEGER T_WRITD1_END       ! Final time step for writing dump
      INTEGER T_WRITD1_INT       ! Interval at which dumps written
!
      COMMON/WTESTP/ WRITD1_TEST_PREV
      COMMON/WRITED1/ T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT
     &              ,L_WRIT_DYN      ,L_WRIT_PHY    ,L_WRIT_INIT
     &              ,L_WRIT_ATMSTEP  ,L_WRIT_OCNSTEP,L_WRIT_WAVSTEP 
!
!
      CHARACTER*80 CMESSAGE  ! Error return message
!
      INTEGER ICODE          ! Indicator for normal or error return
      INTEGER UNIT           ! Namelist input unit no.
!
! --------------------------------------------------------------------
!
      NAMELIST/WRITD1/L_WRIT_ATMSTEP,L_WRIT_DYN  ,L_WRIT_INIT ,
     &                L_WRIT_OCNSTEP,L_WRIT_WAVSTEP,L_WRIT_PHY ,  
     &                T_WRITD1_START,T_WRITD1_END  ,T_WRITD1_INT    
!
!  Read control variables for WRITD1
!
      READ(UNIT,WRITD1,ERR=99)
!
!  Normal return
!
      ICODE=0
      RETURN
!
!  Error return
!
  99  ICODE=1
      CMESSAGE='READWRITD: error reading namelist WRITD1'
      RETURN
      END
