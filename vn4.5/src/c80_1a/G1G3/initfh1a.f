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
CLL  SUBROUTINE INIT_FLH --------------------------------------
CLL
CLL  Written by D. Robinson
CLL
CLL  Model   Date     Modification history
CLL version
CLL   3.4   08/09/94  New routine.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 4  5/2/92
CLL
CLL  System component: R30
CLL
CLL  System task: F3
CLL
CLL  Purpose:
CLL           Initialises the fixed length header to IMDI except
CLL           for all array dimensions which are set to 1.
CLL
CLL  Documentation:
CLL           Unified Model Documentation Paper No F3
CLL           Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE INIT_FLH (FIXHD,LEN_FIXHD)

      IMPLICIT NONE

      INTEGER
     & LEN_FIXHD        ! IN    Length of fixed length header
     &,FIXHD(LEN_FIXHD) ! INOUT Fixed length header

C Local arrays:------------------------------------------------
C None
C -------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C External subroutines called:---------------------------------
C None
C Local variables:---------------------------------------------
      INTEGER J
C -------------------------------------------------------------


! 1.0 Initialise to IMDI
      DO J = 1,LEN_FIXHD
        FIXHD(J) = IMDI
      ENDDO

! 2.0 Set all array dimensions to 1
      FIXHD(101) = 1     !  Integer Constants
      FIXHD(106) = 1     !  Real Constants
      FIXHD(111) = 1     !  1st dim - Level dependent constants
      FIXHD(112) = 1     !  2nd dim - Level dependent constants
      FIXHD(116) = 1     !  1st dim - Row dependent constants
      FIXHD(117) = 1     !  2nd dim - Row dependent constants
      FIXHD(121) = 1     !  1st dim - Column dependent constants
      FIXHD(122) = 1     !  2nd dim - Column dependent constants
      FIXHD(126) = 1     !  1st dim - Field of constants
      FIXHD(127) = 1     !  2nd dim - Field of constants
      FIXHD(131) = 1     !  Extra constants
      FIXHD(136) = 1     !  Temp History file
      FIXHD(141) = 1     !  Compressed field Index 1
      FIXHD(143) = 1     !  Compressed field Index 2
      FIXHD(145) = 1     !  Compressed field Index 3
      FIXHD(151) = 1     !  1st dim - Lookup Table
      FIXHD(152) = 1     !  2nd dim - Lookup Table
      FIXHD(161) = 1     !  Data

      RETURN
      END
