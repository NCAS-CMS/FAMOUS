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
CLL  SUBROUTINE EXPAND32B--------------------------------------
CLL
CLL  Purpose: Expands from 32 to 64 bit for dump reading routines.
CLL
CLL MC          <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.3  08/04/94  Added to avoid problems in readdump. M.Carter
CLL   4.5  28/10/98  Introduce Single Column Model. J-C Thil.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  Logical component: R30
CLL
CLL  System task: F3
CLL
CLL  Documentation: Unified Model Documentation Paper No F3
CLL                 Version No 5 9/2/90
CLLEND---------------------------------------------------------
C
C*L Arguments:-------------------------------------------------
      SUBROUTINE EXPAND32B(LENGTH, ARRAY, VERSION)

      IMPLICIT NONE

      INTEGER
     * LENGTH,       !IN length of the field to be expanded
     * VERSION       !IN model version

      REAL
     * ARRAY(LENGTH)  !IN/OUT array to be expanded in place

C -------------------------------------------------------------
C Local variables: --------------------------------------------
      REAL HOLD(LENGTH)     ! space for expanded array
      INTEGER I             ! Loop index
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
      EXTERNAL EXPAND21,P21BITS
      INTEGER  P21BITS


      CALL EXPAND21(LENGTH,ARRAY,HOLD,
     &              P21BITS(VERSION) )
      DO I=1,LENGTH
        ARRAY(I)=HOLD(I)
      ENDDO

      RETURN
      END
