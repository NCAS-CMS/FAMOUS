C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE SIEVE--------------------------------------------------
CLL
CLL  PURPOSE : TO CALCULATE THE THROUGHFALL OF WATER FALLING
CLL            THROUGH THE SURFACE CANOPY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  WRITTEN FOR CRAY-YMP BY S.ALLEN AND D.GREGORY
CLL  DEC - FEB 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1 18/1/90
CLL
CLL  LOGICAL COMPONENTS COVERED: P252
CLL
CLL  SYSTEM TASK : P252
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
CCL                  SECTION (3B(II)), EQN(P252.9)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE SIEVE (
     & NPNTS,TILE_PTS,TILE_INDEX,AREA,CAN_CPY,R,FRAC,TIMESTEP,
     & CAN_WCNT,TOT_TFALL
     & )

      IMPLICIT NONE

      INTEGER
     & NPNTS                ! IN Total number of land points.
     &,TILE_PTS             ! IN Number of tile points.
     &,TILE_INDEX(NPNTS)    ! IN Index of tile points.

      REAL
     & AREA                 ! IN Fractional area of gridbox over which
                            !    water falls (%).
     &,CAN_CPY(NPNTS)       ! IN Canopy capacity (kg/m2).
     &,R(NPNTS)             ! IN Water fall rate (kg/m2/s).
     &,FRAC(NPNTS)          ! IN Tile fraction.
     &,TIMESTEP             ! IN Timestep (s).

      REAL
     & CAN_WCNT(NPNTS)      ! INOUT Canopy water content (kg/m2).
     &,TOT_TFALL(NPNTS)     ! INOUT Cummulative canopy throughfall
!                           !       (kg/m2/s).

!  Workspace --------------------------------------------------------
      REAL
     & AEXP                 ! Used in calculation of exponential
                            ! in throughfall formula.
     &,CAN_RATIO            ! CAN_WCNT / CAN_CPY
     &,TFALL(NPNTS)         ! Local throughfall (kg/m2/s).

      INTEGER
     & I                    ! Land point index.
     &,J                    ! Counter for loop over tile points.

      DO J=1,TILE_PTS
        I = TILE_INDEX(J)
        IF (CAN_CPY(I) .GT. 0.0 .AND. R(I) .GT. 0.0) THEN
           AEXP = AREA*CAN_CPY(I)/(R(I)*TIMESTEP)
           AEXP = EXP(-AEXP)
           CAN_RATIO = CAN_WCNT(I) / CAN_CPY(I)
           CAN_RATIO = MIN(CAN_RATIO,1.0)
           TFALL(I) = R(I) * ((1.0-CAN_RATIO)*AEXP + CAN_RATIO)
        ELSE
           TFALL(I) = R(I)
        END IF
        CAN_WCNT(I) = CAN_WCNT(I) + (R(I) - TFALL(I))*TIMESTEP
        TOT_TFALL(I) = TOT_TFALL(I) + FRAC(I)*TFALL(I)
      ENDDO

      RETURN
      END
