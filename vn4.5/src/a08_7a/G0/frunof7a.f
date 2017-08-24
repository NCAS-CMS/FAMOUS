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
CLL  SUBROUTINE FRUNOFF------------------------------------------------
CLL
CLL  PURPOSE : TO CALCULATE SURFACE RUNOFF
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
CLL  LOGICAL COMPONENTS : P252
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
CLL                  SECTION (3B(II)), EQN(P252.14)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE FRUNOFF (
     & NPNTS,TILE_PTS,TILE_INDEX,AREA,
     & CAN_CPY,CAN_WCNT,INFIL,R,FRAC,TIMESTEP,
     & SURF_ROFF
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
     &,CAN_WCNT(NPNTS)      ! IN Canopy water content (kg/m2).
     &,INFIL(NPNTS)         ! IN Infiltration rate (kg/m2/s).
     &,R(NPNTS)             ! IN Water fall rate (kg/m2/s).
     &,FRAC(NPNTS)          ! IN Tile fraction.
     &,TIMESTEP             ! IN Timestep (s).

      REAL
     & SURF_ROFF(NPNTS)     ! OUT Cummulative surface runoff (kg/m2/s).

!  Workspace --------------------------------------------------------
      REAL
     & AEXP                 ! Used in the calculation of exponential
     &,AEXP1                ! terms in the surface runoff formula.
     &,AEXP2                !
     &,CM                   ! (CAN_CPY - CAN_WCNT)/TIMESTEP
     &,CAN_RATIO            ! CAN_WCNT / CAN_CPY
     &,RUNOFF               ! Local runoff.

      INTEGER
     & I                    ! Land point index.
     &,J                    ! Counter for loop over tile points.

      DO J=1,TILE_PTS
        I = TILE_INDEX(J)
        RUNOFF = 0.
        IF (R(I).GT.0.0) THEN
          IF ( INFIL(I)*TIMESTEP.LE.CAN_WCNT(I)
     &                                   .AND. CAN_CPY(I).GT.0.0 ) THEN
! Infiltration in timestep < or = canopy water content
             AEXP = AREA*CAN_CPY(I)/R(I)
             IF (CAN_WCNT(I) .GT. 0.0) THEN
               AEXP1 = EXP( -AEXP*INFIL(I)/CAN_WCNT(I))
             ELSE
               AEXP1 = 0.0
             END IF
             AEXP2 = EXP( -AEXP/TIMESTEP)
             CAN_RATIO = CAN_WCNT(I)/CAN_CPY(I)
             CAN_RATIO = MIN(CAN_RATIO,1.0)
             RUNOFF = R(I) * ( CAN_RATIO*AEXP1 +
     *                                       (1. - CAN_RATIO)*AEXP2 )
!                                                        ... P252.14A
          ELSE
! Infiltration in timestep > canopy water content
             CM = (CAN_CPY(I)-CAN_WCNT(I))/TIMESTEP
             CM = MAX(CM,0.0)
             AEXP = EXP( -AREA*(INFIL(I)+CM)/R(I))
             RUNOFF = R(I)*AEXP                    !     ... P252.14B
          ENDIF
        ENDIF
        SURF_ROFF(I) = SURF_ROFF(I) + FRAC(I)*RUNOFF
      ENDDO

      RETURN
      END
