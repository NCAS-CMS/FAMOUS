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
CLL  SUBROUTINE SNOWPR -------------------------------------------------
CLL
CLL  Purpose: This routine calculates the snowfall probabilty
CLL           using an equation based on the 1000-850 TT
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL J.Heming    <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
!LL   4.5    15/04/98 Start-end args added to enable dupicate halo
!LL                   calculations to be avoided. S.D.Mullerworth
CLL
CLL  Logical components covered D432
CLL  Project TASK: D4
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation
CLL
CLLEND------------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE SNOWPR
     1 (P,PSTAR,P_EXNER_HALF,THETA,Q,MODEL_HALF_HEIGHT,
     2  P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,AKH,BKH,
     &  SNPROB,START,END)
C*
C*L---------------------------------------------------------------------
      IMPLICIT NONE
C*
C*L---------------------------------------------------------------------
      INTEGER
     *  P_FIELD          ! IN  No of points in field
     *, P_LEVELS         ! IN  No of pressure levels
     *, Q_LEVELS         ! IN  No of wet levels
     *, Z_REF            ! IN  Level of model used to calculate PMSL
     &, START,END        ! IN  Range of points to calculate
C-----------------------------------------------------------------------
      REAL
     *  P(P_FIELD,P_LEVELS)          ! IN  Pressure array
     *, PSTAR(P_FIELD)               ! IN  Pressure on surface of earth
     *, P_EXNER_HALF(P_FIELD,P_LEVELS+1) !IN Exner press on half levels
     *, THETA(P_FIELD,P_LEVELS)      !IN  Potential temperature
     *, Q(P_FIELD,Q_LEVELS)          !IN  Specific Humidity
     *, MODEL_HALF_HEIGHT(P_FIELD,P_LEVELS+1) !IN Heights on half levels
     *, AKH(P_LEVELS+1)              !IN  A values on half levels
     *, BKH(P_LEVELS+1)              !IN  B values on half levels
     *, SNPROB(P_FIELD)              ! OUT Snow probability in %
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
      REAL
     *  PRESSURE(P_FIELD)  ! Pressure at which height is calculated
     *, Z_100000(P_FIELD)  ! Height at 1000000 Pa
     *, Z_85000(P_FIELD)   ! Height at 850000 Pa
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL V_INT_Z
C*----------------------------------------------------------------------
C
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C
C*L---------------------------------------------------------------------
C   DEFINE LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  I     ! Loop counter
     *, EDGE  ! Highest height at which snow prob will be calculated
C*----------------------------------------------------------------------
CL  1. Calculate height of surface at pressure of 100000Pa
C-----------------------------------------------------------------------
      DO I=START,END
        PRESSURE(I)=100000.
      ENDDO
C-----------------------------------------------------------------------
      CALL V_INT_Z(PRESSURE,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,
     & MODEL_HALF_HEIGHT,Z_100000,P_FIELD,P_LEVELS,Q_LEVELS,
     & Z_REF,AKH,BKH,START,END)
C-----------------------------------------------------------------------
CL  2. Calculate height of surface at pressure of 85000Pa
C-----------------------------------------------------------------------
      DO I=START,END
        PRESSURE(I)=85000.
      ENDDO
C-----------------------------------------------------------------------
      CALL V_INT_Z(PRESSURE,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,
     & MODEL_HALF_HEIGHT,Z_85000,P_FIELD,P_LEVELS,Q_LEVELS,
     & Z_REF,AKH,BKH,START,END)
C-----------------------------------------------------------------------
CL  3. Calculate the snow probability in %
C-----------------------------------------------------------------------
      EDGE=1.E8
      DO I=START,END
        IF (Z_100000(I).LE.EDGE) THEN
          SNPROB(I)=4.*(1305.-(Z_85000(I)-.96666*Z_100000(I)))
        ELSE
          SNPROB(I)=0.
        ENDIF
        IF (SNPROB(I).LE.0.) SNPROB(I)=0.
        IF (SNPROB(I).GE.100.) SNPROB(I)=100.
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
