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
CLL  SUBROUTINEs RH_TO_CC & CC_TO_RH------------------------------------
CLL
CLL  Purpose: calculate cloud cover (fraction) from input rh (%)
CLL         : and vice versa
CLL         : uses eqs P292.19 to P292.21 in UM Doc Paper 29
CLL
CLL  Model            Modification History :
CLL version  Date
CLL   3.3    22/12/93 Coded by Bruce Macpherson and Nigel Richards
CLL
CLL  Programming Standard : UM
CLL
CLL  Project Task : P29
CLL
CLLEND-------------------------------------------------------------
C
C----Arguments:---------------------------------------------------------
      SUBROUTINE RH_TO_CC  (RH,NPTS,RHC,CC)

      IMPLICIT NONE

      INTEGER
     + NPTS                ! IN No. of points on level.

      REAL
     + RHC                 ! IN Critical relative humidity (fraction).
     +,RH(NPTS)            ! IN Rel humidity (%).
     +,CC(NPTS)            ! OUT Cloud cover (fraction).

C----------------------------------------------------------------------
C     External subroutine calls NONE
C-----------------------------------------------------------------------

C Local variables------------------------------------------------------
      REAL WRH              ! Local rh (fraction)

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------


      REAL PC1,PC2,PC3                     ! local constants.
      PARAMETER (
     + PC1=1.060660172                     ! 3/sqrt(8).
     +,PC2=2.0*PC1                         !
     +,PC3=PI/3.0                          ! pi/3
     +)

      INTEGER I     ! Do loop index

        DO I=1,NPTS
C-----------------------------------------------------------------------
CLL Calculate cloud fraction.
C-----------------------------------------------------------------------
C         Work with rh fraction
          WRH=0.01*RH(I)
C         Remove any supersaturation
          IF(WRH.GT.1.0) WRH=1.0
          CC(I)=0.0
C         For WRH<RHC (including WRH<0), CC remains zero.
C         This treats the special MOPS rh=-85% for zero cloud cover.
          IF(WRH.GT.RHC .AND. WRH.LT.(5.+RHC)/6.)THEN
            CC(I)=2.*COS(PC3+ACOS( PC1*(WRH-RHC)/(1.-RHC) )/3.)
            CC(I)=CC(I)*CC(I)
          ENDIF
          IF(WRH.GE.(5.+RHC)/6.)THEN
            CC(I)=PC2*(1.-WRH)/(1.-RHC)
            CC(I)=1.-CC(I)**(2./3.)
          ENDIF
        ENDDO     ! end loop over points

      RETURN
      END

C----Arguments:---------------------------------------------------------
      SUBROUTINE CC_TO_RH  (CC,NPTS,RHC,RH)

      IMPLICIT NONE

      INTEGER
     + NPTS                ! IN No. of points on level.

      REAL
     + RHC                 ! IN Critical relative humidity (fraction).
     +,CC(NPTS)            ! IN Cloud cover (fraction).
     +,RH(NPTS)            ! OUT Rel humidity (%).

C----------------------------------------------------------------------
C     External subroutine calls NONE
C-----------------------------------------------------------------------

C Local variables------------------------------------------------------

      INTEGER I     ! Do loop index

C Code in calling routine restricts CC to range 0-1
C but check to be safe

       DO I=1,NPTS
        IF (CC(I).GT.1.0) CC(I) = 1.0
        IF (CC(I).LT.0.0) CC(I) = 0.0

        IF (CC(I).GT.0.5) THEN
C from eqn p292.21
         RH(I) = 1.0 -
     &           (1.-CC(I))**(3.0/2.0) * SQRT(2.0)/3.0 * (1.-RHC)
        ELSE
C from eqn p292.19
         RH(I) = RHC +
     &           SQRT(2.0*CC(I)) * (1.0-CC(I)/3.0) * (1.-RHC)
        ENDIF
        RH(I) = RH(I) * 100.0
       END DO

      RETURN
      END
