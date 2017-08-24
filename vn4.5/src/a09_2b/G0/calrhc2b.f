C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!
      SUBROUTINE RHCRIT_CALC(AK,BK,AKH,BKH,PSTAR,RHCPT,LEVELS,
     &      POINTS,PFIELD,T,Q,QCF,ROW_LENGTH,LAND,ICE_FRAC,BL_LEVELS)
!
      IMPLICIT NONE
!
!
!     Purpose: To calculate the critical relative humidity in every
!              grid-box.
!
!     Method : The critical relative humidity of a certain grid-box is
!            determined from the variance in a 3*3 set of boxes centred
!            on it. A fit, dependent on pressure, relates the variance
!            of the 3*3 region to the variance within the one grid-box.
!            This variance is converted to a critical relative humidity
!            in a straightforward fashion.
!            Some points in the 3*3 region may be excluded from the
!            variance calculation in the BL layers, if their
!            surfaces do not 'match'. The criterion for matching is that
!            land and sea-ice match, but that open ocean does not match
!            with either of these.
!            In all layers, points in the 3*3 region which lie outside 3
!            std devs of the mean are excluded in a second iteration of
!            the main calculation.
!            Notice that RHc is not the same at all points on the polar
!            row, which gives different values of T, q and qcl at the
!            polar row. However, these polar values are never used by
!            the physics, so this polar discrepancy does not affect
!            model evolution.
!
!     Compatible with version 2B of LS_CLD
!
!
! Current Owner of Code: S. Cusack
!
! History:
! Version   Date     Comment
!  4.5    15/09/98   Original Code     S. Cusack
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   System component covered: P292
!
!   Documentation:
!
!  Global Variables:----------------------------------------------------
!
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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

!
!  Subroutine arguments
!-----------------------------------------------------------------------
! IN variables
!-----------------------------------------------------------------------
      INTEGER LEVELS           ! No. of levels being processed.
!
      INTEGER BL_LEVELS        ! No. of BL levels being processed.
!
      INTEGER POINTS           ! No. of gridpoints being processed.
!
      INTEGER PFIELD           ! No. of points in global field (at one
!                                vertical level).
      INTEGER ROW_LENGTH
!
      REAL PSTAR(PFIELD)       ! Surface pressure (Pa).
      REAL ICE_FRAC(PFIELD)    ! Ice fraction.
!
      REAL AK(LEVELS)          ! Hybrid "A" co-ordinate.
      REAL BK(LEVELS)          ! Hybrid "B" co-ordinate.
      REAL AKH(LEVELS+1)       ! Hybrid "A" co-ordinate.
      REAL BKH(LEVELS+1)       ! Hybrid "B" co-ordinate.
      REAL Q(PFIELD,LEVELS)    ! Total water content (Q+QCL,kg per kg)
      REAL QCF(PFIELD,LEVELS)  ! Ice water content (kg per kg air)
      REAL T(PFIELD,LEVELS)    ! Liquid/frozen water temperature (TL,K)
      LOGICAL LAND(PFIELD)     ! The model land mask
!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
      REAL RHCPT(PFIELD,LEVELS)  ! Critical relative humidity at every
!                                  grid cell.
!
!  Local Parameters
!
!     Comdeck for use with the RHcrit parametrization of the large-scale
!     cloud scheme, A09_2B.
!     Four constants are used to specify the variable (a function of
!     pressure) which relates the variability of the saturation variable
!     in one box to the variability over 9 climate grid boxes. The
!     variability of the saturation variable in one box is required to
!     specify RHcrit.
!     Note that the constants
!       RHC_CON1=0.522, RHC_CON2=0.122, RHC_CON3=2.5E3, RHC_CON4=1.75E4
!     are only suitable for use in the 2.5*3.75 degrees climate model:
!     these constants depend upon the size of a grid-box.
!
!     The fit is of the form:
!          A=RHC_CON1+RHC_CON2*(p-RHC_CON4)/(RHC_CON3+abs(p-RHC_CON4))
!       where p is the pressure at the layer midpoint.
!     Then,
!           sigma(s) = A * sigma(s,9)
!        where sigma(s) is the std dev of the saturation variable s in
!     one grid-box, and sigma(s,9) is the std dev over 9 boxes.
!
!     RHC_MIN and RHC_MAX are user defined limits on the values of RHc.
!
!     S. Cusack   02-09-98
!
      REAL RHC_CON1, RHC_CON2, RHC_CON3, RHC_CON4, RHC_MIN, RHC_MAX
!
      PARAMETER(RHC_CON1 = 0.522
     &        , RHC_CON2 = 0.122
     &        , RHC_CON3 = 2.5E3
     &        , RHC_CON4 = 1.75E4
     &        , RHC_MIN = 0.3
     &        , RHC_MAX = 0.98)

!
      REAL INV9, LS, LSRCP, ERCPR
      PARAMETER ( INV9 = 1./9.
     &          , LS = LC+LF
     &          , LSRCP = (LC+LF)/CP
     &          , ERCPR = EPSILON/(CP*R))
!
!  Local Variables
!
      REAL
     &    MEAN_SUPSAT       ! MEAN RH OF 3*3 REGION
     &   ,P_LEV(POINTS,LEVELS)! Pressure at model levels
     &   ,TL(POINTS,LEVELS) ! Conserved temperature (P292.1, UMDP29)
     &   ,QT(POINTS,LEVELS) ! Conserved WATER(P292.2, UMDP29)
     & ,  TOT_VAR           ! TOTAL VARIANCE OF 3*3 REGION
     & ,  QST(POINTS)       ! SATURATION VAPOUR PRESSURE IN GRID-BOX
     & ,  SUPSAT(POINTS,LEVELS)! 'RELATIVE HUMIDITY' OF GRID-BOX
     & ,  SUPSAT_SD_1       ! STANDARD DEVIATION OF 'R.H.' IN GRID-BOX
     & ,  SUPSAT_SD_3       ! RESOLVED STD DEV OF 'R.H.' IN 3*3 REGION
     & ,  P_GRAD(POINTS,LEVELS)! CONSTANT WHICH RELATES RH_SD_3 TO
!                                 RH_SD_1
     & ,  ROOT_6            ! =sqrt(6.)
     & ,  LATHT             ! =Lc if T>Tm, ELSE = Lc+Lf
     & ,  AL(POINTS)        ! Variable defined in P292.6 in UMDP 29
     & ,  SURF_MULT(POINTS,8)! A multiplier to take into account
!                              surface matching.
     & ,  THREE_SIGMA       ! Three times sigma
     & ,  SUPSAT1           ! Temporary variable
     & ,  SUPSAT2           ! Temporary variable
     & ,  SUPSAT3           ! Temporary variable
     & ,  SUPSAT4           ! Temporary variable
     & ,  SUPSAT5           ! Temporary variable
     & ,  SUPSAT6           ! Temporary variable
     & ,  SUPSAT7           ! Temporary variable
     & ,  SUPSAT8           ! Temporary variable
!
      INTEGER
     &    I, K, J           ! Simple loop variables
     & ,  ICOUNT(POINTS)    ! Counter of points
     & ,  COUNT             ! Counter of points
     & ,  RL, RLM1, RLP1    !
!
      LOGICAL
     &    OCEAN(POINTS)     ! Those points which are not land, and have
!                             a sea-ice fraction less than 0.25.
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT
!- End of Header
!
!
      ROOT_6=SQRT(6.)
      RL=ROW_LENGTH
      RLP1=ROW_LENGTH+1
      RLM1=ROW_LENGTH-1
!
      DO K=1,LEVELS
        DO I=1,POINTS
          P_LEV(I,K)=AK(K)+PSTAR(I)*BK(K)
          P_GRAD(I,K)=RHC_CON1+RHC_CON2*(P_LEV(I,K)-RHC_CON4)/
     &                             (RHC_CON3+ABS(P_LEV(I,K)-RHC_CON4))
! Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
! (Assumes version 3A onwards of Section 4)
          TL(I,K)=T(I,K)-LSRCP*QCF(I,K)
          QT(I,K)=Q(I,K)+QCF(I,K)
        ENDDO
      ENDDO
!
! Ocean points defined now as not land and where ice fraction LT 0.25
      DO I=1,POINTS
        OCEAN(I)=(.NOT.LAND(I)).AND.(ICE_FRAC(I).LT.2.5E-1)
      ENDDO
!
! A real no. is now assigned to every neighbouring point of every point
! on the grid, if their surfaces match it has the value one, else it is
! zero.
      DO J=1,8
        DO I=(ROW_LENGTH+2),(POINTS-ROW_LENGTH)
          SURF_MULT(I,J)=0.
        ENDDO
      ENDDO
      DO I=(ROW_LENGTH+2),(POINTS-ROW_LENGTH)
        ICOUNT(I)=1
        IF ((OCEAN(I).AND.OCEAN(I-RLP1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I-RLP1))) THEN
          SURF_MULT(I,1)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I-RL)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I-RL))) THEN
          SURF_MULT(I,2)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I-RLM1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I-RLM1))) THEN
          SURF_MULT(I,3)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I-1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I-1))) THEN
          SURF_MULT(I,4)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I+1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I+1))) THEN
          SURF_MULT(I,5)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I+RLM1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I+RLM1))) THEN
          SURF_MULT(I,6)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I+RL)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I+RL))) THEN
          SURF_MULT(I,7)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
        IF ((OCEAN(I).AND.OCEAN(I+RLP1)).OR.
     &               (.NOT.OCEAN(I).AND..NOT.OCEAN(I+RLP1))) THEN
          SURF_MULT(I,8)=1.
          ICOUNT(I)=ICOUNT(I)+1
        ENDIF
      ENDDO
!
! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.
!
      DO K=1,BL_LEVELS
        CALL QSAT(QST,TL(1,K),P_LEV(1,K),POINTS)
        DO I=1,POINTS
          IF (TL(I,K).GT.TM) THEN
            LATHT=LC/TL(I,K)
          ELSE
            LATHT=LS/TL(I,K)
          ENDIF
          AL(I)=1./(1.+LATHT*LATHT*ERCPR*QST(I))
! SUPSAT given by P292.3 of UMDP 29.
          SUPSAT(I,K)=AL(I)*(QT(I,K)-QST(I))
        ENDDO
        DO I=(ROW_LENGTH+2),(POINTS-ROW_LENGTH-1)
          SUPSAT1=SURF_MULT(I,1)*SUPSAT(I-RLP1,K)
          SUPSAT2=SURF_MULT(I,2)*SUPSAT(I-RL,K)
          SUPSAT3=SURF_MULT(I,3)*SUPSAT(I-RLM1,K)
          SUPSAT4=SURF_MULT(I,4)*SUPSAT(I-1,K)
          SUPSAT5=SURF_MULT(I,5)*SUPSAT(I+1,K)
          SUPSAT6=SURF_MULT(I,6)*SUPSAT(I+RLM1,K)
          SUPSAT7=SURF_MULT(I,7)*SUPSAT(I+RL,K)
          SUPSAT8=SURF_MULT(I,8)*SUPSAT(I+RLP1,K)
          MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4
     &               + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8
     &               + SUPSAT(I,K)) / ICOUNT(I)
          TOT_VAR=SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2
     &          + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4
     &          + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6
     &          + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8
     &          + SUPSAT(I,K)*SUPSAT(I,K)
     &          - ICOUNT(I)*MEAN_SUPSAT*MEAN_SUPSAT
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 3*sigma of the mean are considered outliers and are
!  rejected.
          IF (ICOUNT(I).GT.1) THEN
            THREE_SIGMA=3.*SQRT(TOT_VAR/ICOUNT(I))
          ELSE
            THREE_SIGMA=QST(I)*0.01
          ENDIF
          COUNT=1
          IF (ABS(SUPSAT(I-RLP1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT1=0.
          ELSE IF (SURF_MULT(I,1).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I-RL,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT2=0.
          ELSE IF (SURF_MULT(I,2).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I-RLM1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT3=0.
          ELSE IF (SURF_MULT(I,3).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I-1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT4=0.
          ELSE IF (SURF_MULT(I,4).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I+1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT5=0.
          ELSE IF (SURF_MULT(I,5).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I+RLM1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT6=0.
          ELSE IF (SURF_MULT(I,6).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I+RL,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT7=0.
          ELSE IF (SURF_MULT(I,7).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (ABS(SUPSAT(I+RLP1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
            SUPSAT8=0.
          ELSE IF (SURF_MULT(I,8).GT.0.5) THEN
            COUNT=COUNT+1
          ENDIF
          IF (COUNT.GT.1) THEN
            MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4
     &               + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8
     &               + SUPSAT(I,K)) / COUNT
            TOT_VAR=SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2
     &          + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4
     &          + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6
     &          + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8
     &          + SUPSAT(I,K)*SUPSAT(I,K)
     &          - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
            SUPSAT_SD_3=SQRT(TOT_VAR/COUNT)
          ELSE
            SUPSAT_SD_3=QST(I)*0.01
          ENDIF
!
! P_GRAD determines the relation between 3*3 and sub-grid variance.
          SUPSAT_SD_1=P_GRAD(I,K)*SUPSAT_SD_3
! RHCPT defined from P292.14 in UMDP 29
          RHCPT(I,K)=1.-ROOT_6*SUPSAT_SD_1/(AL(I)*QST(I))
! RHcrit is now limited to lie between a range defined in RHCCON2B
          RHCPT(I,K)=MAX(RHCPT(I,K),RHC_MIN)
          RHCPT(I,K)=MIN(RHCPT(I,K),RHC_MAX)
        ENDDO
! North and south haloes not determined above, now they are filled in.
! Call to SWAPBOUNDS after this subroutine fills the haloes in properly.
        RHCPT(ROW_LENGTH+1,K)=RHCPT(ROW_LENGTH+2,K)
        DO I=1,ROW_LENGTH
          RHCPT(I,K)=RHCPT(I+ROW_LENGTH,K)
        ENDDO
        RHCPT(POINTS-ROW_LENGTH,K)=RHCPT(POINTS-ROW_LENGTH-1,K)
        DO I=(POINTS-ROW_LENGTH+1),POINTS
          RHCPT(I,K)=RHCPT(I-ROW_LENGTH,K)
        ENDDO
      ENDDO
!
! The same calculations as above are performed, but the 'surface match'
! criterion is now dropped (atmosphere less influenced by surface at
! greater heights).
      IF (LEVELS.GT.BL_LEVELS) THEN
!
        DO K=(BL_LEVELS+1),LEVELS
          CALL QSAT(QST,TL(1,K),P_LEV(1,K),POINTS)
          DO I=1,POINTS
            IF (TL(I,K).GT.TM) THEN
              LATHT=LC/TL(I,K)
            ELSE
              LATHT=LS/TL(I,K)
            ENDIF
            AL(I)=1./(1.+LATHT*LATHT*ERCPR*QST(I))
            SUPSAT(I,K)=AL(I)*(QT(I,K)-QST(I))
          ENDDO
          DO I=(ROW_LENGTH+2),(POINTS-ROW_LENGTH-1)
            SUPSAT1=SUPSAT(I-RLP1,K)
            SUPSAT2=SUPSAT(I-RL,K)
            SUPSAT3=SUPSAT(I-RLM1,K)
            SUPSAT4=SUPSAT(I-1,K)
            SUPSAT5=SUPSAT(I+1,K)
            SUPSAT6=SUPSAT(I+RLM1,K)
            SUPSAT7=SUPSAT(I+RL,K)
            SUPSAT8=SUPSAT(I+RLP1,K)
            MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4
     &               + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8
     &               + SUPSAT(I,K)) * INV9
            TOT_VAR=SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2
     &          + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4
     &          + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6
     &          + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8
     &          + SUPSAT(I,K)*SUPSAT(I,K)
     &          - 9. * MEAN_SUPSAT*MEAN_SUPSAT
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 3*sigma of the mean are considered outliers and are
!  rejected.
            THREE_SIGMA=3.*SQRT(TOT_VAR*INV9)
            COUNT=1
            IF (ABS(SUPSAT(I-RLP1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT1=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I-RL,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT2=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I-RLM1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT3=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I-1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT4=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I+1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT5=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I+RLM1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT6=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I+RL,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT7=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (ABS(SUPSAT(I+RLP1,K)-MEAN_SUPSAT).GT.THREE_SIGMA) THEN
              SUPSAT8=0.
            ELSE
              COUNT=COUNT+1
            ENDIF
            IF (COUNT.GT.1) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4
     &               + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8
     &               + SUPSAT(I,K)) / COUNT
              TOT_VAR=SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2
     &          + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4
     &          + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6
     &          + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8
     &          + SUPSAT(I,K)*SUPSAT(I,K)
     &          - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              SUPSAT_SD_3=SQRT(TOT_VAR/COUNT)
            ELSE
              SUPSAT_SD_3=QST(I)*0.01
            ENDIF
!
            SUPSAT_SD_1=P_GRAD(I,K)*SUPSAT_SD_3
            RHCPT(I,K)=1.-ROOT_6*SUPSAT_SD_1/(AL(I)*QST(I))
! RHcrit is now limited to lie between a range defined in RHCCON2B
            RHCPT(I,K)=MAX(RHCPT(I,K),RHC_MIN)
            RHCPT(I,K)=MIN(RHCPT(I,K),RHC_MAX)
          ENDDO
! North and south haloes not determined above, now they are filled in.
! Call to SWAPBOUNDS after this subroutine fills the haloes in properly.
          RHCPT(ROW_LENGTH+1,K)=RHCPT(ROW_LENGTH+2,K)
          DO I=1,ROW_LENGTH
            RHCPT(I,K)=RHCPT(I+ROW_LENGTH,K)
          ENDDO
          RHCPT(POINTS-ROW_LENGTH,K)=RHCPT(POINTS-ROW_LENGTH-1,K)
          DO I=(POINTS-ROW_LENGTH+1),POINTS
            RHCPT(I,K)=RHCPT(I-ROW_LENGTH,K)
          ENDDO
        ENDDO
      ENDIF  !  LEVELS GT BL_LEVELS
!
      RETURN
      END
