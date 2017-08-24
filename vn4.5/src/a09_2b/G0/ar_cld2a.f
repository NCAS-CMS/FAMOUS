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
      SUBROUTINE AREA_CLD(
     & AK,BK,PSTAR,RHCRIT,LEVELS,RHCPT,POINTS,PFIELD,T,CF_BULK,Q,QCF,
     & QCL,CF_LIQ,CF_ICE,ERROR,CF_AREA,AKH,BKH)
!
      IMPLICIT NONE
!
!
!
!     Purpose: To calculate an area as well as volume fraction of
!              cloud in a layer.
!
!     Method : A vertical profile of (q+qcL), qcF and Tl within a layer
!              is created, from which mean (q+qcL), qcF and Tl values
!              for 3 sub-layers can be found.  These 3 sub-layers are
!              one-third of the total layer thickness.
!              Cloud volume calculations are performed on each sub-
!              layer, and the maximum cloud volume is the cloud area.
!
!     Comments: Compatible with versions 2A and 2B of Section 9.
!
!
! Current Owner of Code: S. Cusack
!
! History:
! Version   Date     Comment
!  4.5    14/05/98   Original Code     S. Cusack
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
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
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

!
!   Subroutine arguments
!----------------------------------------------------------------------
! IN variables
!----------------------------------------------------------------------
      INTEGER LEVELS           ! No. of levels being processed.
!
      INTEGER POINTS           ! No. of gridpoints being processed.
!
      INTEGER PFIELD           ! No. of points in global field (at one
!                                vertical level).
!
      REAL PSTAR(PFIELD)       ! Surface pressure (Pa).
!
      REAL RHCRIT(LEVELS)      ! Critical relative humidity.  See the
!                                the paragraph incorporating eqs P292.11
!                                to P292.14; the values need to be tuned
!                                for the given set of levels.
      REAL QCF(PFIELD,LEVELS)  ! Cloud ice content at processed levels
!                                  (kg per kg air).
      REAL AK(LEVELS)          ! Hybrid "A" co-ordinate.
      REAL BK(LEVELS)          ! Hybrid "B" co-ordinate.
      REAL AKH(LEVELS+1)          ! Hybrid "A" co-ordinate.
      REAL BKH(LEVELS+1)          ! Hybrid "B" co-ordinate.
!-----------------------------------------------------------------------
! INOUT variables
!-----------------------------------------------------------------------
      REAL Q(PFIELD,LEVELS)    ! On input:  Total water content (QW)
!                                           (kg per kg air).
!                                On output: Specific humidity at process
!                                           levels (kg water per kg air)
      REAL T(PFIELD,LEVELS)    ! On input:  Liquid/frozen water
!                                           temperature (TL) (K).
!                                On output: Temperature at processed
!                                           levels (K).
!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
      REAL CF_BULK(PFIELD,LEVELS)! Cloud fraction at processed levels
!                                  (decimal fraction).
      REAL QCL(PFIELD,LEVELS)    ! Cloud liquid water content at
!                                  processed levels (kg per kg air).
      REAL CF_LIQ(PFIELD,LEVELS) ! Grid-box mean cloud
!                                            condensate at processed
!                                            levels (kg per kg air).
      REAL CF_ICE(PFIELD,LEVELS) ! Max moisture fluctuation
!                                          /6*sigma at processed levels
!                                            (kg per kg air).
      REAL CF_AREA(PFIELD,LEVELS)! Area cloud fraction
!                                      (decimal fraction).
      REAL RHCPT(PFIELD,LEVELS)  ! Critical relative humidity for
!                                  every grid cell.
      INTEGER ERROR              ! 0 if OK; 1 if bad arguments.
!
!
! Local parameters
!
      REAL LCRCP, ONETHIRD, DRAT_THRESH
!
      PARAMETER (
     &            LCRCP=LC/CP
!      Latent heat of condensation divided by specific heat at constant
!      pressure
     &          , ONETHIRD=1./3.
!      Threshold for change in (QT-QSAT)/QSAT between 2 layer midpoints
     &          , DRAT_THRESH=3.0E-01
     &           )
!
!
!  Local Variables
!
      REAL
     &     CF_TEMP(PFIELD,LEVELS)
!      Temporary storage of total cloud fraction in 3 sub-layers
     &    ,CF_LIQ_TEMP(PFIELD,LEVELS)
!      Temporary storage of total liquid cloud fraction in 3 sub-layers
     &    ,CF_ICE_TEMP(PFIELD,LEVELS)
!      Temporary storage of total ice cloud fraction in 3 sub-layers
     &    ,QCL_TEMP(PFIELD,LEVELS)
!      Temporary storage of total cloud frozen qc in 3 sub-layers
     &    ,P_LEV(POINTS,LEVELS)
!      Pressure at model levels
     &    ,P_HLEV(POINTS,LEVELS)
!      Pressure at model half-levels
     &    ,T_MID(PFIELD,LEVELS)
!      Mean temperature in middle third segment of layer
     &    ,T_UP(PFIELD,LEVELS)
!      Mean temperature in upper third segment of layer
     &    ,T_LOW(PFIELD,LEVELS)
!      Mean temperature in lower third segment of layer
     &    ,QT_UP(PFIELD,LEVELS)
!      Mean total water in upper third segment of layer
     &    ,QT_LOW(PFIELD,LEVELS)
!      Mean total water in lower third segment of layer
     &    ,QT_MID(PFIELD,LEVELS)
!      Mean total water in middle third segment of layer
     &    ,QCF_UP(PFIELD,LEVELS)
!      Mean total water in upper third segment of layer
     &    ,QCF_LOW(PFIELD,LEVELS)
!      Mean total water in lower third segment of layer
     &    ,QCF_MID(PFIELD,LEVELS)
!      Mean total water in middle third segment of layer
     &    ,QST(POINTS)
!      Saturated mixing ratio at layer temperature
     &    ,DTDP(POINTS,LEVELS)
!      Gradient of temperature w.r.t. pressure
     &    ,DQTDP(POINTS,LEVELS)
!      Gradient of (Q+QCL) w.r.t. pressure
     &    ,DQCFDP(POINTS,LEVELS)
!      Gradient of QCF w.r.t. pressure
     &    ,DRAT(POINTS,LEVELS)
!      Difference in (QT-QSAT)/QSAT between 2 layer midpoints
     &    ,T_LHLEV(POINTS,LEVELS)
!      Temperature at lower half-level in a layer
     &    ,T_UHLEV(POINTS,LEVELS)
!      Temperature at upper half-level in a layer
     &    ,QT_LHLEV(POINTS,LEVELS)
!      (Q+QCL) at lower half-level in a layer
     &    ,QT_UHLEV(POINTS,LEVELS)
!      (Q+QCL) at upper half-level in a layer
     &    ,QCF_LHLEV(POINTS,LEVELS)
!      QCF at lower half-level in a layer
     &    ,QCF_UHLEV(POINTS,LEVELS)
!      QCF at upper half-level in a layer
     &    ,DPLH
!      Difference in pressure between layer upper half-level and layer
!      mean.
     &    ,NUM1
!      Temporary variable
     &    ,T_LEV
!      Temperature at mid-layer calculated so that vertical profile of
!      T in a layer is conservative.
     &    ,QT_LEV
!      (Q+QCL) at mid-layer calculated so that vertical profile of
!      (Q+QCL) in a layer is conservative.
     &    ,QCF_LEV
!      QCF at mid-layer calculated so that vertical profile of QCF in
!      a layer is conservative.
!
      INTEGER
     &     I, LEVEL   ! HORIZ AND VERT LOOP VARIABLES RESPECTIVELY
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL GLUE_CLD, QSAT
!- End of Header
!
!
!  Calculate pressures at model levels and half-levels.
      DO LEVEL=1,LEVELS
        DO I=1,POINTS
          P_LEV(I,LEVEL)=AK(LEVEL)+PSTAR(I)*BK(LEVEL)
          P_HLEV(I,LEVEL)=AKH(LEVEL)+PSTAR(I)*BKH(LEVEL)
        ENDDO
      ENDDO
!
!   Calculate gradients of T, QCF and (Q+QCL) w.r.t. pressure
      DO LEVEL=2,LEVELS
        DO I=1,POINTS
          num1=1./(P_LEV(I,LEVEL-1)-P_LEV(I,LEVEL))
          DTDP(I,LEVEL)=(T(I,LEVEL-1)-T(I,LEVEL))*num1
          DQTDP(I,LEVEL)=(Q(I,LEVEL-1)-Q(I,LEVEL))*num1
          DQCFDP(I,LEVEL)=(QCF(I,LEVEL-1)-QCF(I,LEVEL))*num1
        ENDDO
      ENDDO
      DO I=1,POINTS
        DTDP(I,1)=DTDP(I,2)
        DQTDP(I,1)=DQTDP(I,2)
        DQCFDP(I,1)=DQCFDP(I,2)
      ENDDO
!  Calculate the gradient w.r.t. levels of (QT-QSAT)/QSAT.
      LEVEL=1
      CALL QSAT(QST,T(1,LEVEL),P_LEV(1,LEVEL),POINTS)
      DO LEVEL=2,(LEVELS-1)
        DO I=1,POINTS
          DRAT(I,LEVEL)=(Q(I,LEVEL-1)+QCF(I,LEVEL-1)-QST(I))/QST(I)
        ENDDO
        CALL QSAT(QST,T(1,LEVEL),P_LEV(1,LEVEL),POINTS)
        DO I=1,POINTS
          DRAT(I,LEVEL)=ABS(DRAT(I,LEVEL)-(Q(I,LEVEL)+QCF(I,LEVEL)-
     &                                                 QST(I))/QST(I))
        ENDDO
      ENDDO
!
!  From values at model levels and gradients, calculate values at half-
!  levels.
      DO LEVEL=2,LEVELS
        DO I=1,POINTS
          DPLH=P_HLEV(I,LEVEL)-P_LEV(I,LEVEL)
          T_LHLEV(I,LEVEL)=T(I,LEVEL)+DTDP(I,LEVEL)*DPLH
          QT_LHLEV(I,LEVEL)=Q(I,LEVEL)+DQTDP(I,LEVEL)*DPLH
          QCF_LHLEV(I,LEVEL)=QCF(I,LEVEL)+DQCFDP(I,LEVEL)*DPLH
          T_UHLEV(I,LEVEL-1)=T_LHLEV(I,LEVEL)
          QT_UHLEV(I,LEVEL-1)=QT_LHLEV(I,LEVEL)
          QCF_UHLEV(I,LEVEL-1)=QCF_LHLEV(I,LEVEL)
        ENDDO
      ENDDO
!
!  Now determine whether gradients in (QT-QSAT)/QSAT exceed a threshold,
!  in which case switch to the other vertical interpolation method.
      DO LEVEL=2,(LEVELS-1)
        DO I=1,POINTS
          IF ((DRAT(I,LEVEL).GT.DRAT_THRESH)) THEN
            num1=P_HLEV(I,LEVEL)-P_LEV(I,LEVEL)
            T_LHLEV(I,LEVEL)=T(I,LEVEL)+DTDP(I,LEVEL+1)*num1
            QT_LHLEV(I,LEVEL)=Q(I,LEVEL)+DQTDP(I,LEVEL+1)*num1
            QCF_LHLEV(I,LEVEL)=QCF(I,LEVEL)+DQCFDP(I,LEVEL+1)*num1
            num1=P_LEV(I,LEVEL-1)-P_HLEV(I,LEVEL)
            T_UHLEV(I,LEVEL-1)=T(I,LEVEL-1)-DTDP(I,LEVEL-1)*num1
            QT_UHLEV(I,LEVEL-1)=Q(I,LEVEL-1)-DQTDP(I,LEVEL-1)*num1
            QCF_UHLEV(I,LEVEL-1)=QCF(I,LEVEL-1)-DQCFDP(I,LEVEL-1)*num1
          ENDIF
        ENDDO
      ENDDO
!
!  If switch is activated above and below a layer, do no interpolation
!  at all. Values at upper and lower half levels equal layer mean value.
      DO LEVEL=2,(LEVELS-2)
        DO I=1,POINTS
          IF ((DRAT(I,LEVEL).GT.DRAT_THRESH).AND.
     &              (DRAT(I,LEVEL+1).GT.DRAT_THRESH)) THEN
            T_LHLEV(I,LEVEL)=T(I,LEVEL)
            T_UHLEV(I,LEVEL)=T(I,LEVEL)
            QT_LHLEV(I,LEVEL)=Q(I,LEVEL)
            QT_UHLEV(I,LEVEL)=Q(I,LEVEL)
            QCF_LHLEV(I,LEVEL)=QCF(I,LEVEL)
            QCF_UHLEV(I,LEVEL)=QCF(I,LEVEL)
          ENDIF
        ENDDO
      ENDDO
!
      DO LEVEL=2,(LEVELS-1)
        DO I=1,POINTS
! Re-calculate quantities at midpoint of layer, so as to ensure
! conservation.
          T_LEV=T(I,LEVEL)+T(I,LEVEL)-
     &        0.5*(T_LHLEV(I,LEVEL)+T_UHLEV(I,LEVEL))
          QT_LEV=Q(I,LEVEL)+Q(I,LEVEL)-
     &        0.5*(QT_LHLEV(I,LEVEL)+QT_UHLEV(I,LEVEL))
          QCF_LEV=QCF(I,LEVEL)+QCF(I,LEVEL)-
     &        0.5*(QCF_LHLEV(I,LEVEL)+QCF_UHLEV(I,LEVEL))
! Calculate the mean value of each quantity in the 3 sub-layers
          T_LOW(I,LEVEL)=ONETHIRD*(T_LEV-T_LHLEV(I,LEVEL))+
     &                                           T_LHLEV(I,LEVEL)
          T_UP(I,LEVEL)=0.666*(T_UHLEV(I,LEVEL)-T_LEV)+T_LEV
          T_MID(I,LEVEL)=3.*T(I,LEVEL)-T_LOW(I,LEVEL)-
     &                                           T_UP(I,LEVEL)
          QT_LOW(I,LEVEL)=ONETHIRD*(QT_LEV-QT_LHLEV(I,LEVEL))+
     &                                           QT_LHLEV(I,LEVEL)
          QT_UP(I,LEVEL)=0.666*(QT_UHLEV(I,LEVEL)-QT_LEV)+QT_LEV
          QT_MID(I,LEVEL)=3.*Q(I,LEVEL)-QT_LOW(I,LEVEL)-
     &                                           QT_UP(I,LEVEL)
          QCF_LOW(I,LEVEL)=ONETHIRD*(QCF_LEV-QCF_LHLEV(I,LEVEL))+
     &                                           QCF_LHLEV(I,LEVEL)
          QCF_UP(I,LEVEL)=0.666*(QCF_UHLEV(I,LEVEL)-QCF_LEV)+QCF_LEV
          QCF_MID(I,LEVEL)=3.*QCF(I,LEVEL)-QCF_LOW(I,LEVEL)-
     &                                           QCF_UP(I,LEVEL)
        ENDDO
      ENDDO
!
!  Check for negative values of QT and QCF: set to layer mean value if
!  negative values found.
      DO LEVEL=2,(LEVELS-1)
        DO I=1,POINTS
         IF ((QCF_MID(I,LEVEL).LT.1E-11).OR.(QCF_LOW(I,LEVEL).LT.1E-11)
     &                   .OR.(QCF_UP(I,LEVEL).LT.1E-11)) THEN
           QCF_MID(I,LEVEL)=QCF(I,LEVEL)
           QCF_LOW(I,LEVEL)=QCF(I,LEVEL)
           QCF_UP(I,LEVEL)=QCF(I,LEVEL)
         ENDIF
         IF ((QT_MID(I,LEVEL).LT.1E-11).OR.(QT_LOW(I,LEVEL).LT.1E-11)
     &                       .OR.(QT_UP(I,LEVEL).LT.1E-11)) THEN
           QT_MID(I,LEVEL)=Q(I,LEVEL)
           QT_LOW(I,LEVEL)=Q(I,LEVEL)
           QT_UP(I,LEVEL)=Q(I,LEVEL)
         ENDIF
        ENDDO
      ENDDO
!
!  Set values at top and bottom layer: note gradients are not used here,
!  they are simply set to the adjacent layer's values.
      DO I=1,POINTS
        T_MID(I,LEVELS)=T(I,LEVELS)
        T_UP(I,LEVELS)=T(I,LEVELS)
        T_LOW(I,LEVELS)=T(I,LEVELS)
        QT_MID(I,LEVELS)=Q(I,LEVELS)
        QT_UP(I,LEVELS)=Q(I,LEVELS)
        QT_LOW(I,LEVELS)=Q(I,LEVELS)
        QCF_MID(I,LEVELS)=QCF(I,LEVELS)
        QCF_UP(I,LEVELS)=QCF(I,LEVELS)
        QCF_LOW(I,LEVELS)=QCF(I,LEVELS)
        T_MID(I,1)=T(I,1)
        T_UP(I,1)=T(I,1)
        T_LOW(I,1)=T(I,1)
        QT_MID(I,1)=Q(I,1)
        QT_UP(I,1)=Q(I,1)
        QT_LOW(I,1)=Q(I,1)
        QCF_MID(I,1)=QCF(I,1)
        QCF_UP(I,1)=QCF(I,1)
        QCF_LOW(I,1)=QCF(I,1)
      ENDDO
!
!
      CALL GLUE_CLD(AK,BK,PSTAR,RHCRIT,LEVELS,RHCPT,POINTS,PFIELD,
     &        T_LOW,CF_BULK,QT_LOW,QCF_LOW,QCL,CF_LIQ,CF_ICE,ERROR)
!
      CALL GLUE_CLD(AK,BK,PSTAR,RHCRIT,LEVELS,RHCPT,POINTS,PFIELD,
     &        T_MID,CF_TEMP,QT_MID,QCF_MID,QCL_TEMP,CF_LIQ_TEMP,
     &                                CF_ICE_TEMP,ERROR)
!
!  The outputs of the cloud scheme combined with previous output
      DO LEVEL=1,LEVELS
        DO I=1,POINTS
          CF_AREA(I, LEVEL) = MAX(CF_BULK(I, LEVEL), CF_TEMP(I, LEVEL))
          CF_BULK(I, LEVEL) = CF_BULK(I, LEVEL) + CF_TEMP(I, LEVEL)
          CF_LIQ(I, LEVEL) = CF_LIQ(I, LEVEL) + CF_LIQ_TEMP(I, LEVEL)
          CF_ICE(I, LEVEL) = CF_ICE(I, LEVEL) + CF_ICE_TEMP(I, LEVEL)
          QCL(I, LEVEL) = QCL(I, LEVEL) + QCL_TEMP(I, LEVEL)
        ENDDO
      ENDDO
!
      CALL GLUE_CLD(AK,BK,PSTAR,RHCRIT,LEVELS,RHCPT,POINTS,PFIELD,
     &          T_UP,CF_TEMP,QT_UP,QCF_UP,QCL_TEMP,CF_LIQ_TEMP,
     &                                CF_ICE_TEMP,ERROR)
!
!  The outputs of the cloud scheme combined with previous output
      DO LEVEL=1,LEVELS
        DO I=1,POINTS
          CF_AREA(I, LEVEL) = MAX(CF_AREA(I, LEVEL), CF_TEMP(I, LEVEL))
          CF_BULK(I, LEVEL) = ONETHIRD * (CF_BULK(I, LEVEL) +
     &                                               CF_TEMP(I, LEVEL))
          CF_LIQ(I, LEVEL) = ONETHIRD * (CF_LIQ(I, LEVEL) +
     &                                           CF_LIQ_TEMP(I, LEVEL))
          CF_ICE(I, LEVEL) = ONETHIRD * (CF_ICE(I, LEVEL) +
     &                                           CF_ICE_TEMP(I, LEVEL))
          QCL(I, LEVEL) = ONETHIRD * (QCL(I, LEVEL)+QCL_TEMP(I, LEVEL))
!
!  Grid-box mean temperature and water vapour must now be recalculated,
!  from the grid-box mean QCL calculated above.
          Q(I, LEVEL) = Q(I, LEVEL) - QCL(I, LEVEL)
          T(I, LEVEL) = T(I, LEVEL) + QCL(I, LEVEL) * LCRCP
        ENDDO
      ENDDO
!
      RETURN
      END
