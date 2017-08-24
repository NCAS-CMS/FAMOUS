C *****************************COPYRIGHT******************************
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
      SUBROUTINE BL_LSP( P_FIELD,FIRST_ROW,ROW_LENGTH,ROWS,BL_LEVELS,
     &                   QCF,Q,T )
!
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected AND a full boundary layer treatment is not
!           performed.
!
! D Wilson    <- programmer
!
!  Model            Modification history from model version 4.4:
! version  Date
! 4.4      Sept 97  Originally Coded
!                                                 Damian Wilson
!
        IMPLICIT NONE
!
        INTEGER
     &    P_FIELD               ! IN   Number of points in the field
     &,   FIRST_ROW             ! IN   First row to be processed
     &,   ROW_LENGTH            ! IN   Length of each row
     &,   ROWS                  ! IN   Number of rows
     &,   BL_LEVELS             ! IN   Number of boundary layer levels
!
        REAL
     &    QCF(P_FIELD,BL_LEVELS) ! INOUT Ice water content
     &,   Q(P_FIELD,BL_LEVELS)   ! INOUT
!                                  IN    Vapour+liquid+ice content
!                                  OUT   Vapour+liquid content
     &,   T(P_FIELD,BL_LEVELS)   ! INOUT
!                                  IN    Liquid ice temperature
!                                  OUT   Liquid temperature
! Temporary Space
        INTEGER P1               ! First point
     &,         P2               ! Last point
     &,         I                ! Counter over points
     &,         J                ! Counter over boundary layer levels
        REAL NEWQCF              ! Temporary variable for QCF
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
      REAL LSRCP                 ! IN Latent heat of sublimation / Cp
      PARAMETER( LSRCP=((LC+LF)/CP) )
!
      P1=1+(FIRST_ROW-1)*ROW_LENGTH
      P2=P1+ROWS*ROW_LENGTH-1
!
      DO J=1,BL_LEVELS
        DO I=P1,P2
! Convert Q (vapour+liquid+ice) to (vapour+liquid)
          Q(I,J)=Q(I,J)-QCF(I,J)
! Check that Q is not negative
          IF (Q(I,J) .LT. 0.0) THEN
! Evaporate ice to keep Q positive, but don't let ice go negative
! itself
            NEWQCF=MAX(QCF(I,J)+Q(I,J),0.0)
            Q(I,J)=Q(I,J)+(QCF(I,J)-NEWQCF)
            QCF(I,J)=NEWQCF
          ENDIF
! Adjust T from T liquid ice to T liquid
          T(I,J)=T(I,J)+LSRCP*QCF(I,J)
        END DO
      END DO
! End the subroutine
      RETURN
      END
