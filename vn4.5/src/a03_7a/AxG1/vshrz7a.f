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
!!!   SUBROUTINE VSHR_Z1------------------------------------------------
!!!
!!!  Purpose: Calculate level 1 height and windshear for use in routine
!!!           PHYSIOL.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     7/97     New deck for MOSES II (R. Essery)
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!---------------------------------------------------------------------
      SUBROUTINE VSHR_Z1 (
     & P_FIELD,U_FIELD,LTIMER,
     & N_ROWS,FIRST_ROW,ROW_LENGTH,
     & AKH,BKH,EXNER,PSTAR,Q,QCF,QCL,T,U_1,V_1,U_0,V_0,
     & VSHR,Z1
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid.
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
     &,N_ROWS                      ! IN No. of rows to be treated.
     &,FIRST_ROW                   ! IN First row of data to be treated.
     &,ROW_LENGTH                  ! IN No. of points in one row.

      LOGICAL
     & LTIMER                      ! IN Logical switch for TIMER diags

      REAL
     & AKH(2)                      ! IN Hybrid 'A' for layer 1.
     &,BKH(2)                      ! IN Hybrid 'B' for layer 1.
     &,EXNER(P_FIELD,2)            ! IN Exner function.
     &,PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).
     &,Q(P_FIELD)                  ! IN Specific humidity (kg/kg air).
     &,QCF(P_FIELD)                ! IN Cloud ice (kg/kg air).
     &,QCL(P_FIELD)                ! IN Cloud liquid water (kg/kg air).
     &,T(P_FIELD)                  ! IN Atmospheric temperature (K).
     &,U_1(U_FIELD)                ! IN W'ly wind component (m/s)
     &,V_1(U_FIELD)                ! IN S'ly wind component (m/s)
     &,U_0(U_FIELD)                ! IN W'ly component of surface
!                                  !    current (m/s).
     &,V_0(U_FIELD)                ! IN S'ly component of surface
!                                  !    current (m/s).


      REAL
     & VSHR(P_FIELD)               ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m/s).
     &,Z1(P_FIELD)                 ! OUT Height of lowest level (m).

      REAL
     & DZL(P_FIELD)                ! Depth of layer 1 (m).
     &,RDZ(P_FIELD)                ! Reciprocal of the height of level 1
     &,TV(P_FIELD)                 ! Virtual temperature.
     &,U_1_P(P_FIELD)              ! U_1 on P-grid.
     &,U_0_P(P_FIELD)              ! U_0 on P-grid.
     &,V_1_P(P_FIELD)              ! V_1 on P-grid.
     &,V_0_P(P_FIELD)              ! V_0 on P-grid.
     &,ZLB(P_FIELD,0:1)            ! ZLB(,K) is the height of the upper
!                                  ! boundary of layer K ( = 0 for K=0).

      REAL
     & USHEAR       ! U-component of surface-to-lowest-level wind shear.
     &,VSHEAR       ! V-component of surface-to-lowest-level wind shear.
     &,VSHR2        ! Square of magnitude of surface-to-lowest-level
!                   ! wind shear.

      INTEGER
     & I            ! Loop counter (horizontal field index).
     &,P1           ! First P-point to be processed.
     &,P_POINTS     ! No. of P-points being processed.
     &,U1           ! First UV-point to be processed.
     &,N_P_ROWS     ! No. of P-rows being processed.
     &,N_U_ROWS     ! No. of UV-rows being processed.
     &,U_POINTS     ! No. of UV-points being processed.


      EXTERNAL Z

      N_P_ROWS = N_ROWS
      N_U_ROWS = N_ROWS + 1
      P_POINTS = N_P_ROWS * ROW_LENGTH
      U_POINTS = N_U_ROWS * ROW_LENGTH
      P1 = 1 + (FIRST_ROW-1)*ROW_LENGTH
      U1 = 1 + (FIRST_ROW-2)*ROW_LENGTH

      DO I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
      ENDDO

      CALL Z (P_POINTS,EXNER(P1,1),EXNER(P1,2),PSTAR(P1),
     &        AKH,BKH,Q(P1),QCF(P1),
     &        QCL(P1),T(P1),ZLB(P1,0),TV(P1),
     &        ZLB(P1,1),DZL(P1),RDZ(P1),LTIMER)

      DO I=P1,P1+P_POINTS-1
        Z1(I)=RDZ(I)
      ENDDO

      CALL UV_TO_P(U_1(U1),U_1_P(P1),
     &             U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
      CALL UV_TO_P(V_1(U1),V_1_P(P1),
     &             U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
      CALL UV_TO_P(U_0(U1),U_0_P(P1),
     &             U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
      CALL UV_TO_P(V_0(U1),V_0_P(P1),
     &             U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)

      DO I=P1,P1+P_POINTS-1
        USHEAR = U_1_P(I) - U_0_P(I)
        VSHEAR = V_1_P(I) - V_0_P(I)
        VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
        VSHR(I) = SQRT(VSHR2)
      ENDDO

      RETURN
      END
