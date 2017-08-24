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
CLL   SUBROUTINE MG_CALC_Z
CLL
CLL   PURPOSE:
CLL   -------
CLL   SETS VERTICAL LEVELS FOR ALL GRIDS.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_CALC_Z(Z_Q_GRIDS,Z_MID_GRIDS,START_ADDRESS_Z,
     &                     LAST_ADDRESS_Z,NGRIDS,KMAX,RES_DIRS)

      IMPLICIT NONE

      INTEGER
     &  NGRIDS        !IN NUMBER OF GRIDS

      INTEGER
     &  KMAX(NGRIDS) !IN NUMBER OF LEVELS ON EACH GRID
     &, START_ADDRESS_Z(NGRIDS) !IN START ADDRESS FOR EACH ARRAY OF
     &                          !   LEVELS ON EACH GRID
     &, LAST_ADDRESS_Z          !IN USED TO DIMENSION SPACE REQUIRED
     &, RES_DIRS(NGRIDS)        !IN NUMBER OF DIRECTIONS TO RESTRICT IN

      REAL
     &  Z_Q_GRIDS(LAST_ADDRESS_Z) !INOUT HOLDS ALL Z LEVELS AT Q POINTS
     &, Z_MID_GRIDS(LAST_ADDRESS_Z) !INOUT HOLDS Z VALUES AT MID-POINTS
     &                              !      BETWEEN Q POINTS IN VERTICAL

C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C    LOCAL VARIABLES.
      INTEGER
     & I,J

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. SET LEVELS FOR EACH GRID.
CL----------------------------------------------------------------------

      DO I= NGRIDS-1,1,-1
        IF(MOD(RES_DIRS(I+1),2).EQ.1) THEN
C SET Q LEVELS
          DO J=1,KMAX(I)
            Z_Q_GRIDS(START_ADDRESS_Z(I)+J-1) =
     &                Z_Q_GRIDS(START_ADDRESS_Z(I+1)+2*J-2)
          END DO
C INTERPOLATE MID-LEVELS, CURRENTLY SET TO SIMPLE LINEAR
C FIRST COPY BOTTOM AND TOP LEVEL WHICH REMAIN UNCHANGED
          Z_MID_GRIDS(START_ADDRESS_Z(I)) =
     &                Z_MID_GRIDS(START_ADDRESS_Z(I+1))
          Z_MID_GRIDS(START_ADDRESS_Z(I)+ KMAX(I)) =
     &                Z_MID_GRIDS(START_ADDRESS_Z(I+1)+KMAX(I+1))
          DO J=1,KMAX(I)-1
            Z_MID_GRIDS(START_ADDRESS_Z(I)+J) = .5*
     &               (Z_Q_GRIDS(START_ADDRESS_Z(I)+J-1) +
     &                Z_Q_GRIDS(START_ADDRESS_Z(I)+J))
          END DO
        ELSE
C SAME VERTICAL LEVELS AS LAST GRID SO SIMPLY COPY.
          DO J=1,KMAX(I)
            Z_Q_GRIDS(START_ADDRESS_Z(I)+J-1) =
     &                Z_Q_GRIDS(START_ADDRESS_Z(I+1)+J-1)
          END DO
          DO J=1,KMAX(I)+1
            Z_MID_GRIDS(START_ADDRESS_Z(I)+J-1) =
     &                Z_MID_GRIDS(START_ADDRESS_Z(I+1)+J-1)
          END DO
        END IF
      END DO

CL    END OF ROUTINE MG_CALC_Z

      RETURN
      END
CLL   SUBROUTINE MG_K_MATRIX
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALCULATES CURRENT VALUE OF LEFT-HAND-SIDE OF P.D.E. FOR A GIVEN
CLL   VALUE OF I AND J. SETS MATRIX ENTRIES TO SOLVE FOR CORRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_K_MATRIX(
     &                       LHS,A0,A1,A2,Q,A,B,C1,C2,DEF,D,E,F,G,
     &                       I_DIM,J_DIM,K_DIM,
     &                       COS_P_LATITUDE,SEC_P_LATITUDE,
     &                       COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                       LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,
     &                       I,J,VERSION,Z_Q,Z_MID,I_NT,J_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION

      INTEGER
     &  I_DIM         !IN. NUMBER OF NODES IN THE I-DIRECTION
     &, J_DIM         !IN. NUMBER OF NODES IN THE J-DIRECTION
     &, K_DIM         !IN. NUMBER OF NODES IN THE K-DIRECTION
     &, I             !IN. VALUE OF I FOR WHICH LHS AND COEFFS ARE REQ.
     &, J             !IN. VALUE OF J FOR WHICH LHS AND COEFFS ARE REQ.
     &, VERSION       !IN. VERSION OF MULTIGRID SCHEME.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !IN. SOLUTION
     &, A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)!IN. COSINE OF LATITUDE AT B POINTS
     &, SEC_P_LATITUDE(I_DIM,J_DIM)!IN. 1./COSINE OF LATITUDE AT Q POINT
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  EARTH_RADIUS_INVERSE    !IN. 1/EARTH RADIUS
     &, LATITUDE_STEP_INVERSE   !IN. 1/LATITUDE STEP LENGTH IN RADIANS
     &, LONGITUDE_STEP_INVERSE  !IN. 1/LONGITUDE STEP LENGTH IN RADIANS

      REAL
     &  LHS(K_DIM)         !OUT. LEFT-HAND-SIDE.
     &, A0(K_DIM)          !OUT. DIAGONAL MATRIX COEFFICIENT.
     &, A1(K_DIM)          !OUT. SUPER-DIAGONAL MATRIX COEFFICIENT.
     &, A2(K_DIM)          !OUT. SUB-DIAGONAL MATRIX COEFFICIENT.
C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  K,I2,IC

      REAL
     &  SCALAR
     &, SCALARS
     &, SCALAR1
     &, SCALAR2

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. CALCULATE LEFT-HAND-SIDE OF P.D.E.
CL               AND MATRIX COEFFICIENTS.
CL----------------------------------------------------------------------

C ----------------------------------------------------------------------
CL    SECTION 1.1 CALCULATE K-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

C SET ZERO COEFFICIENTS IN MATRICES.
      A1(K_DIM) = 0.
      A2(1) = 0.

      K=1
      SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
      SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
      LHS(K) = SCALAR*(Q(I,J,K+1)*C2(I,J,K+1)
     &                    -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
     &                      - Q(I,J,K)*G(I,J,K)
      A0(K) = - SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K)
     &            - G(I,J,K)
      A1(K) = SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K+1)
      IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
        LHS(K) = LHS(K)  + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*
     &                     (Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
        A0(K) = A0(K) - .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*C2(I,J,K)
        A1(K) = A1(K) + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*C2(I,J,K+1)
      END IF
      K=K_DIM
      SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
      SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
      LHS(K) = - SCALAR*(Q(I,J,K)*C2(I,J,K)
     &                 -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
     &                      - Q(I,J,K)*G(I,J,K)
      A0(K) = - SCALAR*C1(I,J,K-1)*SCALAR2*C2(I,J,K)
     &            - G(I,J,K)
      A2(K) = SCALAR*C1(I,J,K-1)*SCALAR1*C2(I,J,K-1)
      IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
        LHS(K) = LHS(K)  + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*
     &                     (Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
        A0(K) = A0(K) + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*C2(I,J,K)
        A2(K) = A2(K) - .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*C2(I,J,K-1)
      END IF
      DO K=2,K_DIM-1
        SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
        SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
        SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
        LHS(K) =  SCALAR*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *C1(I,J,K)*SCALAR1
     &               -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *C1(I,J,K-1)*SCALAR2)
     &               +.5*DEF(I,J,K)*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *F(I,J,K)*SCALAR1
     &               +(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *F(I,J,K-1)*SCALAR2)
     &                      - Q(I,J,K)*G(I,J,K)
        A0(K) = - SCALAR*(C1(I,J,K)*SCALAR1
     &                              +C1(I,J,K-1)*SCALAR2)*C2(I,J,K)
     &               +.5*DEF(I,J,K)*(F(I,J,K-1)*SCALAR2
     &                               -F(I,J,K)*SCALAR1)*C2(I,J,K)
     &            - G(I,J,K)
        A1(K) = SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K+1)
     &             + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*C2(I,J,K+1)
        A2(K) = SCALAR*C1(I,J,K-1)*SCALAR1*C2(I,J,K-1)
     &              - .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*C2(I,J,K-1)
      END DO

C ----------------------------------------------------------------------
CL    SECTION 1.2 CALCULATE AND ADD ON I-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

      IF(I_NT) THEN
        SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &           LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
        SCALARS= EARTH_RADIUS_INVERSE*LONGITUDE_STEP_INVERSE

C NOT EXECUTED AT POLES IN GLOBAL VERSION.
        IF(.NOT.(J.EQ.1 .AND. VERSION.LT.3 .AND. J_NT) .AND.
     &     .NOT.(J.EQ.J_DIM .AND. VERSION.LT.3 .AND. J_NT)) THEN
          IF(I.EQ.1) THEN
CL FIRST POINT ON ROW
            IF(VERSION.LT.3) THEN
CL GLOBAL
              SCALAR1= SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
              DO K=1,K_DIM
                A0(K)  = A0(K) - (A(I_DIM,J,K)+A(I,J,K))*SCALAR1
     &                 +.5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                 (D(I_DIM,J,K)-D(I,J,K))
                LHS(K) = LHS(K)+(Q(I_DIM,J,K)*A(I_DIM,J,K) - Q(I,J,K) *
     &              (A(I_DIM,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &               *SCALAR1
     &                   + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                   (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I_DIM,J,K)*
     &                   (Q(I,J,K)-Q(I_DIM,J,K)))
              END DO
            ELSE
CL LIMITED AREA
              SCALAR1= SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
              DO K=1,K_DIM
                A0(K)  = A0(K) - 2.*A(I,J,K)*SCALAR1
                LHS(K) = LHS(K) - 2.*(Q(I,J,K) * A(I,J,K) -
     &                           Q(I+1,J,K)*A(I,J,K))*SCALAR1
              END DO
            END IF
          ELSE IF (I.EQ.I_DIM) THEN
CL LAST POINT ON ROW
            IF(VERSION.LT.3) THEN
CL GLOBAL
              SCALAR1= SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
              DO K=1,K_DIM
                A0(K)  = A0(K) - (A(I-1,J,K)+A(I,J,K))*SCALAR1
     &                 +.5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                 (D(I-1,J,K)-D(I,J,K))
                LHS(K) = LHS(K) + (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &              (A(I-1,J,K) + A(I,J,K)) + Q(1,J,K)*A(I,J,K))
     &               *SCALAR1
     &                   + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                   (D(I,J,K)*(Q(1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &                   (Q(I,J,K)-Q(I-1,J,K)))
              END DO
            ELSE
CL LIMITED AREA
              SCALAR1= SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
              DO K=1,K_DIM
                A0(K)  = A0(K) - 2.*A(I-1,J,K)*SCALAR1
                LHS(K) = LHS(K) + 2.*(Q(I-1,J,K) * A(I-1,J,K) -
     &                           Q(I,J,K)*A(I-1,J,K))*SCALAR1
              END DO
            END IF
          ELSE
CL ALL OTHER POINTS ON ROW
            SCALAR1 = SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
            DO K=1,K_DIM
              A0(K)  = A0(K) - (A(I-1,J,K)+A(I,J,K))*SCALAR1
     &                 +.5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                 (D(I-1,J,K)-D(I,J,K))
              LHS(K) = LHS(K) + (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &              (A(I-1,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &               *SCALAR1
     &                   + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                   (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &                   (Q(I,J,K)-Q(I-1,J,K)))
            END DO
          END IF
        END IF
      END IF

C ----------------------------------------------------------------------
CL    SECTION 1.3 CALCULATE J-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

      IF (J_NT) THEN
        SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &           LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
        SCALARS= EARTH_RADIUS_INVERSE*LATITUDE_STEP_INVERSE
        IF(J.EQ.1) THEN
CL NORTHERN BOUNDARY
          IF(VERSION.LT.3) THEN
CL GLOBAL
            DO K=1,K_DIM
              I2 = I_DIM/2
              SCALAR1 = 0.
              SCALAR2 = 0.

C SUM ALL TERMS AT POLE.
              DO IC= 1,I_DIM/2
                SCALAR1 = SCALAR1 + ((Q(IC,2,K)-Q(IC,1,K))
     &                       *B(IC,1,K)*COS_V_LATITUDE(IC,1)
     &                        -(Q(IC+I2,1,K) -Q(IC+I2,2,K))
     &                        *B(IC+I2,1,K)*COS_V_LATITUDE(IC+I2,1))
     &                        *SEC_P_LATITUDE(IC,1)*SCALAR
     &                    +DEF(1,1,K)*SCALARS*
     &                    (E(IC,1,K)*(Q(IC,1,K)-Q(IC,2,K))+
     &                     E(IC+I2,1,K)*(Q(IC+I2,1,K)-Q(IC+I2,2,K)))
                 SCALAR2 = SCALAR2 - (B(IC,1,K)*COS_V_LATITUDE(IC,1)
     &                    +B(IC+I2,1,K)*COS_V_LATITUDE(IC+I2,1))
     &                     *SEC_P_LATITUDE(IC,1)*SCALAR
     &                    +DEF(1,1,K)*SCALARS*
     &                    (E(IC,1,K)+E(IC+I2,1,K))
              END DO

              LHS(K) = LHS(K) + SCALAR1
              A0(K) = A0(K) + SCALAR2
            END DO

C RESCALE VALUES BY I_DIM SINCE RESIDUAL IS STORED AS R/I_DIM AT POLES
            DO K=1,K_DIM
              LHS(K) = LHS(K) /I_DIM
              A0(K) = A0(K) /I_DIM
              A1(K) = A1(K) /I_DIM
              A2(K) = A2(K) /I_DIM
            END DO

          ELSE

CL LIMITED AREA
            DO K=1,K_DIM
              LHS(K) = LHS(K) -2.*(Q(I,J,K)*B(I,J,K)
     &                        -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
              A0(K)  = A0(K) - 2.*B(I,J,K)*COS_V_LATITUDE(I,J)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
            END DO
          END IF

        ELSE IF (J.EQ.J_DIM) THEN

CL SOUTHERN BOUNDARY
          IF(VERSION.LT.3) THEN
CL GLOBAL
            DO K=1,K_DIM
              I2 = I_DIM/2
              SCALAR1 = 0.
              SCALAR2 = 0.

C SUM ALL TERMS AT POLE.
              DO IC= 1,I_DIM/2
                SCALAR1 = SCALAR1 +((Q(IC,J_DIM-1,K) -Q(IC,J_DIM,K))
     &                      *B(IC,J_DIM-1,K)*COS_V_LATITUDE(IC,J_DIM-1)
     &                      -(Q(IC+I2,J_DIM,K)-Q(IC+I2,J_DIM-1,K))
     &                      *B(IC+I2,J_DIM-1,K)*
     &                      COS_V_LATITUDE(IC+I2,J_DIM-1))
     &                      *SEC_P_LATITUDE(IC,J_DIM)*SCALAR
     &                    +DEF(1,J_DIM,K)*SCALARS*
     &                (E(IC,J_DIM-1,K)*(Q(IC,J_DIM-1,K)-Q(IC,J_DIM,K))+
     &                     E(IC+I2,J_DIM-1,K)*(Q(IC+I2,J_DIM-1,K)
     &                                         -Q(IC+I2,J_DIM,K)))
                SCALAR2 = SCALAR2 - (B(IC,J_DIM-1,K)*
     &                      COS_V_LATITUDE(IC,J_DIM-1)
     &                      +B(IC+I2,J_DIM-1,K)*
     &                      COS_V_LATITUDE(IC+I2,J_DIM-1))
     &                      *SEC_P_LATITUDE(IC,J_DIM)*SCALAR
     &                      -DEF(1,J_DIM,K)*SCALARS*
     &                      (E(IC,J_DIM-1,K)+E(IC+I2,J_DIM-1,K))
               END DO

               LHS(K) = LHS(K) + SCALAR1
               A0(K) = A0(K) + SCALAR2
             END DO

C RESCALE VALUES BY I_DIM SINCE RESIDUAL IS STORED AS R/I_DIM AT POLES
            DO K=1,K_DIM
               LHS(K) = LHS(K) /I_DIM
               A0(K) = A0(K) /I_DIM
               A1(K) = A1(K) /I_DIM
               A2(K) = A2(K) /I_DIM
            END DO

          ELSE
CL LIMITED AREA
            DO K=1,K_DIM
              LHS(K) = LHS(K) + 2.*(Q(I,J-1,K)*B(I,J-1,K)
     &                    -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
              A0(K)  = A0(K) - 2.*B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
             END DO
          END IF
        ELSE
          DO K=1,K_DIM
CL ALL OTHER POINTS
            LHS(K) = LHS(K) + ((Q(I,J-1,K)*B(I,J-1,K)
     &                    -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                    -(Q(I,J,K)*B(I,J,K)
     &                    -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                       E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
            A0(K)  = A0(K) - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                      +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J,K)-E(I,J-1,K))
          END DO
        END IF
      END IF

CL    END OF ROUTINE MG_K_MATRIX

      RETURN
      END
CLL   SUBROUTINE MG_LEFT_HAND_SIDE
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALCULATES CURRENT VALUE OF LEFT-HAND-SIDE OF P.D.E.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_LEFT_HAND_SIDE(
     &                             Q,A,B,C1,C2,DEF,D,E,F,G,
     &                             LHS,I_DIM,J_DIM,K_DIM,
     &                             COS_P_LATITUDE,SEC_P_LATITUDE,
     &                             COS_V_LATITUDE,
     &                             EARTH_RADIUS_INVERSE,
     &                             LATITUDE_STEP_INVERSE,
     &                             LONGITUDE_STEP_INVERSE,VERSION,
     &                             I_NT,J_NT,K_NT,Z_Q,Z_MID)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM         !IN. NUMBER OF NODES IN THE I-DIRECTION
     &, J_DIM         !IN. NUMBER OF NODES IN THE J-DIRECTION
     &, K_DIM         !IN. NUMBER OF NODES IN THE K-DIRECTION
     &, VERSION       !IN. VERSION OF SCHEME BEING USED.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !IN. SOLUTION
     &, A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)!IN. COSINE OF LATITUDE AT B POINTS
     &, SEC_P_LATITUDE(I_DIM,J_DIM)!IN. 1./COSINE OF LATITUDE AT Q POINT
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z AT MID POINTS BETWEEN Q POINTS.

      REAL
     &  EARTH_RADIUS_INVERSE    !IN. 1/EARTH RADIUS
     &, LATITUDE_STEP_INVERSE   !IN. 1/LATITUDE STEP LENGTH IN RADIANS
     &, LONGITUDE_STEP_INVERSE  !IN. 1/LONGITUDE STEP LENGTH IN RADIANS

      REAL
     &  LHS(I_DIM,J_DIM,K_DIM)  !OUT. LEFT-HAND-SIDE.
C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,J,K,I2,JB,JE

      REAL
     &  SCALAR
     &, SCALARS
     &, SCALAR1
     &, SCALAR2

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. CALCULATE LEFT-HAND-SIDE OF P.D.E.
CL----------------------------------------------------------------------

CL    LOOP OVER LEVELS.

      DO K = 1,K_DIM

        DO J=1,J_DIM
          DO I=1,I_DIM
            LHS(I,J,K) = - Q(I,J,K)*G(I,J,K)
          END DO
        END DO

C ----------------------------------------------------------------------
CL    SECTION 1.1 CALCULATE I-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

        IF (I_NT) THEN
          SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &             LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
          SCALARS= EARTH_RADIUS_INVERSE*LONGITUDE_STEP_INVERSE
          IF(VERSION.LT.3) THEN
C I DIRECTION DERIVATIVES ARE MEANINGLESS AT THE BOUNDARIES.
            JB = 2
            JE = J_DIM-1
          ELSE
            JB = 1
            JE = J_DIM
          END IF
          DO J = JB,JE
            DO I= 2,I_DIM-1
              LHS(I,J,K)= LHS(I,J,K) +
     &                   (Q(I-1,J,K)*A(I-1,J,K)- Q(I,J,K) * (A(I-1,J,K)
     &                   +A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &                   *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                   *SCALAR
     &                   + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                   (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &                   (Q(I,J,K)-Q(I-1,J,K)))
            END DO
C END POINTS ON EACH ROW.
            IF (VERSION.LT.3) THEN
              LHS(1,J,K) = LHS(1,J,K) +
     &                    (Q(I_DIM,J,K)*A(I_DIM,J,K) - Q(1,J,K) *
     &                 (A(I_DIM,J,K) + A(1,J,K)) + Q(2,J,K)*A(1,J,K))
     &                   *SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
     &                   *SCALAR
     &                     + .5*DEF(1,J,K)*SCALARS*SEC_P_LATITUDE(1,J)*
     &                     (D(1,J,K)*(Q(2,J,K)-Q(1,J,K))+D(I_DIM,J,K)*
     &                      (Q(1,J,K)-Q(I_DIM,J,K)))
              LHS(I_DIM,J,K)=LHS(I_DIM,J,K) +
     &                       (Q(I_DIM-1,J,K)*A(I_DIM-1,J,K)-Q(I_DIM,J,K)
     &                    *(A(I_DIM-1,J,K) + A(I_DIM,J,K)) + Q(1,J,K)
     &                    *A(I_DIM,J,K))*SEC_P_LATITUDE(I_DIM,J)
     &                    *SEC_P_LATITUDE(I_DIM,J)*SCALAR
     &             + .5*DEF(I_DIM,J,K)*SCALARS*SEC_P_LATITUDE(I_DIM,J)*
     &            (D(I_DIM,J,K)*(Q(1,J,K)-Q(I_DIM,J,K))+D(I_DIM-1,J,K)*
     &             (Q(I_DIM,J,K)-Q(I_DIM-1,J,K)))
            ELSE
              LHS(1,J,K) = LHS(1,J,K)
     &                     +2.*( Q(2,J,K) - Q(1,J,K)) * A(1,J,K)
     &                     *SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
     &                     *SCALAR
              LHS(I_DIM,J,K) = LHS(I_DIM,J,K)
     &                       + 2.*(Q(I_DIM-1,J,K) - Q(I_DIM,J,K))
     &                    * A(I_DIM-1,J,K)*SEC_P_LATITUDE(I_DIM,J)
     &                    *SEC_P_LATITUDE(I_DIM,J)*SCALAR
            END IF
          END DO
        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.2 CALCULATE J-DIRECTION DERIVATIVES.
CL                ADD IN NON-DERIVATIVE TERMS AT POLES.
C ----------------------------------------------------------------------

        IF (J_NT) THEN
          SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &             LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
          SCALARS= EARTH_RADIUS_INVERSE*LATITUDE_STEP_INVERSE
          DO J = 2,J_DIM-1
            DO I= 1,I_DIM
              LHS(I,J,K) = LHS(I,J,K) +
     &                     ((Q(I,J-1,K)-Q(I,J,K))
     &                      *B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                      -(Q(I,J,K)-Q(I,J+1,K))
     &                      *B(I,J,K)*COS_V_LATITUDE(I,J))
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
     &                     +DEF(I,J,K)*.5*SCALARS*
     &                     (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                      E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
            END DO
          END DO

C ----------------------------------------------------------------------
CL    SECTION 1.2.1 BOUNDARY TERMS IN J-DIRECTION.
C ----------------------------------------------------------------------

          IF (VERSION.LT.3.AND.J_NT) THEN
CL    GLOBAL VERSION SO CALCULATE VALUE AT POLES.
CL    NORTH POLE

            I2 = I_DIM/2

C SUM ALL TERMS AT POLE.
            DO I= 1,I_DIM/2
              LHS(1,1,K) = LHS(1,1,K) +
     &                   ((Q(I,2,K)-Q(I,1,K))
     &                    *B(I,1,K)*COS_V_LATITUDE(I,1)
     &                   -(Q(I+I2,1,K)-Q(I+I2,2,K))
     &                    *B(I+I2,1,K)*COS_V_LATITUDE(I+I2,1))
     &                    *SEC_P_LATITUDE(I,1)*SCALAR
     &                    +DEF(1,1,K)*SCALARS*
     &                    (E(I,1,K)*(Q(I,1,K)-Q(I,2,K))+
     &                     E(I+I2,1,K)*(Q(I+I2,1,K)-Q(I+I2,2,K)))
            END DO

CL    SOUTH POLE

C SUM ALL TERMS AT POLE.
            DO I= 1,I_DIM/2
              LHS(1,J_DIM,K) = LHS(1,J_DIM,K) +
     &                       ((Q(I,J_DIM-1,K)-Q(I,J_DIM,K))
     &                        *B(I,J_DIM-1,K)*
     &                        COS_V_LATITUDE(I,J_DIM-1)
     &                        -(Q(I+I2,J_DIM,K)-Q(I+I2,J_DIM-1,K))
     &                        *B(I+I2,J_DIM-1,K)*
     &                        COS_V_LATITUDE(I+I2,J_DIM-1))
     &                        *SEC_P_LATITUDE(I,J_DIM)*SCALAR
     &                    +DEF(1,J_DIM,K)*SCALARS*
     &                (E(I,J_DIM-1,K)*(Q(I,J_DIM-1,K)-Q(I,J_DIM,K))+
     &                     E(I+I2,J_DIM-1,K)*(Q(I+I2,J_DIM-1,K)
     &                            -Q(I2+I,J_DIM,K)))
            END DO

C SET ALL POLAR VALUES TO BE THE SAME
C NOTE VALUES NEED AVERAGING, THIS IS DONE AFTER VERTICAL TERM CODE.
            DO I= 2,I_DIM
              LHS(I,1,K) = LHS(1,1,K)
              LHS(I,J_DIM,K) = LHS(1,J_DIM,K)
            END DO

          ELSE
CL LIMITED AREA CODE, INCLUDE BOUNDARY DERIVATIVES.
            SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &               LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
            J = J_DIM
            DO I= 1,I_DIM
              LHS(I,1,K) = LHS(I,1,K) - 2.*(Q(I,1,K)-Q(I,2,K))
     &                      *B(I,1,K)*COS_V_LATITUDE(I,1)
     &                      *SEC_P_LATITUDE(I,1)*SCALAR
              LHS(I,J,K) = LHS(I,J,K) + 2.*(Q(I,J-1,K)-Q(I,J,K))
     &                    *B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
            END DO

          END IF

        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.3 CALCULATE K-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

        IF (K_NT) THEN
          IF(K.EQ.1) THEN
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            DO J=1,J_DIM
              DO I=1,I_DIM
                LHS(I,J,K) = LHS(I,J,K) + SCALAR*
     &                     (Q(I,J,K+1)*C2(I,J,K+1)
     &                    -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
                IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                  LHS(I,J,K) = LHS(I,J,K) + .5*SCALAR1*DEF(I,J,K)*
     &                            F(I,J,K)*(Q(I,J,K+1)*C2(I,J,K+1)
     &                                      -Q(I,J,K)*C2(I,J,K))
                END IF
              END DO
            END DO
          ELSE IF (K.EQ.K_DIM) THEN
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            DO J=1,J_DIM
              DO I=1,I_DIM
                LHS(I,J,K) = LHS(I,J,K) - SCALAR*
     &                      (Q(I,J,K)*C2(I,J,K)
     &                 -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
                IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                  LHS(I,J,K) = LHS(I,J,K)  + .5*SCALAR2*DEF(I,J,K)*
     &                            F(I,J,K-1)*(Q(I,J,K)*C2(I,J,K)
     &                                      -Q(I,J,K-1)*C2(I,J,K-1))
                END IF
              END DO
            END DO
          ELSE
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            DO J=1,J_DIM
              DO I=1,I_DIM
                LHS(I,J,K) = LHS(I,J,K) + SCALAR*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *C1(I,J,K)*SCALAR1
     &               -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *C1(I,J,K-1)*SCALAR2)
     &               +.5*DEF(I,J,K)*
     &               (F(I,J,K)*SCALAR1*
     &                (Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &               +(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *F(I,J,K-1)*SCALAR2)
              END DO
            END DO
          END IF
        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.4 AVERAGE GLOBAL POLAR VALUES.
C ----------------------------------------------------------------------

        IF (VERSION.LT.3.AND.J_NT) THEN
          DO I=1,I_DIM
            LHS(I,1,K) = LHS(I,1,K) / I_DIM
            LHS(I,J_DIM,K) = LHS(I,J_DIM,K) / I_DIM
          END DO
        END IF

CL    END LOOP OVER LEVELS.
      END DO

CL    END OF ROUTINE MG_LEFT_HAND_SIDE

      RETURN
      END
CLL   SUBROUTINE MG_CNTL
CLL
CLL   PURPOSE:
CLL   -------
CLL   top level multi-grid control routine to solve 1 p.d.e. in 1,2 or 3
CLL   dimensions on a choice of domains.
CLL
CLL   d ( A dQ ) + d ( B dQ ) + C1 d ( C2 dQ )
CLL   --    --     --    --        --     --
CLL   dx    dx     dy    dy        dz     dz
CLL
CLL   + DEF (  D dQ + E dQ + F dQ )  - G Q = RHS
CLL              --     --     --
CLL              dx     dy     dz
CLL
CLL   but with (x,y) replaced by spherical polars
CLL
CLL   See Variable VERSION for currently implemented domain options.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL
CLL
CLL   Documentation: FR SOFTWARE DOCUMENTATION
CLL                  Multigrid Methods for Elliptic Equations.
CLL                  Version 2. December 1993.
CLL                  Mark H. Mawson.
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_CNTL(MXNGRDS,I_LENGTH,J_LENGTH,K_LENGTH,
     &                   MAXITS,TOL_RES,IPRINT,
     &                   KSMOOTH,NPRE,NPOST,NCOARSE,W,
     &                   KREST,NCGC,A,B,C1,C2,DEF,D,E,F,G,Q,RHS,
     &                   COS_P_LATITUDE,SEC_P_LATITUDE,
     &                   COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                   LONGITUDE_STEP_INVERSE,LATITUDE_STEP_INVERSE,
     &                   WORST_SMOOTHING_RATE,VERSION,K_BC,Z_Q,Z_MID)

      IMPLICIT NONE

      INTEGER
     &  VERSION    !IN  1 - Global
     &             !    2 - Global with equatorial solution prescribed
     &             !    3 - Limited Area, Neumann B.C`s in x and y.
     &             !    4 - Limited Area, Dirichlet B.C`s in x and y.
     &, K_BC       !IN  Code for Boundary condition in K direction.
     &             !    1 - Neumann Boundary condition bottom & top.
     &             !    2 - Neumann condition bottom, Dirichlet top.
     &             !    3 - Dirichlet condition bottom, Neumann top.
     &             !    4 - Dirichlet condition bottom & top.

      INTEGER
     &  MXNGRDS    !IN. MAXIMUM NUMBER OF GRIDS ALLOWED.
     &, I_LENGTH   !IN. NUMBER OF POINTS IN I DIRECTION.
     &, J_LENGTH   !IN. NUMBER OF POINTS IN J DIRECTION.
     &, K_LENGTH   !IN. NUMBER OF POINTS IN K DIRECTION.
     &, MAXITS     !IN. MAX NO OF FAS ITERATIONS WITHOUT CONVERGENCE
     &, IPRINT     !IN. PARAMETER CONTROLLING QUANTITY OF OUTPUT:
     &             !    0 - NONE
     &             !    1 - NUMBER OF ITERATIONS REQUIRED AND TIME TAKEN
     &             !    2 - & CONVERGENCE HISTORY
     &             !    3 - & JOB DESCRIPTION

      INTEGER
     & KSMOOTH  !IN. KIND OF ITERATIVE METHOD USED AS A SMOOTHER:
     &          !    1 - I-LINE JACOBI
     &          !    2 - J-LINE JACOBI
     &          !    3 - K-LINE JACOBI
     &          !    4 - I&J-LINE JACOBI
     &          !    5 - I&K LINE JACOBI
     &          !    6 - J&K-LINE JACOBI
     &          !    7 - 3D ALTERNATING LINE JACOBI
     &          !    8 - I-LINE GAUSS-SEIDEL
     &          !    9 - J-LINE GAUSS-SEIDEL
     &          !   10 - K-LINE GAUSS-SEIDEL
     &          !   11 - I&J-LINE GAUSS-SEIDEL
     &          !   12 - I&K-LINE GAUSS-SEIDEL
     &          !   13 - J&K-LINE GAUSS-SEIDEL
     &          !   14 - 3D ALTERNATING LINE GAUSS-SEIDEL
     &          !   15 - 3D ALTERNATING SYMMETRIC LINE GAUSS-SEIDEL
     &          !   16 - I-LINE ZEBRA
     &          !   17 - J-LINE ZEBRA
     &          !   18 - K-LINE ZEBRA
     &          !   19 - I&J-LINE ZEBRA
     &          !   20 - I&K-LINE ZEBRA
     &          !   21 - J&K-LINE ZEBRA
     &          !   22 - 3D ALTERNATING LINE ZEBRA
     &          !   23 - JACOBI POINT SMOOTHER
     &          !   24 - GAUSS-SEIDEL POINT SMOOTHER
     &          !   25 - SYMMETRIC GAUSS-SEIDEL POINT SMOOTHER
     &          !   26 - RED-BLACK POINT SMOOTHER

      INTEGER
     &  NPRE       !IN. NO OF PRE-SMOOTHING SWEEPS
     &, NPOST      !IN. NO OF POST-SMOOTHING SWEEPS
     &, NCOARSE    !IN. NO OF ITERATIONS OF SMOOTHER ON COARSEST MESH
     &, KREST      !IN. KIND OF RESTRICTION USED:
     &             !     1     - INJECTION
     &             !     2     - HALF-INJECTION
     &             !     3     - FULL-WEIGHTING
     &, NCGC       !IN. NO OF COARSE GRID CORRECTIONS:
     &             !    1     - V-CYCLES
     &             !    2     - W-CYCLES
     &             !    >2    - FULL MG-CYCLES

      REAL
     &  EARTH_RADIUS_INVERSE    !IN.
     &, LONGITUDE_STEP_INVERSE  !IN.
     &, LATITUDE_STEP_INVERSE   !IN.
     &, W                       !IN. RELAXATION PARAMETER FOR EACH
     &                          !    VARIABLE IN SYSTEM.
     &, TOL_RES                 !IN. TOLERANCE FOR RESIDUAL NORM
     &                          !    RELATIVE TO INITIAL RESIDUAL
     &, WORST_SMOOTHING_RATE    !IN. WORST PRACTICAL SMOOTHING RATE
     &                          !    ACCEPTABLE BEFORE CONVERGENCE
     &                          !    OR MAXIMUM ITERATIONS REACHED.

      REAL
     &  A(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, B(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, C1(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, C2(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, DEF(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, D(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, E(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, F(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, G(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, RHS(I_LENGTH,J_LENGTH,K_LENGTH)
     &                          !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_LENGTH,J_LENGTH)
     &                           !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_LENGTH,J_LENGTH)
     &                           !IN. 1./COSINE OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_LENGTH,J_LENGTH)
     &                          !IN. COSINE OF LATITUDE AT B POINTS.
     &, Z_Q(K_LENGTH)           !IN. VALUE OF Z AT Q POINTS
     &, Z_MID(K_LENGTH+1)       !IN. VALUE OF Z AT MID POINTS DEFINED
     &                          !    BETWEEN Q POINTS IN VERTICAL
     &                          !    FIRST VALUE IS BELOW FIRST Q POINT

      REAL
     &  Q(I_LENGTH,J_LENGTH,K_LENGTH)    !INOUT. SOLUTION.

C*----------------------------------------------------------------------

C*L   11  LOCAL WORK ARRAYS REQUIRED.

      INTEGER
     &  NNODES(3)     ! NUMBER OF POINTS IN EACH DIMENSION.
     &, IMAX(MXNGRDS) ! NUMBER OF NODES IN THE I-DIRECTION
     &                ! ON EACH GRID.
     &, JMAX(MXNGRDS) ! NUMBER OF NODES IN THE J-DIRECTION
     &                ! ON EACH GRID.
     &, KMAX(MXNGRDS) ! NUMBER OF NODES IN THE K-DIRECTION
     &                ! ON EACH GRID.
     &, START_ADDRESS(MXNGRDS) ! START ADDRESS IN DATA ARRAY FOR
     &                         ! EACH GRID.
     &, START_ADDRESS_2D(MXNGRDS) ! START ADDRESS IN DATA ARRAY FOR
     &                         ! EACH 2-D GRID.
     &, START_ADDRESS_Z(MXNGRDS) ! START ADDRESS IN DATA ARRAY FOR
     &                         ! THE ARRAYS HOLDING Z VALUES
     &, RES_DIRS(MXNGRDS) ! DIRECTIONS IN WHICH RESTRICTION CAN BE
     &                    ! APPLIED.

      REAL
     &  LATITUDE_STEP_GRIDS(MXNGRDS) ! LATITUDE STEP INVERSE ON EACH
     &                               ! GRID.
     &, LONGITUDE_STEP_GRIDS(MXNGRDS) ! LONGITUDE STEP INVERSE ON EACH
     &                               ! GRID.


C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  NGRIDS        ! THE TOTAL NUMBER OF GRIDS IN THE SEQUENCE
     &, NDIM          !IN. NUMBER OF DIMENSIONS
     &, LAST_ADDRESS  ! LAST ADDRESS OF DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION SPACE REQUIRED.
     &, NFB       ! NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID NEEDED
     &, IPAT      ! AN INTEGER USED TO DISTINGUISH RED AND BLACK POINTS
     &, LAST_ADDRESS_2D ! LAST ADDRESS OF 2-D DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION SPACE REQUIRED.
     &, LAST_ADDRESS_Z ! LAST ADDRESS OF Z DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION SPACE REQUIRED.

      LOGICAL
     &  JAC     ! .TRUE. FOR JACOBI METHODS
     &, PAT     ! .TRUE. FOR PATTERN SCHEMES
     &, SYM     ! .TRUE. FOR SYMMETRIC SCHEMES
     &, ILINE   ! .TRUE. FOR I-LINE METHODS AND ALTERNATING SCHEMES
     &, JLINE   ! .TRUE. FOR J-LINE METHODS AND ALTERNATING SCHEMES
     &, KLINE   ! .TRUE. FOR K-LINE METHODS AND ALTERNATING SCHEMES
     &, I_NT    ! .TRUE. If a non-trivial > 1 number of points in I
     &          !        direction.
     &, J_NT    ! .TRUE. If a non-trivial > 1 number of points in J
     &          !        direction.
     &, K_NT    ! .TRUE. If a non-trivial > 1 number of points in K
     &          !        direction.

      REAL
     &  SWJAC   ! A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM   ! A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS

      INTEGER
     &  I,J,K

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_GRIDSEQ,MG_FDMG,MG_STYPE,MG_JOB_INFO

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. FIND THE NUMBER OF GRIDS WHICH MAY BE GENERATED
CL               SET UP GRID INFORMATION AND FIND SPACE REQUIRED.
CL----------------------------------------------------------------------

      I_NT = .FALSE.
      IF (I_LENGTH.GT.1) I_NT = .TRUE.
      J_NT = .FALSE.
      IF (J_LENGTH.GT.1) J_NT = .TRUE.
      K_NT = .FALSE.
      IF (K_LENGTH.GT.1) K_NT = .TRUE.
      IF (.NOT. K_NT) K_BC = 1

CL    CALL MG_GRIDSEQ

      CALL MG_GRIDSEQ (I_NT,J_NT,K_NT,I_LENGTH,J_LENGTH,K_LENGTH,
     &                 MXNGRDS,NGRIDS,IMAX,JMAX,KMAX,
     &                 START_ADDRESS,LAST_ADDRESS,
     &                 LONGITUDE_STEP_INVERSE,LATITUDE_STEP_INVERSE,
     &                 LONGITUDE_STEP_GRIDS,LATITUDE_STEP_GRIDS,
     &                 START_ADDRESS_2D,LAST_ADDRESS_2D,RES_DIRS,
     &                 VERSION,IPRINT,START_ADDRESS_Z,LAST_ADDRESS_Z)

CL----------------------------------------------------------------------
CL    SECTION 2. SET SMOOTHING INFORMATION DEPENDENT ON CHOICE OF
CL               SMOOTH.
CL----------------------------------------------------------------------

CL    CALL STYPE

      CALL MG_STYPE(I_NT,J_NT,K_NT,KSMOOTH,JAC,PAT,SYM,ILINE,JLINE,
     &              KLINE,NFB,IPAT,SWJAC,SWSYM)

CL----------------------------------------------------------------------
CL    SECTION 3. CALL JOB_INFO TO OUTPUT INFORMATION ON CHOICES MADE
CL               IF DESIRED.
CL----------------------------------------------------------------------

      IF(IPRINT.GE.3) THEN
        CALL MG_JOB_INFO(MAXITS,TOL_RES,KSMOOTH,
     &                   NPRE,NPOST,NCOARSE,W,KREST,NCGC,NGRIDS,IMAX,
     &                   JMAX,KMAX,START_ADDRESS,RES_DIRS,
     &                   WORST_SMOOTHING_RATE,VERSION,K_BC)
      END IF

CL----------------------------------------------------------------------
CL    SECTION 4. CALL FDMG TO PERFORM MULTI-GRID SOLVER.
CL----------------------------------------------------------------------

      CALL MG_FDMG(Q,A,B,C1,C2,DEF,D,E,F,G,RHS,I_NT,J_NT,K_NT,
     &             START_ADDRESS,LAST_ADDRESS,IMAX,JMAX,KMAX,JAC,
     &             PAT,SYM,ILINE,JLINE,KLINE,NFB,IPAT,SWJAC,SWSYM,
     &             I_LENGTH,J_LENGTH,K_LENGTH,NGRIDS,MAXITS,
     &             TOL_RES,IPRINT,KSMOOTH,NPRE,NPOST,NCOARSE,W,
     &             KREST,NCGC,COS_P_LATITUDE,SEC_P_LATITUDE,
     &             COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &             LONGITUDE_STEP_GRIDS,LATITUDE_STEP_GRIDS,
     &             START_ADDRESS_2D,LAST_ADDRESS_2D,RES_DIRS,
     &             WORST_SMOOTHING_RATE,VERSION,K_BC,START_ADDRESS_Z,
     &             LAST_ADDRESS_Z,Z_Q,Z_MID)

CL    END OF ROUTINE MG_CNTL

      RETURN
      END
CLL   SUBROUTINE MG_POINTS
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS POINT SMOOTHING.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_POINTS(
     &                     Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                     SEC_P_LATITUDE,COS_V_LATITUDE,I_DIM,J_DIM,
     &                     K_DIM,W,GRID_NUMBER,RMS_RES,RMS_INC,IPRINT,
     &                     SWJAC,SWSYM,JAC,PAT,SYM,
     &                     NFB,IPAT,LATITUDE_STEP_INVERSE,
     &                     LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                     VERSION,k_bc,Z_Q,Z_MID,I_NT,J_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  IPRINT      !IN. PARAMETER CONTROLLING QUANTITY OF OUTPUT.
     &, I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, GRID_NUMBER !IN. NUMBER OF GRID SMOOTHER IS ACTING ON.
     &, NFB         !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &              ! NEEDED.
     &, IPAT        !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &              ! POINTS.
     &, VERSION     !IN. VERSION OF MULTIGRID BEING USED.
     &, K_BC        !IN. BOUNDARY CONDITION IN K DIRECTION.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES

      REAL
     &  SWJAC    !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM    !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
     &, RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.
     &, W        !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM

      REAL
     &  A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, LONGITUDE_STEP_INVERSE !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   1 LOCAL WORK ARRAY REQUIRED.
      REAL
     &  CORRECTION(I_DIM,J_DIM,K_DIM)

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      REAL
     &  SCALAR
     &, SCALAR1
     &, SCALAR2
     &, WEIGHT
     &, R
     &, X
     &, COEFF
     &, LHS

      INTEGER
     &  I,J,K
     &, I_BEGIN
     &, I_END
     &, I_STEP
     &, J_BEGIN(2)
     &, J_END(2)
     &, J_STEP
     &, K_BEGIN
     &, K_END
     &, K_STEP
     &, I_RB
     &, I_FB
     &, SWEEPS
     &, J_S

C*L   NO EXTERNAL ROUTINES CALLED.

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. POINT RELAXATION SCHEMES
CL----------------------------------------------------------------------

C     Initialise Norms
C     ----------------
      RMS_RES=0.0
      RMS_INC=0.0

      WEIGHT=SWJAC*W

C     Define Point Ordering For First Sweep
C     -------------------------------------

      IF(.NOT.I_NT) THEN
        I_BEGIN = 1
        I_END = 1
      ELSE IF(VERSION.EQ.4) THEN
        I_BEGIN=2
        I_END=I_DIM-1
      ELSE
        I_BEGIN=1
        I_END=I_DIM
      END IF
      I_STEP=1
      SWEEPS=1
      IF (.NOT. J_NT) THEN
        J_BEGIN(1) = 1
        J_END(1) = 1
      ELSE IF (VERSION.EQ.3) THEN
        J_BEGIN(1)=1
        J_END(1)=J_DIM
      ELSE IF(VERSION.EQ.2) THEN
        J_BEGIN(1)=2
        J_END(1)=(J_DIM-1)/2
        J_BEGIN(2)= (J_DIM-1)/2+2
        J_END(2)= J_DIM-1
        SWEEPS = 2
      ELSE
        J_BEGIN(1)=2
        J_END(1)=J_DIM-1
      END IF
      J_STEP=1
      K_BEGIN=1
      K_END=K_DIM
      IF (K_BC.EQ.2 ) THEN
        K_END = K_DIM-1
      ELSE IF (K_BC.EQ.3 ) THEN
        K_BEGIN = 2
      ELSE IF (K_BC.EQ. 4) THEN
        K_BEGIN = 2
        K_END = K_DIM-1
      END IF
      K_STEP=1

C     Set Red-Black Switch For First Sweep
C     ------------------------------------
      I_RB=1

      DO I_FB=1,NFB

C     Sweep Over Grid In Prescribed Point Order
C     -----------------------------------------
        RMS_RES=SWSYM*RMS_RES
        RMS_INC=SWSYM*RMS_INC

        DO K=K_BEGIN,K_END,K_STEP
          DO J_S = 1,SWEEPS
            DO J=J_BEGIN(J_S),J_END(J_S),J_STEP
              DO I=I_BEGIN,I_END,I_STEP

C     Skip To Next Point Of The Same Colour If IPAT=2 (Red-Black)
C     ===========================================================
                IF(MOD(I+J+K+I_RB,IPAT).EQ.0) THEN

C     Evaluate Residual At The Point (I,J,K)
C     ======================================

C NON-DERIVATIVE TERM.

                  LHS =  - Q(I,J,K)*G(I,J,K)
                  COEFF = -G(I,J,K)

                  IF (I_NT ) THEN
C I-DIRECTION DERIVATIVES

                    IF(I.EQ.1) THEN
                      IF(VERSION.LT.3) THEN
                        LHS=LHS+(Q(I_DIM,J,K)*A(I_DIM,J,K) - Q(I,J,K) *
     &                   (A(I_DIM,J,K) + A(I,J,K)) + Q(2,J,K)*A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I_DIM,J,K)*
     &                     (Q(I,J,K)-Q(I_DIM,J,K)))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                        COEFF = COEFF -(A(I_DIM,J,K) + A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I_DIM,J,K)-D(I,J,K))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                      ELSE
                        LHS = LHS + 2.*(-Q(I,J,K) *
     &                      A(I,J,K) + Q(2,J,K)*A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
                        COEFF = COEFF -2.* A(I,J,K)*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
                      END IF
                    ELSE IF(I.EQ.I_DIM) THEN
                      IF(VERSION.LT.3) THEN
                        LHS = LHS +(Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &                   (A(I-1,J,K) + A(I,J,K)) + Q(1,J,K)*A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I,J,K)*(Q(1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &                     (Q(I,J,K)-Q(I-1,J,K)))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                        COEFF = COEFF -(A(I-1,J,K) + A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I-1,J,K)-D(I,J,K))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                      ELSE
                        LHS= LHS+2.*(Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &                      A(I-1,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
                        COEFF = COEFF-2.*A(I-1,J,K)*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
                      END IF
                    ELSE
                      LHS = LHS + (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &                   (A(I-1,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &                     (Q(I,J,K)-Q(I-1,J,K)))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                      COEFF = COEFF -(A(I-1,J,K) + A(I,J,K))*
     &                     EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &                    LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
     &                     *SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
     &                     + .5*DEF(I,J,K)*SEC_P_LATITUDE(I,J)*
     &                    (D(I-1,J,K)-D(I,J,K))*
     &                    LONGITUDE_STEP_INVERSE*EARTH_RADIUS_INVERSE
                    END IF
                  END IF

                  IF( J_NT) THEN
C J-DIRECTION DERIVATIVES

                    SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE
     &                     *LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
                    IF (J.EQ.1) THEN
C ONLY POSSIBLE IN LIMITED AREA CODE, VERSION 3.
                      LHS = LHS -2.*(Q(I,J,K)*B(I,J,K)
     &                    -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
                      COEFF = COEFF -2.*B(I,J,K)*COS_V_LATITUDE(I,J)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
                    ELSE IF(J.EQ.J_DIM) THEN
C ONLY POSSIBLE IN LIMITED AREA CODE, VERSION 3.
                      LHS = LHS + 2.*(Q(I,J-1,K)*B(I,J-1,K)
     &                    -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
                      COEFF = COEFF-2.*B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
                    ELSE
                      LHS = LHS + ((Q(I,J-1,K)*B(I,J-1,K)
     &                    -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                    -(Q(I,J,K)*B(I,J,K)
     &                    -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                 +DEF(I,J,K)*.5*LATITUDE_STEP_INVERSE*
     &                  EARTH_RADIUS_INVERSE*
     &                 (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                  E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
                      COEFF = COEFF - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                   +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                 +DEF(I,J,K)*.5*LATITUDE_STEP_INVERSE*
     &                  EARTH_RADIUS_INVERSE*
     &                 (E(I,J,K)-E(I,J-1,K))
                    END IF
                  END IF

C CALCULATE K-DIRECTION DERIVATIVES AND ADD ON.

                  IF (K_NT) THEN
                    IF(K.EQ.1) THEN
                      SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
                      SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
                      LHS = LHS + SCALAR*(Q(I,J,K+1)*C2(I,J,K+1)
     &                    -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
                      COEFF= COEFF - SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K)
                      IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                        LHS = LHS  + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*
     &                               (Q(I,J,K+1)-Q(I,J,K))
                        COEFF = COEFF - .5*SCALAR1*DEF(I,J,K)*F(I,J,K)
                      END IF
                    ELSE IF (K.EQ.K_DIM) THEN
                      SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
                      SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
                      LHS = LHS - SCALAR*(Q(I,J,K)*C2(I,J,K)
     &                     -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
                      COEFF=COEFF- SCALAR*C1(I,J,K-1)*SCALAR2*C2(I,J,K)
                      IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                        LHS = LHS  + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*
     &                                 (Q(I,J,K)-Q(I,J,K-1))
                        COEFF= COEFF + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)
                      END IF
                    ELSE
                      SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
                      SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
                      SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
                      LHS =  LHS + SCALAR*
     &                   ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                   *C1(I,J,K)*SCALAR1
     &                   -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                   *C1(I,J,K-1)*SCALAR2)
     &                   +.5*DEF(I,J,K)*
     &                   ((Q(I,J,K+1)-Q(I,J,K))
     &                    *F(I,J,K)*SCALAR1
     &                    +(Q(I,J,K)-Q(I,J,K-1))
     &                     *F(I,J,K-1)*SCALAR2)
                      COEFF = COEFF - SCALAR*(C1(I,J,K)*SCALAR1
     &                              +C1(I,J,K-1)*SCALAR2)*C2(I,J,K)
     &                           +.5*DEF(I,J,K)*(F(I,J,K-1)*SCALAR2
     &                               -F(I,J,K)*SCALAR1)
                    END IF
                  END IF

                  R = LHS - RHS(I,J,K)
                  RMS_RES = RMS_RES + R*R

C     SOLVE THE LOCAL PROBLEM: COEFF X = R
C     ===================================
                  X = R / COEFF

C     Store Correction In Workspace
C     And Add To Current Solution Unless WT=0 (Jacobi)
C     ================================================
                  CORRECTION(I,J,K) = X
                  Q(I,J,K) = Q(I,J,K) - WEIGHT * X
                  RMS_INC=RMS_INC+X*X
                END IF
              END DO
            END DO
          END DO
        END DO
C
C     End Of Sweep Over Grid
C     ----------------------
C
C     Add Corrections For The Jacobi Method
C     -------------------------------------
        IF(JAC) THEN
          DO K= K_BEGIN,K_END
            DO J_S = 1,SWEEPS
              DO J=J_BEGIN(J_S),J_END(J_S)
                DO I=I_BEGIN,I_END
                  Q(I,J,K) = Q(I,J,K) - W * CORRECTION(I,J,K)
                END DO
              END DO
            END DO
          END DO
        END IF
C
C     Reset Point Ordering For Second Sweep
C     -------------------------------------
        I = I_BEGIN
        I_BEGIN = I_END
        I_END = I
        I_STEP=-1
        DO J_S =1 ,SWEEPS
          J = J_BEGIN(J_S)
          J_BEGIN(J_S) = J_END(J_S)
          J_END(J_S) = J
        END DO
        J_STEP=-1
        K = K_BEGIN
        K_BEGIN = K_END
        K_END = K
        K_STEP=-1
C
C     Reset Red-Black Switch For Second Sweep
C     ---------------------------------------
        I_RB=0

      END DO

CL    END OF ROUTINE MG_POINTS

      RETURN
      END
CLL   SUBROUTINE MG_POLES.
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS POLAR UPDATES.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_POLES(
     &                    Q,RHS,B,C1,C2,DEF,E,F,G,SEC_P_LATITUDE,
     &                    COS_V_LATITUDE,LATITUDE_STEP_INVERSE,
     &                    EARTH_RADIUS_INVERSE,
     &                    I_DIM,J_DIM,K_DIM,RMS_RES,RMS_INC,
     &                    Z_Q,Z_MID,K_NT,K_BC)

      IMPLICIT NONE

      LOGICAL
     &  K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, K_BC        !IN. BOUNDARY CONDITION IN K DIRECTION.

      REAL
     &  RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.

      REAL
     &  B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   NO LOCAL WORK ARRAYS REQUIRED.

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,J,K,I2
     &, K_BEGIN
     &, K_END

      REAL
     &  X_POLE
     &, SCALAR
     &, SCALARS
     &, SCALAR1
     &, SCALAR2
     &, SCALAR3
     &, LHS


C*L   NO EXTERNAL ROUTINES CALLED.

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. FIND CORRECTION AT POLES.
CL----------------------------------------------------------------------

      K_BEGIN = 1
      K_END = K_DIM
      IF (K_BC.EQ.2 ) THEN
        K_END = K_DIM-1
      ELSE IF (K_BC.EQ.3 ) THEN
        K_BEGIN = 2
      ELSE IF (K_BC.EQ. 4) THEN
        K_BEGIN = 2
        K_END = K_DIM-1
      END IF

      DO K= K_BEGIN,K_END
CL    NORTH POLE

        I2 = I_DIM/2
        SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &           LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
        SCALARS= EARTH_RADIUS_INVERSE*LATITUDE_STEP_INVERSE
        LHS = 0.
        X_POLE = 0.

C SUM ALL TERMS AT POLE.
        DO I= 1,I_DIM/2
          LHS = LHS + ((Q(I,2,K)-Q(I,1,K))
     &                       *B(I,1,K)*COS_V_LATITUDE(I,1)
     &                        -(Q(I+I2,1,K) -Q(I+I2,2,K))
     &                        *B(I+I2,1,K)*COS_V_LATITUDE(I+I2,1))
     &                        *SEC_P_LATITUDE(I,1)*SCALAR
     &                    +DEF(1,1,K)*SCALARS*
     &                    (E(I,1,K)*(Q(I,1,K)-Q(I,2,K))+
     &                     E(I+I2,1,K)*(Q(I+I2,1,K)-Q(I+I2,2,K)))
          X_POLE = X_POLE - (B(I,1,K)*COS_V_LATITUDE(I,1)
     &                    +B(I+I2,1,K)*COS_V_LATITUDE(I+I2,1))
     &                     *SEC_P_LATITUDE(I,1)*SCALAR
     &                    +DEF(1,1,K)*SCALARS*
     &                    (E(I,1,K)+E(I+I2,1,K))
        END DO

        LHS = LHS - Q(1,1,K)*G(1,1,K)
        X_POLE = X_POLE - G(1,1,K)

        IF(K_NT) THEN
CL INCLUDE VERTICAL TERMS IN 3-D.
          IF(K.EQ.1) THEN
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            LHS = LHS + SCALAR3*(Q(1,1,K+1)*C2(1,1,K+1)
     &                     -Q(1,1,K)*C2(1,1,K))*C1(1,1,K)*SCALAR1
            X_POLE = X_POLE - SCALAR3*C1(1,1,K)*SCALAR1*C2(1,1,K)
            IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS = LHS  + .5*SCALAR1*DEF(1,1,K)*F(1,1,K)*
     &               (Q(1,1,K+1)*C2(1,1,K+1)
     &                -Q(1,1,K)*C2(1,1,K))
              X_POLE = X_POLE - .5*SCALAR1*DEF(1,1,K)*F(1,1,K)
     &                          *C2(1,1,K)
            END IF
          ELSE IF (K.EQ.K_DIM) THEN
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            LHS = LHS - SCALAR3*(Q(1,1,K)*C2(1,1,K)
     &                   -Q(1,1,K-1)*C2(1,1,K-1))*C1(1,1,K-1)*SCALAR2
            X_POLE = X_POLE - SCALAR3*C1(1,1,K-1)*SCALAR2*C2(1,1,K)
            IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS = LHS  + .5*SCALAR2*DEF(1,1,K)*F(1,1,K-1)*
     &                 (Q(1,1,K)*C2(1,1,K)
     &                  - Q(1,1,K-1)*C2(1,1,K-1))
              X_POLE = X_POLE + .5*SCALAR2*DEF(1,1,K)*F(1,1,K-1)
     &                          *C2(1,1,K)
            END IF
          ELSE
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            LHS = LHS + SCALAR3*
     &               ((Q(1,1,K+1)*C2(1,1,K+1)-Q(1,1,K)*C2(1,1,K))
     &                *C1(1,1,K)*SCALAR1
     &               -(Q(1,1,K)*C2(1,1,K)-Q(1,1,K-1)*C2(1,1,K-1))
     &                *C1(1,1,K-1)*SCALAR2)
     &               +.5*DEF(1,1,K)*
     &               ((Q(1,1,K+1)*C2(1,1,K+1)
     &                -Q(1,1,K)*C2(1,1,K))
     &                *F(1,1,K)*SCALAR1
     &               +(Q(1,1,K)*C2(1,1,K)
     &                 -Q(1,1,K-1)*C2(1,1,K-1))
     &                *F(1,1,K-1)*SCALAR2)
            X_POLE = X_POLE - SCALAR3*(C1(1,1,K)*SCALAR1
     &                              +C1(1,1,K-1)*SCALAR2)*C2(1,1,K)
     &               +.5*DEF(1,1,K)*(F(1,1,K-1)*SCALAR2
     &                               -F(1,1,K)*SCALAR1)*C2(1,1,K)
          END IF
        END IF

C AVERAGE POLAR VALUE AND INCLUDE CORRECTIONS.
        X_POLE = (LHS - RHS(1,1,K)*I_DIM) / X_POLE
        DO I= 1,I_DIM
          Q(I,1,K)=Q(I,1,K)-X_POLE
          RMS_INC=RMS_INC+X_POLE*X_POLE
        END DO

CL    SOUTH POLE

        LHS = 0.
        X_POLE = 0.

C SUM ALL TERMS AT POLE.
        DO I= 1,I_DIM/2
          LHS = LHS +((Q(I,J_DIM-1,K) -Q(I,J_DIM,K))
     &                      *B(I,J_DIM-1,K)*COS_V_LATITUDE(I,J_DIM-1)
     &                      -(Q(I+I2,J_DIM,K)-Q(I+I2,J_DIM-1,K))
     &                      *B(I+I2,J_DIM-1,K)*
     &                      COS_V_LATITUDE(I+I2,J_DIM-1))
     &                      *SEC_P_LATITUDE(I,J_DIM)*SCALAR
     &                    +DEF(1,J_DIM,K)*SCALARS*
     &                (E(I,J_DIM-1,K)*(Q(I,J_DIM-1,K)-Q(I,J_DIM,K))+
     &                     E(I+I2,J_DIM-1,K)*(Q(I+I2,J_DIM-1,K)
     &                                         -Q(I+I2,J_DIM,K)))
          X_POLE = X_POLE - (B(I,J_DIM-1,K)*
     &                      COS_V_LATITUDE(I,J_DIM-1)
     &                      +B(I+I2,J_DIM-1,K)*
     &                      COS_V_LATITUDE(I+I2,J_DIM-1))
     &                      *SEC_P_LATITUDE(I,J_DIM)*SCALAR
     &                      -DEF(1,J_DIM,K)*SCALARS*
     &                      (E(I,J_DIM-1,K)+E(I+I2,J_DIM-1,K))
        END DO
        LHS = LHS - Q(1,J_DIM,K)*G(1,J_DIM,K)
        X_POLE = X_POLE - G(1,J_DIM,K)

        IF(K_NT) THEN
CL INCLUDE VERTICAL TERMS.
          IF(K.EQ.1) THEN
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            LHS = LHS + SCALAR3*(Q(1,J_DIM,K+1)*C2(1,J_DIM,K+1)
     &               -Q(1,J_DIM,K)*C2(1,J_DIM,K))*C1(1,J_DIM,K)*SCALAR1
            X_POLE=X_POLE - SCALAR3*C1(1,J_DIM,K)*SCALAR1*C2(1,J_DIM,K)
            IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS = LHS  + .5*SCALAR1*DEF(1,J_DIM,K)*F(1,J_DIM,K)*
     &                      (Q(1,J_DIM,K+1)*C2(1,J_DIM,K+1)-
     &                       Q(1,J_DIM,K)*C2(1,J_DIM,K))
              X_POLE = X_POLE - .5*SCALAR1*DEF(1,J_DIM,K)*F(1,J_DIM,K)
     &                          *C2(1,J_DIM,K)
            END IF
          ELSE IF (K.EQ.K_DIM) THEN
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            LHS = LHS - SCALAR3*(Q(1,J_DIM,K)*C2(1,J_DIM,K)
     &          -Q(1,J_DIM,K-1)*C2(1,J_DIM,K-1))*C1(1,J_DIM,K-1)*SCALAR2
            X_POLE=X_POLE- SCALAR3*C1(1,J_DIM,K-1)*SCALAR2*C2(1,J_DIM,K)
            IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS = LHS  + .5*SCALAR2*DEF(1,J_DIM,K)*F(1,J_DIM,K-1)*
     &                      (Q(1,J_DIM,K)*C2(1,J_DIM,K)-
     &                       Q(1,J_DIM,K-1)*C2(1,J_DIM,K-1))
              X_POLE = X_POLE + .5*SCALAR2*DEF(1,J_DIM,K)*F(1,J_DIM,K-1)
     &                          *C2(1,J_DIM,K)
            END IF
          ELSE
            SCALAR3 = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            LHS = LHS + SCALAR3*
     &               ((Q(1,J_DIM,K+1)*C2(1,J_DIM,K+1)
     &                 -Q(1,J_DIM,K)*C2(1,J_DIM,K))
     &                *C1(1,J_DIM,K)*SCALAR1
     &               -(Q(1,J_DIM,K)*C2(1,J_DIM,K)
     &                 -Q(1,J_DIM,K-1)*C2(1,J_DIM,K-1))
     &                *C1(1,J_DIM,K-1)*SCALAR2)
     &               +.5*DEF(1,J_DIM,K)*
     &               ((Q(1,J_DIM,K+1)*C2(1,J_DIM,K+1)
     &                -Q(1,J_DIM,K)*C2(1,J_DIM,K))
     &                *F(1,J_DIM,K)*SCALAR1
     &               + (Q(1,J_DIM,K)*C2(1,J_DIM,K)
     &                 -Q(1,J_DIM,K-1)*C2(1,J_DIM,K-1))
     &                *F(1,J_DIM,K-1)*SCALAR2)
            X_POLE = X_POLE - SCALAR3*(C1(1,J_DIM,K)*SCALAR1
     &                        +C1(1,J_DIM,K-1)*SCALAR2)*C2(1,J_DIM,K)
     &               +.5*DEF(1,J_DIM,K)*(F(1,J_DIM,K-1)*SCALAR2
     &                          -F(1,J_DIM,K)*SCALAR1)*C2(1,J_DIM,K)
          END IF
        END IF

C AVERAGE POLAR VALUE AND INCLUDE CORRECTIONS.
        X_POLE = (LHS - RHS(1,J_DIM,K)*I_DIM) / X_POLE
        DO I= 1,I_DIM
          Q(I,J_DIM,K)=Q(I,J_DIM,K)-X_POLE
          RMS_INC=RMS_INC+X_POLE*X_POLE
        END DO
      END DO

CL    END OF ROUTINE MG_POLES

      RETURN
      END
CLL   SUBROUTINE MG_PROLONG
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALCULATES CORRECTION TO SOLUTION AT FINE GRID POINTS
CLL   CORRESPONDING TO THE COARSE GRID ONES. THIS CORRECTION IS THEN
CLL   INTERPOLATED TO ALL THE FINE GRID POINTS.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_PROLONG(
     &                      Q_COARSE,Q_FINE,I_COARSE,J_COARSE,K_COARSE,
     &                      I_FINE,J_FINE,K_FINE,RES_DIRS,VERSION,K_BC,
     &                      Z_Q_FINE,I_NT,J_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_COARSE   !IN. NUMBER OF POINTS IN I DIRECTION ON COARSE GRID.
     &, J_COARSE   !IN. NUMBER OF POINTS IN J DIRECTION ON COARSE GRID.
     &, K_COARSE   !IN. NUMBER OF POINTS IN K DIRECTION ON COARSE GRID.
     &, I_FINE     !IN. NUMBER OF POINTS IN I DIRECTION ON FINE GRID.
     &, J_FINE     !IN. NUMBER OF POINTS IN J DIRECTION ON FINE GRID.
     &, K_FINE     !IN. NUMBER OF POINTS IN K DIRECTION ON FINE GRID.
     &, RES_DIRS   !IN. RESTRICTED DIRECTIONS FROM FINE GRID TO COARSE
     &             !    = 3 NEED TO PROLONG IN ALL DIRECTIONS
     &             !    = 2 NEED TO PROLONG IN I AND J DIRECTIONS
     &             !    = 1 NEED TO PROLONG IN I DIRECTION
     &, VERSION    !IN. VERSION OF CODE TO BE USED.
     &, K_BC       !IN. BOUNDARY CONDITION IN K DIRECTION.

      REAL
     &  Q_COARSE(I_COARSE,J_COARSE,K_COARSE) !IN. COARSE GRID SOLUTION.
     &, Z_Q_FINE(K_FINE)    !IN. Z VALUES AT Q_FINE VERTICAL POINTS.

      REAL
     &  Q_FINE(I_FINE,J_FINE,K_FINE)         !INOUT. FINE GRID SOLUTION.

C*----------------------------------------------------------------------

C*L     LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  Q_CORRECTION(I_FINE,J_FINE,K_COARSE) !USED TO SAVE CALCULATIONS
     &                                       !IN 3-D CODE.
     &, Q_FINE_BOUNDARY(I_FINE*2+J_FINE*2,K_FINE)
     &              !HOLDS BOUNDARY ELEMENTS IN VERSION 4

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,J,K                 ! LOOP COUNTERS
     &, II,JJ,KK,JP,JJP
     &, KC_BEGIN,KC_END

      LOGICAL
     &  J_RES
     &, SWITCH

      REAL
     &  FRAC_ABOVE,
     &  FRAC_BELOW

C*L   NO EXTERNAL ROUTINES CALLED.

C*----------------------------------------------------------------------

CFPP$ NOCONCUR R

CL----------------------------------------------------------------------
CL    SECTION 1. WRITE OUT LIMITED AREA DIRICHLET BOUNDARY CONDITION
CL               IF REQUIRED.
CL               FIND COARSE GRID CORRECTION BY SUBTRACTING
CL               FINE GRID SOLUTION POINTWISE.
CL----------------------------------------------------------------------

C SET KC VARIABLES DEPENDING ON BOUNDARY CONDITIONS IN K DIRECTIONS.
C FOR DIRICHLET BOUNDARY CONDITIONS NO PROLONGATION IS REQUIRED ON
C THAT COARSE LEVEL.
      KC_BEGIN = 1
      KC_END = K_COARSE
      IF (K_BC.EQ.2 ) THEN
        KC_END = K_COARSE - 1
      ELSE IF (K_BC.EQ.3) THEN
        KC_BEGIN = 2
      ELSE IF (K_BC.EQ.4) THEN
        KC_BEGIN = 2
        KC_END = K_COARSE - 1
      END IF

C-----------------------------------------------------------------------
CL    SECTION 1.1 WRITE OUT LIMITED AREA DIRICHLET BOUNDARY CONDITION
CL                IF REQUIRED.
C-----------------------------------------------------------------------

      IF (VERSION.EQ.4) THEN
        DO K=1,K_FINE
          IF(I_NT .AND. J_NT) THEN
            DO I=1,I_FINE
              Q_FINE_BOUNDARY(I,K) = Q_FINE(I,1,K)
            END DO
            II = I_FINE
            DO J=2,J_FINE-1
              II = II+1
              Q_FINE_BOUNDARY(II,K) = Q_FINE(1,J,K)
              II = II+1
              Q_FINE_BOUNDARY(II,K) = Q_FINE(I_FINE,J,K)
            END DO
            DO I=1,I_FINE
              Q_FINE_BOUNDARY(II+I,K) = Q_FINE(I,J_FINE,K)
            END DO
          ELSE IF (I_NT) THEN
            Q_FINE_BOUNDARY(1,K) = Q_FINE(1,1,K)
            Q_FINE_BOUNDARY(2,K) = Q_FINE(I_FINE,1,K)
          ELSE IF (J_NT) THEN
            Q_FINE_BOUNDARY(1,K) = Q_FINE(1,1,K)
            Q_FINE_BOUNDARY(2,K) = Q_FINE(1,J_FINE,K)
          END IF
        END DO
      END IF

C-----------------------------------------------------------------------
CL    SECTION 1.2 FIND COARSE GRID CORRECTION BY SUBTRACTING
CL                FINE GRID SOLUTION POINTWISE.
C-----------------------------------------------------------------------

      J_RES = RES_DIRS.GT.5 .OR. RES_DIRS.EQ.3 .OR. RES_DIRS.EQ.2
      SWITCH = RES_DIRS.EQ.1

      IF(RES_DIRS.GT.5) THEN

        DO K = KC_BEGIN,KC_END
          IF(RES_DIRS.EQ.7) THEN
            KK=2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_COARSE
            JJ=2*J-1
            DO I = 1,I_COARSE
              II=2*I-1
              Q_COARSE(I,J,K) = Q_COARSE(I,J,K) - Q_FINE(II,JJ,KK)
            END DO
          END DO
        END DO

      ELSE

        DO K = KC_BEGIN,KC_END
          IF(MOD(RES_DIRS,2).EQ.1) THEN
            KK = 2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_COARSE
            IF(J_RES) THEN
              JJ = 2*J-1
            ELSE
              JJ = J
            END IF
            DO I = 1,I_COARSE
              IF(RES_DIRS.GT.3) THEN
                II = 2*I-1
              ELSE
                II = I
              END IF
              Q_COARSE(I,J,K) = Q_COARSE(I,J,K) - Q_FINE(II,JJ,KK)
C SET CORRECTION FIELD FOR K ONLY PROLONGATION.
              IF (SWITCH) THEN
                Q_FINE(II,JJ,KK) = Q_COARSE(I,J,K)
                Q_CORRECTION(II,JJ,K) = Q_COARSE(I,J,K)
              END IF
            END DO
          END DO
        END DO

      END IF

CL----------------------------------------------------------------------
CL    SECTION 2. ADD PROLONGED CORRECTION TO FINE GRID SOLUTION
CL----------------------------------------------------------------------

C-----------------------------------------------------------------------
CL    SECTION 2.1 HORIZONTAL PROLONGATION.
CL                CORRECTION ON EACH COARSE GRID LEVEL IS STORED FOR
CL                USE IN VERTICAL PROLONGATION.
C-----------------------------------------------------------------------

      IF(RES_DIRS.GT.5) THEN
C-----------------------------------------------------------------------
CL    I & J PROLONGATION, PLUS POSSIBLY K.
C-----------------------------------------------------------------------

CL    FILL IN ALL POINTS ON FINE GRID (I,J) ON COARSE GRID K LEVELS.

        DO K= KC_BEGIN,KC_END
          IF(RES_DIRS.EQ.7) THEN
            KK=2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_COARSE
            IF(.NOT.(J.EQ.(J_COARSE+1)/2 .AND. VERSION.EQ.2)) THEN
              JJ=2*J-1
CL    DO ALL I POINTS ON J_COARSE ROWS.
              DO I = 1,I_COARSE-1
                II=2*I-1
                Q_CORRECTION(II,JJ,K) = Q_COARSE(I,J,K)
                Q_FINE(II,JJ,KK) = Q_FINE(II,JJ,KK) +
     &                             Q_CORRECTION(II,JJ,K)
                Q_CORRECTION(I*2,JJ,K) = .5*
     &                             (Q_COARSE(I,J,K)+Q_COARSE(I+1,J,K))
                Q_FINE(I*2,JJ,KK) = Q_FINE(I*2,JJ,KK) +
     &                              Q_CORRECTION(I*2,JJ,K)
              END DO
              IF(VERSION.LT.3) THEN
C GLOBAL, 2 EXTRA POINTS TO DO.
                Q_CORRECTION(I_FINE-1,JJ,K) = Q_COARSE(I_COARSE,J,K)
                Q_FINE(I_FINE-1,JJ,KK) = Q_FINE(I_FINE-1,JJ,KK)+
     &                                   Q_CORRECTION(I_FINE-1,JJ,K)
                Q_CORRECTION(I_FINE,JJ,K) = .5*(Q_COARSE(I_COARSE,J,K)+
     &                                          Q_COARSE(1,J,K))
                Q_FINE(I_FINE,JJ,KK) = Q_FINE(I_FINE,JJ,KK) +
     &                                 Q_CORRECTION(I_FINE,JJ,K)
              ELSE
C LIMITED AREA, 1 EXTRA POINT TO DO.
                Q_CORRECTION(I_FINE,JJ,K) = Q_COARSE(I_COARSE,J,K)
                Q_FINE(I_FINE,JJ,KK) = Q_FINE(I_FINE,JJ,KK)+
     &                                 Q_CORRECTION(I_FINE,JJ,K)
              END IF
            END IF
          END DO
CL    FILL IN MISSING ROWS
          DO J = 2,J_FINE-1,2
            JJ = J/2
            DO I = 1,I_FINE,2
              II = (I+1)/2
              Q_CORRECTION(I,J,K) = .5* (Q_COARSE(II,JJ,K) +
     &                                   Q_COARSE(II,JJ+1,K))
              Q_FINE(I,J,KK) = Q_FINE(I,J,KK) + Q_CORRECTION(I,J,K)
            END DO
            DO I = 2,I_FINE-1,2
              II = I/2
              Q_CORRECTION(I,J,K) = .25* (Q_COARSE(II,JJ,K) +
     &                                    Q_COARSE(II,JJ+1,K) +
     &                                    Q_COARSE(II+1,JJ,K) +
     &                                    Q_COARSE(II+1,JJ+1,K))
              Q_FINE(I,J,KK) = Q_FINE(I,J,KK) + Q_CORRECTION(I,J,K)
            END DO
            IF (VERSION.LT.3) THEN
              Q_CORRECTION(I_FINE,J,K) = .25* (Q_COARSE(1,JJ,K) +
     &                                         Q_COARSE(1,JJ+1,K) +
     &                                         Q_COARSE(I_COARSE,JJ,K) +
     &                                        Q_COARSE(I_COARSE,JJ+1,K))
              Q_FINE(I_FINE,J,KK) = Q_FINE(I_FINE,J,KK) +
     &                              Q_CORRECTION(I_FINE,J,K)
            END IF
          END DO
        END DO

      ELSE IF (RES_DIRS.GT.3) THEN
C-----------------------------------------------------------------------
CL    I PROLONGATION, PLUS POSSIBLY K.
C-----------------------------------------------------------------------

        DO K= KC_BEGIN,KC_END
          IF(RES_DIRS.EQ.5) THEN
            KK=2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_COARSE
            IF(.NOT.(J.EQ.(J_COARSE+1)/2 .AND. VERSION.EQ.2)) THEN
              JJ=2*J-1
CL    DO ALL I POINTS ON J_COARSE ROWS.
              DO I = 1,I_COARSE-1
                II=2*I-1
                Q_CORRECTION(II,JJ,K) = Q_COARSE(I,J,K)
                Q_FINE(II,JJ,KK) = Q_FINE(II,JJ,KK) +
     &                             Q_CORRECTION(II,JJ,K)
                Q_CORRECTION(I*2,JJ,K) = .5*
     &                             (Q_COARSE(I,J,K)+Q_COARSE(I+1,J,K))
                Q_FINE(I*2,JJ,KK) = Q_FINE(I*2,JJ,KK) +
     &                              Q_CORRECTION(I*2,JJ,K)
              END DO
              IF(VERSION.LT.3) THEN
C GLOBAL, 2 EXTRA POINTS TO DO.
                Q_CORRECTION(I_FINE-1,JJ,K) = Q_COARSE(I_COARSE,J,K)
                Q_FINE(I_FINE-1,JJ,KK) = Q_FINE(I_FINE-1,JJ,KK)+
     &                                   Q_CORRECTION(I_FINE-1,JJ,K)
                Q_CORRECTION(I_FINE,JJ,K) = .5*(Q_COARSE(I_COARSE,J,K)+
     &                                          Q_COARSE(1,J,K))
                Q_FINE(I_FINE,JJ,KK) = Q_FINE(I_FINE,JJ,KK) +
     &                                 Q_CORRECTION(I_FINE,JJ,K)
              ELSE
C LIMITED AREA, 1 EXTRA POINT TO DO.
                Q_CORRECTION(I_FINE,JJ,K) = Q_COARSE(I_COARSE,J,K)
                Q_FINE(I_FINE,JJ,KK) = Q_FINE(I_FINE,JJ,KK)+
     &                                 Q_CORRECTION(I_FINE,JJ,K)
              END IF
            END IF
          END DO
        END DO

      ELSE IF (RES_DIRS.GT.1) THEN
C-----------------------------------------------------------------------
CL    J PROLONGATION, PLUS POSSIBLY K.
C-----------------------------------------------------------------------

CL    FILL IN ALL POINTS ON FINE GRID (I,J) ON COARSE GRID K LEVELS.

        DO K= KC_BEGIN,KC_END
          IF(RES_DIRS.EQ.3) THEN
            KK=2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_FINE,2
            IF(.NOT.(J.EQ.(J_FINE+1)/2 .AND. VERSION.EQ.2)) THEN
              JJ=(J+1)/2
              DO I = 1,I_FINE
                Q_CORRECTION(I,J,K) = Q_COARSE(I,JJ,K)
                Q_FINE(I,J,KK) = Q_FINE(I,J,KK) + Q_CORRECTION(I,J,K)
              END DO
            END IF
            IF (J.NE.J_FINE) THEN
              JP = J+1
              JJP= JP/2
              DO I = 1,I_FINE
                Q_CORRECTION(I,JP,K) = .5* (Q_COARSE(I,JJP,K) +
     &                                     Q_COARSE(I,JJP+1,K))
                Q_FINE(I,JP,KK) = Q_FINE(I,JP,KK) + Q_CORRECTION(I,JP,K)
              END DO
            END IF
          END DO
        END DO
      END IF

C-----------------------------------------------------------------------
CL    SECTION 2.2 VERTICAL PROLONGATION.
CL                FILLS IN ANY FINE GRID LEVELS NOT ALREADY SET BY
CL                PROLONGATION ON THE COARSE GRID LEVELS.
CL                USES CORRECTION STORED ON EVERY COARSE GRID LEVEL.
C-----------------------------------------------------------------------

      IF (MOD(RES_DIRS,2).EQ.1) THEN
CL    FILL IN ALL POINTS ON FINE GRID (I,J) ON FINE GRID K LEVELS.
CL    USE LINEAR INTERPOLATION BUT WEIGHTED IN CASE NEW FINE LEVEL IS
CL    NOT MID-WAY BETWEEN COARSE LEVELS.

        IF (KC_BEGIN.NE.1) THEN
          DO J= 1,J_FINE
            DO I=1,I_FINE
              Q_CORRECTION(I,J,1) = 0.
            END DO
          END DO
        END IF

        IF (KC_END.NE.K_COARSE) THEN
          DO J= 1,J_FINE
            DO I=1,I_FINE
              Q_CORRECTION(I,J,K_COARSE) = 0.
            END DO
          END DO
        END IF

        DO K= 2,K_FINE-1,2
          KK = K/2
          FRAC_ABOVE = (Z_Q_FINE(K) - Z_Q_FINE(K-1))
     &                 /(Z_Q_FINE(K+1)-Z_Q_FINE(K-1))
          FRAC_BELOW = (Z_Q_FINE(K+1) - Z_Q_FINE(K))
     &                 /(Z_Q_FINE(K+1)-Z_Q_FINE(K-1))
C ALL POINTS.
          DO J = 1,J_FINE
            IF(.NOT.(J.EQ.(J_FINE+1)/2 .AND. VERSION.EQ.2)) THEN
              DO I = 1,I_FINE
                Q_FINE(I,J,K) = Q_FINE(I,J,K) + FRAC_ABOVE *
     &                          Q_CORRECTION(I,J,KK+1) + FRAC_BELOW
     &                          * Q_CORRECTION(I,J,KK)
              END DO
            END IF
          END DO
        END DO
      END IF

C-----------------------------------------------------------------------
CL    SECTION 2.3 LIMITED AREA DIRICHLET BOUNDARY CONDITION ENFORCED.
CL                USES VALUES WRITTEN OUT IN SECTION 1.1
C-----------------------------------------------------------------------

      IF (VERSION.EQ.4) THEN
CL OVER-WRITE BOUNDARY VALUES ON FINE GRID WITH SAVED ONES
        DO K=1,K_FINE
          IF(I_NT .AND. J_NT) THEN
            DO I=1,I_FINE
              Q_FINE(I,1,K) = Q_FINE_BOUNDARY(I,K)
            END DO
            II = I_FINE
            DO J=2,J_FINE-1
              II = II+1
              Q_FINE(1,J,K) = Q_FINE_BOUNDARY(II,K)
              II = II+1
              Q_FINE(I_FINE,J,K) = Q_FINE_BOUNDARY(II,K)
            END DO
            DO I=1,I_FINE
              Q_FINE(I,J_FINE,K) = Q_FINE_BOUNDARY(II+I,K)
            END DO
          ELSE IF (I_NT) THEN
            Q_FINE(1,1,K) = Q_FINE_BOUNDARY(1,K)
            Q_FINE(I_FINE,1,K) = Q_FINE_BOUNDARY(2,K)
          ELSE IF (J_NT) THEN
            Q_FINE(1,1,K) = Q_FINE_BOUNDARY(1,K)
            Q_FINE(1,J_FINE,K) = Q_FINE_BOUNDARY(2,K)
          END IF
        END DO
      END IF

CL    END OF ROUTINE MG_PROLONG

      RETURN
      END
CLL   SUBROUTINE MG_RESTRICT
CLL
CLL   PURPOSE:
CLL   -------
CLL   RESTRICTS SOLUTION AND RESIDUAL FROM FINE GRID TO COARSER ONE.
CLL   IF ON THE FIRST DESCENT OF THE MULTI-GRID ALGORITHM THEN
CLL   RESTRICTS OTHER FIELDS AS WELL.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_RESTRICT(
     &                       Q_COARSE,Q_FINE,RHS_COARSE,RESID_FINE,
     &                       A_COARSE,A_FINE,B_COARSE,B_FINE,
     &                       C1_COARSE,C1_FINE,C2_COARSE,C2_FINE,
     &                       DEF_COARSE,DEF_FINE,D_COARSE,D_FINE,
     &                       E_COARSE,E_FINE,F_COARSE,F_FINE,
     &                       G_COARSE,G_FINE,COS_P_COARSE,COS_P_FINE,
     &                       SEC_P_COARSE,COS_V_COARSE,COS_V_FINE,
     &                       Z_Q_COARSE,Z_Q_FINE,
     &                       Z_MID_COARSE,Z_MID_FINE,
     &                       L_FIRST_DESCENT,
     &                       I_COARSE,J_COARSE,K_COARSE,
     &                       I_FINE,J_FINE,K_FINE,I_NT,J_NT,K_NT,KREST,
     &                       EARTH_RADIUS_INVERSE,LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,RES_DIRS,VERSION)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_COARSE   !IN. NUMBER OF POINTS IN I DIRECTION ON COARSE GRID.
     &, J_COARSE   !IN. NUMBER OF POINTS IN J DIRECTION ON COARSE GRID.
     &, K_COARSE   !IN. NUMBER OF POINTS IN K DIRECTION ON COARSE GRID.
     &, I_FINE     !IN. NUMBER OF POINTS IN I DIRECTION ON FINE GRID.
     &, J_FINE     !IN. NUMBER OF POINTS IN J DIRECTION ON FINE GRID.
     &, K_FINE     !IN. NUMBER OF POINTS IN K DIRECTION ON FINE GRID.
     &, KREST      !IN. KIND OF RESTRICTION USED.
     &, RES_DIRS   !IN. DIRECTIONS TO RESTRICT IN.
     &, VERSION    !IN. VERSION OF SCHEME.

      REAL
     &  EARTH_RADIUS_INVERSE    !IN. 1/EARTH RADIUS
     &, LATITUDE_STEP_INVERSE   !IN. 1/LATITUDE STEP LENGTH IN RADIANS
     &                          !    ON COARSE GRID.
     &, LONGITUDE_STEP_INVERSE  !IN. 1/LONGITUDE STEP LENGTH IN RADIANS
     &                          !    ON COARSE GRID.

      LOGICAL
     &  L_FIRST_DESCENT !IN. IF TRUE THEN ALL FIELDS NEED RESTRICTING.

      REAL
     &  A_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, B_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, C1_FINE(I_FINE,J_FINE,K_FINE)   !IN. COEFFICIENT
     &, C2_FINE(I_FINE,J_FINE,K_FINE)   !IN. COEFFICIENT
     &, DEF_FINE(I_FINE,J_FINE,K_FINE)  !IN. COEFFICIENT
     &, D_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, E_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, F_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, G_FINE(I_FINE,J_FINE,K_FINE)    !IN. COEFFICIENT
     &, COS_P_FINE(I_FINE,J_FINE)       !IN. COSINE OF LATITUDE
     &, COS_V_FINE(I_FINE,J_FINE-1)     !IN. COSINE OF LATITUDE
     &, Q_FINE(I_FINE,J_FINE,K_FINE)    !IN. SOLUTION
     &, RESID_FINE(I_FINE,J_FINE,K_FINE)!IN. RESIDUAL
     &, Z_Q_COARSE(K_COARSE)            !IN. Z AT Q COARSE POINTS.
     &, Z_Q_FINE(K_FINE)                !IN. Z AT Q FINE POINTS.
     &, Z_MID_COARSE(0:K_COARSE)        !IN. Z AT MID Q COARSE POINTS.
     &, Z_MID_FINE(0:K_FINE)            !IN. Z AT MID Q FINE POINTS.

      REAL
     &  RHS_COARSE(I_COARSE,J_COARSE,K_COARSE) !OUT. RIGHT-HAND-SIDE
     &, Q_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT. SOLUTION.

CL    FOLLOWING ARE PASSED OUT ONLY IF L_FIRST_DESCENT=.TRUE. AS
CL    THEY HAVE ALREADY BEEN SET IF IT IS .FALSE.

      REAL
     &  COS_P_COARSE(I_COARSE,J_COARSE)   !OUT. COSINE OF LATITUDE
     &, SEC_P_COARSE(I_COARSE,J_COARSE)   !OUT. 1./COSINE OF LATITUDE
     &, COS_V_COARSE(I_COARSE,J_COARSE-1) !OUT. COSINE OF LATITUDE AT
     &                                    !     MID-POINTS.
     &, A_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT. COEFFICIENT
     &, B_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT.    "
     &, C1_COARSE(I_COARSE,J_COARSE,K_COARSE)  !OUT.    "
     &, C2_COARSE(I_COARSE,J_COARSE,K_COARSE)  !OUT.    "
     &, DEF_COARSE(I_COARSE,J_COARSE,K_COARSE) !OUT.    "
     &, D_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT.    "
     &, E_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT.    "
     &, F_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT.    "
     &, G_COARSE(I_COARSE,J_COARSE,K_COARSE)   !OUT.    "
C*----------------------------------------------------------------------

C*L  1  LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  LEVEL_INCS(I_COARSE,J_COARSE,K_FINE)

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,J,K                 ! LOOP COUNTERS
     &, II,JJ,KK,JB,JE

      REAL
     &  SCALE
     &, SCALAR1
     &, SCALAR2
     &, SCALAR3
     &, SCALAR4
     &, FRAC_ABOVE
     &, FRAC_BELOW

      LOGICAL
     &  J_RES

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_LEFT_HAND_SIDE

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. RESTRICTION OF COEFFICIENTS AND TRIGONOMETRICAL TERMS
CL               IF REQUIRED.
CL----------------------------------------------------------------------

      J_RES = RES_DIRS.GT.5 .OR. RES_DIRS.EQ.2 .OR. RES_DIRS.EQ.3
      IF(.NOT. J_NT) THEN
        JB = 1
        JE = 1
      ELSE IF(VERSION.LT.3) THEN
        JB = 2
        JE = J_COARSE - 1
      ELSE
        JB = 1
        JE = J_COARSE
      END IF

      IF (L_FIRST_DESCENT) THEN

C ----------------------------------------------------------------------
CL    SECTION 1.1 TRIGONOMETRICAL VALUES.
C ----------------------------------------------------------------------

        IF(RES_DIRS.GE.6) THEN
CL    RESTRICT IN BOTH I & J DIRECTIONS.

CL    COSINE AND SECANT AT Q POINTS.

          DO J = JB,JE
            JJ = J*2-1
            DO I = 1,I_COARSE
              II = I*2-1
              COS_P_COARSE(I,J) = COS_P_FINE(II,JJ)
              SEC_P_COARSE(I,J) = 1./COS_P_COARSE(I,J)
            END DO
          END DO

          IF(VERSION.LT.3) THEN
            DO I = 1,I_COARSE
              COS_P_COARSE(I,1) = COS_P_FINE(I,1) * 2.
              SEC_P_COARSE(I,1) = 1./COS_P_COARSE(I,1)
              COS_P_COARSE(I,J_COARSE) = COS_P_FINE(I,J_FINE) * 2.
              SEC_P_COARSE(I,J_COARSE) = 1./COS_P_COARSE(I,J_COARSE)
            END DO
          END IF

CL    COSINE AT B POINTS.

          DO J = 1,J_COARSE-1
            JJ = 2*J
            DO I = 1,I_COARSE
              II = 2*I-1
              COS_V_COARSE(I,J) = COS_P_FINE(II,JJ)
            END DO
          END DO

        ELSE IF (RES_DIRS.GE.4) THEN
CL    RESTRICT IN I DIRECTION ONLY.

CL    COSINE AND SECANT AT Q POINTS.

          DO J = 1,J_COARSE
            DO I = 1,I_COARSE
              II = I*2-1
              COS_P_COARSE(I,J) = COS_P_FINE(II,J)
              SEC_P_COARSE(I,J) = 1./COS_P_COARSE(I,J)
            END DO
          END DO

CL    COSINE AT B POINTS.

          DO J = 1,J_COARSE-1
            DO I = 1,I_COARSE
              II = 2*I-1
              COS_V_COARSE(I,J) = COS_V_FINE(II,J)
           END DO
          END DO

        ELSE IF (RES_DIRS.GE.2) THEN
CL    RESTRICT IN J DIRECTION ONLY.

CL    COSINE AND SECANT AT Q POINTS.

          DO J = JB,JE
            JJ = J*2-1
            DO I = 1,I_COARSE
              COS_P_COARSE(I,J) = COS_P_FINE(I,JJ)
              SEC_P_COARSE(I,J) = 1./COS_P_COARSE(I,J)
            END DO
          END DO

          IF(VERSION.LT.3) THEN
            DO I = 1,I_COARSE
              COS_P_COARSE(I,1) = COS_P_FINE(I,1) * 2.
              SEC_P_COARSE(I,1) = 1./COS_P_COARSE(I,1)
              COS_P_COARSE(I,J_COARSE) = COS_P_FINE(I,J_FINE) * 2.
              SEC_P_COARSE(I,J_COARSE) = 1./COS_P_COARSE(I,J_COARSE)
            END DO
          END IF

CL    COSINE AT B POINTS.

          DO J = 1,J_COARSE-1
            JJ = 2*J
            DO I = 1,I_COARSE
              COS_V_COARSE(I,J) = COS_P_FINE(I,JJ)
            END DO
          END DO

        ELSE
CL    RESTRICT IN K DIRECTION ONLY, NO CHANGE TO TRIG FUNCTIONS.

CL    COSINE AND SECANT AT Q POINTS.

          DO J=1,J_COARSE
            DO I=1,I_COARSE
              COS_P_COARSE(I,J) = COS_P_FINE(I,J)
              SEC_P_COARSE(I,J) = 1./COS_P_COARSE(I,J)
            END DO
          END DO

CL    COSINE AT B POINTS.

          DO J=1,J_COARSE-1
            DO I=1,I_COARSE
              COS_V_COARSE(I,J) = COS_V_FINE(I,J)
            END DO
          END DO

        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.2 COEFFICIENTS.
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
CL    SECTION 1.2.1 COEFFICIENTS HELD AT Q LEVELS.
CL                  NAMELY A,B,C2,DEF,D,E,G.
C ----------------------------------------------------------------------

        IF(RES_DIRS.GT.5) THEN
CL    RESTRICT IN I & J DIRECTIONS, AND POSSIBLY K.
          DO K=1,K_COARSE
            IF(RES_DIRS.EQ.7) THEN
              KK = 2*K-1
            ELSE
              KK = K
            END IF

CL    A & D COEFFICIENTS.
            DO J= JB,JE
              JJ=2*J-1
              DO I= 1,I_COARSE-1
                II=2*I-1
                A_COARSE(I,J,K) = 0.5*(A_FINE(II,JJ,KK)+
     &                                 A_FINE(II+1,JJ,KK))
                D_COARSE(I,J,K) = 0.5*(D_FINE(II,JJ,KK)+
     &                                 D_FINE(II+1,JJ,KK))
              END DO
C SET END POINT ON EACH ROW.
              IF (VERSION.LT.3) THEN
C GLOBAL
                II = 2*I_COARSE
                A_COARSE(I_COARSE,J,K) =  0.5*(A_FINE(II,JJ,KK)+
     &                                         A_FINE(II-1,JJ,KK))
                D_COARSE(I_COARSE,J,K) =  0.5*(D_FINE(II,JJ,KK)+
     &                                         D_FINE(II-1,JJ,KK))
              ELSE
C LIMITED AREA COEFF NOT USED.
                A_COARSE(I_COARSE,J,K) = 0.
                D_COARSE(I_COARSE,J,K) = 0.
              END IF
            END DO

C BOUNDARIES.
            IF (VERSION.LT.3 .AND. J_NT) THEN
C GLOBAL. SET POLAR VALUES.
              DO I= 1,I_COARSE
                A_COARSE(I,1,K) = A_FINE(I,1,KK)
                A_COARSE(I,J_COARSE,K) = A_FINE(I,J_FINE,KK)
                D_COARSE(I,1,K) = D_FINE(I,1,KK)
                D_COARSE(I,J_COARSE,K) = D_FINE(I,J_FINE,KK)
              END DO
            END IF

CL    B & E COEFFICIENTS.

            DO J= 1,J_COARSE-1
              JJ=2*J-1
              DO I= 1,I_COARSE
                II=2*I-1
                B_COARSE(I,J,K) = 0.5*(B_FINE(II,JJ,KK) +
     &                                 B_FINE(II,JJ+1,KK))
                E_COARSE(I,J,K) = 0.5*(E_FINE(II,JJ,KK) +
     &                                 E_FINE(II,JJ+1,KK))
              END DO
            END DO

CL    G,C2 AND DEF  COEFFICIENTS.

            DO J= 1,J_COARSE
              JJ=2*J-1
              DO I= 1,I_COARSE
                II=2*I-1
                C2_COARSE(I,J,K) = C2_FINE(II,JJ,KK)
                DEF_COARSE(I,J,K) = DEF_FINE(II,JJ,KK)
                G_COARSE(I,J,K) = G_FINE(II,JJ,KK)
              END DO
            END DO

          END DO

        ELSE IF (RES_DIRS.GT.3) THEN
C I RESTRICTION, AND POSSIBLY K.

          DO K=1,K_COARSE
            IF(RES_DIRS.EQ.5) THEN
              KK = 2*K-1
            ELSE
              KK = K
            END IF

CL    A & D COEFFICIENTS.
            DO J= JB,JE
              DO I= 1,I_COARSE-1
                II=2*I-1
                A_COARSE(I,J,K) = 0.5*(A_FINE(II,J,KK)+
     &                                 A_FINE(II+1,J,KK))
                D_COARSE(I,J,K) = 0.5*(D_FINE(II,J,KK)+
     &                                 D_FINE(II+1,J,KK))
              END DO
C SET END POINT ON EACH ROW.
              IF (VERSION.LT.3) THEN
C GLOBAL
                II = 2*I_COARSE
                A_COARSE(I_COARSE,J,K) =  0.5*(A_FINE(II,J,KK)+
     &                                         A_FINE(II-1,J,KK))
                D_COARSE(I_COARSE,J,K) =  0.5*(D_FINE(II,J,KK)+
     &                                         D_FINE(II-1,J,KK))
              ELSE
C LIMITED AREA COEFF NOT USED.
                A_COARSE(I_COARSE,J,K) = 0.
                D_COARSE(I_COARSE,J,K) = 0.
              END IF
            END DO

C BOUNDARIES.
            IF (VERSION.LT.3 .AND. J_NT) THEN
C GLOBAL. SET POLAR VALUES.
              DO I= 1,I_COARSE
                A_COARSE(I,1,K) = A_FINE(I,1,KK)
                A_COARSE(I,J_COARSE,K) = A_FINE(I,J_FINE,KK)
                D_COARSE(I,1,K) = D_FINE(I,1,KK)
                D_COARSE(I,J_COARSE,K) = D_FINE(I,J_FINE,KK)
              END DO
            END IF

CL    B & E COEFFICIENTS.

            DO J= 1,J_COARSE-1
              DO I= 1,I_COARSE
                II=2*I-1
                B_COARSE(I,J,K) = B_FINE(II,J,KK)
                E_COARSE(I,J,K) = E_FINE(II,J,KK)
              END DO
            END DO

CL    C2,DEF & G COEFFICIENT.

            DO J= 1,J_COARSE
              DO I= 1,I_COARSE
                II=2*I-1
                C2_COARSE(I,J,K) = C2_FINE(II,J,KK)
                DEF_COARSE(I,J,K) = DEF_FINE(II,J,KK)
                G_COARSE(I,J,K) = G_FINE(II,J,KK)
              END DO
            END DO

          END DO

        ELSE IF (RES_DIRS.GE.2) THEN
CL J RESTRICTION, AND POSSIBLY K.

          DO K=1,K_COARSE
            IF(RES_DIRS.EQ.3) THEN
              KK = 2*K-1
            ELSE
              KK = K
            END IF

CL    A & D COEFFICIENTS.
            DO J= 2,J_COARSE-1
              JJ=2*J-1
              DO I= 1,I_COARSE
                A_COARSE(I,J,K) = A_FINE(I,JJ,KK)
                D_COARSE(I,J,K) = D_FINE(I,JJ,KK)
              END DO
            END DO

C BOUNDARIES.
            DO I= 1,I_COARSE
              A_COARSE(I,1,K) = A_FINE(I,1,KK)
              A_COARSE(I,J_COARSE,K) = A_FINE(I,J_FINE,KK)
              D_COARSE(I,1,K) = D_FINE(I,1,KK)
              D_COARSE(I,J_COARSE,K) = D_FINE(I,J_FINE,KK)
            END DO

CL    B & E COEFFICIENTS.

            DO J= 1,J_COARSE-1
              JJ=2*J-1
              DO I= 1,I_COARSE
                B_COARSE(I,J,K) = 0.5*(B_FINE(I,JJ,KK) +
     &                                 B_FINE(I,JJ+1,KK))
                E_COARSE(I,J,K) = 0.5*(E_FINE(I,JJ,KK) +
     &                                 E_FINE(I,JJ+1,KK))
              END DO
            END DO

CL    G,C2 AND DEF  COEFFICIENTS.

            DO J= 1,J_COARSE
              JJ=2*J-1
              DO I= 1,I_COARSE
                C2_COARSE(I,J,K) = C2_FINE(I,JJ,KK)
                DEF_COARSE(I,J,K) = DEF_FINE(I,JJ,KK)
                G_COARSE(I,J,K) = G_FINE(I,JJ,KK)
              END DO
            END DO

          END DO

        ELSE
CL K RESTRICTION ONLY
          DO K=1,K_COARSE
            KK = 2*K-1

CL    A & D COEFFICIENTS.
            DO J= 1,J_COARSE
              DO I= 1,I_COARSE
                A_COARSE(I,J,K) = A_FINE(I,J,KK)
                D_COARSE(I,J,K) = D_FINE(I,J,KK)
              END DO
            END DO

CL    B & E COEFFICIENTS.

            DO J= 1,J_COARSE-1
              DO I= 1,I_COARSE
                B_COARSE(I,J,K) = B_FINE(I,J,KK)
                E_COARSE(I,J,K) = E_FINE(I,J,KK)
              END DO
            END DO

CL    G,C2 AND DEF  COEFFICIENTS.

            DO J= 1,J_COARSE
              DO I= 1,I_COARSE
                C2_COARSE(I,J,K) = C2_FINE(I,J,KK)
                DEF_COARSE(I,J,K) = DEF_FINE(I,J,KK)
                G_COARSE(I,J,K) = G_FINE(I,J,KK)
              END DO
            END DO

          END DO

        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.2.1 COEFFICIENTS HELD BETWEEN Q LEVELS.
CL                  NAMELY C1 AND F.
C ----------------------------------------------------------------------

        IF (MOD(RES_DIRS,2).EQ.1) THEN
          DO K=1,K_COARSE-1
            KK = 2*K-1
            FRAC_ABOVE = (Z_MID_COARSE(K)-Z_MID_FINE(KK))/
     &                   (Z_MID_FINE(KK+1)-Z_MID_FINE(KK))
            FRAC_BELOW = (Z_MID_FINE(KK+1)-Z_MID_COARSE(K))/
     &                   (Z_MID_FINE(KK+1)-Z_MID_FINE(KK))
            DO J= 1,J_COARSE
              IF (J_RES) THEN
                JJ=2*J-1
              ELSE
                JJ = J
              END IF
              IF(RES_DIRS.GT.3) THEN
                DO I= 1,I_COARSE
                  II=2*I-1
                  C1_COARSE(I,J,K) = FRAC_ABOVE*C1_FINE(II,JJ,KK+1) +
     &                               FRAC_BELOW*C1_FINE(II,JJ,KK)
                  F_COARSE(I,J,K) = FRAC_ABOVE*F_FINE(II,JJ,KK+1) +
     &                              FRAC_BELOW*F_FINE(II,JJ,KK)
                END DO
              ELSE
                DO I= 1,I_COARSE
                  C1_COARSE(I,J,K) = FRAC_ABOVE*C1_FINE(I,JJ,KK+1) +
     &                               FRAC_BELOW*C1_FINE(I,JJ,KK)
                  F_COARSE(I,J,K) = FRAC_ABOVE*F_FINE(I,JJ,KK+1) +
     &                              FRAC_BELOW*F_FINE(I,JJ,KK)
                END DO
              END IF
            END DO
          END DO
        ELSE
          DO K=1,K_COARSE-1
            DO J= 1,J_COARSE
              IF (J_RES) THEN
                JJ=2*J-1
              ELSE
                JJ = J
              END IF
              IF(RES_DIRS.GT.3) THEN
                DO I= 1,I_COARSE
                  II=2*I-1
                  C1_COARSE(I,J,K) = C1_FINE(II,JJ,K)
                  F_COARSE(I,J,K) = F_FINE(II,JJ,K)
                END DO
              ELSE
                DO I= 1,I_COARSE
                  C1_COARSE(I,J,K) = C1_FINE(I,JJ,K)
                  F_COARSE(I,J,K) = F_FINE(I,JJ,K)
                END DO
              END IF
            END DO
          END DO
        END IF

      END IF

CL----------------------------------------------------------------------
CL    SECTION 2. RESTRICTION OF RESIDUAL.
CL----------------------------------------------------------------------

      IF(KREST.LE.2) THEN
C ----------------------------------------------------------------------
CL    SECTION 2.1 INJECTION AND HALF-INJECTION.
C ----------------------------------------------------------------------

C SCALE = 1 FOR FULL INJECTION, .5 FOR HALF-INJECTION.

        SCALE=1.
        IF(KREST.EQ.2) SCALE=0.5
        DO K = 1,K_COARSE
          IF( MOD(RES_DIRS,2).EQ.1) THEN
            KK=2*K-1
          ELSE
            KK = K
          END IF
          DO J = 1,J_COARSE
            IF( J_RES) THEN
              JJ=2*J-1
            ELSE
              JJ = J
            END IF
            IF( RES_DIRS.GT.3) THEN
              DO I = 1,I_COARSE
                II=2*I-1
                RHS_COARSE(I,J,K) = SCALE*RESID_FINE(II,JJ,KK)
              END DO
            ELSE
              DO I = 1,I_COARSE
                RHS_COARSE(I,J,K) = SCALE*RESID_FINE(I,JJ,KK)
              END DO
            END IF
          END DO
        END DO

      ELSE

C ----------------------------------------------------------------------
CL    SECTION 2.2 FULL-WEIGHTING.
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
CL    SECTION 2.2.1 HORIZONTAL RESTRICTION.
CL                  RESTRICTED VALUE STORED IN LEVEL_INCS READY FOR
CL                  VERTICAL RESTRICTION.
C ----------------------------------------------------------------------

        IF (RES_DIRS.GE.6) THEN
C ----------------------------------------------------------------------
CL    RESTRICT IN BOTH I & J DIRECTIONS
C ----------------------------------------------------------------------
C CALCULATE RESTRICTED VALUE ON ALL LEVELS.

          DO K=1,K_FINE
            DO J= 2,J_COARSE-1
              JJ=2*J-1
              DO I= 2,I_COARSE-1
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.125*(RESID_FINE(II+1,JJ,K) +
     &                                   RESID_FINE(II-1,JJ,K))+
     &                           (COS_P_FINE(II,JJ+1)*
     &                            (.125*RESID_FINE(II,JJ+1,K)+.0625*
     &                             (RESID_FINE(II+1,JJ+1,K)+
     &                              RESID_FINE(II-1,JJ+1,K)))+
     &                            COS_P_FINE(II,JJ-1)*
     &                            (.125*RESID_FINE(II,JJ-1,K)+.0625*
     &                             (RESID_FINE(II+1,JJ-1,K)+
     &                              RESID_FINE(II-1,JJ-1,K))))
     &                                  / COS_P_FINE(II,JJ)
              END DO

CL     SET FIRST AND LAST POINT ON EACH ROW.

              IF(VERSION.LT.3) THEN
CL GLOBAL
                LEVEL_INCS(1,J,K) =  0.25*RESID_FINE(1,JJ,K)+
     &                            0.125*(RESID_FINE(2,JJ,K) +
     &                                   RESID_FINE(I_FINE,JJ,K))+
     &                           (COS_P_FINE(1,JJ+1)*
     &                            (.125*RESID_FINE(1,JJ+1,K)+.0625*
     &                             (RESID_FINE(2,JJ+1,K)+
     &                              RESID_FINE(I_FINE,JJ+1,K)))+
     &                            COS_P_FINE(1,JJ-1)*
     &                            (.125*RESID_FINE(1,JJ-1,K)+.0625*
     &                             (RESID_FINE(2,JJ-1,K)+
     &                              RESID_FINE(I_FINE,JJ-1,K))))
     &                                  / COS_P_FINE(1,JJ)
                I = I_COARSE
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.125*(RESID_FINE(II+1,JJ,K) +
     &                                   RESID_FINE(II-1,JJ,K))+
     &                           (COS_P_FINE(II,JJ+1)*
     &                            (.125*RESID_FINE(II,JJ+1,K)+.0625*
     &                             (RESID_FINE(II+1,JJ+1,K)+
     &                              RESID_FINE(II-1,JJ+1,K)))+
     &                            COS_P_FINE(II,JJ-1)*
     &                            (.125*RESID_FINE(II,JJ-1,K)+.0625*
     &                             (RESID_FINE(II+1,JJ-1,K)+
     &                              RESID_FINE(II-1,JJ-1,K))))
     &                                  / COS_P_FINE(II,JJ)
              ELSE
CL LIMITED AREA.
                LEVEL_INCS(1,J,K) =  0.25*RESID_FINE(1,JJ,K)+
     &                            0.25*RESID_FINE(2,JJ,K) +
     &                           (COS_P_FINE(1,JJ+1)*
     &                            (.125*RESID_FINE(1,JJ+1,K)+.125*
     &                             RESID_FINE(2,JJ+1,K))+
     &                            COS_P_FINE(1,JJ-1)*
     &                            (.125*RESID_FINE(1,JJ-1,K)+.125*
     &                              RESID_FINE(2,JJ-1,K)))
     &                                  / COS_P_FINE(1,JJ)
                I = I_COARSE
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.25*RESID_FINE(II-1,JJ,K) +
     &                           (COS_P_FINE(II,JJ+1)*
     &                            (.125*RESID_FINE(II,JJ+1,K)+.125*
     &                             RESID_FINE(II-1,JJ+1,K))+
     &                            COS_P_FINE(II,JJ-1)*
     &                            (.125*RESID_FINE(II,JJ-1,K)+.125*
     &                              RESID_FINE(II-1,JJ-1,K)))
     &                                  / COS_P_FINE(II,JJ)
              END IF
            END DO

CL    SET NORTHERN AND SOUTHERN BOUNDARIES.

            IF(VERSION.LT.3) THEN
CL GLOBAL:  SET POLAR VALUES.
              SCALAR1 = 0.
              SCALAR2 = 0.
              DO I= 1,I_COARSE
                II=  I+I_COARSE
                J = J_FINE-1
                SCALAR1 = SCALAR1 + COS_P_FINE(I,2) *
     &                           (RESID_FINE(I,2,K)+RESID_FINE(II,2,K))
                SCALAR2 = SCALAR2 + COS_P_FINE(I,J) *
     &                           (RESID_FINE(I,J,K)+RESID_FINE(II,J,K))
              END DO
              SCALAR1 = .125*SCALAR1 / (I_FINE * COS_P_FINE(1,1))
              SCALAR2 = .125*SCALAR2 / (I_FINE * COS_P_FINE(1,J_FINE))
              DO I= 1,I_COARSE
                LEVEL_INCS(I,1,K) = (.25*RESID_FINE(I,1,K) + SCALAR1)
                LEVEL_INCS(I,J_COARSE,K) = (.25*RESID_FINE(I,J_FINE,K)
     &                                      + SCALAR2)
              END DO
            ELSE
CL LIMITED AREA
              J= 1
              JJ=2*J-1
              DO I= 2,I_COARSE-1
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.125*(RESID_FINE(II+1,JJ,K) +
     &                                   RESID_FINE(II-1,JJ,K))+
     &                            COS_P_FINE(II,JJ+1)*
     &                            (.25*RESID_FINE(II,JJ+1,K)+.125*
     &                             (RESID_FINE(II+1,JJ+1,K)+
     &                              RESID_FINE(II-1,JJ+1,K)))
     &                                  / COS_P_FINE(II,JJ)
              END DO
              J= J_COARSE
              JJ=2*J-1
              DO I= 2,I_COARSE-1
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.125*(RESID_FINE(II+1,JJ,K) +
     &                                   RESID_FINE(II-1,JJ,K))+
     &                            COS_P_FINE(II,JJ-1)*
     &                            (.25*RESID_FINE(II,JJ-1,K)+.125*
     &                             (RESID_FINE(II+1,JJ-1,K)+
     &                              RESID_FINE(II-1,JJ-1,K)))
     &                                  / COS_P_FINE(II,JJ)
              END DO
C SET FIRST AND LAST POINT ON EACH BOUNDARY ROW.
              J=1
              JJ=2*J-1
              LEVEL_INCS(1,J,K) =  0.25*RESID_FINE(1,JJ,K)+
     &                            0.25*RESID_FINE(2,JJ,K) +
     &                            COS_P_FINE(1,JJ+1)*
     &                            (.25*RESID_FINE(1,JJ+1,K)+.25*
     &                             RESID_FINE(2,JJ+1,K))
     &                                  / COS_P_FINE(1,JJ)
              I = I_COARSE
              II=2*I-1
              LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.25*RESID_FINE(II-1,JJ,K) +
     &                            COS_P_FINE(II,JJ+1)*
     &                            (.25*RESID_FINE(II,JJ+1,K)+.25*
     &                             RESID_FINE(II-1,JJ+1,K))
     &                                  / COS_P_FINE(II,JJ)
              J=J_COARSE
              JJ=2*J-1
              LEVEL_INCS(1,J,K) =  0.25*RESID_FINE(1,JJ,K)+
     &                            0.25*RESID_FINE(2,JJ,K) +
     &                            COS_P_FINE(1,JJ-1)*
     &                            (.25*RESID_FINE(1,JJ-1,K)+.25*
     &                             RESID_FINE(2,JJ-1,K))
     &                                  / COS_P_FINE(1,JJ)
              I = I_COARSE
              II=2*I-1
              LEVEL_INCS(I,J,K) =  0.25*RESID_FINE(II,JJ,K)+
     &                            0.25*RESID_FINE(II-1,JJ,K) +
     &                            COS_P_FINE(II,JJ-1)*
     &                            (.25*RESID_FINE(II,JJ-1,K)+.25*
     &                             RESID_FINE(II-1,JJ-1,K))
     &                                  / COS_P_FINE(II,JJ)
            END IF
          END DO

        ELSE IF (RES_DIRS.GT.3) THEN
C ----------------------------------------------------------------------
CL    RESTRICT IN I DIRECTION ONLY.
C ----------------------------------------------------------------------

          DO K=1,K_FINE
            DO J= JB,JE
              DO I= 2,I_COARSE-1
                II=2*I-1
                LEVEL_INCS(I,J,K) =  0.50*RESID_FINE(II,J,K)+
     &                              0.25*(RESID_FINE(II+1,J,K) +
     &                                    RESID_FINE(II-1,J,K))
              END DO

CL     SET FIRST AND LAST POINT ON EACH ROW.

              IF(VERSION.LT.3) THEN
CL GLOBAL
                LEVEL_INCS(1,J,K) = 0.50*RESID_FINE(1,J,K)+
     &                            0.25*(RESID_FINE(2,J,K) +
     &                                   RESID_FINE(I_FINE,J,K))
                LEVEL_INCS(I_COARSE,J,K) = 0.50*RESID_FINE(I_FINE,J,K)+
     &                            0.25*RESID_FINE(I_FINE-1,J,K)
     &                            + 0.25*RESID_FINE(1,J,K)
              ELSE
CL LIMITED AREA
                LEVEL_INCS(1,J,K) = 0.50*RESID_FINE(1,J,K)+
     &                            0.25*RESID_FINE(2,J,K)
                LEVEL_INCS(I_COARSE,J,K) = 0.50*RESID_FINE(I_FINE,J,K)+
     &                            0.25*RESID_FINE(I_FINE-1,J,K)
              END IF
            END DO

            IF( VERSION.LT.3 .AND. J_NT) THEN
C SET POLAR VALUES.
              DO I= 1,I_COARSE
                LEVEL_INCS(I,1,K) = RESID_FINE(I,1,K)
                LEVEL_INCS(I,J_COARSE,K) = RESID_FINE(I,J_FINE,K)
              END DO
            END IF
          END DO

        ELSE IF (RES_DIRS.GT.1) THEN
C ----------------------------------------------------------------------
CL    RESTRICT IN J DIRECTION ONLY.
C ----------------------------------------------------------------------

          DO K=1,K_FINE
            DO J= 2,J_COARSE-1
              JJ=2*J-1
              DO I= 1,I_COARSE
                LEVEL_INCS(I,J,K) =  0.5*RESID_FINE(I,JJ,K)+
     &                           (COS_P_FINE(I,JJ+1)*
     &                             .25*RESID_FINE(I,JJ+1,K)+
     &                            COS_P_FINE(II,JJ-1)*
     &                            .25*RESID_FINE(I,JJ-1,K))
     &                                  / COS_P_FINE(II,JJ)
              END DO
            END DO

            IF(VERSION.LT.3) THEN
CL GLOBAL
C SET POLAR VALUES.
              SCALAR1 = 0.
              SCALAR2 = 0.
              DO I= 1,I_COARSE
                J = J_FINE-1
                SCALAR1 = SCALAR1 + COS_P_FINE(I,2) *
     &                            RESID_FINE(I,2,K)
                SCALAR2 = SCALAR2 + COS_P_FINE(I,J) *
     &                            RESID_FINE(I,J,K)
              END DO
              SCALAR1 = .125*SCALAR1 / (I_COARSE * COS_P_FINE(1,1))
              SCALAR2 = .125*SCALAR2 / (I_COARSE* COS_P_FINE(1,J_FINE))
              DO I= 1,I_COARSE
                LEVEL_INCS(I,1,K) = (.25*RESID_FINE(I,1,K) + SCALAR1)
                LEVEL_INCS(I,J_COARSE,K) = (.25*RESID_FINE(I,J_FINE,K)
     &                                      + SCALAR2)
              END DO
            ELSE
CL LIMITED AREA
C SET BOUNDARY VALUES
              J= 1
              JJ=2*J-1
              DO I= 1,I_COARSE
                LEVEL_INCS(I,J,K) =  0.5*RESID_FINE(I,JJ,K)+.5*
     &                                   COS_P_FINE(I,JJ+1)*
     &                                   RESID_FINE(I,JJ+1,K)
     &                                  / COS_P_FINE(II,JJ)
              END DO
              J= J_COARSE
              JJ=2*J-1
              DO I= 1,I_COARSE
                LEVEL_INCS(I,J,K) =  0.5*RESID_FINE(I,JJ,K)+.5*
     &                                   COS_P_FINE(I,JJ-1)*
     &                                   RESID_FINE(I,JJ-1,K)
     &                                  / COS_P_FINE(II,JJ)
              END DO
            END IF
          END DO

        ELSE
C ----------------------------------------------------------------------
C K RESTRICTION ONLY, COPY RESID_FINE INTO LEVEL_INCS READY FOR VERTICAL
C RESTRICTION.
C ----------------------------------------------------------------------
          DO K=1,K_FINE
            DO J=1,J_COARSE
              DO I=1,I_COARSE
                LEVEL_INCS(I,J,K) = RESID_FINE(I,J,K)
              END DO
            END DO
          END DO
        END IF

      END IF

C ----------------------------------------------------------------------
CL    SECTION 2.2.2 VERTICAL RESTRICTION.
CL                  HORIZONTALLY RESTRICTED VALUE STORED IN LEVEL_INCS
CL                  READY FOR VERTICAL RESTRICTION.
C ----------------------------------------------------------------------

      IF(MOD(RES_DIRS,2).EQ.1) THEN
CL RESTRICT LEVEL INCS IN VERTICAL USING FULL-WEIGHTING.
C NON-BOUNDARY LEVELS.
        DO K=2,K_COARSE-1
          KK = 2*K-1
          SCALE = Z_MID_COARSE(K) - Z_MID_COARSE(K-1)
          SCALAR1 = (Z_MID_FINE(KK-1) - Z_MID_COARSE(K-1)) / SCALE
          SCALAR2 = (Z_MID_FINE(KK) - Z_MID_FINE(KK-1)) / SCALE
          SCALAR3 = (Z_MID_COARSE(K) - Z_MID_FINE(KK)) / SCALE
          DO J=1,J_COARSE
            DO I=1,I_COARSE
              RHS_COARSE(I,J,K) = SCALAR1*LEVEL_INCS(I,J,KK-1) +
     &                            SCALAR2*LEVEL_INCS(I,J,KK) +
     &                            SCALAR3*LEVEL_INCS(I,J,KK+1)
            END DO
          END DO
        END DO
C BOTTOM.
        K=1
        KK = 2*K-1
        SCALE = Z_MID_COARSE(K) - Z_MID_COARSE(K-1)
        SCALAR2 = (Z_MID_FINE(KK) - Z_MID_COARSE(K-1)) / SCALE
        SCALAR3 = (Z_MID_COARSE(K) - Z_MID_FINE(KK)) / SCALE
        DO J=1,J_COARSE
          DO I=1,I_COARSE
            RHS_COARSE(I,J,K) = SCALAR2*LEVEL_INCS(I,J,KK) +
     &                          SCALAR3*LEVEL_INCS(I,J,KK+1)
          END DO
        END DO
C TOP.
        K=K_COARSE
        KK = 2*K-1
        SCALE = Z_MID_COARSE(K) - Z_MID_COARSE(K-1)
        SCALAR1 = (Z_MID_FINE(KK-1) - Z_MID_COARSE(K-1)) / SCALE
        SCALAR2 = (Z_MID_COARSE(K) - Z_MID_FINE(KK-1)) / SCALE
        DO J=1,J_COARSE
          DO I=1,I_COARSE
            RHS_COARSE(I,J,K) = SCALAR1*LEVEL_INCS(I,J,KK-1) +
     &                          SCALAR2*LEVEL_INCS(I,J,KK)
          END DO
        END DO

      ELSE
C SIMPLY COPY RESTRICTED VALUE INTO RHS_COARSE ARRAY SINCE
C K_FINE = K_COARSE
        DO K=1,K_COARSE
          DO J=1,J_COARSE
            DO I=1,I_COARSE
              RHS_COARSE(I,J,K) = LEVEL_INCS(I,J,K)
            END DO
          END DO
        END DO

      END IF

CL----------------------------------------------------------------------
CL    SECTION 3. RESTRICTION OF SOLUTION.
CL----------------------------------------------------------------------

      DO K = 1,K_COARSE
        IF( MOD(RES_DIRS,2).EQ.1) THEN
          KK=2*K-1
        ELSE
          KK = K
        END IF
        DO J = 1,J_COARSE
          IF( J_RES) THEN
            JJ=2*J-1
          ELSE
            JJ = J
          END IF
          IF(RES_DIRS.GT.3) THEN
            DO I = 1,I_COARSE
              II=2*I-1
              Q_COARSE(I,J,K) = Q_FINE(II,JJ,KK)
            END DO
          ELSE
            DO I = 1,I_COARSE
              Q_COARSE(I,J,K) = Q_FINE(I,JJ,KK)
            END DO
          END IF
        END DO
      END DO

CL----------------------------------------------------------------------
CL    SECTION 4.  ADD COARSE GRID OPERATOR ACTING ON RESTRICTED
CL                SOLUTION TO RESTRICTED RESIDUAL TO OBTAIN RHS ON
CL                COARSE GRID.
CL----------------------------------------------------------------------

CL    CALL LEFT_HAND_SIDE TO GET COARSE GRID OPERATOR ACTING ON
CL    RESTRICTED SOLUTION.

      CALL MG_LEFT_HAND_SIDE(Q_COARSE,A_COARSE,B_COARSE,C1_COARSE,
     &                       C2_COARSE,DEF_COARSE,D_COARSE,E_COARSE,
     &                       F_COARSE,G_COARSE,
     &                       RESID_FINE,I_COARSE,J_COARSE,K_COARSE,
     &                       COS_P_COARSE,SEC_P_COARSE,COS_V_COARSE,
     &                       EARTH_RADIUS_INVERSE,LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,VERSION,I_NT,J_NT,
     &                       K_NT,Z_Q_COARSE,Z_MID_COARSE)

C LEFT_HAND_SIDE STORED IN RESID_FINE TO SAVE WORK-SPACE.

CL    ADD ON LEFT-HAND-SIDE TO RESIDUAL.
      DO K= 1,K_COARSE
        DO J= 1,J_COARSE
          DO I= 1,I_COARSE
            II = ((K-1)*J_COARSE+J-1)*I_COARSE + I
            RHS_COARSE(I,J,K) = RHS_COARSE(I,J,K) + RESID_FINE(II,1,1)
          END DO
        END DO
      END DO

CL    END OF ROUTINE MG_RESTRICT

      RETURN
      END
CLL   SUBROUTINE MG_SMOOTH
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALLS APPROPRIATE SMOOTHING ROUTINES.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_SMOOTH(
     &                     Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                     SEC_P_LATITUDE,COS_V_LATITUDE,I_DIM,J_DIM,
     &                     K_DIM,W,GRID_NUMBER,RMS_RES,IPRINT,
     &                     SWJAC,SWSYM,JAC,PAT,SYM,ILINE,JLINE,KLINE,
     &                     NFB,IPAT,LATITUDE_STEP_INVERSE,
     &                     LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                     VERSION,K_BC,I_NT,J_NT,K_NT,Z_Q,Z_MID)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  IPRINT      !IN. PARAMETER CONTROLLING QUANTITY OF OUTPUT.
     &, I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, GRID_NUMBER !IN. NUMBER OF GRID SMOOTHER IS ACTING ON.
     &, NFB         !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &              ! NEEDED.
     &, IPAT        !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &              ! POINTS.
     &, VERSION     !IN. VERSION OF MULTIGRID BEING USED.
     &, K_BC        !IN. BOUNDARY CONDITIONS IN K DIRECTION.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES
     &, ILINE   !IN. .TRUE. FOR I-LINE METHODS AND ALTERNATING SCHEMES
     &, JLINE   !IN. .TRUE. FOR J-LINE METHODS AND ALTERNATING SCHEMES
     &, KLINE   !IN. .TRUE. FOR K-LINE METHODS AND ALTERNATING SCHEMES

      REAL
     &  SWJAC    !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM    !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
     &, RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, W        !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM

      REAL
     &  A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, LONGITUDE_STEP_INVERSE !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   NO LOCAL WORK ARRAYS REQUIRED.

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      REAL
     &  SCALAR
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_I_LINE,MG_J_LINE,MG_K_LINE,MG_POLES,MG_POINTS

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. LINE RELAXATION SCHEMES
CL----------------------------------------------------------------------

CL    CALL I_LINE IF REQUIRED.

      IF(ILINE) THEN
        CALL MG_I_LINE(Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                 SEC_P_LATITUDE,COS_V_LATITUDE,
     &                 LATITUDE_STEP_INVERSE,
     &                 LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                 I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                 SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                 VERSION,K_BC,Z_Q,Z_MID,J_NT,K_NT)
      END IF

CL    CALL J_LINE IF REQUIRED.

      IF(JLINE) THEN
        CALL MG_J_LINE(Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                 SEC_P_LATITUDE,COS_V_LATITUDE,
     &                 LATITUDE_STEP_INVERSE,
     &                 LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                 I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                 SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                 VERSION,K_BC,Z_Q,Z_MID,I_NT,K_NT)
      END IF

CL    CALL K_LINE IF REQUIRED.

      IF(KLINE) THEN
        CALL MG_K_LINE(Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                 SEC_P_LATITUDE,COS_V_LATITUDE,
     &                 LATITUDE_STEP_INVERSE,
     &                 LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                 I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                 SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                 VERSION,K_BC,Z_Q,Z_MID,I_NT,J_NT)
      END IF

CL----------------------------------------------------------------------
CL    SECTION 2. POINT RELAXATION SCHEMES
CL----------------------------------------------------------------------

      IF(.NOT.ILINE .AND. .NOT.JLINE .AND. .NOT.KLINE) THEN
        CALL MG_POINTS(Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                 SEC_P_LATITUDE,COS_V_LATITUDE,I_DIM,J_DIM,K_DIM,
     &                 W,GRID_NUMBER,RMS_RES,RMS_INC,IPRINT,
     &                 SWJAC,SWSYM,JAC,PAT,SYM,
     &                 NFB,IPAT,LATITUDE_STEP_INVERSE,
     &                 LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                 VERSION,K_BC,Z_Q,Z_MID,I_NT,J_NT,K_NT)
      END IF

      IF (VERSION.LT.3 .AND. J_NT .AND. .NOT.KLINE) THEN
CL    GLOBAL CODE.
CL    CALL POLES TO PERFORM POLAR UPDATES IF KLINE NOT USED.

        CALL MG_POLES(Q,RHS,B,C1,C2,DEF,E,F,G,SEC_P_LATITUDE,
     &                COS_V_LATITUDE,
     &                LATITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                I_DIM,J_DIM,K_DIM,RMS_RES,RMS_INC,
     &                Z_Q,Z_MID,K_NT,K_BC)

      END IF

CL----------------------------------------------------------------------
CL    SECTION 2. CALCULATE NORMS.
CL----------------------------------------------------------------------

      SCALAR = 1./(I_DIM*J_DIM*K_DIM)
      RMS_RES = SQRT(SCALAR*RMS_RES)
      RMS_INC = SQRT(SCALAR*RMS_INC)

      IF(IPRINT.GE.2) WRITE(6,9000) GRID_NUMBER,RMS_RES,RMS_INC

 9000 FORMAT(1X,'GRID',I3,3X,'RESIDUAL NORM=',E9.3,6X,
     &                       'INCREMENT NORM=',E9.3)

CL    END OF ROUTINE MG_SMOOTH

      RETURN
      END
CLL   SUBROUTINE MG_STYPE
CLL
CLL   PURPOSE:
CLL   -------
CLL   SETS SMOOTHING INFORMATION GIVEN A PARTICULAR CHOICE.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_STYPE(I_NT,J_NT,K_NT,KSMOOTH,JAC,PAT,SYM,ILINE,
     &                    JLINE,KLINE,NFB,IPAT,SWJAC,SWSYM)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  KSMOOTH !IN. NUMBER REPRESENTING SMOOTHER SELECTED.
     &          !    1 - I-LINE JACOBI
     &          !    2 - J-LINE JACOBI
     &          !    3 - K-LINE JACOBI
     &          !    4 - I&J-LINE JACOBI
     &          !    5 - I&K LINE JACOBI
     &          !    6 - J&K-LINE JACOBI
     &          !    7 - 3-D ALTERNATING LINE JACOBI
     &          !    8 - I-LINE GAUSS-SEIDEL
     &          !    9 - J-LINE GAUSS-SEIDEL
     &          !   10 - K-LINE GAUSS-SEIDEL
     &          !   11 - I&J-LINE GAUSS-SEIDEL
     &          !   12 - I&K-LINE GAUSS-SEIDEL
     &          !   13 - J&K-LINE GAUSS-SEIDEL
     &          !   14 - 3-D ALTERNATING LINE GAUSS-SEIDEL
     &          !   15 - 3-D ALTERNATING SYMMETRIC LINE GAUSS-SEIDEL
     &          !   16 - I-LINE ZEBRA
     &          !   17 - J-LINE ZEBRA
     &          !   18 - K-LINE ZEBRA
     &          !   19 - I&J-LINE ZEBRA
     &          !   20 - I&K-LINE ZEBRA
     &          !   21 - J&K-LINE ZEBRA
     &          !   22 - 3-D ALTERNATING LINE ZEBRA
     &          !   23 - JACOBI POINT SMOOTHER
     &          !   24 - GAUSS-SEIDEL POINT SMOOTHER
     &          !   25 - SYMMETRIC GAUSS-SEIDEL POINT SMOOTHER
     &          !   26 - RED-BLACK POINT SMOOTHER

      LOGICAL
     &  JAC     !OUT. .TRUE. FOR JACOBI METHODS
     &, PAT     !OUT. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !OUT. .TRUE. FOR SYMMETRIC SCHEMES
     &, ILINE   !OUT. .TRUE. FOR I-LINE METHODS AND ALTERNATING SCHEMES
     &, JLINE   !OUT. .TRUE. FOR J-LINE METHODS AND ALTERNATING SCHEMES
     &, KLINE   !OUT. .TRUE. FOR K-LINE METHODS AND ALTERNATING SCHEMES

      INTEGER
     &  NFB     !OUT.NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID NEEDED
     &, IPAT    !OUT.AN INTEGER USED TO DISTINGUISH RED AND BLACK POINTS

      REAL
     &  SWJAC   !OUT. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM   !OUT. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS

C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C NO LOCAL VARIABLES.

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. SET SWITCHES DEPENDING ON SMOOTHER REQUESTED.
CL----------------------------------------------------------------------

CL    SET JACOBI SWITCH
C     -----------------
      JAC=.FALSE.
      IF(KSMOOTH.LE.7.OR.KSMOOTH.EQ.23) JAC=.TRUE.

CL    SET PATTERN SWITCH
C     ------------------
      PAT=.FALSE.
      IF(KSMOOTH.GE.16.AND.KSMOOTH.LT.23.OR.KSMOOTH.EQ.26) PAT=.TRUE.

CL    SET SYMMETRIC SWITCH
C     --------------------
      SYM=.FALSE.
      IF(KSMOOTH.EQ.15.OR.KSMOOTH.EQ.25) SYM=.TRUE.

CL    LINE METHODS
C     ------------
      ILINE=.FALSE.
      JLINE=.FALSE.
      KLINE=.FALSE.

CL    ALTERNATING
C     ===========
      IF(KSMOOTH.EQ.7.OR.KSMOOTH.EQ.14.OR.
     &   KSMOOTH.EQ.15.OR.KSMOOTH.EQ.22) THEN
        ILINE=.TRUE.
        JLINE=.TRUE.
        KLINE=.TRUE.
      ENDIF

CL    I-LINE
C     ======
      IF(KSMOOTH.EQ.1.OR.KSMOOTH.EQ.4.OR.KSMOOTH.EQ.5.OR.
     &   KSMOOTH.EQ.8.OR.KSMOOTH.EQ.11.OR.KSMOOTH.EQ.12.OR.
     &   KSMOOTH.EQ.16.OR.KSMOOTH.EQ.19.OR.KSMOOTH.EQ.20) ILINE=.TRUE.

CL    J-LINE
C     ======
      IF(KSMOOTH.EQ.2.OR.KSMOOTH.EQ.4.OR.KSMOOTH.EQ.6.OR.
     &   KSMOOTH.EQ.9.OR.KSMOOTH.EQ.11.OR.KSMOOTH.EQ.13.OR.
     &   KSMOOTH.EQ.17.OR.KSMOOTH.EQ.19.OR.KSMOOTH.EQ.21) JLINE=.TRUE.

CL    K-LINE
C     ======
      IF(KSMOOTH.EQ.3.OR.KSMOOTH.EQ.5.OR.KSMOOTH.EQ.6.OR.
     &   KSMOOTH.EQ.10.OR.KSMOOTH.EQ.12.OR.KSMOOTH.EQ.13.OR.
     &   KSMOOTH.EQ.18.OR.KSMOOTH.EQ.20.OR.KSMOOTH.EQ.21) KLINE=.TRUE.

      IF (ILINE .AND. .NOT. I_NT) THEN
        WRITE(6,*)' *********** ERROR ***********'
        WRITE(6,*)' I-LINE SMOOTHING REQUESTED WHEN ONLY 1 POINT IN I',
     &         '-DIRECTION.'
        WRITE(6,*)' I-LINE SMOOTHING NOT ENABLED.'
        WRITE(6,*)' *********** ERROR ***********'
        ILINE=.FALSE.
      END IF

      IF (JLINE .AND. .NOT. J_NT) THEN
        WRITE(6,*)' *********** ERROR ***********'
        WRITE(6,*)' J-LINE SMOOTHING REQUESTED WHEN ONLY 1 POINT IN J',
     &         '-DIRECTION.'
        WRITE(6,*)' J-LINE SMOOTHING NOT ENABLED.'
        WRITE(6,*)' *********** ERROR ***********'
        JLINE=.FALSE.
      END IF

      IF (KLINE .AND. .NOT. K_NT) THEN
        WRITE(6,*)' *********** ERROR ***********'
        WRITE(6,*)' K-LINE SMOOTHING REQUESTED WHEN ONLY 1 POINT IN K',
     &         '-DIRECTION.'
        WRITE(6,*)' K-LINE SMOOTHING NOT ENABLED.'
        WRITE(6,*)' *********** ERROR ***********'
        KLINE=.FALSE.
      END IF

      IF (KSMOOTH.LT.23 .AND. .NOT.ILINE .AND. .NOT.JLINE .AND.
     &                        .NOT.KLINE) THEN
        WRITE(6,*)' *********** ERROR ***********'
        WRITE(6,*)' LINE SMOOTHING REQUESTED BUT NOT POSSIBLE. '
        WRITE(6,*)' CHECK ERROR MESSAGES AND USER SUPPLIED INPUTS.'
      WRITE(6,*)' RED-BLACK POINT SMOOTHING BEING ENABLED FOR THIS RUN.'
        WRITE(6,*)' *********** ERROR ***********'
        KSMOOTH = 26
      END IF

CL    SWITCHES
C     --------
      NFB=1
      IPAT=1
      SWJAC=1.0
      SWSYM=1.0

      IF(SYM.OR.PAT) NFB=2
      IF(PAT) IPAT=2
      IF(JAC) SWJAC=0.0
      IF(SYM) SWSYM=0.0

CL    END OF ROUTINE MG_STYPE

      RETURN
      END
CLL   SUBROUTINE MG_FDMG
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS MAXITS ITERATIONS OF THE FULL APPROXIMATION MULTI-GRID
CLL   SCHEME (FAS-MG).
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_FDMG(
     &                   Q,A,B,C1,C2,DEF,D,E,F,G,RHS,I_NT,J_NT,K_NT,
     &                   START_ADDRESS,LAST_ADDRESS,IMAX,
     &                   JMAX,KMAX,JAC,PAT,SYM,ILINE,JLINE,KLINE,NFB,
     &                   IPAT,SWJAC,SWSYM,I_LENGTH,J_LENGTH,K_LENGTH,
     &                   NGRIDS,MAXITS,TOL_RES,IPRINT,
     &                   KSMOOTH,NPRE,NPOST,NCOARSE,W,KREST,NCGC,
     &                   COS_P_LATITUDE,SEC_P_LATITUDE,
     &                   COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                   LONGITUDE_STEP_GRIDS,LATITUDE_STEP_GRIDS,
     &                   START_ADDRESS_2D,LAST_ADDRESS_2D,
     &                   RES_DIRS,WORST_SMOOTHING_RATE,VERSION,K_BC,
     &                   START_ADDRESS_Z,LAST_ADDRESS_Z,Z_Q,Z_MID)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     & VERSION     !IN. DOMAIN ON WHICH TO SOLVE PROBLEM.
     &,K_BC        !IN. VERTICAL BOUNDARY CONDITIONS.

      INTEGER
     &  NGRIDS     !IN. NUMBER OF GRIDS.
     &, MAXITS     !IN. MAX NO OF FAS ITERATIONS WITHOUT CONVERGENCE
     &, IPRINT     !IN. PARAMETER CONTROLLING QUANTITY OF OUTPUT.
     &, KSMOOTH    !IN. KIND OF ITERATIVE METHOD USED AS A SMOOTHER.

      INTEGER
     &  NPRE       !IN. NO OF PRE-SMOOTHING SWEEPS
     &, NPOST      !IN. NO OF POST-SMOOTHING SWEEPS
     &, NCOARSE    !IN. NO OF ITERATIONS OF SMOOTHER ON COARSEST MESH
     &, KREST      !IN. KIND OF RESTRICTION USED.
     &, NCGC       !IN. NO OF COARSE GRID CORRECTIONS.
     &, I_LENGTH   !IN. NUMBER OF POINTS IN I DIRECTION.
     &, J_LENGTH   !IN. NUMBER OF POINTS IN J DIRECTION.
     &, K_LENGTH   !IN. NUMBER OF POINTS IN K DIRECTION.

      REAL
     &  A(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, B(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, C1(I_LENGTH,J_LENGTH,K_LENGTH)  !IN. COEFFICIENT
     &, C2(I_LENGTH,J_LENGTH,K_LENGTH)  !IN. COEFFICIENT
     &, DEF(I_LENGTH,J_LENGTH,K_LENGTH) !IN. COEFFICIENT
     &, D(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, E(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, F(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, G(I_LENGTH,J_LENGTH,K_LENGTH)   !IN. COEFFICIENT
     &, RHS(I_LENGTH,J_LENGTH,K_LENGTH)
     &                          !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_LENGTH,J_LENGTH)
     &                          !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_LENGTH,J_LENGTH)
     &                          !IN. 1./COSINE OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_LENGTH,J_LENGTH)
     &                          !IN. COSINE OF LATITUDE AT B POINTS.
     &, Z_Q(K_LENGTH)           !IN. VALUE OF Z AT Q POINTS
     &, Z_MID(K_LENGTH+1)       !IN. VALUE OF Z AT MID POINTS DEFINED
     &                          !    BETWEEN Q POINTS IN VERTICAL
     &                          !    FIRST VALUE IS BELOW FIRST Q POINT

      REAL
     &  Q(I_LENGTH,J_LENGTH,K_LENGTH)    !INOUT. SOLUTION.

      INTEGER
     &  IMAX(NGRIDS)  !IN. NUMBER OF NODES IN THE I-DIRECTION
     &                ! ON EACH GRID.
     &, JMAX(NGRIDS)  !IN. NUMBER OF NODES IN THE J-DIRECTION
     &                ! ON EACH GRID.
     &, KMAX(NGRIDS)  !IN. NUMBER OF NODES IN THE K-DIRECTION
     &                ! ON EACH GRID.
     &, START_ADDRESS(NGRIDS) !IN. START ADDRESS IN DATA ARRAY FOR
     &                        !    EACH GRID.
     &, START_ADDRESS_2D(NGRIDS) !IN. START ADDRESS IN DATA ARRAY FOR
     &                        !       EACH 2-D GRID.
     &, START_ADDRESS_Z(NGRIDS) !IN. START ADDRESS IN DATA ARRAY FOR
     &                        !       EACH Z GRID.
     &, RES_DIRS(NGRIDS)      ! RESTRICTED DIRECTIONS.

      REAL
     &  TOL_RES                     !IN. TOLERANCE FOR RESIDUAL NORM
     &                              !    RELATIVE TO INITIAL RESIDUAL.
     &, WORST_SMOOTHING_RATE        !IN. WORST PRACTICAL SMOOTHING RATE
     &                              !    ACCEPTABLE BEFORE CONVERGENCE
     &                              !    OR MAXIMUM ITERATIONS REACHED.
     &, W                           !IN. RELAXATION PARAMETER FOR EACH
     &                              !    VARIABLE IN SYSTEM
     &, LATITUDE_STEP_GRIDS(NGRIDS) !IN.
     &, LONGITUDE_STEP_GRIDS(NGRIDS)!IN.
     &, EARTH_RADIUS_INVERSE        !IN.

      INTEGER
     &  LAST_ADDRESS  !IN. LAST ADDRESS OF DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION SPACE REQUIRED.
     &, LAST_ADDRESS_2D !IN. LAST ADDRESS OF 2-D DATA ARRAY NEEDED. USED
     &                ! TO DIMENSION SPACE REQUIRED.
     &, LAST_ADDRESS_Z !IN. LAST ADDRESS OF Z DATA ARRAY NEEDED. USED
     &                ! TO DIMENSION SPACE REQUIRED.
     &, NFB           !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &                ! NEEDED.
     &, IPAT          !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &                ! POINTS.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES
     &, ILINE   !IN. .TRUE. FOR I-LINE METHODS AND ALTERNATING SCHEMES
     &, JLINE   !IN. .TRUE. FOR J-LINE METHODS AND ALTERNATING SCHEMES
     &, KLINE   !IN. .TRUE. FOR K-LINE METHODS AND ALTERNATING SCHEMES

      REAL
     &  SWJAC   !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM   !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
C*----------------------------------------------------------------------

C*L  18 LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  Q_GRIDS(LAST_ADDRESS)   ! SPACE FOR Q ON ALL GRIDS
     &, RHS_GRIDS(LAST_ADDRESS) ! SPACE FOR RHS ON ALL GRIDS
     &, A_GRIDS(LAST_ADDRESS)   ! SPACE FOR A COEFFICIENT ON ALL GRIDS
     &, B_GRIDS(LAST_ADDRESS)   ! SPACE FOR B COEFFICIENT ON ALL GRIDS
     &, C1_GRIDS(LAST_ADDRESS)  ! SPACE FOR C1 COEFFICIENT ON ALL GRIDS
     &, C2_GRIDS(LAST_ADDRESS)  ! SPACE FOR C2 COEFFICIENT ON ALL GRIDS
     &, D_GRIDS(LAST_ADDRESS)   ! SPACE FOR D COEFFICIENT ON ALL GRIDS
     &, DEF_GRIDS(LAST_ADDRESS) ! SPACE FOR DEF COEFFICIENT ON ALL GRIDS
     &, E_GRIDS(LAST_ADDRESS)   ! SPACE FOR E COEFFICIENT ON ALL GRIDS
     &, F_GRIDS(LAST_ADDRESS)   ! SPACE FOR F COEFFICIENT ON ALL GRIDS
     &, G_GRIDS(LAST_ADDRESS)   ! SPACE FOR G COEFFICIENT ON ALL GRIDS
     &, Z_Q_GRIDS(LAST_ADDRESS_Z)   ! SPACE FOR Z_Q ON ALL GRIDS
     &, Z_MID_GRIDS(LAST_ADDRESS_Z) ! SPACE FOR Z_MID ON ALL GRIDS
     &, WORK_SPACE(I_LENGTH*J_LENGTH*K_LENGTH)
     &, COS_P_GRIDS(LAST_ADDRESS_2D)
     &, SEC_P_GRIDS(LAST_ADDRESS_2D)
     &, COS_V_GRIDS(LAST_ADDRESS_2D)

      INTEGER
     &  IDGC(NGRIDS)            ! COUNTER FOR NUMBER OF TIMES GRID
     &                          ! CORRECTIONED IN DESCENDING PART OF MG.

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,J,K,I2,BASE         ! LOOP COUNTERS
     &, MG_ITERATIONS         ! MULTI-GRID ITERATION NUMBER.
     &, GRID_CURRENT          ! CURRENT GRID NUMBER.
     &, I_CURRENT             ! NUMBER OF I POINTS ON CURRENT GRID
     &, J_CURRENT             ! NUMBER OF J POINTS ON CURRENT GRID
     &, K_CURRENT             ! NUMBER OF K POINTS ON CURRENT GRID
     &, ADDRESS_CURRENT       ! START ADDRESS FOR CURRENT GRID.
     &, ADDRESS_CURRENT_2D    ! START ADDRESS FOR CURRENT 2-D GRID.
     &, ADDRESS_CURRENT_Z     ! START ADDRESS FOR CURRENT Z GRID.
     &, I_FINER               ! NUMBER OF I POINTS ON FINER GRID
     &, J_FINER               ! NUMBER OF J POINTS ON FINER GRID
     &, K_FINER               ! NUMBER OF K POINTS ON FINER GRID
     &, I_COARSER             ! NUMBER OF I POINTS ON COARSER GRID
     &, J_COARSER             ! NUMBER OF J POINTS ON COARSER GRID
     &, K_COARSER             ! NUMBER OF K POINTS ON COARSER GRID
     &, START_ADDRESS_FINER   ! START ADDRESS OF FINER GRID.
     &, START_ADDRESS_FINER_2D! START ADDRESS OF 2-D FINER GRID.
     &, START_ADDRESS_FINER_Z ! START ADDRESS OF Z FINER GRID.
     &, START_ADDRESS_COARSER ! START ADDRESS OF COARSER GRID.
     &, START_GRID            ! START GRID FOR DESCENDING PART OF MG
     &                        ! CYCLE.

      REAL
     &  RMS_RES              ! ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_RES_INIT         ! INITIAL ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_RES_OLD          ! ROOT MEAN SQUARE RESIDUAL NORM FROM
     &                       ! PREVIOUS ITERATION.
     &, CONVERGENCE_RATE     ! MG CONVERGENCE RATE.
     &, PRACTICAL_SMOOTHING_RATE ! SEE CODE FOR DETAILS.

      LOGICAL
     &  L_CNTL               ! CONTROLS DO WHILE LOOPS.
     &, L_DESCEND            ! TRUE IF MG-CYCLE HAS A DESCENT PART TO
     &                       ! DO.
     &, L_ASCEND             ! TRUE IF MG-CYCLE IS IN AN ASCENT STAGE.
     &, L_ITERATE            ! TRUE FOR MULTI-GRID ITERATIONS.
     &, L_FIRST_DESCENT      ! TRUE IF FIRST TIME THROUGH DESCENT LOOP.

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_SMOOTH,MG_RESTRICT,MG_PROLONG,
     &         MG_LEFT_HAND_SIDE,MG_CALC_Z

C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. INITIALISATION.
CL----------------------------------------------------------------------

CL    COPY A,B,C,G,RHS AND Q INTO _GRIDS LOCATIONS.

      DO K= 1,K_LENGTH
        IF(K_NT)Z_Q_GRIDS(START_ADDRESS_Z(NGRIDS)+K-1) = Z_Q(K)
        BASE = START_ADDRESS(NGRIDS) + (K-1)*I_LENGTH*J_LENGTH - 1
        DO J= 1,J_LENGTH
          I2 = (J-1)*I_LENGTH
          DO I= 1,I_LENGTH
            Q_GRIDS(BASE+I2+I)   = Q(I,J,K)
            A_GRIDS(BASE+I2+I)   = A(I,J,K)
            B_GRIDS(BASE+I2+I)   = B(I,J,K)
            C2_GRIDS(BASE+I2+I)  = C2(I,J,K)
            DEF_GRIDS(BASE+I2+I) = DEF(I,J,K)
            D_GRIDS(BASE+I2+I)   = D(I,J,K)
            E_GRIDS(BASE+I2+I)   = E(I,J,K)
            G_GRIDS(BASE+I2+I)   = G(I,J,K)
            RHS_GRIDS(BASE+I2+I) = RHS(I,J,K)
          END DO
        END DO
      END DO
      IF(K_NT) THEN
        DO K= 1,K_LENGTH-1
          BASE = START_ADDRESS(NGRIDS) + (K-1)*I_LENGTH*J_LENGTH - 1
          DO J= 1,J_LENGTH
            I2 = (J-1)*I_LENGTH
            DO I= 1,I_LENGTH
              C1_GRIDS(BASE+I2+I)   = C1(I,J,K)
              F_GRIDS(BASE+I2+I)    = F(I,J,K)
            END DO
          END DO
        END DO
        DO K= 1,K_LENGTH+1
          Z_MID_GRIDS(START_ADDRESS_Z(NGRIDS)+K-1) = Z_MID(K)
        END DO
      END IF
      BASE = START_ADDRESS_2D(NGRIDS) - 1
      DO J= 1,J_LENGTH
        I2 = (J-1)*I_LENGTH
        DO I= 1,I_LENGTH
          COS_P_GRIDS(BASE+I2+I) = COS_P_LATITUDE(I,J)
          SEC_P_GRIDS(BASE+I2+I) = SEC_P_LATITUDE(I,J)
        END DO
      END DO
      DO J= 1,J_LENGTH-1
        I2 = (J-1)*I_LENGTH
        DO I= 1,I_LENGTH
          COS_V_GRIDS(BASE+I2+I) = COS_V_LATITUDE(I,J)
        END DO
      END DO

      RMS_RES_OLD = 1.0

CL----------------------------------------------------------------------
CL    SECTION 2. IF ONLY 1 GRID.
CL----------------------------------------------------------------------

      IF(NGRIDS.EQ.1) THEN
        IF(IPRINT.GE.2) WRITE(6,9000)

        GRID_CURRENT=NGRIDS
        I_CURRENT=IMAX(GRID_CURRENT)
        J_CURRENT=JMAX(GRID_CURRENT)
        K_CURRENT=KMAX(GRID_CURRENT)
        ADDRESS_CURRENT=START_ADDRESS(GRID_CURRENT)
        ADDRESS_CURRENT_2D=START_ADDRESS_2D(GRID_CURRENT)
        ADDRESS_CURRENT_Z=START_ADDRESS_Z(GRID_CURRENT)

CL    CALL SMOOTHER.

        L_CNTL = .TRUE.
        I=1
        DO WHILE (L_CNTL)
          CALL MG_SMOOTH(Q_GRIDS(ADDRESS_CURRENT),
     &                   RHS_GRIDS(ADDRESS_CURRENT),
     &                   A_GRIDS(ADDRESS_CURRENT),
     &                   B_GRIDS(ADDRESS_CURRENT),
     &                   C1_GRIDS(ADDRESS_CURRENT),
     &                   C2_GRIDS(ADDRESS_CURRENT),
     &                   DEF_GRIDS(ADDRESS_CURRENT),
     &                   D_GRIDS(ADDRESS_CURRENT),
     &                   E_GRIDS(ADDRESS_CURRENT),
     &                   F_GRIDS(ADDRESS_CURRENT),
     &                   G_GRIDS(ADDRESS_CURRENT),
     &                   COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                   SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                   COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                   I_CURRENT,J_CURRENT,K_CURRENT,
     &                   W,GRID_CURRENT,RMS_RES,IPRINT,
     &                   SWJAC,SWSYM,JAC,PAT,SYM,ILINE,JLINE,KLINE,
     &                   NFB,IPAT,LATITUDE_STEP_GRIDS(GRID_CURRENT),
     &                   LONGITUDE_STEP_GRIDS(GRID_CURRENT),
     &                   EARTH_RADIUS_INVERSE,VERSION,K_BC,I_NT,J_NT,
     &                   K_NT,Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                   Z_MID_GRIDS(ADDRESS_CURRENT_Z))

          IF(I.EQ.1) THEN
            RMS_RES_INIT = RMS_RES
          END IF

CL    CHECK CONVERGENCE.

          IF(RMS_RES.LT.TOL_RES*RMS_RES_INIT) THEN
            IF(IPRINT.GE.1) WRITE(6,9001) I
            L_CNTL = .FALSE.
          END IF

CL    UPDATE LOOP COUNTER. IF EXCEEDED NUMBER OF ITERATIONS STOP
CL    OUTPUT FAILURE MESSAGE IF ASKED FOR.
          I=I+1
          IF(I.GT.MAXITS) THEN
            L_CNTL = .FALSE.
            IF(IPRINT.GE.1) WRITE(6,9002) MAXITS
          END IF
C END DO LOOP
        END DO

      ELSE

CL----------------------------------------------------------------------
CL    SECTION 3. PERFORM MAXITS OF FAS-MG SCHEME.
CL----------------------------------------------------------------------

C ----------------------------------------------------------------------
CL    SECTION 3.0. SET UP VERTICAL LEVEL INFORMATION FOR ALL GRIDS,
CL                 IF A 3-D RUN.
C ----------------------------------------------------------------------

        IF(K_NT) THEN
          CALL MG_CALC_Z(Z_Q_GRIDS,Z_MID_GRIDS,START_ADDRESS_Z,
     &                   LAST_ADDRESS_Z,NGRIDS,KMAX,RES_DIRS)
          IF(IPRINT.GE.3) THEN
            DO I=1,NGRIDS
              WRITE(6,*)' GRID ',I,' LEVELS= ',KMAX(I)
              DO K=1,KMAX(I)
      WRITE(6,*)' Z_MID = ',Z_MID_GRIDS(START_ADDRESS_Z(I)+K-1),
     &                 ' Z_Q = ',Z_Q_GRIDS(START_ADDRESS_Z(I)+K-1)
              END DO
              K= KMAX(I) + 1
              WRITE(6,*)' Z_MID = ',Z_MID_GRIDS(START_ADDRESS_Z(I)+K-1)
            END DO
          END IF
        END IF

C ----------------------------------------------------------------------
CL    SECTION 3.01.CONTROL LOOP BEGINS.
C ----------------------------------------------------------------------

        IF(IPRINT.GE.2) WRITE(6,9000)

        L_ITERATE = .TRUE.
        L_FIRST_DESCENT = .TRUE.
        MG_ITERATIONS = 0

CL CONTROL LOOP OVER ITERATIONS.
        DO WHILE (L_ITERATE)

          MG_ITERATIONS = MG_ITERATIONS + 1
          IF(IPRINT.GE.2) WRITE(6,9005) MG_ITERATIONS

C SET COUNTER FOR NUMBER OF COARSE GRID CORRECTIONS

          DO I=1,NGRIDS
            IDGC(I)=0
          END DO

CL    SET ADDRESSING AND NUMBER OF POINTS.

          I_CURRENT=IMAX(NGRIDS)
          J_CURRENT=JMAX(NGRIDS)
          K_CURRENT=KMAX(NGRIDS)
          ADDRESS_CURRENT=START_ADDRESS(NGRIDS)
          ADDRESS_CURRENT_2D=START_ADDRESS_2D(NGRIDS)
          ADDRESS_CURRENT_Z=START_ADDRESS_Z(NGRIDS)
          START_GRID = NGRIDS
          L_DESCEND = .TRUE.

C ----------------------------------------------------------------------
CL    SECTION 3.1. PERFORM DESCENT PART OF MULTI-GRID.
CL                 IF A V-CYCLE THEN CODE SIMPLY COMES STRAIGHT BACK UP.
CL                 IF A W-CYCLE THEN CODE KEEPS COMING UP ONE FURTHUR
CL                 GRID EACH TIME THROUGH.
CL    EG: 4 GRIDS.
CL    V-CYCLE GOES 4           4  W-CYCLE 4               4
CL                  \         /            \             /
CL                   3       3              3           3
CL                    \     /                \         /
CL                     2   2                  2   2   2
CL                      \ /                    \ / \ /
CL                       1                      1   1
CL
CL    FULL MULTI-GRID GOES 4                        4
CL                          \                      /
CL                           3            3       3
CL                            \          / \     /
CL                             2   2   2    2   2
CL                              \ / \ /      \ /
CL                               1   1        1
CL
CL                  THUS ON A W-CYCLE THE DESCENT PART OF THE CODE CAN
CL                  BE CALLED MORE THAN ONCE STARTING ON DIFFERENT GRIDS
C ----------------------------------------------------------------------

          DO WHILE (L_DESCEND)
CL    LOOP OVER ALL GRIDS EXCEPT COARSEST.

            DO GRID_CURRENT = START_GRID,2,-1

C ----------------------------------------------------------------------
CL    SECTION 3.1.1 CALL SMOOTH TO CALCULATE SOLUTION ON CURRENT GRID.
C ----------------------------------------------------------------------

CL    SMOOTHING

              DO I=1,NPRE
                CALL MG_SMOOTH(Q_GRIDS(ADDRESS_CURRENT),
     &                         RHS_GRIDS(ADDRESS_CURRENT),
     &                         A_GRIDS(ADDRESS_CURRENT),
     &                         B_GRIDS(ADDRESS_CURRENT),
     &                         C1_GRIDS(ADDRESS_CURRENT),
     &                         C2_GRIDS(ADDRESS_CURRENT),
     &                         DEF_GRIDS(ADDRESS_CURRENT),
     &                         D_GRIDS(ADDRESS_CURRENT),
     &                         E_GRIDS(ADDRESS_CURRENT),
     &                         F_GRIDS(ADDRESS_CURRENT),
     &                         G_GRIDS(ADDRESS_CURRENT),
     &                         COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                         I_CURRENT,J_CURRENT,K_CURRENT,
     &                         W,GRID_CURRENT,RMS_RES,IPRINT,
     &                         SWJAC,SWSYM,JAC,PAT,SYM,ILINE,JLINE,
     &                         KLINE,NFB,IPAT,
     &                         LATITUDE_STEP_GRIDS(GRID_CURRENT),
     &                         LONGITUDE_STEP_GRIDS(GRID_CURRENT),
     &                         EARTH_RADIUS_INVERSE,VERSION,K_BC,
     &                         I_NT,J_NT,K_NT,
     &                         Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                         Z_MID_GRIDS(ADDRESS_CURRENT_Z))

                IF(I.EQ.1 .AND. MG_ITERATIONS.EQ.1) THEN
                  RMS_RES_INIT = RMS_RES
                END IF

              END DO

C UPDATE DESCENDING GRID COUNTER.

              IDGC(GRID_CURRENT) = IDGC(GRID_CURRENT)+1

C ----------------------------------------------------------------------
CL    SECTION 3.1.2 CALCULATE RESIDUAL BY CALLING LEFT_HAND_SIDE AND
CL                  SUBTRACTING FROM RIGHT-HAND-SIDE.
C ----------------------------------------------------------------------

CL    CALL LEFT_HAND_SIDE.

C LHS VALUE IS RETURNED IN WORK_SPACE

              CALL MG_LEFT_HAND_SIDE(Q_GRIDS(ADDRESS_CURRENT),
     &                               A_GRIDS(ADDRESS_CURRENT),
     &                               B_GRIDS(ADDRESS_CURRENT),
     &                               C1_GRIDS(ADDRESS_CURRENT),
     &                               C2_GRIDS(ADDRESS_CURRENT),
     &                               DEF_GRIDS(ADDRESS_CURRENT),
     &                               D_GRIDS(ADDRESS_CURRENT),
     &                               E_GRIDS(ADDRESS_CURRENT),
     &                               F_GRIDS(ADDRESS_CURRENT),
     &                               G_GRIDS(ADDRESS_CURRENT),
     &                               WORK_SPACE,
     &                               I_CURRENT,J_CURRENT,K_CURRENT,
     &                               COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                               SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                               COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                               EARTH_RADIUS_INVERSE,
     &                               LATITUDE_STEP_GRIDS(GRID_CURRENT),
     &                               LONGITUDE_STEP_GRIDS(GRID_CURRENT),
     &                               VERSION,I_NT,J_NT,K_NT,
     &                               Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                               Z_MID_GRIDS(ADDRESS_CURRENT_Z))

CL    CALCULATE RESIDUAL.
C RESIDUAL IS STORED IN WORK_SPACE.

              DO I=1,I_CURRENT*J_CURRENT*K_CURRENT
                WORK_SPACE(I) = RHS_GRIDS(ADDRESS_CURRENT+I-1)
     &                           - WORK_SPACE(I)
              END DO

C ----------------------------------------------------------------------
CL    SECTION 3.1.3 RESTRICT FROM CURRENT GRID TO COARSER ONE.
C ----------------------------------------------------------------------

              I_FINER = I_CURRENT
              J_FINER = J_CURRENT
              K_FINER = K_CURRENT
              START_ADDRESS_FINER = ADDRESS_CURRENT
              START_ADDRESS_FINER_2D = ADDRESS_CURRENT_2D
              START_ADDRESS_FINER_Z = ADDRESS_CURRENT_Z
              I_CURRENT = IMAX(GRID_CURRENT-1)
              J_CURRENT = JMAX(GRID_CURRENT-1)
              K_CURRENT = KMAX(GRID_CURRENT-1)
              ADDRESS_CURRENT = START_ADDRESS(GRID_CURRENT-1)
              ADDRESS_CURRENT_2D = START_ADDRESS_2D(GRID_CURRENT-1)
              ADDRESS_CURRENT_Z = START_ADDRESS_Z(GRID_CURRENT-1)

CL    CALL RESTRICT.

              CALL MG_RESTRICT(Q_GRIDS(ADDRESS_CURRENT),
     &                         Q_GRIDS(START_ADDRESS_FINER),
     &                         RHS_GRIDS(ADDRESS_CURRENT),WORK_SPACE,
     &                         A_GRIDS(ADDRESS_CURRENT),
     &                         A_GRIDS(START_ADDRESS_FINER),
     &                         B_GRIDS(ADDRESS_CURRENT),
     &                         B_GRIDS(START_ADDRESS_FINER),
     &                         C1_GRIDS(ADDRESS_CURRENT),
     &                         C1_GRIDS(START_ADDRESS_FINER),
     &                         C2_GRIDS(ADDRESS_CURRENT),
     &                         C2_GRIDS(START_ADDRESS_FINER),
     &                         DEF_GRIDS(ADDRESS_CURRENT),
     &                         DEF_GRIDS(START_ADDRESS_FINER),
     &                         D_GRIDS(ADDRESS_CURRENT),
     &                         D_GRIDS(START_ADDRESS_FINER),
     &                         E_GRIDS(ADDRESS_CURRENT),
     &                         E_GRIDS(START_ADDRESS_FINER),
     &                         F_GRIDS(ADDRESS_CURRENT),
     &                         F_GRIDS(START_ADDRESS_FINER),
     &                         G_GRIDS(ADDRESS_CURRENT),
     &                         G_GRIDS(START_ADDRESS_FINER),
     &                         COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         COS_P_GRIDS(START_ADDRESS_FINER_2D),
     &                         SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                         COS_V_GRIDS(START_ADDRESS_FINER_2D),
     &                         Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                         Z_Q_GRIDS(START_ADDRESS_FINER_Z),
     &                         Z_MID_GRIDS(ADDRESS_CURRENT_Z),
     &                         Z_MID_GRIDS(START_ADDRESS_FINER_Z),
     &                         L_FIRST_DESCENT,
     &                         I_CURRENT,J_CURRENT,K_CURRENT,
     &                         I_FINER,J_FINER,K_FINER,I_NT,J_NT,K_NT,
     &                         KREST,EARTH_RADIUS_INVERSE,
     &                         LATITUDE_STEP_GRIDS(GRID_CURRENT-1),
     &                         LONGITUDE_STEP_GRIDS(GRID_CURRENT-1),
     &                         RES_DIRS(GRID_CURRENT),VERSION)

CL    END LOOP OVER GRIDS
            END DO

            L_FIRST_DESCENT = .FALSE.

C ----------------------------------------------------------------------
CL    SECTION 3.1.4 COARSEST GRID SOLUTION
C ----------------------------------------------------------------------

            DO I=1,NCOARSE
              CALL MG_SMOOTH(Q_GRIDS(ADDRESS_CURRENT),
     &                       RHS_GRIDS(ADDRESS_CURRENT),
     &                       A_GRIDS(ADDRESS_CURRENT),
     &                       B_GRIDS(ADDRESS_CURRENT),
     &                       C1_GRIDS(ADDRESS_CURRENT),
     &                       C2_GRIDS(ADDRESS_CURRENT),
     &                       DEF_GRIDS(ADDRESS_CURRENT),
     &                       D_GRIDS(ADDRESS_CURRENT),
     &                       E_GRIDS(ADDRESS_CURRENT),
     &                       F_GRIDS(ADDRESS_CURRENT),
     &                       G_GRIDS(ADDRESS_CURRENT),
     &                       COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                       SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                       COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                       I_CURRENT,J_CURRENT,K_CURRENT,
     &                       W,GRID_CURRENT,RMS_RES,IPRINT,
     &                       SWJAC,SWSYM,JAC,PAT,SYM,ILINE,JLINE,
     &                       KLINE,NFB,IPAT,
     &                       LATITUDE_STEP_GRIDS(1),
     &                       LONGITUDE_STEP_GRIDS(1),
     &                       EARTH_RADIUS_INVERSE,VERSION,K_BC,
     &                       I_NT,J_NT,K_NT,
     &                       Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                       Z_MID_GRIDS(ADDRESS_CURRENT_Z))

            END DO

C ----------------------------------------------------------------------
CL    SECTION 3.2 ASCENDING PART OF MULTI-GRID CYCLE.
CL                IF A W-CYCLE THEN ASCENDING PART TERMINATES AT
CL                A NEW GRID LEVEL AND CONTROL RETURNS TO DESCENDING
CL                SECTION 3.1
C ----------------------------------------------------------------------

            L_ASCEND = .TRUE.
            START_GRID = 1
            DO WHILE (L_ASCEND)

C ----------------------------------------------------------------------
CL    SECTION 3.2.1 PROLONG FROM COARSE GRID TO FINER GRID.
C ----------------------------------------------------------------------

              I_COARSER = I_CURRENT
              J_COARSER = J_CURRENT
              K_COARSER = K_CURRENT
              START_ADDRESS_COARSER = ADDRESS_CURRENT
              START_GRID = START_GRID+1
              I_CURRENT = IMAX(START_GRID)
              J_CURRENT = JMAX(START_GRID)
              K_CURRENT = KMAX(START_GRID)
              ADDRESS_CURRENT = START_ADDRESS(START_GRID)
              ADDRESS_CURRENT_2D = START_ADDRESS_2D(START_GRID)
              ADDRESS_CURRENT_Z = START_ADDRESS_Z(START_GRID)

CL    CALL PROLONG

              CALL MG_PROLONG(Q_GRIDS(START_ADDRESS_COARSER),
     &                        Q_GRIDS(ADDRESS_CURRENT),
     &                        I_COARSER,J_COARSER,K_COARSER,
     &                        I_CURRENT,J_CURRENT,K_CURRENT,
     &                        RES_DIRS(START_GRID),VERSION,K_BC,
     &                        Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                        I_NT,J_NT,K_NT)

C ----------------------------------------------------------------------
CL    SECTION 3.2.2 POST-SMOOTHING.
C ----------------------------------------------------------------------

CL    CALL SMOOTH

              DO I = 1,NPOST
                CALL MG_SMOOTH(Q_GRIDS(ADDRESS_CURRENT),
     &                         RHS_GRIDS(ADDRESS_CURRENT),
     &                         A_GRIDS(ADDRESS_CURRENT),
     &                         B_GRIDS(ADDRESS_CURRENT),
     &                         C1_GRIDS(ADDRESS_CURRENT),
     &                         C2_GRIDS(ADDRESS_CURRENT),
     &                         DEF_GRIDS(ADDRESS_CURRENT),
     &                         D_GRIDS(ADDRESS_CURRENT),
     &                         E_GRIDS(ADDRESS_CURRENT),
     &                         F_GRIDS(ADDRESS_CURRENT),
     &                         G_GRIDS(ADDRESS_CURRENT),
     &                         COS_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         SEC_P_GRIDS(ADDRESS_CURRENT_2D),
     &                         COS_V_GRIDS(ADDRESS_CURRENT_2D),
     &                         I_CURRENT,J_CURRENT,K_CURRENT,
     &                         W,START_GRID,RMS_RES,IPRINT,
     &                         SWJAC,SWSYM,JAC,PAT,SYM,ILINE,JLINE,
     &                         KLINE,NFB,IPAT,
     &                         LATITUDE_STEP_GRIDS(START_GRID),
     &                         LONGITUDE_STEP_GRIDS(START_GRID),
     &                         EARTH_RADIUS_INVERSE,VERSION,K_BC,
     &                         I_NT,J_NT,K_NT,
     &                         Z_Q_GRIDS(ADDRESS_CURRENT_Z),
     &                         Z_MID_GRIDS(ADDRESS_CURRENT_Z))
              END DO

C ----------------------------------------------------------------------
CL    SECTION 3.2.3 DECIDE WHETHER :
CL                  1. FINEST GRID REACHED SO END MULTI-GRID ITERATION.
CL                  2. W-CYCLE REACHED NEW GRID-LEVEL SO GO BACK TO
CL                     DESCENT PHASE (SECTION 3.1)
CL                  3. IF NOT 1 OR 2 THEN CONTINUE WITH ASCENT PHASE.
C ----------------------------------------------------------------------

CL    CHECK TO SEE IF ON FINEST GRID.
              IF (START_GRID.EQ.NGRIDS) THEN
C NO MORE ASCENT STEPS.
                L_ASCEND  = .FALSE.
C NO MORE DESCENT STEPS.
                L_DESCEND = .FALSE.
CL    CHECK TO SEE IF A W CYCLE REACHING GRID 2 FOR THE
CL    FIRST TIME.
              ELSE IF (START_GRID.EQ.2.AND.NCGC.EQ.2.AND.
     &                 IDGC(2).LT.NCGC)THEN
C NO MORE ASCENT STEPS.
                L_ASCEND  = .FALSE.
C BUT DESCENT STEPS WILL RESTART AT GRID 2
CL    CHECK TO SEE IF FULL MULTI-GRID AND REACHING A NEW GRID-LEVEL.
              ELSE IF (NCGC.GT.2.AND.IDGC(START_GRID).EQ.1)THEN
C NO MORE ASCENT STEPS.
                L_ASCEND  = .FALSE.
C BUT DESCENT STEPS WILL RESTART AT START_GRID.
              END IF

CL    END ASCENDING LOOP.
            END DO

CL    END DESCENDING LOOP.
          END DO

C ----------------------------------------------------------------------
CL    SECTION 3.3.  ANALYSE CONVERGENCE BEHAVIOUR FOR LAST ITERATION
CL                  OF MULTI-GRID.
C ----------------------------------------------------------------------

CL    CALCULATE CONVERGENCE FIGURES AND OUTPUT AS REQUESTED.

          CONVERGENCE_RATE = RMS_RES / RMS_RES_OLD
          RMS_RES_OLD = RMS_RES
          PRACTICAL_SMOOTHING_RATE = CONVERGENCE_RATE**(1./(NPRE+NPOST))
          IF(MG_ITERATIONS.NE.1.AND.IPRINT.GE.2)
     &                  WRITE(6,9006) CONVERGENCE_RATE,
     &                                PRACTICAL_SMOOTHING_RATE

CL    CHECK FOR CONVERGENCE TO REQUIRED TOLERANCE.

          IF(RMS_RES.LT.TOL_RES*RMS_RES_INIT) THEN
            IF(IPRINT.GE.1) THEN
              WRITE(6,9003) MG_ITERATIONS
            ENDIF
CL    TERMINATE ITERATIONS IF CONVERGENCE ACHIEVED.
            L_ITERATE = .FALSE.
          ENDIF

CL    CHECK TO SEE IF MAXIMUM NUMBER OF ITERATIONS PERFORMED.

          IF (MG_ITERATIONS.EQ.MAXITS.AND.L_ITERATE) THEN
            L_ITERATE = .FALSE.
            IF(IPRINT.GE.1) WRITE(6,9004) MAXITS
          ENDIF

CL    CHECK TO SEE IF SMOOTHING FACTOR LESS THAN TOLERANCE.
          IF(MG_ITERATIONS.NE.1.AND.
     &          PRACTICAL_SMOOTHING_RATE.GT.WORST_SMOOTHING_RATE) THEN
            L_ITERATE = .FALSE.
            IF(IPRINT.GE.1) WRITE(6,9007) MG_ITERATIONS
            IF(IPRINT.GE.1) WRITE(6,9008) PRACTICAL_SMOOTHING_RATE,
     &                                    WORST_SMOOTHING_RATE
          END IF

CL    END LOOP OVER MULTI-GRID ITERATIONS.
        END DO

CL END CHECK ON NUMBER OF GRIDS.
      END IF

CL----------------------------------------------------------------------
CL    SECTION 4. FINAL SOLUTION OUTPUT IF REQUIRED.
CL               COPY SOLUTION BACK INTO SOLUTION ARRAY.
CL----------------------------------------------------------------------

CL    COPY SOLUTION.

      DO K= 1,K_LENGTH
        BASE = START_ADDRESS(NGRIDS) + (K-1)*I_LENGTH*J_LENGTH - 1
        DO J= 1,J_LENGTH
          I2 = (J-1)*I_LENGTH
          DO I= 1,I_LENGTH
            Q(I,J,K) = Q_GRIDS(BASE+I+I2)
          END DO
        END DO
      END DO

 9000 FORMAT(1X,//,' Convergence History',/,20('-'))
 9001 FORMAT(1X,//,' Converged Within Tolerance In ',I3,' Relaxations')
 9002 FORMAT(1X,//,' Failed To Converge In ',I3,' Relaxations')
 9003 FORMAT(1X,//,' Converged Within Tolerance In ',I3,' MG Itns')
 9004 FORMAT(1X,//,' Failed To Converge In ',I3,' MG Itns')
 9005 FORMAT(1X,/,' Multigrid Iteration No :',I3,/,28('~'))
 9006 FORMAT(1X,/,' MG Convergence Rate =',F5.3,/,
     &            ' Practical Smoothing Rate =',F5.3)
 9007 FORMAT(1X,/,' Multigrid stopped on iteration number ',I3)
 9008 FORMAT(1X,/,' Practical Smoothing Rate =',F5.3,
     &            ' worse than allowed rate =',F5.3)

CL    END OF ROUTINE MG_FDMG

      RETURN
      END
CLL   SUBROUTINE MG_GRIDSEQ
CLL
CLL   PURPOSE:
CLL   -------
CLL   DECIDES HOW MANY COARSE GRIDS MAY BE GENERATED FROM THE FINE
CLL   GRID SUPPLIED BY THE USER. DEFINES THE NUMBER OF NODES IN EACH
CLL   DIRECTION AND SETS STORAGE POINTERS FOR EACH GRID IN THE SEQUENCE
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_GRIDSEQ(I_NT,J_NT,K_NT,I_LENGTH,J_LENGTH,K_LENGTH,
     &                      MXNGRDS,NGRIDS,IMAX,JMAX,KMAX,
     &                      START_ADDRESS,LAST_ADDRESS,
     &                      LONGITUDE_STEP_INVERSE,
     &                      LATITUDE_STEP_INVERSE,
     &                      LONGITUDE_STEP_GRIDS,LATITUDE_STEP_GRIDS,
     &                      START_ADDRESS_2D,LAST_ADDRESS_2D,RES_DIRS,
     &                      VERSION,IPRINT,START_ADDRESS_Z,
     &                      LAST_ADDRESS_Z)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_LENGTH   !IN. NUMBER OF POINTS IN I DIRECTION.
     &, J_LENGTH   !IN. NUMBER OF POINTS IN J DIRECTION.
     &, K_LENGTH   !IN. NUMBER OF POINTS IN K DIRECTION.

      INTEGER
     &  MXNGRDS  !IN. MAXIMUM NUMBER OF GRIDS ALLOWED.
     &, IPRINT   !IN. CONTROLS PRINTING OPTIONS
     &, VERSION  !IN. VERSION OF MULTIGRID WANTED.

      REAL
     &  LONGITUDE_STEP_INVERSE !IN. 1./LONGITUDE STEP ON FINEST GRID.
     &, LATITUDE_STEP_INVERSE  !IN. 1./LATITUDE STEP ON FINEST GRID.

      REAL
     &  LONGITUDE_STEP_GRIDS(MXNGRDS) !OUT. LONGITUDE STEP INVERSE ON
     &                                !     EACH GRID.
     &, LATITUDE_STEP_GRIDS(MXNGRDS)  !OUT. LATITUDE STEP INVERSE ON
     &                                !     EACH GRID.

      INTEGER
     &  NGRIDS        !OUT. THE TOTAL NUMBER OF GRIDS IN THE SEQUENCE
     &, IMAX(MXNGRDS) !OUT. NUMBER OF NODES IN THE I-DIRECTION
     &                !     ON EACH GRID.
     &, JMAX(MXNGRDS) !OUT. NUMBER OF NODES IN THE J-DIRECTION
     &                !     ON EACH GRID.
     &, KMAX(MXNGRDS) !OUT. NUMBER OF NODES IN THE K-DIRECTION
     &                !     ON EACH GRID.
     &, START_ADDRESS(MXNGRDS) !OUT. START ADDRESS IN DATA ARRAY FOR
     &                         !     EACH GRID.
     &, RES_DIRS(MXNGRDS)    !OUT. CODE GIVING NUMBER OF DIRECTIONS
     &                       !   TO RESTRICT IN.
     &, LAST_ADDRESS  !OUT. LAST ADDRESS OF DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION SPACE REQUIRED.
     &, START_ADDRESS_2D(MXNGRDS) !OUT. START ADDRESS IN DATA ARRAY FOR
     &                         !     EACH 2-D GRID.
     &, LAST_ADDRESS_2D !OUT. LAST ADDRESS OF DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION 2-D SPACE REQUIRED.
     &, START_ADDRESS_Z(MXNGRDS) !OUT. START ADDRESS IN DATA ARRAY FOR
     &                         !     EACH Z GRID.
     &, LAST_ADDRESS_Z !OUT. LAST ADDRESS OF DATA ARRAY NEEDED. USED TO
     &                ! DIMENSION Z SPACE REQUIRED.

C*----------------------------------------------------------------------

C*L  LOCAL WORK ARRAYS REQUIRED.
      INTEGER
     & WORK(MXNGRDS)
C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  IMAXW    ! NUMBER OF NODES IN I DIMENSION.
     &, JMAXW    ! NUMBER OF NODES IN J DIMENSION.
     &, KMAXW    ! NUMBER OF NODES IN K DIMENSION.
     &, INEXT    ! NUMBER OF NODES ON NEXT GRID DOWN
     &, JNEXT    !             "
     &, KNEXT    !             "
     &, NNODESW  ! TOTAL NUMBER OF NODES.
     &, I,J,K
     &, EVAL     ! A NUMBER DESCRIBING THE DIRECTIONS TO RESTRICT IN
     &           ! CALCULATED AS I+J+K WHERE I=4 IF TRYING TO RESTRICT
     &           ! IN I DIRECTION, J=2, K=1 LIKEWISE.

CMHM THESE VARIABLES MIGHT BE BETTER AS PASSED IN.
      INTEGER
     &  MIN_I    !MINIMUM NUMBER OF POINTS ALLOWED IN I DIRECTION.
     &, MIN_J    !MINIMUM NUMBER OF POINTS ALLOWED IN J DIRECTION.
     &, MIN_K    !MINIMUM NUMBER OF POINTS ALLOWED IN K DIRECTION.

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. FIND THE NUMBER OF GRIDS WHICH MAY BE GENERATED
CL----------------------------------------------------------------------

C SET MAXIMUM NUMBER OF NODES IN EACH DIMENSION.

      IMAXW = I_LENGTH
      JMAXW = J_LENGTH
      KMAXW = K_LENGTH

      IF(VERSION.LT.3) THEN
CL GLOBAL VERSIONS.
C SET MINIMUM NUMBER OF POINTS ALLOWED ON EACH GRID.
        MIN_I = 3
        MIN_J = 5
        IF(VERSION .EQ. 2) MIN_J=9
        MIN_K = 3
C CALCULATE MAXIMUM NUMBER OF GRIDS POSSIBLE.
C CALCULATION DOES NOT CHECK FOR MORE LEGAL GRIDS THAN USER REQUESTED.
C CALCULATION ASSUMES PERIODICITY IN I DIRECTION AND EVEN NUMBER.

        NGRIDS=1
        I = 0
        J = 0
        K = 0
        IF (I_NT) I=4
        IF (J_NT) J=2
        IF (K_NT) K=1

C CHECK INPUT VALUES ARE LEGAL.

        IF(MOD(IMAXW,2).NE.0 .AND. I_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
          WRITE(6,*)' INPUT NUMBER OF I NODES WAS ODD, VALUE WAS ',IMAXW
          WRITE(6,*)' VALUE MUST BE EVEN FOR PERIODIC MULTIGRID '
          WRITE(6,*)' COARSENING IN I DIRECTION NOT ATTEMPTED.'
          I=0
        END IF
        IF(MOD(JMAXW,2).EQ.0 .AND. J_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
      WRITE(6,*)' INPUT NUMBER OF J NODES WAS EVEN, VALUE WAS ',JMAXW
          WRITE(6,*)' COARSENING IN J DIRECTION NOT ATTEMPTED.'
          J=0
        END IF
        IF(MOD(KMAXW,2).EQ.0 .AND. K_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
      WRITE(6,*)' INPUT NUMBER OF K NODES WAS EVEN, VALUE WAS ',KMAXW
          WRITE(6,*)' COARSENING IN K DIRECTION NOT ATTEMPTED.'
          K=0
        END IF

        EVAL = I + J + K

CL TEST ALL DIRECTIONS.
        DO WHILE (EVAL.EQ.7 .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          JNEXT=(JMAXW+1)/2
          KNEXT=(KMAXW+1)/2
          NNODESW=(INEXT+1)*JNEXT*KNEXT
C NEXT GRID ILLEGAL
          IF(MOD(INEXT,2).NE.0 .OR. INEXT.LT.MIN_I ) THEN
            EVAL = 3
          ELSE IF (JNEXT.LT.MIN_J .OR.
     &              MOD(JNEXT,2).EQ.0 ) THEN
            EVAL = 5
          ELSE IF  (KNEXT.LT.MIN_K) THEN
            EVAL = 6
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF (MOD(INEXT+1,2).EQ.0) I = 0
            IF (MOD(JNEXT,2).EQ.0) J = 0
            IF (MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
            KMAXW = KNEXT
          END IF
        END DO

C RESTRICT IN I AND J DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.6  .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          JNEXT=(JMAXW+1)/2
          NNODESW=(INEXT+1)*JNEXT*KMAXW
C NEXT GRID ILLEGAL
          IF(MOD(INEXT,2).NE.0 .OR. INEXT.LT.MIN_I ) THEN
            EVAL = 2
          ELSE IF(JNEXT.LT.MIN_J .OR.
     &             MOD(JNEXT,2).EQ.0 ) THEN
            EVAL = 4
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(INEXT+1,2).EQ.0) I = 0
            IF(MOD(JNEXT,2).EQ.0) J = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
          END IF
        END DO

C RESTRICT IN I AND K DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.5  .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          KNEXT = (KMAXW+1)/2
          NNODESW=(INEXT+1)*JMAXW*KNEXT
C NEXT GRID ILLEGAL
          IF(MOD(INEXT,2).NE.0 .OR. INEXT.LT.MIN_I ) THEN
            EVAL = 1
          ELSE IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 4
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(INEXT+1,2).EQ.0) I = 0
            IF(MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            KMAXW = KNEXT
          END IF
        END DO

C TRY AND RESTRICT IN I DIRECTION ONLY.

        DO WHILE (EVAL.EQ.4 .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          NNODESW=(INEXT+1)*JMAXW*KMAXW
C NEXT GRID ILLEGAL
          IF(MOD(INEXT,2).NE.0.OR.INEXT.LT.MIN_I) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
          END IF
        END DO

C RESTRICT IN J AND K DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.3  .AND. NGRIDS.LT.MXNGRDS)
          JNEXT=(JMAXW+1)/2
          KNEXT=(KMAXW+1)/2
          NNODESW=(IMAXW+1)*JNEXT*KNEXT
C NEXT GRID ILLEGAL
          IF(JNEXT.LT.MIN_J .OR.
     &        MOD(JNEXT,2).EQ.0 ) THEN
            EVAL = 1
          ELSE IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 2
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(JNEXT,2).EQ.0) J = 0
            IF(MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
            KMAXW = KNEXT
          END IF
        END DO

C TRY AND RESTRICT IN J DIRECTION ONLY.

        DO WHILE (EVAL.EQ.2 .AND. NGRIDS.LT.MXNGRDS)
          JNEXT = (JMAXW+1)/2
          NNODESW=(IMAXW+1)*JNEXT*KMAXW
C NEXT GRID ILLEGAL
          IF(JNEXT.LT.MIN_J .OR.
     &        MOD(JNEXT,2).EQ.0 ) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
          END IF
        END DO

C TRY AND RESTRICT IN K DIRECTION ONLY.

        DO WHILE (EVAL.EQ.1 .AND. NGRIDS.LT.MXNGRDS)
          KNEXT = (KMAXW+1)/2
          NNODESW=(IMAXW+1)*JMAXW*KNEXT
C NEXT GRID ILLEGAL
          IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            KMAXW = KNEXT
          END IF
        END DO

      ELSE
CL LAM VERSIONS
C SET MINIMUM NUMBER OF POINTS ALLOWED ON EACH GRID.
        MIN_I = 3
        MIN_J = 3
        MIN_K = 3
C CALCULATE MAXIMUM NUMBER OF GRIDS POSSIBLE.
C CALCULATION DOES NOT CHECK FOR MORE LEGAL GRIDS THAN USER REQUESTED.
C CALCULATION ASSUMES NON-PERIODICITY IN I DIRECTION.

        NGRIDS=1
        I = 0
        J = 0
        K = 0
        IF (I_NT) I=4
        IF (J_NT) J=2
        IF (K_NT) K=1

C CHECK INPUT VALUES ARE LEGAL.

        IF(MOD(IMAXW,2).EQ.0 .AND. I_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
      WRITE(6,*)' INPUT NUMBER OF I NODES WAS EVEN, VALUE WAS ',IMAXW
          WRITE(6,*)' COARSENING IN I DIRECTION NOT ATTEMPTED.'
          I=0
        END IF
        IF(MOD(JMAXW,2).EQ.0 .AND. J_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
      WRITE(6,*)' INPUT NUMBER OF J NODES WAS EVEN, VALUE WAS ',JMAXW
          WRITE(6,*)' COARSENING IN J DIRECTION NOT ATTEMPTED.'
          J=0
        END IF
        IF(MOD(KMAXW,2).EQ.0 .AND. K_NT) THEN
          WRITE(6,*)' GRIDSEQ: ILLEGAL INPUT DATA FOR MULTI-GRID '
      WRITE(6,*)' INPUT NUMBER OF K NODES WAS EVEN, VALUE WAS ',KMAXW
          WRITE(6,*)' COARSENING IN K DIRECTION NOT ATTEMPTED.'
          K=0
        END IF

        EVAL = I + J + K

CL TEST ALL DIRECTIONS.
        DO WHILE (EVAL.EQ.7 .AND. NGRIDS.LT.MXNGRDS)
          INEXT=(IMAXW+1)/2
          JNEXT=(JMAXW+1)/2
          KNEXT=(KMAXW+1)/2
          NNODESW=INEXT*JNEXT*KNEXT
C NEXT GRID ILLEGAL
          IF( INEXT.LT.MIN_I ) THEN
            EVAL = 3
          ELSE IF ( JNEXT.LT.MIN_J .OR. MOD(JNEXT,2).EQ.0) THEN
            EVAL = 5
          ELSE IF  (KNEXT.LT.MIN_K) THEN
            EVAL = 6
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF (MOD(INEXT,2).EQ.0) I = 0
            IF (MOD(JNEXT,2).EQ.0) J = 0
            IF (MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
            KMAXW = KNEXT
          END IF
        END DO

C RESTRICT IN I AND J DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.6  .AND. NGRIDS.LT.MXNGRDS)
          INEXT=(IMAXW+1)/2
          JNEXT=(JMAXW+1)/2
          NNODESW=INEXT*JNEXT*KMAXW
C NEXT GRID ILLEGAL
          IF(INEXT.LT.MIN_I ) THEN
            EVAL = 2
          ELSE IF ( JNEXT.LT.MIN_J .OR. MOD(JNEXT,2).EQ.0) THEN
            EVAL = 4
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(INEXT,2).EQ.0) I = 0
            IF(MOD(JNEXT,2).EQ.0) J = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            JMAXW = JNEXT
          END IF
        END DO

C RESTRICT IN I AND K DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.5  .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          KNEXT = (KMAXW+1)/2
          NNODESW= INEXT*JMAXW*KNEXT
C NEXT GRID ILLEGAL
          IF(INEXT.LT.MIN_I ) THEN
            EVAL = 1
          ELSE IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 4
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(INEXT,2).EQ.0) I = 0
            IF(MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
            KMAXW = KNEXT
          END IF
        END DO

C TRY AND RESTRICT IN I DIRECTION ONLY.

        DO WHILE (EVAL.EQ.4 .AND. NGRIDS.LT.MXNGRDS)
          INEXT = IMAXW/2
          NNODESW=INEXT*JMAXW*KMAXW
C NEXT GRID ILLEGAL
          IF(INEXT.LT.MIN_I) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            IMAXW = INEXT
          END IF
        END DO

C RESTRICT IN J AND K DIRECTIONS ONLY.

        DO WHILE (EVAL.EQ.3  .AND. NGRIDS.LT.MXNGRDS)
          JNEXT=(JMAXW+1)/2
          KNEXT=(KMAXW+1)/2
          NNODESW=IMAXW*JNEXT*KNEXT
C NEXT GRID ILLEGAL
          IF ( JNEXT.LT.MIN_J .OR. MOD(JNEXT,2).EQ.0) THEN
            EVAL = 1
          ELSE IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 2
C NEXT GRID LEGAL BUT THAT'S THE LAST
          ELSE IF(MOD(NNODESW,2).EQ.0) THEN
            WORK( NGRIDS) = EVAL
            IF(MOD(JNEXT,2).EQ.0) J = 0
            IF(MOD(KNEXT,2).EQ.0) K = 0
            EVAL = I + J + K
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
            KMAXW = KNEXT
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
            KMAXW = KNEXT
          END IF
        END DO

C TRY AND RESTRICT IN J DIRECTION ONLY.

        DO WHILE (EVAL.EQ.2 .AND. NGRIDS.LT.MXNGRDS)
          JNEXT = (JMAXW+1)/2
          NNODESW=IMAXW*JNEXT*KMAXW
C NEXT GRID ILLEGAL
          IF ( JNEXT.LT.MIN_J .OR. MOD(JNEXT,2).EQ.0) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            JMAXW = JNEXT
          END IF
        END DO

C TRY AND RESTRICT IN K DIRECTION ONLY.

        DO WHILE (EVAL.EQ.1 .AND. NGRIDS.LT.MXNGRDS)
          KNEXT = (KMAXW+1)/2
          NNODESW=IMAXW*JMAXW*KNEXT
C NEXT GRID ILLEGAL
          IF(KNEXT.LT.MIN_K ) THEN
            EVAL = 0
C NEXT GRID LEGAL AND POSSIBLE OTHERS TO FOLLOW
          ELSE
            WORK( NGRIDS) = EVAL
            NGRIDS = NGRIDS + 1
            KMAXW = KNEXT
          END IF
        END DO

      END IF

CL--- ------------------------------------------------------------------
CL    SECTION 2. CALCULATE GRID-SIZES AND ADDRESS INFORMATION.
CL--- ------------------------------------------------------------------

CFPP$ NOCONCUR
      DO I = NGRIDS,2,-1
        RES_DIRS(I) = WORK(NGRIDS+1-I)
      END DO
      RES_DIRS(1) = 0.

CL    DEFINE THE PARAMETERS OF THE COARSEST GRID: IGRID=1

      START_ADDRESS(1)=1
      START_ADDRESS_2D(1)=1
      START_ADDRESS_Z(1)=1
      IMAX(1)=IMAXW
      JMAX(1)=JMAXW
      KMAX(1)=KMAXW

CL    FIND NUMBER OF NODES AND STORAGE POINTERS FOR THE REMAINING GRIDS

      DO I=2,NGRIDS
        J= I-1
        IMAX(I) = IMAX(J)
        JMAX(I) = JMAX(J)
        KMAX(I) = KMAX(J)
        IF(VERSION.LT.3) THEN
CL GLOBAL
          IF(RES_DIRS(I).GT.3) IMAX(I)=2*IMAX(J)
        ELSE
CL LAM
          IF(RES_DIRS(I).GT.3) IMAX(I)=2*IMAX(J)-1
        END IF
        IF(MOD(RES_DIRS(I),2).EQ.1) KMAX(I)=2*KMAX(J)-1
        IF(RES_DIRS(I).EQ.2 .OR.
     &     RES_DIRS(I).EQ.3 .OR.
     &     RES_DIRS(I).EQ.6 .OR.
     &     RES_DIRS(I).EQ.7 ) JMAX(I) = 2*JMAX(J)-1
        START_ADDRESS(I)= START_ADDRESS(J)+IMAX(J)*JMAX(J)*KMAX(J)
        START_ADDRESS_2D(I)= START_ADDRESS_2D(J)+IMAX(J)*JMAX(J)
        START_ADDRESS_Z(I)= START_ADDRESS_Z(J)+KMAX(J)+1
      END DO
      LAST_ADDRESS = START_ADDRESS(NGRIDS) + IMAX(NGRIDS)*JMAX(NGRIDS)
     &               *KMAX(NGRIDS) - 1
      LAST_ADDRESS_2D = START_ADDRESS_2D(NGRIDS) +
     &                  IMAX(NGRIDS)*JMAX(NGRIDS) - 1
      LAST_ADDRESS_Z = START_ADDRESS_Z(NGRIDS) +
     &                  KMAX(NGRIDS)

CL    SET HORIZONTAL GRID SPACING ON EACH GRID.

      LATITUDE_STEP_GRIDS(NGRIDS) = LATITUDE_STEP_INVERSE
      LONGITUDE_STEP_GRIDS(NGRIDS) = LONGITUDE_STEP_INVERSE
      DO I = NGRIDS-1,1,-1
        LATITUDE_STEP_GRIDS(I) = LATITUDE_STEP_GRIDS(I+1)
        LONGITUDE_STEP_GRIDS(I) = LONGITUDE_STEP_GRIDS(I+1)
        IF(RES_DIRS(I+1).GT.3)
     &     LONGITUDE_STEP_GRIDS(I) = LONGITUDE_STEP_GRIDS(I+1)*.5
        IF(RES_DIRS(I+1).EQ.2 .OR.
     &     RES_DIRS(I+1).EQ.3 .OR.
     &     RES_DIRS(I+1).EQ.6 .OR.
     &     RES_DIRS(I+1).EQ.7 )
     &     LATITUDE_STEP_GRIDS(I) = LATITUDE_STEP_GRIDS(I+1)*.5
      END DO

CL    END OF ROUTINE MG_GRIDSEQ

      RETURN
      END
CLL   SUBROUTINE MG_I_LINE
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS LINE SMOOTHING IN I DIRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_I_LINE(
     &                     Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                     SEC_P_LATITUDE,
     &                     COS_V_LATITUDE,LATITUDE_STEP_INVERSE,
     &                     LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                     I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                     SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                     VERSION,K_BC,Z_Q,Z_MID,J_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, NFB         !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &              ! NEEDED.
     &, IPAT        !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &              ! POINTS.
     &, VERSION     !IN. VERSION OF MULTIGRID SCHEME.
     &, K_BC        !IN. BOUNDARY CONDITION IN K DIRECTION.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES

      REAL
     &  SWJAC    !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM    !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
     &, RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.
     &, W        !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM

      REAL
     &  A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, LONGITUDE_STEP_INVERSE !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   9 LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  LHS(I_DIM)
     &, A0(I_DIM,J_DIM)     !\
     &, A1(I_DIM,J_DIM)     !  \
     &, A2(I_DIM,J_DIM)     !   > ARRAYS USED IN TRI-DIAGONAL SOLVER
     &, X(I_DIM,J_DIM)      !  /
     &, R(I_DIM,J_DIM)      !/
     &, F_VECTOR(I_DIM-1,J_DIM) !USED IN PERIODIC TRIDIAGONAL SOLVER
     &, G_VECTOR(2,J_DIM)       !USED IN PERIODIC TRIDIAGONAL SOLVER
     &, CORRECTIONS(I_DIM,J_DIM,K_DIM)
     &, AN(J_DIM)           !USED IN PERIODIC TRIDIAGONAL SOLVER

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I_BEGIN
     &, I_END
     &, J_BEGIN
     &, JK_BEGIN
     &, J_END
     &, J_STEP
     &, K_BEGIN
     &, K_END
     &, K_STEP
     &, I,J,K
     &, I_FB
     &, I_ZEBRA

      REAL
     &  WEIGHT
     &, FACTOR


C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_I_MATRIX

C*----------------------------------------------------------------------

CFPP$ NOCONCUR R

CL----------------------------------------------------------------------
CL    SECTION 1. INITIALISATION.
CL----------------------------------------------------------------------

CL    INITIALISE NORMS

      RMS_RES=0.0
      RMS_INC=0.0

CL    SET TEMPORARY RELAXATION PARAMETER TO ZERO FOR JACOBI:

      WEIGHT = SWJAC*W

CL    DEFINE POINT ORDERING IN JK PLANE FOR FIRST SWEEP

      IF (.NOT.J_NT) THEN
        J_BEGIN = 1
        J_END = 1
      ELSE IF (VERSION.EQ.3) THEN
        J_BEGIN = 1
        J_END   = J_DIM
      ELSE
C DO NOT DO BOUNDARY IN OTHER VERSIONS.
        J_BEGIN = 2
        J_END   = J_DIM-1
      END IF
      IF (IPAT.EQ.2) THEN
        J_STEP  = 2
      ELSE
        J_STEP  = 1
      END IF
      K_BEGIN = 1
      K_END   = K_DIM
      IF (K_BC.EQ.2 ) THEN
        K_END = K_DIM-1
      ELSE IF (K_BC.EQ.3 ) THEN
        K_BEGIN = 2
      ELSE IF (K_BC.EQ. 4) THEN
        K_BEGIN = 2
        K_END = K_DIM-1
      END IF
      K_STEP  = 1
      I_BEGIN=1
      I_END=I_DIM
      IF (VERSION.EQ.4) THEN
C DO NOT USE BOUNDARIES WHEN THEY ARE DIRICHLET BOUNDARY CONDITIONS
        I_BEGIN = 2
        I_END = I_DIM-1
      END IF

CL    SET ZEBRA SWITCH FOR FIRST SWEEP

      I_ZEBRA=0

CL----------------------------------------------------------------------
CL    SECTION 2. PERFORM APPROPRIATE LINE SMOOTHING.
CL----------------------------------------------------------------------

CL    FORWARD / BACKWARD SWEEP

      DO I_FB= 1,NFB

        RMS_RES = SWSYM*RMS_RES
        RMS_INC = SWSYM*RMS_INC

CL    SWEEP OVER JK PLANE IN PRESCRIBED POINT ORDER

        DO K = K_BEGIN,K_END,K_STEP

          IF (IPAT.EQ.2) THEN
CL    START EACH LINE ON A DIFFERENT COLOUR.
            IF(MOD(J_BEGIN+K+I_ZEBRA,IPAT).NE.0 .AND. I_FB.EQ.1) THEN
              JK_BEGIN = J_BEGIN + 1
            ELSE IF(MOD(J_BEGIN+K+I_ZEBRA,IPAT).NE.0 .AND.
     &              I_FB.EQ.2) THEN
              JK_BEGIN = J_BEGIN - 1
            ELSE
              JK_BEGIN = J_BEGIN
            END IF
          ELSE
            JK_BEGIN = J_BEGIN
          END IF

          IF (VERSION .EQ. 2) THEN
            DO J = JK_BEGIN,J_END,J_STEP

              IF(J.EQ.(J_DIM+1)/2) THEN
C SET ZERO VALUES EXCEPT FOR MATRIX DIAGONAL TO ENABLE COMPLETE
C VECTORISATION OF MATRIX SOLVERS.
                DO I=1,I_DIM
                  R(I,J) = 0.
                  A0(I,J) = 1.
                  A1(I,J) = 0.
                  A2(I,J) = 0.
                END DO
                DO I=1,I_DIM-1
                  F_VECTOR(I,J) = 0.
                END DO
                G_VECTOR(1,J) = 0.
                G_VECTOR(2,J) = 0.
                AN(J) = 1

              ELSE
CL    ASSEMBLE TRIDIAGONAL SYSTEM TO INVERT FOR LINE VALUES

                CALL MG_I_MATRIX(LHS,A0(1,J),A1(1,J),A2(1,J),
     &                           AN(J),F_VECTOR(1,J),G_VECTOR(1,J),
     &                           Q,A,B,C1,C2,DEF,D,E,F,G,
     &                           I_DIM,J_DIM,K_DIM,
     &                           COS_P_LATITUDE,SEC_P_LATITUDE,
     &                           COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                           LATITUDE_STEP_INVERSE,
     &                           LONGITUDE_STEP_INVERSE,J,K,VERSION,
     &                           Z_Q,Z_MID,J_NT,K_NT)

                DO I=I_BEGIN,I_END
                  R(I,J) = LHS(I)-RHS(I,J,K)
                  RMS_RES = RMS_RES+R(I,J)*R(I,J)
                END DO

              END IF
            END DO
          ELSE
            DO J = JK_BEGIN,J_END,J_STEP
CL    ASSEMBLE TRIDIAGONAL SYSTEM TO INVERT FOR LINE VALUES

                CALL MG_I_MATRIX(LHS,A0(1,J),A1(1,J),A2(1,J),
     &                           AN(J),F_VECTOR(1,J),G_VECTOR(1,J),
     &                           Q,A,B,C1,C2,DEF,D,E,F,G,
     &                           I_DIM,J_DIM,K_DIM,
     &                           COS_P_LATITUDE,SEC_P_LATITUDE,
     &                           COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                           LATITUDE_STEP_INVERSE,
     &                           LONGITUDE_STEP_INVERSE,J,K,VERSION,
     &                           Z_Q,Z_MID,J_NT,K_NT)

                DO I=I_BEGIN,I_END
                  R(I,J) = LHS(I)-RHS(I,J,K)
                  RMS_RES = RMS_RES+R(I,J)*R(I,J)
                END DO

            END DO
          END IF

          IF (VERSION.LT.3) THEN
CL    GLOBAL VERSION.  SOLVE CYCLIC TRIDIAGONAL SYSTEM FOR CORRECTIONS

C CREATE V_VECTOR IN R(I)
C AND CREATE U_VECTOR FOR IN F_VECTOR
C REDUCE MATRIX TO UPPER DIAGONAL FORM.

            DO J = JK_BEGIN,J_END,J_STEP
              A0(1,J) = 1./A0(1,J)
            END DO
            DO I=2,I_DIM-1
              DO J = JK_BEGIN,J_END,J_STEP
                FACTOR = A2(I,J) * A0(I-1,J)
                A0(I,J) = 1./(A0(I,J) - FACTOR*A1(I-1,J))
                R(I,J) = R(I,J) - FACTOR*R(I-1,J)
                F_VECTOR(I,J) = F_VECTOR(I,J) - FACTOR*F_VECTOR(I-1,J)
              END DO
            END DO

C BACK SUBSTITUTE TO GET SOLUTION.

            DO J = JK_BEGIN,J_END,J_STEP
              R(I_DIM-1,J) = A0(I_DIM-1,J)*R(I_DIM-1,J)
              F_VECTOR(I_DIM-1,J) = A0(I_DIM-1,J)*F_VECTOR(I_DIM-1,J)
            END DO
            DO I= I_DIM-2,1,-1
              DO J = JK_BEGIN,J_END,J_STEP
                R(I,J) = A0(I,J)*(R(I,J)-A1(I,J)*R(I+1,J))
                F_VECTOR(I,J) = A0(I,J)*(F_VECTOR(I,J)-
     &                                   A1(I,J)*F_VECTOR(I+1,J))
              END DO
            END DO

C SOLVE FOR LAST ELEMENT OF FIELD.

            DO J = JK_BEGIN,J_END,J_STEP
              X(I_DIM,J) = (R(I_DIM,J)-G_VECTOR(1,J)*R(1,J)-
     &                        G_VECTOR(2,J)*R(I_DIM-1,J))/
     &                        (AN(J)-G_VECTOR(1,J)*F_VECTOR(1,J)-
     &                         G_VECTOR(2,J)*F_VECTOR(I_DIM-1,J))
            END DO

C SOLVE FOR ALL OTHER ELEMENTS OF FIELD.

            DO J = JK_BEGIN,J_END,J_STEP
              DO I=1,I_DIM-1
                X(I,J) = R(I,J) - F_VECTOR(I,J)* X(I_DIM,J)
              END DO
            END DO

          ELSE
CL    LIMITED AREA CODE. SOLVES NON-CYCLIC TRI-DIAGONAL SYSTEM.
CL    SOLVE TRIDIAGONAL SYSTEM FOR CORRECTIONS X
C REDUCE MATRIX TO UPPER DIAGONAL FORM.

            DO J = JK_BEGIN,J_END,J_STEP
              A0(I_BEGIN,J) = 1./A0(I_BEGIN,J)
            END DO
            DO I=I_BEGIN+1,I_END
              DO J = JK_BEGIN,J_END,J_STEP
                FACTOR = A2(I,J) * A0(I-1,J)
                A0(I,J) = 1./(A0(I,J) - FACTOR*A1(I-1,J))
                R(I,J) = R(I,J) - FACTOR*R(I-1,J)
              END DO
            END DO

C BACK SUBSTITUTE TO GET SOLUTION.

            DO J = JK_BEGIN,J_END,J_STEP
              X(I_END,J) = A0(I_END,J)*R(I_END,J)
            END DO
            DO I= I_END-1,1,-1
              DO J = JK_BEGIN,J_END,J_STEP
                X(I,J) = A0(I,J)*(R(I,J)-A1(I,J)*X(I+1,J))
              END DO
            END DO

          END IF

CL    STORE CORRECTIONS IN WORKSPACE
CL    AND ADD TO CURRENT SOLUTION UNLESS WT=0 (JACOBI)

          DO J = JK_BEGIN,J_END,J_STEP
            DO I=I_BEGIN,I_END
              Q(I,J,K)=Q(I,J,K)-WEIGHT*X(I,J)
              RMS_INC=RMS_INC+X(I,J)*X(I,J)
            END DO
          END DO
          IF (JAC) THEN
            DO J = JK_BEGIN,J_END,J_STEP
              DO I=I_BEGIN,I_END
                CORRECTIONS(I,J,K)=X(I,J)
              END DO
            END DO
          END IF

CL    END OF JK PLANE SWEEP

        END DO

CL    ADD CORRECTIONS FOR THE JACOBI METHOD

        IF(JAC) THEN
          DO K=K_BEGIN,K_END
            DO J=J_BEGIN,J_END
              IF(.NOT. (J.EQ.(J_DIM+1)/2 .AND. VERSION.EQ.2)) THEN
                DO I=I_BEGIN,I_END
                  Q(I,J,K)=Q(I,J,K)-W*CORRECTIONS(I,J,K)
                END DO
              END IF
            END DO
          END DO
        ENDIF

CL    RESET POINT ORDERING OF JK PLANE FOR SECOND SWEEP

        J = J_BEGIN
        J_BEGIN = J_END
        J_END   =  J
        J_STEP  = -J_STEP
        K = K_BEGIN
        K_BEGIN = K_END
        K_END   =  K
        K_STEP  = -1

CL    RESET ZEBRA SWITCH FOR SECOND SWEEP

        I_ZEBRA=1

      END DO

CL    END OF ROUTINE MG_I_LINE

      RETURN
      END
CLL   SUBROUTINE MG_I_MATRIX
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALCULATES CURRENT VALUE OF LEFT-HAND-SIDE OF P.D.E. FOR A GIVEN
CLL   VALUE OF J AND K. SETS MATRIX ENTRIES TO SOLVE FOR CORRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_I_MATRIX(
     &                       LHS,A0,A1,A2,AN,F_VECTOR,G_VECTOR,
     &                       Q,A,B,C1,C2,DEF,D,E,F,G,I_DIM,J_DIM,K_DIM,
     &                       COS_P_LATITUDE,SEC_P_LATITUDE,
     &                       COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                       LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,
     &                       J,K,VERSION,Z_Q,Z_MID,J_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM         !IN. NUMBER OF NODES IN THE I-DIRECTION
     &, J_DIM         !IN. NUMBER OF NODES IN THE J-DIRECTION
     &, K_DIM         !IN. NUMBER OF NODES IN THE K-DIRECTION
     &, J             !IN. VALUE OF J FOR WHICH LHS AND COEFFS ARE REQ.
     &, K             !IN. VALUE OF K FOR WHICH LHS AND COEFFS ARE REQ.
     &, VERSION       !IN. VERSION OF MULTIGRID SCHEME.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !IN. SOLUTION
     &, A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)!IN. COSINE OF LATITUDE AT B POINTS
     &, SEC_P_LATITUDE(I_DIM,J_DIM)!IN. 1./COSINE OF LATITUDE AT Q POINT
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  EARTH_RADIUS_INVERSE    !IN. 1/EARTH RADIUS
     &, LATITUDE_STEP_INVERSE   !IN. 1/LATITUDE STEP LENGTH IN RADIANS
     &, LONGITUDE_STEP_INVERSE  !IN. 1/LONGITUDE STEP LENGTH IN RADIANS

      REAL
     &  LHS(I_DIM)         !OUT. LEFT-HAND-SIDE.
     &, A0(I_DIM)          !OUT. DIAGONAL MATRIX COEFFICIENT.
     &, A1(I_DIM)          !OUT. SUPER-DIAGONAL MATRIX COEFFICIENT.
     &, A2(I_DIM)          !OUT. SUB-DIAGONAL MATRIX COEFFICIENT.
     &, AN                 !OUT. LAST VALUE OF DIAGONAL
     &, F_VECTOR(I_DIM-1)  !OUT. LAST COLUMN VECTOR
     &, G_VECTOR(2)        !OUT. LAST ROW VECTOR.
C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I,I2
     &, I_BEGIN
     &, I_END

      REAL
     &  SCALAR
     &, SCALARS
     &, SCALAR1
     &, SCALAR2

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. CALCULATE LEFT-HAND-SIDE OF P.D.E.
CL               AND MATRIX COEFFICIENTS.
CL----------------------------------------------------------------------

C ----------------------------------------------------------------------
CL    SECTION 1.1 CALCULATE I-DIRECTION DERIVATIVES.
C ----------------------------------------------------------------------

      SCALARS= EARTH_RADIUS_INVERSE*LONGITUDE_STEP_INVERSE
      SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &         LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE

C SET LAST POINT IN A0 WHICH IS USED, WHICH IS I_DIM-1 IN CYCLIC CODE.
      IF(VERSION.LT.3) THEN
        I_END = I_DIM-1
        DO I=1,I_DIM-1
          F_VECTOR(I) = 0.
        END DO
      ELSE
        I_END = I_DIM
      END IF

C SET ZERO COEFFICIENTS IN MATRICES.
      A1(I_END) = 0.
      A2(1) = 0.

      DO I= 2,I_END-1
        SCALAR1 = SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
        A1(I)  =  A(I,J,K) *SCALAR1 + .5*DEF(I,J,K)*D(I,J,K)*SCALARS
     &            *SEC_P_LATITUDE(I,J)
        A2(I)  =  A(I-1,J,K) *SCALAR1 - .5*DEF(I,J,K)*D(I-1,J,K)
     &            *SCALARS*SEC_P_LATITUDE(I,J)
        A0(I)  = - A1(I) - A2(I) - G(I,J,K)
        LHS(I) = (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &            (A(I-1,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &            *SCALAR1 + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &            (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &             (Q(I,J,K)-Q(I-1,J,K)))
     &            - Q(I,J,K)*G(I,J,K)
      END DO

      IF (VERSION.LT.3) THEN
CL GLOBAL VERSION.
CL    FIRST POINT.

        SCALAR1 = SCALAR*SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
        A1(1)  =  A(1,J,K) *SCALAR1 + .5*DEF(1,J,K)*D(1,J,K)*SCALARS
     &            *SEC_P_LATITUDE(1,J)
        F_VECTOR(1) = A(I_DIM,J,K)*SCALAR1
     &                - .5*DEF(1,J,K)*D(I_DIM,J,K)
     &                    *SCALARS*SEC_P_LATITUDE(1,J)
        A0(1)       = - A1(1) - F_VECTOR(1) - G(1,J,K)
        LHS(1)      =  (Q(I_DIM,J,K)*A(I_DIM,J,K) - Q(1,J,K) *
     &                 (A(I_DIM,J,K) + A(1,J,K)) + Q(2,J,K)*A(1,J,K))
     &            *SCALAR1 + .5*DEF(1,J,K)*SCALARS*SEC_P_LATITUDE(1,J)*
     &            (D(1,J,K)*(Q(2,J,K)-Q(1,J,K))+D(I_DIM,J,K)*
     &             (Q(1,J,K)-Q(I_DIM,J,K)))
     &            - Q(1,J,K)*G(1,J,K)

CL    PENULTIMATE POINT.

        I = I_DIM-1
        SCALAR1 = SCALAR*SEC_P_LATITUDE(I,J)*SEC_P_LATITUDE(I,J)
        F_VECTOR(I)= A(I,J,K) *SCALAR1 + .5*DEF(I,J,K)*D(I,J,K)*SCALARS
     &            *SEC_P_LATITUDE(I,J)
        A2(I)  =  A(I-1,J,K) *SCALAR1 - .5*DEF(I,J,K)*D(I-1,J,K)
     &            *SCALARS*SEC_P_LATITUDE(I,J)
        A0(I)       = - F_VECTOR(I) - A2(I) - G(I,J,K)
        LHS(I)      = (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &               (A(I-1,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &            *SCALAR1 + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &            (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &             (Q(I,J,K)-Q(I-1,J,K)))
     &            - Q(I,J,K)*G(I,J,K)

CL    LAST POINT.

        SCALAR1 = SCALAR*SEC_P_LATITUDE(I_DIM,J)*SEC_P_LATITUDE(I_DIM,J)
        G_VECTOR(1) = A(I_DIM,J,K) *SCALAR1 +
     &                .5*DEF(I_DIM,J,K)*D(I_DIM,J,K)*SCALARS
     &            *SEC_P_LATITUDE(I_DIM,J)
        G_VECTOR(2) = A(I_DIM-1,J,K)*SCALAR1
     &                - .5*DEF(I_DIM,J,K)*D(I_DIM-1,J,K)
     &            *SCALARS*SEC_P_LATITUDE(I_DIM,J)
        AN          = - G_VECTOR(1) - G_VECTOR(2) - G(I_DIM,J,K)
        LHS(I_DIM) =(Q(I_DIM-1,J,K)*A(I_DIM-1,J,K) - Q(I_DIM,J,K)
     &                 *(A(I_DIM-1,J,K) + A(I_DIM,J,K)) + Q(1,J,K)
     &                 *A(I_DIM,J,K))*SCALAR1
     &             + .5*DEF(I_DIM,J,K)*SCALARS*SEC_P_LATITUDE(I_DIM,J)*
     &            (D(I_DIM,J,K)*(Q(1,J,K)-Q(I_DIM,J,K))+D(I_DIM-1,J,K)*
     &             (Q(I_DIM,J,K)-Q(I_DIM-1,J,K)))
     &            - Q(I_DIM,J,K)*G(I_DIM,J,K)

      ELSE
CL LIMITED AREA VERSION.
CL    FIRST POINT.

        SCALAR1 = SCALAR*SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
        A1(1)       =  2.*A(1,J,K)*SCALAR1
        A0(1)       = - A1(1) - G(1,J,K)
        LHS(1)      =  2.*( - Q(1,J,K) *
     &                   A(1,J,K) + Q(2,J,K)*A(1,J,K))
     &                  *SCALAR1 - G(1,J,K) * Q(1,J,K)

CL    LAST POINT.

        SCALAR1 = SCALAR*SEC_P_LATITUDE(I_DIM,J)*SEC_P_LATITUDE(I_DIM,J)
        A2(I_DIM) = 2.*A(I_DIM-1,J,K)*SCALAR1
        A0(I_DIM)   = - A2(I_DIM) - G(I_DIM,J,K)
        LHS(I_DIM) = 2.*(Q(I_DIM-1,J,K)*A(I_DIM-1,J,K) - Q(I_DIM,J,K)
     &                 *A(I_DIM-1,J,K))*SCALAR1
     &                  - G(I_DIM,J,K)*Q(I_DIM,J,K)

      END IF

C ----------------------------------------------------------------------
CL    SECTION 1.2 CALCULATE J-DIRECTION DERIVATIVES AND ADD ON.
C ----------------------------------------------------------------------

      IF(J_NT) THEN
        SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &           LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
        SCALARS= EARTH_RADIUS_INVERSE*LATITUDE_STEP_INVERSE
        IF (J.EQ.1) THEN
C ONLY POSSIBLE IN LIMITED AREA CODE, AS NO I-LINE SOLVE CALLED AT POLE
C IN GLOBAL VERSION.
          DO I= 1,I_END
            LHS(I) = LHS(I) -2.*(Q(I,J,K)*B(I,J,K)
     &                      -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J)
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
            A0(I)  = A0(I) -2.*B(I,J,K)*COS_V_LATITUDE(I,J)
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
          END DO

        ELSE IF(J.EQ.J_DIM) THEN
C ONLY POSSIBLE IN LIMITED AREA CODE, AS NO I-LINE SOLVE CALLED AT POLE
C IN GLOBAL VERSION.
          DO I= 1,I_END
            LHS(I) = LHS(I) + 2.*(Q(I,J-1,K)*B(I,J-1,K)
     &                      -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
            A0(I)  = A0(I) - 2.*B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
          END DO

        ELSE
          DO I= 1,I_END
            LHS(I) = LHS(I) + ((Q(I,J-1,K)*B(I,J-1,K)
     &                      -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                      -(Q(I,J,K)*B(I,J,K)
     &                      -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &                      *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                       E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
            A0(I)  = A0(I) - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                   +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J,K)-E(I,J-1,K))
          END DO

          IF (VERSION.LT.3) THEN
CL    LAST POINT.

            I = I_DIM
            LHS(I) = LHS(I) +  ((Q(I,J-1,K)*B(I,J-1,K)
     &                    -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &                    -(Q(I,J,K)*B(I,J,K)
     &                    -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                       E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
            AN   = AN - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                   +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                    *SEC_P_LATITUDE(I,J)*SCALAR
     &                      +DEF(I,J,K)*.5*SCALARS*
     &                      (E(I,J,K)-E(I,J-1,K))
          END IF
        END IF
      END IF

C ----------------------------------------------------------------------
CL    SECTION 1.3 CALCULATE K-DIRECTION DERIVATIVES AND ADD ON.
C ----------------------------------------------------------------------

      IF (K_NT) THEN
        IF(K.EQ.1) THEN
          SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
          SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
          DO I=1,I_END
            LHS(I) = LHS(I) + SCALAR*(Q(I,J,K+1)*C2(I,J,K+1)
     &                   -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
            A0(I) = A0(I) - SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K)
            IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS(I) = LHS(I)  + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*
     &                        (Q(I,J,K+1)*C2(I,J,K+1)-
     &                         Q(I,J,K)*C2(I,J,K))
              A0(I) = A0(I) - .5*SCALAR1*DEF(I,J,K)*F(I,J,K)
     &                         *C2(I,J,K)
            END IF
          END DO
          IF (VERSION.LT.3) THEN
            I = I_DIM
            LHS(I) = LHS(I) + SCALAR*(Q(I,J,K+1)*C2(I,J,K+1)
     &                   -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
            AN = AN - SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K)
            IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS(I) = LHS(I) + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*
     &                        (Q(I,J,K+1)*C2(I,J,K+1)-
     &                         Q(I,J,K)*C2(I,J,K))
              AN = AN - .5*SCALAR1*DEF(I,J,K)*F(I,J,K)
     &                         *C2(I,J,K)
            END IF
          END IF
        ELSE IF (K.EQ.K_DIM) THEN
          SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
          SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
          DO I=1,I_END
            LHS(I) = LHS(I) - SCALAR*(Q(I,J,K)*C2(I,J,K)
     &                   -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
            A0(I) = A0(I) - SCALAR*C1(I,J,K-1)*SCALAR2*C2(I,J,K)
            IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS(I) = LHS(I)  + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*
     &                        (Q(I,J,K)*C2(I,J,K)-
     &                         Q(I,J,K-1)*C2(I,J,K-1))
              A0(I) = A0(I) + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)
     &                         *C2(I,J,K)
            END IF
          END DO
          IF (VERSION.LT.3) THEN
            I = I_DIM
            LHS(I) = LHS(I) - SCALAR*(Q(I,J,K)*C2(I,J,K)
     &                   -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
            AN = AN - SCALAR*C1(I,J,K-1)*SCALAR2*C2(I,J,K)
            IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
              LHS(I) = LHS(I)  + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*
     &                        (Q(I,J,K)*C2(I,J,K)-
     &                         Q(I,J,K-1)*C2(I,J,K-1))
              AN = AN + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)
     &                         *C2(I,J,K)
            END IF
          END IF
        ELSE
          SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
          SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
          SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
          DO I=1,I_END
            LHS(I) = LHS(I) + SCALAR*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *C1(I,J,K)*SCALAR1
     &               -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *C1(I,J,K-1)*SCALAR2)
     &               +.5*DEF(I,J,K)*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *F(I,J,K)*SCALAR1
     &               +(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *F(I,J,K-1)*SCALAR2)
            A0(I) = A0(I) - SCALAR*(C1(I,J,K)*SCALAR1
     &                              +C1(I,J,K-1)*SCALAR2)*C2(I,J,K)
     &               +.5*DEF(I,J,K)*(F(I,J,K-1)*SCALAR2
     &                               -F(I,J,K)*SCALAR1)*C2(I,J,K)
          END DO
          IF (VERSION.LT.3) THEN
            I = I_DIM
            LHS(I) = LHS(I) + SCALAR*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *C1(I,J,K)*SCALAR1
     &               -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *C1(I,J,K-1)*SCALAR2)
     &               +.5*DEF(I,J,K)*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *F(I,J,K)*SCALAR1
     &               +(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *F(I,J,K-1)*SCALAR2)
            AN = AN - SCALAR*(C1(I,J,K)*SCALAR1
     &                              +C1(I,J,K-1)*SCALAR2)*C2(I,J,K)
     &               +.5*DEF(I,J,K)*(F(I,J,K-1)*SCALAR2
     &                               -F(I,J,K)*SCALAR1)*C2(I,J,K)
          END IF
        END IF

      END IF

CL    END OF ROUTINE MG_I_MATRIX

      RETURN
      END
CLL   SUBROUTINE MG_J_LINE
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS LINE SMOOTHING IN J DIRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_J_LINE
     &                    (Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                     SEC_P_LATITUDE,
     &                     COS_V_LATITUDE,LATITUDE_STEP_INVERSE,
     &                     LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                     I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                     SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                     VERSION,K_BC,Z_Q,Z_MID,I_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, NFB         !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &              ! NEEDED.
     &, IPAT        !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &              ! POINTS.
     &, VERSION     !IN. VERSION OF MULTIGRID SCHEME.
     &, K_BC        !IN. BOUNDARY CONDITION IN K DIRECTION.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES

      REAL
     &  SWJAC    !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM    !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
     &, RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.
     &, W        !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM

      REAL
     &  A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, LONGITUDE_STEP_INVERSE !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   6 LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  LHS(J_DIM)
     &, A0(J_DIM,I_DIM)
     &, A1(J_DIM,I_DIM)
     &, A2(J_DIM,I_DIM)
     &, R(J_DIM,I_DIM)
     &, CORRECTIONS(I_DIM,J_DIM,K_DIM)

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I_BEGIN
     &, IK_BEGIN
     &, I_END
     &, I_STEP
     &, K_BEGIN
     &, K_END
     &, K_STEP
     &, I,J,K,J2
     &, J_FB
     &, I_ZEBRA
     &, J_EQ
     &, J_BEGIN(2) !TWO VALUES. 2ND ONLY USED IN VERSION 2
     &, J_END(2)   !  "           "              "
     &, SWEEPS

      REAL
     &  WEIGHT
     &, SCALAR
     &, FACTOR

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_J_MATRIX

C*----------------------------------------------------------------------

CFPP$ NOCONCUR R

CL----------------------------------------------------------------------
CL    SECTION 1. INITIALISATION.
CL----------------------------------------------------------------------

CL    INITIALISE NORMS

      RMS_RES=0.0
      RMS_INC=0.0

CL    SET TEMPORARY RELAXATION PARAMETER TO ZERO FOR JACOBI:

      WEIGHT = SWJAC*W

CL    DEFINE POINT ORDERING IN JK PLANE FOR FIRST SWEEP

      IF(VERSION.EQ.4 .AND. I_NT) THEN
        I_BEGIN = 2
        I_END   = I_DIM-1
      ELSE
        I_BEGIN = 1
        I_END   = I_DIM
      END IF
      IF (IPAT.EQ.2) THEN
        I_STEP  = 2
      ELSE
        I_STEP  = 1
      END IF
      K_BEGIN = 1
      K_END   = K_DIM
      IF (K_BC.EQ.2 ) THEN
        K_END = K_DIM-1
      ELSE IF (K_BC.EQ.3 ) THEN
        K_BEGIN = 2
      ELSE IF (K_BC.EQ. 4) THEN
        K_BEGIN = 2
        K_END = K_DIM-1
      END IF
      K_STEP  = 1
      J_EQ = (J_DIM+1)/2
      SWEEPS = 1
      IF (VERSION. EQ .3) THEN
CL LIMITED AREA VERSION. BOUNDARIES INCLUDED.
        J_BEGIN(1) = 1
        J_END(1) = J_DIM
      ELSE
CL BOUNDARIES NOT INCLUDED.
        J_BEGIN(1) = 2
        J_END(1) = J_DIM-1
        IF(VERSION .EQ.2) THEN
          J_END(1) = J_EQ - 1
          J_BEGIN(2) = J_EQ+1
          J_END(2) = J_DIM-1
          SWEEPS = 2
        END IF
      END IF

CL    SET ZEBRA SWITCH FOR FIRST SWEEP

      I_ZEBRA=0

CL----------------------------------------------------------------------
CL    SECTION 2. PERFORM APPROPRIATE LINE SMOOTHING.
CL----------------------------------------------------------------------

CL    FORWARD / BACKWARD SWEEP

      DO J_FB= 1,NFB

        RMS_RES = SWSYM*RMS_RES
        RMS_INC = SWSYM*RMS_INC

CL    SWEEP OVER IK PLANE IN PRESCRIBED POINT ORDER

        DO K = K_BEGIN,K_END,K_STEP

          IF (IPAT.EQ.2) THEN
CL    START EACH LINE ON A DIFFERENT COLOUR.
            IF(MOD(I_BEGIN+K+I_ZEBRA,IPAT).NE.0 .AND. J_FB.EQ.1) THEN
              IK_BEGIN = I_BEGIN + 1
            ELSE IF(MOD(I_BEGIN+K+I_ZEBRA,IPAT).NE.0 .AND.
     &              J_FB.EQ.2) THEN
              IK_BEGIN = I_BEGIN - 1
            ELSE
              IK_BEGIN = I_BEGIN
            END IF
          ELSE
            IK_BEGIN = I_BEGIN
          END IF

          DO I = IK_BEGIN,I_END,I_STEP

CL    ASSEMBLE TRIDIAGONAL SYSTEM TO INVERT FOR LINE VALUES

            CALL MG_J_MATRIX(LHS,A0(1,I),A1(1,I),A2(1,I),Q,A,B,
     &                       C1,C2,DEF,D,E,F,G,I_DIM,J_DIM,K_DIM,
     &                       COS_P_LATITUDE,SEC_P_LATITUDE,
     &                       COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                       LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,I,K,VERSION,
     &                       Z_Q,Z_MID,I_NT,K_NT)

            DO J2 = 1,SWEEPS
              DO J=J_BEGIN(J2),J_END(J2)
                R(J,I) = LHS(J)-RHS(I,J,K)
                RMS_RES = RMS_RES+R(J,I)*R(J,I)
              END DO
            END DO
          END DO

CL    SOLVE TRIDIAGONAL SYSTEM FOR CORRECTIONS X

C SOLUTION IS IN R(J)
C REDUCE MATRIX TO UPPER DIAGONAL FORM.

          DO J2 = 1,SWEEPS
            DO I=IK_BEGIN,I_END,I_STEP
              A0(J_BEGIN(J2),I) = 1./A0(J_BEGIN(J2),I)
            END DO
            DO J=J_BEGIN(J2)+1,J_END(J2)
              DO I=IK_BEGIN,I_END,I_STEP
                FACTOR = A2(J,I) * A0(J-1,I)
                A0(J,I) = 1./(A0(J,I) - FACTOR*A1(J-1,I))
                R(J,I) = R(J,I) - FACTOR*R(J-1,I)
              END DO
            END DO

C BACK SUBSTITUTE TO GET SOLUTION.

            DO I=IK_BEGIN,I_END,I_STEP
              R(J_END(J2),I)  = A0(J_END(J2),I)*R(J_END(J2),I)
            END DO
            DO J= J_END(J2)-1,J_BEGIN(J2),-1
              DO I=IK_BEGIN,I_END,I_STEP
                  R(J,I) = A0(J,I)*(R(J,I)-A1(J,I)*R(J+1,I))
              END DO
            END DO

CL    STORE CORRECTIONS IN WORKSPACE
CL    AND ADD TO CURRENT SOLUTION UNLESS WT=0 (JACOBI)

            DO J=J_BEGIN(J2),J_END(J2)
              DO I=IK_BEGIN,I_END,I_STEP
                Q(I,J,K)=Q(I,J,K)-WEIGHT*R(J,I)
                RMS_INC=RMS_INC+R(J,I)*R(J,I)
              END DO
            END DO

            IF (JAC) THEN
              DO J=J_BEGIN(J2),J_END(J2)
                DO I=IK_BEGIN,I_END,I_STEP
                  CORRECTIONS(I,J,K)=R(J,I)
                END DO
              END DO
            END IF

          END DO
CL    END OF IK PLANE SWEEP
        END DO

CL    ADD CORRECTIONS FOR THE JACOBI METHOD

        IF(JAC) THEN
          DO K=K_BEGIN,K_END
            DO J2 = 1,SWEEPS
              DO J=J_BEGIN(J2),J_END(J2)
                DO I=1,I_DIM
                  Q(I,J,K)=Q(I,J,K)-W*CORRECTIONS(I,J,K)
                END DO
              END DO
            END DO
          END DO
        ENDIF

CL    RESET POINT ORDERING OF JK PLANE FOR SECOND SWEEP

        I = I_BEGIN
        I_BEGIN = I_END
        I_END = I
        I_STEP  = -I_STEP
        K = K_BEGIN
        K_BEGIN = K_END
        K_END   =  K
        K_STEP  = -1

CL    RESET ZEBRA SWITCH FOR SECOND SWEEP

        I_ZEBRA=1

      END DO

CL    END OF ROUTINE MG_J_LINE

      RETURN
      END
CLL   SUBROUTINE MG_J_MATRIX
CLL
CLL   PURPOSE:
CLL   -------
CLL   CALCULATES CURRENT VALUE OF LEFT-HAND-SIDE OF P.D.E. FOR A GIVEN
CLL   VALUE OF I AND K. SETS MATRIX ENTRIES TO SOLVE FOR CORRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3             WRITTEN BY M. H. MAWSON.
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_J_MATRIX(
     &                       LHS,A0,A1,A2,Q,A,B,C1,C2,DEF,D,E,F,G,
     &                       I_DIM,J_DIM,K_DIM,
     &                       COS_P_LATITUDE,SEC_P_LATITUDE,
     &                       COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                       LATITUDE_STEP_INVERSE,
     &                       LONGITUDE_STEP_INVERSE,
     &                       I,K,VERSION,Z_Q,Z_MID,I_NT,K_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, K_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN K DIRECTION

      INTEGER
     &  I_DIM         !IN. NUMBER OF NODES IN THE I-DIRECTION
     &, J_DIM         !IN. NUMBER OF NODES IN THE J-DIRECTION
     &, K_DIM         !IN. NUMBER OF NODES IN THE K-DIRECTION
     &, I             !IN. VALUE OF I FOR WHICH LHS AND COEFFS ARE REQ.
     &, K             !IN. VALUE OF K FOR WHICH LHS AND COEFFS ARE REQ.
     &, VERSION       !IN. VERSION OF MULTIGRID SCHEME.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !IN. SOLUTION
     &, A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)!IN. COSINE OF LATITUDE AT B POINTS
     &, SEC_P_LATITUDE(I_DIM,J_DIM)!IN. 1./COSINE OF LATITUDE AT Q POINT
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  EARTH_RADIUS_INVERSE    !IN. 1/EARTH RADIUS
     &, LATITUDE_STEP_INVERSE   !IN. 1/LATITUDE STEP LENGTH IN RADIANS
     &, LONGITUDE_STEP_INVERSE  !IN. 1/LONGITUDE STEP LENGTH IN RADIANS

      REAL
     &  LHS(J_DIM)       !OUT. LEFT-HAND-SIDE.
     &, A0(J_DIM)        !OUT. DIAGONAL MATRIX COEFFICIENT.
     &, A1(J_DIM)        !OUT. SUPER-DIAGONAL MATRIX COEFFICIENT.
     &, A2(J_DIM)        !OUT. SUB-DIAGONAL MATRIX COEFFICIENT.
C*----------------------------------------------------------------------

C*L NO LOCAL WORK ARRAYS REQUIRED.
C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  J,J2
     &, J_BEGIN(2) !2ND VALUE USED ONLY IN VERSION 2
     &, J_END(2)   !  "        "           "
     &, SWEEPS

      REAL
     &  SCALAR
     &, SCALARS
     &, SCALAR1
     &, SCALAR2

C*L NO EXTERNAL ROUTINES CALLED.
C*----------------------------------------------------------------------

CL----------------------------------------------------------------------
CL    SECTION 1. CALCULATE LEFT-HAND-SIDE OF P.D.E.
CL               AND MATRIX COEFFICIENTS.
CL----------------------------------------------------------------------

      SWEEPS = 1
      IF (VERSION. LT.3) THEN
CL GLOBAL VERSION. POLES NOT INCLUDED.
        J_BEGIN(1) = 2
        J_END(1) = J_DIM-1
        IF(VERSION .EQ.2) THEN
          J_END(1) = (J_DIM+1)/2 - 1
          J_BEGIN(2) = (J_DIM+1)/2+1
          J_END(2) = J_DIM-1
          SWEEPS = 2
        END IF
      ELSE
CL LIMITED AREA VERSION. BOUNDARIES INCLUDED.
        J_BEGIN(1) = 1
        J_END(1) = J_DIM
      END IF

C ----------------------------------------------------------------------
CL    SECTION 1.1 CALCULATE J-DIRECTION DERIVATIVES AND MATRIX COEFFS.
C ----------------------------------------------------------------------

      DO J2 =1,SWEEPS

        A1(J_END(J2)) = 0.
        A2(J_BEGIN(J2)) = 0.
        SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &           LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
        SCALARS=EARTH_RADIUS_INVERSE*LATITUDE_STEP_INVERSE

CL    ALL POINTS EXCEPT ENDS.

        DO J= J_BEGIN(J2)+1,J_END(J2)-1
          SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
          LHS(J) = ((Q(I,J-1,K)*B(I,J-1,K)
     &               -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &               -(Q(I,J,K)*B(I,J,K)
     &               -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &                *SCALAR1 - Q(I,J,K)*G(I,J,K)
     &                 +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                  E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
          A0(J)  = - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &              +B(I,J,K)*COS_V_LATITUDE(I,J))
     &              *SCALAR1 - G(I,J,K)
     &                +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J,K)-E(I,J-1,K))
          A1(J)  = B(I,J,K)*COS_V_LATITUDE(I,J)*SCALAR1
     &             -DEF(I,J,K)*.5*SCALARS*E(I,J,K)
          A2(J)  =  B(I,J-1,K)*COS_V_LATITUDE(I,J-1)*SCALAR1
     &             +DEF(I,J,K)*.5*SCALARS*E(I,J-1,K)
        END DO

CL    END POINTS.

        IF(VERSION.NE.3) THEN
CL GLOBAL
          J= J_BEGIN(J2)
          SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
          LHS(J) = ((Q(I,J-1,K)*B(I,J-1,K)
     &             -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &             -(Q(I,J,K)*B(I,J,K)
     &             -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &             *SCALAR1 - Q(I,J,K)*G(I,J,K)
     &                 +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                  E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
          A0(J)  = - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                *SCALAR1 - G(I,J,K)
     &                +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J,K)-E(I,J-1,K))
          A1(J)  = B(I,J,K)*COS_V_LATITUDE(I,J) *SCALAR1
     &             -DEF(I,J,K)*.5*SCALARS*E(I,J,K)

          J= J_END(J2)
          SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
          LHS(J) = ((Q(I,J-1,K)*B(I,J-1,K)
     &               -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &               -(Q(I,J,K)*B(I,J,K)
     &               -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J))
     &              *SCALAR1 - Q(I,J,K)*G(I,J,K)
     &                 +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J-1,K)*(Q(I,J-1,K)-Q(I,J,K))+
     &                  E(I,J,K)*(Q(I,J,K)-Q(I,J+1,K)))
          A0(J)  = - (B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                +B(I,J,K)*COS_V_LATITUDE(I,J))
     &                *SCALAR1 - G(I,J,K)
     &                +DEF(I,J,K)*.5*SCALARS*
     &                 (E(I,J,K)-E(I,J-1,K))
          A2(J)  =  B(I,J-1,K)*COS_V_LATITUDE(I,J-1) *SCALAR1
     &             +DEF(I,J,K)*.5*SCALARS*E(I,J-1,K)

        ELSE IF(VERSION.EQ.3) THEN
CL LIMITED AREA
          J= J_BEGIN(J2)
          SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
          LHS(J) = -2.*(Q(I,J,K)*B(I,J,K)
     &                -Q(I,J+1,K)*B(I,J,K))*COS_V_LATITUDE(I,J)
     &             *SCALAR1 - Q(I,J,K)*G(I,J,K)
          A0(J)  = -2.*B(I,J,K)*COS_V_LATITUDE(I,J)
     &                *SCALAR1 - G(I,J,K)
          A1(J)  = 2.*B(I,J,K)*COS_V_LATITUDE(I,J) *SCALAR1

          J= J_END(J2)
          SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
          LHS(J) = 2.*(Q(I,J-1,K)*B(I,J-1,K)
     &               -Q(I,J,K)*B(I,J-1,K))*COS_V_LATITUDE(I,J-1)
     &              *SCALAR1 - Q(I,J,K)*G(I,J,K)
          A0(J)  = - 2.*B(I,J-1,K)*COS_V_LATITUDE(I,J-1)
     &                *SCALAR1 - G(I,J,K)
          A2(J)  =  2.*B(I,J-1,K)*COS_V_LATITUDE(I,J-1) *SCALAR1

        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.2 CALCULATE I-DIRECTION DERIVATIVES
CL                AND ADD ON TO J-DIRECTION TERMS.
C ----------------------------------------------------------------------

        IF(I_NT) THEN

          SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &             LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE
          SCALARS = EARTH_RADIUS_INVERSE*LONGITUDE_STEP_INVERSE

          IF(I.EQ.1) THEN

CL    FIRST POINT ON A ROW.

            IF(VERSION.LT.3) THEN
CL GLOBAL
              DO J=J_BEGIN(J2),J_END(J2)
                SCALAR1= SCALAR*SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
                A0(J) = A0(J) -(A(1,J,K)+ A(I_DIM,J,K))*SCALAR1
     &                 +.5*DEF(1,J,K)*SCALARS*SEC_P_LATITUDE(1,J)*
     &                 (D(I_DIM,J,K)-D(1,J,K))
                LHS(J) =  LHS(J) +(Q(I_DIM,J,K)*A(I_DIM,J,K)- Q(1,J,K) *
     &                   (A(I_DIM,J,K) + A(1,J,K)) + Q(2,J,K)*A(1,J,K))
     &                    *SCALAR1
     &                   + .5*DEF(1,J,K)*SCALARS*SEC_P_LATITUDE(1,J)*
     &                   (D(1,J,K)*(Q(2,J,K)-Q(1,J,K))+D(I_DIM,J,K)*
     &                   (Q(1,J,K)-Q(I_DIM,J,K)))
              END DO
            ELSE
CL LIMITED AREA
              DO J=J_BEGIN(J2),J_END(J2)
                SCALAR1= SCALAR*SEC_P_LATITUDE(1,J)*SEC_P_LATITUDE(1,J)
                A0(J) = A0(J) - 2.*A(1,J,K)*SCALAR1
                LHS(J)    =  LHS(J) - 2.*(Q(1,J,K) * A(1,J,K)
     &                                 - Q(2,J,K)*A(1,J,K))*SCALAR1
              END DO
            END IF

          ELSE IF(I.EQ.I_DIM) THEN

CL    LAST POINT ON A ROW.

            IF(VERSION.LT.3) THEN
CL GLOBAL
              DO J=J_BEGIN(J2),J_END(J2)
                SCALAR1 = SCALAR * SEC_P_LATITUDE(I_DIM,J)
     &                           * SEC_P_LATITUDE(I_DIM,J)
                A0(J) = A0(J) - (A(I_DIM,J,K)+
     &                          A(I_DIM-1,J,K))*SCALAR1
     &              +.5*DEF(I_DIM,J,K)*SCALARS*SEC_P_LATITUDE(I_DIM,J)*
     &                 (D(1,J,K)-D(I_DIM,J,K))
                LHS(J) = LHS(J) +
     &                    (Q(I_DIM-1,J,K)*A(I_DIM-1,J,K) - Q(I_DIM,J,K)
     &                    *(A(I_DIM-1,J,K) + A(I_DIM,J,K)) + Q(1,J,K)
     &                    *A(I_DIM,J,K))*SCALAR1
     &                   + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &            (D(I_DIM,J,K)*(Q(1,J,K)-Q(I_DIM,J,K))+D(I_DIM-1,J,K)*
     &             (Q(I_DIM,J,K)-Q(I_DIM-1,J,K)))
              END DO
            ELSE
CL LIMITED AREA
              DO J=J_BEGIN(J2),J_END(J2)
                SCALAR1 = SCALAR * SEC_P_LATITUDE(I_DIM,J)
     &                           * SEC_P_LATITUDE(I_DIM,J)
                A0(J) = A0(J) -  2.*A(I_DIM-1,J,K)*SCALAR1
                LHS(J) = LHS(J) + 2.*
     &                    (Q(I_DIM-1,J,K)*A(I_DIM-1,J,K) - Q(I_DIM,J,K)
     &                    *A(I_DIM-1,J,K))*SCALAR1
              END DO
            END IF

          ELSE
CL ALL OTHER POINTS.

            DO J=J_BEGIN(J2),J_END(J2)
              SCALAR1 = SCALAR * SEC_P_LATITUDE(I,J)
     &                         * SEC_P_LATITUDE(I,J)
              A0(J) = A0(J) - ( A(I,J,K) + A(I-1,J,K))*SCALAR1
     &                 +.5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &                 (D(I-1,J,K)-D(I,J,K))
              LHS(J)    = LHS(J) + (Q(I-1,J,K)*A(I-1,J,K) - Q(I,J,K) *
     &               (A(I-1,J,K) + A(I,J,K)) + Q(I+1,J,K)*A(I,J,K))
     &                *SCALAR1
     &                + .5*DEF(I,J,K)*SCALARS*SEC_P_LATITUDE(I,J)*
     &               (D(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))+D(I-1,J,K)*
     &               (Q(I,J,K)-Q(I-1,J,K)))
            END DO

          END IF

        END IF

C ----------------------------------------------------------------------
CL    SECTION 1.3 CALCULATE K-DIRECTION DERIVATIVES AND ADD ON.
C ----------------------------------------------------------------------

        IF (K_NT) THEN
          IF(K.EQ.1) THEN
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            DO J=J_BEGIN(J2),J_END(J2)
              LHS(J) = LHS(J) + SCALAR*(Q(I,J,K+1)*C2(I,J,K+1)
     &                     -Q(I,J,K)*C2(I,J,K))*C1(I,J,K)*SCALAR1
              A0(J) = A0(J) - SCALAR*C1(I,J,K)*SCALAR1*C2(I,J,K)
              IF(Z_Q(1).NE.Z_MID(1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                LHS(J) = LHS(J)  + .5*SCALAR1*DEF(I,J,K)*F(I,J,K)*
     &                               (Q(I,J,K+1)*C2(I,J,K+1)
     &                                -Q(I,J,K)*C2(I,J,K))
                A0(J) = A0(J) - .5*SCALAR1*DEF(I,J,K)*F(I,J,K)
     &                          *C2(I,J,K)
              END IF
            END DO
          ELSE IF (K.EQ.K_DIM) THEN
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            DO J=J_BEGIN(J2),J_END(J2)
              LHS(J) = LHS(J) - SCALAR*(Q(I,J,K)*C2(I,J,K)
     &                   -Q(I,J,K-1)*C2(I,J,K-1))*C1(I,J,K-1)*SCALAR2
              A0(J) = A0(J) - SCALAR*C1(I,J,K-1)*SCALAR2*C2(I,J,K)
              IF(Z_Q(K_DIM).NE.Z_MID(K_DIM+1)) THEN
C INCLUDE FIRST DERIVATIVE TERM WHEN NON-ZERO
                LHS(J) = LHS(J)  + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)*
     &                               (Q(I,J,K)*C2(I,J,K)
     &                                -Q(I,J,K-1)*C2(I,J,K-1))
                A0(J) = A0(J) + .5*SCALAR2*DEF(I,J,K)*F(I,J,K-1)
     &                          *C2(I,J,K)
              END IF
            END DO
          ELSE
            SCALAR = 1./(Z_MID(K) - Z_MID(K-1))
            SCALAR1 = 1./(Z_Q(K+1) - Z_Q(K))
            SCALAR2 = 1./(Z_Q(K) - Z_Q(K-1))
            DO J=J_BEGIN(J2),J_END(J2)
              LHS(J) = LHS(J) + SCALAR*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *C1(I,J,K)*SCALAR1
     &               -(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *C1(I,J,K-1)*SCALAR2)
     &               +.5*DEF(I,J,K)*
     &               ((Q(I,J,K+1)*C2(I,J,K+1)-Q(I,J,K)*C2(I,J,K))
     &                *F(I,J,K)*SCALAR1
     &               +(Q(I,J,K)*C2(I,J,K)-Q(I,J,K-1)*C2(I,J,K-1))
     &                *F(I,J,K-1)*SCALAR2)
              A0(J) = A0(J) - SCALAR*(C1(I,J,K)*SCALAR1
     &                              +C1(I,J,K-1)*SCALAR2)*C2(I,J,K)
     &               +.5*DEF(I,J,K)*(F(I,J,K-1)*SCALAR2
     &                               -F(I,J,K)*SCALAR1)*C2(I,J,K)
            END DO
          END IF

        END IF

CL END LOOP OVER SWEEPS
      END DO

CL    END OF ROUTINE MG_J_MATRIX

      RETURN
      END
CLL   SUBROUTINE MG_JOB_INFO
CLL
CLL   PURPOSE:
CLL   -------
CLL   PROVIDES INFORMATION ON CHOICES OF SCHEME MADE.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_JOB_INFO(MAXITS,TOL_RES,KSMOOTH,
     1                       NPRE,NPOST,NCOARSE,W,KREST,NCGC,NGRIDS,
     2                       IMAX,JMAX,KMAX,START_ADDRESS,RES_DIRS,
     3                       WORST_SMOOTHING_RATE,VERSION,K_BC)

      IMPLICIT NONE

      INTEGER
     &  VERSION    !IN. VERSION OF SCHEME BEING USED.
     &, K_BC       !IN. BOUNDARY CONDITIONS TO USE IN K DIRECTION.

      INTEGER
     &  NGRIDS     !IN. NUMBER OF GRIDS.
     &, MAXITS     !IN. MAX NO OF FAS ITERATIONS WITHOUT CONVERGENCE
     &, KSMOOTH    !IN. KIND OF ITERATIVE METHOD USED AS A SMOOTHER.
     &, NPRE       !IN. NO OF PRE-SMOOTHING SWEEPS
     &, NPOST      !IN. NO OF POST-SMOOTHING SWEEPS
     &, NCOARSE    !IN. NO OF ITERATIONS OF SMOOTHER ON COARSEST MESH
     &, KREST      !IN. KIND OF RESTRICTION USED.
     &, NCGC       !IN. NO OF COARSE GRID CORRECTIONS.

      INTEGER
     &  IMAX(NGRIDS) !IN. NUMBER OF NODES IN THE I-DIRECTION
     &               !    ON EACH GRID.
     &, JMAX(NGRIDS) !IN. NUMBER OF NODES IN THE J-DIRECTION
     &               !    ON EACH GRID.
     &, KMAX(NGRIDS) !IN. NUMBER OF NODES IN THE K-DIRECTION
     &               !    ON EACH GRID.
     &, START_ADDRESS(NGRIDS) !IN. START ADDRESS IN DATA ARRAY FOR
     &                        !    EACH GRID.
     &, RES_DIRS(NGRIDS)      !IN. RESTRICTED DIRECTIONS.

      REAL
     &  W          !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM
     &, WORST_SMOOTHING_RATE !IN WORST PRACTICAL SMOOTHING RATE ALLOWED
     &, TOL_RES    !IN. TOLERANCE FOR RESIDUAL NORM
C*----------------------------------------------------------------------

C*L  2 LOCAL WORK ARRAYS REQUIRED.

      CHARACTER*12 SMOOTHER(26)
      CHARACTER*14 RESTN(3)

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I

CL----------------------------------------------------------------------
CL    SECTION 1. SET DATA VALUES.
CL----------------------------------------------------------------------

      DATA SMOOTHER /'    I-JACOBI','    J-JACOBI','    K-JACOBI',
     &               '  I&J-JACOBI','  I&K-JACOBI','  J&K-JACOBI',
     &               '    A-JACOBI','    I-L.G.S.','    J-L.G.S.',
     &               '    K-L.G.S.','  I&J-L.G.S.','  I&K-L.G.S.',
     &               '  J&K-L.G.S.','    A-L.G.S.','  A.S.L.G.S.',
     &               '     I-ZEBRA','     J-ZEBRA','     K-ZEBRA',
     &               '   I&J-ZEBRA','   I&K-ZEBRA','   J&K-ZEBRA',
     &               '     A-ZEBRA','   JACOBI-PT','     G.S.-PT',
     &               ' SYM G.S.-PT','  RED/BLK-PT'/

      DATA RESTN /'Injection','Half-Injection','Full-Weighting'/

CL----------------------------------------------------------------------
CL    SECTION 2. OUTPUT INFORMATION ON CHOICES.
CL----------------------------------------------------------------------

CL    VERSION INFORMATION.

      IF (VERSION.EQ.1) THEN
        WRITE(6,9000)
      ELSE IF (VERSION.EQ.2) THEN
        WRITE(6,8998)
        WRITE(6,8999)
      ELSE IF (VERSION.EQ.3) THEN
        WRITE(6,8997)
        WRITE(6,8996)
      ELSE IF (VERSION.EQ.4) THEN
        WRITE(6,8997)
        WRITE(6,8995)
      END IF

      IF (K_BC .EQ. 1) THEN
        WRITE(6,8980)
      ELSE IF (K_BC .EQ. 2) THEN
        WRITE(6,8981)
      ELSE IF (K_BC .EQ. 3) THEN
        WRITE(6,8982)
      ELSE IF (K_BC .EQ. 4) THEN
        WRITE(6,8983)
      END IF
CL    NAME OF CHOSEN SMOOTHING METHOD

      WRITE(6,9001) SMOOTHER(KSMOOTH)
      WRITE(6,9002)

CL    OUTPUT DETAILS OF THE SEQUENCE OF GRIDS TO BE USED

      DO I=1,NGRIDS
        WRITE(6,9003) I,IMAX(I),JMAX(I),KMAX(I),START_ADDRESS(I),
     &                RES_DIRS(I)
      END DO

C     PARAMETERS OF THE MULTIGRID ALGORITHM

      WRITE(6,9004) MAXITS,NPRE,NPOST,NCOARSE,NCGC,RESTN(KREST)

CL    RELAXATION PARAMETERS FOR EACH VARIABLE

      WRITE(6,9006)
      WRITE(6,9007) W

CL    TOLERANCES TO BE ACHIEVED

      WRITE(6,9008) TOL_RES
      WRITE(6,9009) WORST_SMOOTHING_RATE

 8980 FORMAT(1X,/,' Using Neumann vertical boundary conditions.')
 8981 FORMAT(1X,/,' Using Neumann bottom and Dirichlet top'
     &           ,' boundary conditions.')
 8982 FORMAT(1X,/,' Using Dirichlet bottom and Neumann top'
     &           ,' boundary conditions.')
 8983 FORMAT(1X,/,' Using Dirichlet vertical boundary conditions.')
 8995 FORMAT(1X,/,' Using Dirichlet lateral boundary conditions.')
 8996 FORMAT(1X,/,' Using Neumann lateral boundary conditions.')
 8997 FORMAT(1X,/,' Limited Area version of Multi-Grid requested.')
 8998 FORMAT(1X,/,' Global version of Multi-Grid requested with given'
     &           ,' equatorial solution.')
 8999 FORMAT(1X,/,' NB: Equatorial solution MUST have been in solution '
     &            ,'array when call to MGCNTL was made.')
 9000 FORMAT(1X,/,' Global version of Multi-Grid requested.')
 9001 FORMAT(1X,/,' Chosen smoother is ',A12,/,32('~'))
 9002 FORMAT(1X,/,' Grid Structure',/,15('-'),//,
     &       '  GRID INDEX   IMAX    JMAX    KMAX     IST     ',
     &       'RES DIRS',/,31('-'))
 9003 FORMAT(7X,I2,2X,5I8)
 9004 FORMAT(1X,/,' MG Strategy',/,12('-'),/,
     &            ' No Of FAS Iterations =',I6,/,
     &            ' No Of Pre-Smoothing Iterations =',I2,/,
     &            ' No Of Post-Smoothing Iterations =',I2,/,
     &            ' No Of Coarsest Grid Iterations =',I2,/,
     &            ' No Of Coarse Grid Corrections =',I2,//,
     &            ' Restriction Is By ',A14)
 9006 FORMAT(1X,/,' RELAXATION PARAMETER',/,37('-'))
 9007 FORMAT(5X,F4.2)
 9008 FORMAT(1X,/,' Tolerances : Residual Norm   =',E12.6,' times ',
     &            'initial Residual.')
 9009 FORMAT(1x,/,' Worst Practical smoothing rate allowed = ',F10.6)

      RETURN
      END
CLL   SUBROUTINE MG_K_LINE
CLL
CLL   PURPOSE:
CLL   -------
CLL   PERFORMS LINE SMOOTHING IN K DIRECTION.
CLL
CLL   VERSION NUMBER 2.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.3
CLL
CLL
CLLEND -----------------------------------------------------------------

C*L ARGUMENT LIST.

      SUBROUTINE MG_K_LINE
     &                    (Q,RHS,A,B,C1,C2,DEF,D,E,F,G,COS_P_LATITUDE,
     &                     SEC_P_LATITUDE,
     &                     COS_V_LATITUDE,LATITUDE_STEP_INVERSE,
     &                     LONGITUDE_STEP_INVERSE,EARTH_RADIUS_INVERSE,
     &                     I_DIM,J_DIM,K_DIM,W,RMS_RES,RMS_INC,
     &                     SWJAC,SWSYM,JAC,PAT,SYM,NFB,IPAT,
     &                     VERSION,K_BC,Z_Q,Z_MID,I_NT,J_NT)

      IMPLICIT NONE

      LOGICAL
     &  I_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN I DIRECTION
     &, J_NT    !IN. TRUE IF NON-TRIVIAL NUMBER OF POINTS IN J DIRECTION

      INTEGER
     &  I_DIM       !IN. NUMBER OF POINTS IN I-DIRECTION.
     &, J_DIM       !IN. NUMBER OF POINTS IN J-DIRECTION.
     &, K_DIM       !IN. NUMBER OF POINTS IN K-DIRECTION.
     &, NFB         !IN. NUMBER OF FORWARD / BACKWARD SWEEPS OF GRID
     &              ! NEEDED.
     &, IPAT        !IN. USED TO DISTINGUISH BETWEEN RED AND BLACK
     &              ! POINTS.
     &, NDIM        !IN. NUMBER OF DIMENSIONS.
     &, VERSION     !IN. VERSION OF MULTIGRID SCHEME.
     &, K_BC        !IN. BOUNDARY CONDITION IN K DIRECTION.

      LOGICAL
     &  JAC     !IN. .TRUE. FOR JACOBI METHODS
     &, PAT     !IN. .TRUE. FOR PATTERN SCHEMES
     &, SYM     !IN. .TRUE. FOR SYMMETRIC SCHEMES

      REAL
     &  SWJAC    !IN. A SWITCH WHICH IS ZERO FOR JACOBI METHODS
     &, SWSYM    !IN. A SWITCH WHICH IS ZERO FOR SYMMETRIC METHODS
     &, RMS_RES  !IN. ROOT MEAN SQUARE RESIDUAL NORM.
     &, RMS_INC  !IN. ROOT MEAN SQUARE INCREMENT NORM.
     &, W        !IN. RELAXATION PARAMETER FOR EACH VARIABLE IN SYSTEM

      REAL
     &  A(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, B(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, C1(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, C2(I_DIM,J_DIM,K_DIM)   !IN. COEFFICIENT
     &, DEF(I_DIM,J_DIM,K_DIM)  !IN. COEFFICIENT
     &, D(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, E(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, F(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, G(I_DIM,J_DIM,K_DIM)    !IN. COEFFICIENT
     &, RHS(I_DIM,J_DIM,K_DIM)  !IN. RIGHT-HAND-SIDE OF EQUATION.
     &, COS_P_LATITUDE(I_DIM,J_DIM) !IN. COSINE OF LATITUDE AT Q POINTS.
     &, SEC_P_LATITUDE(I_DIM,J_DIM) !IN. 1./COS OF LATITUDE AT Q POINTS.
     &, COS_V_LATITUDE(I_DIM,J_DIM-1)
     &                          !IN. COSINE OF LATITUDE AT V POINTS.
     &, Z_Q(K_DIM)              !IN. Z AT Q POINTS.
     &, Z_MID(0:K_DIM)          !IN. Z MIDWAY BETWEEN Q POINTS.

      REAL
     &  Q(I_DIM,J_DIM,K_DIM)    !OUT. SOLUTION.

      REAL
     &  LATITUDE_STEP_INVERSE  !IN.
     &, LONGITUDE_STEP_INVERSE !IN.
     &, EARTH_RADIUS_INVERSE   !IN

C*----------------------------------------------------------------------

C*L   6 LOCAL WORK ARRAYS REQUIRED.

      REAL
     &  LHS(K_DIM)
     &, A0(K_DIM,I_DIM)
     &, A1(K_DIM,I_DIM)
     &, A2(K_DIM,I_DIM)
     &, R(K_DIM,I_DIM)
     &, CORRECTIONS(I_DIM,J_DIM,K_DIM)

C*----------------------------------------------------------------------
C LOCAL VARIABLES.

      INTEGER
     &  I_BEGIN
     &, IK_BEGIN
     &, I_END
     &, I_STEP
     &, J_BEGIN
     &, J_END
     &, J_STEP
     &, K_BEGIN
     &, K_END
     &, I,J,K
     &, K_FB
     &, I_ZEBRA

      REAL
     &  WEIGHT
     &, SCALAR
     &, FACTOR

C*L   EXTERNAL ROUTINES CALLED.
      EXTERNAL MG_K_MATRIX

C*----------------------------------------------------------------------

CFPP$ NOCONCUR R

CL----------------------------------------------------------------------
CL    SECTION 1. INITIALISATION.
CL----------------------------------------------------------------------

CL    INITIALISE NORMS

      RMS_RES=0.0
      RMS_INC=0.0

CL    SET TEMPORARY RELAXATION PARAMETER TO ZERO FOR JACOBI:

      WEIGHT = SWJAC*W

CL    DEFINE POINT ORDERING IN IJ PLANE FOR FIRST SWEEP

      I_BEGIN = 1
      I_END   = I_DIM
      IF( IPAT.EQ.2) THEN
        I_STEP  = 2
      ELSE
        I_STEP  = 1
      END IF
      J_BEGIN = 1
      J_END   = J_DIM
      J_STEP  = 1
      IF(VERSION.EQ.4) THEN
CL DO NOT DO BOUNDARIES.
        IF(I_NT) THEN
          I_BEGIN = 2
          I_END = I_DIM-1
        END IF
        IF(J_NT) THEN
          J_BEGIN = 2
          J_END = J_DIM-1
        END IF
      END IF
      K_BEGIN = 1
      K_END = K_DIM
      IF (K_BC.EQ.2 ) THEN
        K_END = K_DIM-1
      ELSE IF (K_BC.EQ.3 ) THEN
        K_BEGIN = 2
      ELSE IF (K_BC.EQ. 4) THEN
        K_BEGIN = 2
        K_END = K_DIM-1
      END IF

CL    SET ZEBRA SWITCH FOR FIRST SWEEP

      I_ZEBRA=0

CL----------------------------------------------------------------------
CL    SECTION 2. PERFORM APPROPRIATE LINE SMOOTHING.
CL----------------------------------------------------------------------

CL    FORWARD / BACKWARD SWEEP

      DO K_FB= 1,NFB

        RMS_RES = SWSYM*RMS_RES
        RMS_INC = SWSYM*RMS_INC

CL    SWEEP OVER IK PLANE IN PRESCRIBED POINT ORDER

        DO J = J_BEGIN,J_END,J_STEP
          IF(.NOT.(VERSION.EQ.2.AND.J.EQ.(J_DIM+1)/2)) THEN

            IF (IPAT.EQ.2) THEN
CL    START EACH LINE ON A DIFFERENT COLOUR.
              IF(MOD(I_BEGIN+J+I_ZEBRA,IPAT).NE.0 .AND. K_FB.EQ.1) THEN
                IK_BEGIN = I_BEGIN + 1
              ELSE IF(MOD(I_BEGIN+J+I_ZEBRA,IPAT).NE.0 .AND.
     &                K_FB.EQ.2) THEN
                IK_BEGIN = I_BEGIN - 1
              ELSE
                IK_BEGIN = I_BEGIN
              END IF
            ELSE
              IK_BEGIN = I_BEGIN
            END IF

            DO I = IK_BEGIN,I_END,I_STEP

CL    ASSEMBLE TRIDIAGONAL SYSTEM TO INVERT FOR LINE VALUES

              CALL MG_K_MATRIX(LHS,A0(1,I),A1(1,I),A2(1,I),Q,A,B,
     &                         C1,C2,DEF,D,E,F,G,I_DIM,J_DIM,K_DIM,
     &                         COS_P_LATITUDE,SEC_P_LATITUDE,
     &                         COS_V_LATITUDE,EARTH_RADIUS_INVERSE,
     &                         LATITUDE_STEP_INVERSE,
     &                         LONGITUDE_STEP_INVERSE,I,J,VERSION,
     &                         Z_Q,Z_MID,I_NT,J_NT)

              DO K=K_BEGIN,K_END
                R(K,I) = LHS(K)-RHS(I,J,K)
                RMS_RES = RMS_RES+R(K,I)*R(K,I)
              END DO
            END DO

CL    SOLVE TRIDIAGONAL SYSTEM FOR CORRECTIONS X

C SOLUTION IS IN R(K)
C REDUCE MATRIX TO UPPER DIAGONAL FORM.

            DO I=IK_BEGIN,I_END,I_STEP
              A0(K_BEGIN,I) = 1./A0(K_BEGIN,I)
            END DO
            DO K=K_BEGIN+1,K_END
              DO I=IK_BEGIN,I_END,I_STEP
                FACTOR = A2(K,I) * A0(K-1,I)
                A0(K,I) = 1./(A0(K,I) - FACTOR*A1(K-1,I))
                R(K,I) = R(K,I) - FACTOR*R(K-1,I)
              END DO
            END DO

C BACK SUBSTITUTE TO GET SOLUTION.

            DO I=IK_BEGIN,I_END,I_STEP
              R(K_END,I)  = A0(K_END,I)*R(K_END,I)
            END DO
            DO K= K_END-1,K_BEGIN,-1
              DO I=IK_BEGIN,I_END,I_STEP
                R(K,I) = A0(K,I)*(R(K,I)-A1(K,I)*R(K+1,I))
              END DO
            END DO

CL    STORE CORRECTIONS IN WORKSPACE
CL    AND ADD TO CURRENT SOLUTION UNLESS WT=0 (JACOBI)

            DO K=K_BEGIN,K_END
              DO I=IK_BEGIN,I_END,I_STEP
                Q(I,J,K)=Q(I,J,K)-WEIGHT*R(K,I)
                RMS_INC=RMS_INC+R(K,I)*R(K,I)
              END DO
            END DO

            IF (JAC) THEN
              DO K=K_BEGIN,K_END
                DO I=IK_BEGIN,I_END,I_STEP
                  CORRECTIONS(I,J,K)=R(K,I)
                END DO
              END DO
            END IF

          END IF
CL    END OF IJ PLANE SWEEP
        END DO

CL    ADD CORRECTIONS FOR THE JACOBI METHOD

        IF(JAC) THEN
          DO K=K_BEGIN,K_END
            DO J=1,J_DIM
              IF(.NOT.(VERSION.EQ.2.AND.J.EQ.(J_DIM+1)/2)) THEN
                DO I=1,I_DIM
                  Q(I,J,K)=Q(I,J,K)-W*CORRECTIONS(I,J,K)
                END DO
              END IF
            END DO
          END DO
        END IF

CL    RESET POINT ORDERING OF IJ PLANE FOR SECOND SWEEP

        I = I_BEGIN
        I_BEGIN = I_END
        I_END   =  I
        I_STEP  = -I_STEP
        J = J_BEGIN
        J_BEGIN = J_END
        J_END   =  J
        J_STEP  = -1

CL    RESET ZEBRA SWITCH FOR SECOND SWEEP

        I_ZEBRA=1

      END DO

CL    END OF ROUTINE MG_K_LINE

      RETURN
      END
