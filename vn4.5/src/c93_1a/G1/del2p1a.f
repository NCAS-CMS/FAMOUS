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
CLL  SUBROUTINE DEL SQUARED FFT P -------------------------------------
CLL
CLL  PURPOSE: Solves DEL SQUARED Q = RHS on surface of sphere.
CLL           Uses Fast Fourier Fransforms to decompose right-hand-
CLL           side in x-direction and then solves in y-direction
CLL           for the fourier coefficients of the solution.
CLL
CLL  N.B.  1. This version is for a problem where the solution Q and
CLL           the right-hand-side RHS are held at pressure points.
CLL        2. A solution to this equation exists uniquely, upto
CLL           an arbitrary constant if and only if the sum of the
CLL           right-hand-side values over the sphere is equal to zero.
CLL           This is known as the Compatability Condition.
CLL           This routine assumes that this condition is satisfied.
CLL           This routine chooses the arbitrary constant as follows;
CLL           If A(y) is the coefficient of wave number 0 then the
CLL           constant is chosen by setting the mean of A to zero.
CLL
CLL  NOT SUITABLE FOR SINGLE CLOUMN USE.
CLL  VERSION FOR CRAY Y-MP
CLL
CLL M.Mawson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.1  6/1/92   FPP directives to fix bug   (R.S.Bell)
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL                        STANDARD B.
CLL
CLL  SYSTEM COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK:
CLL
CLL  DOCUMENTATION: U.M.D.P. No 14 by M.H. Mawson
CLL
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE DEL_SQUARED_FFT_P
     1                            (Q,RHS,SEC_P_LATITUDE,COS_U_LATITUDE,
     2                             TRIGS,IFAX,LATITUDE_STEP_INVERSE,
     3                             EARTH_RADIUS_INVERSE,P_FIELD,
     4                             ROW_LENGTH)

      IMPLICIT NONE

      INTEGER
     &  P_FIELD            !IN Horizontal dimension of pressure field
     &, ROW_LENGTH         !IN Number of points on a row.
     &, IFAX(10)           !IN Holds factors of row_length used in FFT's

      REAL
     & RHS(P_FIELD)            !IN Holds right-hand-side.
     &,SEC_P_LATITUDE(P_FIELD) !IN  1./cos(lat) at pressure points
     &,COS_U_LATITUDE(P_FIELD-ROW_LENGTH) !IN cos(lat) at velocity point
     &,LATITUDE_STEP_INVERSE   !IN. 1/(delta phi)
     &,EARTH_RADIUS_INVERSE    !IN  1./(radius of earth).

      REAL
     & TRIGS(ROW_LENGTH)   !IN Holds trigonometric terms used in FFT's

      REAL
     & Q(P_FIELD)          !OUT Holds solution.

C*---------------------------------------------------------------------

C*L  DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C   DEFINE LOCAL ARRAYS: 5 ARE REQUIRED

      REAL
     &  RHS_DATA(ROW_LENGTH+2,P_FIELD/ROW_LENGTH) ! Fourier modes of
     &, Q_DATA(ROW_LENGTH+2,P_FIELD/ROW_LENGTH)   ! Q and RHS.
     &, A_DIAG(P_FIELD/ROW_LENGTH,ROW_LENGTH+2)   ! Matrix diagonal
     &, A_SUB_DIAG(P_FIELD/ROW_LENGTH,ROW_LENGTH+2) ! Sub diagonal
     &, A_SUP_DIAG(P_FIELD/ROW_LENGTH,ROW_LENGTH+2) ! Super diagonal

C*---------------------------------------------------------------------
C   DEFINE LOCAL VARIABLES

C   COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J,K,IK
     &, ROWS                ! Number of rows in field.
     &, LOT                 ! Number of data vectors passed to FFT's.
     &, JUMP                ! Number of storage locations between the
     &                      ! start of consecutive data vectors.
     &, INCREMENT           ! Number of storage locations between each
     &                      ! element of the same data vector, 1, if
     &                      ! consecutive.
     &, FFT_ISIGN            ! Parameter determining whether spectral
     &                       ! to grid-point (1) or grid-point to
     &                       ! spectral (-1) FFT's are required.

      REAL
     & SCALAR               ! Generic real work variable.
     &,FACTOR               ! Holds factor in matrix gaussian
     &                      ! elimination.
     &,WAVE_NUMBER          ! Holds wave number for which fourier
     &                      ! coefficients are being calculated.

C*L  EXTERNAL SUBROUTINE CALLS:------------------------------------
      EXTERNAL FOURIER
C*---------------------------------------------------------------------

CL  MAXIMUM VECTOR LENGTH ASSUMED IS ROW_LENGTH+2
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL

CL---------------------------------------------------------------------
CL    SECTION 1.     Put all RHS data in bigger array before calling
CL                   fourier decomposition. Two extra addresses per row
CL                   required.
CL---------------------------------------------------------------------

      ROWS = P_FIELD/ROW_LENGTH

CFPP$ NODEPCHK
! Fujitsu vectorization directive
!OCL NOVREC
      DO 100 J=1,ROWS
        IK = (J-1)*ROW_LENGTH
        DO 110 I=1,ROW_LENGTH
          RHS_DATA(I,J) = RHS(IK+I)
 110    CONTINUE
 100  CONTINUE

CL---------------------------------------------------------------------
CL    SECTION 2.     Call FOURIER to get fourier decomposition of data.
CL---------------------------------------------------------------------

      INCREMENT = 1
      JUMP = ROW_LENGTH+2
      FFT_ISIGN = -1
      LOT = ROWS
      CALL FOURIER(RHS_DATA,TRIGS,IFAX,INCREMENT,JUMP,ROW_LENGTH,
     &             LOT,FFT_ISIGN)

CL---------------------------------------------------------------------
CL    SECTION 3.     Solve equation.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 3.1 Solve for real constant term.
C ---------------------------------------------------------------------

CL i) Set up tri-diagonal matrix left-hand-side.

      SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &         LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
      J=1
      WAVE_NUMBER = 0.
      A_SUB_DIAG(1,J) = 0.
      A_SUP_DIAG(ROWS,J) = 0.

C set up A matrix.
C Northern boundary condition.
      I=1
      K = (I-1)*ROW_LENGTH+1
      A_DIAG(I,J) = -2.*SEC_P_LATITUDE(K)*SCALAR*
     &                     COS_U_LATITUDE(K)
      A_SUP_DIAG(I,J) = 2.*SEC_P_LATITUDE(K)*SCALAR*
     &                        COS_U_LATITUDE(K)

C Inner points.
      DO 310 I=2,ROWS-1
        K = (I-1)*ROW_LENGTH+1
        A_DIAG(I,J) = -SEC_P_LATITUDE(K)*SCALAR*
     &                  (COS_U_LATITUDE(K)+COS_U_LATITUDE(K-ROW_LENGTH))
        A_SUP_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*
     &                      COS_U_LATITUDE(K)
        A_SUB_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*
     &                      COS_U_LATITUDE(K-ROW_LENGTH)
 310  CONTINUE

C Southern boundary condition.
      I=ROWS
      K = (I-1)*ROW_LENGTH+1
      A_DIAG(I,J) = -2.*SEC_P_LATITUDE(K)*SCALAR*
     &                   COS_U_LATITUDE(K-ROW_LENGTH)
      A_SUB_DIAG(I,J) = 2.*SEC_P_LATITUDE(K)*SCALAR*
     &                      COS_U_LATITUDE(K-ROW_LENGTH)

CL ii) Solve matrix system.

      A_DIAG(1,J) = 1./A_DIAG(1,J)
      DO 312 I=2,ROWS
        FACTOR = A_SUB_DIAG(I,J) * A_DIAG(I-1,J)
        A_DIAG(I,J) = 1./(A_DIAG(I,J) - FACTOR*A_SUP_DIAG(I-1,J))
        RHS_DATA(J,I) = RHS_DATA(J,I) - FACTOR*RHS_DATA(J,I-1)
 312  CONTINUE

C Back substitute to get solution.

      Q_DATA(J,ROWS) = A_DIAG(ROWS,J)*RHS_DATA(J,ROWS)
      DO 314 I= ROWS-1,1,-1
        Q_DATA(J,I) = A_DIAG(I,J)*(RHS_DATA(J,I)-
     &                             A_SUP_DIAG(I,J)*Q_DATA(J,I+1))
 314  CONTINUE

C Set constant imaginery mode to zero.
C Remove arbitrary constant from real constant mode.
C Remove mean value as guess to arbitrary constant.

      DO I=1,ROWS
        Q_DATA(2,I) = 0.
      END DO
      SCALAR = 0.
      DO I=1,ROWS
        SCALAR = SCALAR + Q_DATA(1,I)
      END DO
      SCALAR = SCALAR / ROWS
CFPP$ NOINNER
      DO I=1,ROWS
        Q_DATA(1,I) = Q_DATA(1,I) - SCALAR
      END DO

C ---------------------------------------------------------------------
CL    SECTION 3.2 Solve for wave-number > 0 modes.
CL                Solution at poles is zero and this is thus the
CL                boundary condition for the solver.
C ---------------------------------------------------------------------

CL Set solution at poles = 0 for wave numbers > 0

CFPP$ NOINNER
      DO I=3,ROW_LENGTH+2
        Q_DATA(I,1)    = 0.
        Q_DATA(I,ROWS) = 0.
      END DO

CL i) Set up tri-diagonal matrix system left-hand-side.
C     Real coefficients only.

      SCALAR = EARTH_RADIUS_INVERSE*EARTH_RADIUS_INVERSE*
     &         LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
      DO J=3,ROW_LENGTH+2,2
        WAVE_NUMBER = (J+1)/2-1
        A_SUB_DIAG(2,J) = 0.
        A_SUP_DIAG(ROWS-1,J) = 0.

C Set up A matrix.
C Northern boundary condition.
        I=2
        K = (I-1)*ROW_LENGTH+1

        A_SUP_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*COS_U_LATITUDE(K)
        A_DIAG(I,J) = -SEC_P_LATITUDE(K)*SCALAR*
     &                  (COS_U_LATITUDE(K)+COS_U_LATITUDE(K-ROW_LENGTH))
     &                  - WAVE_NUMBER*WAVE_NUMBER*EARTH_RADIUS_INVERSE
     &                    *EARTH_RADIUS_INVERSE*SEC_P_LATITUDE(K)
     &                    *SEC_P_LATITUDE(K)

C Inner points.
        DO 320 I=3,ROWS-2
          K = (I-1)*ROW_LENGTH+1
          A_SUP_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*
     &                        COS_U_LATITUDE(K)
          A_SUB_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*
     &                        COS_U_LATITUDE(K-ROW_LENGTH)
          A_DIAG(I,J) = -SEC_P_LATITUDE(K)*SCALAR*
     &                  (COS_U_LATITUDE(K)+COS_U_LATITUDE(K-ROW_LENGTH))
     &                  - WAVE_NUMBER*WAVE_NUMBER*EARTH_RADIUS_INVERSE
     &                    *EARTH_RADIUS_INVERSE*SEC_P_LATITUDE(K)
     &                    *SEC_P_LATITUDE(K)
 320    CONTINUE

C Southern boundary condition.
        I=ROWS-1
        K = (I-1)*ROW_LENGTH+1
        A_SUB_DIAG(I,J) = SEC_P_LATITUDE(K)*SCALAR*
     &                        COS_U_LATITUDE(K-ROW_LENGTH)
        A_DIAG(I,J) = -SEC_P_LATITUDE(K)*SCALAR*
     &                  (COS_U_LATITUDE(K)+COS_U_LATITUDE(K-ROW_LENGTH))
     &                  - WAVE_NUMBER*WAVE_NUMBER*EARTH_RADIUS_INVERSE
     &                    *EARTH_RADIUS_INVERSE*SEC_P_LATITUDE(K)
     &                    *SEC_P_LATITUDE(K)
      END DO

C Matrix for imaginery coefficients is copy of real one for same wave
C number.
      DO 322 J=4,ROW_LENGTH+2,2
        DO I=2,ROWS-1
          A_DIAG(I,J) = A_DIAG(I,J-1)
          A_SUB_DIAG(I,J) = A_SUB_DIAG(I,J-1)
          A_SUP_DIAG(I,J) = A_SUP_DIAG(I,J-1)
        END DO
 322  CONTINUE

CL ii) Solve matrix system for all right-hand-sides.

CFPP$ NOINNER
      DO 324 J=3,ROW_LENGTH+2
        A_DIAG(2,J) = 1./A_DIAG(2,J)
 324  CONTINUE
      DO 326 I=3,ROWS-1
        DO J=3,ROW_LENGTH+2
          FACTOR = A_SUB_DIAG(I,J) * A_DIAG(I-1,J)
          A_DIAG(I,J) = 1./(A_DIAG(I,J) - FACTOR*A_SUP_DIAG(I-1,J))
          RHS_DATA(J,I) = RHS_DATA(J,I) - FACTOR*RHS_DATA(J,I-1)
        END DO
 326  CONTINUE

C Back substitute to get solution.

      DO 328 J=3,ROW_LENGTH+2
        Q_DATA(J,ROWS-1) = A_DIAG(ROWS-1,J)*RHS_DATA(J,ROWS-1)
 328  CONTINUE
      DO 329 I= ROWS-2,2,-1
        DO J=3,ROW_LENGTH+2
          Q_DATA(J,I) = A_DIAG(I,J)*(RHS_DATA(J,I)-
     &                              A_SUP_DIAG(I,J)*Q_DATA(J,I+1))
        END DO
 329  CONTINUE

CL---------------------------------------------------------------------
CL    SECTION 4.     Call FOURIER to create grid-point fields from
CL                   fourier solution modes.
CL---------------------------------------------------------------------

      INCREMENT = 1
      JUMP = ROW_LENGTH+2
      FFT_ISIGN = 1
      CALL FOURIER(Q_DATA,TRIGS,IFAX,INCREMENT,JUMP,ROW_LENGTH,
     &             LOT,FFT_ISIGN)

CL---------------------------------------------------------------------
CL    SECTION 5.     Copy solution into output array.
CL---------------------------------------------------------------------

      DO 500 J=1,ROWS
        IK = (J-1)*ROW_LENGTH
        DO 510 I=1,ROW_LENGTH
            Q(IK+I) = Q_DATA(I,J)
 510    CONTINUE
 500  CONTINUE

CL END OF ROUTINE DEL_SQUARED_FFT_P

      RETURN
      END
