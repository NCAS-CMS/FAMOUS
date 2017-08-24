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
CLL -------------- SUBROUTINE MERGEFLD ---------------------------------
CLL
CLL Purpose:  To merge an array held for points in boundary zone of a
CLL         rectangular area with an array covering the full area.
CLL Not suitable for single column use
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL  3.1  30/03/93 Correct rimweights at corners when J<I et sim. RTHB
CLL
CLL Programing standard: UM Documentation paper No4,
CLL                      Version No 2, dated 18/01/90
CLL
CLL System components covered: C72 (part)
CLL
CLL System task: C7
CLL
CLL Documentation: UM Documentation paper No C7,
CLL                draft version No 6, Dated 22/01/90
CLL                UM Documentation paper No 10,
CLL                draft version No 7, Dated 05/02/90
CLL
CLLEND

C*L Arguments

      SUBROUTINE MERGEFLD (ROW_LENGTH,ROW_SIZE,ROWS,RIMWIDTH,
     & RIMWEIGHTS,RIM,FIELD)
C*
      IMPLICIT NONE

C*L
      INTEGER
     &       ROW_LENGTH,  ! In Row length of input boundary data
     &       ROW_SIZE,    ! In Row length of data field
     &       ROWS,        ! In
     &       RIMWIDTH     ! In Width of boundary zone

      REAL
     &       RIMWEIGHTS(RIMWIDTH),  ! In Weights to be given to
C                                   !    boundary zone values.
     &       RIM((ROW_LENGTH+ROWS-2)*RIMWIDTH*2),
C                                   ! In Input boundary data.
     &       FIELD(ROWS*ROW_SIZE)   ! In/Out Output field

C*
C Local variables

      INTEGER
     &       I,           ! Loop over rim gridpoints.
     &       J,           ! Loop over N & S rows or E & W columns.
     &       IRIM         ! Position in RIM data array.

      REAL
     &       RWT          ! Modified rimweight for N & S rows.

CL Internal Structure

CL 1.0  Copy N rows into final positions

      IRIM = 1
      DO 10 I=1,RIMWIDTH
        DO 11 J=1,ROW_LENGTH
          IF (J .LT. I) THEN
            RWT = RIMWEIGHTS(J)
          ELSE IF (J .GT. ROW_LENGTH+1-I) THEN
            RWT = RIMWEIGHTS(ROW_LENGTH+1-J)
          ELSE
            RWT = RIMWEIGHTS(I)
          END IF
          FIELD(J+(I-1)*ROW_SIZE) = RIM(IRIM)*RWT
     &                              +FIELD(J+(I-1)*ROW_SIZE)*(1.0-RWT)
          IRIM=IRIM+1
 11     CONTINUE
 10   CONTINUE

CL 2.0  Copy E rows into final positions

      DO 20 I=RIMWIDTH+1,ROWS-RIMWIDTH
        DO 21 J=ROW_LENGTH-RIMWIDTH+1,ROW_LENGTH
          FIELD(J+(I-1)*ROW_SIZE)=RIM(IRIM)*RIMWEIGHTS(ROW_LENGTH+1-J)
     &                              +FIELD(J+(I-1)*ROW_SIZE)*(1.0-
     &                              RIMWEIGHTS(ROW_LENGTH+1-J))
          IRIM=IRIM+1
 21     CONTINUE
 20   CONTINUE

CL 3.0  Copy S rows into final positions

      DO 30 I=ROWS-RIMWIDTH+1,ROWS
        DO 31 J=1,ROW_LENGTH
          IF (J .LT. ROWS+1-I) THEN
            RWT = RIMWEIGHTS(J)
          ELSE IF (J .GT. ROW_LENGTH-ROWS+I) THEN
            RWT = RIMWEIGHTS(ROW_LENGTH+1-J)
          ELSE
            RWT = RIMWEIGHTS(ROWS+1-I)
          END IF
          FIELD(J+(I-1)*ROW_SIZE) = RIM(IRIM)*RWT
     &                              +FIELD(J+(I-1)*ROW_SIZE)*(1.0-RWT)
          IRIM=IRIM+1
 31     CONTINUE
 30   CONTINUE

CL 4.0  Copy W rows into final positions

      DO 40 I=RIMWIDTH+1,ROWS-RIMWIDTH
        DO 41 J=1,RIMWIDTH
          FIELD(J+(I-1)*ROW_SIZE)=RIM(IRIM)*RIMWEIGHTS(J)
     &                              +FIELD(J+(I-1)*ROW_SIZE)*(1.0-
     &                              RIMWEIGHTS(J))
          IRIM=IRIM+1
 41     CONTINUE
 40   CONTINUE

      RETURN
      END


