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
C*LL
CLL   SUBROUTINE SLBHCADJ
CLL   -------------------
CLL
CLL   THIS ROUTINE IS FOR USE WITH THE 'SLAB' OCEAN MODEL ONLY.
CLL
CLL   MODIFIED VERSION:  2/8/1993
CLL
CLL   THE HEAT CONVERGENCE IS ADJUSTED TO PREVENT PROBLEMS
CLL   DUE TO LARGE NEGATIVE VALUES. SPECIFICALLY:-
CLL
CLL   IF THE HEAT CONVERGENCE AT A SEA-ICE POINT IS LESS
CLL   THAN HCLIMIT,THEN IT IS SET TO BE HCLIMIT, AND THE
CLL   HEAT IS REDISTRIBUTED OVER ALL THE SEA POINTS IN THE
CLL   HEMISPHERE. HCLIMIT = -40 WM2
CLL
CLL   THIS ROUTINE FORMS PART OF SYSTEM COMPONENT P40.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS
CLL   AND ENDDO STATEMENTS.
CLL
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   CALLED BY: SLABCNTL
CLL
CLL   WRITTEN BY A.B.KEEN (12/03/92)
CLL   MODIFIED BY A.B.KEEN (27/04/93)
CLL   MODIFIED BY A.B.KEEN (17/06/93)
CLL   MODIFIED BY C.A.SENIOR (28/02/94)
CLL   VERSION NUMBER 1.2
CLL   REVIEWER: W.INGRAM (01/03/93)
CLL
CLLEND---------------------------------------------------------------
C*L
      SUBROUTINE SLBHCADJ(L1,L2,
     +                    ADJHCONV,
     +                    WEIGHTS,
     +                    ICY,
     +                    HCLIMIT,
     +                    OPENSEA)
C
C
      INTEGER L1   ! IN SIZE OF INPUT DATA ARRAY
     +,L2          ! IN AMOUNT OF DATA TO BE PROCESSED
C
      REAL
     + ADJHCONV(L1)       ! INOUT ADJUSTED HEAT CONVERGENCE RATE (WM-2)
     +,WEIGHTS(L1)        ! IN WEIGHTS (COS LATITUDE) FOR AREA SUMS
     +,HCLIMIT      ! LIMIT FOR USING HEAT CONVERGENCES AT ICE POINTS
C
      LOGICAL
     + ICY(L1)            ! IN TRUE IF BOX CONTAINS ICE.
     +,OPENSEA(L1)        ! IN TRUE FOR OPEN SEA POINTS
C
C     VARIABLES LOCAL TO THIS ROUTINE ARE NOW DEFINED.
C
      REAL
     + SEASUM      ! SUM OF WEIGHTS OVER SEA POINTS
     +,HCONVSUM    ! SUM OF WEIGHTED HEAT CONVERGENCES
     +,HC_CORR     ! HEAT CONVERGENCE CORRECTION
C
      INTEGER
     + J           ! LOOP COUNTER
C
C    1. COMPUTE HEAT CONVERGENCE ADJUSTEMENTS
C
      SEASUM   = 0.0
      HCONVSUM = 0.0
      DO J=1,L2
          IF ( OPENSEA(J) ) THEN
            SEASUM = SEASUM + WEIGHTS(J)
          ENDIF
          IF (ICY(J) .AND. ( ADJHCONV(J) .LT. HCLIMIT ) ) THEN
            HCONVSUM    = HCONVSUM + ( ADJHCONV(J) - HCLIMIT )
     &                    * WEIGHTS(J)
            ADJHCONV(J) = HCLIMIT
          ENDIF
      END DO
      HC_CORR = HCONVSUM / SEASUM
      DO J=1,L2
          IF ( OPENSEA(J) ) THEN
            ADJHCONV(J) = ADJHCONV(J) + HC_CORR
          ENDIF
      END DO
C
C
C
      RETURN
      END
