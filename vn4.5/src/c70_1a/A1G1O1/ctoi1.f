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
C
!+Changes a character into an integer
! Subroutine Interface:
      SUBROUTINE CTOI(CHAR,LEN,INTEG,ERROR)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.5     Jul. 98    Change call to INTRFACE with call to C_MDI
!                      (A Van der Wal)
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------

! Subroutine arguments:

!   Scalar arguments with intent(in):

      INTEGER LEN

!   Scalar arguments with intent(out):

      CHARACTER*1 CHAR(LEN)
      CHARACTER*8 ERROR
      INTEGER INTEG
      INTEGER IFIRST
      INTEGER ILAST
      INTEGER IMULT
      INTEGER IPOS

! Local Scalars:

      LOGICAL LBLANK
      INTEGER I
      INTEGER IDIGIT

!- End of Header ------------------------------------------------------

      INTEG=IMDI
      ERROR='        '
      ILAST=LEN
C
      IFIRST=0
      IPOS=1
      LBLANK=.TRUE.
      DO 100 I=1,LEN
        IF((CHAR(I).NE.'0').AND.
     &  (CHAR(I).NE.'1').AND.
     &  (CHAR(I).NE.'2').AND.
     &  (CHAR(I).NE.'3').AND.
     &  (CHAR(I).NE.'4').AND.
     &  (CHAR(I).NE.'5').AND.
     &  (CHAR(I).NE.'6').AND.
     &  (CHAR(I).NE.'7').AND.
     &  (CHAR(I).NE.'8').AND.
     &  (CHAR(I).NE.'9').AND.
     &  (CHAR(I).NE.' ')) THEN
C         ILLEGAL CHARACTER IN INTEGER
          ERROR='M109    '
          RETURN
        ELSE
          IF(CHAR(I).EQ.' ') THEN
            IF(IPOS.EQ.2) ILAST=I-1
            IF(.NOT.LBLANK) THEN
              IPOS=IPOS+1
              LBLANK=.TRUE.
              IF(IPOS.GT.3) THEN
C               ILLEGAL INTEGER
                 ERROR='M109    '
                 RETURN
              END IF
            END IF
          ELSE
            IF(IPOS.EQ.1) IFIRST=I
            IF(LBLANK) THEN
              IPOS=IPOS+1
              LBLANK=.FALSE.
              IF(IPOS.GT.3) THEN
C               ILLEGAL INTEGER
                 ERROR='M109    '
                 RETURN
              END IF
            END IF
          END IF
        END IF
  100 CONTINUE
C
      INTEG=0
      IMULT=1
      DO 200 I=ILAST,IFIRST,-1
        READ(CHAR(I),210) IDIGIT
  210   FORMAT(I1)
        INTEG=INTEG+IDIGIT*IMULT
        IMULT=IMULT*10
  200 CONTINUE
      RETURN
      END
