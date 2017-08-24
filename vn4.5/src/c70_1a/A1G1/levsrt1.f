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
!+
! Subroutine Interface:
      SUBROUTINE LEVSRT(TYPE,NLEVS,IL,RL)
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
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

! Subroutine arguments:

!   Scalar arguments with intent(in):

      CHARACTER*1 TYPE
      INTEGER     NLEVS

!   Array arguments with intent(inout):

      REAL        RL(NLEVS)
      INTEGER     IL(NLEVS)

! Local variables:

      LOGICAL     LSWAP
      INTEGER     I
      INTEGER     J
      INTEGER     ILT
      REAL        RLT

!- End of Header ----------------------------------------------------

      DO 100 I=1,NLEVS
        LSWAP=.FALSE.
        DO 200 J=1,NLEVS-1
          IF(TYPE.EQ.'I') THEN
            IF(IL(J).GT.IL(J+1)) THEN
              LSWAP=.TRUE.
              ILT=IL(J)
              IL(J)=IL(J+1)
              IL(J+1)=ILT
            END IF
          ELSE
            IF(RL(J).LT.RL(J+1)) THEN
              LSWAP=.TRUE.
              RLT=RL(J)
              RL(J)=RL(J+1)
              RL(J+1)=RLT
            END IF
          END IF
  200   CONTINUE
        IF(.NOT.LSWAP) RETURN
  100 CONTINUE
      RETURN
      END
