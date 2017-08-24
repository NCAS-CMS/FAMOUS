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
CLL Subroutine SET_LEVELS_LIST
CLL
CLL Purpose : To set up a list of levels at which a diagnostic is
CLL           required, using information in the STASH list.
CLL Service routine  version for Cray YMP
CLL
CLL W.Ingram    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL
CLL Programming Standard : Unified Model Documentation paper number 4
CLL                      : Version no 2, dated 18/01/90
CLL
CLL System components covered : D3
CLL
CLL System task : P0
CLL
CLL Documentation: U.M. Documentation paper number P0,C4
CLL
CLLEND

C*L Arguments

      SUBROUTINE SET_LEVELS_LIST(LEVELS,
     &                    LEN_STLIST,STLIST,LIST,STASH_LEVELS,
     &      LEN_STASHLEVELS,ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     &       LEVELS,      ! IN Number of levels in input data
     &       LEN_STLIST,  ! IN
     &       STLIST(LEN_STLIST), ! IN STASH list
     &  LEN_STASHLEVELS,
     &  STASH_LEVELS(LEN_STASHLEVELS,*),  ! IN - list of levels required
     &       ICODE        ! OUT Return code =0 Normal exit
C                                       >1 Error message

      LOGICAL
     &       LIST(LEVELS) ! OUT List of levels required.

      CHARACTER*(80) CMESSAGE ! Error message

CL Local variables

      INTEGER
     &      K,
     &      KOUT


CL Initialise levels list to false

      DO K=1,LEVELS
        LIST(K)= .FALSE.
      END DO

CL Check for method of levels selection
CL Levels list must be present.

      IF(STLIST(10).LT.0) THEN

C Set logical array list to identify levels required.

        DO KOUT=2,STASH_LEVELS(1,-STLIST(10))+1
          IF((STASH_LEVELS(KOUT,-STLIST(10)).GE.1).AND.
     &    (STASH_LEVELS(KOUT,-STLIST(10)).LE.LEVELS)) THEN
C         LEVEL IS IN THE RANGE OF LIST.
              LIST(STASH_LEVELS(KOUT,-STLIST(10))) =.TRUE.
          ELSE
C         LEVEL IS OUT OF THE RANGE OF LIST.
              CMESSAGE=  ' SET_LEVELS_LIST: level out of range'
              WRITE(6,*) ' SET_LEVELS_LIST: level out of range'
              WRITE(6,*) ' level=',STASH_LEVELS(KOUT,-STLIST(10))
              WRITE(6,*) ' Section, Item =',STLIST(2),STLIST(1)
              ICODE=2
          END IF
        END DO


      ELSE

CL Illegal control data

        ICODE=1
        CMESSAGE='SET_LEVELS_LIST: Illegal control data'
      WRITE(6,*) 'Illegal control data SET_LEVELS_LIST,STLIST(10,11)=',
     &         STLIST(10) ,STLIST(11)
      WRITE(6,*) 'Section and item numbers ',STLIST(2),STLIST(1)
        RETURN

      END IF

      RETURN
      END
