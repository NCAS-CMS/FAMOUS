C ******************************COPYRIGHT******************************
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
C
      SUBROUTINE CHK1BOX(LENL,INDEX,NX,NY,ICODE,CMESSAGE)
CLL   Subroutine CHK1BOX ---------------------------------------------
CLL
CLL Purpose:
CLL
CLL   Check that each gridbox appears no more than once in a
CLL   supplied list of box indices. This routine is used to check
CLL   that it is valid to use the ADJUST mode of DO_AREAVER.
CLL
CLL   Programming Standard, paper 4 version 4 (14.12.90)
CLL
CLL Modification history:
CLL
CLL   Original version: J.M.Gregory 5.9.96 for HADCM3
CLL
CLL Logical components covered :
CLL
CLL Project task :
CLL
CLL External documentation: Unified Model documentation paper No:
CLL                         Version:
CLL
CLLEND -----------------------------------------------------------------
C
      IMPLICIT NONE
C*L
      INTEGER
     & LENL                    !IN length of index list
     &,NX,NY                   !IN dimensions of grid
     &,INDEX(LENL)             !IN indices of gridboxes
     &,ICODE                   !OUT return code
      CHARACTER
     & CMESSAGE*(*)            !OUT error message
C*
      INTEGER
     & IP                      ! pointer into list
     &,IX,IY                   ! working indices
     &,COUNT(NX,NY)            ! counts of occurrences of indexes
C
      DO IY=1,NY
      DO IX=1,NX
        COUNT(IX,IY)=0
      ENDDO
      ENDDO
C
      DO IP=1,LENL
        IX=MOD(INDEX(IP)-1,NX)+1
        IY=(INDEX(IP)-1)/NX+1
        COUNT(IX,IY)=COUNT(IX,IY)+1
      ENDDO
C
      ICODE=0
      DO IY=1,NY
      DO IX=1,NX
        IF (COUNT(IX,IY).GT.1) ICODE=1
      ENDDO
      ENDDO
C
      IF (ICODE.NE.0)
     &CMESSAGE='ADJUST mode of DO_AREAVER should not be used'
C
      RETURN
      END
