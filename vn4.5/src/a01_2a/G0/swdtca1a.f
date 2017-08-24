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
CLL Subroutine SWDTCA   ----------------------------------------------
CLL
CLL Purpose :
CLL  It calculates a total cloud amount diagnostic, the fraction of the
CLL  grid-box where cloud exists at some level(s), consistent with the
CLL  random overlap assumption used in the SW and whichever cloud
CLL  distribution is being used there.
CLL  Suitable for single column model use.
CLL
CLL                      Author: William Ingram 15 Oct 1992
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL  It conforms to programming standard A of UMDP 4 (version 2 18/1/90)
CLL  and has no features deprecated by 8X.
CLL  The code is standard FORTRAN 77 except for having ! comments
CLL
CLL Logical components covered : D23 (radiation diagnostics)
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
      SUBROUTINE SWDTCA (LCA, CCA, NCLDS, L1, L2, TCA)
C
      INTEGER!, INTENT (IN)
     &     L1,                       ! First dimension of input arrays
     &     L2,                       ! Number of points to be treated
     &     NCLDS                     ! Number of layers with cloud
      REAL!, INTENT (IN)
     &     LCA(L1,NCLDS), CCA(L1)    ! Layer & convective cloud fraction
      REAL!, INTENT (OUT) ::
     &     TCA(L1)                   ! Total cloud amount
C
CL    !  SWDTCA has no EXTERNAL calls
CL    !  and no dynamically allocated workspace.
C*
      INTEGER LEVEL, J               ! Loopers over level & point
C
CL    ! Since LCA is the fractional cover by layer cloud outside the
CL    ! convective cloud, can just multiply the cloud-free fractions
CL    ! together for each cloud present:
C     !
C
      DO J=1, L2
        TCA(J) = 1. - CCA(J)
      ENDDO
      DO 100 LEVEL=1, NCLDS
        DO J=1, L2
          TCA(J) = TCA(J) * ( 1. - LCA(J,LEVEL) )
        ENDDO
  100 CONTINUE
C
CL    !  Finally, convert clear fraction into total cloud amount:
C
      DO J=1, L2
        TCA(J) = 1. - TCA(J)
      ENDDO
C
      RETURN
      END
