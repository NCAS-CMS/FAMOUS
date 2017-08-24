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
CLL Subroutine LWDCSF  ------------------------------------------------
CLL
CLL        Purpose :
CLL  It calculates the clear-sky fraction (i.e. the fraction of the
CLL  grid-box where no cloud exists at any level) for use in
CLL  constructing diagnostics - the "total cloud amount" diagnostic (the
CLL  fraction of the grid-box where there is some cloud at some level)
CLL  is one minus this, and the clear-sky (type I) diagnostics are the
CLL  clear-sky (type II) diagnostics multiplied by it.
CLL  Suitable for single column model use.
CLL
CLL       Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.3  15/12/93   Corrected to allow for zero convective cloud cover
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled.
CLL                       S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms to programming standard A of UMDP 4 (version 2 18/1/90)
CLL  and has no features deprecated by 8X.
CLL  If *DEF IBM or *DEF RANDOVER or both are set, the code is standard
CLL  FORTRAN 77 except for having ! comments (it then sets the "vector
CLL  length" to be 1) but otherwise it includes CRAY automatic arrays
CLL  also.
CLL
CLL Logical components covered : D23 (radiation diagnostics)
CLL
CLL Project task : P23
CLL
CLL External documentation:
CLL  The cloud overlap assumptions are documented in UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE LWDCSF (LCA, CCA, CCB, CCT, NCLDS, L1,
     &     L2,                                                          
     &     CSF)
C
      INTEGER!, INTENT (IN)
     &     L1,                       ! First dimension of input arrays
     &     L2,                       ! Number of points to be treated   
     &     NCLDS,                    ! Number of layers with cloud
     &     CCB(L1),                  ! Convective cloud base & top,
     &     CCT(L1)                   !  counting upward & the surface=1
      REAL!, INTENT (IN)
     &     LCA(L1,NCLDS), CCA(L1)    ! Layer & convective cloud fraction
C*IF -DEF,RANDOVER                                                      
C     !  Array dimensions must be constants in FORTRAN:
C*CALL L2VAL                                                            
C*ENDIF -RANDOVER                                                       
      REAL!, INTENT (OUT) ::
     &     CSF(L2)                   ! Clear-sky fraction returned
C
CL    !  LWDCSF has no EXTERNAL calls
CL    !  and no dynamically allocated workspace.
C*
      INTEGER LEVEL, J               ! Loopers over level & point
C
CL    !  First an initialization loop:
C
      DO J=1, L2
        CSF(J) = 1. - CCA(J)
      ENDDO
C
CL    ! Then work up through the cloudy layers, remembering that LCA is
CL    ! the fractional cover by layer cloud outside the convective cloud
C
      DO 100 LEVEL=1, NCLDS
        DO 200 J=1, L2
C         !  Thus we can just multiply the cloud-free fractions together
C         !  for each layer of layer cloud, putting in the CCA term at
C         !  the beginning - if we did do it layer-by-layer the term
C         !  would be (1-(LCA*(1-CCA)+CCA))) in the layer where the
C         !  convective cloud base is, & (1-LCA) elsewhere.
          CSF(J) = CSF(J) * ( 1. - LCA(J,LEVEL) )
  200   CONTINUE
  100 CONTINUE
C
      RETURN
      END
