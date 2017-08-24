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
CL    !  but one dynamically allocated array MAXCON:                    
      REAL MAXCON(L2)                ! Maximum total cloud cover in the
C     ! layer currently being considered and those below it through
C     ! which cloud extends contiguously.
C*
      REAL TOCLE                     ! Total cloud in this layer
C*
      INTEGER LEVEL, J               ! Loopers over level & point
C
CL    !  First an initialization loop:
C
      DO J=1, L2
        CSF(J) = 1.
        MAXCON(J) = 0.
      ENDDO
C
CL    ! Then work up through the cloudy layers, remembering that LCA is
CL    ! the fractional cover by layer cloud outside the convective cloud
C
      DO 100 LEVEL=1, NCLDS
        DO 200 J=1, L2
C         !  So total cloud amount in a layer is (1-(LCA*(1-CCA)+CCA)))
C         !  if the convective cloud base extends through it, & (1-LCA)
C         !  if not.  We want the product of 1 - the maxima of this for
C         !  each group of contiguous (in the vertical) cloudy layers.
C         !  So accumulate this maximum through each such group, and at
C         !  each cloud-free layer multiply it in & re-zero it.
          IF ( LCA(J,LEVEL) .EQ. 0. .AND. ( CCA(J) .EQ. 0. .OR.
     &        LEVEL .LT. CCB(J) .OR. LEVEL .GE. CCT(J) )  ) THEN
             CSF(J) = CSF(J) * ( 1. - MAXCON(J) )
             MAXCON(J) = 0.
           ELSE
             TOCLE = LCA(J,LEVEL)
             IF ( LEVEL .GE. CCB(J)  .AND.  LEVEL .LT. CCT(J) )
     &          TOCLE = TOCLE + CCA(J) * ( 1. - TOCLE )
             MAXCON(J) = MAX(MAXCON(J),TOCLE)
          ENDIF
  200   CONTINUE
  100 CONTINUE
C
CL    !  The term from the highest cloud block has still to be put in if
CL    !    it extends into layer NCLDS.
      DO J=1, L2
        CSF(J) = CSF(J) * ( 1. - MAXCON(J) )
      ENDDO
C
      RETURN
      END
