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
CLL Subroutine SWMSAL  -----------------------------------------------
CLL
CLL Purpose :
CLL  It is part of component P234 (interaction of shortwave radiation
CLL  with the atmosphere)
CLL  It modifies the surface albedo to allow crudely for multiple
CLL  reflections.
CLL     Release 2.8 of the UM modified to allow for direct &
CLL  diffuse surface albedos being different.
CLL  It is suitable for single column model use.
CLL
CLL                      Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled.
CLL                       S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms to programming standard A of UMDP 4, version 3 (07/9/90
CLL  Except for containing ! comments, it conforms to the FORTRAN 77
CLL  standard with no features deprecated by 8X if *DEF CRAY is off:
CLL  otherwise it contains automatic arrays.
CLL
CLL Logical components covered : P234
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation: UMDP23 sub-section "Modifications to the
CLL  surface albedo".
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SWMSAL (TSA, LCDDR, LCA, CCDDR, CCA, CCB, OFFSET,
     &     L2,
     &     L1, NBANDS, NCLDS,                  MSA)
      INTEGER!, INTENT (IN)
     &     L1,                       ! First dimension of input arrays
     &     L2,                       ! Number of points to be treated   
     &     NCLDS,                    ! Number of layers with cloud
     &     NBANDS                    ! Number of bands
      REAL!, INTENT (IN)
     &     TSA(L1,NBANDS,2),         ! True surface albedo - mean over
C     !  the whole grid-box, direct-beam value followed by diffuse-beam
     &     LCA(L1,NCLDS),            ! Layer cloud amount and
     &     LCDDR(L2,NBANDS,NCLDS),   ! diffuse/diffuse reflectivity a5
     &     CCA(L1), CCDDR(L2,NBANDS) ! Same for convective cloud
C     !  - except that LCDDR has been multiplied by LCA and so is the
C     !  mean value over the grid-box, or at least that part of the
C     !  grid-box not occupied by any convective cloud at that level,
C     !  while CCDDR is a mean over the convective cloud only.
      INTEGER!, INTENT (IN)
     &     CCB(L1),                  ! Convective cloud base
     &     OFFSET                    ! Allows for CCB being numbered
C     ! from the top of the model down, while LCA & LCDDR begin in the
C     !                          first layer where cloud is allowed
      REAL!, INTENT (OUT) ::
     &     MSA(L2,NBANDS,2)          ! Modified surface albedo
C
CL    !  SWMSAL has no EXTERNAL calls and no significant structure
CL    !  but two dynamically allocated arrays, REFRAC & VISFRC.         
C
C*
      REAL REFRAC(L2),               ! The sky's fractional reflectivity
     &     VISFRC(L2)                ! The fraction of the sky at the
C     ! current level (and so, given random overlap, of the cloud in
C     ! that level) visible from the surface.
      REAL MODF                      ! Modification factor which
C                                    !           converts TSA into MSA
C     !
      INTEGER BAND, LEVEL, J         ! Loopers over band, level & point
C     !
      DO 100 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
       DO 110 J=1, L2
C       !
C       !  First, accumulate through the "DO 101" loop mean cloud
C       !  reflectivity over all the box, with weighting by the area of
C       !  each cloud visible at the surface.
C       !
        REFRAC(J) = 0.
        VISFRC(J) = 1.
  110  CONTINUE
       DO 101 LEVEL=NCLDS, 1, -1
Cfpp$   Select(CONCUR)
        DO 111 J=1, L2
C        !  Since LCA is in fact the fractional cover by layer cloud
C        !  outside the convective cloud, we can do the calculations
C        !  just by working up, allowing for the effects of each cloud
C        !  on VISFRC and REFRAC as we reach its base, and treating
C        !  convective cloud as if it were just below where it actually
C        !  is.
         IF ( CCB(J) .EQ. (LEVEL+OFFSET) ) THEN
           REFRAC(J) = REFRAC(J) + VISFRC(J) * CCA(J) * CCDDR(J,BAND)
           VISFRC(J) = VISFRC(J) * ( 1. - CCA(J) )
         ENDIF
         REFRAC(J) = REFRAC(J) + VISFRC(J) * LCDDR(J,BAND,LEVEL)
         VISFRC(J) = VISFRC(J) * ( 1. - LCA(J,LEVEL) )
  111   CONTINUE
  101  CONTINUE
C      !
Cfpp$  Select(CONCUR)
       DO 100 J=1, L2
        MODF = ( 1. - REFRAC(J) ) / ( 1. - TSA(J,BAND,2) * REFRAC(J) )
        MSA(J,BAND,1) = TSA(J,BAND,1) * MODF
        MSA(J,BAND,2) = TSA(J,BAND,2) * MODF
C
  100 CONTINUE
C
      RETURN
      END
