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
CLL Subroutine SWPTSC   ----------------------------------------------
CLL
CLL Purpose :
CLL  It calculates scaled pathlengths of each gaseous absorber for each
CLL  layer and returns them in DPATH for use by SWMAST, which sums them
CLL  to get the total scaled pathlengths for each beam considered, so
CLL  that the gaseous transmissivities can be calculated.
CLL
CLL Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled.
CLL                       S.J.Swarbrick
CLL   4.3    Feb. 97  T3E optimisation: code restructured, cray vector
CLL                    library functions introduced.
CLL                       D.Salmond & S.J.Swarbrick
CLL
CLL   4.5    Jan. 98  T3E optimisation:  rtor_v replaced by powr_v 
CLL                                      D.Salmond
CLL
CLL Programming standard :
CLL  It conforms to standard A of UMDP 4 (version 3, 07/9/90).
CLL  If UPDATE *DEF CRAY is off, a version is produced which except
CLL  for the addition of ! comments is standard FORTRAN 77 with no
CLL  8X-deprecated features (and which sets the "vector length" to 1)
CLL  but the standard version includes automatic arrays also.
CLL
CLL Logical components covered : P234
CLL  (interaction of shortwave radiation with the atmosphere)
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation:   UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SWPTSC (H2O, CO2, O3, PSTAR, AB, BB,
     &    L2,                                                           
     &    NLEVS, NWET, NOZONE, L1, DPATH)
C*
      INTEGER NGASES
C Number of absorbing gases treated in the shortwave
      PARAMETER (NGASES=3)    !  Standard set is water vapour, ozone
C                             !  and carbon dioxide.
C*L
      INTEGER!, INTENT (IN) ::
     &     L2,                       ! Number of points to be treated   
     &     NLEVS,                    ! Number of levels
     &     NWET,                     ! Number of levels with moisture -
C     ! above them a small value H2OMN is used (zero would give trouble)
     &     NOZONE,                   ! Number of levels with ozone data
C     ! provided - below them the lowest layer's is used.
     &     L1                        ! First dimension of input arrays
      REAL!, INTENT(IN) ::
     &     H2O(L1,NWET), CO2,        ! Mass mixing ratio (mK in UMDP 23)
     &     O3(L1,NOZONE),            !             of each absorbing gas
     &     PSTAR(L1),                ! Surface pressure
     &     AB(NLEVS+1), BB(NLEVS+1)  ! As & Bs at layer boundaries
      REAL!, INTENT(OUT) ::
     &     DPATH(L2,NGASES,NLEVS)
C     !  The scaled pathlengths are returned in DPATH, indexed by NGASES
C     !  1 is H2O, 2 is O3 & 3 is CO2
C*
CL    !  SWPTSC has no EXTERNAL calls and no significant structure
CL    !     but it has one dynamically allocated array, WORK.           
C
      REAL WORK(L2,2,2)
C     !  WORK is used to hold powers of layer boundary pressures used
C     !  in 2.3.1 and passed from one level to the next to save
C     !  re-calculation.  (This does prevent autotasking over levels.)
      REAL PSNH2O,                   ! Pressure scaling normalization
     &     PSNCO2,                   ! constants for water vapour & CO2
     &     PSXH2O,                   ! Pressure scaling exponents for
     &     PSXCO2,                   !               water vapour & CO2
     &     PX1H2O, PX1CO2,           ! 1+PSXH2O, 1+PSXCO2
     &     PRFH2O,                   ! Reference pressures for scaling
     &     PRFCO2,                   !               water vapour & CO2
     &     PSTRO3,                   ! Standard surface pressure for O3
     &     H2OMN                     ! Minimum water vapour pathlength
      REAL                           ! Pressure at top of current layer 
     &  power,pbot(l2,nlevs+1),pbot_h2o(l2,nlevs+1),      
     &     DPOBYG                    ! Pressure difference for ozone, /g
      INTEGER LEVEL, J,              ! Loopers over levels & points
     &     ONETWO,                   ! Flipper
     &     NDRY,                     ! Number of levels without moisture
     &     OLEVEL                    ! Index for the ozone to be used in
C                                    !                 the current level
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      PARAMETER ( PSTRO3 = 101325. )
      PARAMETER ( H2OMN = 1.E-10 )
      PARAMETER ( PSXH2O = 0.9, PX1H2O = 1. + PSXH2O, PRFH2O = 50000.,
     &            PSXCO2 = 0.7, PX1CO2 = 1. + PSXCO2, PRFCO2 = 25000. )
C     !  FORTRAN 77 will not allow the next two constants to be
C     !  defined in a PARAMETER statement, but the CRAY compiler will
C     !  give the same effect as if they were.
      PSNH2O = PRFH2O**(-PSXH2O) / (G*PX1H2O)
      PSNCO2 = PRFCO2**(-PSXCO2) / (G*PX1CO2)
C
      NDRY = NLEVS - NWET
C
C     ! Initialize the WORK term for CO2:
C
      power=px1co2                                                   
      do level=1, nlevs + 1
      DO 1 J=1, L2
       pbot(j,level) = ( PSTAR(J) * BB(level) + AB(level) ) 
    1 CONTINUE
      do j=1,L2
        pbot(j,level)=pbot(j,level)**power
      end do
      enddo
C
C     !   The next loop deals with H2O and CO2 for levels where there is
C     !   no moisture.  We just put a minimum value in for moisture,
C     !   treating CO2 as below.
C
      DO 2 LEVEL=1, NDRY
       DO 20 J=1, L2
        DPATH(J,1,LEVEL) = H2OMN
        DPATH(J,3,LEVEL) =
     &        CO2 * PSNCO2 * ( pbot(j,level+1) - pbot(j,level) )
   20  CONTINUE
    2 CONTINUE
       power=(PX1H2O/PX1CO2)                                         
C
C     !   This is the more general loop, calculating scaled pathlengths
C     !   for H2O and CO2.
C
      do LEVEL=NDRY+1, NLEVS + 1
      do j=1,L2
        pbot_h2o(j,level)=pbot(j,level)**power
      end do
      enddo
      DO 4 LEVEL=NDRY+1, NLEVS
       DO 40 J=1, L2
        IF (H2O(J,LEVEL-NDRY) .NE. 0.)  THEN
        DPATH(J,1,LEVEL) = H2O(J,LEVEL-NDRY) * PSNH2O *
     &                       ( pbot_h2o(j,level+1) - pbot_h2o(j,level) )
        ELSE
        DPATH(J,1,LEVEL) = H2OMN        
        ENDIF
        DPATH(J,3,LEVEL) =
     &        CO2 * PSNCO2 * ( pbot(j,level+1) - pbot(j,level) )  
   40  CONTINUE
    4 CONTINUE
C
C     !  Ozone has no pressure scaling, and to calculate the pathlengths
C     !  from the mass mixing ratios we use a "standard" surface
C     !  pressure, so that the climatology can be used without
C     !  interpolation but preserving total column ozone.  There are
C     !  thus no calculations in common with those for H2O and CO2, and
C     !  it is most conveniently treated quite separately, with no
C     !  repetition of code for wet and dry levels.
C
      DO 5 LEVEL=1, NLEVS
       DPOBYG = ( ( AB(LEVEL+1) - AB(LEVEL) ) + PSTRO3 *
     &                            ( BB(LEVEL+1) - BB(LEVEL) ) ) / G
       OLEVEL = MIN (LEVEL, NOZONE)
       DO 50 J=1, L2
        DPATH(J,2,LEVEL) = DPOBYG * O3(J,OLEVEL)
   50  CONTINUE
    5 CONTINUE
C
      RETURN
      END
