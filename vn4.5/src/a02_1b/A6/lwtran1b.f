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
CLL Subroutine LWTRAN   ----------------------------------------------
CLL
CLL Purpose :
CLL  It calculates clear-sky transmissivities in each of the NBANDS
CLL  longwave spectral bands (and, optionally, additional diagnostic
CLL  ones) from the pathlengths for each effective absorbing gas in
CLL  each band.  (Where the absorption by a gas includes terms with
CLL  different pathlength scaling, like water vapour line & continuum,
CLL  they are treated as two gases.)
CLL     The version of routine LWTRAN used in
CLL  version 1B (gaseous effects treated as Morcrette et al, 1986) of
CLL  the UM LW code, and a dummy version of routine LWLKIN for
CLL  compatibility with other versions.
CLL  Version 3, part of the alternative code giving "ECMWF-like"
CLL  treatment of LW gaseous transmissivities, following Morcrette et al
CLL  (J.-J. Morcrette, L.D. Smith & Y. Fouquart, 1986, Beitr. Phys.
CLL  Atmosph., 59, 455-469).  All the calculations are changed: instead
CLL  of look-up tables for line absorption, ozone uses a Malkmus model,
CLL  and CO2 and water vapour use Horner's algorithm, while the
CLL  continuum terms, though still exponential, are calculated rather
CLL  indirectly to save evaluation of exponentials.  Also, different
CLL  pathlengths are generally used for different bands, and the
CLL  pathlengths are assumed to contain a diffusivity factor.  Version 3
CLL  of LWTRAN was set up from version 2.1 to be part of version 1B
CLL  (ECMWF-like gaseous transmissivities) of the LW from release 2.7 of
CLL  the UM.                            William Ingram 22 June 1992
CLL
CLL           Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: routine optimised using cray vector
CLL                   library functions & appropriately restructured.
CLL                              D.Salmond & S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL  If UPDATE *DEF CRAY is off, the code is standard FORTRAN 77 except
CLL  for having ! comments (it then sets the "vector length" to 1) but
CLL  otherwise it includes an automatic array also.
CLL
CLL Logical components covered : P232
CLL  P232 (longwave radiation),
CLL  It is also intended to be easily extended to perform
CLL  some of the functions of D23 (radiation diagnostics), by diagnosing
CLL  additional transmissivities.
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation:      UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE LWTRAN (PATH, DUMMY1, DUMMY2, DUMMY3, DUMMY4,
     &     L,                                                           
     &     TRANS)
C*
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
!  "ECMWF" bands as described by Morcrette et al (J.-J. Morcrette,
!   L.D. Smith & Y. Fouquart, 1986, Beitr. Phys. Atmosph., 59, 455-469)
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=12)
C     ! This set is for the "ECMWF-like" code.  Gases mostly have
C     ! different pressure and temperature scaling of their pathlengths
C     ! in each band, so that there are 6 absorber amounts for water
C     ! vapour line absorption, 1 for each of the foreign-broadened &
C     ! self-broadened water vapour continua, 2 for CO2 (for 3 different
C     ! bands) and 2 for ozone (with and without pressure broadening,
C     ! though only in one band).
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
C*L
      INTEGER!, INTENT(IN) ::
     &     L,                    ! Number of points                     
     &     DUMMY3                ! Not used in this version, but left
C                                !    as an argument for compatibility
      REAL!, INTENT(IN) ::
     &     PATH (L,NGASES),      ! Scaled pathlengths for each gas
     &     DUMMY1, DUMMY2,       ! Not used in this version, but left
     &     DUMMY4                !    as arguments for compatibility
      REAL!, INTENT(OUT) ::
     &     TRANS(L,NTRANS)       ! Transmissivities
C*
CL    !  No EXTERNAL routines called
C     ! *COMDECK LWGAFD holds numbers used to calculate LW gaseous
C     !   transmissivities from scaled pathlengths similarly to
C     !   the ECMWF LW code.
      REAL O3M1, O3M2, PIALF, O3WT1, O3WT2
CL    ! Coefficients for Malkmus model of ozone absorption
      PARAMETER ( O3M1 = 2230., O3M2 = 100., PIALF = 2., O3WT1 = 0.7554,
     &     O3WT2 = 1. - O3WT1 )
C
      REAL PADE(3,8)          ! Pade/ approximant coefficients for water
C     !  vapour line and CO2 absorption - second subscript indexes water
C     !  vapour in each band in standard order, & then CO2 in band 2 and
C     !  bands-3-&-4.  (N.B. This matches the order of filling of GA &
C     !  GB in JJM's original code, but means that 3-5 are a permutation
C     !  of ZTT/PTT(3-5) there.)  In the original code there are
C     !  different Pade/ fits for each temperature from 187.5K to
C     !  312.5K: here a single one is used for each band, for 237.5K,
C     !  275K, 262.5K, 287.5K, 287.5K & 275K respectively for water
C     !  vapour.  These temperatures were chosen subjectively by
C     !  considering the position of the cooling peak in each band - not
C     !  always a clear choice, but there is in fact little sensitivity.
C     !  For CO2, the 225K sets are used.
      DATA PADE /
     &    0.73585943E-02,  -0.10847662E-02,   0.10475952E+00,
     &    0.12167192E+01,   0.52341830E+00,   0.10800762E+01,
     &    0.79978921E+01,   0.71929934E+01,   0.73878952E+01,
     &    0.24063614E+02,   0.10392022E+02,   0.10509317E+02,
     &    0.18097099E+00,  -0.25423873E-01,   0.42353379E+00,
     &   -0.37133165E+01,   0.44809588E+00,  -0.81329826E+01,
     &    0.77659686E-01,   0.12191543E+01,   0.20855896E+01,
     &    0.13213894E+02,   0.22259478E+02,   0.22268925E+02 /
C
      REAL SRAP,                           ! SQRT ( scaled pathlength )
     &     TCO22,                          ! CO2 transmissivities for
     &     TCO234,                         ! band 2, & for both 3 & 4.
     &     FBWCP1,                         !   Foreign- and Self-
     &     SBWCP1,                         ! Broadened Water vapour
     &     FBWCP2, FBWCP4, FBWCP8,         ! Continuum terms, being two
     &     FBWC16, FBWC32, FBWCPQ, FBWCPH, ! particular exponentials to 
     &     SBWCP2, SBWCP4, SBWCP8,         ! the Powers 1, 2, 4, 8, 16,
     &     SBWC16, SBWC32, SBWCPQ, SBWCPH, ! 32, a Half & a Quarter.    
     &     FBWCB2, FBWCB3,                 !   Foreign- and Self-
     &     FBWCB4, FBWCB5,                 ! Broadened Water vapour
     &     SBWCB2, SBWCB3,                 ! Continuum transmissivities,
     &     SBWCB4, SBWCB5,                 ! in Bands 2-5
     &     TO31,                           ! Ozone transmissivities in
     &     TO32,                           !  2 fractions of the band.
     &     UXY(2*L),                       ! Used for the Malkmus       
     &     VXY(2*L)                        ! calculation of TO31,2 -    
!  VXY is the ozone-pathlength-weighted mean pressure over the          
!    pathlength, & UXY is twice the unscaled path divided by VXY.       
      INTEGER J, BAND,K                    ! Loopers over point & band  
C
! Local workspace used for t3e optimisation
      REAL sqrt_uxy(2*L)
      REAL exp_vxy(2*L)
      REAL utemp,vtemp
      REAL EXPPATH(L,4)  
C
      DO JTRANS=1, NTRANS                                             
        DO J=1, L                                                    
CL        !   Horner's algorithm for H2O transmission
          SRAP  = sqrt(PATH(j,jtrans))
          TRANS(J,JTRANS) = ( PADE(1,JTRANS) + SRAP * PADE(2,JTRANS) )
     &          /( PADE(1,JTRANS)+SRAP*(PADE(3,JTRANS)+SRAP) )          
        END DO                                                          
      END DO        
C
      DO J=1,L
        EXPPATH(J,1)=-0.002*PATH(J, 9)
        EXPPATH(J,2)=0.25*EXPPATH(J,1)
        EXPPATH(J,3)=      -PATH(J,10)                                  
        EXPPATH(J,4)=0.25*EXPPATH(J,3)
      END DO
      DO K=1,4
        DO J=1,L
          EXPPATH(J,K)=EXP(EXPPATH(J,K))
        END DO
      END DO                              
      DO J = 1, L
      utemp = 4. * PATH(J,11) * PATH(J,11)
     &          /( PIALF * PATH(J,12) )
      uxy(j)=1.+O3M1*utemp
      uxy(j+l)=1.+O3M2*utemp
      ENDDO
      DO J=1,2*L
      sqrt_uxy(j)=sqrt(uxy(j))
      ENDDO
      DO J = 1, L
      vtemp = PIALF * PATH(J,12)
     &         /( 2. * PATH(J,11) )
      vxy(j)=-(sqrt_uxy(j)-1.)*vtemp
      vxy(j+l)=-(sqrt_uxy(j+l)-1.)*vtemp
      ENDDO
      DO J=1,2*L
      exp_vxy(j)=exp(vxy(j))
      ENDDO
      DO J = 1, L                                                       
CL      !   Horner's algorithm for CO2 transmission
        SRAP  = sqrt(PATH(j,7))                    
        TCO22 = ( PADE(1,7) + SRAP*PADE(2,7) ) / 
     &              ( PADE(1,7)+ SRAP * ( PADE(3,7) + SRAP ) )
C
        SRAP  = sqrt(PATH(j,8))                                         

        TCO234 = ( PADE(1,8) + SRAP*PADE(2,8) ) /
     &              ( PADE(1,8) + SRAP * ( PADE(3,8) + SRAP ) )
C
        FBWCP1 = EXPPATH(J,1)  

        FBWCP2 = FBWCP1 * FBWCP1
        FBWCP4 = FBWCP2 * FBWCP2
        FBWCP8 = FBWCP4 * FBWCP4
        FBWC16 = FBWCP8 * FBWCP8
        FBWC32 = FBWC16 * FBWC16

        FBWCPQ = EXPPATH(J,2)                                           
C
        FBWCB5 = FBWC32 * FBWC32 * FBWC16                               
        FBWCB2 = FBWC32 * FBWCB5
        FBWCB3 = FBWCP4 * FBWCP2 * FBWCPQ
        FBWCB4 = FBWCP4 * FBWCP1
C
        SBWCP1 = EXPPATH(J,3)                                           
C
        SBWCP2 = SBWCP1 * SBWCP1
        SBWCP4 = SBWCP2 * SBWCP2
        SBWCP8 = SBWCP4 * SBWCP4
        SBWC16 = SBWCP8 * SBWCP8
        SBWC32 = SBWC16 * SBWC16

        SBWCPQ = EXPPATH(J,4)                                           
C
        SBWCB2 = SBWCP8 * SBWCP4
        SBWCB3 = SBWCP4 * SBWCP2 * SBWCPQ
        SBWCB4 = SBWCP4 * SBWCP1
        SBWCB5 = SBWC32 * SBWC32 * SBWC16
C
        TRANS(J,2) = TRANS(J,2) * FBWCB2 * SBWCB2 * TCO22               
        TRANS(J,3) = TRANS(J,3) * FBWCB3 * SBWCB3 * TCO234              
        TRANS(J,4) = TRANS(J,4) * FBWCB4 * SBWCB4 * TCO234        
     1   * ( O3WT1 * exp_vxy(j) +  O3WT2 * exp_vxy(j+l))
        TRANS(J,5) = TRANS(J,5) * FBWCB5 * SBWCB5
C
      ENDDO
C
      RETURN
      END
      SUBROUTINE LWLKIN (DUMMY)
CLL   ! Dummy routine provided in version 1B of the LW from release 2.7
CLL   !   of the UM for compatibility with version 1A and the
CLL   !   control-level routines.
      REAL DUMMY                       ! A dummy argument in all senses.
      RETURN
      END
