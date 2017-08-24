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
CLL Subroutines SWTRAN, SWLKIN ---------------------------------------
CLL
CLL Purpose :
CLL  It calculates shortwave transmissivities for each band for a
CLL  particular set of pathlengths of the absorbing gases.
CLL
CLL  Before SWTRAN is CALLed (normally via SWRAD and SWMAST), SWLKIN
CLL  must be CALLed to initialize LUT, which contains the transmissivity
CLL  and difference-of-transmissivity look-up tables.  (Some other
CLL  numbers used to access them are set by SWMAST and passed in in
CLL  TTEC: a single array is used for what are logically 3 types of
CLL  quantity to reduce CALLing overheads.)
CLL
CLL  It is intended to be easily modified to perform also some of the
CLL  functions of D23 (radiation diagnostics).
CLL  Suitable for single column model use.
CLL
CLL    Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                    CRAY HF functions replaced by T3E fast vector 
CLL                    functions.     S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms to standard A of UMDP 4 (version 3, 07/9/90), and
CLL  includes no features deprecated in 8X.
CLL  If *DEF CRAY is off, the code is standard FORTRAN 77 except for
CLL  having ! comments (it then sets the "vector length" to be 1) but
CLL  otherwise it includes CRAY automatic arrays also.
CLL
CLL Logical components covered : P234
CLL  (interaction of shortwave radiation with the atmosphere)
CLL
CLL Project task :
CLL
CLL External documentation:
CLL  Offline documentation in UMDP 23, particularly Appendix 2.
CLL
CLLEND ---------------------s-------------------------------------------
C*L
      SUBROUTINE SWTRAN (PATH, TTEC, TRTAB, DTRTAB,
     &     L,                                                           
     &     TRANS)
C*
      INTEGER NGASES
C Number of absorbing gases treated in the shortwave
      PARAMETER (NGASES=3)    !  Standard set is water vapour, ozone
C                             !  and carbon dioxide.
      INTEGER NBANDS         ! Number of spectral bands in the shortwave
      PARAMETER (NBANDS=4)   ! This run uses the standard set of
C shortwave bands as described by Slingo (May 1989, J.Atmos.Sci., 46,
C 10, 1419-1427) or UMDP 23.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER NLKUPS          ! DIMENSION OF LOOK-UP TABLES
      PARAMETER (NLKUPS=50)
      INTEGER!, INTENT (IN)
     &     L                           !  Number of points
C*L
      REAL!, INTENT (IN)
     &     PATH(L,NGASES),             !  Total pathlength for each gas
     &     TRTAB(NLKUPS,NTRANS,NGASES),
C     !  Look-up tables for transmissivities for each gas and of
C     ! differences of their successive elements.
     &     DTRTAB(NLKUPS,NTRANS,NGASES),
     &     TTEC(NGASES,NTRANS+2)
      REAL!, INTENT (OUT)
     &     TRANS(L,NBANDS)             !  Transmissivities in each band
C     !  Note that the transmissivities are the fraction of the total
C     !  incoming solar to be transmitted in each band, i.e. the
C     !  transmissivity for the band alone multiplied by the fraction
C     !  of the solar constant in the band.  This simplifies SWMAST.
C     !
C     ! SWTRAN has 2 dynamically allocated workspace arrays, no EXTERNAL
C     ! calls and no significant structure - just nested loops.
C*
      REAL FSCIEB(NBANDS)     ! fraction of solar constant in each band
      DATA FSCIEB / .429760, .326158, .180608, .033474 /
C     !  First value has 3% knocked off for Rayleigh scattering.
      INTEGER FSTBAND(NGASES),         !  First and last band in which
     &     LSTBAND(NGASES)             !               each gas absorbs
      DATA FSTBAND / 2, 1, 3 /, LSTBAND / 4, 1, 4 /
C                                                  
      INTEGER I(L)              !  Index for which element of the
C                               !  transmissivity look-up table is used
      REAL TR1GAS,              !  Transmissivity considering only 1 gas
     &     Y(L)                 !  Scaled log(pathlength): its integer
C     !  part is I and its fractional part gives the fraction to move
C     !  towards the next entry in the look-up table.
      INTEGER BAND, GAS, J ,K   !  Loopers over bands, gases and points 
! Local workspace
      REAL LOGPATH(L,NGASES)
C
C     !  Initialize TRANS from FSCIEB:
C
      DO 100 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
       DO 100 J=1, L
        TRANS(J,BAND) = FSCIEB(BAND)
  100 CONTINUE
C
C     !  Find and combine the TR1GAS terms:
C
      DO   J=1,L
        DO K=1,NGASES
          logpath(j,k)=log(path(j,k))
        END DO
      END DO 

      DO 1000 GAS=1, NGASES
Cfpp$  Select(CONCUR)
       DO 101 J=1, L
        Y(J) = TTEC(GAS,NTRANS+1)
     &       + TTEC(GAS,NTRANS+2) * LOGPATH(J,GAS)
        I(J) = INT(Y(J))
        Y(J) = Y(J) - REAL(I(J))
C       !  For very large pathlengths, use maximum values in the table:
        I(J) = MIN(I(J),NLKUPS)
  101  CONTINUE
       DO 1000 BAND=FSTBAND(GAS), LSTBAND(GAS)
Cfpp$   Select(CONCUR)
        DO 1000 J=1, L
         IF ( I(J) .GT. 0 ) THEN
C (Equivalent to IF ( PATH(J,GAS) .GT. RMNPTH(GAS) ) but safer.)
            TR1GAS = TRTAB(I(J),BAND,GAS) + Y(J) * DTRTAB(I(J),BAND,GAS)
          ELSE
C           !  For very small pathlengths, absorption goes linearly to 0
            TR1GAS = 1. - PATH(J,GAS) * TTEC(GAS,BAND)
         ENDIF
C        !  We assume random overlap of different gases' absorption
C        !  lines, so that their transmissivities just multiply:
         TRANS(J,BAND) = TRANS(J,BAND) * TR1GAS
 1000 CONTINUE
C
      RETURN
      END
      SUBROUTINE SWLKIN (SWLUT)
      INTEGER NBANDS         ! Number of spectral bands in the shortwave
      PARAMETER (NBANDS=4)   ! This run uses the standard set of
C shortwave bands as described by Slingo (May 1989, J.Atmos.Sci., 46,
C 10, 1419-1427) or UMDP 23.
      INTEGER NGASES
C Number of absorbing gases treated in the shortwave
      PARAMETER (NGASES=3)    !  Standard set is water vapour, ozone
C                             !  and carbon dioxide.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER NLKUPS          ! DIMENSION OF LOOK-UP TABLES
      PARAMETER (NLKUPS=50)
      REAL!, INTENT(OUT)
     &     SWLUT(NLKUPS,NTRANS,NGASES,2)
      REAL TRTAB(NLKUPS,NTRANS,NGASES)
C
      INTEGER JTRANS, GAS, J     ! Loop over transmissivity, gas & ...
CLL   *COMDECK SWLKUPNU contains the basic data for the transmissivity
CLL look-up tables for the shortwave.  These data were obtained by Tony
CLL Slingo and John Edwards in August and September 1991, by fitting to
CLL LOWTRAN 7.  They also checked that the previous scalings still
CLL seemed appropriate.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
C
C     !  The first gas is water vapour.  Data were derived for 50000 Pa
C     !    and 250 K, and will be scaled by the 0.9th power of pressure,
C     !    which fits bands 2, 3 and 4 reasonably well.
      DATA  ((TRTAB(J,JTRANS,1),J=1, NLKUPS), JTRANS=1, NTRANS) /
     &   NLKUPS*0.,
     &   .999925, .999902, .999873, .999834, .999782, .999712, .999627,
     &   .999515, .999373, .999182, .998941, .998628, .998221, .997700,
     &   .997019, .996140, .995004, .993543, .991673, .989243, .986138,
     &   .982166, .977118, .970772, .962696, .952648, .940233, .925080,
     &   .906964, .885235, .860017, .831196, .798890, .763523, .724730,
     &   .683486, .640193, .595601, .550869, .506163, .463130, .422344,
     &   8*.384329,
     &   .999351, .999162, .998924, .998615, .998220, .997710, .997056,
     &   .996219, .995155, .993775, .992024, .989785, .986943, .983364,
     &   .978792, .973070, .965924, .957079, .946315, .933128, .917455,
     &   .899109, .878105, .854785, .829092, .801972, .773946, .745562,
     &   .717345, .688908, .660691, .632481, .604137, .575597, .546026,
     &   .515551, .483807, .450740, .416722, .381512, .346202, .311246,
     &   8*.277160,
     &   .996628, .995673, .994443, .992886, .990863, .988295, .985024,
     &   .980869, .975646, .968990, .960674, .950325, .937567, .922114,
     &   .903300, .881119, .855402, .826302, .794441, .759906, .724076,
     &   .687667, .651233, .615066, .578117, .540298, .500889, .459635,
     &   .416955, .372655, .328404, .285283, .244560, .207373, .173516,
     &   .143549, .117077,.0938588,.0738897,.0568386,.0430195, .0322597,
     &   8*.0242632 /
C     !  The second gas is ozone.  Data were derived for 25000 Pa and
C     !   225 K, but there is little pressure or temperature dependence.
      DATA  ((TRTAB(J,JTRANS,2),J=1, NLKUPS), JTRANS=1, NTRANS) /
     &   .999721, .999649, .999560, .999449, .999310, .999137, .998911,
     &   .998639, .998304, .997905, .997421, .996827, .996140, .995327,
     &   .994397, .993338, .992178, .990903, .989540, .988132, .986673,
     &   .985179, .983627, .982006, .980291, .978461, .976524, .974414,
     &   .972120, .969655, .966930, .963926, .960547, .956720, .952319,
     &   .947216, .941332, .934339, .926017, .916281, .904617, .890829,
     &   .874387, .854918, .831874, .804833, .773879, .738109, .697762,
     &   .654226,
     &   NLKUPS*0.,  NLKUPS*0.,  NLKUPS*0. /
C     !  The third gas is carbon dioxide.  Data were derived for
C     !    25000 Pa and 225 K, and will be scaled by the 0.7th power
C     !    of pressure, which fits both bands 3 and 4 fairly well.
      DATA  ((TRTAB(J,JTRANS,3),J=1, NLKUPS), JTRANS=1, NTRANS) /
     &   NLKUPS*0.,  NLKUPS*0.,
     &   .999997, .999994, .999986, .999969, .999937, .999869, .999734,
     &   .999460, .998918, .997831, .995705, .991615, .984135, .971311,
     &   36*.951543,
     &   .999893, .999783, .999560, .999106, .998194, .996355, .992726,
     &   .985687, .972763, .950973, .920079, .885743, .855990, .831226,
     &   36*.809535 /
C
      DO 1 GAS=1, NGASES
       DO 1 JTRANS=1, NTRANS
        DO 1 J=1, NLKUPS
         SWLUT(J,JTRANS,GAS,1) = TRTAB(J,JTRANS,GAS)
    1 CONTINUE
C
      DO 2 GAS=1, NGASES
       DO 2 JTRANS=1, NTRANS
        DO 2 J=1, NLKUPS - 1
         SWLUT(J,JTRANS,GAS,2) =
     &    TRTAB(J+1,JTRANS,GAS) - TRTAB(J,JTRANS,GAS)
    2 CONTINUE
C
C     ! Set the last element for each gas and band to zero, so that the
C     ! extrapolation done for any pathlength greater than the maximum
C     ! catered for just gives the greatest value in TRTAB.
C
      DO 3 GAS=1, NGASES
       DO 3 JTRANS=1, NTRANS
        SWLUT(NLKUPS,JTRANS,GAS,2) = 0.
    3 CONTINUE
C
      RETURN
      END
