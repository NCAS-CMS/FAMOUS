*IDENT CALC_SOLAR
*/
*/  calculate secular variations in orbital parameters for 
*/  paleo runs based on the routines of Berger78 (JAS 35).
*/  basically a backport of SOLPOS and ORBPRM from UM6.1
*/
*COMDECK CSOLAR
c
c common block to hold the orbital parameters so that they 
c can be read in, and then only have to be redone once a year
c
      REAL    SC,
     &        GAMMA,E,TAU0,SINOBL,
     &        E1, E2, E3, E4, DINY,
     &        SEC_VAR_FACTOR
      LOGICAL L_SEC_VAR,L_SEC_VAR_ONLINE,L_SEC_VAR_FILE
      CHARACTER*80 SEC_VAR_FILE
      INTEGER SEC_VAR_YEAR

      COMMON /COMSOLAR/ SC,
     &                  GAMMA,E,TAU0,SINOBL,
     &                  E1, E2, E3, E4, DINY,
     &                  L_SEC_VAR,L_SEC_VAR_ONLINE,
     &                  L_SEC_VAR_FILE,SEC_VAR_FILE,
     &                  SEC_VAR_YEAR, SEC_VAR_FACTOR

      NAMELIST /NLSTSOLAR/ SC,SEC_VAR_YEAR,
     &               L_SEC_VAR,L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &               SEC_VAR_FILE,SEC_VAR_FACTOR
*/
*DECLARE SWSC
*/  replace declaration of SC
*D SWSC.2,SWSC.3
*CALL CSOLAR
*/
*DECLARE READLSA1
*/ read in namelist values for orbital params, set up param. values
*/ for the start of the run
*/
*B READLSA1.34
*CALL CSOLAR
      REAL OBLQ                            ! Obliquity (local)
*B READLSA1.145
c set defaults for solar param namelist
      SC=1365.
      SEC_VAR_YEAR=2000
      L_SEC_VAR=.FALSE.
      L_SEC_VAR_ONLINE=.FALSE.
      L_SEC_VAR_FILE=.FALSE.
      SEC_VAR_FILE="dummy_filename"
      SEC_VAR_FACTOR=1.
c read in any user specified values
      READ(5,NLSTSOLAR)

c if using secular variations, take account of the acceleration factor
      if (L_SEC_VAR) SEC_VAR_YEAR=
     &   ( A_FIXHD(28)-MODEL_BASIS_TIME(1) )*SEC_VAR_FACTOR  + 
     &     MODEL_BASIS_TIME(1)

c initialise orbital params, whether fixed or varying
      CALL ORBPRM(L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &            SEC_VAR_FILE,SEC_VAR_YEAR,LCAL360,
     &            E,GAMMA,OBLQ,TAU0,DINY)

      SINOBL=SIN(OBLQ)

      write(6,*)"solar constants: year,varying?",SEC_VAR_YEAR,L_SEC_VAR
      write(6,*)"               : table of variations?",L_SEC_VAR_FILE
      write(6,*)"               : online variations?",L_SEC_VAR_ONLINE
      write(6,'(a,f8.1,f10.6,f10.6,f10.5,f10.6)')
     &          "                : values",SC,GAMMA,E,TAU0,SINOBL

      E1 = E * (2.-.25*E*E)
      E2 = 1.25 * E*E                ! Coefficients for 3.1.2        
      E3 = E*E*E * 13./12.
      E4=( (1.+E*E*.5)/(1.-E*E) )**2 ! Constant for 3.1.4         
*/
*DECLARE SOLPOS1A
*/ get rid of fixed parameter versions of gamma, e, tau and sinobl
*/ if we want varying parameters and it's a new year, re-call the
*/ orbprm routine. Use the 6.1 calculation to get new SINDEC/SCS
*/ 
*D GSS1F304.678
      SUBROUTINE SOLPOS (DAY, YEAR, SINDEC, SCS, LCAL360_IN)
*D GSS1F304.681        
      LOGICAL LCAL360_IN    !In, true if 360 day calendar in use.
*D SOLPOS1A.43,SOLPOS1A.79
      REAL M, V                            ! Mean & true anomaly
      REAL OBLQ                            ! Obliquity (local)
*CALL C_PI                                                    
      REAL TWOPI
      PARAMETER ( TWOPI = 2. * PI ) 
*CALL CSUBMODL
*CALL CMAXSIZE
*CALL CHSUNITS
*CALL CNTLALL
*CALL CTIME
*CALL CSOLAR

C============================================================
C
c get new versions of the orbital params if needed
c

      IF (L_SEC_VAR .AND. PREVIOUS_TIME(1).lt.I_YEAR) THEN

c work out "orbital year" using the acceleration factor
        SEC_VAR_YEAR= (I_YEAR-MODEL_BASIS_TIME(1))*SEC_VAR_FACTOR  + 
     &            MODEL_BASIS_TIME(1)

        CALL ORBPRM(L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &              SEC_VAR_FILE,SEC_VAR_YEAR,LCAL360,
     &              E,GAMMA,OBLQ,TAU0,DINY)

        SINOBL=SIN(OBLQ)

        write(6,'(a,i8,f10.6,f10.6,f10.5,f10.6)')
     &       "new solar constants for",SEC_VAR_YEAR, GAMMA,E,TAU0,SINOBL

        E1 = E * (2.-.25*E*E)
        E2 = 1.25 * E*E                ! Coefficients for 3.1.2        
        E3 = E*E*E * 13./12.
        E4=( (1.+E*E*.5)/(1.-E*E) )**2 ! Constant for 3.1.4         


      END IF

c     Calculate the mean anomaly at 12Z on the current day.
c     The 0.5 accounts for the time in days since mid-night.            
c     The references are to Smart 1944 (and UMDP23)                     
c     Eq 67 p. 113 and n=2pi/orbital period     (Eq 3.1.1)        
      IF (LCAL360) THEN
        M = (TWOPI / DINY)           * (FLOAT(DAY) - TAU0 - .5)
      ELSE
        M = (TWOPI / 365.2424)       * (FLOAT(DAY) - TAU0 - .5)
      END IF

c       True anomaly, equation 87 in Smart on p. 120 (UMDP23 Eq 3.1.2)  
      V  = M + E1*SIN(M) + E2*SIN(2.*M) + E3*SIN(3.*M)

c       Solar constant scaling factor (UMDP23 Eq 3.1.4)                 
      SCS = E4 * ( 1. + E * COS(V) ) **2

c       sin(solar declination) (UMDP23 Eq 3.1.5)                        
c       The solar declination is related to                             
c        the true longitude of the earth (lambda) by:                   
c        sindec = sin(obliquity) * sin(lambda)                          
c       Lambda is counted counterclockwise from the vernal equinox      
c        and is related to v (the true anomaly) through                 
c        lambda = v + (longitude of perihelion)                         

      SINDEC = SINOBL * SIN (V - GAMMA)
C
*DECK ORBPRM
c                                                                         
c+ Subroutine to calculate the parameters of the Earth's orbit.           
c                                                                         
c (slightly adapted to fit in FAMOUS UM4.5 solpos by r.s.smith 
c  can get E,OBLQ and LPH from precalculated table if req'd
c  NB: mod(TAU0,DINY), mod(GAMMA,TWOPI) loop added at end
c  to match Julia/Michel's offline lookup numbers. I think is
c  OK. The date of the vernal equinox is fixed for the use of
c  LCAL360 - it's more useful in paleoruns to have this as a
c  fixed point rather than use the "real" Gregorian DOY for it)
c
c  Purpose:                                                               
c  This routine returns the parameters of the Earth's orbit               
c   (the eccentricity, obliquity, supplement of the longitude of          
c    perihelion) and the time of the perihelion passage in days
c                                                                          
c  Method:                                                                 
c  For long runs there may be an interest in running with secular          
c   variations in the astronomy. The orbital constants have                
c   been derived from A. L. Berger 1978, J. Atm. Sci, Volume 35            
c   2362-2367. A copy of which can be found in the Met Office library.     
c  For short current runs, or long control runs it is preferrable          
c   not to allow the astronomy to vary, so fixed values are used.          
c                                                                          
c Current Owner of Code: J. M. Edwards                                     
c                                                                          
c History:                                                                 
c       Version         Date                    Comment                    
c       5.2             15/11/00                Original Code              
c                                               E. Ostrom                  
c                                                                          
c Description of Code:                                                     
c   FORTRAN90 complying with UMDP3 Version 7.2 from 5/2/98                 
c                                                                          
c- ---------------------------------------------------------------------   

      SUBROUTINE ORBPRM(L_SEC_VAR_ONLINE, L_SEC_VAR_FILE
     &                , SEC_VAR_FILE, YEAR, LCAL360
     &                , E, GAMMA, OBLQ, TAU0, DINY)                       
                                                                           
      IMPLICIT NONE                                                        
                                                                           
      INTEGER  YEAR       ! Calendar year                   
      LOGICAL  L_SEC_VAR_ONLINE ! Use a table of precalulated values
      LOGICAL  L_SEC_VAR_FILE   ! Calculate orbital values online

      LOGICAL  LCAL360    ! Use a calendar of 360 days      

      CHARACTER*80 SEC_VAR_FILE ! Location of orbital lookup table
                                                                           
c     Parameters of the Earth's orbit:                                     
c                                                                          
      REAL  E             ! Eccentricity of the orbit       
      REAL  GAMMA         ! Supplement of the longitude     
c                                        !  of the perihelion              
      REAL  OBLQ          ! Obliquity of the orbit          
      REAL  TAU0          ! Time of the perihelion          
c                                        !  passage in days                
      REAL  DINY          ! Length of the calendar year     
c                                        !  (in whole days)                
                                                                           
c     Local Variables for use within ORBPRM                                
c                                                                          
      REAL  YEAR_OFFSET                ! Offset of the year from the     
c                                        !  reference year when default    
c                                        !  values apply                   
      REAL  ECN_SN                     ! Eccentricity multiplied by      
c                                        !  the sine of the longitude      
c                                        !  of the perihelion              
      REAL  ECN_CN                     ! Eccentricity multiplied by      
c                                        !  the cosine of the longitude    
c                                        !  of the perihelion              
      REAL  LPH_FIXED_VE               ! Longitude of the perihelion     
c                                        !  relative to a fixed vernal     
c                                        !  equinox                        
      REAL  GN_PRCS                    ! General precession              
      REAL  DATE_VE                    ! Date of the vernal equinox      
c                                        !  in days into the year          
      REAL  NO_LEAP_DAYS               ! The number of leap days,        
c                                        !  used to calculate DATE_VE      
      REAL  MEAN_ANOM_VE               ! Mean anomaly at the vernal      
c                                        !  equinox                        
                                                                           
c     Synthetic constants                                                  
c                                                                          
      REAL  BETA                                                         
      REAL  EE1                                                          
      REAL  EE2                                                          
      REAL  EE3                                                          
                                                                           
      INTEGER  I                       ! Loop variable                   
                                                                           
c     Mathematical constants:                                              
*CALL C_PI
      REAL TWOPI 
      PARAMETER ( TWOPI = 2. * PI )
                                                                           
      REAL TropYearLength
      PARAMETER ( TropYearLength=365.2424)

c     Astronomical Parameters:                                             
c     Default values of the orbital elements
c      (currently those for the epoch J2000 which is 1.5d Jan. 2000):
c     The Eccentricity and Longitue of perhelion are recommended by NAS
c      see (http://ssd.jpl.nasa.gov/elem_planets.html)                     
c     The Obliquity value comes from the Astronomical Almanac for 1984
c      page S26 and is used on several webpages e.g.
c      nedwww.ipac.caltech.edu/help/calc_doc.txt
c      www.stargazing.net/kepler/astrovba2.html
c      http://edhs1.gsfc.nasa.gov/waisdata/docsw/txt/tp4450505.txt         
c
c     The data in the series expansions are adapted from Berger 1978.
c      Andre' L. Berger Journ. of Atm. Sci. Volume 35 p. 2362-2367,
c      and is available from the Met Office library.
c
c     ! Eccentricity of the orbit
      Real E_DFLT
      Parameter ( E_DFLT         = 1.6710222E-02)
c
c     ! Longitude of the perihelion in radians
      Real LPH_DFLT
      Parameter ( LPH_DFLT       = 102.94719*PI/180.0)
c
c     ! Obliquity of the orbit - corresponds to 23.43929111 degrees
      Real OBLQ_DFLT
      Parameter ( OBLQ_DFLT      = 0.409092804)
c
c     ! Reference year for setting the date of the vernal equinox
      Integer YEAR_REF_VE
      Parameter ( YEAR_REF_VE    = 2000)
c
c     ! Date of the vernal equinox in days after the start of the year
c     !  This date is for the year 2000.
      Real DATE_VE_DFLT
      Parameter ( DATE_VE_DFLT   = 79.3159)
c
c     The final parameter required is the time of the perihelion
c     passage, TAU0. For a pure Keplerian orbit, with a specified
c     eccentricity and longitude of the perihelion, this can be
c     deduced from the date of the vernal equinox (as is specified
c     in AMIP-2, for example). In practice it is somewhat more
c     complicated to calculate the time of the perihelion.
c     For simplicity, a mean value for the years 1995-2005 is used
c     here: note that the range of TAU0 in this period is from      
c     1.0 to 3.75 and that there is no simple relationship with leap
c     years.
c                                             
c     ! Time of the perihelion passage in days     
      Real TAU0_DFLT
      Parameter ( TAU0_DFLT      = 2.667)
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
c     The parameters used to calculate secular variations of the orbital
c     elements are taken from A. L. Berger (1978), J. Atm. Sci., vol 35,
c     p. 2362.These have been converted so that:
c     amplitudes are in radians,                  
c     angular frequencies in radians per year and
c     phases in radians                                      
c     Phases have also been converted to be taken relative to
c     J2000 for consistency with the time for the default values of the
c     orbital parameters.
c     Berger's numbers (with time correction) differ slightly from the
c     default values above.
c
c     The obliquity and longitude of the perihelion have been adjusted
c     to agree with the NASA values for J2000, but adjustment of the
c     eccentricity to NASA values is not so easy and has not been done.
c
c
c     ! Reference year     YEAR_REF
      Integer  YEAR_REF
      Parameter  ( YEAR_REF       = 2000)
c
c   -----------------------------------------------------------------
c     Obliquity (Table 1) from the first 24 terms:
c     (enough for deviations of less than 0.002 degrees)
c
c     ! Constant term in the obliquity: from the Astr. Almanac for 1984
c     !  The following value corresponds to 23.320870 degrees at J2000
      Real OBLQ_CNST
      Parameter ( OBLQ_CNST      = 0.40702597)
c
c     ! Number of terms retained in the series for the obliquity
      Integer N_TERM_OBQ 
      Parameter (  N_TERM_OBQ     = 24)
c
c     ! Amplitude
      Real  A(N_TERM_OBQ)
c     ! Angular frequency                                                  
      Real  F(N_TERM_OBQ)                                                
c     ! Phase in the series                                                
      Real  D(N_TERM_OBQ)                                                
c                                                                          
c   -----------------------------------------------------------------      
c     Eccentricity and longitude of the fixed perihelion (Table 4):        
c                                                                          
c     ! Number of terms retained in the series for the                     
c     !  eccentricty and longitude of the perihelion                       
      Integer N_TERM_ECN_LPH
      Parameter  ( N_TERM_ECN_LPH = 19  )                          
c                                                                          
c     ! Amplitude                                                          
      Real  M(N_TERM_ECN_LPH)                                            
c     ! Angular frequency                                                  
      Real  G(N_TERM_ECN_LPH)                                            
c     ! Phase in the series                                                
      Real  B(N_TERM_ECN_LPH)                                            

c ---------------------------------------------------------------
c   some variables for using a lookup table
      INTEGER SEC_VARUNIT,ICODE
      CHARACTER DUMMY
      REAL TABLEYEAR,TABLEYEAR_REF
      REAL FILEYEAR1,FILEYEAR2,FILEYEAR3
      REAL FE1,FOBLQ1,FLPH1
      REAL FE2,FOBLQ2,FLPH2,FLPH3
      REAL LPH_MOVING_VE
      REAL DY1,DY2
c                                                                          
c   ------------------------------------------------------------------72   
c     General Precession (Table 5):                                        
c                                                                          
c     ! Linear rate of precession!                                         
c     ! The value corresponds to 50.439273 seconds per year -Berger 1979   
      Real LIN_RATE_GN_PRCS
      Parameter ( LIN_RATE_GN_PRCS  = 2.44536496E-04  )
c                                                                          
c     ! Constant offset to general precession (in seconds pre year),       
c     ! corrected for 50 years difference in reference time.               
      Real GN_PRCS_CNST
      Parameter ( GN_PRCS_CNST      = 7.14372244E-02  )
c                                                                          
c     ! Number of terms kept in the series for the general precession      
      Integer N_TERM_GN_PRCS
      Parameter ( N_TERM_GN_PRCS = 10  )
c                                                                          
c     ! Amplitude                                                          
      Real  C(N_TERM_GN_PRCS)                                            
c     ! Angular frequency                                                  
      Real  H(N_TERM_GN_PRCS)                                            
c     ! Phase in the series                                                
      Real  R(N_TERM_GN_PRCS)                                            
c                                                                          
c   -----------------------------------------------------------------      
c    Table 1                                                               
c                                                                          
      DATA A/                                                           &  
     &    -1.19372E-02, -4.15640E-03, -3.05103E-03, -2.00849E-03        &  
     &  , -1.51146E-03,  1.49778E-03, -7.88065E-04, -5.62917E-04        &  
     &  ,  4.90244E-04, -3.28170E-04,  1.20767E-04,  1.09471E-04        &  
     &  , -1.02587E-04, -7.58733E-05,  7.46128E-05,  7.11222E-05        &  
     &  , -5.68686E-05,  4.97904E-05,  3.14644E-05,  2.83616E-05        &  
     &  , -2.66163E-05, -2.63254E-05,  2.50164E-05,  2.46285E-05/          
      DATA F/                                                           &  
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04                &  
     &  ,  1.5506174E-04,  2.1733392E-04,  1.5016256E-04                &  
     &  ,  2.1170962E-04,  1.5633636E-04,  1.4835028E-04                &  
     &  ,  2.0692488E-04,  2.1252514E-04,  2.2999289E-04                &  
     &  ,  3.0649899E-04,  3.1139817E-04,  4.8991877E-06                &  
     &  ,  3.6059331E-05,  2.7043965E-04,  1.8122966E-06                &  
     &  ,  6.4084427E-05,  3.0341210E-04,  3.0831127E-04                &  
     &  ,  3.7058338E-04,  2.2211866E-04,  4.0958519E-05/                  
      DATA D/                                                           &  
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  5.1167E+00            &  
     &  ,  2.7912E-01,  4.6115E+00,  5.3935E+00,  4.1966E+00            &  
     &  ,  3.8990E+00,  4.7014E+00,  5.5397E+00,  5.5896E+00            &  
     &  ,  2.5251E+00,  3.0303E+00,  5.0517E-01,  2.1589E+00            &  
     &  ,  3.6608E-01,  7.1253E-01,  2.1582E+00,  2.7325E+00            &  
     &  ,  3.2376E+00,  4.6833E+00,  9.7121E-01,  2.6640E+00/              
c                                                                          
c   -----------------------------------------------------------------      
c    Table 4                                                               
c                                                                          
      DATA M/                                                           &  
     &     1.8607980E-02,  1.6275220E-02, -1.3006600E-02                &  
     &  ,  9.8882900E-03, -3.3670000E-03,  3.3307700E-03                &  
     &  , -2.3540000E-03,  1.4001500E-03,  1.0070000E-03                &  
     &  ,  8.5700000E-04,  6.4990000E-04,  5.9900000E-04                &  
     &  ,  3.7800000E-04, -3.3700000E-04,  2.7600000E-04                &  
     &  ,  1.8200000E-04, -1.7400000E-04, -1.2400000E-04                &  
     &  ,  1.2500000E-05/                                                  
      DATA G/                                                           &  
     &     2.0397105E-05,  3.5614854E-05,  8.6574454E-05                &  
     &  ,  8.3487563E-05,  8.1675266E-05,  2.5205846E-05                &  
     &  ,  8.8386751E-05,  1.2710243E-04,  3.0830121E-05                &  
     &  ,  7.8588375E-05,  1.4860417E-05,  8.0400672E-05                &  
     &  ,  8.9661345E-05,  3.0014587E-05,  9.1473642E-05                &  
     &  ,  8.4481533E-05,  2.9990579E-05,  8.9290274E-05                &  
     &  ,  3.2378912E-06/                                                  
      DATA B/                                                           &  
     &     5.0053E-01,  3.3839E+00,  5.3852E+00,  5.5925E+00            &  
     &  ,  4.8800E+00,  1.5230E+00,  6.0977E+00,  2.2481E+00            &  
     &  ,  2.6918E+00,  5.0874E+00,  2.0054E+00,  5.8001E+00            &  
     &  ,  5.1778E+00,  2.5455E+00,  5.8903E+00,  2.6587E+00            &  
     &  ,  2.2151E+00,  3.6812E+00,  1.2585E+00/                           
c                                                                          
c   -----------------------------------------------------------------      
c    Table 5                                                               
c                                                                          
      DATA C/                                                           &  
     &     3.58327E-02,  1.23877E-02,  9.80662E-03, -9.56853E-03        &  
     &  ,  6.01280E-03,  4.62449E-03, -4.51725E-03,  4.22942E-03        &  
     &  ,  2.93967E-03, -2.40482E-03/                                      
      DATA H/                                                           &  
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04,  3.0868911E-06&  
     &  ,  1.5506174E-04,  1.5217749E-05,  1.5016256E-04,  2.1733392E-04&  
     &  ,  4.8087409E-06,  1.8122966E-06/                                  
      DATA R/                                                           &  
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  6.0756E+00            &  
     &  ,  5.1167E+00,  2.8833E+00,  4.6115E+00,  2.7912E-01            &  
     &  ,  1.0225E+00,  7.1253E-01/
                                                                           
c     The length of the calendar year may be set for a 360-day calendar    
c      (as is often used in climate runs),                                 
c      or for a real Gregorian calendar which has 365 days in              
c      non-leap years and 366 in leap years.                               

                                                                           
      IF (LCAL360) THEN                                                    
                                                                           
        DINY=360.0                                                         
                                                                           
      ELSE                                                                 
c      Is this a leap year?                                                
        IF (mod(year,4)   .eq. 0 .AND.                                     
     &     (mod(year,400) .eq. 0 .OR. mod(year,100) .ne. 0)) then          
                                                                           
          DINY = 366.0                                                     
                                                                           c      Is this a normal year?                                              
        ELSE                                                               
                                                                           
          DINY = 365.0                                                     
                                                                           
        END IF                                                             
      END IF                                                               
                                                                           
      IF (L_SEC_VAR_FILE .OR. L_SEC_VAR_ONLINE) THEN
c     The orbital elements are normally set to default values, but         
c     secular variations may be required in some longer climate runs.      
c  
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     r.s.smith
c     Eccentricity, obliquity and the longitude of perihelion can
c     be calculated online for +/- 1Myr via Berger78 formula. For
c     older periods, precalculated tables such as Laskar04 can be used.
c     Assumes a certain format for the tables (note: dates must
c     go back in time, as in Laskar paleo tables):
c      KYEAR FROM
c      <offset year>        ECC   OBLQ     PERH_MOVING_VE
c       ---------------------------------------------------------
c      <kyear from offset> <ecc> <oblq>  <long.perih. rel. to v.e.>
ccccccccccccccccccccccccccccccccccccccccccccccccc

      IF (L_SEC_VAR_FILE) THEN !L_SEC_VAR_FILE, use precalculated values

       SEC_VARUNIT=615

       OPEN(SEC_VARUNIT,FILE=SEC_VAR_FILE,IOSTAT=ICODE)
       REWIND(SEC_VARUNIT)

       READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY
       READ(SEC_VARUNIT,*,IOSTAT=ICODE)TABLEYEAR_REF
       READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY

       TABLEYEAR = REAL( YEAR - TABLEYEAR_REF )                              

       do while (icode.eq.0)
         READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR1,FE1,FOBLQ1,FLPH1
         READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR2,FE2,FOBLQ2,FLPH2

         FILEYEAR1=FILEYEAR1*1000.
         FILEYEAR2=FILEYEAR2*1000.

         if (tableyear.le.FILEYEAR1 .AND. tableyear.ge.FILEYEAR2) 
     &       icode=9999

         BACKSPACE(SEC_VARUNIT)
       end do

       if (ICODE.ne.9999) then
         write(6,*)"orbital lookup: year not found/file error",
     &             TABLEYEAR,icode
         stop
       end if


c linearly interpolate between the years we've got 

         DY2=1-( (TABLEYEAR-FILEYEAR2)/(FILEYEAR1-FILEYEAR2) )
         DY1=1-( (FILEYEAR1-TABLEYEAR)/(FILEYEAR1-FILEYEAR2) )

         E            =   FE1*DY1 + FE2*DY2
         OBLQ         =FOBLQ1*DY1 + FOBLQ2*DY2

c CAREFUL HERE - doing mod(GAMMA,2PI) later goes wrong when LPH
c uses a linear interpolation as above due to LPH wrapping at 0,2PI
c in the table
c LPH decreases as we go back in time. If 2 is greater than 1,
c then we must have wrapped. Use the general gradient to get a
c -ve value to interpolate towards - the wrap on GAMMA later will
c sort the signs out in the end
         if (FLPH2 .gt. FLPH1) then
           READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY
           READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR3,FE2,FOBLQ2,FLPH3
           FILEYEAR3=FILEYEAR3*1000.

           FLPH2=(FLPH1-(FLPH2-FLPH3))*
     &           (FILEYEAR1-FILEYEAR2)/(FILEYEAR2-FILEYEAR3)

         end if

         LPH_MOVING_VE= FLPH1*DY1 + FLPH2*DY2

ccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSE ! L_SEC_VAR_ONLINE, calculate online, not from table

                                                                           
        YEAR_OFFSET = REAL( YEAR - YEAR_REF )                              
        if (ABS(year_offset).GT.1e6) then
          write(6,*) year_offset,
     &  "ONLINE SOLUTION FOR SOLAR PARAMS ONLY VALID FOR +/- 1MYr BP"
          stop
        end if
                                                                           
c       Obliquity: (Equation 1 from Berger 1978)                           
                                                                           
        OBLQ = OBLQ_CNST                                                   
        DO I=1, N_TERM_OBQ                                                 
          OBLQ = OBLQ+A(I)*COS(F(I)*YEAR_OFFSET+D(I))                      
        END DO                                                             
                                                                           
c       Eccentricity: this is better computed from its components          
c       than directly from the series.(Equation (4) of Berger 1978).       
                                                                           
        ECN_SN = M(1) * SIN (G(1) * YEAR_OFFSET + B(1))                    
        ECN_CN = M(1) * COS (G(1) * YEAR_OFFSET + B(1))                    
                                                                           
        DO I=2, N_TERM_ECN_LPH                                             
          ECN_SN = ECN_SN + M(I) * SIN (G(I) * YEAR_OFFSET + B(I))         
          ECN_CN = ECN_CN + M(I) * COS (G(I) * YEAR_OFFSET + B(I))         
        END DO                                                             
        E = SQRT(ECN_SN*ECN_SN+ECN_CN*ECN_CN)                              
                                                                           
c       We now obtain the longitude of the perihelion relative to the      
c       fixed equinox.                                                     
                                                                           
        LPH_FIXED_VE = ATAN2 (ECN_SN,ECN_CN)                               
                                                                           
c       The longitude of perihelion and                                    
c        the supplement of the longitude of the perihelion relative to     
c        the actual vernal equinox requires the general precession.        
                                                                           
c      General Precession.                                                 
        GN_PRCS = LIN_RATE_GN_PRCS * YEAR_OFFSET + GN_PRCS_CNST            
        DO I=1, N_TERM_GN_PRCS                                             
          GN_PRCS = GN_PRCS + C(I) * SIN (H(I) * YEAR_OFFSET + R(I))       
        END DO                                                             
                                                                           
        LPH_MOVING_VE = LPH_FIXED_VE + GN_PRCS

      END IF  !L_SEC_VAR_FILE, lookup or online calc identical from here
ccccccccccccccccccccccccccccccccccccccccccccccccc

c      Supplement of the longitude of the perihelion                       
        GAMMA = PI - LPH_MOVING_VE
                                                                           
c       Time of perihelion: The time at which an object is at perihelion   
c        (its closest distance to the sun).                                
c       The time of perihelion is inferred from the date of                
c        the vernal equinox using the Gregorian calendar.                  
c                                                                          
c      Calculate the date of the vernal equinox.                           
c       First we need to:                                                  
c        Calculate the no of leap days between year & year_ref_ve.         
c        This needs to be corrected when using the Gregorian calendar.     
c         by adding (DINY-366.0) when the year_ref_ve is a leap year or    
c         by adding (DINY-365.0) when the year_ref_ve is a normal year.    
c        This correction is done when the DATE_VE is calculated below!     
c                                                                          
c        In the calculation of NO_LEAP_DAYS below, the divisions of type   
c         'YEAR'/x (where x is 4, 100 or 400) are integer computations.    
c         These integers are then subtracted and the resulting integer     
c         is then converted to a real.                                     
                                                                           
        NO_LEAP_DAYS = ( TropYearLength - 365.0)                           
     &    * REAL( YEAR     - YEAR_REF_VE     )                             
     &    - REAL( YEAR/4   - YEAR_REF_VE/4   )                             
     &    + REAL( YEAR/100 - YEAR_REF_VE/100 )                             
     &    - REAL( YEAR/400 - YEAR_REF_VE/400 )                             
                                                                           
c      Now we can calculate the date of the vernal equinox!                
c      Because the date of the vernal equinox is varying with the year,    
c      we have to keep track of its position in the sky.                   
c      In order to accomodate a time varying vernal equinox when using     
c      a 360-day year, we still have to calculate the difference in        
c      the vernal equinox depending on leap years, normal years and the    
c      difference between the length of the tropical year and the          
c      "normal" year and then we adjust this by multiplying the            
c      DATE_VE by 360/(length of tropical year).                           
c                                                                          
c      Is a 360 day calendar being used?                                   
                                                                           
        IF (LCAL360) THEN                                                  
c         DATE_VE = DATE_VE_DFLT + NO_LEAP_DAYS                            
c         DATE_VE = DATE_VE * DINY / TropYearLength                        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for paleo runs it's far more useful to use the VE as a fixed point
c to define your year's forcing around
          DATE_VE = DATE_VE_DFLT * DINY / TropYearLength                        
                                                                           
                                                                           
c      Is a 365 day calendar being used?                                   
        ELSE                                                               
                                                                           
c        Is the epoch reference year a leap year?                          
                                                                           
          IF (mod(YEAR_REF_VE,4)   .eq. 0 .AND.                            
     &       (mod(YEAR_REF_VE,400) .eq. 0 .OR.                             
     &        mod(YEAR_REF_VE,100) .ne. 0)) THEN                           
                                                                           
            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 366.0))       
                                                                           
c        Is the epoch reference year a normal year?                        
                                                                           
          ELSE                                                             
                                                                           
            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 365.0))       
                                                                           
          END IF                                                           
        END IF                                                             
                                                                           
        BETA = SQRT(1.0E+00-E*E)                                           
        EE1  = (0.5*E + 0.125*E*E*E)*(1.0 + BETA)                          
        EE2  = -0.25*E*E* (0.5 + BETA)                                     
        EE3  = 0.125*E*E*E*((1.0/3.0) + BETA)                              
        MEAN_ANOM_VE = GAMMA - 2.0E+00 * (                                 
     &      EE1 * SIN (GAMMA)                                              
     &    + EE2 * SIN (2.0 * GAMMA)                                        
     &    + EE3 * SIN (3.0 * GAMMA)                                        
     &    )                                                                
                                                                           
        TAU0 = DATE_VE - MEAN_ANOM_VE * TropYearLength/(TWOPI)             

      ELSE ! neither L_SEC_VAR_ONLINE or _FILE, just use default values

        E     = E_DFLT
        OBLQ  = OBLQ_DFLT
        GAMMA = PI - LPH_DFLT
        TAU0  = TAU0_DFLT

      END IF 

                                                                           
c     If using a 360-day calendar the time of the perihelion is            
c     adjusted.                                                            
      IF (LCAL360) THEN                                                    
        TAU0 = TAU0*(360.0/TropYearLength)+0.71                            
      ENDIF                                                               

ccccccccccccccccccccccccccccccccccccccccccccccccc
c r.s.smith - added to get numbers similar to those
c in Julia/Michel's offline lookup table for paleodates
ccccccccccccccccccccccccccccccccccccccccccccccccc
      do while (tau0 .lt.0)
        tau0=tau0+DINY
      end do
      do while (tau0 .gt.DINY)
        tau0=tau0-DINY
      end do
      do while (gamma .lt.0)
        gamma=gamma+2*PI
      end do
      do while (gamma .gt.2*PI)
        gamma=gamma-2*PI
      end do
c                                                                          
      RETURN                                                              
      END                                                                 
