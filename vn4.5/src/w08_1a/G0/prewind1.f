C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! subroutine prewind (for wave model)
!
! Description:
!   input values from wave model dump of wind components
!   are converted to speed and direction, and are then rearranged
!   into the blocked data arrays required by WAM
!   the call to airsea will calculate values for ustar (array us_out)
!   and also roughness length Z0 (array Z0_out)
!
! Method:
!
!
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! UM4.1    June 1996 Original code.  M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header

      SUBROUTINE PREWIND(UWIND,VWIND,TAUW_X,TAUW_Y,mdata,
     & U10_OUT,
     & THW_OUT,
     & US_OUT,
     & Z0_OUT,
     & TAUW,

C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,


     & icode,cmessage)


C*    *PARAMETER* OF GLOBAL CONSTANTS.
C
      PARAMETER (G = 9.806, PI = 3.14159265358978, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *G*         REAL      ACCELLERATION OF GRAVITY.
C      *PI*        REAL      PI.
C      *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
C      *RAD*       REAL      PI / 180.
C      *DEG*       REAL      180. / PI.
C      *ZPI*       REAL      2. * PI.
C      *R*         REAL      EARTH RADIUS        (METRES).
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C     ! table dimensions !
      INTEGER    ITAUMAX, JUMAX, IUSTAR, IALPHA
      PARAMETER (ITAUMAX=100, JUMAX=100, IUSTAR=100, IALPHA=100)
C

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
      real DELPHI  ! grid increment for latitude (in metres)
      real DELLAM  ! grid increment for long. at equator (in metres)
      real SINPH(NGY), COSPH(NGY) ! sin / cos of latitude

      integer IGL   ! number of blocks
      integer IJS(NBLO)  ! index of first point of second lat
      integer IJL2(NBLO) ! index of last  point of second lat
      integer IJLS(NBLO) ! index of first point of lat before last
      integer IJL(NBLO)  ! index of last  point of lat before last
      integer IJLT(NBLO) ! total number of gridpoints in a block
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
       INTEGER
     & NANG,       ! number of direction components
     & NFRE,       ! number of frequency components
     & NGX,        ! number of cols in LS mask grid
     & NGY,        ! number of rows in LS mask grid
     & NBLO,       ! max number of blocks
     & NIBLO,      ! max number datapoints per block
     & NOVER,      ! max number datapoints in overlap row
     & NIBLD, NBLD, NIBLC, NBLC
C
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
      real BETAMAX      ! parameter for wind input
      real ZALP         ! shifts growth curve
      real ALPHA        ! charnock constant
      real XKAPPA       ! von karman constant
      real XNLEV        ! assumed height of input winds
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
      real TAUT(0:ITAUMAX,0:JUMAX)   ! stress table
      real DELTAUW                   ! wave stress increment
      real DELU                      ! wind increment
      real TAUHFT(0:IUSTAR,0:IALPHA) ! high freq. stress table
      real DELUST                    ! ustar increment
      real DELALP                    ! alpha increment


      INTEGER ICODE            ! OUT return code
      CHARACTER*80 CMESSAGE    ! OUT message accompanying return code

      integer mdata

      real
     &     UWIND(mdata), ! IN u component of wind on sea-points grid
     &     VWIND(mdata), ! IN v component of wind on sea-points grid
     &     TAUW_X(MDATA),! wave stress tauw on sea points grid x comp
     &     TAUW_Y(MDATA) ! wave stress tauw on sea points grid y comp


C     ARRAYS OUT

      real
     & U10_OUT(niblo,nblo), ! blocked wind speeds
     & THW_OUT(niblo,nblo), ! blocked wind direction
     & US_OUT(niblo,nblo),  ! blocked ustar
     & Z0_OUT(niblo,nblo),  ! blocked roughness length
     & TAUW(niblo,nblo)     ! blocked wave induced stress

c     local arrays:

      real windspeed(mdata), windir(mdata), tw_fld(mdata)

      integer
     &       i    ! loop counter
     &      ,ig   ! loop counter over blocks


CLL 1.0 convert UWIND and VWIND into speed and direction


      do i=1,mdata

       windspeed(i)=sqrt(uwind(i)*uwind(i) + vwind(i)*vwind(i) )

       tw_fld(i)=sqrt(tauw_x(i)*tauw_x(i) + tauw_y(i)*tauw_y(i) )

       windir(i) = atan2(uwind(i),vwind(i))

      enddo


C   NOW call blokfld for windspeed wind direction and wave stress

      call blokfld(windspeed,U10_OUT,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

      if(icode.eq.1)then
       CMESSAGE='error in PREWIND: BLOKFLD U10'
       return
      endif

      call blokfld(windir,THW_OUT,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

      if(icode.eq.1)then
       CMESSAGE='error in PREWIND: BLOKFLD THW'
       return
      endif

      call blokfld(tw_fld,tauw,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)
C
      if(icode.eq.1)then
       CMESSAGE='error in WAVSTEP: BLOKFLD TAUW'
       return
      endif


C     now call AIRSEA for each block in turn NOTE IGL is in ARGWVGD

      DO IG = 1,IGL

C      this call is for u10 to convert to ustar

         CALL AIRSEA (U10_OUT(1,ig), TAUW(1,IG), US_OUT(1,ig),
     & Z0_OUT(1,ig), IJS(ig), IJL(ig),
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

     & icode)

      if(icode.eq.1)then
       CMESSAGE='error in PREWIND: AIRSEA '
       return
      endif
C
      ENDDO

      RETURN
      END
