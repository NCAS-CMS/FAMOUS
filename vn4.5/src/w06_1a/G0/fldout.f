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
C*LLL
C     12. SUBROUTINE FLDOUT
C
!
! Description:
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
!
C     DOCUMENTATION
C
C     SEE WAVE MODEL DOCUMENTATION PAPER.
C
C     DESCRIPTION
C
C     THIS ROUTINE OUTPUTS INTEGRATED WAVE PARAMETERS CALCULATED
C     BY CALLING WAVEH. Arrays are prepared for use by STASH
C
C *********************************************************************

      SUBROUTINE FLDOUT(ia1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     +       l_wvtra,len1,energy,uwind,vwind,mdata,ndata,nfldmax
     +  ,icode)



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
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
      real FR(nfre)    ! frequencies (Hz)
      real DFIM(nfre)  ! frequency interval * direction interval
      real GOM(nfre)   ! deep water group velocity
      real C(nfre)     ! deep water phase velocity
      real DELTH       ! angular increment of spectrum (radians)
      real DELTR       ! delth times radius of earth (m)
      real TH(nang)    ! directions in radians
      real COSTH(nang), SINTH(nang)
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

C
      LOGICAL IA1(ngx,ngy)
      REAL WH(ngx*ngy,nfldmax)
      real windsp(ndata),windir(ndata)
      REAL uwind(mdata),vwind(mdata) ! input U10 wind components

      LOGICAL l_wvtra ! IN    T if wavetrains required

      INTEGER ICODE
C
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C        VALUES STORED IN WH                       FIELD NUMBER  PP
C      1       TOTAL WAVE HEIGHT                       94       387
C      2       MEAN DIRECTION                           -       ?
C      3       PRINCIPAL DIRECTION                     96       394
C      4       ZERO UP CROSSING PERIOD                 95       393
C      5       MEAN PERIOD                              -       ?
C      6       PEAK PERIOD                            116       392
C      7       TOTAL WINDSEA WAVE HEIGHT               97       385
C      8       MEAN DIRECTION                           -       ?
C      9       PRINCIPAL DIRECTION                     99       389
C     10       ZERO UP CROSSING PERIOD                 98       388
C     11       MEAN PERIOD                              -       ?
C     12       PEAK PERIOD                              -       ?
C     13       TOTAL SWELL WAVE HEIGHT                100       386
C     14       MEAN DIRECTION                           -       ?
C     15       PRINCIPAL DIRECTION                    102       391
C     16       ZERO UP CROSSING PERIOD                101       390
C     17       MEAN PERIOD                              -       ?
C     18       PEAK PERIOD                              -       ?

C     ??       WAVE STRESS - ABSOLUTE VALUE           211       366
C     ??       ROUGHNESS LENGTH WAVE DEPENDENT        212       367
C
c     19  swh      )  first
c     20  swh      )  second wave train
c     21  swh      )  third
c     22  swh      )  fourth
c
c     23  period   )  first
c     24  period   )  second wave train
c     25  period   )  third
c     26  period   )  fourth
c
c     27  dir      )  first
c     28  dir      )  second wave train
c     29  dir      )  third
c     30  dir      )  fourth
c
c     31  number of wave trains present
C ----------------------------------------------------------------------
C
C still need somewhere to put windsp windir and other data arrays into
C WH for passing to stash eg stress and zo

C
C   need to have all diagnostic arrays on the full grid for STASH
C

C calculate wind speed and direction from components

      do i=1,ndata

       windsp(i) = sqrt(uwind(i)*uwind(i) + vwind(i)*vwind(i) )
       windir(i) = atan2(uwind(i),vwind(i))

      enddo

C
C *** CALCULATE WAVE INTEGRATED FIELDS
C
       ngrid=ngx*ngy

       call waveh(ia1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
     &   l_wvtra,ndata,ngrid,nfldmax,len1,energy,windsp,windir,wh,icode)

      RETURN
      END
