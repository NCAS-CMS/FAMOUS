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

      SUBROUTINE PROPDOT(ishallo, irefra,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,

C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
     & THDD, THDC, SIDC,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

     & icode)

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
C

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
      real U(0:NIBLC,NBLC)  ! u-component of current
      real V(0:NIBLC,NBLC)  ! v component of current
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
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.
C
      integer IXLG(NIBLO,NBLO) ! longitude grid index
      integer KXLT(NIBLO,NBLO) ! latitude grid index
      integer NY, NX           ! number of lats / longs in grid
      integer IPER             ! set to 1 for a periodic grid
      integer ICASE            ! set to 1 for lat-long grid

      real AMOWEP    ! most western longitude (deg)
      real AMOSOP    ! most southern latitude (deg)
      real AMOEAP    ! most eastern longitude (deg)
      real AMONOP    ! most northern latitdue (deg)
      real XDELLA    ! grid increment for latitude (deg)
      real XDELLO    ! grid increment for longitude (deg)
C
C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
      real THDD(NIBLD,NANG,nblo) ! depth gradient part of theta dot
      real THDC(NIBLC,NANG,nblo) ! current gradient part of theta dot
c
CCNEW MH array SIDC to replace use of SL as workspace in WAM
cccc later make nangc / nfrec etc and  either set to nang or 1
      real SIDC(NIBLC,NANG,nfre,nblo)
ccc                        ! current gradient part of SIGMA dot
C
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
      real DEPTH(NIBLO, NBLO)  ! water depth (metres)
      real DEPTHA, DEPTHD      ! min depth and increment for tables (m)
      real TCGOND(NDEPTH,NFRE) ! shallow water group velocity table
      real TFAK(NDEPTH,NFRE)   ! wave number table
      real TSIHKD(NDEPTH,NFRE) ! table for omega /sinh(2kd)

      integer INDEP(NIBLO)     ! depth index for gridpoint :one block

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
      integer KLAT(NIBLO,2,nblo) ! index of gridpoint south and north
      integer KLON(NIBLO,2,nblo) ! index of gridpoint west and east
c                                ! land points marked by zero
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

C ----------------------------------------------------------------------
C
C**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.
C
C     H. GUNTHER   GKSS/ECMWF   17/02/91
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PROPDOT*
C
C     METHOD.
C     -------
C
C       IN A LOOP OVER THE BLOCKS THE COMMON UBUF IS READ,
C       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED,
C       COMMON REFDOT (DEPTH AND CURRENT REFRACTION FOR THETA DOT)
C       IS COMPUTED AND WRITTEN TO MASS STORAGE (IU16).
C       IN CASE OF CURRENT REFRACTION THE COMPLETE SIGMA DOT TERM
C       IS COMPUTED AND WRITTEN TO IU16 ADDITIONALLY.
C       WRITE OPERATIONS ARE NOT DONE FOR COMMON UBUF AND REFDOT
C       IF THIS IS A ONE BLOCK MODEL.
C
C     EXTERNALS.
C     ----------
C
C       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
ccc*CALL PARALL
C
ccc*CALL COMCURR
C
ccc*CALL COMFRED
C
ccc*CALL COMGRID
C
ccc*CALL COMMAP
C
ccc*CALL COMREFD
C
ccc*CALL COMSHAL
C
ccc*CALL COMSOUR
C
ccc*CALL COMSTAT
C
ccc*CALL COMUBUF
C
ccc*CALL COMUNIT
C
C ----------------------------------------------------------------------
C
C       local arrays for one block only
C
      DIMENSION DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC),
     1          DVPHI(NIBLC), DVLAM(NIBLC), DCO(NIBLD), OMDD(NIBLC)
C
C ----------------------------------------------------------------------
C
ccmh notes - array sl from comSR is used as work space only
ccmh for shallow water current refraction need array (ij,k,l)
ccmh for deep water current refraction use array (ij,k)

C*    1. IF CARTESIAN PROPAGATION SET COSINE OF LAT TO 1.
C         -----------------------------------------------
C
      WRITE(6,*)'in propdot'
 1000 CONTINUE
      IF (ICASE.NE.1) THEN
         DO 1001 IJ = 1,NIBLD
            DCO(IJ) = 1.
 1001    CONTINUE
      ENDIF
C
C*    2. LOOP OVER BLOCKS.
C        -----------------
C
      DO 2000 IG = 1,IGL
C
C*    2.1 IF MULTI BLOCK VERSION.
C         -----------------------
C
cccc     IF (IGL.NE.1) THEN ! UMwave needs to do block one also
C
C*    2.1.2 COMPUTE SHALLOW WATER TABLE INDICES. for this block
c      but not for   block 1 ?? where is indep filled for block one??
C           ------------------------------------
C
            IF (ISHALLO.NE.1) THEN
               DO 2121 IJ=1,IJLT(IG)
                if(depth(ij,ig).gt.0.) then
                  XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
                  ID = NINT(XD)
                  ID = MAX(ID,1)
                  INDEP(IJ) = MIN(ID,NDEPTH)
                else
      WRITE(6,*)'fatal error in propags: zero depth encountered'
                 icode=1
                 goto 999
                endif
 2121          CONTINUE
            ENDIF
CSHALLOW
cccc     ENDIF ! commented out as UMwave needs to do block one also
C
C*    2.2 DEPTH AND CURRENT GRADIENTS. for this block / every block
C         ----------------------------
C
        WRITE(6,*)'calling gradi from propdot'
         CALL GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI,
     &               DULAM, DVPHI, DVLAM,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

     & icode)
C
C*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
C         ---------------------------------------------
C
         IF (ICASE.EQ.1) THEN
            DO 2301 IJ = IJS(IG),IJL(IG)
               JH = KXLT(IJ,IG)
               DCO(IJ) = 1./COSPH(JH)
 2301       CONTINUE
         ENDIF
C
C*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
C         ---------------------------------
C
         IF (ISHALLO.NE.1 .AND. IREFRA.EQ.2) THEN
            DO 2401 IJ = IJS(IG),IJL(IG)
               OMDD(IJ) = V(IJ,IG)*DDPHI(IJ) +
     1                    U(IJ,IG)*DDLAM(IJ)*DCO(IJ)
 2401       CONTINUE
         ENDIF
C
C*    2.5. LOOP OVER DIRECTIONS.
C          ---------------------
C
         DO 2501 K=1,NANG
            SD = SINTH(K)
            CD = COSTH(K)
C
C*    2.5.1. DEPTH GRADIENT OF THETA DOT.
C            ----------------------------
C
            IF (ISHALLO.NE.1) THEN
               DO 2511 IJ = IJS(IG),IJL(IG)
                  THDD(IJ,K,ig) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
 2511          CONTINUE
            ENDIF
C
C*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
C           ---------------------------------------------------
C
            IF (IREFRA.EQ.2) THEN
               SS  = SD**2
               SC  = SD*CD
               CC  = CD**2
               DO 2521 IJ = IJS(IG),IJL(IG)
ccc old line      SL(IJ,K,NFRE) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)
              SIDC(IJ,K,NFRE,ig) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)
     1                          - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
                  THDC(IJ,K,ig) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)
     1                          - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
 2521          CONTINUE
C
C*    2.5.3 LOOP OVER FREQUENCIES. if shallow water + currents
C           ----------------------
C
               IF (ISHALLO.NE.1) THEN
               DO 2530 M=1,NFRE
                  DO 2531 IJ=IJS(IG),IJL(IG)
ccc old line         SL(IJ,K,M) = (SL(IJ,K,NFRE)*TCGOND(INDEP(IJ),M)
           SIDC(IJ,K,M,ig) = (SIDC(IJ,K,NFRE,ig)*TCGOND(INDEP(IJ),M)
     1                          + OMDD(IJ)*TSIHKD(INDEP(IJ),M))
     2                          * TFAK(INDEP(IJ),M)
 2531             CONTINUE
C
C*    BRANCH BACK TO 2.5.3 FOR NEXT FREQUENCY.
C
 2530          CONTINUE
               ENDIF
            ENDIF
C
C*    BRANCH BACK TO 2.5 FOR NEXT DIRECTION.
C
 2501    CONTINUE
C
C*    BRANCH BACK TO 2. FOR NEXT BLOCK.
C
 2000 CONTINUE

  999 continue
      RETURN
      END
