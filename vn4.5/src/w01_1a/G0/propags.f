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

      SUBROUTINE PROPAGS (F1, F3, IG, irefra, ishallo, idelpro,
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
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
     & THDD, THDC, SIDC,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

c* argument list for source term diagnostic arrays for PROPAGS
c len_p2 (source diagnostics) =1 or nang*nfre*niblo as required
c
c THIS set for one block only - pass from WAMODEL to PROPAGS
c
     & sadv2, len_p2,
c
     & icode)

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
C
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

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
      real U(0:NIBLC,NBLC)  ! u-component of current
      real V(0:NIBLC,NBLC)  ! v component of current
C
C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
      integer KLAT(NIBLO,2,nblo) ! index of gridpoint south and north
      integer KLON(NIBLO,2,nblo) ! index of gridpoint west and east
c                                ! land points marked by zero
C
c       local diagnostic arrays passed through argument list
c THIS set for one block only - pass down to PROPAGS
c
       real sadv2 (len_p2)
c
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
C**** *PROPAGS* - COMPUTATION OF A PROPAGATION TIME STEP.
C
C     S.D. HASSELMANN.
C     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER
C
C     MODIFIED BY   H. GUNTHER   01/06/90    -   LAND POINTS ARE TAKEN
C                             OUT OF BLOCKS AND REFRACTION INTEGRATION
C                             CORRECTED FOR N-S AND S-N PROPAGATION.
C
C     K.P. HUBBERT                /07/89    -   DEPTH AND CURRENT
C     S. HASSELMANN   MPIFM       /04/90        REFRACTION SHALLOW
C
C     H. GUNTHER   GKSS/ECMWF   17/01/91    -   MODIFIED FOR CYCLE_4
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF A PROPAGATION TIME STEP.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PROPAGS(F1, F3, IG)*
C          *F1* - SPECTRUM AT TIME T.
C          *F3* - SPECTRUM AT TIME T+DELT.
C          *IG* - BLOCK NUMBER.
C
C     METHOD.
C     -------
C
C       FIRST ORDER FLUX SCHEME.
C
C     EXTERNALS.
C     ----------
C
C       *DOTDC*     - READ DOT TERMS FOR REFRACTION AND SCATTER TABLE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C local arrays
c
      DIMENSION F1(0:NIBLO,NANG,NFRE), F3(0:NIBLO,NANG,NFRE)

      DIMENSION DCO(NIBLO), DP1(NIBLO), DP2(NIBLO)
      DIMENSION DPN(NIBLO), DPS(NIBLO), DPH(0:NIBLO)
      DIMENSION DLE(NIBLO), DLW(NIBLO), DLA(0:NIBLO)
      DIMENSION DOP(NIBLC), DOM(NIBLC)
      DIMENSION DTP(NIBLO), DTM(NIBLO), DRGP(NIBLO), DRGM(NIBLO),
     1          DRDP(NIBLD), DRDM(NIBLD), DRCP(NIBLC), DRCM(NIBLC)
      DIMENSION DTC(NIBLO), CGOND(0:NIBLO)

      real fconst(nibld,nfre)
C
C ----------------------------------------------------------------------
C
c      initialise diagnostics array for this block
c
      do i=1,len_p2
       sadv2(i)=0.
      enddo

C*    0. SPECTRUM AT LAND TO ZERO.
C        -------------------------
C
      DO 1 M=1,NFRE
      DO 1 K=1,NANG
         F3(0,K,M) = 0.
         F1(0,K,M) = 0.
   1  CONTINUE
      CGOND(0) = 0.0
      DPH(0) = 0.0
      DLA(0) = 0.0
C
C*    0.1 READ REFRACTION DOT TERMS.
C         --------------------------
C
c       call to dotdc deleted from here.
c       replaced with arrays held in memory calculated in propdot
CCUM
cc array thdd is already filled for all blocks in setupwv / propdot
cc this following bit of code copied out from dotdc
cc
cc array indep is filled in wamodel for each block
cc and is carried in common argwvsh
cc
cc WAM use of FCONST is a memory saving trick as FCONST is used
cc elsewhere for a different role. The array is redefined in UM wam
cc as a local array in this routine and reference to ARGWVSR removed
cc from propags argument list
ccUM
      if(irefra.ne.0) then
      if(ishallo.ne.1) then
       DO  M=1,NFRE
         DO IJ=IJS(ig),IJL(ig)
            FCONST(IJ,M) = TSIHKD(INDEP(IJ),M)
         ENDDO
       ENDDO
      endif
      endif
C
C*    0.2 SPHERICAL OR CARTESIAN GRID?
C         ----------------------------
C
      IF (ICASE.EQ.1) THEN
C
C*    0.2.1 SPHERICAL GRID.
C           ---------------
C
C*    0.2.1.1 COSINE OF LATITUDE.
C             -------------------
C
         DO 211 IJ = 1,IJLT(IG)
            JH = KXLT(IJ,IG)
            DCO(IJ) = 1./COSPH(JH)
  211    CONTINUE
C
C*    0.2.1.2 COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
C             ------------------------------------------------
C
         DO 212 IJ = IJS(IG),IJL(IG)
            JH = KLAT(IJ,1,ig)
            IF (JH.LE.0) THEN
               DP1(IJ) = 1.
            ELSE
               DP1(IJ) = DCO(IJ)/DCO(JH)
            ENDIF
            JH = KLAT(IJ,2,ig)
            IF (JH.LE.0) THEN
               DP2(IJ) = 1.
            ELSE
               DP2(IJ) = DCO(IJ)/DCO(JH)
            ENDIF
  212    CONTINUE
         IF (IREFRA.NE.2) THEN
C
C*       BRANCH TO 3. IF WITHOUT REFRACTION OR DEPTH.
C        --------------------------------------------
C
            GOTO 3000
         ELSE
C
C*       BRANCH TO 4. IF DEPTH AND CURRENT REFRACTION.
C        ---------------------------------------------
C
            GOTO 4000
         ENDIF
      ELSE
C
C*    0.2.2 CARTESIAN GRID.
C           ---------------
C
C*    0.2.2.1 BRANCH TO 2. IF DEPTH AND CURRENT REFRACTION.
C             ---------------------------------------------
C
         IF (IREFRA.EQ.2) GOTO 2000
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    1. PROPAGATION FOR CARTESIAN GRID
C*       WITHOUT REFRACTION OR DEPTH REFRATION.
C        --------------------------------------
C
 1000 CONTINUE
C
      DELPRO = FLOAT(IDELPRO)
      DELPH0 = DELPRO/DELPHI
      DELLA0 = DELPRO/DELLAM
      DELTH0 = 0.25*DELPRO/DELTH
C
C*    1.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 1100 K=1,NANG
         SD = SINTH(K)*DELLA0
         CD = COSTH(K)*DELPH0
C
C*    1.1.1 INDEX FOR ADJOINING POINTS.
C           ---------------------------
C
         IF (SD.LT.0) THEN
            IJLA = 2
         ELSE
            IJLA = 1
         ENDIF
         IF (CD.LT.0) THEN
            IJPH = 2
         ELSE
            IJPH = 1
         ENDIF
C
         IF (ISHALLO.EQ.1) THEN
C
C*    1.1.2 DEEP WATER.
C           -----------
C
            SD = ABS(SD)
            CD = ABS(CD)
            DTH = SD + CD
C
C*    1.1.2.1 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 1120 M=1,NFRE
C
C*    1.1.2.1.1 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DTT = 1.- DTH*GOM(M)
               DNO = CD*GOM(M)
               DEA = SD*GOM(M)
               DO 1121 IJ = IJS(IG),IJL(IG)
                  F3(IJ,K,M) = DTT * F1(IJ,K,M )
     1                       + DNO * F1(KLAT(IJ,IJPH,ig),K  ,M)
     2                       + DEA * F1(KLON(IJ,IJLA,ig),K  ,M)
 1121          CONTINUE
C
C*    BRANCH BACK TO 1.1.2.1 FOR NEXT FREQUENCY.
C
 1120       CONTINUE
         ELSE
CSHALLOW
C
C*    1.1.3 SHALLOW WATER.
C           --------------
C
            SD = 0.5*SD
            CD = 0.5*CD
C
C*    1.1.3.1 DEPTH REFRACTION.
C             -----------------
C
            IF(IREFRA.EQ.1) THEN
               KP1 = K+1
               IF (KP1.GT.NANG) KP1 = 1
               KM1 = K-1
               IF (KM1.LT.1) KM1 = NANG
               DO 1131 IJ = IJS(IG),IJL(IG)
                  DRDP(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KP1,ig))*DELTH0
                  DRDM(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KM1,ig))*DELTH0
 1131          CONTINUE
            ENDIF
C
C*    1.1.3.2 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 1130 M=1,NFRE
C
C*    1.1.3.2.1 GROUP VELOCITIES.
C               -----------------
C
               CGOND(0) = TCGOND(NDEPTH,M)
               DO 1132 IJ=1,IJLT(IG)
                  CGOND(IJ) = TCGOND(INDEP(IJ),M)
 1132          CONTINUE
C
C*    1.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               IF (SD.GE.0.) THEN
                  DO 1133 IJ=IJS(IG),IJL(IG)
                     DLA(IJ) = SD*(CGOND(KLON(IJ,1,ig)) + CGOND(IJ))
                     DTC(IJ) = SD*(CGOND(KLON(IJ,2,ig)) + CGOND(IJ))
 1133             CONTINUE
               ELSE
                  DO 1134 IJ=IJS(IG),IJL(IG)
                     DLA(IJ) =-SD*(CGOND(KLON(IJ,2,ig)) + CGOND(IJ))
                     DTC(IJ) =-SD*(CGOND(KLON(IJ,1,ig)) + CGOND(IJ))
 1134             CONTINUE
               ENDIF

               IF (CD.GE.0.) THEN
                  DO 1135 IJ=IJS(IG),IJL(IG)
                     DPH(IJ) = CD*(CGOND(KLAT(IJ,1,ig)) + CGOND(IJ))
                     DTC(IJ) = DTC(IJ)
     1                       + CD*(CGOND(KLAT(IJ,2,ig)) + CGOND(IJ))
 1135             CONTINUE
               ELSE
                  DO 1136 IJ=IJS(IG),IJL(IG)
                     DPH(IJ) =-CD*(CGOND(KLAT(IJ,2,ig)) + CGOND(IJ))
                     DTC(IJ) = DTC(IJ)
     1                        -CD*(CGOND(KLAT(IJ,1,ig)) + CGOND(IJ))
 1136             CONTINUE
               ENDIF
               IF (IREFRA.EQ.1) THEN
                  DO 1137 IJ = IJS(IG),IJL(IG)
                     DTHP = FCONST(IJ,M)*DRDP(IJ)
                     DTHM = FCONST(IJ,M)*DRDM(IJ)
                     DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 1137             CONTINUE
               ENDIF
C
C*    1.1.3.2.3 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 1138 IJ = IJS(IG),IJL(IG)
                  F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                       + DPH(IJ) * F1(KLAT(IJ,IJPH,ig),K  ,M)
     2                       + DLA(IJ) * F1(KLON(IJ,IJLA,ig),K  ,M)
 1138          CONTINUE
               IF (IREFRA.EQ.1) THEN
                  DO 1139 IJ = IJS(IG),IJL(IG)
                     F3(IJ,K,M) = F3(IJ,K,M )
     1                          + DTP(IJ) * F1(IJ,KP1,M)
     2                          + DTM(IJ) * F1(IJ,KM1,M)
 1139             CONTINUE
               ENDIF
C
C*    BRANCH BACK TO 1.1.3.2 FOR NEXT FREQUENCY.
C
 1130       CONTINUE
CSHALLOW
         ENDIF
C
C*    BRANCH BACK TO 1.1 FOR NEXT DIRECTION.
C
 1100 CONTINUE
C
C*    1.2 END OF PROPAGATION FOR CARTESIAN GRID
C*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
C         -----------------------------------------------
cc
cc     here extract propagation source term diagnostics:
cc
      if(len_p2.eq.nang*nfre*niblo) then
       WRITE(6,*)'extracting diagnostics Sadv'
       do l=1,nfre
        do m=1,nang
         nstart=((l-1)*nang + m-1)*niblo
         do ip=ijs(ig),ijl(ig)
          sadv2(nstart+ip)=(F3(ip,m,l) - F1(ip,m,l))
         enddo
        enddo
       enddo
      endif


      RETURN
C
C ----------------------------------------------------------------------
C
C*    2. PROPAGATION FOR CARTESIAN GRID
C*       WITH DEPTH AND CURRENT REFRACTION.
C        ----------------------------------
C
 2000 CONTINUE
C
      DELPRO = FLOAT(IDELPRO)
      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELPRO/DELLAM
      DELFR0 = 0.25*DELPRO/(0.1*ZPI)
C
C*    2.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 2100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)*DELLA0
         CD = COSTH(K)*DELPH0
C
C*    2.1.1 DEPTH REFRACTION IF SHALLOW WATER.
C           ----------------------------------
C
         IF (ISHALLO.NE.1) THEN
            DO 2101 IJ = IJS(IG),IJL(IG)
               DRDP(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KP1,ig))*DELTH0
               DRDM(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KM1,ig))*DELTH0
 2101       CONTINUE
         ENDIF
C
C*    2.1.2 CURRENT REFRACTION.
C           -------------------
C
         DO 2102 IJ = IJS(IG),IJL(IG)
            DRCP(IJ) = (THDC(IJ,K,ig) + THDC(IJ,KP1,ig))*DELTH0
            DRCM(IJ) = (THDC(IJ,K,ig) + THDC(IJ,KM1,ig))*DELTH0
 2102    CONTINUE
C
C*    2.1.3 LOOP OVER FREQUENCIES.
C           ----------------------
C
         DO 2130 M=1,NFRE
            IF (ISHALLO.EQ.1) THEN
C
C*    2.1.3.1 DEEP WATER.
C             -----------
C
               MP1 = MIN(NFRE,M+1)
               MM1 = MAX(1,M-1)
               DFP = PI*2.1*DELFR0
C
C*    2.1.3.1.1 GROUP VELOCITIES.
C               -----------------
C
               CGS = GOM(M)*SD
               CGC = GOM(M)*CD
C
C*    2.1.3.1.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA(0) = CGS
               DPH(0) = CGC
               DO 2131 IJ=1,IJLT(IG)
                  DLA(IJ) = U(IJ,IG)*DELLA0 + CGS
                  DPH(IJ) = V(IJ,IG)*DELPH0 + CGC
 2131          CONTINUE
               DO 2132 IJ=IJS(IG),IJL(IG)
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1,ig))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2,ig))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1,ig))
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2,ig))
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) =  DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = DRCP(IJ)
                  DTHM = DRCM(IJ)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)

                  DTHP    = sidc(IJ,K,NFRE,ig) * DFP
                  DTC(IJ) = DTC(IJ) + 2.* ABS(DTHP)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHP+ABS(DTHP))*1.1
 2132          CONTINUE
            ELSE
CSHALLOW
C
C*    2.1.3.2 SHALLOW WATER.
C             --------------
C
               MP1 = MIN(NFRE,M+1)
               MM1 = MAX(1,M-1)
               DFP = DELFR0/FR(M)
               DFM = DELFR0/FR(MM1)
C
C*    2.1.3.2.1 GROUP VELOCITIES.
C               -----------------
C
               CGOND(0) = TCGOND(NDEPTH,M)
               DO 2133 IJ=1,IJLT(IG)
                  CGOND(IJ) = TCGOND(INDEP(IJ),M)
 2133          CONTINUE
C
C*    2.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA(0) = SD*CGOND(0)
               DPH(0) = CD*CGOND(0)
               DO 2134 IJ=1,IJLT(IG)
                  DLA(IJ) = U(IJ,IG)*DELLA0 + SD*CGOND(IJ)
                  DPH(IJ) = V(IJ,IG)*DELPH0 + CD*CGOND(IJ)
 2134          CONTINUE
               DO 2135 IJ=IJS(IG),IJL(IG)
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1,ig))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2,ig))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) = DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1,ig))
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2,ig))
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = FCONST(IJ,M)*DRDP(IJ) + DRCP(IJ)
                  DTHM = FCONST(IJ,M)*DRDM(IJ) + DRCM(IJ)
                  DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP = (sidc(IJ,K,M,ig) + sidc(IJ,K,MP1,ig))*DFP
                  DTHM = (sidc(IJ,K,M,ig) + sidc(IJ,K,MM1,ig))*DFM
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHM+ABS(DTHM))*1.1
 2135          CONTINUE
CSHALLOW
            ENDIF
C
C*    2.1.3.3 LOOP OVER GRIDPOINTS.
C             ---------------------
C
            DO 2136 IJ = IJS(IG),IJL(IG)
               F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                    + DPN(IJ) * F1(KLAT(IJ,2,ig),K  ,M)
     2                    + DPS(IJ) * F1(KLAT(IJ,1,ig),K  ,M)
     3                    + DLE(IJ) * F1(KLON(IJ,2,ig),K  ,M)
     4                    + DLW(IJ) * F1(KLON(IJ,1,ig),K  ,M)
     5                    + DTP(IJ) * F1(IJ        ,KP1,M)
     6                    + DTM(IJ) * F1(IJ        ,KM1,M)
     7                    + DOP(IJ) * F1(IJ        ,K  ,MP1)
     8                    + DOM(IJ) * F1(IJ        ,K  ,MM1)
 2136       CONTINUE
C
C*    BRANCH BACK TO 2.1.3 FOR NEXT FREQUENCY.
C
 2130    CONTINUE
C
C*    BRANCH BACK TO 2.1 FOR NEXT DIRECTION.
C
 2100 CONTINUE
C
C*    2.2 END OF PROPAGATION FOR CARTESIAN GRID
C*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
C         ------------------------------------------
cc
cc     here extract propagation source term diagnostics:
cc
      if(len_p2.eq.nang*nfre*niblo) then
       WRITE(6,*)'extracting diagnostics Sadv'
       do l=1,nfre
        do m=1,nang
         nstart=((l-1)*nang + m-1)*niblo
         do ip=ijs(ig),ijl(ig)
          sadv2(nstart+ip)=(F3(ip,m,l) - F1(ip,m,l))
         enddo
        enddo
       enddo
      endif


      RETURN
C
C ----------------------------------------------------------------------
C
C*    3. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
C*       WITHOUT OR DEPTH REFRACTION.
C        -------------------------------------------------
C
 3000 CONTINUE
C
      DELPRO = FLOAT(IDELPRO)
      DELTH0 = 0.25*DELPRO/DELTH
      DELPH0 = 0.5*DELPRO/DELPHI
      IF (ISHALLO.EQ.1) THEN
         DELLA0 = DELPRO/DELLAM
      ELSE
         DELLA0 = 0.5*DELPRO/DELLAM
      ENDIF
C
C*    3.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 3100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)*DELLA0
         CD = COSTH(K)*DELPH0
         SDA = ABS(SD)
         CDA = ABS(CD)
C
C*    3.1.1 COMPUTE GRID REFRACTION.
C           ------------------------
C
         SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
         SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R
         DO 3101 IJ = IJS(IG),IJL(IG)
            JH = KXLT(IJ,IG)
            TANPH = SINPH(JH)*DCO(IJ)
            DRGP(IJ) = TANPH*SP
            DRGM(IJ) = TANPH*SM
 3101    CONTINUE
C
C*    3.1.2 INDEX FOR ADJOINING POINTS.
C           ---------------------------
C
         IF (SD.LT.0) THEN
            IJLA = 2
         ELSE
            IJLA = 1
         ENDIF
         IF (CD.LT.0) THEN
            IJPH = 2
         ELSE
            IJPH = 1
         ENDIF
C
         IF (ISHALLO.EQ.1) THEN
C
C*    3.1.3 DEEP WATER.
C           -----------
C
C*    3.1.3.1 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
C             -----------------------------------------
C
            DO 3131 IJ=IJS(IG),IJL(IG)
               DLE(IJ) = DCO(IJ)*SDA
 3131       CONTINUE
            IF (CD.GT.0.) THEN
               DO 3132 IJ=IJS(IG),IJL(IG)
                  DTC(IJ) = DLE(IJ) + CDA*(DP2(IJ) + 1.)
                  DPN(IJ) = CDA*(DP1(IJ) + 1.)
 3132          CONTINUE
            ELSE
               DO 3133 IJ=IJS(IG),IJL(IG)
                  DTC(IJ) = DLE(IJ) + CDA*(DP1(IJ) + 1.)
                  DPN(IJ) = CDA*(DP2(IJ) + 1.)
 3133          CONTINUE
            ENDIF
C
C*    3.1.3.2 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
C             -----------------------------------------
C
            DO 3134 IJ=IJS(IG),IJL(IG)
               DTHP = DRGP(IJ)
               DTHM = DRGM(IJ)
               DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
               DTP(IJ) = -DTHP+ABS(DTHP)
               DTM(IJ) =  DTHM+ABS(DTHM)
 3134       CONTINUE
C
C*    3.1.3.3 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 3135 M=1,NFRE
C
C*    3.1.3.3.1 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 3136 IJ = IJS(IG),IJL(IG)
                  DTT = 1. - DTC(IJ)*GOM(M)
                  F3(IJ,K,M) = DTT*F1(IJ,K,M ) + GOM(M) *
     1                        (DPN(IJ) * F1(KLAT(IJ,IJPH,ig),K  ,M)
     2                       + DLE(IJ) * F1(KLON(IJ,IJLA,ig),K  ,M)
     3                       + DTP(IJ) * F1(IJ           ,KP1,M)
     4                       + DTM(IJ) * F1(IJ           ,KM1,M))
 3136          CONTINUE
C
C*    BRANCH BACK TO 3.1.3.3 FOR NEXT FREQUENCY.
C
 3135       CONTINUE
         ELSE
CSHALLOW
C
C*    3.1.4 SHALLOW WATER.
C           --------------
C
C
C*    3.1.4.1 COMPUTE DEPTH REFRACTION.
C             -------------------------
C
         IF (IREFRA.EQ.1) THEN
            DO 3141 IJ = IJS(IG),IJL(IG)
               DRDP(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KP1,ig))*DELTH0
               DRDM(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KM1,ig))*DELTH0
 3141       CONTINUE
         ENDIF
C
C*    3.1.4.2 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 3142 M=1,NFRE
C
C*    3.1.4.2.1 GROUP VELOCITIES.
C               -----------------
C
               CGOND(0) = TCGOND(NDEPTH,M)
               DO 3143 IJ=1,IJLT(IG)
                  CGOND(IJ) = TCGOND(INDEP(IJ),M)
 3143          CONTINUE
C
C*    3.1.4.3.2 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
C               -----------------------------------------
C
               IF (SD.GT.0.) THEN
                  DO 3144 IJ=IJS(IG),IJL(IG)
                     DTC(IJ) = 1. - DCO(IJ)*SDA*
     1                         (CGOND(KLON(IJ,2,ig)) + CGOND(IJ))
                     DLE(IJ) = DCO(IJ)*SDA*
     1                         (CGOND(KLON(IJ,1,ig)) + CGOND(IJ))
 3144             CONTINUE
               ELSE
                  DO 3145 IJ=IJS(IG),IJL(IG)
                     DTC(IJ) = 1. - DCO(IJ)*SDA*
     1                         (CGOND(KLON(IJ,1,ig)) + CGOND(IJ))
                     DLE(IJ) = DCO(IJ)*SDA*
     1                         (CGOND(KLON(IJ,2,ig)) + CGOND(IJ))
 3145             CONTINUE
               ENDIF
               IF (CD.GT.0.) THEN
                  DO 3146 IJ=IJS(IG),IJL(IG)
                     DTC(IJ) = DTC(IJ) - CDA*
     1                      (CGOND(KLAT(IJ,2,ig))*DP2(IJ) + CGOND(IJ))
                     DPN(IJ) = CDA*
     1                      (CGOND(KLAT(IJ,1,ig))*DP1(IJ) + CGOND(IJ))
 3146             CONTINUE
               ELSE
                  DO 3147 IJ=IJS(IG),IJL(IG)
                     DTC(IJ) = DTC(IJ) - CDA*
     1                      (CGOND(KLAT(IJ,1,ig))*DP1(IJ) + CGOND(IJ))
                     DPN(IJ) = CDA*
     1                      (CGOND(KLAT(IJ,2,ig))*DP2(IJ) + CGOND(IJ))
 3147             CONTINUE
               ENDIF
C
C*    3.1.4.2.3 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
C               -----------------------------------------
C
               IF (IREFRA.EQ.0) THEN
                  DO 3148 IJ=IJS(IG),IJL(IG)
                     DTHP = DRGP(IJ)*CGOND(IJ)
                     DTHM = DRGM(IJ)*CGOND(IJ)
                     DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 3148             CONTINUE
               ELSE
                  DO 3149 IJ=IJS(IG),IJL(IG)
                     DTHP = DRGP(IJ)*CGOND(IJ)+FCONST(IJ,M)*DRDP(IJ)
                     DTHM = DRGM(IJ)*CGOND(IJ)+FCONST(IJ,M)*DRDM(IJ)
                     DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 3149             CONTINUE
                ENDIF
C
C*    3.1.4.2.4 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 3150 IJ = IJS(IG),IJL(IG)
                  F3(IJ,K,M) = DTC(IJ)*F1(IJ,K,M )
     1                       + DPN(IJ) * F1(KLAT(IJ,IJPH,ig),K  ,M)
     2                       + DLE(IJ) * F1(KLON(IJ,IJLA,ig),K  ,M)
     3                       + DTP(IJ) * F1(IJ           ,KP1,M)
     4                       + DTM(IJ) * F1(IJ           ,KM1,M)
 3150          CONTINUE
C
C*    BRANCH BACK TO 3.1.4.2 FOR NEXT FREQUENCY.
C
 3142       CONTINUE
CSHALLOW
         ENDIF
C
C*    BRANCH BACK TO 3.1 FOR NEXT DIRECTION.
C
 3100 CONTINUE
C
C*    3.2 END OF PROPAGATION FOR SPHERICAL GRID
C*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
C         -----------------------------------------------
C
cc
cc     here extract propagation source term diagnostics:
cc
      if(len_p2.eq.nang*nfre*niblo) then
       WRITE(6,*)'extracting diagnostics Sadv'
       do l=1,nfre
        do m=1,nang
         nstart=((l-1)*nang + m-1)*niblo
         do ip=ijs(ig),ijl(ig)
          sadv2(nstart+ip)=(F3(ip,m,l) - F1(ip,m,l))
         enddo
        enddo
       enddo
      endif

      RETURN
C
C ----------------------------------------------------------------------
C
C*    4. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
C*       WITH DEPTH AND CURRENT REFRACTION.
C        -------------------------------------------------
C
 4000 CONTINUE
C
      DELPRO = FLOAT(IDELPRO)
      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELPRO/DELLAM
      DELFR0 = 0.25*DELPRO/(0.1*ZPI)
C
C*    4.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 4100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)*DELLA0
         CD = COSTH(K)*DELPH0
C
C*    4.1.1 COMPUTE GRID REFRACTION.
C           ------------------------
C
         SP = DELTH0*(SINTH(K)+SINTH(KP1))/R
         SM = DELTH0*(SINTH(K)+SINTH(KM1))/R
         DO 4111 IJ = IJS(IG),IJL(IG)
            JH = KXLT(IJ,IG)
            TANPH = SINPH(JH)*DCO(IJ)
            DRGP(IJ) = TANPH*SP
            DRGM(IJ) = TANPH*SM
 4111    CONTINUE
C
C*    4.1.2 COMPUTE DEPTH REFRACTION.
C           -------------------------
C
         IF (ISHALLO.NE.1) THEN
            DO 4121 IJ = IJS(IG),IJL(IG)
               DRDP(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KP1,ig))*DELTH0
               DRDM(IJ) = (THDD(IJ,K,ig) + THDD(IJ,KM1,ig))*DELTH0
 4121       CONTINUE
         ENDIF
C
C*    4.1.3 COMPUTE CURRENT REFRACTION.
C           ---------------------------
C
         DO 4131 IJ = IJS(IG),IJL(IG)
            DRCP(IJ) = (THDC(IJ,K,ig) + THDC(IJ,KP1,ig))*DELTH0
            DRCM(IJ) = (THDC(IJ,K,ig) + THDC(IJ,KM1,ig))*DELTH0
 4131    CONTINUE
C
C*    4.1.4 LOOP OVER FREQUENCIES.
C           ----------------------
C
         DO 4140 M=1,NFRE
            MP1 = MIN(NFRE,M+1)
            MM1 = MAX(1,M-1)
            IF (ISHALLO.EQ.1) THEN
C
C*    4.1.4.1 DEEP WATER.
C             -----------
C
C*    4.1.4.1.1 GROUP VELOCITIES.
C               -----------------
C
               DFP = PI*2.1*DELFR0
               CGS = GOM(M)*SD
               CGC = GOM(M)*CD
C
C*    4.1.4.1.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA( 0) = CGS
               DPH( 0) = CGC
               DO 4141 IJ=1,IJLT(IG)
                  DLA(IJ) = (U(IJ,IG)*DELLA0 + CGS)*DCO(IJ)
                  DPH(IJ) =  V(IJ,IG)*DELPH0 + CGC
 4141          CONTINUE
               DO 4142 IJ=IJS(IG),IJL(IG)
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1,ig))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2,ig))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1,ig))*DP1(IJ)
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2,ig))*DP2(IJ)
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = DRGP(IJ)*GOM(M) + DRCP(IJ)
                  DTHM = DRGM(IJ)*GOM(M) + DRCM(IJ)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP =  sidc(IJ,K,NFRE,ig) * DFP
                  DTC(IJ) =  DTC(IJ) + 2. * ABS(DTHP)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHP+ABS(DTHP))*1.1
 4142          CONTINUE
            ELSE
CSHALLOW
C
C*    4.1.4.2 SHALLOW WATER.
C             --------------
C
C*    4.1.4.2.1 GROUP VELOCITIES.
C               -----------------
C
               DFP = DELFR0/FR(M)
               DFM = DELFR0/FR(MM1)
               CGOND(0) = TCGOND(NDEPTH,M)
               DO 4143 IJ=1,IJLT(IG)
                  CGOND(IJ) = TCGOND(INDEP(IJ),M)
 4143          CONTINUE
C
C*    4.1.4.2.2 LON/LAT/DIR WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------------------
C
               DLA( 0) = SD*CGOND(0)
               DPH( 0) = CD*CGOND(0)
               DO 4144 IJ=1,IJLT(IG)
                  DLA(IJ) = (U(IJ,IG)*DELLA0 + SD*CGOND(IJ))*DCO(IJ)
                  DPH(IJ) =  V(IJ,IG)*DELPH0 + CD*CGOND(IJ)
 4144          CONTINUE
               DO 4145 IJ=IJS(IG),IJL(IG)
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1,ig))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2,ig))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1,ig))*DP1(IJ)
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2,ig))*DP2(IJ)
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP=DRGP(IJ)*CGOND(IJ)+FCONST(IJ,M)*DRDP(IJ)+DRCP(IJ)
                  DTHM=DRGM(IJ)*CGOND(IJ)+FCONST(IJ,M)*DRDM(IJ)+DRCM(IJ)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP = (sidc(IJ,K,M,ig) + sidc(IJ,K,MP1,ig))*DFP
                  DTHM = (sidc(IJ,K,M,ig) + sidc(IJ,K,MM1,ig))*DFM
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHM+ABS(DTHM))*1.1
 4145          CONTINUE
CSHALLOW
            ENDIF
C
C*    4.1.4.3 LOOP OVER GRIDPOINTS.
C             ---------------------
C
            DO 4146 IJ = IJS(IG),IJL(IG)
               F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                    + DPN(IJ) * F1(KLAT(IJ,2,ig),K  ,M)
     2                    + DPS(IJ) * F1(KLAT(IJ,1,ig),K  ,M)
     3                    + DLE(IJ) * F1(KLON(IJ,2,ig),K  ,M)
     4                    + DLW(IJ) * F1(KLON(IJ,1,ig),K  ,M)
     5                    + DTP(IJ) * F1(IJ        ,KP1,M)
     6                    + DTM(IJ) * F1(IJ        ,KM1,M)
     7                    + DOP(IJ) * F1(IJ        ,K  ,MP1)
     8                    + DOM(IJ) * F1(IJ        ,K  ,MM1)
 4146          CONTINUE
C
C*    BRANCH BACK TO 4.1.4 FOR NEXT FREQUENCY.
C
 4140    CONTINUE
C
C*    BRANCH BACK TO 4.2 FOR NEXT DIRECTION.
C
 4100 CONTINUE
C
C*    4.4 END OF PROPAGATION FOR SPHERICAL GRID
C*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
C         ------------------------------------------
C

cc
cc     here extract propagation source term diagnostics:
cc
      if(len_p2.eq.nang*nfre*niblo) then
       WRITE(6,*)'extracting diagnostics Sadv'
       do l=1,nfre
        do m=1,nang
         nstart=((l-1)*nang + m-1)*niblo
         do ip=ijs(ig),ijl(ig)
          sadv2(nstart+ip)=(F3(ip,m,l) - F1(ip,m,l))
         enddo
        enddo
       enddo
      endif


      RETURN
      END
