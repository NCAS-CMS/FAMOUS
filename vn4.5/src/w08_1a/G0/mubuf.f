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

      SUBROUTINE MUBUF (IA1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

     & icode)

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
ccmh
cc revised for UKMO UMwave implementation to do all blocks.
cc arrays klat klon given third dimension of nblo
ccmh
C**** *MUBUF* - ROUTINE TO ARRANGE COMMON UBUF FOR ONE BLOCK.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MUBUF (IA1, IG, IU08, IU18, IFORM)*
C          *IA1*     - TOPOGRAPHIC DATA.
C          *IG*      - BLOCK NUMBER.
C          *IU08*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
C                      COMMON UBUF (UNFORMATED)
C          *IU18*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
C                      COMMON UBUF (FORMATED)
C          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
C                                           = 2 FORMATED
C                                           OTHERWISE BOTH
C
C     METHOD.
C     -------
C
C       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE
C       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.
C       THE FINAL COMMON UBUF IS WRITTEN OUT.
C
C     EXTERNALS.
C     ----------
C
C       *OUTUBUF*   - WRITE OUT COMMON UBUF.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
        logical ia1(ngx,ngy)   ! in land-sea mask land=T
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE ARRAYS.
C        ------------------
C
 1000 CONTINUE
      do ig=1,igl
      DO 1001 IJ=1,NIBLO
         KLAT(IJ,1,ig) = 0
         KLAT(IJ,2,ig) = 0
         KLON(IJ,1,ig) = 0
         KLON(IJ,2,ig) = 0
 1001 CONTINUE
      enddo
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.
C        ----------------------------------------
C
 2000 CONTINUE
C
C*    2.1 LONGITUDE NEIGHBOURS.
C         ---------------------
C
      do ig=1,igl
      DO 2100 IP = 1,IJLT(IG)
         I = IXLG(IP,IG)
         K = KXLT(IP,IG)
c
C       * now using the UM logical landsea mask with land=T
c
         IF (I.GT.1) THEN
           IF (.not.IA1(I-1,K)) KLON(IP,1,ig) = IP-1
         ELSE
            IF (IPER.EQ.1 .AND. .not.IA1(NX,K)) THEN
               KLON(IP,1,ig) = IP
               DO 2101 IH=2,NX
                IF (.not.IA1(IH,K)) KLON(IP,1,ig) = KLON(IP,1,ig)+1
 2101          CONTINUE
            ENDIF
         ENDIF
         IF (I.LT.NX) THEN
          IF (.not. IA1(I+1,K) ) KLON(IP,2,ig) = IP+1
         ELSE
           IF (IPER.EQ.1 .AND. .not. IA1(1,K)) THEN
              KLON(IP,2,ig) = IP
              DO 2102 IH=NX-1,1,-1
               IF (.not. IA1(IH,K)) KLON(IP,2,ig) = KLON(IP,2,ig)-1
 2102         CONTINUE
            ENDIF
         ENDIF
 2100 CONTINUE
C
C*    2.2 LATITUDE NEIGHBOURS.
C         --------------------
C
      DO 2200 IP = 1,IJLT(IG)
         I = IXLG(IP,IG)
         K = KXLT(IP,IG)
         IF (K.GT.1) THEN
            IF (.not. IA1(I,K-1)) THEN
               DO 2201 IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.I .AND.
     1                KXLT(IH,IG).EQ.K-1) GOTO 2202
 2201          CONTINUE
 2202          CONTINUE
               KLAT(IP,1,ig) = IH
            ENDIF
         ENDIF
         IF (K.LT.NY) THEN
            IF (.not. IA1(I,K+1)) THEN
               DO 2203 IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.I .AND.
     1                KXLT(IH,IG).EQ.K+1) GOTO 2204
 2203          CONTINUE
 2204          CONTINUE
               KLAT(IP,2,ig) = IH
           ENDIF
         ENDIF
 2200 CONTINUE
      enddo
C
      RETURN
      END
