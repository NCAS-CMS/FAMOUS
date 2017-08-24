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

      SUBROUTINE MGRID (IA1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

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

C ----------------------------------------------------------------------
C
C**** *MGRID* - ROUTINE TO ARRANGE WAMODEL GRID.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO ARRANGE WAMODEL GRID FOR A GIVEN AREA AND COMPUTE VARIOUS
C       MODEL CONSTANTS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MGRID (IA1)*
C          *IA1*     - TOPOGRAPHIC DATA OF PART
C
C     METHOD.
C     -------
C
C       THE NUMBER OF SEA POINTS PER LATITUDE IS COUNTED AND MODEL
C       BLOCKS OF MAXIMUM LENGTH OF NIBLO ARE CONSTRUCTED.
C
C     EXTERNALS.
C     ----------
C
C       *MBLOCK*    - SUB. TO GENERATE A BLOCK.
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
      integer   IPP(NGY)
C
C          *IPP*   INTEGER   NUMBER OF SEA POINTS PER LATITUDE.
C
C ----------------------------------------------------------------------
C
C*    1. COUNT NUMBER OF SEA POINTS PER LATITUDE.
C        ----------------------------------------
C
 1000 CONTINUE

      DO 1001 K=1,NY
         IPP(K) = 0
         DO 1002 I=1,NX
            IF (.not.IA1(I,K)) THEN  ! a sea point
              IPP(K) = IPP(K) + 1
            ENDIF
 1002    CONTINUE
 1001 CONTINUE

C ----------------------------------------------------------------------
C
C*    2. MAKE BLOCKS.
C        ------------
C
      IGL=0
 2000 CONTINUE
      IL = 0
      KA = 1

      DO 2001 K = 1,NY
         IL = IL + IPP(K)
         IF (IL.GT.NIBLO) THEN
            KE = K-1

            CALL MBLOCK (KA, KE, IPP, ia1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

         if(icode.ne.0) then
          WRITE(6,*)'calling abort in setupwv : mgrid'
          WRITE(6,*)'icode ',icode,' returned from mblock'
          call abort
         endif

            KA = KE-1
ccmh     the line following starts counting from line ke-1
ccmh                                             (=ka for next block)
ccmh     line ke+1 is the current line k in the loop

            IL = IPP(KA)+IPP(KE)+IPP(KE+1)
         ENDIF
 2001 CONTINUE

      CALL MBLOCK (KA, NY, IPP, ia1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     &icode)

      if(icode.ne.0) then
       WRITE(6,*)'calling abort in setupwv : mgrid'
       WRITE(6,*)'icode ',icode,' returned from mblock'
       call abort
      endif

      RETURN
      END
