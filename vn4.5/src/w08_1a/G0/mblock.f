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

      SUBROUTINE MBLOCK (ka,ke,ipp,ia1,
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

      integer   ka           ! in
      integer   ke           ! in

      integer   ipp(ngy)     ! in
        logical ia1(ngx,ngy)   ! in land-sea mask

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
C**** *MBLOCK* - ROUTINE TO ARRANGE WAMODEL GRID FOR ONE BLOCK.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       *MBLOCK* ARRANGES WAMODEL GRID FOR A BLOCK AND
C                COMPUTES VARIOUS MODEL CONSTANTS
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MBLOCK (IA1, KA, KE, IPP)*
C          *IA1*     - TOPOGRAPHIC DATA.
C          *KA*      - NUMBER OF FIRST LAT IN BLOCK.
C          *KE*      - NUMBER OF LAST LAT IN BLOCK.
C          *IPP*     - NUMBER OF SEA POINTS PER LAT.
C
C     METHOD.
C     -------
C
C       THE LAND POINTS ARE REMOVED. ALL MODEL CONSTANTS WHICH ARE
C       GRID DEPENDENT ARE COMPUTED AND STORED IN THE COMMON BLOCKS.
C
C     EXTERNALS.
C     ----------
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C*    1. UPDATE BLOCK NUMBER AND INITIALIZES ARRAYS.
C        -------------------------------------------
C
       iu06=6

 1000 CONTINUE
      IGL = IGL + 1
      IF (IGL.GT.NBLO) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* MORE BLOCKS THAN DIMENSION ALLOWS.         *'
         WRITE (IU06,*) '* BLOCK NUMBER IS                 IGL = ', IGL
         WRITE (IU06,*) '* DIMENSION IS                   NBLO = ', NBLO
         WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS      KA = ', KA
         WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS      KE = ', KE
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         icode=-1
         call abort
      ENDIF
      DO 1001 IJ=1,NIBLO
         IXLG(IJ,IGL) = 0
         KXLT(IJ,IGL) = 0
 1001 CONTINUE

C ----------------------------------------------------------------------
C
C*    2. THE FIRST AND LAST BLOCK MUST CONTAIN MORE THAN 2
C*       ALL OTHER BLOCKS MORE  THAN 3 LATITUDES.
C        -------------------------------------------------
C
 2000 CONTINUE
cc
cc
ccUM the if test below was originally :
ccUM 1    ((KA.NE.1) .AND. (KE.EQ.NY) .AND. (KE-KA.LT.2))) THEN
ccUM                          ~~~~
ccUM  corrected to .NE. after UM Review
ccUM  S Kelsall / M Holt 26/7/95
cc
cc    this now matches what is described in the comment block above
cc
      IF ((KE.EQ.1) .OR. (KA.EQ.NY) .OR.
     1    ((KA.NE.1) .AND. (KE.NE.NY) .AND. (KE-KA.LT.2))) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* BLOCK LENGTH IS TO SHORT.                  *'
         WRITE (IU06,*) '* LESS THAN 2 LATITUDES IN FIRST OR LAST, OR *'
         WRITE (IU06,*) '* LESS THAN 3 LATITUDES IN OTHER BLOCKS.     *'
         WRITE (IU06,*) '* BLOCK NUMBER IS               IGL = ', IGL
         WRITE (IU06,*) '* BLOCK LENGTH IS             NIBLO = ', NIBLO
         WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS    KA = ', KA
         WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS    KE = ', KE
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         icode=-1
         call abort
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. COMPUTE INDICES OF FIRST, SECOND, BEFORE LAST, AND LAST LAT.
C        -----------------------------------------------------------
C
 3000 CONTINUE
      IF (KA.EQ.1) THEN
         IJS (IGL) = 1
         IJL2(IGL) = IPP(1)
      ELSE
        IJS (IGL) = IPP(KA)+1
        IJL2(IGL) = IPP(KA)+IPP(KA+1)
      ENDIF
      IJLT(IGL) = 0
      DO 3001 K=KA,KE
         IJLT(IGL) = IJLT(IGL)+IPP(K)
 3001 CONTINUE
      IF (KE.EQ.NY) THEN
         IJL (IGL) = IJLT(IGL)
      ELSE
         IJL (IGL) = IJLT(IGL)-IPP(KE)
      ENDIF
      IJLS(IGL) = IJL(IGL)-IPP(KE-1)+1
C
C ----------------------------------------------------------------------
C
C*    4. REMOVE LAND POINTS AND STORE COS AND SIN OF LAT.
C        ------------------------------------------------
C
 4000 CONTINUE
      IP = 0
      DO 4001 K=KA,KE
         DO 4002 I=1,NX
            IF (.not.IA1(I,K)) THEN ! a sea point
              IP = IP+1
              IXLG(IP,IGL) = I
              KXLT(IP,IGL) = K
            ENDIF
 4002    CONTINUE
 4001 CONTINUE
      IF (IP.NE.IJLT(IGL)) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* TOTAL NUMBER OF SEAPOINTS DO NOT MATCH.    *'
         WRITE (IU06,*) '* BLOCK NUMBER                    IGL = ', IGL
         WRITE (IU06,*) '* NO. OF SEAPOINTS COUNTED         IP = ', IP
         WRITE (IU06,*) '* NO. OF SEAPOINTS EXPECTED IJLT(IGL) = ',
     1                                            IJLT(IGL)
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         icode=-1
         call abort
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. PRINTER PROTOCOL OF BLOCK.
C        --------------------------
C
 5000 CONTINUE
      IF (IGL.EQ.1) THEN
         WRITE (IU06,'(1H0,'' BLOCKING INFORMATION:'')')
         WRITE (IU06,'(1H ,''            LATITUDES   '',
     1                     ''   SECOND LAT. INDEX '',
     2                     '' SECOND TO LAST LAT  '',
     3                     ''   TOTAL'')')
         WRITE (IU06,'(1H ,''  NO     SOUTH     NORTH'',
     1                     ''     START       END'',
     2                     ''     START       END'',
     3                     ''    POINTS'')')
      ENDIF
      WRITE (IU06,'(1X,I4,2F10.2,5I10)')
     1        IGL, AMOSOP+(KA-1)*XDELLA, AMOSOP+(KE-1)*XDELLA,
     2        IJS(IGL), IJL2(IGL), IJLS(IGL), IJL(IGL), IJLT(IGL)

      RETURN
      END
