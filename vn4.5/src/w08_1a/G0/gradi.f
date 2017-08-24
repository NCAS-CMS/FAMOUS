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

      SUBROUTINE GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
     &                  DVPHI, DVLAM,
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
C**** *GRADI* - CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS.
C
C     K.P. HUBBERT              AUGUST   1988
C     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C       CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS OF A BLOCK.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
C                      DVPHI, DVLAM)*
C          *IG*     - BLOCK NUMBER.
C          *IREFRA* - REFRACTION OPTION.
C          *DDPHI*  - LATITUDE DEPTH GRADIENT.
C          *DDLAM*  - LONGITUDE DEPTH GRADIENT.
C          *DUPHI*  - LATITUDE  U-COMPONENT GRADIENT.
C          *DULAM*  - LONGITUDE U-COMPONENT GRADIENT.
C          *DVPHI*  - LATITUDE  V-COMPONENT GRADIENT.
C          *DVLAM*  - LONGITUDE V-COMPONENT GRADIENT.
C
C     METHOD.
C     ------
C
C       CENTRAL DIFFERENCING FOR DEPTH AND CURRENT VELOCITY GRADIENTS.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
cccc*CALL PARALL
C
ccc*CALL COMCURR
C      PARAMETER (CO = 1.1, ALN00 = 24.1589,
C    1           FRE0 = CO-1., PI2G = ZPI/G)
ccc*CALL COMGRID
C
ccc*CALL COMSHAL
C
ccc*CALL COMUBUF
C
C ----------------------------------------------------------------------
C
      DIMENSION DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC),
     1          DVPHI(NIBLC), DVLAM(NIBLC)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE.
C        -----------
C
 1000 CONTINUE
      DELPHI2 = 1./(DELPHI*2.0)
      DELLAM2 = 1./(DELLAM*2.0)
C
C ----------------------------------------------------------------------
C
C*    2. CALCULATE DEPTH GRADIENTS.
C        --------------------------
C
 2000 CONTINUE
      DO 2001 IJ=1,NIBLD
         DDPHI(IJ) = 0.0
         DDLAM(IJ) = 0.0
 2001 CONTINUE
      DO 2002 IJ=IJS(IG),IJL(IG)
         IPP = KLAT(IJ,2,ig)
         IPM = KLAT(IJ,1,ig)
         IF (IPP.GT.0 .AND. IPM.GT.0) THEN
            DDPHI(IJ) = (DEPTH(IPP,IG)-DEPTH(IPM,IG))*DELPHI2
         ENDIF
         ILP = KLON(IJ,2,ig)
         ILM = KLON(IJ,1,ig)
         IF (ILP.GT.0 .AND. ILM.GT.0) THEN
            DDLAM(IJ) = (DEPTH(ILP,IG)-DEPTH(ILM,IG))*DELLAM2
         ENDIF
 2002 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
C        -------------------------------------
C
 3000 CONTINUE
      IF (IREFRA.EQ.2) THEN
         DO 3001 IJ=1,NIBLC
            DUPHI(IJ) = 0.0
            DULAM(IJ) = 0.0
            DVPHI(IJ) = 0.0
            DVLAM(IJ) = 0.0
 3001    CONTINUE
         DO 3002 IJ=IJS(IG),IJL(IG)
            IPP = KLAT(IJ,2,ig)
            IPM = KLAT(IJ,1,ig)
            IF (IPP.GT.0 .AND. IPM.GT.0) THEN
               DUPHI(IJ) = (U(IPP,IG)-U(IPM,IG))*DELPHI2
               DVPHI(IJ) = (V(IPP,IG)-V(IPM,IG))*DELPHI2
            ENDIF
            ILP = KLON(IJ,2,ig)
            ILM = KLON(IJ,1,ig)
            IF (ILP.GT.0 .AND. ILM.GT.0) THEN
               DULAM(IJ) = (U(ILP,IG)-U(ILM,IG))*DELLAM2
               DVLAM(IJ) = (V(ILP,IG)-V(ILM,IG))*DELLAM2
            ENDIF
 3002    CONTINUE
      ENDIF

      RETURN
      END
