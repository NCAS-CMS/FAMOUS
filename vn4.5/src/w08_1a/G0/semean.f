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

      SUBROUTINE SEMEAN (F3, IJS, IJL,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
     & EMEAN, FMEAN, THQ, AKMEAN,

     & icode)

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
C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
      real EMEAN(NIBLO)  ! total energy
      real FMEAN(NIBLO)  ! mean frequency
      real THQ(NIBLO)    ! mean wave direction (radians)
      real AKMEAN(NIBLO) ! mean wave number
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
C**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C
C     S.D. HASSELMANN.
C     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER
C
C*    PURPOSE.
C     --------
C
C       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SEMEAN(F3, IJS, IJL)*
C          *F3*  - SPECTRUM.
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C
C     METHOD.
C     -------
C
C       NONE.
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
      DIMENSION F3(0:NIBLO,NANG,NFRE), TEMP(NIBLO)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE ENERGY ARRAY.
C        ------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
         EMEAN(IJ) = 0.
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
C        -----------------------------------------
C
 2000 CONTINUE
      DO 2001 M=1,NFRE
         DO 2002 IJ=IJS,IJL
            TEMP(IJ) = 0.
 2002   CONTINUE
        DO 2003 K=1,NANG
           DO 2004 IJ=IJS,IJL
              TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)
 2004      CONTINUE
 2003   CONTINUE
        DO 2005 IJ=IJS,IJL
           EMEAN(IJ) = EMEAN(IJ)+DFIM(M)*TEMP(IJ)
 2005   CONTINUE
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. ADD TAIL ENERGY.
C        ----------------
C
 3000 CONTINUE
      DELT25 = FR(NFRE)/4.*DELTH
      DO 3001 IJ=IJS,IJL
         EMEAN(IJ) = EMEAN(IJ)+DELT25*TEMP(IJ)
 3001 CONTINUE

      RETURN
      END
