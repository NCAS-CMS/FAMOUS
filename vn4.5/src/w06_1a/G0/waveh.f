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
C     13. SUBROUTINE WAVEH
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

C
C     DOCUMENTATION
C
C     SEE WAVE MODEL DOCUMENTATION PAPER.
C
C     DESCRIPTION
C
C     THIS ROUTINE CALCULATES WAVE HEIGHTS ETC FROM 2D SPECTRA
C
C *********************************************************************
       subroutine waveh(ia1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
     &   l_wvtra,ndata,ngrid,nfldmax,len1,energy,windsp,windir,wh,icode)
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

C
c      * set max number of wavetrains to search for *
       parameter(kwtmax=4)

C UM constants

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
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------


       INTEGER ICODE
       real wh(ngrid,nfldmax) ! array for output fields ! OUT
       real energy(len1),windir(ndata),windsp(ndata) ! IN
C
C      # LOCAL arrays for use by wavetrain partitioning
C
CC note can rearange and not declare these arrays here by calling
CC with array WH - have kwtmax consec wavehts then kwtmax periods etc.
CC       real pswh(ndata,kwtmax)  ! height of wavetrains
CC       real perio(ndata,kwtmax) ! period of wavetrain
CC       real pdir(ndata,kwtmax)  ! direction of wavetrain
       integer kwtot(ndata)     ! number of wavetrains present
C
C     LOCAL ARRAYS
C
      REAL      ECOS(ndata),ESIN(ndata),EMAX(ndata),
     -          ESIN1(ndata),ECOS1(ndata),EMSW(ndata),EMWS(ndata),
     -          ESIN2(ndata),ECOS2(ndata),EWS(ndata)


C

CC local work arrays replacing WH(ip,??) for development of output
CC on full grid
CC work7(ip) equates to WH(ip,7) in the original WAVEH routine

      REAL    work1(ndata),work4(ndata),work5(ndata)
      REAL    work7(ndata),work10(ndata),work11(ndata)

C
      INTEGER   IMAX(ndata),IMSW(ndata),IMWS(ndata)
C
      LOGICAL   LFLIM,LTRUE,LWS(ndata)
      LOGICAL   l_wvtra

      LOGICAL IA1(ngx,ngy)
C
C ----------------------------------------------------------------------
C        VALUES STORED IN WH(IP,I) IN THIS SUBROUTINE
C  I=  1       TOTAL WAVE HEIGHT
C      2       MEAN DIRECTION
C      3       PRINCIPAL DIRECTION
C      4       ZERO UP CROSSING PERIOD
C      5       MEAN PERIOD
C      6       PEAK PERIOD
C      7       TOTAL WINDSEA WAVE HEIGHT
C      8       MEAN DIRECTION
C      9       PRINCIPAL DIRECTION
C     10       ZERO UP CROSSING PERIOD
C     11       MEAN PERIOD
C     12       PEAK PERIOD
C     13       TOTAL SWELL WAVE HEIGHT
C     14       MEAN DIRECTION
C     15       PRINCIPAL DIRECTION
C     16       ZERO UP CROSSING PERIOD
C     17       MEAN PERIOD
C     18       PEAK PERIOD
C     19
C     20
C ----------------------------------------------------------------------
C LOCALLY USED CONSTANTS

      PI2  = 2.0*PI
      PIO2 = PI/2.0

CCC UKMO value      G14  = 0.14*G : use 0.13*G for 10m wind -
CCC see wavetrain code of Anne Guillaume

      G14  = 0.13*G

CL 13.0 ZERO DATA ARRAY
C ---------------------
C
c      * this loop initialises array WH

      do i=1,nfldmax
       DO II=1,ngrid
        WH(II,I) = RMDI
       ENDDO
      ENDDO

      do ip=1,ndata
       if(windsp(ip).lt.0.5) then
        windsp(ip)=0.5
       endif
      enddo
C
      EMIN = 1.0E-7
      DO IP=1,NDATA
       ECOS(IP) = 0.0
       ESIN(IP) = 0.0
       ECOS1(IP) = 0.0
       ESIN1(IP) = 0.0
       EWS (IP) = 0.0
       EMAX(IP) = EMIN
       EMWS(IP) = EMIN
       EMSW(IP) = EMIN
       IMAX(IP) = 0
       IMWS(IP) = 0
       IMSW(IP) = 0
       work1(ip)=0.
       work4(ip)=0.
       work5(ip)=0.
       work7(ip)=0.
       work10(ip)=0.
       work11(ip)=0.
      ENDDO
C
CL 13.1 FIRST GUESS AT WINDSEA ENERGY - EWS
C --------------------------------------------
C
      DO 10 L=1,nfre
CCC    DTHDF = DTH*DF(L)
       DTHDF = DFIM(L)
       DO 11 M=1,nang
CC NOTE need to ensure this indexing is ok for energy
        NST = ((L-1)*nang+M-1)*NDATA
        DO 12 IP=1,NDATA
         FPMC = 0.8*G14/WINDSP(IP)
         LFLIM = FR(L).GT.FPMC
         IF (L.EQ.nfre) LFLIM=.TRUE.
         ANG = ABS(TH(M)-WINDIR(IP))
         IF (ANG.GT.PI) ANG = PI2-ANG
         ANG = PIO2-ANG
         LTRUE= ANG.GT.-0.01
         IF (LTRUE.AND.LFLIM) THEN                   ! WINDSEA
          EWS(IP) = EWS(IP) + ENERGY(NST+IP) * DTHDF
         END IF
   12   CONTINUE
   11  CONTINUE
   10 CONTINUE
C
CL 13.2 CALCULATE VALUES USING NEW CUT OFF FREQUENCY
C --------------------------------------------------
      DO 20 L=1,nfre
ccc    DTHDF = DTH*DF(L)
       DTHDF = DFIM(L)
       DO 21 M=1,nang
        NST = ((L-1)*nang+M-1)*NDATA
        DO 22 IP=1,NDATA
         FPM = G14/WINDSP(IP)
         EPM = (FPM**4)*0.0001
         IF (EWS(IP).GT.EPM) EWS(IP)=EPM
         IF (EWS(IP).GT.0.0) THEN

CC CHECK THIS COEFFICIENT 0.31    elsewhere is 0.33 ??

          FPMC = 0.8*FPM*((EPM/EWS(IP))**0.31)
         ELSE
          FPMC = 0.8*G14
         END IF
         LFLIM = FR(L).GT.FPMC
         IF (L.EQ.nfre) LFLIM=.TRUE.
         ANG = ABS(TH(M)-WINDIR(IP))
         IF (ANG.GT.PI) ANG = PI2-ANG
         ANG = PIO2-ANG
         LTRUE= ANG.GT.-0.01
         LWS(IP) = LTRUE.AND.LFLIM
C
C *** WINDSEA
C
         IF (LWS(IP)) THEN
          Work7(IP) =  Work7 (IP)+ENERGY(NST+IP) * DTHDF
          Work10(IP)=  Work10(IP)+ENERGY(NST+IP)*FR(L)*FR(L)*DTHDF
          Work11(IP)=  Work11(IP)+ENERGY(NST+IP)*FR(L)*DTHDF
          ECOS1(IP) = ECOS1(IP) + ENERGY(NST+IP)*COSTH(M)*DTHDF
          ESIN1(IP) = ESIN1(IP) + ENERGY(NST+IP)*SINTH(M)*DTHDF
         END IF
C
C *** TOTAL  VALUES
C
cccc   if(energy(nst+ip).gt.0.0001) then
         Work1(IP) = Work1(IP)+ENERGY(NST+IP) * DTHDF
         Work4(IP) = Work4(IP)+ENERGY(NST+IP)*FR(L)*FR(L)*DTHDF
         Work5(IP) = Work5(IP)+ENERGY(NST+IP)*FR(L)*DTHDF
         ECOS(IP) = ECOS(IP) + ENERGY(NST+IP)*COSTH(M)*DTHDF
         ESIN(IP) = ESIN(IP) + ENERGY(NST+IP)*SINTH(M)*DTHDF
cccc    endif

   22   CONTINUE
C
        DO 24 IP=1,NDATA
C
C *** WINDSEA
C
         IF (LWS(IP)) THEN
          IF (ENERGY(NST+IP).GT.EMWS(IP)) THEN
           EMWS(IP) = ENERGY(NST+IP)
           IMWS(IP) = L
          END IF
         ELSE
C
C *** SWELL
C
          IF (ENERGY(NST+IP).GT.EMSW(IP)) THEN
           EMSW(IP) = ENERGY(NST+IP)
           IMSW(IP) = L
          END IF
         END IF
C
C *** TOTAL  VALUES
C
          IF (ENERGY(NST+IP).GT.EMAX(IP)) THEN
           EMAX(IP) = ENERGY(NST+IP)
           IMAX(IP) = L
          END IF

   24   CONTINUE
C
   21  CONTINUE
   20 CONTINUE


C
C 13.3 EVALUATE SWELL VALUES FROM TOTAL AND WINDSEA
C -------------------------------------------------
C

cc make this a loop over ngrid ?? with test on IA1 lsmask value ??

cc    DO 30 IP=1,Ngrid

      idata=0

      do j=1,ngy

       do i=1,ngx

        ip=(j-1)*ngx + i

        if(.not.ia1(i,j)) then

         idata=idata+1

         if(idata.gt.ndata)then
          WRITE(6,*)'error idata gt ndata ',idata,ndata
          WRITE(6,*)'RETURNING FROM WAVEH ERROR CODE  1'
          ICODE=1
          return
         endif

         wh(ip,1) = work1(idata)
         wh(ip,4) = work4(idata)
         wh(ip,5) = work5(idata)
         wh(ip,7) = work7(idata)
         wh(ip,10)= work10(idata)
         wh(ip,11)= work11(idata)

C
C *** SWELL ENERGY
C
        WH(IP,13) = WH(IP,1) - WH(IP,7)
        IF (WH(IP,13).LT.0.0) WH(IP,13) = 0.0

CC re-use work1 to hold WH(ip,13) for later use
        work1(idata)=wh(ip,13)

C
C *** SWELL INTEGRATED VALUES FOR PERIODS
C
        WH(IP,16) = WH(IP,4) - WH(IP,10)
        IF (WH(IP,16).LT.0.0) WH(IP,16) = 0.0
        WH(IP,17) = WH(IP,5) - WH(IP,11)
        IF (WH(IP,17).LT.0.0) WH(IP,17) = 0.0
C
C *** TOTAL SEA DIRECTIONS - MEAN
C
C note for UM wave that WAM directions are from zero=North increasing
C anticlockwise ( ie wind direction is atan2(u,v) which is tan -1 (u/v)
C
C SO here we have that ECOS is the northward component and
C                      ESIN is the eastward  component
C
C SO taking ATAN2(ESIN, ECOS) gives tan -1 (esin/ecos) = tan -1 (u/v)
C SO the direction evaluated here is with ZERO=NORTH
C                                     and increasing clockwise
C
C (in radians) so only need to convert to degrees
C
         IF(WH(IP,1).GT.0.00001) THEN
          WH(IP,2) = ATAN2(ESIN(Idata),ECOS(Idata))
         END IF
C
C *** WINDSEA DIRECTIONS
C
c             ?? 7 not 8
         IF(WH(IP,7).GT.0.00001) THEN
          WH(IP,8) = ATAN2(ESIN1(Idata),ECOS1(Idata))
         END IF
C
C *** NOTE CHECKS TOTAL SWELL ENERGY NONZERO
C
         AA = WH(IP,13)
         IF (AA.GT.0.00001) THEN
C
C *** SWELL DIRECTIONS
C
         ESIN2(Idata) = ESIN(Idata)-ESIN1(Idata)
         ECOS2(Idata) = ECOS(Idata)-ECOS1(Idata)
         WH(IP,14) = ATAN2(ESIN2(Idata),ECOS2(Idata))
        END IF
C
C *** RESET WORK ARRAYS TO ZERO
C
       ECOS(idata) = 0.0
       ESIN(idata) = 0.0
       ECOS1(idata) = 0.0
       ESIN1(idata) = 0.0
       ECOS2(idata) = 0.0
       ESIN2(idata) = 0.0
C

       endif
       enddo
       enddo
cc 30 CONTINUE
C
C 13.4 CALCULATE PRINCIPAL DIRECTIONS
C ------------------------------------
C
      DO 45 M=1,nang
       DO 42 IP=1,NDATA
        LL=IMAX(IP)
        IF (LL.GT.0) THEN
         NST = ((LL-1)*nang+M-1)*NDATA
         ECOS(IP) = ECOS(IP) + ENERGY(NST+IP)*COSTH(M)
         ESIN(IP) = ESIN(IP) + ENERGY(NST+IP)*SINTH(M)
        END IF
        L1=IMWS(IP)
        IF (L1.GT.0) THEN
         NS1 = ((L1-1)*nang+M-1)*NDATA
         ECOS1(IP) = ECOS1(IP) + ENERGY(NS1+IP)*COSTH(M)
         ESIN1(IP) = ESIN1(IP) + ENERGY(NS1+IP)*SINTH(M)
        END IF
CC need to sort out this use of WH(,13)
CC WH(IP,13) is copied into data work array work1 in previous loop
CC      IF (WH(IP,13).GT.0.00001) THEN
        IF (Work1(IP).GT.0.00001) THEN
         L2=IMSW(IP)
         NS2 = ((L2-1)*nang+M-1)*NDATA
         ECOS2(IP) = ECOS2(IP) + ENERGY(NS2+IP)*COSTH(M)
         ESIN2(IP) = ESIN2(IP) + ENERGY(NS2+IP)*SINTH(M)
        END IF
   42  CONTINUE
   45 CONTINUE
C

CC MAKE THIS LOOP over ngx  ngy

      idata=0 ! must reset idata before next use

cc    DO IP=1,NDATA
      do j=1,ngy
       do i=1,ngx

        ip=(j-1)*ngx + i    ! counter for points on full grid

        if(.not.ia1(i,j)) then

         idata=idata+1      ! counter for data points only

         if(idata.gt.ndata)then
          WRITE(6,*)'error idata gt ndata ',idata,ndata
          WRITE(6,*)'setting error code and returning'
          icode=99
          return
         endif

       IF(IMAX(idata).GT.0)WH(IP,3) = ATAN2(ESIN(idata),ECOS(idata))
       IF(IMWS(idata).GT.0)WH(IP,9) = ATAN2(ESIN1(idata),ECOS1(idata))
       IF (WH(IP,13).GT.0.00001.AND.IMSW(idata).GT.0) THEN
         WH(IP,15) = ATAN2(ESIN2(idata),ECOS2(idata))
       ENDIF
C
C
C *** TOTAL SPECTRUM  - ZERO UP CROSSING PERIOD
C                     - & MEAN PERIOD
C                     - & PEAK PERIOD
C
        IF(WH(IP,4).GT.0.00001) THEN
         WH(IP,4) = SQRT(WH(IP,1)/WH(IP,4))
        ENDIF

        IF(WH(IP,5).GT.0.00001) THEN
         WH(IP,5) = WH(IP,1)/WH(IP,5)
        ENDIF

        IF(IMAX(idata).GT.0) WH(IP,6) = 1./FR(IMAX(idata))
C
C *** WIND SEA PERIODS
C
        IF (WH(IP,10).GT.0.0) THEN
         WH(IP,10) = SQRT(WH(IP,7)/WH(IP,10))
        ELSE
         WH(IP,10) = 0.0
        END IF
        IF (WH(IP,11).GT.0.0) THEN
         WH(IP,11) = WH(IP,7)/WH(IP,11)
        ELSE
         WH(IP,11) = 0.0
        ENDIF
        IF(IMWS(idata).GT.0) WH(IP,12) = 1./FR(IMWS(idata))
C
C *** NOTE CHECKS TOTAL SWELL ENERGY NONZERO
C *** SWELL PERIODS
C
         IF (WH(IP,16).GT.0.000001) THEN
          WH(IP,16) = SQRT(WH(IP,13)/WH(IP,16))
         ELSE
          WH(IP,16) = 0.0
         END IF
         IF (WH(IP,17).GT.0.000001) THEN
          WH(IP,17) = WH(IP,13)/WH(IP,17)
         ELSE
          WH(IP,17) = 0.0
         END IF
        AA = WH(IP,13)
        IF (AA.GT.0.00001) THEN
         IF(IMSW(idata).GT.0) WH(IP,18) = 1./FR(IMSW(idata))
        ELSE
         WH(IP,18) = 0.0
        END IF
C
C *** TOTAL WAVE HEIGHT
C
       WH(IP,1) = 4.*SQRT(WH(IP,1))
C
C *** WINDSEA  WAVE HEIGHT
C
       WH(IP,7) = 4.*SQRT(WH(IP,7))
C
C *** SWELL WAVE HEIGHT
C
       WH(IP,13) = 4.*SQRT(WH(IP,13))
C

c     convert directions to degrees in range 0-360

         WH(ip,2)=wh(ip,2)*RECIP_PI_OVER_180
         WH(ip,3)=wh(ip,3)*RECIP_PI_OVER_180
         WH(ip,8)=wh(ip,8)*RECIP_PI_OVER_180
         WH(ip,9)=wh(ip,9)*RECIP_PI_OVER_180
         WH(ip,14)=wh(ip,14)*RECIP_PI_OVER_180
         WH(ip,15)=wh(ip,15)*RECIP_PI_OVER_180

         if(wh(ip,2).lt.0.) wh(ip,2)=wh(ip,2)+360.
         if(wh(ip,3).lt.0.) wh(ip,3)=wh(ip,3)+360.
         if(wh(ip,8).lt.0.) wh(ip,8)=wh(ip,8)+360.
         if(wh(ip,9).lt.0.) wh(ip,9)=wh(ip,9)+360.
         if(wh(ip,14).lt.0.) wh(ip,14)=wh(ip,14)+360.
         if(wh(ip,15).lt.0.) wh(ip,15)=wh(ip,15)+360.

         if(wh(ip,2).gt.360.) wh(ip,2)=wh(ip,2)-360.
         if(wh(ip,3).gt.360.) wh(ip,3)=wh(ip,3)-360.
         if(wh(ip,8).gt.360.) wh(ip,8)=wh(ip,8)-360.
         if(wh(ip,9).gt.360.) wh(ip,9)=wh(ip,9)-360.
         if(wh(ip,14).gt.360.) wh(ip,14)=wh(ip,14)-360.
         if(wh(ip,15).gt.360.) wh(ip,15)=wh(ip,15)-360.


        endif ! check for sea points if(ia1....
       enddo
      enddo

C
       if(l_wvtra) then

cccccc   len1=nfre*nang*ndata

ccc      call wavetr(energy,pswh,perio,pdir,kwtot,fr,dfim,th,
ccc  +len1,ndata,kwtmax,nfre,nang,rmdi,icode)

CC pointers to wavetrain data in array WH
CC usage altered from UKMO model, so that now array WH has
CC all wvtrain heights together then all periods together then
CC all directions together
CC (means we don't need array pswh / perio / pdir in this routine)

       i_SWH=19
       i_PER=i_SWH + kwtmax
       i_DIR=i_PER + kwtmax
       i_NUM=i_DIR + kwtmax

       call wavetr(energy,WH(1,i_swh),WH(1,i_per),WH(1,i_dir),
     +kwtot,fr,dfim,th,
     +len1,ndata,kwtmax,nfre,nang,rmdi,icode)

C       * here add  test return code if ne 0 then error trap  *
C
C
        do ip=1,ndata
C          wave heights
ccc      wh(ip,25)=pswh(ip,1)
ccc      wh(ip,28)=pswh(ip,2)
ccc      wh(ip,31)=pswh(ip,3)

C          wave periods
ccc      wh(ip,26)=perio(ip,1)
ccc      wh(ip,29)=perio(ip,2)
ccc      wh(ip,32)=perio(ip,3)

C          wave directions
ccc      wh(ip,27)=pdir(ip,1)
ccc      wh(ip,30)=pdir(ip,2)
ccc      wh(ip,33)=pdir(ip,3)

C          number of trains
         wh(ip,i_num) = 1.0*kwtot(ip)
C
        enddo
c
       endif

      RETURN
      END
