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
!+ subroutine SHAPE
!
! Description:
!   subroutine to calculate JONSWAP spectral shape array for use in
!   unified wave model    wave height assimilation routines
!
! Method:
!   As in UKMO 2G operational wave model
!
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! UM4.1    June 1996 Original code.  Sophie Kelsall
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header


      subroutine shape(nfre,nang,df,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
     & jonswap_shape,
     &                 icode)
C
C
C subroutine to calculate JONSWAP spectral shape array for use is
C UKMO wave height assimilation routines
C
C
ccc   REAL jonswap_shape(24,220,nfre) ! out
      real jonswap_shape(24,220,nfre) ! wave model JONSWAP spectral
C                                     ! shape used in assimilation



      REAL shnor(24,220)              ! local
      REAL df(nfre)                   ! local

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

      integer l,ifj,igam
      real a,b,d

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


      WRITE(6,*)'subroutine shape called to set array jonswap_shape'

      mfj=220
      sigma=0.08

CCC      pi=3.1415926
      rpio2=2.0/pi

C ----------------------------------------------------------------------
CL 2.5 CREATE A LOOKUP TABLE OF JONSWAP SHAPE FUNCTIONS. A THREE
C *** PARAMETER FORM OF THE JONSWAP SPECTRUM IS USED: EJ(F,FJ,GAMMA)
C *** AND THE SHAPE FUNCTION IS NORMALISED OVER FREQUENCY
C
C *** ZERO LOOKUP TABLE
      DO L=1,nfre
       DO IFJ=1,mfj
        DO IGAM=1,24
          SHNOR(IGAM,IFJ)=0.0
          jonswap_shape(IGAM,IFJ,L)=0.0
        enddo
       enddo
      enddo


C
C *** FILL LOOKUP TABLE AND INTEGRATE OVER FREQUENCY
C *** NEW IMPROVED LOOKUP TABLE IGAM FROM 1.0 TO 3.3,220 FREQUENCIES.

      F1L = ALOG10(fr(1))
      ffac=1.8
      FN1L=ALOG10(ffac*fr(nfre-1))
      DFLOOK = (F1L-FN1L)/(MFJ-1.)
      RDFLK = 1./DFLOOK

C
C
C
C *** SET UP ARRAY OF MODEL FREQUENCY INTERVALS
      DO L=3,NFRE
       DF(L-1)=(fr(L)-fr(L-2))*0.5
      ENDDO
      DF(1)=(fr(2)-fr(1))*0.5+0.005
C *** TOP INTERVAL CONTAINS CORRECTION FACTOR BASED ON INTEGRATION OF
C *** PHILLIPS SPECTRUM ABOVE TOP FREQUENCY
      DF(NFRE)=fr(NFRE)-fr(NFRE-1)-DF(NFRE-1)*0.5+fr(NFRE)*0.25


      DO L=1,nfre
       DO IFJ=1,MFJ
        FJ = 10.**(FN1L+(IFJ-1)*DFLOOK)

        IF (fr(L).GE.FJ*0.8) THEN
         DO IGAM=1,24
          IGG = IGAM+9
          GAM = FLOAT(IGG)*0.1
          A = -1.25*(FJ/fr(L))**4
          B = 0.0
          d = ((fr(L)-FJ)/(FJ*SIGMA))**2
          IF (d.LT.50.) then
            B = ALOG(GAM)*EXP(-d*0.5)
          endif
          jonswap_shape(IGAM,IFJ,L) = EXP(A+B)/fr(L)**5
          SHNOR(IGAM,IFJ) = SHNOR(IGAM,IFJ)
     &                    + jonswap_shape(IGAM,IFJ,L)*DF(L)
         enddo
        END IF
       enddo
      enddo

C
C *** NORMALIZE OVER FREQUENCY AND MULTIPLY BY 2/PI TO ALLOW FOR LATER
C *** COS SQUARED DISTRIBUTION

      DO L=1,nfre
       DO IFJ=1,MFJ
        DO IGAM=1,24
         IF (SHNOR(IGAM,IFJ).GE.0.0000001) THEN
          jonswap_shape(IGAM,IFJ,L) = RPIO2 * jonswap_shape(IGAM,IFJ,L)
     &                                / SHNOR(IGAM,IFJ)
         END IF
        enddo
       enddo
      enddo


      return
      end
