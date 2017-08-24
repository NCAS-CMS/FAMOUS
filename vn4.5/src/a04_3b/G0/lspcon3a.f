C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!  SUBROUTINE LSPCON-----------------------------------------
!   PURPOSE: CALCULATES CONSTANTS USED IN PRECIP
!
!    Modification History from Version 4.4
!      Version    Date
!       4.4       Sept 97          New Deck         Damian Wilson
!                                                   Sue Ballard
! ----------------------------------------------------------------
           SUBROUTINE LSPCON(CX,CONSTP)
             IMPLICIT NONE
!
! LOCAL VARIABLES
             INTEGER I
! Counter to print out the contents of CX and CONSTP
             REAL TEMP
! Forms input to the GAMMAF routine which calculates gamma functions.
     &,           GBD1,GB1,GD3,GD52,GDR3,GDR4,GDR52,G1,G2,G3
! Represents the gamma function of BI+DI+1 etc.
!
! PROCEDURE CALL
             EXTERNAL GAMMAF
!
! ------------------COMDECK C_LSPDRP-----------------------------------
!
! input variables
         REAL                 !, INTENT(IN)
     &     X1R,X2R
!            Raindrops distn coeff
     &,    X1I,X2I
!            Ice distn coeff
     &,    X3I
!            Whether there is a temperature dependence to N0snow
     &,    X4I
!            m parameter in gamma distribution
     &,    AI,BI
!            Ice mass-diameter coeff
     &,    CI,DI
!            Ice fallspeed-diameter coeff
     &,    CR,DR
!            Raindrop fallspeed-diameter coeff
!
! COX/GOLDING values - 1.5 x Heymsfield fallspeed
!
      PARAMETER(X1R=8.E6,
     &          X2R=0.0,
! Drop size distribution for rain: N(D) = N0 exp(-lambda D)
! where N0 = X1R lambda^X2R
!
     &          X1I=2.E6,
     &          X2I=0.0,
     &          X3I=1.0,
     &          X4I=0.0,
! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
! where N0 = X1I exp( -X3I T[deg C]/8.18 ) lambda^X2I   and   m=X4I
!
     &          AI=6.9E-2,
     &          BI=2.0,
! Mass diameter relationship for ice:  m(D) = AI D^BI
!
     &          CI=25.2,
     &          DI=.527,
! Fall speed diameter relationship for ice: vt(D) = CI D^DI
!
     &          CR=386.8,
     &          DR=0.67)
! Fall speed diameter relationship for rain: vt(D) = CR D^DR
!
! Obtain the size of CONSTP and CX
! Sets up the size of arrays for CX and CONSTP
      REAL CX(16),CONSTP(16)
!
! OUTPUT VARIABLES
!
             CX(1)=NINT(1000.*(2.+X4I-X2I)/(BI+1-X2I+X4I))/1000.
! Used in deposition, evaporation of snow and melting calculations.
!
             CX(2)=NINT(1000.*(5.+DI-2.*X2I+2.*X4I)/
     &                  (2.*(BI+1.-X2I+X4I)))/1000.
! Used in deposition, evaporation of snow and melting calculations.
!
             CX(3)=NINT(1000.*DI/(BI+1-X2I+X4I))/1000.
! Used in fall speed and capture calculations.
!
             CX(4)=NINT(1000.*(3.+DI-X2I+X4I)/(BI+1.-X2I+X4I))/1000.
! Used in riming calculations.
!
             CX(5)=NINT(1000.*DR/(4.+DR-X2R))/1000.
! Used in capture calculations.
!
             CX(6)=NINT(1000./(X2I-X4I-1.-BI))/1000.
! Used in capture calculations.
!
             CX(7)=NINT(1000.*(3.0+DR-X2R)/(4.0+DR-X2R))/1000.
! Used in accretion calculations.
!
             CX(8)=NINT(1000.*X2I)/1000.
! Used in capture calculations.
!
             CX(9)=NINT(1000.*X2R)/1000.
! Used in capture calculations.
!
             CX(10)=NINT(1000./(4.0+DR-X2R))/1000.
! Used in capture and evaporation of rain calculations.
!
             CX(11)=NINT(1000.*((DR+5.0)/2.0-X2R))/1000.
! Used in evaporation of rain calculations.
!
             CX(12)=NINT(1000.*(2.0-X2R))/1000.
! Used in evaporation of rain calculations.
!
             CX(13)=X3I
! Used to define temperature dependence of ice particle distribution.
!
             CX(14)=3.+X4I
! Used in capture calculations.
!
             CX(15)=2.+X4I
! Used in capture calculations.
!
             CX(16)=1.+X4I
! Used in capture calculations.
!
! Define values of GBD1 etc.
             TEMP=BI+DI+1.+X4I
             CALL GAMMAF(TEMP,GBD1)
             TEMP=BI+1.+X4I
             CALL GAMMAF(TEMP,GB1)
             TEMP=3.+DI+X4I
             CALL GAMMAF(TEMP,GD3)
             TEMP=2.5+DI/2.+X4I
             CALL GAMMAF(TEMP,GD52)
             TEMP=4.+DR
             CALL GAMMAF(TEMP,GDR4)
             TEMP=3.+DR
             CALL GAMMAF(TEMP,GDR3)
             TEMP=2.5+DR/2.
             CALL GAMMAF(TEMP,GDR52)
             TEMP=1.+X4I
             CALL GAMMAF(TEMP,G1)
             TEMP=2.+X4I
             CALL GAMMAF(TEMP,G2)
             TEMP=3.+X4I
             CALL GAMMAF(TEMP,G3)
!
! Define values of CONSTP
!
             CONSTP(1)=1.0/(AI*X1I*GB1)
! Used in fallspeed, deposition, riming, capture, evap of snow and
! melting of snow calculations.
!
             CONSTP(2)=6.2832*X1R
!              6.2832 = 2 * pi
! Used in evaporation of rain calculations.
!
             CONSTP(3)=CI*GBD1/GB1
! Used in fallspeed and capture calculations.
!
             CONSTP(4)=0.7854*X1I*CI*GD3
!              0.7854 = pi / 4
! Used in riming calculations.
!
!
             CONSTP(5)=1.0*6.2832*X1I
! 6.2832 = 2 pi   1.0 represents the relative capacitance of spheres.
! Used in deposition and evaporation of snow calculations.
!
             CONSTP(6)=89.48*SQRT(CI)*GD52
!              89.48 = 0.44Sc**0.333 /dyn viscosity**0.5
! Used in deposition, evaporation of snow and melting calculations.
!
             CONSTP(7)=4.57E-7*X1I
!              4.57E-7 = 2 pi ka / Lf
! Used in melting of snow calculations.
!
             CONSTP(8)=523.6*X1R*GDR4*CR
!              523.6 = pi rho(water)/ 6
! Used in capture, evap of rain and melting of snow calculations.
!
             CONSTP(9)=9869.6*X1I*X1R
!              9869.6 = pi**2 rho(water)
! Used in capture calculations.
!
             CONSTP(10)=0.7854*X1R*CR*GDR3
!              0.7854 = pi / 4
! Used in accretion calculations.
!
             CONSTP(11)=CR*GDR4
! Used in capture calculations.
!
             CONSTP(12)=63.100*GDR52*SQRT(CR)
!              63.1 = 0.31Sc**0.333/dyn viscosity**0.5
! Used in evaporation of rain calculations.
!
             CONSTP(13)=G2
! Used in deposition, evap of snow and melting of snow calculations.
!
             CONSTP(14)=G3*0.25
! Used in capture calulations.
!
             CONSTP(15)=G2*2.
! Used in capture calulations.
!
             CONSTP(16)=G1*5.
! Used in capture calulations.
!
! End the subroutine
             RETURN
           END
!
!  SUBROUTINE GAMMAF-------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
!   A POLYNOMIAL APPROXIMATION
! ----------------------------------------------------------------
         SUBROUTINE GAMMAF(Y,GAM)
           IMPLICIT NONE
           REAL               !, INTENT(IN)
     &       Y
           REAL               !, INTENT(OUT)
     &       GAM
! Gamma function of Y
!
! LOCAL VARIABLE
           INTEGER I,M
           REAL GG,G,PARE,X
! --------------------------------------------------------------------
           GG=1.
           M=Y
           X=Y-M
           IF (M.NE.1) THEN
             DO I=1,M-1
               G=Y-I
               GG=GG*G
             END DO
           END IF
           PARE=-0.5748646*X+0.9512363*X*X-0.6998588*X*X*X
     &     +0.4245549*X*X*X*X-0.1010678*X*X*X*X*X+1.
           GAM=PARE*GG
           RETURN
         END
!
