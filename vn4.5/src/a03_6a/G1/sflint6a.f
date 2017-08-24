C ******************************COPYRIGHT******************************
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
C
!!!  SUBROUTINE SFL_INT------------------------------------------------
!!!
!!!  Purpose: To calculate interpolation coefficients for 10m winds
!!!           and 1.5m temperature/specific humidity diagnostics.
!!!
!!!  Suitable for single column use (via *IF definition IBM).
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.4  09/05/97   New exact formulation based on Monin-Obukhov
!!!                   stability functions.
!!!                                                    R.N.B.Smith
!!!
!!!  Programming standard:
!!!
!!!  Logical component covered: Part of P243.
!!!
!!!  System Task:
!!!
!!!  External Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!*L  Arguments :-
      SUBROUTINE SFL_INT (
     & P_POINTS,P_FIELD,P1
     &,Z0M,Z0H,CD,CH
     &,Z0M_STD,CD_STD
     &,RESFT,RECIP_L_MO,V_S,V_S_STD
     &,CDR10M,CHR1P5M,CER1P5M
     +,SU10,SV10,ST1P5,SQ1P5,LTIMER
     +)
      IMPLICIT NONE

      INTEGER
     & P_POINTS          ! IN No. of P-grid points to be processed.
     &,P_FIELD           ! IN Total No. of P-grid points.
     &,P1                ! IN First P-grid point to be processed.

      REAL
     + Z0M(P_FIELD)      ! IN Roughness length for momentum (m).
     +,Z0H(P_FIELD)      ! IN Roughness length for heat and
!                        !    moisture (m).
     &,Z0M_STD(P_FIELD)  ! IN Roughness length for momentum without 
!                        !    orographic component (m).
     &,CD(P_FIELD)       ! IN Surface drag coefficient.
     &,CH(P_FIELD)       ! IN Surface transfer coefficient for heat and
!                        !    moisture.
     &,CD_STD(P_FIELD)   ! IN Surface drag coefficient excluding 
!                        !    orographic from drag.
     +,RESFT(P_FIELD)    ! IN Total resistance factor for moisture
!                        !    transfer from the surface.
     &,RECIP_L_MO(P_FIELD)
!                        ! IN Reciprocal of the Monin-Obukhov length (m)
     &,V_S(P_FIELD)      ! IN Surface layer scaling velocity including
!                        !    orographic form drag (m/s).
     &,V_S_STD(P_FIELD)  ! IN Surface layer scaling velocity excluding
!                        !    orographic form drag (m/s).

      LOGICAL
     + SU10                      ! IN 10m U-wind diagnostic flag
     +,SV10                      ! IN 10m V-wind diagnostic flag
     +,ST1P5                     ! IN screen temp diagnostic flag
     +,SQ1P5                     ! IN screen specific humidity
!                                !    diagnostic flag
     +,LTIMER                    ! IN TIMER diagnostics flag
! Output variables
!
      REAL
     + CDR10M(P_FIELD)   ! OUT interpolation coefficicent for 10m wind
     +,CHR1P5M(P_FIELD)  ! OUT Interpolation coefficient for 1.5m
!                        !     temperature
     +,CER1P5M(P_FIELD)  ! OUT Interpolation coefficient for 1.5m
!                        !     specific humidity
!*
!*L---------------------------------------------------------------------
      EXTERNAL TIMER , PHI_M_H
!*
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL Z_OBS_TQ,Z_OBS_WIND
      PARAMETER (
     + Z_OBS_TQ = 1.5    ! Height of screen observations of temperature
!                        ! and humidity.
     +,Z_OBS_WIND = 10.0 ! Height of surface wind observations.
     +)
      LOGICAL EFF_INT
      PARAMETER (EFF_INT = .FALSE.)
!
!  Define local storage.
!
!  (a) Local work arrays.
!
      REAL
     & Z_WIND(P_FIELD)     ! Height of wind observations.
     &,Z_TEMP(P_FIELD)     ! Height of temperature and humidity
!                          ! observations.
     &,PHI_M_OBS(P_FIELD)  ! Monin-Obukhov stability function for
!                          ! momentum integrated to the wind observation
!                          ! height.
     &,PHI_H_OBS(P_FIELD)  ! Monin-Obukhov stability function for
!                          ! scalars integrated to their observation
!                          ! height.
      LOGICAL
     & L_D_ARRAY(P_FIELD)
!
!  (b) Scalars.
!
      INTEGER
     + I       ! Loop counter (horizontal field index).
      LOGICAL
     & L_DUMMY
!*
      IF (LTIMER) THEN
        CALL TIMER('SFL_INT   ',3)
      ENDIF
!
!-----------------------------------------------------------------------
!! 1. If diagnostics required calculate M-O stability functions at
!!    observation heights.
!-----------------------------------------------------------------------

      IF (SU10 .OR. SV10 .OR. ST1P5 .OR. SQ1P5) THEN
        L_DUMMY = .FALSE.
        DO I=P1,P1+P_POINTS-1
          Z_WIND(I) = Z_OBS_WIND
          Z_TEMP(I) = Z_OBS_TQ
          L_D_ARRAY(I) = .TRUE.
        ENDDO
        CALL PHI_M_H (P_POINTS,P_FIELD,P1,L_DUMMY,L_D_ARRAY,
     &                RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,
     &                PHI_M_OBS,PHI_H_OBS,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
!! 2. If diagnostics required calculate interpolation coefficient
!!    for 1.5m screen temperature and specific humidity.                
!-----------------------------------------------------------------------
!
      IF (ST1P5 .OR. SQ1P5) THEN                                        
        DO I=P1,P1+P_POINTS-1
          CHR1P5M(I) = CH(I) * PHI_H_OBS(I)/(VKMAN*V_S_STD(I))          
          CER1P5M(I) = ( CHR1P5M(I) - 1.0 ) * RESFT(I)                  
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!! 3. If diagnostics required calculate interpolation coefficient
!!    for 10m winds.                                                    
!-----------------------------------------------------------------------
!
      IF ( (SU10 .OR. SV10) .AND. EFF_INT ) THEN
        DO I=P1,P1+P_POINTS-1
          CDR10M(I) = CD(I) * PHI_M_OBS(I)/(VKMAN*V_S(I))               
        ENDDO                                                           
      ELSEIF ( (SU10 .OR. SV10) .AND. .NOT.EFF_INT ) THEN
        CALL PHI_M_H (P_POINTS,P_FIELD,P1,L_DUMMY,L_D_ARRAY,
     &                RECIP_L_MO,Z_WIND,Z_TEMP,Z0M_STD,Z0H,   
     &                PHI_M_OBS,PHI_H_OBS,LTIMER)
        DO I=P1,P1+P_POINTS-1
          CDR10M(I) = CD_STD(I) * PHI_M_OBS(I)/(VKMAN*V_S_STD(I))
        ENDDO
      ENDIF
!
      IF (LTIMER) THEN
        CALL TIMER('SFL_INT ',4)
      ENDIF
      RETURN
      END
