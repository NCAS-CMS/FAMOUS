C ******************************COPYRIGHT******************************
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
!    SUBROUTINE FRUNOFF------------------------------------------------ 
!
! Subroutine Interface:
      SUBROUTINE FRUNOFF (NPNTS,R,TFALL,CAN_WCNT,CAN_CPY,INFIL,AREA
     &,                   RUNOFF,TIMESTEP)                              
 
      IMPLICIT NONE
!
! Description:
!     Calculates the fast (surface) runoff at the soil surface    
!
! Documentation:                      
!
! Current Code Owner: Peter Cox
!
! History:
! Version  Date    Comment
! -------  ----    -------
! 4.5      6/98    Options to interpret canopy moisture as bimodally
!                  rather than uniformily distributed, with either a
!                  random or maximum overlap of the wet canopy and 
!                  the precipitating area (P.M.Cox)
! 4.5  01/10/98    Removed old section-version defs. K Rogers
!
! Code Description:                                                     
!   Language: FORTRAN 77 + common extensions.                           
!                                                                       
! System component covered: P25                                         
! System Task: P25                                                      
!                                                                       
! Subroutine arguments
! Scalar arguments with intent(IN):
      INTEGER                                                           
     & NPNTS                ! IN Number of gridpoints.                  
                                                                        
      REAL                                                              
     & TIMESTEP             ! IN Model timestep (s).                    
     &,AREA                 ! IN Fractional area of the gridbox over 
C                           !    which water falls.

! Array arguments with intent(IN):
      REAL
     & R(NPNTS)             ! IN Flux of water incident on the 
C                           !    canopy (kg/m2/s).    
     &,TFALL(NPNTS)         ! IN Throughfall (kg/m2/s).    
     &,CAN_WCNT(NPNTS)      ! IN Canopy water content (kg/m2).   
     &,CAN_CPY(NPNTS)       ! IN Canopy capacity (kg/m2).    
     &,INFIL(NPNTS)         ! IN Maximum infiltration rate (kg/m2/s).

! Array arguments with intent(OUT):
      REAL
     & RUNOFF(NPNTS)        ! OUT Surface runoff (kg/m2/s).

! Local scalars:
      INTEGER
     & I                    ! WORK Loop counter.

      REAL
     & AEXP,AEXP1,AEXP2     ! WORK Exponential terms.
     &,CAN_RATIO            ! WORK Fractional saturation of the
C                           !      canopy.                   
     &,CM                   ! WORK (CAN_CPY - CAN_WCNT)/TIMESTEP        
C                           !      (kg/m2/s).    

! Local parameters:
      INTEGER
     & CAN_MODEL            ! 1 for no thermal canopy (pre 4.5 UM).
C                           ! 2 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature.
C                           ! 3 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature, plus canopy thermal
C                           !   capacity.
     &,REX_MODEL            ! 1 for uniform root density profile
C                           !   (pre 4.5 UM) with MOSES I
C                           !   rootdepths.
C                           ! 2 for exponential root density
C                           !   profile with "old" rootdepths.
     &,TF_MODEL             ! 1 for uniformily distributed canopy
C                           !   water (pre 4.5 UM).
C                           ! 2 for bimodally distributed canopy
C                           !   water with random overlap.
C                           ! 3 for bimodally distributed canopy
C                           !   water with maximum overlap.
C
C For pre 4.5 UM MOSES I choose:
C     PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)
C
      PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)

C-----------------------------------------------------------------------
C Uniform canopy water
C-----------------------------------------------------------------------
      IF (TF_MODEL.EQ.1) THEN
        DO I=1,NPNTS                                
          IF (R(I).GT.0.0) THEN                  
            IF (INFIL(I)*TIMESTEP.LE.CAN_WCNT(I).AND.                   
     &               CAN_CPY(I).GT.0.0) THEN                            
C-----------------------------------------------------------------------
C Equation (P252.14A)
C-----------------------------------------------------------------------
              AEXP=AREA*CAN_CPY(I)/R(I)                                 
              IF (CAN_WCNT(I).GT.0.0) THEN                              
                AEXP1=EXP(-AEXP*INFIL(I)/CAN_WCNT(I))                   
              ELSE                                                    
                AEXP1=0.0         
              ENDIF                    
              AEXP2=EXP(-AEXP/TIMESTEP)                             
              CAN_RATIO=CAN_WCNT(I)/CAN_CPY(I)        
              RUNOFF(I)=R(I)*(CAN_RATIO*AEXP1+(1.0-CAN_RATIO)*AEXP2)    
            ELSE                                        
C-----------------------------------------------------------------------
C Equation (P254.14B)
C-----------------------------------------------------------------------
              CM=(CAN_CPY(I)-CAN_WCNT(I))/TIMESTEP                      
              AEXP=EXP(-AREA*(INFIL(I)+CM)/R(I))                        
              RUNOFF(I)=R(I)*AEXP                             
            ENDIF                                                       
          ELSE                                                          
            RUNOFF(I)=R(I)                       
          ENDIF        
        ENDDO

C-----------------------------------------------------------------------
C Bimodal canopy water
C-----------------------------------------------------------------------
      ELSEIF (TF_MODEL.EQ.2.OR.TF_MODEL.EQ.3) THEN
        DO I=1,NPNTS                                
          IF (R(I).GT.0.0) THEN
            RUNOFF(I)=TFALL(I)*EXP(-AREA*INFIL(I)/R(I)) 
          ELSE
            RUNOFF(I)=R(I)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
