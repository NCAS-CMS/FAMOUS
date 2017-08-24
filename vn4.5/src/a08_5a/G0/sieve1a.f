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
!    SUBROUTINE SIEVE-------------------------------------------------- 
!
! Subroutine Interface:
      SUBROUTINE SIEVE (NPNTS,R,CAN_WCNT,CAN_CPY,AREA,TFALL    
     &,                 TIMESTEP)                
 
      IMPLICIT NONE
!
! Description:
!     Calculates the flux of water passing through the canopy    
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
     &,CAN_WCNT(NPNTS)      ! IN Canopy water content (kg/m2).   
     &,CAN_CPY(NPNTS)       ! IN Canopy capacity (kg/m2).    

! Array arguments with intent(OUT):
      REAL
     & TFALL(NPNTS)         ! OUT Throughfall (kg/m2/s).    

! Local scalars:
      INTEGER
     & I                    ! WORK Loop counter.

      REAL
     & AEXP                 ! WORK Exponential term.
     &,CAN_RATIO            ! WORK Fractional saturation of the
C                           !      canopy.                   
     &,FDT_TERM             ! WORK Finite timestep term.

! Local parameters:
      REAL
     & GAMMA                ! Forward timestep weighting
      PARAMETER(GAMMA=1.0)

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
          IF (CAN_CPY(I).GT.0.0.AND.R(I).GT.0.0) THEN        
            AEXP=AREA*CAN_CPY(I)/(R(I)*TIMESTEP)                        
            AEXP=EXP(-AEXP)                                             
            CAN_RATIO=CAN_WCNT(I)/CAN_CPY(I)                      
            TFALL(I)=R(I)*((1.0-CAN_RATIO)*AEXP+CAN_RATIO)            
          ELSE                       
           TFALL(I)=R(I)          
          END IF              
        ENDDO

C-----------------------------------------------------------------------
C Bimodel canopy water, random overlap
C-----------------------------------------------------------------------
      ELSEIF (TF_MODEL.EQ.2) THEN
        DO I=1,NPNTS                                                    
          IF (CAN_CPY(I).GT.0.0.AND.R(I).GT.0.0) THEN                   
            CAN_RATIO=CAN_WCNT(I)/CAN_CPY(I)                            
            FDT_TERM=GAMMA*R(I)*TIMESTEP/CAN_CPY(I) 
            TFALL(I)=R(I)*((CAN_RATIO+FDT_TERM)/(1.0+FDT_TERM))
          ELSE                                                          
            TFALL(I)=R(I)                                               
          END IF                                                        
        ENDDO

C-----------------------------------------------------------------------
C Bimodel canopy water, maximum overlap
C-----------------------------------------------------------------------
      ELSEIF (TF_MODEL.EQ.3) THEN
        DO I=1,NPNTS                                                    
          IF (CAN_CPY(I).GT.0.0.AND.R(I).GT.0.0) THEN                   
            CAN_RATIO=CAN_WCNT(I)/(AREA*CAN_CPY(I))                     
            IF (CAN_RATIO.LT.1.0) THEN
              FDT_TERM=GAMMA*R(I)*TIMESTEP/(AREA*CAN_CPY(I)) 
              TFALL(I)=R(I)*((CAN_RATIO+FDT_TERM)/(1.0+FDT_TERM))
            ELSE
              TFALL(I)=R(I)
            ENDIF
          ELSE
            TFALL(I)=R(I)
          ENDIF
        ENDDO
      ENDIF 
      RETURN
      END
