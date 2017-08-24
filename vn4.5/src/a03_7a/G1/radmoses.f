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
!!!   SUBROUTINE RAD_MOSES--------------------------------------------- 
!!!                                                                     
!!!  Purpose: Calculate surface radiation fluxes over snow-free land    
!!!           and snow, and adjust atmospheric heating rates between    
!!!           radiation timesteps for MOSES II land surface scheme.     
!!!                                                                     
!!!  Model            Modification history:                             
!!! version  Date                                                       
!!!  4.4     8/97   New deck    Richard Essery                          
!!!                                                                     
!!!---------------------------------------------------------------------
                                                                        
      SUBROUTINE RAD_MOSES (                                            
     & P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,NLEVS,BL_LEVELS,
     & AKH,BKH,COS_Z,ALBEDO,ALBSNF,LW_SURF,LW_DT,PSTAR,SW_SURF_NCZ,     
     & SNOW_FRAC,FRAC,TSTAR_RAD,TSTAR_TILE,TIMESTEP,
     & T,RAD_NO_SNOW,RAD_SNOW                                           
     & )                                                                
                                                                        
      IMPLICIT NONE                                                     
                                                                        
      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice
                                                                        
      INTEGER                                                           
     & P_FIELD                     ! IN Total number of P-grid points.  
     &,LAND_FIELD                  ! IN Total number of land points.    
     &,LAND_PTS                    ! IN Number of land points processed.
     &,LAND1                       ! IN First land point to be processed
     &,LAND_INDEX(LAND_FIELD)      ! IN Index of land points.
     &,NLEVS                       ! IN Number of atmospheric levels.   
     &,BL_LEVELS                   ! IN Number of boundary layer levels.
                                                                        
      REAL                                                              
     & AKH(NLEVS+1)                ! IN Hybrid 'A' for layer interfaces.
     &,BKH(NLEVS+1)                ! IN Hybrid 'B' for layer interfaces.
     &,COS_Z(P_FIELD)              ! IN cos ( zenith angle )            
     &,ALBEDO(P_FIELD)             ! IN Gridbox-mean surface albedo.    
     &,ALBSNF(LAND_FIELD)          ! IN Snow-free surface albedo.       
     &,LW_SURF(P_FIELD)            ! IN Net surface LW flux (W/m2).     
     &,LW_DT(P_FIELD,NLEVS)        ! IN LW atmospheric heating rates    
!                                  !    (K/timestep).                   
     &,PSTAR(P_FIELD)              ! IN Surface pressure (Pa).          
     &,SW_SURF_NCZ(P_FIELD)        ! IN Net surface SW flux divided by  
!                                  !    COS_Z (W/m2).                   
     &,SNOW_FRAC(LAND_FIELD)       ! IN Snow-cover fraction.            
     &,FRAC(LAND_FIELD,NTYPE)      ! IN Tile fractions.
     &,TSTAR_RAD(P_FIELD)          ! IN Effective radiative surface     
!                                  !    temperature on the last         
!                                  !    radiation timestep (K).         
     &,TSTAR_TILE(LAND_FIELD,NTYPE)! IN Surface tile temperatures (K).  
     &,TIMESTEP                    ! IN Timesetep (s).                  
                                                                        
      REAL                                                              
     & T(P_FIELD,NLEVS)            ! INOUT Atmospheric temperatures (K).
                                                                        
      REAL                                                              
     & RAD_NO_SNOW(P_FIELD)        ! OUT Net surface radiation over     
!                                  !     snow-free land (W/m2).         
     &,RAD_SNOW(P_FIELD)           ! OUT Net surface radiation over     
!                                  !     snow or land-ice (W/m2).       
                                                                        
! Workspace                                                             
      REAL                                                              
     & SDT(LAND_FIELD)             ! Sum of absolute radiative heating  
!                                  ! rates.                             
     &,TSTAR_RAD_NOW(LAND_FIELD)   ! TSTAR_RAD on this timestep.
                                                                        
      REAL                                                              
     & DACON                       ! Used in calculation of heating     
     &,DBCON                       ! rates.                             
                                                                        
      INTEGER                                                           
     & I                           ! Horizontal field index.            
     &,K                           ! Vertical level index.              
     &,L                           ! Land field index.                  
     &,N                           ! Tile index.                        
                                                                        
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

       REAL
     & SBCON                       ! Stefan-Boltzmann constant
!                                  ! (W/m**2/K**4).
      PARAMETER ( SBCON=5.67E-8 )
                                                                        
!---------------------------------------------------------------------- 
! Add shortwave contribution to net surface radiation fluxes            
!---------------------------------------------------------------------- 
! Snow-free land tiles                                                  
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)                                               
        RAD_NO_SNOW(I) = (1. - ALBSNF(L)) * SW_SURF_NCZ(I) * COS_Z(I) / 
     &                                                 (1. - ALBEDO(I)) 
      ENDDO                                                             
                                                                        
! Snow and land-ice                                                     
      N = NTYPE                                                         
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)                                               
        IF ( SNOW_FRAC(L) .GT. 0. ) THEN
          RAD_SNOW(I) = ( SW_SURF_NCZ(I) * COS_Z(I) -
     &                      (1. - SNOW_FRAC(L))*RAD_NO_SNOW(I) )        
     &                                                   / SNOW_FRAC(L) 
        ENDIF
      ENDDO                                                             
                                                                        
!---------------------------------------------------------------------- 
! Add longwave contribution to net surface radiation fluxes             
!---------------------------------------------------------------------- 
                                                                        
! Snow-free land tiles                                                  
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)                                               
        RAD_NO_SNOW(I) = RAD_NO_SNOW(I) + LW_SURF(I)                    
     &                                          + SBCON*TSTAR_RAD(I)**4 
        TSTAR_RAD_NOW(L) = 0.                                           
      ENDDO                                                             
      DO N=1,NTYPE-1                                                    
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)                                             
          RAD_NO_SNOW(I) = RAD_NO_SNOW(I) -                             
     &                             FRAC(L,N)*SBCON*TSTAR_TILE(L,N)**4
          TSTAR_RAD_NOW(L) = TSTAR_RAD_NOW(L) + (1. - SNOW_FRAC(L)) *   
     &                                   FRAC(L,N)*TSTAR_TILE(L,N)**4
        ENDDO                                                           
      ENDDO                                                             
                                                                        
! Snow and land-ice                                                     
      N = NTYPE                                                         
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)                                               
        IF ( SNOW_FRAC(L) .GT. 0. ) THEN
          RAD_SNOW(I) = RAD_SNOW(I) + LW_SURF(I) +
     &                  SBCON*(TSTAR_RAD(I)**4 - TSTAR_TILE(L,N)**4)
        ENDIF
        TSTAR_RAD_NOW(L) = TSTAR_RAD_NOW(L) +                           
     &                                  SNOW_FRAC(L)*TSTAR_TILE(L,N)**4 
      ENDDO                                                             
                                                                        
!---------------------------------------------------------------------- 
! Adjust temperatures on boundary layer levels for changes in surface   
! temperature and LW heating rates between radiation timesteps          
!---------------------------------------------------------------------- 
                                                                        
      DO L=LAND1,LAND1+LAND_PTS-1
        TSTAR_RAD_NOW(L) = TSTAR_RAD_NOW(L)**0.25                       
      ENDDO                                                             
                                                                        
      DO K=1,BL_LEVELS                                                  
        DACON = (AKH(K) - AKH(K+1))*CP / (G*TIMESTEP)                   
        DBCON = (BKH(K) - BKH(K+1))*CP / (G*TIMESTEP)                   
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)                                             
          T(I,K) = T(I,K) -                                             
     &                SBCON*(TSTAR_RAD(I)**4 - TSTAR_RAD_NOW(L)**4) 
     &                       / (BL_LEVELS*(DACON + PSTAR(I)*DBCON))
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
