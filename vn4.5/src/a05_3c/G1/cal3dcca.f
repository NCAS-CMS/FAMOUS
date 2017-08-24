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
!  Subroutine CALC_3D_CCA: Calculates a conv. cld amt on model levels.
!
!  Subroutine Interface:
      SUBROUTINE CALC_3D_CCA(NP_FIELD,NPNTS,NLEV,NBL
     &                      ,ANVIL_FACTOR,TOWER_FACTOR
     &                      ,AKM12,BKM12,CLOUD_BASE,CLOUD_TOP
     &                     ,FREEZE_LEV,PSTAR,CCA_2D,CCA_3D,L_CLOUD_DEEP)
!
      IMPLICIT NONE
!
! Description: Calculates a 3D convective cloud amount (i.e. on model
!             levels) from the 2D convective cloud amount array
!             according to parameters specified in the umui and the
!             position of cloud base, cloud top and freezing level.
!
! Method: The 2D convective cloud amount is expanded into the vertical
!         by applying it between the cloud base and top with the
!         additional constraints that
!         (i)   If the cloud base is in the boundary layer,
!         (ii)  cloud top is above the freezing level and
!         (iii) the cloud is more than 500mb deep
!         then the cloud below the freezing level will be multiplied
!         by TOWER_FACTOR, and the cloud above the freezing level
!         will be linearly (with model level) increased to cloud top
!         where it will be equal to the 2D fraction * ANVIL_FACTOR.
!
! Current Code Owner: Julie M. Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.4     18/9/97   Original code. J.Gregory.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Global variables (*CALLed COMDECKs etc...):
!------------------------------------------------------------------
!   Scalar arguments with intent(in):
!------------------------------------------------------------------
      INTEGER NPNTS               ! IN Number of points
     &       ,NP_FIELD            ! IN Full size of data
     &       ,NLEV                ! IN Number of levels
     &       ,NBL                 ! IN Number of Boundary layer levels
      REAL ANVIL_FACTOR           ! IN Needed in calculation of vertical
     &    ,TOWER_FACTOR           ! IN cloud amount distribution
      LOGICAL L_CLOUD_DEEP        ! IN Apply depth criterion if true
!------------------------------------------------------------------
!   Array  arguments with intent(in):
!------------------------------------------------------------------
      INTEGER CLOUD_TOP(NP_FIELD) ! IN Convective cloud top level
     &       ,CLOUD_BASE(NP_FIELD)! IN Convective cloud base level
     &       ,FREEZE_LEV(NPNTS)   ! IN Freezing level
!
      REAL PSTAR(NP_FIELD)        ! IN Surface pressure
     &    ,AKM12(NLEV+1)          ! IN Hybrid co-ord coeffs to define
     &    ,BKM12(NLEV+1)          !    pressure at level k-1/2
     &    ,CCA_2D(NPNTS)          ! IN 2D convective cloud amount
!------------------------------------------------------------------
!   Array  arguments with intent(out):
!------------------------------------------------------------------
      REAL CCA_3D(NP_FIELD,NLEV)  ! OUT Convective cloud amount on
!                                 !     model levels
!------------------------------------------------------------------
! Local parameters:
!------------------------------------------------------------------
      REAL DEEP_DP                ! Depth cloud must reach to be 'deep'
!
      PARAMETER (DEEP_DP = 50000) ! Critical depth of clouds = 500hPa
!------------------------------------------------------------------
! Local scalars:
!------------------------------------------------------------------
      INTEGER ANVIL_LEV           ! Base level of 'anvil' if it is to
!                                 ! be applied.
      INTEGER I,K                 ! Loop counters
!
      REAL ANVIL_DEPTH
     &    ,P_CLOUD_BASE
     &    ,P_CLOUD_TOP
!
      LOGICAL DEEP
!
!======================================================================
!  ANVIL CLOUD CALCULATION:
!  If cloud base is in the PBL, and cloud top is above (or at)
!  the freezing level, then add an anvil cloud by increasing the
!  cloud fraction linearly from freezing lev to cloud top. Also
!  decrease the cloud fraction below this level to represent the
!  'tower'.
!======================================================================
!
      IF (L_CLOUD_DEEP) THEN
      DO I = 1,NPNTS
        ANVIL_DEPTH = 0
        ANVIL_LEV = 0
        DEEP = .FALSE.
C----------------------------------------------------------------------
C Calculate cloud depth:
C----------------------------------------------------------------------
        IF (CCA_2D(I).GT.0.0) THEN
          P_CLOUD_BASE = AKM12(CLOUD_BASE(I)) +
     &                   BKM12(CLOUD_BASE(I))*PSTAR(I)
          P_CLOUD_TOP  = AKM12(CLOUD_TOP(I)+1) +
     &                   BKM12(CLOUD_TOP(I)+1)*PSTAR(I)
          DEEP   = (P_CLOUD_BASE - P_CLOUD_TOP) .GE. DEEP_DP
C----------------------------------------------------------------------
C Check to see if cloud is deep and above freezing level:
C----------------------------------------------------------------------
          IF ( ( CLOUD_BASE(I) .LT. NBL )  .AND.
     &         ( CLOUD_TOP(I) .GT. FREEZE_LEV(I) ) .AND.
     &         ( DEEP ) ) THEN
C----------------------------------------------------------------------
C Define anvil base level as freezing level or cloud base if above FL:
C----------------------------------------------------------------------
            ANVIL_DEPTH = ( CLOUD_TOP(I) - FREEZE_LEV(I) )
            ANVIL_LEV   = FREEZE_LEV(I)
            IF ( ANVIL_DEPTH .GT. (CLOUD_TOP(I)-CLOUD_BASE(I)) ) THEN
              ANVIL_DEPTH = CLOUD_TOP(I)-CLOUD_BASE(I)
              ANVIL_LEV   = CLOUD_BASE(I)
            ENDIF
C----------------------------------------------------------------------
C Apply wedge-shaped anvil from anvil base to cloud top:
C----------------------------------------------------------------------
            DO K = ANVIL_LEV,(CLOUD_TOP(I) - 1)
              CCA_3D(I,K) = (ANVIL_FACTOR - TOWER_FACTOR)
     &                    * CCA_2D(I)
     &                    * (K - ANVIL_LEV + 1)/ANVIL_DEPTH
     &                    + (CCA_2D(I) * TOWER_FACTOR)
              IF (CCA_3D(I,K) .GE. 1.0) THEN
                CCA_3D(I,K) = 0.99
              ENDIF
            ENDDO
C----------------------------------------------------------------------
C ...and tower below (i.e. from cloud base to anvil base):
C----------------------------------------------------------------------
            DO K = CLOUD_BASE(I),ANVIL_LEV-1
              CCA_3D(I,K) = TOWER_FACTOR * CCA_2D(I)
            ENDDO
          ELSE
C----------------------------------------------------------------------
C If cloud is not 'deep' keep old fraction, but put on model levels:
C----------------------------------------------------------------------
            DO K = CLOUD_BASE(I),(CLOUD_TOP(I) - 1)
              CCA_3D(I,K) = CCA_2D(I)
            ENDDO
          ENDIF
C----------------------------------------------------------------------
C Finally check there is no cloud below cloud base or above cloud top!
C----------------------------------------------------------------------
          DO K = 1,(CLOUD_BASE(I)-1)
            CCA_3D(I,K) = 0.0
          END DO
          DO K = CLOUD_TOP(I),NLEV
            CCA_3D(I,K) = 0.0
          END DO
        ENDIF
      ENDDO
      ELSE
        DO I = 1,NPNTS                                                  
          ANVIL_DEPTH = 0                                               
          ANVIL_LEV = 0                                                 
!---------------------------------------------------------------------- 
! Calculate cloud depth:                                                
!---------------------------------------------------------------------- 
          IF (CCA_2D(I).GT.0.0) THEN                                    
!---------------------------------------------------------------------- 
! Check to see if cloud is deep and above freezing level:               
!---------------------------------------------------------------------- 
            IF ( ( CLOUD_BASE(I) .LT. NBL )  .AND.                      
     &           ( CLOUD_TOP(I) .GT. FREEZE_LEV(I) ) ) THEN 
!---------------------------------------------------------------------- 
! Define anvil base level as freezing level or cloud base if above FL:  
!---------------------------------------------------------------------- 
              ANVIL_DEPTH = ( CLOUD_TOP(I) - FREEZE_LEV(I) )            
              ANVIL_LEV   = FREEZE_LEV(I)                               
              IF ( ANVIL_DEPTH .GT. (CLOUD_TOP(I)-CLOUD_BASE(I)) ) THEN 
                ANVIL_DEPTH = CLOUD_TOP(I)-CLOUD_BASE(I)                
                ANVIL_LEV   = CLOUD_BASE(I)                             
              ENDIF                                                     
!---------------------------------------------------------------------- 
! Apply wedge-shaped anvil from anvil base to cloud top:                
!---------------------------------------------------------------------- 
              DO K = ANVIL_LEV,(CLOUD_TOP(I) - 1)                       
                CCA_3D(I,K) = (ANVIL_FACTOR - TOWER_FACTOR)             
     &                      * CCA_2D(I)                                 
     &                      * (K - ANVIL_LEV + 1)/ANVIL_DEPTH           
     &                      + (CCA_2D(I) * TOWER_FACTOR)                
                IF (CCA_3D(I,K) .GE. 1.0) THEN                          
                  CCA_3D(I,K) = 0.99                                    
                ENDIF                                                   
              ENDDO                                                     
!---------------------------------------------------------------------- 
! ...and tower below (i.e. from cloud base to anvil base):              
!---------------------------------------------------------------------- 
              DO K = CLOUD_BASE(I),ANVIL_LEV-1                          
                CCA_3D(I,K) = TOWER_FACTOR * CCA_2D(I)                  
              ENDDO                                                     
            ELSE                                                        
!---------------------------------------------------------------------- 
! If cloud does not satisfy anvil criteria, keep old fraction, but put 
! on model levels:    
!---------------------------------------------------------------------- 
              DO K = CLOUD_BASE(I),(CLOUD_TOP(I) - 1)                   
                CCA_3D(I,K) = CCA_2D(I)                                 
              ENDDO                                                     
            ENDIF                                                       
!---------------------------------------------------------------------- 
! Finally check there is no cloud below cloud base or above cloud top!  
!---------------------------------------------------------------------- 
            DO K = 1,(CLOUD_BASE(I)-1)                                  
              CCA_3D(I,K) = 0.0                                         
            END DO                                                      
            DO K = CLOUD_TOP(I),NLEV                                    
              CCA_3D(I,K) = 0.0                                         
            END DO                                                      
          ENDIF                                                         
        ENDDO                                                           
      ENDIF
C
C
C======================================================================
C  END OF ANVIL CALCULATION
C======================================================================
C
      RETURN
      END
