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
C
CLL    Subroutine LWCLD.
CLL
CLL      Purpose
CLL      ~~~~~~~
CLL  It calculates, from the fractional cloud cover, cloud water paths,
CLL  ice and water bulk absorption co-efficients and the fraction of the
CLL  cloud to be frozen, the effective cloud amount (cloud amount times
CLL  emissivity), for each layer and each longwave band, and returns
CLL  1-this, the "effective clear fraction", in ECA
CLL  for use in the overlap calculations when LWMAST constructs the
CLL  longwave fluxes.  It is a separate routine to make it as easy as
CLL  possible to change the cloud generation scheme & the way clouds are
CLL  passed into P232.
CLL
CLL          Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
!LL   4.0   28/9/95  FOCWWIL COMDECK now subroutine CALL. (A Bushell)
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions replaced by T3E lib functions.
CLL                       S.J.Swarbrick
!     4.4   09/4/97  Allow replacement of FOCWWIL parametrization by
!                    direct ratio of prognostic cloud ice to liquid
!                    in layer cloud calculations.   (A C Bushell)
CLL
CLL     Standard
CLL     ~~~~~~~~
CLL  If UPDATE *DEF CRAY is off, the code is standard FORTRAN 77 except
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL  for having ! comments (it then sets the "vector length" to 1) but
CLL  otherwise it includes automatic arrays also.
CLL
CLL  It is part of component P232 (longwave radiation), in task P23
CLL  (radiation).
CLL
CLL  Offline documentation is in UMDP 23.
C*L
      SUBROUTINE LWCLD (LCA, LCCWC1, LCCWC2, CCA, CCCWP, CCB, CCT,
     &     TAC, PSTAR, AB, BB,
     &     L_CLOUD_WATER_PARTITION,
     &     L1, NLEVS, NCLDS,
     &     L2,                                                          
     &     ECA)
C*
      EXTERNAL LSP_FOCWWIL
!     EITHER
!       Use temperature dependent focwwil for convection but calculate
!       ratio in layer cloud from prognostic cloud ice produced as part
!       of large-scale precipitation scheme 3A, OR
!       Use the subroutine LSP_FOCWWIL (from Section 4) consistently to
C     ! derive cloud radiative properties and precipitation amount,
C     ! taking into account that cloud does not freeze as soon as it is
C     ! cooled below the freezing point of bulk water.  The release of
C     ! latent heat of fusion (not a major term) is done differently in
C     ! order to allow energy conservation (UMDP 29).  This is the
C     ! reason for two layer cloud water contents being passed in and
C     ! then combined and differently split.
C     !  Array dimensions must be constants in FORTRAN:
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
!  "ECMWF" bands as described by Morcrette et al (J.-J. Morcrette,
!   L.D. Smith & Y. Fouquart, 1986, Beitr. Phys. Atmosph., 59, 455-469)
C*L
      INTEGER!, INTENT(IN) ::
     &     L1,                    ! Full field dimension
     &     L2,                    ! Number of points to be treated      
     &     NLEVS,                 ! Number of model levels
     &     NCLDS                  ! Number of possibly cloudy levels
      REAL!, INTENT(IN) ::
     &     LCA(L1,NCLDS),         ! Layer cloud fractional cover
     &     LCCWC1(L1,NCLDS),      ! layer cloud condensed water content
     &     LCCWC2(L1,NCLDS),      ! layer cloud condensed ice content
C     !   These are specific cloud water contents, mass per unit mass,
C     !               and, as explained above, only their sum is used.
     &     CCA(L1),               ! Convective cloud fractional cover
     &     CCCWP(L1),             !             and condensed water path
     &     TAC(L1,NLEVS),         ! Mid-layer atmospheric temperature
     &     PSTAR(L1),             ! Surface pressure
     &     AB(NLEVS+1),BB(NLEVS+1)! As & Bs at layer boundaries
C     ! Note that the fractional cover is that given by P29 without
C     ! any knowledge of convective cloud, and in the layers though
C     ! which the convective cloud extends it is taken to be the
C     ! fractional cover by layer cloud not over the whole grid-box but
C     ! the parts outside the convective cloud.
C     !  The LCCWC are averages over the whole grid-box, while CCCWP is
C     ! the in-cloud value.
      INTEGER!, INTENT(IN) ::
     &     CCB(L1), CCT(L1)       ! Convective cloud base and top
      LOGICAL!, INTENT(IN) ::
     &    L_CLOUD_WATER_PARTITION ! True if prognostic cloud ice used
      REAL!, INTENT(OUT) ::
     &     ECA(L2,NCLDS,NBANDS)   ! "effective clear amount"
C*
! It has one array CECC, L2 in size, of dynamic storage.           
! Its structure consists of one nested set of loops to find the    
! effective cloud contributions from layer cloud, and one for      
! convective cloud.                                                
C
C
      REAL CECL,                  ! Contribution to effective cloud
     &     CECC(L2),              !     from layer and convective cloud
     &     ABSCCL,                ! Mean absorption coeff. for a cloud
     &     LQFR,                  !  Liquid fraction of layer
     &     CCLQFR(L2),            !               & convective cloud
     &     DPBYG,                 ! Converts mixing ratio to pathlength
     &     DAB, DBB,              ! Differences of As & Bs
     &     EXPONC(L1,NCLDS),      ! Exponent calculating emissivity     
     &     EXPONB(L1,NBANDS)
      INTEGER BAND, LEVEL, J      ! Loopers over band, level and points
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL ABSIW(NBANDS,2)   ! Absorption co-efficients for ice & water
      DATA ABSIW / NBANDS*-65., NBANDS*-130. /
C                                                                       
      REAL expEXPONC(L2,NCLDS)
      REAL expEXPONB(L2,NBANDS)
C
C
CL    ! Section 1
CL    ! ~~~~~~~~~
CL    !  First find contributions to ECA, ECTA & ECBA from layer cloud:
C
      IF (.NOT. L_CLOUD_WATER_PARTITION) THEN
!        Use end of ECA as workspace for liquid fraction for each layer:
      DO LEVEL=1, NCLDS
        CALL LSP_FOCWWIL (TAC(1,LEVEL), L2, ECA(1,LEVEL,NBANDS))
      END DO
      END IF
C
      DO BAND=1, NBANDS                                                
       DO LEVEL=1, NCLDS                                               
        DAB = AB(LEVEL) - AB(LEVEL+1)
        DBB = BB(LEVEL) - BB(LEVEL+1)
Cfpp$   Select(CONCUR)
        DO J=1, L2                                                      
C        !  From the liquid fraction find the average absorption coefft:
          IF (L_CLOUD_WATER_PARTITION) THEN
!         calculate liquid fraction focwwil as ratio qcl/(qcl+qcf)
            IF (LCA(J,LEVEL).GT.0.) THEN
              LQFR = LCCWC1(J,LEVEL)/ (LCCWC1(J,LEVEL)+LCCWC2(J,LEVEL))
            ELSE
!           Arbitrary number: makes it safe & vectorizable
              LQFR = 0.
            ENDIF
          ELSE
!         set proportion of liquid water focwwil from lsp_focwwil
         LQFR = ECA(J,LEVEL,NBANDS)
          ENDIF
!
         ABSCCL = ( 1. - LQFR ) * ABSIW(BAND,1) + LQFR * ABSIW(BAND,2)
C        !  Calculate cloud water path and convert to in-cloud mean:
         DPBYG = ( DAB + PSTAR(J) * DBB ) / G
         EXPONC(J,LEVEL) =
     &     ABSCCL * ( LCCWC1(J,LEVEL) + LCCWC2(J,LEVEL) ) * DPBYG    
         IF ( LCA(J,LEVEL) .NE. 0. ) 
     &        EXPONC(J,LEVEL) = EXPONC(J,LEVEL) / LCA(J,LEVEL)      
        end do                                                          
       end do                                                           
      end do   
      DO J=1,L2
        DO LEVEL=1, NCLDS
          expEXPONC(J,LEVEL)=exp(EXPONC(J,LEVEL))
        END DO
      END DO
C        ! Equation 2.3.1:
      DO BAND=1, NBANDS                                             
       DO LEVEL=1, NCLDS                                            
        DO J=1, L2                                                  
         CECL = LCA(J,LEVEL) * (1. - expEXPONC(J,LEVEL))                
         IF (  LEVEL .GE. CCB(J)  .AND.  LEVEL .LT. CCT(J)  )
     &             CECL = CECL * ( 1. - CCA(J) )
         ECA(J,LEVEL,BAND) = CECL
        end do                                                          
       end do                                                           
      end do                                                    
C
CL    ! Section 2
CL    ! ~~~~~~~~~
CL    !  And then convective cloud contributes similarly,
C     !  except that the temperature from which the ice/water fraction
C     !  is calculated is that of the top layer into which it extends,
C     !  CCCWP does not need to be converted into an in-cloud value,
C     !  and that the effective cloud cover then has to be put into a
C     !  range of levels.
C     !
      DO J=1, L2
        CCLQFR(J) = TAC(J,CCT(J))
      END DO
      CALL LSP_FOCWWIL (CCLQFR, L2, CCLQFR)
C
      DO BAND=1, NBANDS           
Cfpp$  Select(CONCUR)
       DO J=1, L2                
        ABSCCL = (1.-CCLQFR(J))*ABSIW(BAND,1) + CCLQFR(J)*ABSIW(BAND,2)
        EXPONB(J,BAND) = ABSCCL * CCCWP(J)                              
       end do                                                           
      end do
C        ! Equation 2.3.1 again:
      DO J=1,L2
        DO BAND=1, NBANDS  
          expEXPONB(J,BAND)=exp(EXPONB(J,BAND))
        END DO
      END DO
C
      DO BAND=1, NBANDS                                             
       DO J=1, L2                                                   
        CECC(J) = CCA(J) * (1. - expEXPONB(j,band)) 
       end do                
C      ! The asymmetry beween CCB and CCT is because the indexing of the
C      !  effective cloud arrays is set up to simplify the layer cloud
C      !  loop, with a top being indexed like a bottom a layer below.
       DO LEVEL=1, NCLDS        
Cfpp$   Select(CONCUR)
        DO J=1, L2          
         IF (  LEVEL .GE. CCB(J)  .AND.  LEVEL .LT. CCT(J)  )
     &             ECA(J,LEVEL,BAND) = ECA(J,LEVEL,BAND) + CECC(J)
        end do
       end do
      end do    
C     !
CL    ! Section 3
CL    ! ~~~~~~~~~
CL    !  Finally change ECA from effective cloud amount to effective
CL    !                                                    clear amount
      DO BAND=1, NBANDS                                                 
       DO LEVEL=1, NCLDS                                              
        DO J=1, L2                                                     
         ECA(J,LEVEL,BAND) = 1. - ECA(J,LEVEL,BAND)
        end do                                                          
       end do                                                           
      end do  
C
      RETURN
      END
