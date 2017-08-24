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
CLL  SUBROUTINE  LWPTSC
CLL
CLL      PURPOSE
CLL  It calculates scaled pathlengths of each gaseous absorber for each
CLL  layer and returns them in DPATH for use by LWMAST, which sums them
CLL  to get the total scaled pathlengths between each pair of layers so
CLL  that the gaseous transmissivities can be calculated.
CLL  Used in version 1B (gaseous effects treated as Morcrette et al,
CLL  1986) of the UM LW code.
CLL  If UPDATE *DEF CRAY is off, a version is produced which except
CLL  for the addition of ! comments is standard FORTRAN 77 (and which
CLL  sets the "vector length" to 1) but the standard version includes
CLL  CRAY automatic arrays also.
CLL  Version 3, part of the alternative code giving ECMWF-like treatment
CLL  of LW gaseous transmissivities.  Almost all the calculations are
CLL  changed: gases are generally scaled differently in each band they
CLL  have an effect in, pressure scaling is no longer by fractional
CLL  powers, quite elaborate pathlength-dependent temperature scaling
CLL  is used, a diffusivity factor is included, zero pathlengths are
CLL  permissible, and the indentation is changed.
CLL  The increased complexity of the scaling means that the loop finding
CLL  the water vapour scaled pathlengths does not (with current CRAY
CLL  compilers) vectorize without compiler option "-o aggress".
CLL  Version 3 of LWPTSC was set up from version 2.2 to be part of
CLL  version 1B (ECMWF-like gaseous transmissivities) of the LW from
CLL  release 2.7 of the UM.                William Ingram 22 June 1992
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions removed.
CLL                       S.J.Swarbrick
CLL   4.3    Feb. 97  T3E optimisation: code restructured, cray vector
CLL                    library functions introduced.
CLL                       D.Salmond & S.J.Swarbrick
CLL  4.4  20/06/97  Add missing array indices to pbypr, dco2
CLL                 & do3 in dry levels loop.   RTHBarnes.
CLL
CLL
CLL  It conforms to standard A of UMDP 4 (version 2, 18/1/90), and
CLL  includes no 8X-deprecated features.
CLL
CLL  It is part of component P232 (longwave radiation) which is in task
CLL  P23 (radiation).
CLL
CLL  External documentation is in UMDP 23.
C*L
      SUBROUTINE LWPTSC (H2O, CO2, O3, PSTAR, AC, BC, AB, BB, TAC,
     &     L2,
     &     NWET, NOZONE, NLEVS, L1,DPATH)   
C*
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=12)
C     ! This set is for the "ECMWF-like" code.  Gases mostly have
C     ! different pressure and temperature scaling of their pathlengths
C     ! in each band, so that there are 6 absorber amounts for water
C     ! vapour line absorption, 1 for each of the foreign-broadened &
C     ! self-broadened water vapour continua, 2 for CO2 (for 3 different
C     ! bands) and 2 for ozone (with and without pressure broadening,
C     ! though only in one band).
C*L
      INTEGER!, INTENT (IN) ::
     &     L2,                       ! Number of points to be treated   
     &     NWET,                     ! Number of levels with moisture -
C                                    ! above these zero is used.
     &     NOZONE,                   ! Number of levels with ozone data
C     ! provided - below these the value in the lowest of them is used
     &     NLEVS,                    ! Number of levels
     &     L1                        ! First dimension of input arrays
C     !  (The different physical assumptions about water vapour and
C     !  ozone in levels where no data is provided means that separate
C     !  loops are used for levels with and without water vapour but
C     !  only the indexing needs changing for levels with and without
C     !  their own ozone data.)
      REAL!, INTENT(IN) ::
     &     H2O(L1,NWET), CO2,        ! Mass mixing ratio (mK in UMDP 23)
     &     O3(L1,NOZONE),            !             of each absorbing gas
     &     TAC(L1,NLEVS),            ! Mid-layer temperatures
     &     PSTAR(L1),                ! Surface pressure
     &     AC(NLEVS), BC(NLEVS),     ! A & B for layer centres and
     &     AB(NLEVS+1), BB(NLEVS+1)  !                       boundaries
      REAL!, INTENT(OUT) ::
     &     DPATH(L2,NGASES,NLEVS)
C     !  The scaled pathlengths are returned in DPATH, indexed by NGASES
C     ! 1-6 are H2O line absorption in bands 1-6 respectively, 7 is CO2
C     ! scaled for band 2, 8 is CO2 scaled for bands 3 & 4, 9 & 10 are
C     ! foreign & self-broadened water vapour continuum, and 11 and 12
C     ! are ozone without and with pressure scaling.
CL    !  LWPTSC has no EXTERNAL calls and no significant structure
C     ! *COMDECK LWABTSAA holds numbers used to apply temperature
C     !   scaling to gaseous pathlengths similarly to the ECMWF LW code.
      REAL TRTSAA                    !  Reference temperature for
      PARAMETER ( TRTSAA = 250. )    ! temperature scaling of pathlength
      REAL O3T1, O3T2, O3T3, O3T4    !  Coefficients giving the
      PARAMETER ( O3T1 = -.326E-03,  ! temperature dependence of ozone
     &     O3T2 = -.102E-05,         !                       absorption
     &     O3T3 = .274E-02, O3T4 = -.107E-04 )
      REAL ABTSAA(3,8,2)             !  Polynomials giving the
C     ! temperature dependence of the absorption for water vapour & CO2.
C     !  The second subscript indexes water vapour in each band in std
C     !  order, & then CO2 in band 2 & 3/4.
C     !  (N.B. This means 3-5 permute onto ZTT/PTT(3-5) in the original)
C
      DATA ABTSAA /
     &   0.298199E-02,  -.394023E-03,  0.319566E-04,
     &   0.143676E-01,  0.366501E-02,  -.160822E-02,
     &   0.197861E-01,  0.315541E-02,  -.174547E-02,
     &   0.289560E-01,  -.208807E-02,  -.121943E-02,
     &   0.103800E-01,  0.436296E-02,  -.161431E-02,
     &   0.868859E-02,  -.972752E-03,  0.000000E-00,
     &   0.250073E-03,  0.455875E-03,  0.109242E-03,
     &   0.307423E-01,  0.110879E-02,  -.322172E-03,
C
     &  -0.106432E-04,  0.660324E-06,  0.174356E-06,
     &  -0.553979E-04,  -.101701E-04,  0.920868E-05,
     &  -0.877012E-04,  0.513302E-04,  0.523138E-06,
     &  -0.165960E-03,  0.157704E-03,  -.146427E-04,
     &   -.276744E-04,  -.327381E-04,  0.127646E-04,
     &   -.278412E-04,  -.713940E-06,  0.117469E-05,
     &   0.199846E-05,  -.216313E-05,  0.175991E-06,
     &  -0.108482E-03,  0.258096E-05,  -.814575E-06 /
      REAL RLNR10,                   !  lg(sqrt(e))
     &     EPSP1,                    !  1 + epsilon
     &     DIFFAC                    !  Diffusivity factor
C
      REAL DABBYG,                   !  Difference of As & Bs across
     &     DBBBYG,                   !  model layer, divided by 10 g.
     &     DABMBP,                   !  Mean As & Bs across model layer,
     &     DBBMBP,                   !   divided by a reference pressure
C     ! These four are used to calculate the next two quantities :
     &     DPBYGA,                   !  Pressure difference across model
C                                    !   layer, divided by 10 g.
     &     pbypr(l2),                !  Mean of layer-boundary pressure,
C                                    !   divided by a reference pressure
C     ! These two together give the pressure scaled pathlength - the
C     !   integral across the layer with respect to pressure of the
C     !   local pressure over the reference one.
     &     TN,                       !  Temperature diffce from TRTSAA
     &     dh2o(l2),                 !  Water vapour pathlength for a   
C                                    !   single layer, pressure-scaled
     &     UPH2O,                    !  Logarithmic function of DH2O
C     !                     used to calculate the temperature scaling
     &     dco2(l2),                 !  Equivalents of DH2O & UPH2O for 
     &     UPCO2,                    !                             CO2.
     &     do3(l2),                  ! Unscaled 1-layer ozone pathlength
     &     TSCCON,                   ! Temperature-scaling term for the
C     !                self-broadened water vapour continuum pathlength
     &     SBWV,                     ! Water-vapour fraction of the
C     !   atmosphere for calculating water vapour continuum pathlengths
     &     X,                        !  Dummy argument for statement    
     & exp_d(l2*10)
C                                    !                       functions
      INTEGER LEVEL, J,              ! Loopers over levels & points
     &     ONETWO,                   ! Flipper
     &     OLEVEL                    ! Index for the ozone data to be
C                                    !        used in the current level
C*
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      PARAMETER ( EPSP1 = 1. + EPSILON )
      PARAMETER ( DIFFAC = 1.66 )
C     !  Simplify the code for temperature scaling of water vapour line
C     !    and CO2 pathlengths by defining statement functions:
      REAL POLYTS,
     &     TSCAL
      POLYTS(X,J,ONETWO) = ABTSAA(1,J,ONETWO) +
     &     X * ( ABTSAA(2,J,ONETWO) + X * ABTSAA(3,J,ONETWO) )
      TSCAL(TN,X,J) = TN * ( POLYTS(X,J,1) + TN * POLYTS(X,J,2) )
C     !  FORTRAN 77 will not allow the following constants to be
C     !  defined in a PARAMETER statement, but the CRAY compiler will
C     !  give the same effect as if they were.
      RLNR10 = .5 / LOG (10.)
C
      DO 2 LEVEL=1, NWET
        DABBYG = ( AB(LEVEL) - AB(LEVEL+1) ) / ( G * 10. )
        DBBBYG = ( BB(LEVEL) - BB(LEVEL+1) ) / ( G * 10. )
        DABMBP = ( AB(LEVEL) + AB(LEVEL+1) ) * 0.5 / 101325.
        DBBMBP = ( BB(LEVEL) + BB(LEVEL+1) ) * 0.5 / 101325.
        OLEVEL = MAX (1, LEVEL+NOZONE-NLEVS)


        DO 20 J=1, L2
          DPBYGA = DABBYG + PSTAR(J) * DBBBYG
          pbypr(j) = DABMBP + PSTAR(J) * DBBMBP           
          TN = TAC(J,LEVEL) - TRTSAA
C
          dh2o(j) = DPBYGA * H2O(J,LEVEL) * pbypr(j)
          z=1.0
          if(dh2o(j).eq.0) dh2o(j)=sign(dh2o(j),z)
                             
          IF ( dh2o(j) .NE. 0. ) THEN                                   
             UPH2O =
     &         AMIN1 ( AMAX1( RLNR10 * alog(dh2o(j))+5.,0.), 6.0)   
           ELSE
             UPH2O = 0.
          ENDIF
          dh2o(j) = dh2o(j) * DIFFAC                                    
          exp_d(j+0*l2) = TSCAL (TN, UPH2O, 1)
          exp_d(j+1*l2) = TSCAL (TN, UPH2O, 2)
          exp_d(j+2*l2) = TSCAL (TN, UPH2O, 3)
          exp_d(j+3*l2) = TSCAL (TN, UPH2O, 4)
          exp_d(j+4*l2) = TSCAL (TN, UPH2O, 5)
          exp_d(j+5*l2) = TSCAL (TN, UPH2O, 6)
C
          dco2(j) = DPBYGA * CO2 * pbypr(j)                  
          UPCO2 = AMAX1 ( RLNR10 * alog ( dco2(j) ) + 5., 0.)  
          dco2(j) = dco2(j) * DIFFAC                                    
          exp_d(j+6*l2) = TSCAL (TN, UPCO2, 7)
          exp_d(j+7*l2) = TSCAL (TN, UPCO2, 8)
C
          TSCCON = EXP ( 6.08 * ( 296. / TAC(J,LEVEL) - 1. ) )
          SBWV = EPSP1 * H2O(J,LEVEL) / ( 1. + EPSILON*H2O(J,LEVEL) )
          DPATH(J,9,LEVEL) = (1.-SBWV) * dh2o(j)    
          DPATH(J,10,LEVEL) = SBWV * dh2o(j) * TSCCON                   
C
          do3(j) = DIFFAC * DPBYGA * O3(J,OLEVEL)                       
          exp_d(j+8*l2) = TN * ( O3T1 + TN * O3T2 )
          exp_d(j+9*l2) = TN * ( O3T3 + TN * O3T4 )
   20   CONTINUE
      do j=1,l2*10
        exp_d(j)=exp(exp_d(j))
      end do
      do j=1,l2
          DPATH(J,1,LEVEL) = DH2O(j) * exp_d(j+0*l2)
          DPATH(J,2,LEVEL) = DH2O(j) * exp_d(j+1*l2)
      enddo
      do j=1,l2
          DPATH(J,3,LEVEL) = DH2O(j) * exp_d(j+2*l2)
          DPATH(J,4,LEVEL) = DH2O(j) * exp_d(j+3*l2)
      enddo
      do j=1,l2
          DPATH(J,5,LEVEL) = DH2O(j) * exp_d(j+4*l2)
          DPATH(J,6,LEVEL) = DH2O(j) * exp_d(j+5*l2)
      enddo
      do j=1,l2
          DPATH(J,7,LEVEL) = DCO2(j) * exp_d(j+6*l2)
          DPATH(J,8,LEVEL) = DCO2(j) * exp_d(j+7*l2)
      enddo
      do j=1,l2
          DPATH(J,11,LEVEL) = DO3(j) * exp_d(j+8*l2)
          DPATH(J,12,LEVEL) = PBYPR(j) * DO3(j)* exp_d(j+9*l2)
      enddo
    2 CONTINUE
C
C     !  for the H2O pathlength but treat CO2 and O3 the same:
C
      DO 3 LEVEL=NWET+1, NLEVS
        DABBYG = ( AB(LEVEL) - AB(LEVEL+1) ) / ( G * 10. )
        DBBBYG = ( BB(LEVEL) - BB(LEVEL+1) ) / ( G * 10. )
        DABMBP = ( AB(LEVEL) + AB(LEVEL+1) ) * 0.5 / 101325.
        DBBMBP = ( BB(LEVEL) + BB(LEVEL+1) ) * 0.5 / 101325.
        OLEVEL = MAX (1, LEVEL+NOZONE-NLEVS)
        DO 30 J=1, L2
          DPBYGA = DABBYG + PSTAR(J) * DBBBYG
          PBYPR(J) = DABMBP + PSTAR(J) * DBBMBP
          TN = TAC(J,LEVEL) - TRTSAA
C
          DCO2(J) = DPBYGA * CO2 * PBYPR(J)
          UPCO2 = AMAX1 ( RLNR10 * alog ( DCO2(j) ) + 5., 0.)
          DCO2(J) = DCO2(J) * DIFFAC
          exp_d(j+0*l2) = TSCAL (TN, UPCO2, 7)
          exp_d(j+1*l2) = TSCAL (TN, UPCO2, 8)
C                                                                       
          DO3(J) = DIFFAC * DPBYGA * O3(J,OLEVEL)
          exp_d(j+2*l2) = TN * ( O3T1 + TN * O3T2 )
          exp_d(j+3*l2) = TN * ( O3T3 + TN * O3T4 )
   30   CONTINUE                                                        
      do j=1,l2*4
        exp_d(j)=exp(exp_d(j))
      end do
      do j=1,l2
          DPATH(J,1,LEVEL) = 0.
          DPATH(J,2,LEVEL) = 0.
      enddo
      do j=1,l2
          DPATH(J,3,LEVEL) = 0.
          DPATH(J,4,LEVEL) = 0.
      enddo
      do j=1,l2
          DPATH(J,5,LEVEL) = 0.
          DPATH(J,6,LEVEL) = 0.
      enddo
      do j=1,l2
          DPATH(J,7,LEVEL) = DCO2(j) * exp_d(j+0*l2)
          DPATH(J,8,LEVEL) = DCO2(j) * exp_d(j+1*l2)
      enddo
      do j=1,l2
          DPATH(J,9,LEVEL) = 0.
          DPATH(J,10,LEVEL) = 0.
      enddo
      do j=1,l2
          DPATH(J,11,LEVEL) = DO3(j) * exp_d(j+2*l2)
          DPATH(J,12,LEVEL) = PBYPR(j) * DO3(j) * exp_d(j+3*l2)
      enddo
    3 CONTINUE
C
      RETURN
      END
