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
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_GAS_MIX_RATIO(IERR
     &   , N_PROFILE, NLEVS, NWET, NOZONE
     &   , I_GATHER
     &   , N_ABSORB, TYPE_ABSORB
     &   , L_N2O, L_CH4, L_CFC11, L_CFC12, L_O2
     &   , L_CFC113, L_HCFC22, L_HFC125, L_HFC134A
     &   , H2O, CO2, O3, N2O_MIX_RATIO, CH4_MIX_RATIO
     &   , C11_MIX_RATIO, C12_MIX_RATIO, O2_MIX_RATIO
     &   , C113_MIX_RATIO, HCFC22_MIX_RATIO
     &   , HFC125_MIX_RATIO, HFC134A_MIX_RATIO
     &   , GAS_MIX_RATIO
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_SPECIES
     &   )
!
!
!     COMDECKS INCLUDED
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INDEXING NUMBERS OF GASEOUS ABSORBING SPECIES.
!     THE NUMBERING 1-12 CORRESPONDS TO LOWTRAN 7.
!
      INTEGER
     &     NPD_GASES
!             NUMBER OF INDEXED GASES
      PARAMETER(NPD_GASES=19)
!
      INTEGER
     &     IP_H2O
!             INDEX NUMBER OF WATER VAPOUR
     &   , IP_CO2
!             INDEX NUMBER OF CARBON DIOXIDE
     &   , IP_O3
!             INDEX NUMBER OF OZONE
     &   , IP_N2O
!             INDEX NUMBER OF DINITROGEN OXIDE
     &   , IP_CO
!             INDEX NUMBER OF CARBON MONOXIDE
     &   , IP_CH4
!             INDEX NUMBER OF METHANE
     &   , IP_O2
!             INDEX NUMBER OF OXYGEN
     &   , IP_NO
!             INDEX NUMBER OF NITROGEN MONOXIDE
     &   , IP_SO2
!             INDEX NUMBER OF SULPHUR DIOXIDE
     &   , IP_NO2
!             INDEX NUMBER OF NITROGEN DIOXIDE
     &   , IP_NH3
!             INDEX NUMBER OF AMMONIA
     &   , IP_HNO3
!             INDEX NUMBER OF NITRIC ACID
     &   , IP_N2
!             INDEX NUMBER OF NITROGEN
     &   , IP_CFC11
!             INDEX NUMBER OF CFC11
     &   , IP_CFC12
!             INDEX NUMBER OF CFC12
     &   , IP_CFC113
!             INDEX NUMBER OF CFC113
     &   , IP_HCFC22
!             INDEX NUMBER OF HCFC22
     &   , IP_HFC125
!             INDEX NUMBER OF HFC125
     &   , IP_HFC134A
!             INDEX NUMBER OF HCF134A
!
      PARAMETER(
     &     IP_H2O=1
     &   , IP_CO2=2
     &   , IP_O3=3
     &   , IP_N2O=4
     &   , IP_CO=5
     &   , IP_CH4=6
     &   , IP_O2=7
     &   , IP_NO=8
     &   , IP_SO2=9
     &   , IP_NO2=10
     &   , IP_NH3=11
     &   , IP_HNO3=12
     &   , IP_N2=13
     &   , IP_CFC11=14
     &   , IP_CFC12=15
     &   , IP_CFC113=16
     &   , IP_HCFC22=17
     &   , IP_HFC125=18
     &   , IP_HFC134A=19
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF ARRAY FROM UM
     &   , NPD_PROFILE
!             SIZE OF ARRAY
     &   , NPD_LAYER
!             SIZE OF ARRAY
     &   , NPD_SPECIES
!             SIZE OF ARRAY
!
!     SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF LEVELS
     &   , NWET
!             NUMBER OF WET LEVELS
     &   , NOZONE
!             NUMBER OF OZONE LEVELS
!
!     GATHERING ARRAY:
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
!
!     TYPES OF GASES:
      INTEGER   !, INTENT(IN)
     &     N_ABSORB
!             NUMBER OF ABSORBERS
     &   , TYPE_ABSORB(NPD_SPECIES)
!             TYPES OF ABSORBERS
!
!     FLAGS FOR MINOR GASES:
      LOGICAL   !,INTENT(IN)
     &     L_N2O
!             FLAG FOR NITROUS OXIDE
     &   , L_CH4
!             FLAG FOR METHANE
     &   , L_CFC11
!             FLAG FOR CFC11
     &   , L_CFC12
!             FLAG FOR CFC12
     &   , L_O2
!             FLAG FOR O2
     &   , L_CFC113
!             FLAG FOR CFC113
     &   , L_HCFC22
!             FLAG FOR HCFC22
     &   , L_HFC125
!             FLAG FOR HFC125
     &   , L_HFC134A
!             FLAG FOR HFC134A
!
!     MIXING RATIOS SUPPLIED:
      INTEGER  CO2_DIM1, CO2_DIM2   ! dimensions of CO2_3D field
      LOGICAL  L_CO2_3D    !  controls use of 3D co2 field
      REAL      !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)
!             MASS MIXING RATIO OF WATER VAPOUR
     &   , CO2
!             MASS MIXING RATIO OF CARBON DIOXIDE
     &   , CO2_3D(CO2_DIM1, CO2_DIM2)
!             3D MASS MIXING RATIO OF CO2 (full field)
     &   , O3(NPD_FIELD, NOZONE)
!             MASS MIXING RATIO OF OZONE
     &   , N2O_MIX_RATIO
!             MASS MIXING RATIO OF NITROUS OXIDE
     &   , CH4_MIX_RATIO
!             MASS MIXING RATIO OF METHANE
     &   , C11_MIX_RATIO
!             MASS MIXING RATIO OF CFC11
     &   , C12_MIX_RATIO
!             MASS MIXING RATIO OF CFC12
     &   , O2_MIX_RATIO
!             MASS MIXING RATIO OF O2
     &   , C113_MIX_RATIO
!             MASS MIXING RATIO OF CFC113
     &   , HCFC22_MIX_RATIO
!             MASS MIXING RATIO OF HCFC22
     &   , HFC125_MIX_RATIO
!             MASS MIXING RATIO OF HFC125
     &   , HFC134A_MIX_RATIO
!             MASS MIXING RATIO OF HFC134A
!
!     ARRAY OF ASSIGNED MXING RATIOS:
      REAL      !, INTENT(OUT)
     &     GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             MIXING RATIOS
!
!     LOCAL VARIABLES.
!
!     POINTERS TO GASES:
      INTEGER
     &     IUMP_H2O
!             POINTER TO WATER VAPOUR
     &   , IUMP_CO2
!             POINTER TO CARBON DIOXIDE
     &   , IUMP_O3
!             POINTER TO OZONE
     &   , IUMP_N2O
!             POINTER TO NITOUS OXIDE
     &   , IUMP_CH4
!             POINTER TO METHANE
     &   , IUMP_CFC11
!             POINTER TO CFC11
     &   , IUMP_CFC12
!             POINTER TO CFC12
     &   , IUMP_O2
!             POINTER TO O2
     &   , IUMP_CFC113
!             POINTER TO CFC113
     &   , IUMP_HCFC22
!             POINTER TO HCFC22
     &   , IUMP_HFC125
!             POINTER TO HFC125
     &   , IUMP_HFC134A
!             POINTER TO HFC134A
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , LG
!             CORRESPONDING UNGATHERED INDEX
!
      REAL
     &     H2OLMN
!             LIMITING CONCENTRATION OF WATER VAPOUR
      PARAMETER(H2OLMN=1.E-8)
!
!
!
!
!     MATCH THE INDEXING NUMBERS OF GASEOUS SPECIES IN THE SPECTRAL
!     FILE WITH ACTUAL TYPES OF GASES KNOWN TO THE UM.
!
!     SET ALL POINTERS TO 0 INITIALLY TO FLAG MISSING GASES.
      IUMP_H2O=0
      IUMP_CO2=0
      IUMP_O3=0
      IUMP_N2O=0
      IUMP_CH4=0
      IUMP_CFC11=0
      IUMP_CFC12=0
      IUMP_O2=0
      IUMP_CFC113=0
      IUMP_HCFC22=0
      IUMP_HFC125=0
      IUMP_HFC134A=0
!
!
      DO I=1, N_ABSORB
!
         IF (TYPE_ABSORB(I).EQ.IP_H2O) THEN
            IUMP_H2O=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_CO2) THEN
            IUMP_CO2=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_O3) THEN
            IUMP_O3=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_N2O) THEN
            IUMP_N2O=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_CH4) THEN
            IUMP_CH4=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_CFC11) THEN
            IUMP_CFC11=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_CFC12) THEN
            IUMP_CFC12=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_O2) THEN
            IUMP_O2=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_CFC113) THEN
            IUMP_CFC113=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_HCFC22) THEN
            IUMP_HCFC22=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_HFC125) THEN
            IUMP_HFC125=I
         ELSE IF (TYPE_ABSORB(I).EQ.IP_HFC134A) THEN
            IUMP_HFC134A=I
         ENDIF
!
      ENDDO
!
!
!     ASSIGN MIXING RATIOS OF THE GASES TO THE MAIN ARRAYS.
!
!     WATER VAPOUR:
!
      IF (IUMP_H2O.GT.0) THEN
!        THE UPPER LEVELS RECEIVE A CONSTANT SMALL VALUE.
         DO I=1, NLEVS-NWET
            DO L=1, N_PROFILE
               GAS_MIX_RATIO(L, I, IUMP_H2O)=H2OLMN
            ENDDO
         ENDDO
         DO I=NLEVS-NWET+1, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_H2O)=H2O(LG, NLEVS-I+1)
            ENDDO
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: WATER VAPOUR IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     CARBON DIOXIDE:
!
      IF (IUMP_CO2.GT.0) THEN
         DO I=1, NLEVS
           IF (L_CO2_3D) THEN
             DO L=1, N_PROFILE   
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_CO2)=CO2_3D(LG, NLEVS-I+1)
             ENDDO
           ELSE
             DO L=1, N_PROFILE   
               GAS_MIX_RATIO(L, I, IUMP_CO2)=CO2
             ENDDO
           ENDIF
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: CARBON DIOXIDE IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     OZONE:
!
      IF (IUMP_O3.GT.0) THEN
!        THE CLIMATOLOGY OF OZONE IS GIVEN ON NOZONE LEVELS,
!        THE LOWEST VALUE SUPPLYING THE MIXING RATIO ON
!        ALL LOWER LEVELS.
         DO I=1, NOZONE
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_O3)=O3(LG, NOZONE+1-I)
            ENDDO
         ENDDO
         DO I=NOZONE+1, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_O3)=O3(LG, 1)
            ENDDO
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: OZONE IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     OTHER TRACE GASES:
!
!     THESE GASES ARE NOT ALWAYS INCLUDED IN THE CALCULATION.
!     TESTING IS THEREFORE MORE INTRICATE.
!
      IF (IUMP_N2O.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_N2O) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_N2O)=N2O_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_N2O)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_N2O) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: NITROUS OXIDE IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CH4.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CH4) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CH4)=CH4_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CH4)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CH4) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: METHANE IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC11.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC11) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC11)=C11_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC11)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC11) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: CFC11 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC12.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC12) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC12)=C12_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC12)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC12) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: CFC12 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_O2.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_O2) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_O2)=O2_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_O2)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_O2) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: O2 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC113.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC113) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC113)=C113_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC113)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC113) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: CFC113 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HCFC22.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HCFC22) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HCFC22)=HCFC22_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HCFC22)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HCFC22) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: HCFC22 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HFC125.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HFC125) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC125)=HFC125_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC125)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HFC125) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: HFC125 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HFC134A.GT.0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HFC134A) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC134A)=HFC134A_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC134A)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HFC134A) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: HFC134A IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
!
      RETURN
      END
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_THERMODYNAMIC(
     &     N_PROFILE, NLEVS, I_GATHER, L_BOUNDARY_TEMPERATURE
     &   , PSTAR, TSTAR, AB, BB, AC, BC, PEXNER, TAC
     &   , P, T, T_BDY, T_SURFACE, D_MASS
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDED COMDECKS
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

!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF ARRAY FROM UM
     &   , NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF LEVELS
     &   , I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
      REAL      !, INTENT(IN)
     &     PSTAR(NPD_FIELD)
!             SURFACE PRESSURES
     &   , TSTAR(NPD_FIELD)
!             SURFACE TEMPERSTURES
     &   , AB(NLEVS+1)
!             A AT EDGES OF LAYERS
     &   , BB(NLEVS+1)
!             B AT EDGES OF LAYERS
     &   , AC(NLEVS)
!             A AT CENTRES OF LAYERS
     &   , BC(NLEVS)
!             B AT CENTRES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
     &   , PEXNER(NPD_FIELD, NLEVS+1)
!             EXNER FUNCTION AT BOUNDARIES
!
      LOGICAL   !, INTENT(IN)
     &     L_BOUNDARY_TEMPERATURE
!             FLAG TO CALCULATE TEMPERATURES AT BOUNADRIES OF LAYERS.
!
!
      REAL      !, INTENT(OUT)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE FIELD
     &   , T_BDY(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURES AT EDGES OF LAYERS
     &   , T_SURFACE(NPD_PROFILE)
!             GATHERED TEMPERATURE OF SURFACE
!
!
!     LOCAL VARIABLES.
!
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , II
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , LG
!             INDEX TO GATHER
!
      REAL
     &     PU
!             PRESSURE FOR UPPER LAYER
     &   , PL
!             PRESSURE FOR LOWER LAYER
     &   , PML1
!             PRESSURE FOR INTERPOLATION
     &   , WTL
!             WEIGHT FOR LOWER LAYER
     &   , WTU
!             WEIGHT FOR UPPER LAYER
!
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------

!
!
!
!     CALCULATE PROPERTIES AT THE CENTRES OF LAYERS.
      DO I=1, NLEVS
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
            P(L, I)=AC(NLEVS+1-I)+BC(NLEVS+1-I)*PSTAR(LG)
            T(L, I)=TAC(LG, NLEVS+1-I)
            D_MASS(L, I)=(AB(NLEVS+1-I)-AB(NLEVS+2-I)
     &         +PSTAR(LG)*(BB(NLEVS+1-I)-BB(NLEVS+2-I)))
     &         /G
         ENDDO
      ENDDO
!
!
      IF (L_BOUNDARY_TEMPERATURE) THEN
!
!        GATHER THE SURFACE TEMPERATURE.
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
            T_SURFACE(L)=TSTAR(LG)
         ENDDO
!
!        INTERPOLATE TEMPERATURES AT THE BOUNDARIES OF LAYERS
!        FROM THE EXNER FUNCTION.
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
!
!           TAKE THE TEMPERATURE OF THE AIR JUST ABOVE THE SURFACE AS
!           THE TEMPERATURE AT THE MIDDLE OF THE BOTTOM LAYER.
            T_BDY(L, NLEVS)=TAC(LG, 1)
!           TAKE THE TEMPERATURE AS CONSTANT ACROSS THE TOP HALF-LAYER.
            T_BDY(L, 0)=TAC(LG, NLEVS)
!
         ENDDO
!
         DO I=1, NLEVS-1
            II=NLEVS-I
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               PU=PSTAR(LG)*BB(II+2)+AB(II+2)
               PL=PSTAR(LG)*BB(II+1)+AB(II+1)
               PML1=PSTAR(LG)*BB(II)+AB(II)
               WTU=TAC(LG, II+1)*(PEXNER(LG, II+1)
     &            /P_EXNER_C(PEXNER(LG, II+2), PEXNER(LG, II+1)
     &            , PU, PL, KAPPA)-1.0E+00)
               WTL=TAC(LG, II)*(PEXNER(LG, II)
     &            /P_EXNER_C(PEXNER(LG, II+1), PEXNER(LG, II)
     &            , PL, PML1, KAPPA)-1.0E+00)
               T_BDY(L, I)=(WTU*TAC(LG, NLEVS+1-I)
     &            +WTL*TAC(LG, NLEVS-I))/(WTL+WTU)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_FIELD(N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , P, T, D_MASS
     &   , CCB, CCT, CCA, CCCWP
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK
     &   , L_MICROPHYSICS, L_AEROSOL_CCN
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE
     &   , L_CLOUD_WATER_PARTITION, LAND_G
     &   , I_CLOUD_REPRESENTATION, I_CONDENSED_PARAM
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , W_CLOUD, FRAC_CLOUD, L_LOCAL_CNV_PARTITION
     &   , CONDENSED_MIX_RAT_AREA, CONDENSED_DIM_CHAR
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG
     &   , LWP_STRAT, LWP_STRAT_FLAG
     &   , NTOT_DIAG, NTOT_DIAG_FLAG
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &   , N_CCA_LEV, L_3D_CCA
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INTERNAL DIMENSIONS TIED TO ALGORITHMS,
!     MOSTLY FOR CLOUDS.
!
      INTEGER
     &     NPD_CLOUD_COMPONENT
!             NUMBER OF COMPONENTS OF CLOUDS
     &   , NPD_CLOUD_TYPE
!             NUMBER OF PERMITTED TYPES OF CLOUDS.
     &   , NPD_CLOUD_REPRESENTATION
!             NUMBER OF PERMITTED REPRESENTATIONS OF CLOUDS.
     &   , NPD_OVERLAP_COEFF
!             NUMBER OF OVERLAP COEFFICIENTS FOR CLOUDS
     &   , NPD_SOURCE_COEFF
!             NUMBER OF COEFFICIENTS FOR TWO-STREAM SOURCES
     &   , NPD_REGION
!             NUMBER OF REGIONS IN A LAYER
!
      PARAMETER(
     &     NPD_CLOUD_COMPONENT=4
     &   , NPD_CLOUD_TYPE=4
     &   , NPD_CLOUD_REPRESENTATION=4
     &   , NPD_OVERLAP_COEFF=18
     &   , NPD_SOURCE_COEFF=2
     &   , NPD_REGION=3
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET COMPONENTS OF CLOUDS.
!
      INTEGER
     &     IP_CLCMP_ST_WATER
!             STRATIFORM WATER DROPLETS
     &   , IP_CLCMP_ST_ICE
!             STRATIFORM ICE CRYSTALS
     &   , IP_CLCMP_CNV_WATER
!             CONVECTIVE WATER DROPLETS
     &   , IP_CLCMP_CNV_ICE
!             CONVECTIVE ICE CRYSTALS
!
      PARAMETER(
     &     IP_CLCMP_ST_WATER=1
     &   , IP_CLCMP_ST_ICE=2
     &   , IP_CLCMP_CNV_WATER=3
     &   , IP_CLCMP_CNV_ICE=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET TYPES OF CLOUDS.
!
      INTEGER
     &     IP_CLOUD_TYPE_HOMOGEN
!             CLOUD COMPOSED OF MIXED WATER AND ICE
     &   , IP_CLOUD_TYPE_WATER
!             CLOUD COMPOSED ONLY OF WATER
     &   , IP_CLOUD_TYPE_ICE
!             CLOUD COMPOSED ONLY OF ICE
     &   , IP_CLOUD_TYPE_STRAT
!             MIXED-PHASE STRATIFORM CLOUD
     &   , IP_CLOUD_TYPE_CONV
!             MIXED-PHASE CONVECTIVE CLOUD
     &   , IP_CLOUD_TYPE_SW
!             STRATIFORM WATER CLOUD
     &   , IP_CLOUD_TYPE_SI
!             STRATIFORM ICE CLOUD
     &   , IP_CLOUD_TYPE_CW
!             CONVECTIVE WATER CLOUD
     &   , IP_CLOUD_TYPE_CI
!             CONVECTIVE ICE CLOUD
!
      PARAMETER(
     &     IP_CLOUD_TYPE_HOMOGEN=1
     &   , IP_CLOUD_TYPE_WATER=1
     &   , IP_CLOUD_TYPE_ICE=2
     &   , IP_CLOUD_TYPE_STRAT=1
     &   , IP_CLOUD_TYPE_CONV=2
     &   , IP_CLOUD_TYPE_SW=1
     &   , IP_CLOUD_TYPE_SI=2
     &   , IP_CLOUD_TYPE_CW=3
     &   , IP_CLOUD_TYPE_CI=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET REPRESENTATIONS OF CLOUDS.
!
!     REPRESENTATIONS
      INTEGER
     &     IP_CLOUD_HOMOGEN
!             ALL COMPONENTS ARE MIXED HOMOGENEOUSLY
     &   , IP_CLOUD_ICE_WATER
!             ICE AND WATER CLOUDS ARE TREATED SEPARATELY
     &   , IP_CLOUD_CONV_STRAT
!             CLOUDS ARE DIVIDED INTO HOMOGENEOUSLY MIXED
!             STRATIFORM AND CONVECTIVE PARTS
     &   , IP_CLOUD_CSIW
!             CLOUDS DIVIDED INTO ICE AND WATER PHASES AND
!             INTO STRATIFORM AND CONVECTIVE COMPONENTS.
!
      PARAMETER(
     &     IP_CLOUD_HOMOGEN=1
     &   , IP_CLOUD_ICE_WATER=2
     &   , IP_CLOUD_CONV_STRAT=3
     &   , IP_CLOUD_CSIW=4
     &   )
!
!     TYPES OF CLOUDS:
      INTEGER
     &     NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)
!             NUMBER OF TYPE OF CLOUDS IN REPRESENTATION
!
      INTEGER
     &     IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT
     &       , NPD_CLOUD_REPRESENTATION)
!            MAP OF COMPONENTS CONTRIBUTING TO TYPES OF CLOUDS
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET NUMBERS FOR ICE CLOUD SCHEMES.
!
      INTEGER
     &     NPD_ICE_CLOUD_FIT
!             NUMBER OF CLOUD FITTING SCHEMES
     &   , IP_SLINGO_SCHRECKER_ICE
!             PARAMETRIZATION OF SLINGO AND SCHRECKER.
     &   , IP_ICE_UNPARAMETRIZED
!             UNPARAMETRIZED ICE CRYSTAL DATA
     &   , IP_SUN_SHINE_VN2_VIS
!             SUN AND SHINE'S PARAMETRIZATION IN THE VISIBLE (VERSION 2)
     &   , IP_SUN_SHINE_VN2_IR
!             SUN AND SHINE'S PARAMETRIZATION IN THE IR (VERSION 2)
     &   , IP_ICE_ADT
!             SCHEME BASED ON ANOMALOUS DIFFRACTION THEORY
!             FOR ICE CRYSTALS

!
      PARAMETER(
     &     NPD_ICE_CLOUD_FIT=6
     &   , IP_SLINGO_SCHRECKER_ICE=1
     &   , IP_ICE_UNPARAMETRIZED=3
     &   , IP_SUN_SHINE_VN2_VIS=4
     &   , IP_SUN_SHINE_VN2_IR=5
     &   , IP_ICE_ADT=6
     &   )
!
!     ------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
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

!
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL_SPECIES
     &   , SULP_DIM1
!             1ST DIMENSION OF ARRAYS OF SULPHATE
     &   , SULP_DIM2
!             2ND DIMENSION OF ARRAYS OF SULPHATE
     &   , N_CCA_LEV
!             NUMBER OF LEVELS FOR CONVECTIVE CLOUD AMOUNT
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
!
!     GATHERING ARRAY:
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
!
!     THERMODYNAMIC FIELDS:
      REAL      !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURES
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURES
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESSES OF LAYERS
!
!     CONVECTIVE CLOUDS:
      INTEGER   !, INTENT(IN)
     &     CCB(NPD_FIELD)
!             BASE OF CONVECTIVE CLOUD
     &   , CCT(NPD_FIELD)
!             TOP OF CONVECTIVE CLOUD
      REAL      !, INTENT(IN)
     &     CCA(NPD_FIELD,N_CCA_LEV)
!             FRACTION OF CONVECTIVE CLOUD
     &   , CCCWP(NPD_FIELD)
!             WATER PATH OF CONVECTIVE CLOUD
      LOGICAL   !, INTENT(IN)
     &     L_3D_CCA
     &   , L_LOCAL_CNV_PARTITION
!             FLAG TO CARRY OUT THE PARTITIONING BETWEEN ICE
!             AND WATER IN CONVECTIVE CLOUDS AS A FUNCTION OF
!             THE LOCAL TEMPERATURE
!
!     LAYER CLOUDS:
      REAL      !, INTENT(IN)
     &     LCCWC1(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             LIQUID WATER CONTENTS
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             ICE WATER CONTENTS
     &   , LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             AREA COVERAGE FRACTIONS OF LAYER CLOUDS
     &   , LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK COVERAGE FRACTIONS OF LAYER CLOUDS
!
!     ARRAYS FOR MICROPHYSICS:
      LOGICAL   !, INTENT(IN)
     &     L_MICROPHYSICS
!             MICROPHYSICAL FLAG
     &   , L_AEROSOL_CCN
!             FLAG TO USE AEROSOLS TO FIND CCN
     &   , L_CLOUD_WATER_PARTITION
!             FLAG TO USE PROGNOSTIC CLOUD ICE CONTENTS
     &   , LAND_G(NPD_PROFILE)
!             FLAG FOR LAND POINTS
      REAL      !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF ACCUMULATION-MODE SULPHATE
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF DISSOLVED SULPHATE
!
!     REPRESENTATION OF CLOUDS
      INTEGER   !, INTENT(IN)
     &     I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
!
!     PARAMETRIZATIONS FOR CLOUDS:
      INTEGER   !, INTENT(IN)
     &     I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
!     LIMITS ON SIZES OF PARTICLES
      REAL      !, INTENT(IN)
     &     CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
!     ASSIGNED CLOUD FIELDS:
      INTEGER   !, INTENT(OUT)
     &     N_CONDENSED
!             NUMBER OF CONDENSED COMPONENTS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF CONDENSED COMPONENTS
      REAL      !, INTENT(OUT)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             TOTAL AMOUNTS OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTION OF EACH TYPE OF CLOUD
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER
     &      , NPD_CLOUD_COMPONENT)
!             CHARACTERISTIC DIMENSIONS OF CLOUDY COMPONENTS
     &   , CONDENSED_MIX_RAT_AREA(NPD_PROFILE, 0: NPD_LAYER
     &      , NPD_CLOUD_COMPONENT)
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING AREA CLD
     &   , NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
     &   , STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
     &   , SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
!
!
!     MICROPHYSICAL DIAGNOSTICS:
      LOGICAL
     &     RE_CONV_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV_FLAG
!             DIAGNOSE WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT_FLAG
!             DIAGNOSE WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT_FLAG
!             DIAGNOSE LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG_FLAG
!             DIAGNOSE DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG_FLAG
!             DIAGNOSE STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG_FLAG
!             DIAGNOSE SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT_FLAG
!             DIAGNOSE CONDITIONAL SAMPLING WEIGHT
!
      REAL
     &     RE_CONV(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV(NPD_FIELD, NCLDS)
!             WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT(NPD_FIELD, NCLDS)
!             WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT(NPD_FIELD, NCLDS)
!             LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG(NPD_FIELD, NCLDS)
!             DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG(NPD_FIELD, NCLDS)
!             STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG(NPD_FIELD, NCLDS)
!             SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT(NPD_FIELD, NCLDS)
!             CONDITIONAL SAMPLING WEIGHT
!
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , LG
!             INDEX TO GATHER
      LOGICAL
     &     L_GLACIATED_TOP(NPD_PROFILE)
!             LOGICAL FOR GLACIATED TOPS IN CONVECTIVE CLOUD.

!
      REAL
     &     LIQ_FRAC(NPD_PROFILE)
!             FRACTION OF LIQUID CLOUD WATER
     &   , LIQ_FRAC_CONV(NPD_PROFILE)
!             FRACTION OF LIQUID WATER IN CONVECTIVE CLOUD
     &   , T_GATHER(NPD_PROFILE)
!             GATHERED TEMPERATURE FOR LSP_FOCWWIL
     &   , T_LIST(NPD_PROFILE)
!             LIST OF TEMPERATURES
     &   , TOTAL_MASS(NPD_PROFILE)
!             TOTAL MASS IN CONVECTIVE CLOUD
     &   , CC_DEPTH(NPD_PROFILE)
!             DEPTH OF CONVECTIVE CLOUD
     &   , CONDENSED_MIX_RAT_BULK(NPD_PROFILE, 0: NPD_LAYER
     &      , NPD_CLOUD_COMPONENT)
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING BULK CLD
     &   , DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             DENSITY OF AIR
     &   , CONVECTIVE_CLOUD_LAYER(NPD_PROFILE)
!             AMOUNT OF CONVECTIVE CLOUD IN TH CURRENT LAYER
     &   , MEAN_WATER_CONTENT
!             MEAN WATER CONTENT
     &   , MEAN_ICE_CONTENT
!             MEAN ICE CONTENT
     &   , CONDENSED_LIMIT
!             LOWER LIMIT ON WATER CONTENTS
!
      PARAMETER(CONDENSED_LIMIT=1.E-8)
!
!
!
!     CHECK THE LIMITS FOR CONVECTIVE CLOUD.
      DO L=1, N_PROFILE
         LG=I_GATHER(L)
         IF ( (CCB(LG).GT.NCLDS).OR.(CCB(LG).LT.1) ) CCB(LG)=1
         IF ( (CCT(LG).GT.NCLDS+1).OR.(CCT(LG).LT.2) ) CCT(LG)=NCLDS+1
         IF (L_3D_CCA) THEN
           IF (CCA(LG,CCB(LG)).LT.TOL_TEST) CCCWP(LG)=0.0E+00
         ELSE
           IF (CCA(LG,1).LT.TOL_TEST) CCCWP(LG)=0.0E+00
         ENDIF
      ENDDO
!
!
!     SET THE COMPONENTS WITHIN THE CLOUDS. IN THE UNIFIED MODEL WE
!     HAVE FOUR COMPONENTS: STRATIFORM ICE AND WATER AND CONVECTIVE
!     ICE AND WATER.
      N_CONDENSED=4
      TYPE_CONDENSED(1)=IP_CLCMP_ST_WATER
      TYPE_CONDENSED(2)=IP_CLCMP_ST_ICE
      TYPE_CONDENSED(3)=IP_CLCMP_CNV_WATER
      TYPE_CONDENSED(4)=IP_CLCMP_CNV_ICE
!
!
!
!     SET THE TOTAL AMOUNTS OF CLOUD AND THE FRACTIONS COMPRISED BY
!     CONVECTIVE AND STRATIFORM COMPONENTS.
!
!     ZERO THE AMOUNTS OF CLOUD IN THE UPPER LAYERS.
      DO I=1, NLEVS-NCLDS
         DO L=1, N_PROFILE
            W_CLOUD(L, I)=0.0E+00
         ENDDO
      ENDDO
!
      IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT .AND.
     &    .NOT. L_CLOUD_WATER_PARTITION) THEN
!  This cloud representation not available with new cloud microphysics
!
!        THE CLOUDS ARE DIVIDED INTO MIXED-PHASE STRATIFORM AND
!        CONVECTIVE CLOUDS: LSP_FOCWWIL GIVES THE PARTITIONING BETWEEN
!        ICE AND WATER IN STRATIFORM CLOUDS AND IN CONVECTIVE CLOUD,
!        UNLESS THE OPTION TO PARTITION AS A FUNCTION OF THE LOCAL
!        TEMPERATURE IS SELECTED. WITHIN CONVECTIVE CLOUD THE LIQUID
!        WATER CONTENT IS DISTRIBUTED UNIFORMLY THROUGHOUT THE CLOUD.
!
!        CONVECTIVE CLOUD:
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
            ENDDO
         ENDDO
!
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT 
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!
            DO L=1, N_PROFILE
               IF (T(L, NLEVS+2-CCT(I_GATHER(L))).LT.TM) THEN
                  L_GLACIATED_TOP(L)=.TRUE.
               ELSE
                  L_GLACIATED_TOP(L)=.FALSE.
               ENDIF
            ENDDO

         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!
            DO L=1, N_PROFILE
               T_GATHER(L)=T(L, NLEVS+2-CCT(I_GATHER(L)))
            ENDDO
            CALL LSP_FOCWWIL(T_GATHER, N_PROFILE, LIQ_FRAC_CONV)
!
         ENDIF
!
!
         DO L=1, N_PROFILE
            TOTAL_MASS(L)=0.0E+00
         ENDDO
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG).GE.NLEVS+2-I).AND.
     &              (CCB(LG).LE.NLEVS+1-I) ) THEN
                  TOTAL_MASS(L)=TOTAL_MASS(L)+D_MASS(L, I)
               ENDIF
            ENDDO
         ENDDO
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG).GE.NLEVS+2-I).AND.
     &              (CCB(LG).LE.NLEVS+1-I) ) THEN
                  IF (L_LOCAL_CNV_PARTITION) THEN
!                    THE PARTITIONING IS RECALCULATED FOR EACH LAYER
!                    OTHERWISE A GENERIC VALUE IS USED.
                     LIQ_FRAC_CONV(L)=MAX(0.0E+00, MIN(1.0E+00
     &                  , 1.61E-02*(T(L, I)-TM)+8.9E-01))
                     IF ((T(L, I).GT.TM).AND.(.NOT.L_GLACIATED_TOP(L)))
     &                  LIQ_FRAC_CONV(L)=1.0E+00
                  ENDIF
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
     &               =CCCWP(LG)*LIQ_FRAC_CONV(L)
     &               /(TOTAL_MASS(L)+TOL_MACHINE)
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)
     &               =CCCWP(LG)*(1.0E+00-LIQ_FRAC_CONV(L))
     &               /(TOTAL_MASS(L)+TOL_MACHINE)
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)
     &               =CCCWP(LG)*LIQ_FRAC_CONV(L)
     &               /(TOTAL_MASS(L)+TOL_MACHINE)
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)
     &               =CCCWP(LG)*(1.0-LIQ_FRAC_CONV(L))
     &               /(TOTAL_MASS(L)+TOL_MACHINE)
               ENDIF
            ENDDO
         ENDDO
!
!
!        STRATIFORM CLOUDS:
!
!        PARTITION BETWEEN ICE AND WATER DEPENDING ON THE
!        LOCAL TEMPERATURE.
!
         DO I=1, NCLDS
            CALL LSP_FOCWWIL(T(L, NLEVS+1-I), N_PROFILE, LIQ_FRAC)
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I).GT.TOL_TEST) THEN
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))
     &               *LIQ_FRAC(L)/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))
     &               *(1.0E+00-LIQ_FRAC(L))/LCA_AREA(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I).GT.TOL_TEST) THEN
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))
     &               *LIQ_FRAC(L)/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))
     &               *(1.0E+00-LIQ_FRAC(L))/LCA_BULK(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!        CLOUD FRACTIONS:
!
       IF (L_3D_CCA) THEN
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               W_CLOUD(L, NLEVS+1-I)
     &            =CCA(LG,I)+(1.0E+00-CCA(LG,I))*LCA_AREA(LG, I)
               FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)
     &            =CCA(LG,I)/(W_CLOUD(L, NLEVS+1-I)+TOL_MACHINE)
               FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_STRAT)
     &            =1.0E+00-FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ELSE
         DO I=1, NCLDS
            DO L=1, N_PROFILE
              LG=I_GATHER(L)
               IF ( (I.LE.CCT(LG)-1).AND.(I.GE.CCB(LG)) ) THEN
                  W_CLOUD(L, NLEVS+1-I)
     &               =CCA(LG,1)+(1.0E+00-CCA(LG,1))*LCA_AREA(LG, I)
                  FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)
     &               =CCA(LG,1)/(W_CLOUD(L, NLEVS+1-I)+TOL_MACHINE)
               ELSE
                  W_CLOUD(L, NLEVS+1-I)=LCA_AREA(LG, I)
                  FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)=0.0E+00
               ENDIF
               FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_STRAT)
     &            =1.0E+00-FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ENDIF
!
!        REMOVE VERY THIN CLOUDS TO PREVENT
!        PROBLEMS OF ILL-CONDITIONING.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               MEAN_WATER_CONTENT
     &            =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER)
     &            +(CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
     &            -CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER))
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
               MEAN_ICE_CONTENT
     &            =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE)
     &            +(CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)
     &            -CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE))
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
               IF ( (MEAN_WATER_CONTENT.LT.CONDENSED_LIMIT)
     &              .AND.(MEAN_ICE_CONTENT.LT.CONDENSED_LIMIT) ) THEN
                  W_CLOUD(L, I)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!
!
      ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
!
!        HERE THE CLOUDS ARE SPLIT INTO FOUR SEPARATE TYPES.
!        THE PARTITIONING BETWEEN ICE AND WATER IS REGARDED AS
!        DETERMINING THE AREAS WITHIN THE GRID_BOX COVERED BY
!        ICE OR WATER CLOUD, RATHER THAN AS DETERMINING THE IN-CLOUD
!        MIXING RATIOS. THE GRID-BOX MEAN ICE WATER CONTENTS IN 
!        STRATIFORM CLOUDS MAY BE PREDICTED BY THE ICE MICROPHYSICS 
!        SCHEME OR MAY BE DETERMINED AS A FUNCTION OF THE TEMPERATURE
!        (LSP_FOCWWIL). IN CONVECTIVE CLOUDS THE PARTITIONING MAY BE
!        DONE USING THE SAME FUNCTION, LSP_FOCWWIL, BASED ON A SINGLE
!        TEMPERATURE, OR USING A PARTITION BASED ON THE LOCAL 
!        TEMPERATURE. 
!
!        CONVECTIVE CLOUD:
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
            TOTAL_MASS(L)=0.0E+00
         ENDDO
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG).GE.NLEVS+2-I).AND.
     &              (CCB(LG).LE.NLEVS+1-I) ) THEN
                  TOTAL_MASS(L)=TOTAL_MASS(L)+D_MASS(L, I)
               ENDIF
            ENDDO
         ENDDO
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG).GE.NLEVS+2-I).AND.
     &              (CCB(LG).LE.NLEVS+1-I) ) THEN
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
     &               =CCCWP(LG)/(TOTAL_MASS(L)+TOL_MACHINE)
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)
     &               =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)
     &               =CCCWP(LG)/(TOTAL_MASS(L)+TOL_MACHINE)
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)
     &               =CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)
               ENDIF
            ENDDO
         ENDDO
!
!        STRATIFORM CLOUDS:
!
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I).GT.TOL_TEST) THEN
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &          =CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I).GT.TOL_TEST) THEN
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &          =CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_WATER)
     &               =0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, NLEVS+1-I, IP_CLCMP_ST_ICE)
     &               =0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!        CLOUD FRACTIONS:
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT 
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!
            DO L=1, N_PROFILE
               IF (T(L, NLEVS+2-CCT(I_GATHER(L))).LT.TM) THEN
                  L_GLACIATED_TOP(L)=.TRUE.
               ELSE
                  L_GLACIATED_TOP(L)=.FALSE.
               ENDIF
            ENDDO

         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!
            DO L=1, N_PROFILE
               T_GATHER(L)=T(L, NLEVS+2-CCT(I_GATHER(L)))
            ENDDO
            CALL LSP_FOCWWIL(T_GATHER, N_PROFILE, LIQ_FRAC_CONV)
!
         ENDIF
!
!
         DO I=NLEVS+1-NCLDS, NLEVS
!
            IF (.NOT. L_CLOUD_WATER_PARTITION)
!           PARTITION STRATIFORM CLOUDS USING THE LOCAL TEMPERATURE.
     &        CALL LSP_FOCWWIL(T(1, I), N_PROFILE, LIQ_FRAC)
!
          IF (L_3D_CCA) THEN
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
              CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,NLEVS+1-I)
            ENDDO
          ELSE
            DO L=1, N_PROFILE
            LG=I_GATHER(L)
               IF ( (CCT(LG).GE.NLEVS+2-I).AND.
     &              (CCB(LG).LE.NLEVS+1-I) ) THEN
                CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,1)
               ELSE
                  CONVECTIVE_CLOUD_LAYER(L)=0.0E+00
               ENDIF
            ENDDO
          ENDIF
!
            DO L=1, N_PROFILE
            LG=I_GATHER(L)
               W_CLOUD(L, I)
     &            =CONVECTIVE_CLOUD_LAYER(L)
     &            +(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))
     &            *LCA_AREA(LG, NLEVS+1-I)
!
               IF (L_CLOUD_WATER_PARTITION) THEN
!  PARTITION STRATIFORM CLOUDS USING THE RATIO OF CLOUD WATER CONTENTS.
                 IF (LCA_AREA(LG, NLEVS+1-I).GT.TOL_TEST) THEN
                   LIQ_FRAC(L) = LCCWC1(LG, NLEVS+1-I) /
     &              (LCCWC1(LG, NLEVS+1-I) + LCCWC2(LG, NLEVS+1-I))
                 ELSE
                   LIQ_FRAC(L) = 0.0E+00
                 ENDIF
               ENDIF
!
               IF (L_LOCAL_CNV_PARTITION) THEN
!
!                THE PARTITIONING BETWEEN ICE AND WATER MUST BE
!                RECALCULATED FOR THIS LAYER AS A FUNCTION OF THE
!                LOCAL TEMPERATURE, BUT ICE IS ALLOWED ABOVE THE
!                FREEZING POINT ONLY IF THE TOP OF THE CLOUD IS i
!                GLACIATED.
                 LIQ_FRAC_CONV(L)=MAX(0.0E+00, MIN(1.0E+00
     &              , 1.61E-02*(T(L, I)-TM)+8.9E-01))
                 IF ( (T(L, I).GT.TM).AND.
     &              .NOT.L_GLACIATED_TOP(L) ) THEN
                    LIQ_FRAC_CONV(L)=1.0E+00
                 ENDIF

               ENDIF
!
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)
     &            =LIQ_FRAC(L)*(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))
     &            *LCA_AREA(LG, NLEVS+1-I)
     &            /(W_CLOUD(L, I)+TOL_MACHINE)
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
     &            =(1.0E+00-LIQ_FRAC(L))
     &            *(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))
     &            *LCA_AREA(LG, NLEVS+1-I)/(W_CLOUD(L, I)+TOL_MACHINE)
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)
     &            =LIQ_FRAC_CONV(L)*CONVECTIVE_CLOUD_LAYER(L)
     &            /(W_CLOUD(L, I)+TOL_MACHINE)
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
     &            =(1.0E+00-LIQ_FRAC_CONV(L))*CONVECTIVE_CLOUD_LAYER(L)
     &            /(W_CLOUD(L, I)+TOL_MACHINE)
!
            ENDDO
         ENDDO
!
!
!        REMOVE VERY THIN CLOUDS TO PREVENT
!        PROBLEMS OF ILL-CONDITIONING.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               MEAN_WATER_CONTENT
     &            =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER)
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)
     &            +CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)
               MEAN_ICE_CONTENT
     &            =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE)
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
     &            +CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)
     &            *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
               IF ( (MEAN_WATER_CONTENT.LT.CONDENSED_LIMIT)
     &              .AND.(MEAN_ICE_CONTENT.LT.CONDENSED_LIMIT) ) THEN
                  W_CLOUD(L, I)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
!     EFFECTIVE RADII OF WATER CLOUDS: A MICROPHYSICAL PARAMETRIZATION
!     IS AVAILABLE; OTHERWISE STANDARD VALUES ARE USED.
!
      IF (L_MICROPHYSICS) THEN
!
!        STANDARD VALUES ARE USED FOR ICE CRYSTALS, BUT
!        A PARAMETRIZATION PROVIDED BY UMIST AND MRF
!        IS USED FOR DROPLETS.
!
!        CALCULATE THE DENSITY OF AIR.
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               DENSITY_AIR(L, I)=P(L, I)/(R*T(L, I))
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
            CC_DEPTH(L)=0.0E+00
         ENDDO
!
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
            DO I=NLEVS+2-CCT(LG), NLEVS+1-CCB(LG)
               CC_DEPTH(L)=CC_DEPTH(L)+D_MASS(L, I)/DENSITY_AIR(L, I)
            ENDDO
         ENDDO
!
         CALL R2_RE_MRF_UMIST(N_PROFILE, NLEVS, NCLDS
     &      , I_GATHER
     &      , L_AEROSOL_CCN
     &      , ACCUM_SULPHATE, DISS_SULPHATE
     &      , I_CLOUD_REPRESENTATION
     &      , LAND_G, DENSITY_AIR, CONDENSED_MIX_RAT_BULK, CC_DEPTH
     &      , CONDENSED_DIM_CHAR
     &      , NTOT_DIAG_G
     &      , STRAT_LWC_DIAG_G
     &      , SO4_CCN_DIAG_G
     &      , SULP_DIM1, SULP_DIM2
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &      )
!
!        CONSTRAIN THE SIZES OF DROPLETS TO LIE WITHIN THE RANGE OF
!        VALIDITY OF THE PARAMETRIZATION SCHEME.
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)
     &            =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)
     &            , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)
     &            , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)))
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)
     &            =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)
     &            , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)
     &            , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)))
            ENDDO
         ENDDO
!
!
!        SET MICROPHYSICAL DIAGNOSTICS. WEIGHTS FOR CLOUD CALCULATED
!        HERE ARE USED SOLELY FOR THE MICROPHYSICS AND DO NOT HAVE
!        AN INDEPENDENT MEANING.
!
         IF (WGT_CONV_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     LG=I_GATHER(L)
                     WGT_CONV(LG, I)=W_CLOUD(L, NLEVS+1-I)
     &                  *FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CONV)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     LG=I_GATHER(L)
                     WGT_CONV(LG, I)=W_CLOUD(L, NLEVS+1-I)
     &                  *FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_CW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_CONV_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_CONV(LG, I)
     &               =CONDENSED_DIM_CHAR(L, NLEVS+1-I
     &               , IP_CLCMP_CNV_WATER)
     &               *WGT_CONV(LG, I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF
!
         IF (WGT_STRAT_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     LG=I_GATHER(L)
                     WGT_STRAT(LG, I)=W_CLOUD(L, NLEVS+1-I)
     &                  *FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_STRAT)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     LG=I_GATHER(L)
                     WGT_STRAT(LG, I)=W_CLOUD(L, NLEVS+1-I)
     &                  *FRAC_CLOUD(L, NLEVS+1-I, IP_CLOUD_TYPE_SW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_STRAT(LG, I)
     &               =CONDENSED_DIM_CHAR(L, NLEVS+1-I
     &               , IP_CLCMP_ST_WATER)
     &               *WGT_STRAT(LG, I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF

         IF (LWP_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  LWP_STRAT(LG, I)
     &               =CONDENSED_MIX_RAT_AREA(L, NLEVS+1-I
     &               , IP_CLCMP_ST_WATER)*D_MASS(L, NLEVS+1-I)
     &               *WGT_STRAT(LG, I)
               ENDDO
            ENDDO
         ENDIF

         IF (NTOT_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  NTOT_DIAG(LG, I)
     &               =NTOT_DIAG_G(L, NLEVS+1-I)*WGT_STRAT(LG, I)
               ENDDO
            ENDDO
         ENDIF

         IF (STRAT_LWC_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  STRAT_LWC_DIAG(LG, I)
     &               =STRAT_LWC_DIAG_G(L, NLEVS+1-I)*WGT_STRAT(LG, I)
               ENDDO
            ENDDO
         ENDIF

! Non-cloud diagnostics are "weighted" by the conditional sampling
! weight COND_SAMP_WGT, but as this is 1.0 if the SW radiation is
! active, and 0.0 if it is not, there is no need to actually
! multiply by it.

         IF (COND_SAMP_WGT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  COND_SAMP_WGT(LG, I)=1.0
               ENDDO
            ENDDO
         ENDIF

         IF (SO4_CCN_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  SO4_CCN_DIAG(LG, I)
     &                    =SO4_CCN_DIAG_G(L, NLEVS+1-I)
               ENDDO
            ENDDO
         ENDIF
!
!
      ELSE
!
!        ALL EFFECTIVE RADII ARE SET TO STANDARD VALUES.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)=7.E-6
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)=7.E-6
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     SET THE CHARACTERISTIC DIMENSIONS OF ICE CRYSTALS:
!
!     ICE CRYSTALS IN STRATIFORM CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE).EQ.
     &   IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE).EQ.
     &   IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)
     &            =MIN(7.198755E-04
     &            , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     ICE CRYSTALS IN CONVECTIVE CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE).EQ.
     &   IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE).EQ.
     &   IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)
     &            =MIN(7.198755E-04
     &            , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     CONSTRAIN THE SIZES OF ICE CRYSTALS TO LIE WITHIN THE RANGE
!     OF VALIDITY OF THE PARAMETRIZATION SCHEME.
      DO I=NLEVS+1-NCLDS, NLEVS
         DO L=1, N_PROFILE
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)
     &         =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)
     &         , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)
     &         , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)))
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)
     &         =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)
     &         , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)
     &         , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)))
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND
     &   , I_ST_WATER, I_CNV_WATER, I_ST_ICE, I_CNV_ICE
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , NPD_BAND, NPD_DROP_TYPE, NPD_ICE_TYPE, NPD_CLOUD_PARAMETER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INTERNAL DIMENSIONS TIED TO ALGORITHMS,
!     MOSTLY FOR CLOUDS.
!
      INTEGER
     &     NPD_CLOUD_COMPONENT
!             NUMBER OF COMPONENTS OF CLOUDS
     &   , NPD_CLOUD_TYPE
!             NUMBER OF PERMITTED TYPES OF CLOUDS.
     &   , NPD_CLOUD_REPRESENTATION
!             NUMBER OF PERMITTED REPRESENTATIONS OF CLOUDS.
     &   , NPD_OVERLAP_COEFF
!             NUMBER OF OVERLAP COEFFICIENTS FOR CLOUDS
     &   , NPD_SOURCE_COEFF
!             NUMBER OF COEFFICIENTS FOR TWO-STREAM SOURCES
     &   , NPD_REGION
!             NUMBER OF REGIONS IN A LAYER
!
      PARAMETER(
     &     NPD_CLOUD_COMPONENT=4
     &   , NPD_CLOUD_TYPE=4
     &   , NPD_CLOUD_REPRESENTATION=4
     &   , NPD_OVERLAP_COEFF=18
     &   , NPD_SOURCE_COEFF=2
     &   , NPD_REGION=3
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET COMPONENTS OF CLOUDS.
!
      INTEGER
     &     IP_CLCMP_ST_WATER
!             STRATIFORM WATER DROPLETS
     &   , IP_CLCMP_ST_ICE
!             STRATIFORM ICE CRYSTALS
     &   , IP_CLCMP_CNV_WATER
!             CONVECTIVE WATER DROPLETS
     &   , IP_CLCMP_CNV_ICE
!             CONVECTIVE ICE CRYSTALS
!
      PARAMETER(
     &     IP_CLCMP_ST_WATER=1
     &   , IP_CLCMP_ST_ICE=2
     &   , IP_CLCMP_CNV_WATER=3
     &   , IP_CLCMP_CNV_ICE=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS:
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_BAND
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_DROP_TYPE
!             MAXIMUM NUMBER OF TYPES OF DROPLETS
     &   , NPD_ICE_TYPE
!             MAXIMUM NUMBER OF TYPES OF ICE CRYSTALS
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF PARAMETERS FOR CLOUDS
!
      INTEGER   !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
!     TYPES OF DROPLETS AND CRYSTALS:
      INTEGER   !, INTENT(IN)
     &     I_ST_WATER
!             TYPE OF WATER DROPLETS IN STRATIFORM CLOUDS
     &   , I_CNV_WATER
!             TYPE OF WATER DROPLETS IN CONVECTIVE CLOUDS
     &   , I_ST_ICE
!             TYPE OF ICE CRYSTALS IN STRATIFORM CLOUDS
     &   , I_CNV_ICE
!             TYPE OF ICE CRYSTALS IN CONVECTIVE CLOUDS
!
      LOGICAL   !, INTENT(IN)
     &     L_DROP_TYPE(NPD_DROP_TYPE)
!             FLAGS FOR TYPES OF DROPLET PRESENT
     &   , L_ICE_TYPE(NPD_ICE_TYPE)
!             FLAGS FOR TYPES OF ICE CRYSTAL PRESENT
      INTEGER   !, INTENT(IN)
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
!             PARAMETRIZATIONS OF TYPES OF DROPLETS
     &   , I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             PARAMETRIZATIONS OF TYPES OF ICE CRYSTALS
      REAL      !, INTENT(IN)
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF DROPLETS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)
!             MINIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)
!             MAXIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)
!             MINIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
!
      INTEGER   !, INTENT(OUT)
     &     I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
      REAL      !, INTENT(OUT)
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_CLOUD_COMPONENT, NPD_BAND)
!             COEFFICIENTS FOR PARAMETRIZATION OF CONDENSED PHASES
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , I_SCHEME
!             PARAMETRIZATION SCHEME
!
!     FUNCTIONS CALLED:
      INTEGER
     &     SET_N_CLOUD_PARAMETER
!             FUNCTION TO FIND NUMBER OF PARAMETERS FOR CLOUDS
      EXTERNAL
     &     SET_N_CLOUD_PARAMETER
!
!
!
!     SELECT PARAMETRIZATION FOR WATER IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_WATER.LE.NPD_DROP_TYPE).AND.
     &     (L_DROP_TYPE(I_ST_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_ST_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_ST_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)
     &      =DROP_PARM_MIN_DIM(I_ST_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)
     &      =DROP_PARM_MAX_DIM(I_ST_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE '
     &      , 'OF DROPLET SELECTED IN STRATIFORM WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_WATER, I)
     &         =DROP_PARAMETER_LIST(J, I, I_ST_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR WATER IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_WATER.LE.NPD_DROP_TYPE).AND.
     &     (L_DROP_TYPE(I_CNV_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_CNV_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)
     &      =DROP_PARM_MIN_DIM(I_CNV_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)
     &      =DROP_PARM_MAX_DIM(I_CNV_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE '
     &      , 'OF DROPLET SELECTED IN CONVECTIVE WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_WATER, I)
     &         =DROP_PARAMETER_LIST(J, I, I_CNV_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_ICE.LE.NPD_ICE_TYPE).AND.
     &     (L_ICE_TYPE(I_ST_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_ST_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_ST_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)
     &      =ICE_PARM_MIN_DIM(I_ST_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)
     &      =ICE_PARM_MAX_DIM(I_ST_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE '
     &      , 'OF CRYSTAL SELECTED IN STRATIFORM ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_ICE, I)
     &         =ICE_PARAMETER_LIST(J, I, I_ST_ICE)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_ICE.LE.NPD_ICE_TYPE).AND.
     &     (L_ICE_TYPE(I_CNV_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_CNV_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)
     &      =ICE_PARM_MIN_DIM(I_CNV_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)
     &      =ICE_PARM_MAX_DIM(I_CNV_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE '
     &      , 'OF CRYSTAL SELECTED IN CONVECTIVE ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_ICE, I)
     &         =ICE_PARAMETER_LIST(J, I, I_CNV_ICE)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_AEROSOL_FIELD(IERR
     &   , N_PROFILE, NLEVS, N_AEROSOL, TYPE_AEROSOL
     &   , I_GATHER
     &   , L_CLIMAT_AEROSOL, N_LEVELS_BL
     &   , L_USE_SULPC_DIRECT
     &   , SULP_DIM1, SULP_DIM2
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT
     &   , LAND, LYING_SNOW, PSTAR, AB, BB, TRINDX
     &   , AEROSOL_MIX_RATIO
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
      INTEGER
     &     NPD_AEROSOL_COMPONENT
!             MAXIMUM NUMBER OF AEROSOL COMPONENTS
      INTEGER
     &     IP_WATER_SOLUBLE
!             WATER SOLUBLE AEROSOL
     &   , IP_DUST_LIKE
!             DUST-LIKE AEROSOL
     &   , IP_OCEANIC
!             OCEANIC AEROSOL
     &   , IP_SOOT
!             SOOT AEROSOL
     &   , IP_ASH
!             VOLCANIC ASH
     &   , IP_SULPHURIC
!             SULPHURIC ACID
     &   , IP_ACCUM_SULPHATE
!             ACCUMULATION MODE SULPHATE
     &   , IP_AITKEN_SULPHATE
!             AITKEN MODE SULPHATE
     &   , IP_FRESH_SOOT
     &   , IP_AGED_SOOT
!
!
      PARAMETER(
     &     NPD_AEROSOL_COMPONENT=13
     &   )
      PARAMETER(
     &     IP_WATER_SOLUBLE=1
     &   , IP_DUST_LIKE=2
     &   , IP_OCEANIC=3
     &   , IP_SOOT=4
     &   , IP_ASH=5
     &   , IP_SULPHURIC=6
     &   , IP_ACCUM_SULPHATE=10
     &   , IP_AITKEN_SULPHATE=11
     &   , IP_FRESH_SOOT=12
     &   , IP_AGED_SOOT=13

     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL SPECIES
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , N_LEVELS_BL
!             NUMBER OF LEVELS IN THE BOUNDARY LAYER
     &   , N_AEROSOL
!             NUMBER OF AEROSOLS IN SPECTRAL FILE
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             ACTUAL TYPES OF AEROSOLS
!
!     GATHERING ARRAY:
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
!     FLAG FOR THE CLIMATOLOGICAL AEROSOL DISTRIBUTION.
      LOGICAL      !, INTENT(IN)
     &     L_CLIMAT_AEROSOL
!             FLAG FOR CLIMATOLOGICAL AEROSOL DISTRIBUTION
!
!     VARIABLES FOR THE SULPHUR CYCLE:
      LOGICAL      !, INTENT(IN)
     &     L_USE_SULPC_DIRECT
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
      INTEGER      !, INTENT(IN)
     &     SULP_DIM1,SULP_DIM2
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL      !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIOS OF ACCUMULATION MODE AEROSOL
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIOS OF AITKEN MODE AEROSOL
!
! Declare soot variables:
      LOGICAL L_USE_SOOT_DIRECT !USE DIRECT RAD. EFFECT OF SOOT AEROSOL
      INTEGER SOOT_DIM1,SOOT_DIM2
                !DIMENSIONS FOR SOOT ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_SOOT(SOOT_DIM1, SOOT_DIM2)      ! MMR OF FRESH SOOT
     &   , AGED_SOOT(SOOT_DIM1, SOOT_DIM2)       ! MMR OF AGED SOOT
!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER   !, INTENT(IN)
     &     TRINDX(NPD_FIELD)
!             LAYER BOUNDARY OF TROPOPAUSE
      REAL      !, INTENT(IN)
     &     PSTAR(NPD_FIELD)
!             SURFACE PRESSURES
     &   , AB(NLEVS+1)
!             A AT BOUNDARIES OF LAYERS
     &   , BB(NLEVS+1)
!             B AT BOUNDARIES OF LAYERS
!
!     SURFACE FIELDS
      LOGICAL   !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND SEA MASK
      REAL      !, INTENT(IN)
     &     LYING_SNOW(NPD_FIELD)
!             DEPTH OF LYING SNOW
!
      REAL      !, INTENT(OUT)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_AEROSOL_SPECIES)
!             MIXING RATIOS OF AEROSOLS
!
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , LG
!             INDEX FOR GATHERING
     &   , BLTOP
!             INDEX OF UPPER BOUNDARY OF PLANETARY BOUNDARY LAYER
     &   , I_AEROSOL
!             ACTUAL TYPE OF AEROSOL BEING CONSIDERED
!
!
!     ARRAYS FOR THE CLIMATOLOGICAL AEROSOL MODEL
      LOGICAL
     &     L_IN_CLIMAT(NPD_AEROSOL_COMPONENT)
!             FLAGS TO INDICATE WHICH AEROSOLS ARE INCLUDED IN
!             THE CLIMATOLOGY: THIS MAY BE USED
!             TO ENABLE VARIOUS COMPONENTS TO BE REPLACED BY
!             FULLY PROGNOSTIC SCHEMES.
      INTEGER
     &     I_CLIM_POINTER(NPD_AEROSOL_COMPONENT)
!             POINTERS TO HARD-WIRED INDICES OF THE ORIGINAL
!             CLIMATOLOGICAL AEROSOL MODEL
      REAL
     &     AEROSOL_MIX_RATIO_CLIM(NPD_PROFILE, 0: NPD_LAYER, 5)
!             MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     R2_SET_AERO_CLIM_HADCM3
!
!
!     INITIALIZATION FOR THE CLIMATOLOGICAL AEROSOL MODEL:
!
      DATA L_IN_CLIMAT/NPD_AEROSOL_COMPONENT*.FALSE./
      DATA L_IN_CLIMAT(IP_WATER_SOLUBLE)/.TRUE./
      DATA L_IN_CLIMAT(IP_DUST_LIKE)/.TRUE./
      DATA L_IN_CLIMAT(IP_OCEANIC)/.TRUE./
      DATA L_IN_CLIMAT(IP_SOOT)/.TRUE./
      DATA L_IN_CLIMAT(IP_SULPHURIC)/.TRUE./
!
!     MATCHING OF COMPONENTS TO ORIGINAL HARD-WIRED SETTINGS:
      DATA I_CLIM_POINTER(IP_WATER_SOLUBLE)/1/
      DATA I_CLIM_POINTER(IP_DUST_LIKE)/2/
      DATA I_CLIM_POINTER(IP_OCEANIC)/3/
      DATA I_CLIM_POINTER(IP_SOOT)/4/
      DATA I_CLIM_POINTER(IP_SULPHURIC)/5/
!
!  Use climatological soot if climatological aerosols are on and not
!  using interactive soot.
       L_IN_CLIMAT(IP_SOOT) = L_IN_CLIMAT(IP_SOOT)
     &                     .AND.(.NOT.L_USE_SOOT_DIRECT)
!
!
      IF (L_CLIMAT_AEROSOL) THEN
!
!        SET THE MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!        USED IN THE CLIMATOLOGY OF HADCM3. A SEPARATE SUBROUTINE
!        IS USED TO ENSURE BIT-REPRODUCIBLE RESULTS BY USING
!        EARLIER CODE. THIS COULD BE ALTERED IF A NEW CLIMATOLOGY WERE
!        USED.
!
         CALL R2_SET_AERO_CLIM_HADCM3(N_PROFILE, NLEVS
     &      , I_GATHER
     &      , N_LEVELS_BL
     &      , LAND, LYING_SNOW, PSTAR, AB, BB, TRINDX
     &      , AEROSOL_MIX_RATIO_CLIM
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &      )
!
      ENDIF
!
!
!     THE AEROSOLS REQUIRED BY FOR THE CALCULATION SHOULD HAVE BEEN
!     SELECTED WHEN THE SPECTRAL FILE WAS READ IN. EACH TYPE SHOULD
!     BE SET APPROPRIATELY.
!
      DO J=1, N_AEROSOL
!
         I_AEROSOL=TYPE_AEROSOL(J)
!
         IF (L_CLIMAT_AEROSOL.AND.L_IN_CLIMAT(I_AEROSOL)) THEN
!
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  AEROSOL_MIX_RATIO(L, I, J)
     &               =AEROSOL_MIX_RATIO_CLIM(L, I
     &               , I_CLIM_POINTER(I_AEROSOL))
               ENDDO
            ENDDO
!
         ELSE IF ( (I_AEROSOL.EQ.IP_ACCUM_SULPHATE).AND.
     &      L_USE_SULPC_DIRECT) THEN
!
!           Aerosols related to the sulphur cycle (note that dissolved
!           sulphate does not contribute to the direct effect):
!
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)
     &               =ACCUM_SULPHATE(LG, NLEVS+1-I)
               ENDDO
            ENDDO
!
         ELSE IF ( (I_AEROSOL.EQ.IP_AITKEN_SULPHATE).AND.
     &      L_USE_SULPC_DIRECT) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)
     &               =AITKEN_SULPHATE(LG, NLEVS+1-I)
               ENDDO
            ENDDO
!
         ELSE IF ((I_AEROSOL.EQ.IP_FRESH_SOOT)
     &        .AND.L_USE_SOOT_DIRECT) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=FRESH_SOOT(LG, NLEVS+1-I)
               ENDDO
            ENDDO
!
         ELSE IF ((I_AEROSOL.EQ.IP_AGED_SOOT)
     &        .AND.L_USE_SOOT_DIRECT) THEN
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=AGED_SOOT(LG, NLEVS+1-I)
               ENDDO
            ENDDO
!
!
         ELSE
!
!           The options to the radiation code do not require this
!           aerosol to be considered: its mixing ratio is set to 0.
!           This block of code should not normally be executed,
!           but may be required for ease of including modifications.
!
            DO I=1, NLEVS
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=0.0E+00
               ENDDO
            ENDDO
!
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_AERO_CLIM_HADCM3(N_PROFILE, NLEVS
     &   , I_GATHER
     &   , N_LEVELS_BL
     &   , LAND, LYING_SNOW, PSTAR, AB, BB, TRINDX
     &   , AEROSOL_MIX_RATIO_CLIM
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , N_LEVELS_BL
!             NUMBER OF LEVELS IN THE BOUNDARY LAYER
!
!     GATHERING ARRAY:
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER   !, INTENT(IN)
     &     TRINDX(NPD_FIELD)
!             LAYER BOUNDARY OF TROPOPAUSE
      REAL      !, INTENT(IN)
     &     PSTAR(NPD_FIELD)
!             SURFACE PRESSURES
     &   , AB(NLEVS+1)
!             A AT BOUNDARIES OF LAYERS
     &   , BB(NLEVS+1)
!             B AT BOUNDARIES OF LAYERS
!
!     SURFACE FIELDS
      LOGICAL   !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND-SEA MASK
      REAL      !, INTENT(IN)
     &     LYING_SNOW(NPD_FIELD)
!             DEPTH OF LYING SNOW
!
      REAL      !, INTENT(OUT)
     &     AEROSOL_MIX_RATIO_CLIM(NPD_PROFILE, 0: NPD_LAYER, 5)
!             MIXING RATIOS OF CLIMATOLOGICAL AEROSOLS
!
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , LG
!             INDEX FOR GATHERING
     &   , BLTOP
!             INDEX OF UPPER BOUNDARY OF PLANETARY BOUNDARY LAYER
      REAL
     &     PRESSURE_WT(NPD_FIELD)
!             ARRAY FOR SCALING AEROSOL AMOUNTS FOR DIFFERENT SURFACE
!             PRESSURES
!
!     TOTAL COLUMN MASS (KG M-2) OF EACH AEROSOL SPECIES IN
!     THE BOUNDARY LAYER, THE FREE TROPOSPHERE AND THE STRATOSPHERE
!     RESPECTIVELY. THIS MODEL ASSUMES THAT THERE ARE FIVE AEROSOLS.
      REAL
     &     BL_OCEANMASS(5)
     &   , BL_LANDMASS(5)
     &   , FREETROP_MASS(5)
     &   , STRAT_MASS(5)
!
!     INITIALIZATION FOR THE CLIMATOLOGICAL AEROSOL MODEL
      DATA BL_LANDMASS/2.77579E-5, 6.70018E-5, 0.0, 9.57169E-7, 0.0/
      DATA BL_OCEANMASS/1.07535E-5, 0.0, 2.043167E-4, 0.0, 0.0/
      DATA FREETROP_MASS/3.46974E-6, 8.37523E-6, 0.0, 1.19646E-7, 0.0/
      DATA STRAT_MASS/0.0, 0.0, 0.0, 0.0, 1.86604E-6/
!
!
!
!     TROPOSPHERIC AEROSOL LOADING IS A SIMPLE FUNCTION OF SURFACE
!     PRESSURE: HALVING PSTAR HALVES THE TROPOSPHERIC AEROSOL BURDEN.
!     THE STRATOSPHERIC BURDEN IS INDEPENDENT OF PSTAR.  NOTE THE
!     FACTOR MULTIPLING AEROSOL AMOUNTS USES A REFERENCE PRESSURE
!     OF 1013 mbars.
      DO L=1, N_PROFILE
        PRESSURE_WT(L)=PSTAR(I_GATHER(L))*(1.0/1.013E5)
      END DO
!
!     For each of the 5 aerosol species, the column amount in the
!     boundary layer, free troposphere and stratosphere are known for
!     a standard atmosphere over ocean and land. These can be used
!     to find mixing ratios for the UM by dividing total aerosol by
!     total air mass (and using pressure weighting in the
!     troposphere).
!
!     Firstly, mixing ratios are set for the 5 aerosol species in the
!     stratosphere.
      DO I=1,5
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-TRINDX(LG),I)
     &      =STRAT_MASS(I)*G/
     &      ((AB(TRINDX(LG))+BB(TRINDX(LG))*PSTAR(LG))
     &      -(AB(NLEVS+1)+BB(NLEVS+1)*PSTAR(LG)))
        END DO
      END DO
      DO I=1,5
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
            DO J=(TRINDX(LG)+1),NLEVS
              AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-J,I)=
     &          AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-TRINDX(LG),I)
            END DO
         END DO
       END DO
!      Now, the mixing ratios are set for the 5 aerosol species
!      in the free troposphere.
!      The half-level at the top of the boundary layer is BLTOP
       BLTOP=N_LEVELS_BL+1
       DO I=1,5
         DO L=1, N_PROFILE
           LG=I_GATHER(L)
           AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-BLTOP,I)
     &       =FREETROP_MASS(I)*G*
     &       PRESSURE_WT(L)/((AB(BLTOP)+BB(BLTOP)*PSTAR(LG))-
     &       (AB(TRINDX(LG))+BB(TRINDX(LG))*PSTAR(LG)))
         END DO
       END DO
       DO L=1, N_PROFILE
         LG=I_GATHER(L)
         IF ((BLTOP+1).LE.(TRINDX(LG)-1)) THEN
           DO I=1,5
             DO J=(BLTOP+1),(TRINDX(LG)-1)
               AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-J,I)=
     &           AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-BLTOP,I)
             END DO
           END DO
         END IF
       END DO
!
!      Now, the boundary layer mixing ratios are set for the
!      5 aerosol species. A continental aerosol is used over most land
!      areas, but not over ice sheets, which are identified by the
!      criterion used in the boundary layer scheme that the mass of
!      lying snow exceeds 5000 kgm-2. Over ice sheets a maritime
!      aerosol is used.
       DO I=1,5
         DO L=1, N_PROFILE
           LG=I_GATHER(L)
           IF ( LAND(LG).AND.(LYING_SNOW(LG).LT.5.0E+03) ) THEN
            AEROSOL_MIX_RATIO_CLIM(L,NLEVS+2-BLTOP,I)
     &        =BL_LANDMASS(I)*G*PRESSURE_WT(L)
     &        /(PSTAR(LG)-(AB(BLTOP)+BB(BLTOP)*PSTAR(LG)))
           ELSE
            AEROSOL_MIX_RATIO_CLIM(L,NLEVS+2-BLTOP,I)
     &        =BL_OCEANMASS(I)*G*PRESSURE_WT(L)
     &        /(PSTAR(LG)-(AB(BLTOP)+BB(BLTOP)*PSTAR(LG)))
           END IF
         END DO
       END DO
       DO I=1,5
         DO J=1,(BLTOP-2)
           DO L=1, N_PROFILE
             AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-J,I)=
     &         AEROSOL_MIX_RATIO_CLIM(L,NLEVS+2-BLTOP,I)
           END DO
         END DO
       END DO
!
!
!
      RETURN
      END
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, NLEVS, NCLDS
     &   , I_CLOUD, W_CLOUD, TOTAL_CLOUD_COVER
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DECLARATION OF ARRAY SIZES.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     COMDECKS INCLUDED
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO DEFINE REFERENCE NUMBERS FOR CLOUD SCHEMES.
!
      INTEGER
     &     IP_CLOUD_MIX_MAX
!             MAXIMUM/RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_MIX_RANDOM
!             RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_COLUMN_MAX
!             MAXIMUM OVERLAP IN A COLUMN MODEL
     &   , IP_CLOUD_CLEAR
!             CLEAR COLUMN
     &   , IP_CLOUD_TRIPLE
!             MIXED COLUMN WITH SPLIT BETWEEN 
!             CONVECTIVE AND LAYER CLOUD.
!
      PARAMETER(
     &     IP_CLOUD_MIX_MAX=2
     &   , IP_CLOUD_MIX_RANDOM=4
     &   , IP_CLOUD_COLUMN_MAX=3
     &   , IP_CLOUD_CLEAR=5
     &   , IP_CLOUD_TRIPLE=6
     &   )
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF LAYERS
     &   , NCLDS
!             NUMBER OF CLOUDY LAYERS
     &   , I_CLOUD
!             CLOUD SCHEME EMPLOYED
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
      REAL      !, INTENT(OUT)
     &     TOTAL_CLOUD_COVER(NPD_PROFILE)
!             TOTAL CLOUD COVER
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
!
!
!
!     DIFFERENT OVERLAP ASSUMPTIONS ARE CODED INTO EACH SOLVER.
!
      IF (I_CLOUD.EQ.IP_CLOUD_MIX_MAX) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, NLEVS+1-NCLDS)
         ENDDO
         DO I=NLEVS+1-NCLDS, NLEVS-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1).GT.W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)
     &               *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00
         ENDDO
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)
     &            *(1.0E+00-W_CLOUD(L, I))
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_COLUMN_MAX) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=MAX(TOTAL_CLOUD_COVER(L)
     &            , W_CLOUD(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, NLEVS+1-NCLDS)
         ENDDO
         DO I=NLEVS+1-NCLDS, NLEVS-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1).GT.W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)
     &               *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_CLEAR) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius 
!                                               removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_RE_MRF_UMIST(N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , L_AEROSOL_CCN, ACCUM_SULPHATE, DISS_SULPHATE
     &   , I_CLOUD_REPRESENTATION
     &   , LAND_G, DENSITY_AIR, CONDENSED_MIX_RATIO, CC_DEPTH
     &   , CONDENSED_RE
     &   , NTOT_DIAG_G
     &   , STRAT_LWC_DIAG_G
     &   , SO4_CCN_DIAG_G
     &   , SULP_DIM1, SULP_DIM2
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED:
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

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
!
! Description:
!
!  Contains various parameters for the SW effective radius
!  parametrisation, defined for land and sea areas.
!
!  NTOT is the number concentration (m-3) of activated CCN;
!  KPARAM is the ratio of volume mean radius to effective radius;
!  DCONRE is the effective radius (m) for deep convective clouds.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!
      REAL NTOT_LAND,
     &     NTOT_SEA,
     &     KPARAM_LAND,
     &     KPARAM_SEA,
     &     DCONRE_LAND,
     &     DCONRE_SEA

      PARAMETER ( NTOT_LAND = 6.0E08,
     &            NTOT_SEA = 1.5E08,
     &            KPARAM_LAND = 0.67,
     &            KPARAM_SEA = 0.80,
     &            DCONRE_LAND = 9.5E-06,
     &            DCONRE_SEA = 13.5E-06 )
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INTERNAL DIMENSIONS TIED TO ALGORITHMS,
!     MOSTLY FOR CLOUDS.
!
      INTEGER
     &     NPD_CLOUD_COMPONENT
!             NUMBER OF COMPONENTS OF CLOUDS
     &   , NPD_CLOUD_TYPE
!             NUMBER OF PERMITTED TYPES OF CLOUDS.
     &   , NPD_CLOUD_REPRESENTATION
!             NUMBER OF PERMITTED REPRESENTATIONS OF CLOUDS.
     &   , NPD_OVERLAP_COEFF
!             NUMBER OF OVERLAP COEFFICIENTS FOR CLOUDS
     &   , NPD_SOURCE_COEFF
!             NUMBER OF COEFFICIENTS FOR TWO-STREAM SOURCES
     &   , NPD_REGION
!             NUMBER OF REGIONS IN A LAYER
!
      PARAMETER(
     &     NPD_CLOUD_COMPONENT=4
     &   , NPD_CLOUD_TYPE=4
     &   , NPD_CLOUD_REPRESENTATION=4
     &   , NPD_OVERLAP_COEFF=18
     &   , NPD_SOURCE_COEFF=2
     &   , NPD_REGION=3
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET COMPONENTS OF CLOUDS.
!
      INTEGER
     &     IP_CLCMP_ST_WATER
!             STRATIFORM WATER DROPLETS
     &   , IP_CLCMP_ST_ICE
!             STRATIFORM ICE CRYSTALS
     &   , IP_CLCMP_CNV_WATER
!             CONVECTIVE WATER DROPLETS
     &   , IP_CLCMP_CNV_ICE
!             CONVECTIVE ICE CRYSTALS
!
      PARAMETER(
     &     IP_CLCMP_ST_WATER=1
     &   , IP_CLCMP_ST_ICE=2
     &   , IP_CLCMP_CNV_WATER=3
     &   , IP_CLCMP_CNV_ICE=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET REPRESENTATIONS OF CLOUDS.
!
!     REPRESENTATIONS
      INTEGER
     &     IP_CLOUD_HOMOGEN
!             ALL COMPONENTS ARE MIXED HOMOGENEOUSLY
     &   , IP_CLOUD_ICE_WATER
!             ICE AND WATER CLOUDS ARE TREATED SEPARATELY
     &   , IP_CLOUD_CONV_STRAT
!             CLOUDS ARE DIVIDED INTO HOMOGENEOUSLY MIXED
!             STRATIFORM AND CONVECTIVE PARTS
     &   , IP_CLOUD_CSIW
!             CLOUDS DIVIDED INTO ICE AND WATER PHASES AND
!             INTO STRATIFORM AND CONVECTIVE COMPONENTS.
!
      PARAMETER(
     &     IP_CLOUD_HOMOGEN=1
     &   , IP_CLOUD_ICE_WATER=2
     &   , IP_CLOUD_CONV_STRAT=3
     &   , IP_CLOUD_CSIW=4
     &   )
!
!     TYPES OF CLOUDS:
      INTEGER
     &     NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)
!             NUMBER OF TYPE OF CLOUDS IN REPRESENTATION
!
      INTEGER
     &     IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT
     &       , NPD_CLOUD_REPRESENTATION)
!            MAP OF COMPONENTS CONTRIBUTING TO TYPES OF CLOUDS
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS:
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF INPUT FIELDS TO THE RADIATION
     &   , NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL SPECIES
     &   , SULP_DIM1
!             1ST DIMENSION OF ARRAYS OF SULPHATE
     &   , SULP_DIM2
!             2ND DIMENSION OF ARRAYS OF SULPHATE
!
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , NLEVS
!             NUMBER OF LEVELS
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
!
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
      LOGICAL   !, INTENT(IN)
     &     LAND_G(NPD_PROFILE)
!             GATHERED MASK FOR LAND POINTS
      INTEGER   !, INTENT(IN)
     &     I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
!
!     VARIABLES FOR AEROSOLS
      LOGICAL   !, INTENT(IN)
     &     L_AEROSOL_CCN
!             FLAG TO USE AEROSOLS TO FIND CCN.
      REAL      !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF ACCUMULATION MODE SULPHATE
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF DISSOLVED SULPHATE
!
      REAL      !, INTENT(IN)
     &     DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             DENSITY OF AIR
!
      REAL      !, INTENT(IN)
     &     CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             MIXING RATIOS OF CONDENSED SPECIES
     &   , CC_DEPTH(NPD_PROFILE)
!             DEPTH OF CONVECTIVE CLOUD
!
      REAL      !, INTENT(OUT)
     &     CONDENSED_RE(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CONDENSED COMPONENTS OF CLOUDS
!
      REAL      !, INTENT(OUT)
     &     NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
     &   , STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
     &   , SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     TOTAL_MIX_RATIO_ST(NPD_PROFILE)
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
     &   , TOTAL_MIX_RATIO_CNV(NPD_PROFILE)
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
!
      REAL
     &     N_DROP(NPD_PROFILE, NPD_LAYER)
!             NUMBER DENSITY OF DROPLETS
     &   , KPARAM
!             RATIO OF CUBES OF VOLUME RADIUS TO EFFECTIVE RADIUS
!
!     FIXED CONSTANTS OF THE PARAMETRIZATION:
      REAL
     &     DEEP_CONVECTIVE_CLOUD
!             THRESHOLD VALUE FOR DEEP CONVECTIVE CLOUD
      PARAMETER(
     &     DEEP_CONVECTIVE_CLOUD=5.0E+02
     &   )
!
!
!
!     CALCULATE THE NUMBER DENSITY OF DROPLETS
      CALL R2_FIND_NUMBER_DROP(N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , DENSITY_AIR, L_AEROSOL_CCN
     &   , ACCUM_SULPHATE, DISS_SULPHATE
     &   , LAND_G
     &   , N_DROP
     &   , SO4_CCN_DIAG_G
     &   , SULP_DIM1, SULP_DIM2
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &   )
!
      DO I=NLEVS+1-NCLDS, NLEVS
!
!        FIND THE TOTAL MIXING RATIO OF WATER SUBSTANCE IN THE CLOUD
!        AS IMPLIED BY THE REPRESENTATION.
         IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)
     &            +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_ICE)
               TOTAL_MIX_RATIO_CNV(L)
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)
     &            +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_ICE)
            ENDDO
         ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)
               TOTAL_MIX_RATIO_CNV(L)
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)
            ENDDO
         ENDIF
         DO L=1, N_PROFILE
            IF (LAND_G(L)) THEN
               KPARAM=KPARAM_LAND
            ELSE
               KPARAM=KPARAM_SEA
            ENDIF
            CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)
     &         =(3.0E+00*TOTAL_MIX_RATIO_CNV(L)*DENSITY_AIR(L, I)
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
     &         **(1.0E+00/3.0E+00)
            CONDENSED_RE(L, I, IP_CLCMP_ST_WATER)
     &         =(3.0E+00*TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
     &         **(1.0E+00/3.0E+00)
         ENDDO
         DO L=1, N_PROFILE
            NTOT_DIAG_G(L, I)=N_DROP(L, I)*1.0E-06
            STRAT_LWC_DIAG_G(L, I)
     &         =TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)*1.0E03
         ENDDO
      ENDDO
!
!     RESET THE EFFECTIVE RADII FOR DEEP CONVECTIVE CLOUDS.
      DO I=NLEVS+1-NCLDS, NLEVS
         DO L=1, N_PROFILE
            IF (CC_DEPTH(L).GT.DEEP_CONVECTIVE_CLOUD) THEN
               IF (LAND_G(L)) THEN
                  CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_LAND
               ELSE
                  CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_SEA
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to calculate the number density of droplets.
!
! Purpose:
!   The number density of cloud droplets is calculated.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               The number of CCN now
!                                               depends on the
!                                               dissolved sulphate as
!                                               well as the accumulation
!                                               mode sulphate.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_FIND_NUMBER_DROP(N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , DENSITY_AIR, L_AEROSOL_CCN
     &   , ACCUM_SULPHATE, DISS_SULPHATE
     &   , LAND_G
     &   , N_DROP
     &   , SO4_CCN_DIAG_G
     &   , SULP_DIM1, SULP_DIM2
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED:
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

!
! Description:
!
!  Contains various parameters for the SW effective radius
!  parametrisation, defined for land and sea areas.
!
!  NTOT is the number concentration (m-3) of activated CCN;
!  KPARAM is the ratio of volume mean radius to effective radius;
!  DCONRE is the effective radius (m) for deep convective clouds.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!
      REAL NTOT_LAND,
     &     NTOT_SEA,
     &     KPARAM_LAND,
     &     KPARAM_SEA,
     &     DCONRE_LAND,
     &     DCONRE_SEA

      PARAMETER ( NTOT_LAND = 6.0E08,
     &            NTOT_SEA = 1.5E08,
     &            KPARAM_LAND = 0.67,
     &            KPARAM_SEA = 0.80,
     &            DCONRE_LAND = 9.5E-06,
     &            DCONRE_SEA = 13.5E-06 )
!
!
!     DUMMY ARGUMENTS:
!
!     SIZES OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF INPUT FIELDS
     &   , NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL SPECIES
     &   , SULP_DIM1
!             1ST DIMENSION OF ARRAYS OF SULPHATE
     &   , SULP_DIM2
!             2ND DIMENSION OF ARRAYS OF SULPHATE
      INTEGER   !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , NLEVS
!             NUMBER OF LEVELS
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
      LOGICAL   !, INTENT(IN)
     &     L_AEROSOL_CCN
!             FLAG TO USE AEROSOLS TO FIND CCN
     &   , LAND_G(NPD_PROFILE)
!             GATHERED MASK FOR LAND POINTS
      REAL      !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF ACCUMULATION-MODE SULPHATE
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF DISSOLVED SULPHATE
     &   , DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             DENSITY OF AIR
!
      REAL      !, INTENT(OUT)
     &     N_DROP(NPD_PROFILE, NPD_LAYER)
!             NUMBER DENSITY OF DROPLETS
!
      REAL
     &     SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             SO4 CCN MASS CONC DIAGNOSTIC ARRAY (GATHERED)
!
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     PARTICLE_VOLUME
!             MEAN VOLUME OF A PARTICLE
     &   , N_CCN
!             NUMBER DENSITY OF CCN
!
      REAL
     &     RADIUS_0
!             MEDIAN RADIUS OF LOG-NORMAL DISTRIBUTION
     &   , SIGMA_0
!             GEOMETRIC STANDARD DEVIATION
     &   , DENSITY_SULPHATE
!             DENSITY OF SULPHATE AEROSOL
      PARAMETER(
     &     RADIUS_0=5.0E-08
     &   , SIGMA_0=2.0
     &   , DENSITY_SULPHATE=1.769E+03
     &   )
!
!
!
      IF (L_AEROSOL_CCN) THEN
!
!        IF AEROSOLS ARE INCLUDED THE NUMBER OF CCN IS FOUND FROM THE
!        CONCENTRATION OF ACCUMULATION-MODE AND DISSOLVED SULPHATE.
!        NOTE THAT IN PRINCIPLE EACH MODE MIGHT HAVE A DIFFERENT
!        DENSITY AND SIZE DISTRIBUTION.
!        THE DROPLET NUMBER CONCENTRATION IS HELD TO A MINIMUM
!        VALUE OF 5.0E+06 (5cm-3).
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               PARTICLE_VOLUME=(4.0E+00*PI/3.0E+00)*RADIUS_0**3
     &            *EXP(4.5E+00*(LOG(SIGMA_0))**2)
               N_CCN=(ACCUM_SULPHATE(I_GATHER(L), NLEVS+1-I)
     &            +DISS_SULPHATE(I_GATHER(L), NLEVS+1-I))
     &            *DENSITY_AIR(L, I)
     &            /(DENSITY_SULPHATE*PARTICLE_VOLUME)
               N_DROP(L, I)=3.75E+08*(1.0E+00-EXP(-2.5E-9*N_CCN))
               IF (N_DROP(L, I) .LT. 5.0E+06) N_DROP(L, I)=5.0E+06
!              CONVERT THE MASS MIXING RATIOS FROM (NH4)2SO4
!              TO MASS PER UNIT VOLUME OF SO4(IN MICROGRAMMES
!              PER CUBIC METRE) FOR DIAGNOSTIC PURPOSES.
               SO4_CCN_DIAG_G(L, I)=
     &                  (ACCUM_SULPHATE(I_GATHER(L), NLEVS+1-I)
     &                  +DISS_SULPHATE(I_GATHER(L), NLEVS+1-I))
     &                  * DENSITY_AIR(L, I) * (96./132.) * 1.0E+09
            ENDDO
         ENDDO
!
      ELSE
!
!        WITHOUT AEROSOLS THE NUMBERS OF DROPLETS ARE FIXED.
!
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               IF (LAND_G(L)) THEN
                  N_DROP(L, I)=NTOT_LAND
               ELSE
                  N_DROP(L, I)=NTOT_SEA
               ENDIF
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COMPARE_PROC(IERR, L_PRESENT
     &   , L_RAYLEIGH_PERMITTED, L_GAS_PERMITTED, L_CONTINUUM_PERMITTED
     &   , L_DROP_PERMITTED, L_AEROSOL_PERMITTED
     &   , L_AEROSOL_CCN_PERMITTED, L_ICE_PERMITTED
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT
     &    ,L_USE_SOOT_DIRECT
     &   , L_CLIMAT_AEROSOL
     &   , L_RAYLEIGH, L_GAS, L_CONTINUUM
     &   , L_DROP, L_AEROSOL, L_AEROSOL_CCN, L_ICE
     &   , NPD_TYPE
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS:
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     NPD_TYPE
!             NUMBER OF TYPES OF SPECTRAL DATA
!
      LOGICAL   !, INTENT(IN)
     &     L_PRESENT(0: NPD_TYPE)
!             ARRAY INDICATING BLOCKS OF DATA PRESENT
!             IN THE SPECTRAL FILE.
!
!     PROCESSES PERMITTED WITHIN THE UNIFIED MODEL.
      LOGICAL   !, INTENT(IN)
     &     L_RAYLEIGH_PERMITTED
!             RAYLEIGH SCATTERING PERMITTED IN THE MODEL
     &   , L_GAS_PERMITTED
!             GASEOUS ABSORPTION PERMITTED IN THE MODEL
     &   , L_CONTINUUM_PERMITTED
!             CONTINUUM ABSORPTION PERMITTED IN THE MODEL
     &   , L_DROP_PERMITTED
!             CLOUD DROPLET EXTINCTION PERMITTED IN THE MODEL
     &   , L_AEROSOL_PERMITTED
!             AEROSOL EXTINCTION PERMITTED IN THE MODEL
     &   , L_AEROSOL_CCN_PERMITTED
!             DETERMINATION OF CCN FROM AEROSOLS PERMITTED IN THE MODEL
     &   , L_ICE_PERMITTED
!             ICE EXTINCTION PERMITTED IN THE MODEL
!
!     OPTIONS PASSED IN
      LOGICAL
     &     L_USE_SULPC_DIRECT
!             LOGICAL TO USE SULPHUR CYCLE FOR THE DIRECT EFFECT
     &   , L_USE_SULPC_INDIRECT
!             LOGICAL TO USE SULPHUR CYCLE FOR THE INDIRECT EFFECT
     &    ,L_USE_SOOT_DIRECT
!             LOGICAL TO USE DIRECT RADIATIVE EFFECT DUE TO SOOT
     &   , L_CLIMAT_AEROSOL
!             LOGICAL TO USE CLIMATOLOGICAL AEROSOL MODEL
!
!     PROCESSES TO BE ENABLED IN THE RUN.
      LOGICAL   !, INTENT(OUT)
     &     L_RAYLEIGH
!             RAYLEIGH SCATTERING TO BE ENABLED IN THE RUN
     &   , L_GAS
!             GASEOUS ABSORPTION TO BE ENABLED IN THE RUN
     &   , L_CONTINUUM
!             CONTINUUM ABSORPTION TO BE ENABLED IN THE RUN
     &   , L_DROP
!             CLOUD DROPLET EXTINCTION TO BE ENABLED IN THE RUN
     &   , L_AEROSOL
!             AEROSOL EXTINCTION TO BE ENABLED IN THE RUN
     &   , L_AEROSOL_CCN
!             DETERMINATION OF CCN FROM AEROSOL TO BE ENABLED IN THE RUN
     &   , L_ICE
!             ICE EXTINCTION TO BE ENABLED IN THE RUN
!
!
!
!     EACH OPTICAL PROCESS INCLUDED IN THE RADIATION CODE MAY BE
!     PERMITTED OR DENIED IN THE UNIFIED MODEL, DEPENDING ON THE
!     PRESENCE OF SUPPORTING CODE. TO BE ENABLED IN A RUN AN OPTICAL
!     PROCESS MUST BE PERMITTED IN THE UNIFIED MODEL AND HAVE
!     SUITABLE SPECTRAL DATA.
      L_RAYLEIGH=L_RAYLEIGH_PERMITTED.AND.L_PRESENT(3)
      L_GAS=L_GAS_PERMITTED.AND.L_PRESENT(5)
      L_CONTINUUM=L_CONTINUUM_PERMITTED.AND.L_PRESENT(9)
      L_DROP=L_DROP_PERMITTED.AND.L_PRESENT(10)
      L_ICE=L_ICE_PERMITTED.AND.L_PRESENT(12)
!
!     SET THE CONTROLLING FLAG FOR THE DIRECT RADIATIVE EFFECTS OF
!     AEROSOLS.
      IF (L_AEROSOL_PERMITTED) THEN
!        SET THE FLAG AND THEN CHECK THE SPECTRAL FILE.
         L_AEROSOL=L_USE_SULPC_DIRECT.OR.L_CLIMAT_AEROSOL
     &    .OR. L_USE_SOOT_DIRECT
         IF (L_AEROSOL.AND.(.NOT.L_PRESENT(11))) THEN
            WRITE(IU_ERR, '(/A, /A)')
     &         '*** ERROR: THE SPECTRAL FILE CONTAINS NO DATA '
     &         //'FOR AEROSOLS.', 'SUCH DATA ARE REQUIRED FOR THE '
     &         //'DIRECT EFFECT.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ELSE
!        CHECK THAT AEROSOLS HAVE NOT BEEN REQUESTED
!        WHEN NOT PERMITTED.
         IF (L_USE_SULPC_DIRECT
     &   .OR.L_CLIMAT_AEROSOL
     &   .OR.L_USE_SOOT_DIRECT) THEN
            WRITE(IU_ERR, '(/A, /A)')
     &         '*** ERROR: THE DIRECT EFFECTS AEROSOLS ARE NOT '
     &         , 'PERMITTED IN THIS CONFIGURATION OF THE '
     &         //'RADIATION CODE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!     SET THE CONTROLLING FLAG FOR THE INDIRECT EFFECTS OF AEROSOLS.
!     AT PRESENT THIS DEPENDS SOLELY ON THE SULPHUR CYCLE.
      L_AEROSOL_CCN=L_USE_SULPC_INDIRECT
!
!
!
      RETURN
      END
