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
!!! Subroutine SOILCARB -----------------------------------------------
!!!
!!! Purpose : Updates carbon contents of the soil.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE SOILCARB (LAND_FIELD,TRIF_PTS,TRIF_INDEX
     &,                    FORW,GAMMA,LIT_C_T,RESP_S,CS)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,TRIF_PTS                   ! IN Number of points on which 
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_FIELD)     ! IN Indices of land points on 
!                                 !    which TRIFFID may operate
     &,L,T                        ! WORK Loop counters

      REAL
     & FORW                       ! IN Forward timestep weighting.
     &,GAMMA                      ! IN Inverse timestep (/360days).
     &,LIT_C_T(LAND_FIELD)        ! IN Total carbon litter 
C                                 !    (kg C/m2/360days).
     &,RESP_S(LAND_FIELD)         ! INOUT Soil respiration 
C                                 !    (kg C/m2/360days).
     &,CS(LAND_FIELD)             ! INOUT Soil carbon (kg C/m2).
     &,DCS(LAND_FIELD)            ! WORK Increment to the soil carbon
C                                 !      (kg C/m2).
     &,DPC_DCS(LAND_FIELD)        ! WORK Rate of change of PC with
C                                 !      soil carbon (/360days).
     &,PC(LAND_FIELD)             ! WORK Net carbon accumulation in
C                                 !      the soil (kg C/m2/360days).



      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 

C----------------------------------------------------------------------
C Diagnose the net local carbon flux into the soil
C----------------------------------------------------------------------
        PC(L) = LIT_C_T(L)-RESP_S(L)

C----------------------------------------------------------------------
C Variables required for the implicit and equilibrium calculations
C----------------------------------------------------------------------
        DPC_DCS(L) = RESP_S(L)/CS(L)

C----------------------------------------------------------------------
C Save current value of soil carbon
C----------------------------------------------------------------------
        DCS(L) = CS(L)

      ENDDO


C----------------------------------------------------------------------
C Update soil carbon
C----------------------------------------------------------------------
      CALL DECAY (LAND_FIELD,TRIF_PTS,TRIF_INDEX
     &,           DPC_DCS,FORW,GAMMA,PC,CS)

C----------------------------------------------------------------------
C Apply implicit correction to the soil respiration rate.
C----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 

        DCS(L) = CS(L) - DCS(L)
        RESP_S(L) = RESP_S(L) + FORW*DPC_DCS(L)*DCS(L)

      ENDDO

      RETURN
      END
