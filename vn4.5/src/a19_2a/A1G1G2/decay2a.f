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
!!! Subroutine DECAY --------------------------------------------------
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
      SUBROUTINE DECAY (LAND_FIELD,TRIF_PTS,TRIF_INDEX
     &,                 DPC_DCS,FORW,GAMMA,PC,CS)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,TRIF_PTS                   ! IN Number of points on which 
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_FIELD)     ! IN Indices of land points on 
!                                 !    which TRIFFID may operate
     &,L,T                        ! WORK Loop counters

      REAL
     & DPC_DCS(LAND_FIELD)        ! IN Rate of change of PC with
C                                 !    soil carbon (yr).
     &,FORW                       ! IN Forward timestep weighting.
     &,GAMMA                      ! IN Inverse timestep (/360days).
     &,PC(LAND_FIELD)             ! IN Net carbon flux into the
C                                 !    soil (kg C/m2/360days).
     &,CS(LAND_FIELD)             ! INOUT Soil carbon (kg C/m2).
     &,DENOM                      ! WORK Denominator of update
C                                 !      equation.
     &,DENOM_MIN                  ! WORK Minimum value for the
C                                 !      denominator of the update
C                                 !      equation. Ensures that
C                                 !      gradient descent does not
C                                 !      lead to an unstable solution.
     &,NUMER                      ! WORK Numerator of the update
C                                 !      equation.
C----------------------------------------------------------------------
C Local parameters
C----------------------------------------------------------------------
      REAL
     + CS_MIN                     ! Minimum soil carbon (kg C/m2).
      PARAMETER(CS_MIN = 1.0E-6)
      INTEGER
     + ITER_EQ                    ! Number of TRIFFID iterations for
C                                 ! gradient descent to equilibrium.
      REAL
     + GAMMA_EQ                   ! Inverse timestep for gradient
C                                 ! descent to equilibrium (/360days).
      PARAMETER(GAMMA_EQ = 1.0E-4, ITER_EQ = 10)


      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 

        NUMER = PC(L)
        DENOM = GAMMA+FORW*DPC_DCS(L)
        DENOM_MIN = GAMMA_EQ
        DENOM = MAX(DENOM,DENOM_MIN)

        CS(L) = CS(L)+NUMER/DENOM

        CS(L) = MAX(CS_MIN,CS(L))

      ENDDO

      RETURN
      END
