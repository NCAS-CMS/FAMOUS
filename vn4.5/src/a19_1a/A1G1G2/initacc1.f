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
! Initialises accumulated carbon fluxes to zero if new calling period
!
! Subroutine Interface:
      SUBROUTINE INIT_ACC(LAND_FIELD,
     &       NPP_PFT_ACC,G_LEAF_PHEN_PFT_ACC,
     &       RESP_W_PFT_ACC,RESP_S_ACC,ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Resets accumulation prognostics to zero if a new TRIFFID calling
!   period is starting.  This routine is needed when starting an NRUN
!   from an initial dump created in either of the following situations:
!
!   i)  Initial dump created from a non-TRIFFID run
!
!   ii) Initial dump created in a TRIFFID run mid-way through a TRIFFID
!       calling period.  The NRUN may re-start at the same point within
!       this calling period and continue with the accumulation already
!       part-completed in this dump; in this case this routine will not
!       be used.  Alternatively, the NRUN may start a new calling
!       period, in which case the accumulation must begin; this routine
!       allows this by re-setting the relevant prognostics to zero.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4   10/10/97   Original code.  Richard Betts
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
C Arguments

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
     & LAND_FIELD                          ! IN number of land points


      REAL
     & NPP_PFT_ACC(LAND_FIELD,NPFT)        !INOUT Accumulated NPP on
C                                        !      Plant Functional Types
     &,G_LEAF_PHEN_PFT_ACC(LAND_FIELD,NPFT)!INOUT Accum. phenological
C                                        !      leaf turnover rate PFTs
     &,RESP_W_PFT_ACC(LAND_FIELD,NPFT)     !INOUT Accumulated wood
C                                        !      respiration on PFTs
     &,RESP_S_ACC(LAND_FIELD)              !INOUT Accumulated soil resp.


      INTEGER
     & L                       ! Loop counter for land points
     &,N                       ! Loop counter for plant functional types

      INTEGER ICODE            ! Work - Internal return code
      CHARACTER*80 CMESSAGE    ! Work - Internal error message


      WRITE (6,*)
     & 'INIT_ACC: setting accumulation prognostics to zero'

      DO L=1,LAND_FIELD
        DO N=1,NPFT
          NPP_PFT_ACC(L,N) = 0.0
          G_LEAF_PHEN_PFT_ACC(L,N) = 0.0
          RESP_W_PFT_ACC(L,N) = 0.0
        ENDDO
        RESP_S_ACC(L) = 0.0
      ENDDO

      RETURN
      END
