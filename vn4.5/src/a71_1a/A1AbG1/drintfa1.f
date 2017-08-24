C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Subroutine DERV_INTF_A : Calculates Interface array dimensions.
!
! Subroutine Interface :
!
      SUBROUTINE DERV_INTF_A (TOT_LEN_INTFA_P,TOT_LEN_INTFA_U,
     &           MAX_INTF_P_LEVELS,N_INTF_A,U_FIELD,U_FIELD_INTFA)

      IMPLICIT NONE
!
! Description : Calculate array dimensions for boundary data output.
!
! Method : Reads in INTFCNSTA namelist to get grid dimensions of
!          interface area. Calculates array dimensions for boundary
!          data. Also sets dimensions to 1 if no interface areas
!          required to prevent zero dynamic allocation.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.5    03/08/98  Original Code
!
! Code Description :
! Language : FORTRAN 77 + common extensions
!
! Declarations :

!     Arguments
      Integer TOT_LEN_INTFA_P   ! OUT  Total length of interface data
                                !      on P grid
      Integer TOT_LEN_INTFA_U   ! OUT  Total length of interface data
                                !      on U grid
      Integer MAX_INTF_P_LEVELS ! OUT  Max no of levels for all areas
      Integer N_INTF_A          ! IN   No of interface areas
      Integer U_FIELD           ! IN   Dimension of U_FIELD
      Integer U_FIELD_INTFA     ! OUT  Dimension of U_FIELD for dynamic
                                !      allocation.

C*L================ COMDECK CMAXSIZE ==========================
C   Description:
C     This COMDECK contains maximum sizes for dimensioning arrays
C   of model constants whose sizes are configuration dependent. This
C   allows constants to be read in from a NAMELIST file and maintain
C   the flexibility of dynamic allocation for primary variables. The
C   maximum sizes should agree with the maximum sizes implicit in the
C   front-end User Interface.
C
CLL
CLL  Model            Modification history:
CLL version  Date
CLL 3.2   26/03/93  New COMDECK. Author R.Rawlins
CLL  3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
CLL                 in macro-tasked calls to SWRAD, LWRAD & CONVECT.
CLL                 Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
CLL  3.5  22/05/95  Add MAX_N_INTF. D. Robinson
CLL  4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.

CLL
C
C

C Define Parameters:
      INTEGER  MAX_P_LEVELS     ! Maximum no. of p levels
        PARAMETER (MAX_P_LEVELS = 99  )
      INTEGER  MAX_REQ_THPV_LEVS  ! Max no. of levels for pvort output
        PARAMETER (MAX_REQ_THPV_LEVS = MAX_P_LEVELS )
      INTEGER  MAX_ADJ_TSL      ! Max A_ADJSTEPS
        PARAMETER (MAX_ADJ_TSL  = 10  )
      INTEGER  MAX_N_INTF_A     ! Max no. of atmos interface areas
        PARAMETER (MAX_N_INTF_A =  8  )
      INTEGER  MAX_INTF_LEVELS  ! Max no. of atmos interface levels
        PARAMETER (MAX_INTF_LEVELS = MAX_P_LEVELS )
      INTEGER  MAX_NO_OF_SEGS   ! Maximum number of physics segments
        PARAMETER (MAX_NO_OF_SEGS = 200  )
C     MAX_N_INTF/MAX_N_INTF_A to be sorted out in next version
      INTEGER  MAX_N_INTF     ! Max no. of interface areas
        PARAMETER (MAX_N_INTF =  8  )


C*L------------------ COMDECK CINTFA ----------------------------------
CL CMAXSIZE should be called first.
C
C    Contains Variables, Headers and Index blocks for control of
C    generation of boundary information for the limited area model.
C
C    Interfaces to all other models are handled by STASH, and there is
C    no explicit coding written for them in the model.
C
C Interface variables initialised through INTFCNSTA
C namelist read in the interface control routine INTF_CTL.
CL
CL 29/07/98  CINTF comdeck renamed to CINTFA. New arrays LBC_STREAM_A
CL           and LBC_UNIT_NO_A added. INTF_AK/BK/AKH/BKH removed - now
CL           in ARGINFA/TYPINFA. D. Robinson.
CL
      INTEGER
     &  INTF_ROW_LENGTH  ! Interface field row length
     & ,INTF_P_ROWS      ! Interface field no of rows
     & ,INTF_P_LEVELS    ! Interface field no of levels
     & ,INTF_Q_LEVELS    ! Interface field no of wet levels
     & ,INTF_TR_LEVELS   ! Interface field no of tracer levels
     & ,INTFWIDTHA       ! Width of interface zone (atmosphere)
     & ,A_INTF_START_HR  ! ) Start, Frequency and End time in
     & ,A_INTF_FREQ_HR   ! ) hours for which atmosphere interface
     & ,A_INTF_END_HR    ! ) data is to be generated.
     & ,LEN_INTFA_P      ! Length of interface p field
     & ,LEN_INTFA_U      ! Length of interface u field
     & ,LEN_INTFA_DATA   ! Length of interface data
     & ,INTF_PACK        ! Packing Indicator for boundary data
     & ,LBC_STREAM_A     ! Output streams in UMUI
     & ,LBC_UNIT_NO_A    ! Unit Nos for Atmos Boundary Dataset
!
! Following 3 variables not in common ; in namelist
     & ,INTF_METH_LEV_CALC(MAX_N_INTF_A)
!                              !Method of calculating Eta level (ETAK)
!                              !from layers (ETAH)
     & ,INTF_MAX_SIG_HLEV(MAX_N_INTF_A)
!                              !level below which sigma coordinates used
     & ,INTF_MIN_PRS_HLEV(MAX_N_INTF_A)
!                              !level above which pressure coordinates

      REAL
     *  INTF_EWSPACE     ! E-W grid spacing (degrees)
     * ,INTF_NSSPACE     ! N-S grid spacing (degrees)
     * ,INTF_FIRSTLAT    ! Latitude of first row (degrees)
     * ,INTF_FIRSTLONG   ! Longitude of first row (degrees)
     * ,INTF_POLELAT     ! Real latitude of coordinate pole (degrees)
     * ,INTF_POLELONG    ! Real longitude of coordinate pole (degrees)

! Following variable not in common ; in namelist
      REAL INTF_ETAH(MAX_INTF_LEVELS+1,MAX_N_INTF_A)
C                           !Eta values at model layer boundaries ETAKH

      LOGICAL
     +  INTF_VERT_INTERP ! Switch to request vertical interpolation
     + ,LNEWBND          ! True for initialising new boundary data file

C*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/
     &  INTF_EWSPACE(MAX_N_INTF_A)    ,INTF_NSSPACE(MAX_N_INTF_A)
     & ,INTF_FIRSTLAT(MAX_N_INTF_A)   ,INTF_FIRSTLONG(MAX_N_INTF_A)
     & ,INTF_POLELAT(MAX_N_INTF_A)    ,INTF_POLELONG(MAX_N_INTF_A)
     & ,INTF_ROW_LENGTH(MAX_N_INTF_A) ,INTF_P_ROWS(MAX_N_INTF_A)
     & ,INTF_P_LEVELS(MAX_N_INTF_A)   ,INTF_Q_LEVELS(MAX_N_INTF_A)
     & ,INTF_TR_LEVELS(MAX_N_INTF_A)  ,INTFWIDTHA(MAX_N_INTF_A)
     & ,A_INTF_START_HR(MAX_N_INTF_A) ,A_INTF_FREQ_HR(MAX_N_INTF_A)
     & ,A_INTF_END_HR(MAX_N_INTF_A)   ,LEN_INTFA_P(MAX_N_INTF_A)
     & ,LEN_INTFA_U(MAX_N_INTF_A)     ,LEN_INTFA_DATA(MAX_N_INTF_A)
     & ,LNEWBND(MAX_N_INTF_A)         ,INTF_VERT_INTERP(MAX_N_INTF_A)
     & ,INTF_PACK(MAX_N_INTF_A)       ,LBC_STREAM_A(MAX_N_INTF_A)
     & ,LBC_UNIT_NO_A(MAX_N_INTF_A)
C---------------------------------------------------------------------

C  Namelist for atmos interface constants
!+ COMDECK CNAMINFA
!
!    Description:
!       This COMDECK contains the INTFCNSTA namelist which
!       contains all the variables required to define the grids
!       of the interface areas for Atmosphere Boundary data.
!
!       All variables are set up in the UMUI.
!       All variables are declared in comdeck CINTFA
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.5    03/08/98   New COMDECK created.
!
      NAMELIST/INTFCNSTA/
     &         INTF_EWSPACE,INTF_NSSPACE,INTF_FIRSTLAT,INTF_FIRSTLONG,
     &         INTF_POLELAT,INTF_POLELONG,INTF_VERT_INTERP,
     &         INTF_METH_LEV_CALC,INTF_ETAH,
     &         INTF_MAX_SIG_HLEV,INTF_MIN_PRS_HLEV,
     &         INTF_ROW_LENGTH,INTF_P_ROWS,INTF_P_LEVELS,INTF_Q_LEVELS,
     &         INTF_TR_LEVELS,INTFWIDTHA,INTF_PACK,
     &         A_INTF_FREQ_HR,A_INTF_START_HR,A_INTF_END_HR
     &        ,LBC_STREAM_A

!- End of comdeck CNAMINFA

!     Local variables
      INTEGER JINTF  !  loop index

!     Read in INTFCNSTA namelist to get output grids for
!     generating boundary data.

      REWIND 5
      READ (5,INTFCNSTA)
      REWIND 5

      IF (N_INTF_A.GT.0) THEN

!       Boundary data to be generated in this run.

        TOT_LEN_INTFA_P=0
        TOT_LEN_INTFA_U=0
        MAX_INTF_P_LEVELS=0

        DO JINTF=1,N_INTF_A

!         Calculate lengths for interface area JINTF

          LEN_INTFA_P(JINTF) = ( INTF_ROW_LENGTH(JINTF) +
     &    INTF_P_ROWS(JINTF) - 2*INTFWIDTHA(JINTF) )
     &    * 2 * INTFWIDTHA(JINTF)
          LEN_INTFA_U(JINTF) = LEN_INTFA_P(JINTF) - 4*INTFWIDTHA(JINTF)

!         Add on to total length

          TOT_LEN_INTFA_P = TOT_LEN_INTFA_P + LEN_INTFA_P(JINTF)
          TOT_LEN_INTFA_U = TOT_LEN_INTFA_U + LEN_INTFA_U(JINTF)
          MAX_INTF_P_LEVELS =
     &    MAX ( MAX_INTF_P_LEVELS , INTF_P_LEVELS(JINTF) )

        ENDDO

!       U_FIELD_INTFA Dimensions COEFF3 & COEFF4 in TYPINFA

        U_FIELD_INTFA = U_FIELD

        write (6,*) ' '
        write (6,*) ' Data lengths calculated in DERV_INTF_A.'
        do jintf=1,n_intf_a
        write (6,*) ' Area no ',jintf,
     &              ' len_intfa_p ',len_intfa_p(jintf),
     &              ' len_intfa_u ',len_intfa_u(jintf)
        enddo
        write (6,*) ' n_intf_a ',n_intf_a
        write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
        write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
        write (6,*) ' max_intf_p_levels ',max_intf_p_levels
        write (6,*) ' u_field_intfa ',u_field_intfa

      ELSE

!       No boundary conditions to be generated.
!       Initialise to prevent zero length dynamic allocation.

        write (6,*) ' n_intf_a ',n_intf_a

        N_INTF_A = 1
        TOT_LEN_INTFA_P = 1
        TOT_LEN_INTFA_U = 1
        MAX_INTF_P_LEVELS = 1
        U_FIELD_INTFA = 1

      write (6,*) ' n_intf_a ',n_intf_a
      write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
      write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
      write (6,*) ' max_intf_p_levels ',max_intf_p_levels
      write (6,*) ' u_field_intfa ',u_field_intfa

      ENDIF

      RETURN
      END
