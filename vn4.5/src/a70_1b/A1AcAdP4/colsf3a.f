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
!+ Subroutine to assemble lists of points with the same surfaces.
!
! Method:
!       The surfaces at the bottom of each profile are examined
!       and lists of points with the same type of surface are made.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENEQ replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE COLLECT_SURFACE(N_PROFILE
     &   , I_SURFACE
     &   , N_POINT_TYPE, INDEX_SURFACE
     &   , NPD_PROFILE, NPD_SURFACE
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_SURFACE
!             MAXIMUM NUMBER OF SURFACES
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_SURFACE(NPD_PROFILE)
!             SURFACE SPECIFICATIONS
      INTEGER   !, INTENT(OUT)
     &     N_POINT_TYPE(NPD_SURFACE)
!             NUMBER OF POINTS OF EEACH TYPE
     &   , INDEX_SURFACE(NPD_PROFILE, NPD_SURFACE)
!             LIST OF POINTS OF EACH TYPE
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , J
!
!
!     PASS THROUGH ALL THE COLUMNS COLLECTING LISTS OF POINTS WHICH
!     HAVE THE SAME SURFACE TYPE.
      DO K=1, NPD_SURFACE
         N_POINT_TYPE(K)=0
      ENDDO
C
C
      DO   K=1,NPD_SURFACE        
           J=1
        DO L=1,N_PROFILE
          IF (I_SURFACE(L).EQ.K) THEN
            INDEX_SURFACE(J,K)=L
            J=J+1
            N_POINT_TYPE(K)=N_POINT_TYPE(K)+1
          END IF
        END DO
      ENDDO
!
!
!
      RETURN
      END
