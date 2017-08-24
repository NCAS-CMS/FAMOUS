C      Subroutine WAV_TO_D1
!+ Subroutine to convert blocked stress fields onto full grid and
!+ place into D1
! Description:
!   input arrays are blocked ustar and wave induced stresses for the
!   timestep just completed
!   these are converted onto the full grid and placied into D1
!
! Method:
!   use wind directions in thwnew to get components
!   call unblok to go to full grid
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! UM4.1    July 1996 Original code.  M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header
      Subroutine WAV_TO_D1(mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C

     & ustar_in,
     & windir_in,
     & tauw_in,

     & wvstressx_out,
     & wvstressy_out,
     & taux_out,
     & tauy_out,
     & icode,cmessage)

      implicit none

C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
       INTEGER
     & NANG,       ! number of direction components
     & NFRE,       ! number of frequency components
     & NGX,        ! number of cols in LS mask grid
     & NGY,        ! number of rows in LS mask grid
     & NBLO,       ! max number of blocks
     & NIBLO,      ! max number datapoints per block
     & NOVER,      ! max number datapoints in overlap row
     & NIBLD, NBLD, NIBLC, NBLC
C
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
      real DELPHI  ! grid increment for latitude (in metres)
      real DELLAM  ! grid increment for long. at equator (in metres)
      real SINPH(NGY), COSPH(NGY) ! sin / cos of latitude

      integer IGL   ! number of blocks
      integer IJS(NBLO)  ! index of first point of second lat
      integer IJL2(NBLO) ! index of last  point of second lat
      integer IJLS(NBLO) ! index of first point of lat before last
      integer IJL(NBLO)  ! index of last  point of lat before last
      integer IJLT(NBLO) ! total number of gridpoints in a block
C

C ARRAYS IN

      INTEGER mdata ! total number of data points
      INTEGER ICODE

      CHARACTER*80 CMESSAGE

      REAL
     &  ustar_in(niblo,nblo)   ! blocked u* values (usnew)
     & ,windir_in(niblo,nblo)  ! blocked dir values (thwnew)
     & ,tauw_in(niblo,nblo)    ! blocked wv stress values (tauw)


C ARRAYS OUT

       REAL
     &  wvstressx_out(mdata) ! x cmpt tauw  D1(jwv_stressx)
     & ,wvstressy_out(mdata) ! y cmpt tauw  D1(jwv_stressy)
     & ,taux_out(mdata)      ! x cmpt ustar D1(jwv_taux)
     & ,tauy_out(mdata)      ! y cmpt ustar D1(jwv_tauy)

C LOCAL ARRAYS

      INTEGER ig,ii         ! loop counters

      REAL
     &   ustarx(niblo,nblo) ! to hold x comp ustar before unblocking
     &  ,ustary(niblo,nblo) !         y            back into D1
     &  ,tauwx(niblo,nblo)  !      "         tauw
     &  ,tauwy(niblo,nblo)  !      "


C first get arrays of components

       do ig=1,nblo
         do ii=1,niblo
            ustarx(ii,ig)=0.
            ustary(ii,ig)=0.

            tauwy(ii,ig)=0.
            tauwx(ii,ig)=0.
         enddo
      enddo

      do ig=1,igl

        do ii=ijs(ig),ijl(ig)

C note WAM direction convention has zero north increasing clockwise

          ustarx(ii,ig)=ustar_in(ii,ig)*sin(windir_in(ii,ig))
          ustary(ii,ig)=ustar_in(ii,ig)*cos(windir_in(ii,ig))

          tauwx(ii,ig)=tauw_in(ii,ig)*sin(windir_in(ii,ig))
          tauwy(ii,ig)=tauw_in(ii,ig)*cos(windir_in(ii,ig))

        enddo
       enddo

       CALL unblok(wvstressx_out
     & ,tauwx
     & ,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     &  icode)

       CALL unblok(wvstressy_out
     & ,tauwy
     & ,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     &  icode)

       CALL unblok(taux_out
     & ,ustarx
     & ,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     &  icode)

       CALL unblok(tauy_out
     & ,ustary
     & ,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     &  icode)

      return
      end
