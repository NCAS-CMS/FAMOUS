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

!+  subroutine to convert fields from full grid to blocks
!
! Subroutine Interface:
      subroutine blokfld(field,blkfld,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

      IMPLICIT NONE
!
! Description:
! This subroutine converts the input array of data on the full grid
! (mdata points) to blocked data for use by WAM : fldout(niblo,nblo)
!
! Method:
!   The WAM model grid and block indexing arrays are used
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! <1.0>    July 1995 Original code.    M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <UM wave model>
! System Task:              <UM wave model>
!
!- End of header

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


       real blkfld(niblo,nblo) ! out - field arranged by blocks
       integer icode           ! out - return code

       integer nstart ! index on full grid of first point in block
       integer nend   ! index on full grid of last  point in block
       integer mdata  ! array size - max number of datapoints

       integer ig       ! index to loop over block number
       integer ii,i     ! index looping over points in block
       integer ifill    ! counting index to fill blocks
       real field(mdata)       !in field of values at data points

C      initialise output array with zeros

       do ig=1,nblo
        do ii=1,niblo
         blkfld(ii,ig)=0.
        enddo
       enddo

       nstart=1
       do ig=1,igl
        nend  = nstart + ijlt(ig) - 1

        ifill=1
        do i=nstart,nend
         blkfld(ifill,ig)=field(i)
         ifill=ifill+1
        enddo
        nstart = nend - (ijlt(ig)-ijls(ig))

       enddo


       return
      end                                                              
