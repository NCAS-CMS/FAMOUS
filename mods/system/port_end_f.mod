*ID GQZ1F505
*/--------------------------------------------------------
*/
*/ Backport vn5.5 mod gqz1505 and vn6.0 mods gqz3600, gel6600 
*/ to vn4.5. Original author frqz (Paul Dando)
*/
*/ Jeff Cole 18/11/03
*/
*/ Modset created by frqz on 26/02/03
*/
*/ Portable I/O changes.  Allows for input and output
*/ of UM start dumps, LBCs, Fields Files, PP files etc
*/ in big_endian format on LITTLE_ENDIAN platforms.
*/--------------------------------------------------------
*/
*DECLARE RDMULT1A
*B RDMULT1A.28
!    5.5    05/02/03  Portability changes allowing for big_endian
!                     I/O on little_endian platforms.        P.Dando
*D GPB4F403.70,GPB4F403.71
*D GBC5F404.225,GPB4F403.76
! LBC data is packed using CRAY 32 bit method - note that we need to
! read in ISIZE 32 bit words using BUFFIN32
          CALL BUFFIN32_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
*D GPB4F403.78
          IF (MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*ISIZE 32 bit words using BUFFIN32
            CALL BUFFIN32_SINGLE(NFT,buf,2*ISIZE,LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against ISIZE
            LEN_IO = LEN_IO/2
          ELSE
! For non-packed data
            CALL BUFFIN_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
          ENDIF
*/----------------------------------------------------------------------
*DECLARE READDM1A
*D READDM1A.167
        IF (MOD((LOOKUP(LBPACK,K)),10) .EQ. 2) THEN
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*IPTS 32 bit words using BUFFIN32
          CALL BUFFIN32(NFTIN,D1(LOOKUP(NADDR,K)),2*IPTS,LEN_IO,A)
! And then halve LEN_IO to satisfy tests against IPTS
          LEN_IO = LEN_IO/2
        ELSE
! For non-packed data
          CALL BUFFIN(NFTIN,D1(LOOKUP(NADDR,K)),IPTS,LEN_IO,A)
        ENDIF
*/----------------------------------------------------------------------
*DECLARE READFL1A
*I UBC1F402.27
         else if(pack_type.eq.2) then
           call buffin32(nftin,buf(1),2*ipts,len_io,a_io)
           len_io = len_io/2
*D READFL1A.99
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*IPTS 32 bit words using BUFFIN32
          CALL BUFFIN32(NFTIN,BUF(1),2*IPTS,LEN_IO,A_IO)
          LEN_IO = LEN_IO/2
*/----------------------------------------------------------------------
*DECLARE WTMULT1A
*B WTMULT1A.29
!     5.5  05/02/03  Portability changes allowing for big_endian
!                    I/O on little_endian platforms.        P.Dando
*D GPB0F401.957
        IF(MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
          IF(LOOKUP(DATA_TYPE) .EQ. 1) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*ISIZE 32 bit words using BUFFO32
            CALL buffo32_single(NFT,COMPBUF,2*ISIZE,LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against ISIZE
            LEN_IO = LEN_IO/2
          ENDIF
        ELSE
! For non-packed data
          CALL buffout_single(NFT,COMPBUF,ISIZE,LEN_IO,IOSTAT)
        ENDIF
*/----------------------------------------------------------------------
*DECLARE INTFOUT1
*B INTFOUT1.30
!    5.5  05/02/03    Portability changes allowing for big_endian
!                     I/O on little_endian platforms.        P.Dando
*D INTFOUT1.417

        IF (MOD(LOOKUP_INTF(LBPACK,1,JINTF),10) .EQ. 2) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*LEN_DATA 32 bit words using BUFFO32
          CALL BUFFO32(NFTOUT,INTF_DATA(1),2*LEN_DATA,LEN_IO,A_IO)
! And then halve LEN_IO to satisfy tests against LEN_DATA
          LEN_IO = LEN_IO/2
        ELSE
! For non-packed data
          CALL BUFFOUT(NFTOUT,INTF_DATA(1),LEN_DATA,LEN_IO,A_IO)
        ENDIF

*/----------------------------------------------------------------------
*DECLARE GENINTF1
*D GDR1F400.191
        IF (MOD(LOOKUP_INTFA(LBPACK,1,JINTF),10) .EQ. 2) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*LEN_DATA 32 bit words using BUFFO32
          CALL BUFFO32(NFTOUT,INTF_DATA(1),2*LEN_DATA,LEN_IO,A_IO)
! And then halve LEN_IO to satisfy tests against LEN_DATA
          LEN_IO = LEN_IO/2
        ELSE
! For non-packed data
          CALL BUFFOUT(NFTOUT,INTF_DATA(1),LEN_DATA,LEN_IO,A_IO)
        ENDIF
*/----------------------------------------------------------------------
*DECLARE WRITFL1A
*D GPB0F405.99
        if (mod(lookup(lbpack,k),10) .eq. 2) then
! Data is packed using cray 32 bit method - note that we need to write
! out 2*um_sector_ipts 32 bit words using buffo32
          call buffo32(nftout,buf(1),2*um_sector_ipts,len_io,a_io)
! And then halve len_io to satisfy tests against ipts_write
          len_io = len_io/2
        else
! For non-packed data
          call buffout(nftout,buf(1),um_sector_ipts,len_io,a_io)
        endif
*/----------------------------------------------------------------------
*DECK PBFIN32
*IF DEF,C96_1A,OR,DEF,C96_1B
*IF DEF,MPP
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2002, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
C
!+ Parallel UM version of BUFFIN32
!
! Subroutine Interface:
      SUBROUTINE BUFFIN32(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFIN32 routine for the Parallel Unified
!  Model. It is used where the same data which is 32-bit packed has
!  to be read in to all processors. (Contrast with READ_MULTI where
!  each processor reads its own local data from a larger global
!  field.)
!
! Method:
!  The C BUFFIN32 is renamed BUFFIN32_SINGLE under *DEF,MPP. This
!  routine causes BUFFIN32_SINGLE to be called by PE 0 only, and
!  then the data is broadcast to all other processors.
!
! Current Code Owner: Paul Dando
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    5.5    25/03/03 Code copied and modified from BUFFIN during
!                    development of portable I/O changes allowing
!                    big-endian data to be input and output on
!                    little-endian platforms.          P.Dando
!    6.0    13/10/03 Declare IOSTAT as REAL rather than INTEGER.
!                                                      P.Dando
!
! Code Description:
!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine Arguments with INTENT IN:

      INTEGER, INTENT(IN)  :: NFT          ! FORTRAN unit number
      INTEGER, INTENT(IN)  :: ISIZE        ! no. of words to write out

      REAL,    INTENT(IN)  :: ARRAY(ISIZE) ! Array to write out

! Subroutine Arguments with INTENT OUT:

      INTEGER, INTENT(OUT) :: LEN_IO       ! no. of words written out
      REAL,    INTENT(OUT) :: IOSTAT       ! Return code

! Parameters and Common blocks

*CALL PARVARS

! Local variables

      INTEGER :: info                      ! Broadcast information code
      REAL    :: return_codes(2)           ! Array to hold BUFFO32 code


! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE


      IF (mype .EQ. 0) THEN
        CALL BUFFIN32_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        return_codes(1)=LEN_IO
        return_codes(2)=IOSTAT
      ENDIF

c--get the broadcast flag
      call find_unit_bcast_flag(nft, info)
c--skip the broadcasts if the flag is set
      IF (info .EQ. 0) THEN
        CALL GC_RBCAST(1,ISIZE,0,nproc,info,ARRAY)    ! data
        CALL GC_RBCAST(2,2,0,nproc,info,return_codes) ! return codes
      ENDIF

      LEN_IO=return_codes(1)
      IOSTAT=return_codes(2)

      RETURN
      END

*ENDIF
*ENDIF
*/----------------------------------------------------------------------
*DECK PBFOUT32
*IF DEF,C96_1A,OR,DEF,C96_1B
*IF DEF,MPP
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2002, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
C
!+ Parallel UM version of BUFFO32
!
! Subroutine Interface:
      SUBROUTINE BUFFO32(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)

      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFO32 routine for the Parallel Unified
!  Model. It is used where the same data which is 32-bit packed has
!  to be written out from all processors. (Contrast with WRITE_MULTI
!  where each processor reads its own local data from a larger global
!  field.)
!
! Method:
!  The C BUFFO32 is renamed BUFFO32_SINGLE under *DEF,MPP. This
!  routine causes BUFFO32_SINGLE to be called by PE 0 only.
!  Note : No check is made that all processors are attempting
!         to write out identical data. It is assumed that the
!         data on PE 0 is the same as that on all other processors.
!
! Current Code Owner: Paul Dando
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    5.5    25/03/03 Code copied and modified from BUFFOUT as
!                    part of the portable I/O changes allowing
!                    big-endian data to be input and output on
!                    little-endian platforms.          P.Dando
!    6.0    13/10/03 Declare IOSTAT as REAL rather than INTEGER.
!                                                      P.Dando
!
! Code Description:
!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine Arguments with INTENT IN:

      INTEGER, INTENT(IN)  :: NFT          ! FORTRAN unit number
      INTEGER, INTENT(IN)  :: ISIZE        ! no. of words to write out

      REAL,    INTENT(IN)  :: ARRAY(ISIZE) ! Array to write out

! Subroutine Arguments with INTENT OUT:

      INTEGER, INTENT(OUT) :: LEN_IO       ! no. of words written out
      REAL,    INTENT(OUT) :: IOSTAT       ! Return code

! Parameters and Common blocks

*CALL PARVARS

! Local variables

      INTEGER :: info                      ! Broadcast information code
      REAL    :: stats(2)                  ! Array to hold BUFFO32 code

! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE

      IF (mype .EQ. 0) THEN
        CALL BUFFO32_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        stats(1)=LEN_IO
        stats(2)=IOSTAT
      ENDIF

      CALL GC_RBCAST(521,2,0,nproc,info,stats)
      LEN_IO=stats(1)
      IOSTAT=stats(2)

      RETURN
      END

*ENDIF
*ENDIF
