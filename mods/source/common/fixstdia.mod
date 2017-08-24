*ID FIXSTDIA
*/
*/  Change the way that deck ST_DIA11 does the selection of levels
*/  for diagnostics.  This method will choose the *closest* of the
*/  available levels.
*/
*/  The code being replaced insists on an *exact* match; this involves
*/  a floating point equality test -- and if it evaluates false, the
*/  level index is uninitialised and is later used as an array offset,
*/  segfaulting the program.
*/
*DECLARE ST_DIA11      
*D ST_DIA11.338,ST_DIA11.345
      Call level_indices(STASH_LEVELS(1,NI),
     $     T_P_LEVS,T_PRESS,T2_IND)
*D ST_DIA11.366,ST_DIA11.373
      Call level_indices(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,U2_IND)
*D ST_DIA11.394,ST_DIA11.401
      Call level_indices(STASH_LEVELS(1,NI),
     $     VCOMP_P_LEVS,VCOMP_PRESS,V2_IND)
*D ST_DIA11.226,ST_DIA11.238
      Call level_indices2(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     UV_IND)
*D ST_DIA11.272,ST_DIA11.284
      Call level_indices2(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,T_P_LEVS,T_PRESS,
     $     UT_IND)
*D ST_DIA11.305,ST_DIA11.317
      Call level_indices2(STASH_LEVELS(1,NI),
     $     VCOMP_P_LEVS,VCOMP_PRESS,T_P_LEVS,T_PRESS,
     $     VT_IND)
*D ST_DIA11.435,ST_DIA11.447
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,T_P_LEVS,T_PRESS,
     $     WT_IND)
*D ST_DIA11.468,ST_DIA11.480
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     WU_IND)
*D ST_DIA11.501,ST_DIA11.513
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     WV_IND)
*D ST_DIA11.548,ST_DIA11.560
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     QU_IND)
*D ST_DIA11.581,ST_DIA11.593
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     QV_IND)
*D ARS1F404.210,ARS1F404.222
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,W_P_LEVS,W_PRESS,
     $     QW_IND)
*D ARS1F404.264,ARS1F404.276
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Z_P_LEVS,Z_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     UZ_IND)
*D ARS1F404.293,ARS1F404.305
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Z_P_LEVS,Z_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     VZ_IND)
*B ST_DIA11.736
!
      CONTAINS
!
!---------------------------------------
      Subroutine level_indices2
     $     (stlev, npres1, press1, npres2, press2, levels)
!
! wrapper for level_indices where indices from two sets of levels 
! are to be extracted on the stash levels
!      
      implicit none
      integer,intent(in)::stlev(*),npres1,npres2
      real,intent(in)::press1(npres1),press2(npres2)
      integer,intent(out)::levels(stlev(1)*2)
!
      Call level_indices(stlev,npres1,press1,levels(1))
      Call level_indices(stlev,npres2,press2,levels(1+stlev(1)))
      End subroutine level_indices2
!     
!---------------------------------------
      Subroutine level_indices(stlev, npres, press, levels)
      implicit none
      integer,intent(in)::npres
      real,intent(in)::press(npres) ! pressures we have; sorted descending
!
      integer,intent(in)::stlev(*)
             ! pressures we want; sorted descending
             !  STLEV(1) is number of pressures
             !  STLEV(2:STLEV(1)+1) is integers representing pressures
!
      integer,intent(out)::levels(stlev(1))
          ! returns indices of pressures in PRESS array closest to 
          ! those specified by STLEV
!---
      integer::i,index
      real::required
      real,parameter::multiplier=2e-3
          ! factor of 1e-3 is to convert STLEV integers into pascals
          ! factor of 2 is so we can compare directly with sum of 2
          !     pressures rather than explicitly computing the mean.
!
      index=1
      do i=1,stlev(1)
         required=stlev(i+1)*multiplier
         findmatch: do while (index < npres)
           if (required > press(index)+press(index+1)) exit findmatch
           index=index+1
         enddo findmatch
         levels(i)=index
      enddo
!
!      PRINT *,'DEBUG level_indices...'
!      PRINT *,'INPUTS:'
!      PRINT *,npres,press
!      PRINT *,stlev(1:1+stlev(1))
!      PRINT *,'OUTPUT:'
!      PRINT *,levels
!      PRINT *,'============'
!      CALL FLUSH(6)
!      
      End subroutine level_indices
