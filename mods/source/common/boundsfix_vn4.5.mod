*IDENT BOUNDSFIX
*/
*/ General fixes necessary for bounds checking vn4.5 code.
*/ Tested on HadCM3 and FAMOUS
*/ 
*/ These fixes are all from Paul Valdes's coupled_bugs4.
*/ 	Annette Osprey 11-Sep-2006
*/
*/
*/   For some unknown reason, with some combinations of stash/input 
*/   dumps, the reconfiguration gives a variable which has undefined
*/   (i.e. 0) stashcode. This causes a bounds error. Not too sure
*/   about solution since the new if statements means that no space
*/   is allocated for these mysterious stash variables which may be
*/   OK but ...... (PJV)
*/
*DECLARE INANCA1A
*I GRB4F305.223
        if (stashancil(i).gt.0) then
*I GRB4F305.232
        end if
*I INANCA1A.655
        if (stashancil(i).gt.0) then
*I GDR8F400.95
        end if
*I ADR1F304.115
        if (fileancil(i).gt.0) then
*I INANCA1A.701
        end if
*/
*/
*/   This is a real bounds error problem but don't fully
*/   understand. In PP files, level dependent constants are
*/   defined with a second dimension 4 within subroutine
*/   initpp but the calling argument suggests that they 
*/   should be dimensioned 1 for ocean pp files. Simple
*/   correction stops bounds error. (PJV)
*/
*DECLARE INITPP1A
*D INITPP1A.98
      PP_FIXHD(112)=
     : MIN(LEN2_LEVDEPC,PP_LEN2_LEVDEPC)
*D INITPP1A.136
      DO 5 II=1,LEN1_LEVDEPC*
     :    MIN(LEN2_LEVDEPC,PP_LEN2_LEVDEPC)
*/
*/
*/   The FFT routine FOURIER passes array arguments to FTRANS in the 
*/   f77 style by the starting array element only. FTRANS accesses 
*/   beyond the dimensions of the dummy array (which causes bounds 
*/   error), but refers to element of the actual array in FOURIER 
*/   (presumably at the appropriate place). 
*/   So just define all as assumed-size arrays in FTRANS. (AO)
*/
*DECLARE FOURIE3A
*D PXORDER.17,21
      REAL A(*),     ! First real input vector
     &     B(*),
     &     C(*),     ! First real output vector
     &     D(*),
     &     TRIGS(*)  ! Precalculated list of sines & cosines
