*IDENT GNUML45
*/
*/ Update to set gnu to FKAPB (i.e. to set vert tracer diffusivity 
*/ to its background value) everywhere.
*/
*/ Written for UM vn 3.1              RAW 09/11/94
*/ Converted for vn 3.2               RAW 3/3/95
*/ Should work OK for vn4.3           JRP 27/6/97
*/
*/ Re-written to replace all the executable code in VERTCOFT.
*/	Annette Osprey 15/09/06
*/ 
*/ This is a FAMOUS mod. 
*/ 
*/
*DECLARE VERTCOFT
*D VERTCOFT.60,VERTCOFT.157
C
      DO K=1,KM
        DO I=1,IMT
          gnu(I,K)=KAPPA_B_SI(K)*1.E4*FM(I,K) ! use background value
        END DO 
      END DO
C