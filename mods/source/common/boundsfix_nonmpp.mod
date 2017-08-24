*IDENT BOUNDSNONMPP
*/
*/ Fixes necessary for bounds checking non-MPP code. 
*/
*/ These fixes are all from Paul Valdes's coupled_bugs4.
*/ 	Annette Osprey 11-Sep-2006
*/
*/
*/ First problem was that someone forgot to caveat last call 
*/ for when not running in MPP. (PJV)
*/
*DECLARE CNVSTOP
*B CNVSTOP.173
*IF DEF,MPP
*I CNVSTOP.174
*ENDIF
*/
*/
*/   This is a bounds problem assoicated with the calculation of 
*/   mountain torque. It only effects the non-MPP version. Fairly
*/   confident about the solution. (PJV)
*/
*DECLARE DYNDIA1A
*D GSM1F405.732
        DO i=FIRST_FLD_PT,LAST_U_FLD_PT-1
