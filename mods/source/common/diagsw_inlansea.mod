*/
*/ if L_INLANSEA is used to cap SSS values to be between
*/ SALLOW and SALUP, the diagnostic DIAGSW (2,30,280) is 
*/ modified to reflect the virtual salinity flux that would
*/ be required to make this change. DIAGSW has the same sense
*/ as P-E however, such that a negative change in SSS (fresher)
*/ requires a positive change in DIAGSW (more water). This
*/ appears to have been originally coded with the opposite sign
*/
*/ r.s.smith@reading.ac.uk 24/09/15
*/ 
*DECLARE TRACER
*D OJL1F405.31
     &          -1e9*dz(1)*(sallow-TA(I,1,2))/(C2DTTS*100.0)
*D OJL1F405.36
     &          -1e9*dz(1)*(salup-TA(I,1,2))/(C2DTTS*100.0)
