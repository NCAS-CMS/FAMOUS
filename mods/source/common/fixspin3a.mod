*ID FIXSPIN3A
*/
*/  fix a multiple initialisation in deck SPIN3A
*/
*/  Most compilers treat this as a warning, but the Intel compiler
*/  treats it as an error.
*/
*/  Just set the array elements using assignments instead...
*/
*/  Lois Steenman-Clark 18.08.05
*/
*DECLARE SPIN3A
*D SPIN3A.1129
*I SPIN3A.1130
      N_SCALE_VARIABLE(IP_SCALE_POWER_LAW)=2
      N_SCALE_VARIABLE(IP_SCALE_FNC_NULL)=0
      N_SCALE_VARIABLE(IP_SCALE_POWER_QUAD)=3
      N_SCALE_VARIABLE(IP_SCALE_DOPPLER_QUAD)=4
*/
