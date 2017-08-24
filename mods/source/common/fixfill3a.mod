*ID FIXFILL3A
*/
*/  fix a multiple initialisation in deck FILL3A
*/
*/  Having initialised the whole array to .FALSE. we are not
*/  allowed to initialise individual elements to .TRUE. using
*/  DATA statements.
*/
*/  Most compilers treat this as a warning, but the Intel compiler
*/  treats it as an error.
*/
*/  Just set the array elements using assignments instead...
*/
*/  Alan Iwi 26/11/01
*/
*DECLARE FILL3A
*D ADB2F404.235,239
*I ADB2F404.247
      L_IN_CLIMAT(IP_WATER_SOLUBLE) = .TRUE.
      L_IN_CLIMAT(IP_DUST_LIKE) = .TRUE.
      L_IN_CLIMAT(IP_OCEANIC) = .TRUE.
      L_IN_CLIMAT(IP_SOOT) = .TRUE.
      L_IN_CLIMAT(IP_SULPHURIC) = .TRUE.
