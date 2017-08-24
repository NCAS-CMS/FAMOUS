*ID MODELFIX
*/
*/ Delete barrier call in TIMER3A routine, can get deadlock otherwise.
*/
*DECLARE TIMER3A
*D GPB8F405.15
C   Skipping sync
C        CALL GC_GSYNC(nproc,info)
*/
*/ Initialise O_ATMCO2 array correctly
*/
*DECLARE SWAPA2O2
*D CCN1F405.160
*I GRR0F403.229
      DO I=1,CO2_DIMA
      O_ATMCO2(I) = RMDI
      ENDDO ! I
*/
*/ This code will work on any MPP machine, therefore change 
*/ *IF  DEF,MPP,AND,DEF,T3E to *IF  DEF,MPP.
*/
*DECLARE SWPLND1A
*D SWPLND1A.3
*IF  DEF,MPP
*/
*/ If running coupled model then J will be declared twice
*/
*DECLARE INTFCTL1
*I GMB1F405.146
      INTEGER J                   ! Ocean: For setting LBC_UNIT_NO_O
*D UDG1F305.178
*D GMB1F405.154
*/
*/ If using VECTLIB delete non-fortran line
*/
*DECLARE TRSFC3A
*DELETE PXVECTLB.149
*/
*/ On some machines (e.g. hp,ibm) abort doesn't take any arguments
*/ this causes a compilation error. If you declare abort as external
*/ compiler doesn't give an error but abort is still called, without
*/ writing the error message.
*/
*DECLARE FIELDOP1
*D FIELDOP1.105
      External readff,setpos,ioerror,fieldop_main,abort
*/
*DECLARE FIELDCOS
*D FIELDCOS.35
      EXTERNAL READFF,SETPOS,IOERROR,READ_WRITE,ABORT
*/
*DECLARE DUMMYVEG
*I DUMMYVEG.30

      EXTERNAL ABORT
*/
*DECLARE RDLSM1A
*B RDLSM1A.60

      EXTERNAL ABORT
*/
*DECLARE SOLANG1A
*B SOLANG1A.71

      EXTERNAL ABORT
