*ID ABORTMOD
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
