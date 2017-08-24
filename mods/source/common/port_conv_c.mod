*ID PORTCONVC
*/--------------------------------------------------------
*/
*/ Backport vn5.5 mod gqz3505 and vn6.0 mod gqz1600 to vn4.5
*/ Original author frqz (Paul Dando)
*/
*/ Jeff Cole 14/11/03
*/
*/ Provides portable data conversion routines to replace
*/ the Cray-specific routines on non-Cray platforms
*/--------------------------------------------------------
*/
*/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*DECLARE PORTIO2A
*/---------------------------------------------------------
*B GBC3C405.5
  /* 5.5  28/02/03  Add portable data conversion routines:      */
  /*                IEEE2IBM, IBM2IEEE, IEEE2IEG, IEG2IEEE      */
  /*                MOVEBITS and MOVEBYTES.                     */
  /*                These replace the Cray-specific routines    */
  /*                CRI2IBM, IBM2CRI, CRI2IEG, IEG2CRI, MOVBIT  */
  /*                and STRMOV respectively.                    */
  /*                (Author: S.Carroll)               P.Dando   */
  /* 6.0  23/07/03  Additional functionality added to portable  */
  /*                data conversion routines.                   */
  /*                Some improvement to efficiency of MOVEBITS  */
  /*                and MOVEBYTES.                    P.Dando   */
*B PORTIO2A.16

/****************************************************************/
/* The following global constants are required by the portable  */
/* data conversion routines.  These constants provide masks for */
/* for selectively extracting the sign, exponent and mantissa   */
/* parts from IEEE and IBM 32 and 64 bit reals and integers.    */
/* The constants are expressed here in hex format.  The         */
/* corresponding binary representations are given by the usual  */
/* relationships:                                               */
/*    0=0000, 1=0001, 2=0010,....., 9=1001, A=1010,....F=1111   */
/****************************************************************/

/* Format 1: IBM 32 bit floating point numbers                  */
#define ibm_sign_32   0x80000000
#define ibm_expo_32   0x7F000000
#define ibm_mant_32   0x00FFFFFF
#define ibm_mant_bits_32 24
#define ibm_expo_bits_32 7
#define ibm_expo_bias_32 64

*IF DEF,FRL8,OR,DEF,CRAY
/* Format 2: IBM 64 bit floating point numbers                  */
#define ibm_sign_64   0x8000000000000000LL
#define ibm_expo_64   0x7F00000000000000LL
#define ibm_mant_64   0x00FFFFFFFFFFFFFFLL
#define ibm_mant_bits_64 56
#define ibm_expo_bits_64 7
#define ibm_expo_bias_64 64
*ENDIF

/* Format 3: IEEE 32 bit floating point numbers                 */
#define ieee_sign_32  0x80000000
#define ieee_expo_32  0x7F800000
#define ieee_mant_32  0x007FFFFF
#define ieee_mant_bits_32 23
#define ieee_expo_bits_32 8
#define ieee_expo_bias_32 127

*IF DEF,FRL8,OR,DEF,CRAY
/* Format 4: IEEE 64 bit floating point numbers                 */
#define ieee_sign_64  0x8000000000000000LL
#define ieee_expo_64  0x7FF0000000000000LL
#define ieee_mant_64  0x000FFFFFFFFFFFFFLL
#define ieee_mant_bits_64 52
#define ieee_expo_bits_64 11
#define ieee_expo_bias_64 1023
*ENDIF

/* Format 6: IEEE 16 bit integers (two's complement)            */
#define ieee_int_sign_16  0x8000
#define ieee_int_mant_16  0x7FFF
#define ieee_int_mant_bits_16 15

/* Format 7: IEEE 32 bit integers (two's complement)            */
#define ieee_int_sign_32  0x80000000
#define ieee_int_mant_32  0x7FFFFFFF
#define ieee_int_mant_bits_32 31

*IF DEF,FRL8,OR,DEF,CRAY
/* Format 8: IEEE 64 bit integers (two's complement)            */
#define ieee_int_sign_64  0x8000000000000000LL
#define ieee_int_mant_64  0x7FFFFFFFFFFFFFFFLL
#define ieee_int_mant_bits_64 63
*ENDIF

*B PORTIO2A.279

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine CRI2IEG -        */
/* interface to read_number/write_number                        */
/****************************************************************/
integer
*IF DEF,C_LOW
ieee2ieg
*ELSEIF DEF,C_LOW_U
ieee2ieg_
*ELSE
IEEE2IEG
*ENDIF
(type, num, ieg_num_out, offset_out, cri_num_in, stride,
    size_num_in, size_num_out)
integer *type, *num, *offset_out, *stride, *size_num_in, *size_num_out;
integer ieg_num_out[], cri_num_in[];
{
/* Prototype functions */
void read_number();
void write_number();
/* Local variables */
integer errorcode=0 ;   /* Return error code (success=0) */
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer); /* Length of integer in bits    */
int out_offset ;

/* Variables in calls to read/write_number */
integer sign=0, expo=0, mant=0, mant_bits_read=0;

/* Check for valid input */

if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid num = %d. Return code = %d\n",
     *num, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid stride = %d. Return code = %d\n",
     *stride, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Convert Cray offset to my offset */

 offset = word_length - (int) *offset_out - (int) *size_num_out;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%( (int) *size_num_out ) != 0
               || offset > word_length
               || (int) *offset_out < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid bitoff = %d. Return code = %d\n",
     *offset_out, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if ((int) *type == 2){

  /* INTEGER data types */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 3){

  /* REAL data types */

  switch ((int) *size_num_in){               /* real IEEE in */
    case 32:  type_num_in = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 4;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* real IEEE out */
    case 32:  type_num_out = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 4;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid forlen %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 5){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
          "IEEE2IEG: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else {
  /* ERROR - unsupported type */
  errorcode = -2 ;
  sprintf(message,
   "IEEE2IEG: Error - unsupported data type = %d. Return code = %d\n",
     *type, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Loop over all values converting them as we go along */

out_offset = offset ;

for (i=0; i < *num; i++){
  read_number(type_num_in, (int) *size_num_in, 0, *cri_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, out_offset,
               ieg_num_out, sign,  expo,  mant,  mant_bits_read);

  if ( out_offset == 0 ) {
    out_offset = word_length  - (int) *size_num_out ;
    ieg_num_out++ ;
  }
  else {
    out_offset = out_offset - (int) *size_num_out ;
  }

  cri_num_in++ ;

}
 return errorcode ;
}

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine IEG2CRI -        */
/* interface to read_number/write_number                        */
/****************************************************************/
integer
*IF DEF,C_LOW
ieg2ieee
*ELSEIF DEF,C_LOW_U
ieg2ieee_
*ELSE
IEG2IEEE
*ENDIF
(type, num, ieg_num_in, offset_in, cri_num_out, stride,
    size_num_out, size_num_in)
integer *type, *num, *offset_in,  *stride, *size_num_in, *size_num_out;
integer ieg_num_in[], cri_num_out[];
{
/* Prototype functions */
void read_number();
void write_number();
/* Local variables */
integer errorcode=0 ;  /* Return error code (success=0) */
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer); /* Length of integer in bits    */
int in_offset ;

/* Variables in calls to read/write_number */
integer sign=0, expo=0, mant=0, mant_bits_read=0;

/* Check for valid input */

if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid num = %d. Return code = %d\n",
     *num, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid stride = %d. Return code = %d\n",
     *stride, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Convert Cray offset to my offset */

 offset = word_length - (int) *offset_in - (int) *size_num_in;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%( (int) *size_num_in ) != 0
               || offset > word_length
               || (int) *offset_in < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid bitoff = %d. Return code = %d\n",
     *offset_in, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}


/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if ((int) *type == 2){

  /* INTEGER data types */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 3){

  /* REAL data types  */

  switch ((int) *size_num_in){               /* real IEEE in */
    case 32:  type_num_in = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 4;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* real IEEE out */
    case 32:  type_num_out = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 4;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
           "IEG2IEEE: Error - Invalid forlen %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 5){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
          "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else {

  /* ERROR - unsupported data type */

  errorcode = -2 ;
  sprintf(message,
   "IEG2IEEE: Error - unsupported data type = %d. Return code = %d\n",
     *type, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Loop over all values converting them as we go along */

in_offset = offset ;

for (i=0; i<*num; i++){
  read_number(type_num_in, (int) *size_num_in, in_offset, *ieg_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, 0, cri_num_out,
                      sign,  expo,  mant,  mant_bits_read);

  if ( in_offset == 0 ) {
    in_offset = word_length - (int) *size_num_in ;
    ieg_num_in++ ;
  }
  else {
    in_offset = in_offset - (int) *size_num_in ;
  }

  cri_num_out++ ;

}
 return errorcode ;
}

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine CRI2IBM -        */
/* interface to read_number/write_number                        */
/****************************************************************/
integer
*IF DEF,C_LOW
ieee2ibm
*ELSEIF DEF,C_LOW_U
ieee2ibm_
*ELSE
IEEE2IBM
*ENDIF
(type, num, ibm_num_out, offset_out, cri_num_in, stride,
    size_num_in, size_num_out)
integer *type, *num, *offset_out, *stride, *size_num_in, *size_num_out;
integer ibm_num_out[], cri_num_in[];
{
/* Prototype functions */
void read_number();
void write_number();
/* local variables */
integer errorcode=0 ;  /* Return error code (success=0) */
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer); /* Length of integer in bits    */
int out_offset ;

/* Variables in calls to read/write_number */
integer sign=0, expo=0, mant=0, mant_bits_read=0;

/* Check for valid input */

if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid num = %d. Return code = %d\n",
     *num, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid stride = %d. Return code = %d\n",
     *stride, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Convert Cray offset to my offset */
 offset = word_length - (int) *offset_out - (int) *size_num_out;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%( (int) *size_num_out ) != 0
               || offset > word_length
               || (int) *offset_out < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid bitoff = %d. Return code = %d\n",
     *offset_out, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Cray fortran uses IEEE data types */
/* decide which data types are used */

if ((int) *type == 2){

  /* INTEGER data types */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 3){

  /* REAL data types */

  switch ((int) *size_num_in){               /* real IEEE in */
    case 32:  type_num_in = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 4;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* real IBM out */
    case 32:  type_num_out = 1;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 2;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
           "IEEE2IBM: Error - Invalid forlen %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 5){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d . Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else {

  /* ERROR - unsupported type */

  errorcode = -2 ;
  sprintf(message,
   "IEEE2IBM: Error - unsupported data type = %d. Return code = %d\n",
     *type, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Loop over all values converting them as we go along */

out_offset = offset ;

for (i=0; i < *num; i++){

  read_number(type_num_in, (int) *size_num_in, 0, *cri_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, out_offset,
               ibm_num_out, sign,  expo,  mant,  mant_bits_read);

  if ( out_offset == 0 ) {
    out_offset = word_length - (int) *size_num_out ;
    ibm_num_out++ ;
  }
  else {
    out_offset = out_offset - (int) *size_num_out ;
  }

  cri_num_in++ ;

}
 return errorcode ;
}

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine IBM2CRI -        */
/* interface to read_number/write_number                        */
/****************************************************************/
integer
*IF DEF,C_LOW
ibm2ieee
*ELSEIF DEF,C_LOW_U
ibm2ieee_
*ELSE
IBM2IEEE
*ENDIF
(type, num, ibm_num_in, offset_in, cri_num_out, stride,
size_num_out, size_num_in)
integer *type, *num, *offset_in, *stride, *size_num_in, *size_num_out;
integer ibm_num_in[], cri_num_out[];
{
/* Prototype functions */
void read_number();
void write_number();
/* Local variables */
integer errorcode=0 ;   /* Return error code (success=0) */
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer); /* Length of integer in bits    */
int in_offset ;

/* Variables in calls to read/write_number */
integer sign=0, expo=0, mant=0, mant_bits_read=0;

/* Check for valid input */

if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid num = %d. Return code = %d\n",
     *num, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid stride = %d. Return code = %d\n",
     *stride, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Convert Cray offset to my offset */
 offset = word_length - (int) *offset_in - (int) *size_num_in;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%( (int) *size_num_in ) != 0
               || offset > word_length
               || (int) *offset_in < 0){
  errorcode = -4 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid bitoff = %d. Return code = %d\n",
     *offset_in, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if ((int) *type == 2){

  /* INTEGER data types */

  switch ((int) *size_num_in){                  /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){                 /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 3){

  /* REAL data types */

  switch ((int) *size_num_in){                  /* real IBM in */
    case 32:  type_num_in = 1;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 2;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){                 /* real IEEE out */
    case 32:  type_num_out = 3;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 4;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else if ((int) *type == 5){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch ((int) *size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = 6;  break;
    case 32:  type_num_in = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_in = 8;  break;
*ENDIF
    default:  errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           *size_num_in, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
  switch ((int) *size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = 6;  break;
    case 32:  type_num_out = 7;  break;
*IF DEF,FRL8,OR,DEF,CRAY
    case 64:  type_num_out = 8;  break;
*ENDIF
    default:  errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           *size_num_out, errorcode);
        CALL_MESSAGE_PRINT(message);
        return errorcode ;
  }
}
else {
  /* ERROR - unsupported type */
  errorcode = -2 ;
  sprintf(message,
   "IBM2IEEE: Error - unsupported data type = %d. Return code = %d\n",
     *type, errorcode);
  CALL_MESSAGE_PRINT(message);
  return errorcode ;
}

/* Loop over all values converting them as we go along */

in_offset = offset ;

for (i=0; i<*num; i++){
  read_number(type_num_in, (int) *size_num_in, in_offset, *ibm_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, 0, cri_num_out,
                      sign,  expo,  mant,  mant_bits_read);

  if ( in_offset == 0 ) {
    in_offset = word_length - (int) *size_num_in ;
    ibm_num_in++ ;
  }
  else {
    in_offset = in_offset - (int) *size_num_in ;
  }

  cri_num_out++ ;

}
 return errorcode ;
}

/****************************************************************/
/* The following function is to be called from C.  It           */
/* reads in numbers of different types:                         */
/* real/integer, IBM or IEEE  with offsets/packed etc           */
/* Called from the portable data conversion routines            */
/****************************************************************/
void read_number
(format,size,offset,in_number,out_sign,out_expo,out_mant,mant_bits_read)
int format, size, offset;
integer in_number;
integer *out_sign, *out_expo, *out_mant, *mant_bits_read;
{

integer expo_bits, expo_bias;
integer sign_mask, expo_mask, mant_mask, in_number_mask;
integer one ;
int i, first_bit;
int word_length=8*sizeof(integer) ; /* Length of word in bits */

switch (format){
  case 1:
    /* format 1 IBM 32 bits */
    sign_mask =      ibm_sign_32;
    expo_mask =      ibm_expo_32;
    expo_bits =      ibm_expo_bits_32;
    expo_bias =      ibm_expo_bias_32;
    mant_mask =      ibm_mant_32;
    *mant_bits_read = ibm_mant_bits_32;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 2:
    /* format 2 IBM 64 bits */
    sign_mask =      ibm_sign_64;
    expo_mask =      ibm_expo_64;
    expo_bits =      ibm_expo_bits_64;
    expo_bias =      ibm_expo_bias_64;
    mant_mask =      ibm_mant_64;
    *mant_bits_read = ibm_mant_bits_64;
  break;
*ENDIF

  case 3:
    /* format 3  IEEE 32 bits */
    sign_mask =      ieee_sign_32;
    expo_mask =      ieee_expo_32;
    expo_bits =      ieee_expo_bits_32;
    expo_bias =      ieee_expo_bias_32;
    mant_mask =      ieee_mant_32;
    *mant_bits_read = ieee_mant_bits_32;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 4:
    /* format 4  IEEE 64 bits */
    sign_mask =      ieee_sign_64;
    expo_mask =      ieee_expo_64;
    expo_bits =      ieee_expo_bits_64;
    expo_bias =      ieee_expo_bias_64;
    mant_mask =      ieee_mant_64;
    *mant_bits_read = ieee_mant_bits_64;
  break;
*ENDIF

  case 6:
    /* format 6 IEEE 16 bit integers */
    sign_mask =      ieee_int_sign_16;
    mant_mask =      ieee_int_mant_16;
    *mant_bits_read = ieee_int_mant_bits_16;
    expo_bits = 0;
    expo_bias = 0;
  break;

  case 7:
    /* format 7 IEEE 32 bit integers */
    sign_mask =      ieee_int_sign_32;
    mant_mask =      ieee_int_mant_32;
    *mant_bits_read = ieee_int_mant_bits_32;
    expo_bits = 0;
    expo_bias = 0;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 8:
    /* format 8 IEEE 64 bit integers */
    sign_mask =      ieee_int_sign_64;
    mant_mask =      ieee_int_mant_64;
    *mant_bits_read = ieee_int_mant_bits_64;
    expo_bits = 0;
    expo_bias = 0;
  break;
*ENDIF
}

if ( offset != 0 || size != (*mant_bits_read + expo_bits + 1)
     || (word_length - offset -size != 0) ){
  /* For offsets and size differences */

  /* shift number to right and clear off any other data in string */
  in_number = in_number >> offset;
  one = 1 ;
  in_number_mask = ( (one << size) - 1 );
  in_number = in_number & in_number_mask;

  /* sign_mask */
  one = 1 ;
  sign_mask = one << (size - 1);
}


/* For all types of nos */

/* Extract 1 bit sign and remove trailing 0s*/
  *out_sign= (in_number & sign_mask) >> (*mant_bits_read + expo_bits);

/* mantissa and expo */

if ((in_number & mant_mask) == 0 && (in_number & expo_mask) == 0) {
  *out_mant = 0 ;
  *out_expo = 0 ;
}
else{
if( format == 1 || format ==2 ){

  /* IBM type numbers */

  /*  find first non-zero bit
  (up to first three leading bits could be 0 from IBM)*/
  first_bit = 0;
  for (i = 1; first_bit == 0; i++){
    if(((in_number & mant_mask) >> (*mant_bits_read - i)) & 1  == 1){
      first_bit = i;
    }
    if(i == *mant_bits_read) first_bit = *mant_bits_read;
  }

  /* shift mantissa to be of the form 1.fraction
     and scale expo appropriately */
  *out_mant = (in_number & mant_mask) << first_bit;

  /*  IBM type - scale base 16 exponent to decimal and then correct
  for normalised mantissa */

  *out_expo =((((in_number & expo_mask) >> *mant_bits_read)
                - expo_bias) * 4) - first_bit;

  /* added extra bit to mantissa forming 1.frac */
  *mant_bits_read = *mant_bits_read + 1;

}
else if(format == 3 || format == 4){

  /* IEEE floating point number */

  /* For IEEE number, add in 1 above fractional mantissa */
  one = 1 ;
  one = one << *mant_bits_read ;

  *out_mant= (in_number & mant_mask) | one ;

  /* just expo, removes trailing 0s, removes offset*/
  *out_expo= ((in_number & expo_mask) >> *mant_bits_read ) - expo_bias;

  /* mant has grown by one bit */
  *mant_bits_read = *mant_bits_read + 1;
}
else {
/* For integers */
  *out_mant = mant_mask & in_number;
  *out_expo = 0;
}
}
/* NB - DATA FROM READ TO WRITE
sign:  1 bit representing sign of mantissa 0=+, 1=-
expo:  without any offset
mantissa: explicit 1 before fraction, normalised
mant_bits_read:  no. of bits in mantissa
offset:POSITIVE, no. of bits to left integer number is offset in string
size:  bit size of number read / wrote */

}

/****************************************************************/
/* The following function is to be called from C.  It           */
/* write out numbers of different types:                        */
/* real/integer, IBM or IEEE with offsets/packed etc            */
/* Called from the portable data conversion routines            */
/****************************************************************/
void write_number
(format,size,offset,out_number,in_sign,in_expo,in_mant,mant_bits_read)
int format, size, offset;
integer *out_number;
integer in_sign, in_expo, in_mant, mant_bits_read;
{
integer mant_bits_write, expo_bits_write;
integer sign_mask, expo_mask, mant_mask;
integer sign, mant, expo, conv_number, conv_number_mask;
integer first_bit, temp;
integer one, add_one ;
integer mant_bits_diff, expo_bias;
int i, expo_diff;
int word_length=8*sizeof(integer) ; /* Length of word in bits */

switch(format){
  case 1:
    /* format 1 IBM 32 bit */
    sign= in_sign << 31;
    mant_bits_write= ibm_mant_bits_32;
    expo_bits_write= ibm_expo_bits_32;
    sign_mask= ibm_sign_32;
    expo_mask= ibm_expo_32;
    mant_mask= ibm_mant_32;
    expo_bias= ibm_expo_bias_32;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 2:
    /* format 2 IBM 64 bit */
    sign= in_sign << 63;
    mant_bits_write= ibm_mant_bits_64;
    expo_bits_write= ibm_expo_bits_64;
    sign_mask= ibm_sign_64;
    expo_mask= ibm_expo_64;
    mant_mask= ibm_mant_64;
    expo_bias= ibm_expo_bias_64;
  break;
*ENDIF

  case 3:
    /* format 3 IEEE 32 bit */
    sign= in_sign << 31;
    mant_bits_write= ieee_mant_bits_32;
    expo_bits_write= ieee_expo_bits_32;
    sign_mask= ieee_sign_32;
    expo_mask= ieee_expo_32;
    mant_mask= ieee_mant_32;
    expo_bias= ieee_expo_bias_32;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 4:
    /* format 4 IEEE 64 bit */
    sign= in_sign << 63;
    mant_bits_write= ieee_mant_bits_64;
    expo_bits_write= ieee_expo_bits_64;
    sign_mask= ieee_sign_64;
    expo_mask= ieee_expo_64;
    mant_mask= ieee_mant_64;
    expo_bias= ieee_expo_bias_64;
  break;
*ENDIF

  case 6:
    /* format 6 IEEE 16 bit integers */
    sign= in_sign << ieee_int_mant_bits_16;
    mant_bits_write= ieee_int_mant_bits_16;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_16;
    mant_mask= ieee_int_mant_16;
  break;

  case 7:
    /* format 7 IEEE 32 bit integers */
    sign= in_sign << ieee_int_mant_bits_32;
    mant_bits_write= ieee_int_mant_bits_32;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_32;
    mant_mask= ieee_int_mant_32;
  break;
*IF DEF,FRL8,OR,DEF,CRAY

  case 8:
    /* format 8 IEEE 64 bit integers */
    sign= in_sign << ieee_int_mant_bits_64;
    mant_bits_write= ieee_int_mant_bits_64;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_64;
    mant_mask= ieee_int_mant_64;
  break;
*ENDIF

}

if ( in_mant == 0 && in_expo == 0 ) {
  mant = 0 ;
  expo = 0 ;
}
else {
/* write for real types - formats 1-4 inclusive */

if(format == 1 || format == 2){
  /* IBM numbers - reposition mantissa and re-scale expo to base 16 */

  /*  Find base 16 expo closest to but just larger than in_expo */

  if ( in_expo > -4 ) {
    expo = in_expo/4 + 1 ;

    if (in_expo < 0) {
      /* Correct exponent for negative powers of 2 */
      expo-- ;
    }

    expo_diff = in_expo - 4*expo ;
  }
  else {
    expo = in_expo/4 ;
    expo_diff = in_expo - 4*expo ;
    if (expo_diff == 0){
      expo++ ;
      expo_diff = in_expo - 4*expo ;
    }
  }
  mant = in_mant >> ( 0 - expo_diff ) ;

  /* mant has shrunk back again, lose extra bit */

  mant_bits_read--;
  expo = expo + expo_bias;

  /* Truncate/grow mantissa as appropriate to write_number format */
  mant_bits_diff = mant_bits_write - mant_bits_read;

  /* check for rounding */
  add_one = 0 ;
  if (mant_bits_diff < 0) {
    mant = (mant >> (0 - mant_bits_diff - 1) ) ;
    if (  mant & 1  == 1 ) {
      add_one = 1 ;
    }
    mant = (mant >> 1) + add_one ;
  }
  else if(mant_bits_diff > 0) mant = ( mant << mant_bits_diff );
  else                        mant = mant;

}
else if (format == 3 || format == 4 ){

  /* For IEEE numbers knock off the 1 before the decimal -
  mantissa is only the fraction and take one from mant_bits_read, as
  it effectively shrinks one bit */

  one = 1 ;
  mant = in_mant ^ (one << (mant_bits_read - 1) );
  mant_bits_read--;

  expo = in_expo + expo_bias;

  /* truncate/grow mantissa as appropriate to write_number format */
  mant_bits_diff = mant_bits_write - mant_bits_read;

  /* check for rounding */
  add_one = 0 ;
  if (mant_bits_diff < 0) {
    mant = (mant >> (0 - mant_bits_diff - 1) ) ;
    if (  mant & 1  == 1 ) {
      add_one = 1 ;
    }
    mant = (mant >> 1) + add_one ;
  }
  else if(mant_bits_diff > 0) mant = (mant << mant_bits_diff );
  else                        mant = mant;

}

else if(format>5){

  /* Write for integer types - formats 6-8 inclusive */
  /* NB integers are written in two's complement form.
     The Most Sig Bit contains sign info (0=+).
     Neg numbers are converse of their pos opposite + 1  */

  /* Positive integers */
  first_bit=0;
  if (in_sign == 0){
    /* Find first data bit */
    for (i=1;first_bit==0;i++){
      if (( in_mant >> (mant_bits_read - i) ) & 1  == 1) first_bit = i;
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }

  /* Negetive integers */
  else{
    /* Find first data bit */
    temp = 0;
    for (i=1;first_bit==0;i++){
      if ( (in_mant >> (mant_bits_read - i)) < ( (temp<<1) + 1) )
          first_bit = i;
      temp = in_mant >> (mant_bits_read - i);
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }

  /* Shrink/truncate as appropriate */
  mant = in_mant & mant_mask;

  /* If no is negative add ones above the Most Significant Bit
     to fill out two's complement form */

  /* Only need to fill out bits before bit 0 and bit mant_bits_read */

  if(in_sign != 0  &&  mant_bits_write > mant_bits_read){
    one = 1 ;
    mant = mant | ( ( (one << (mant_bits_write - mant_bits_read)) - 1 )
                                << mant_bits_read ) ;
    /* replace sign bit on mant */
    /* mant = mant | ( one << (mant_bits_read - 1) ); */
  }

  /* Warn if siginificant bits are lost
  mant_bits_diff = mant_bits_write - mant_bits_read;
  if ( (mant_bits_diff + first_bit) > 1) WARN
  */

  /* Null unused expo for build up number */
  expo = 0;
}
}

/* Build up converted number */
sign= sign_mask & sign;
mant= mant_mask & mant;
expo= expo_mask & (expo << (mant_bits_write) );
conv_number= sign | expo | mant;

/* Dealing with different sizes and offsets within the longer
   data type */
/* i.e.,  format = 32 bit int, but size = 8, offset = 8 */
/*  ----------------11111111--------  */
/*  32 bit data     8b size  8b offest
    ^ this kind of packing is used in the UM */

  if (offset != 0 || size != (mant_bits_write + expo_bits_write + 1)
      || (word_length - offset - size) != 0) {
  /* create mask for converted number */
    one = 1;
    conv_number_mask = (one << size) - 1;
  }
  else {
    conv_number_mask = -1 ;
  }

  /* apply offset to mask and number */
  conv_number_mask = conv_number_mask << offset;
  conv_number = conv_number << offset;

  /* write in the correct place within the original string */
  *out_number = conv_number | ( *out_number & (~conv_number_mask) );

}

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine STRMOV -         */
/* mirrored parameters and functionality: shifts bytes from     */
/* one variable to another.                                     */
/****************************************************************/
void
*IF DEF,C_LOW
movebytes
*ELSEIF DEF,C_LOW_U
movebytes_
*ELSE
MOVEBYTES
*ENDIF
(src,isb,num,dest,idb)
integer src[], *isb, *num, dest[], *idb;
/* src:  source variable
   isb:  starting byte of moved byte, bytes numbered from 1 from left
   num:  number of bytes moved
   dest: destianation variable
   idb:  starting byte in destination variable */
/* MOVEBYTES assumes kind=8 ie 64bit data */
{
  /* Local variables */
  integer  bytes, bytes_mask, dest_temp, one;
  integer  nbytes_moved, nbytes_to_move, num_left ;
  integer  src_start, dest_start ;
  integer  word_length=sizeof(integer) ; /* Length of word in bytes */
  int      i ;

  nbytes_moved = 0 ;
  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;

  /* Check which element we're reading from and adjust accordingly */
  for (src_start = *isb ; src_start>word_length ;
                          src_start=src_start-word_length) {
    src++ ;
  }

  /* Check which element we're writing to and adjust accordingly */
  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {

    /* First, set the nbytes_to_move to the number of bytes in
       src[src_index] between byte=src_start and byte=word_length */

    nbytes_to_move = word_length - src_start + 1;

    /* Check whether nbytes_to_move is greater than the number
       left to be moved.  If it is, then set equal to num_left */

    if (nbytes_to_move > num_left) nbytes_to_move = num_left ;

    /* Check whether there is enough room in dest[dest_index]
       between byte=dest_start and byte=8*word_length to insert
       nbytes_to_move and, if not, then set equal to the number
       of bytes available */

    if (nbytes_to_move > (word_length - dest_start + 1) ) {
       nbytes_to_move = word_length - dest_start + 1 ;
    }

    /* if (nbytes_to_move > (word_length - src_start + 1) ) {
       nbytes_to_move = word_length - src_start + 1;
    } */

    /* Move nbytes_to_move bytes from src_start in *src to dest_start
       in *dest  */

    /* make mask of same length as bits */

    if (nbytes_to_move == word_length) {
    /* This is effectively a copy dest[dest_index] = src[src_index] */
      bytes_mask = -1 ;
    }
    else{
      bytes_mask = 0;
      one = 1 ;
      bytes_mask = (one << (word_length * nbytes_to_move) ) - 1;
    }

    /* grab bits from src */

    bytes =
         (*src >>
             word_length*(8 - (src_start-1) - nbytes_to_move ) )
          & bytes_mask;

    /* put bits in dest */
    dest_temp = *dest  &
                ( ~(bytes_mask << word_length*( 8 -
                   (dest_start-1) - nbytes_to_move ) ) );

    *dest =  dest_temp |
        (bytes << word_length * ( 8 -
                    (dest_start-1) - nbytes_to_move ) );

    /* Update number of bytes moved and number left to move */

    nbytes_moved = nbytes_moved + nbytes_to_move ;
    num_left = num_left - nbytes_to_move ;

    /* Update src_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    src_start = src_start + nbytes_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;
      /* Increment to next data element of src */
      src++ ;
    }

    /* Update dest_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    dest_start = dest_start + nbytes_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;
      /* Increment to next data element of dest */
      dest++ ;
    }
  }
}

/****************************************************************/
/* The following function is to be called from Fortran.  It     */
/* provides a portable version of Cray routine MOVBIT -         */
/* mirrored parameters and functionality: shifts bytes from     */
/* one variable to another.                                     */
/****************************************************************/
void
*IF DEF,C_LOW
movebits
*ELSEIF DEF,C_LOW_U
movebits_
*ELSE
MOVEBITS
*ENDIF
(src,isb,num,dest,idb)
integer src[], *isb, *num, dest[], *idb;
/* src:  source variable
   isb:  starting bit of moved bits, bits numbered from 1 from left
   num:  number of bits moved
   dest: destination variable
   idb:  starting bit in destination variable */

{
  /* Local variables */

  integer  bits, bits_mask, dest_temp, one;
  integer  nbits_moved, nbits_to_move, num_left ;
  integer  src_start, dest_start ;
  integer  word_length=8*sizeof(integer) ; /* Word length in bits */
  int      i ;

  nbits_moved = 0 ;
  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;

  /* Check which element we're reading from and adjust accordingly */
  for (src_start = *isb; src_start>word_length ;
                         src_start=src_start-word_length) {
    src++ ;
  }

  /* Check which element we're writing to and adjust accordingly */
  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {

    /* First, set the nbits_to_move to the number of bits in
       src[src_index] between bit=src_start and bit=word_length */

    nbits_to_move = word_length - src_start + 1;

    /* Check whether nbits_to_move is greater than the number
       left to be moved.  If it is, then set equal to num_left */

    if (nbits_to_move > num_left) nbits_to_move = num_left ;

    /* Check whether there is enough room in dest[dest_index]
       between bit=dest_start and bit=word_length to insert
       nbits_to_move and, if not, then set equal to the number
       of bits available */

    if (nbits_to_move > (word_length - dest_start + 1) ) {
       nbits_to_move = word_length - dest_start + 1 ;
    }

    /* Move nbits_to_move bits from src_start in *src to dest_start
       in *dest  */

    /* make mask of same length as bits */

    if (nbits_to_move == word_length) {
    /*  This is effectively a copy: dest[dest_index]=src[src_index] */
      bits_mask = -1 ;
    }
    else{
      bits_mask = 0;
      one = 1 ;
      bits_mask = (one << nbits_to_move) - 1;
    }

    /* grab bits from src */
    bits = (*src >> ( word_length - src_start - nbits_to_move + 1) )
         & bits_mask;

    /* put bits in dest */
    dest_temp = *dest
         & ( ~(bits_mask <<
              ( word_length - dest_start - nbits_to_move +1)) );
    *dest =  dest_temp |
                ( bits <<
                ( word_length - dest_start - nbits_to_move +1) );

    /* Update number of bits moved and number left to move */

    nbits_moved = nbits_moved + nbits_to_move ;
    num_left = num_left - nbits_to_move ;

    /* Update src_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    src_start = src_start + nbits_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;
      /* Increment to next data element of src */
      src++ ;
    }

    /* Update dest_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    dest_start = dest_start + nbits_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;
      /* Increment to next data element of dest */
      dest++ ;
    }
  }
}
