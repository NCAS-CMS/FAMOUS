*ID GQZ1C505
*/--------------------------------------------------------
*/
*/ Backport vn5.5 mod gqz1505 and vn6.0 mod gel6600 to vn4.5
*/ Original author frqz (Paul Dando)
*/
*/ Jeff Cole 17/11/03
*/
*/ Modset created by frqz on 26/02/03
*/
*/ Portable I/O changes.  Allows for input and output
*/ of UM start dumps, LBCs, Fields Files, PP files etc
*/ on LITTLE_ENDIAN machines.
*/--------------------------------------------------------
*/
*/-----
*DECLARE PORTIO2A
*/-----
*/
*B GBC3C405.5
  /* 5.5  25/02/03  Portability changes to allow for input and    */
  /*                output of BIGENDIAN files on LITTLEENDIAN     */
  /*                machines. Required byte-swap only performed   */
  /*                *IF DEF,LITTLE_END. Special case included        */
  /*                to handle 32-bit packed data.       P.Dando   */
  /* 6.0  18/09/03  Fix for BUFFO32 - NB for IEEE *only*        */
  /*                                             E.Leung        */
*/
*/ buffin_single
*/
*I PORTIO2A.23
*/*IF DEF,LITTLE_END
  int i;                               /* array counter             */
  integer *ptr_integer = 0;            /* temporary integer pointer */
  void change_endian(integer *, int);  /* function to swap endian   */
*/*ENDIF
*I PORTIO2A.28
*/*IF DEF,LITTLE_END
    for (i = 0; i < *maxlen; i++) {
      ptr_integer = (integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
*/*ENDIF

*/
*/ buffout_single
*/
*I PORTIO2A.46
*/*IF DEF,LITTLE_END
  int i;                              /* array counter             */
  integer *ptr_integer = 0;           /* temporary integer pointer */
  void change_endian(integer *,int);  /* function to swap endian   */
*/*ENDIF
*I PORTIO2A.49
*/*IF DEF,LITTLE_END
    for (i = 0; i < *maxlen; i++) {
      ptr_integer=(integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
*/*ENDIF

*I PORTIO2A.51
*/*IF DEF,LITTLE_END
    for (i = 0; i < *maxlen; i++) {
      ptr_integer=(integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
*/*ENDIF

*/
*/ buffo32
*/
*I GBC2C405.83
*IF DEF,MPP
*IF DEF,C_LOW
buffo32_single
*ELSEIF DEF,C_LOW_U
buffo32_single_
*ELSE
BUFFO32_SINGLE
*ENDIF
*ELSE
*I TS140793.69
*ENDIF
*I PORTIO2A.267
*/*IF DEF,LITTLE_END
    int i;                              /* array counter             */
    integer *ptr_integer = 0;           /* temporary integer pointer */
    void change_endian(integer *,int);  /* function to swap endian   */
*/*ENDIF
*I GBC0C402.179

*/*IF DEF,LITTLE_END
*IF DEF,FRL8,OR,DEF,CRAY,OR,DEF,IEEE
     for (i = 0; i < (*maxlen + 1)/2; i++) {
*ELSE
     for (i = 0; i < *maxlen; i++) {
*ENDIF
        ptr_integer = (integer *)&array[i] ;
*IF DEF,IEEE
        change_endian(ptr_integer, sizeof(integer));
*ELSE
        change_endian(ptr_integer, 4);
*ENDIF
      }
*/*ENDIF

*I PORTIO2A.271
*/*IF DEF,LITTLE_END
*IF DEF,FRL8,OR,DEF,CRAY,OR,DEF,IEEE
      for (i = 0; i < (*maxlen + 1)/2; i++) {
*ELSE
      for (i = 0; i < *maxlen; i++) {
*ENDIF
        ptr_integer = (integer *)&array[i] ;
*IF DEF,IEEE
        change_endian(ptr_integer, sizeof(integer));
*ELSE
        change_endian(ptr_integer, 4);
*ENDIF
      }
*/*ENDIF

*/
*/ buffin32
*/
*I GPB0C405.12
*IF DEF,MPP
*IF DEF,C_LOW
buffin32_single
*ELSEIF DEF,C_LOW_U
buffin32_single_
*ELSE
BUFFIN32_SINGLE
*ENDIF
*ELSE
*I GPB0C405.19
*ENDIF
*I GPB0C405.30
*/*IF DEF,LITTLE_END
    int i;                              /* array counter             */
    integer *ptr_integer = 0;           /* temporary integer pointer */
    void change_endian(integer *,int);  /* function to swap endian   */
*/*ENDIF
*I GPB0C405.59
*/*IF DEF,LITTLE_END
*IF DEF,FRL8,OR,DEF,CRAY
      for (i = 0; i < (*maxlen + 1)/2; i++) {
*ELSE
      for (i = 0; i < *maxlen; i++) {
*ENDIF
        ptr_integer = (integer *)&array[i];
        change_endian(ptr_integer, 4);
      }
*/*ENDIF

*I GBC0C402.296
*/*IF DEF,LITTLE_END

/* In order to be consistent with present systems all data files      */
/* will be saved in big endian format.                                */
/*                                                                    */
/* The following functions are used to convert byte order (endian)    */
/* of data                                                            */
/*                                                                    */
/* Use only when reading (writing) data on little endian machines     */

void change_endian ( integer *ptr_array, int Nbytes ) {

/* Swap byte order of *ptr_array */

  int i;
  unsigned char *ptr_IVal=0;
  unsigned char *ptr_OVal=0;

  integer IVal=*(integer *)ptr_array ;
  ptr_IVal=(unsigned char *)&IVal ;

  /* unsigned char is one byte */

  ptr_OVal=(unsigned char *)ptr_array;

  for (i=0; i<Nbytes; i++) {
    /* reverse byte ordering */
    ptr_OVal[Nbytes-1-i]=ptr_IVal[i];
  }

  /* Packed data (Nbytes=4) needs the second 4 bytes swapping too */

  if ( Nbytes == 4 && Nbytes != sizeof(integer) ) {
    for (i=0; i<Nbytes; i++){
      /* reverse byte ordering */
      ptr_OVal[2*Nbytes-1-i]=ptr_IVal[Nbytes+i];
    }
  }
  return;
}

*/*ENDIF
