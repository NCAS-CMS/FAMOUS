#include <stdio.h>
#include <stdlib.h>


void ctest(float *);
void ctest_(float *);
void CTEST(float *);
void ctest_main(float [2],double [2]);


void ctest(float array[])
{
printf("Fortran/C Interface:\n");
printf("Lowercase, no underscore - use nupdate *DEF C_LOW\n");
printf("FC_LINK=C_LOW\n");
ctest_main(array,(double *)array);
}

void ctest_(float array[])
{
printf("Fortran/C Interface:\n");
printf("Lowercase with underscore - use nupdate *DEF C_LOW_U\n");
printf("FC_LINK=C_LOW_U\n");
ctest_main(array,(double *)array);
}

void CTEST(float array[])

{
printf("Fortran/C Interface:\n");
printf("Uppercase, no underscore - no nupdate *DEF required\n");
printf("FC_LINK=C_UP\n");
ctest_main(array,(double *)array);
}

void ctest_main(float array1[2],double array2[2])

{
int word_size;
long long_integer;

unsigned char *one_byte;
int full_int=1,first_val,last_val;
  
printf("\nFortran/C REAL types:\n");

/*
printf("array1: value1: %f and value2: %f \n",array1[0],array1[1]);
printf("array2: value1: %f and value2: %f \n",array2[0],array2[1]);
*/

if (array1[0]==array1[1]) {
  printf("Fortran type \"REAL\" is equivalent to C type \"float\"\n");
  printf("FRL8=false\n");
  word_size=sizeof(float)*8;}
else {
  if (array2[0]==array2[1]) {
    printf("Fortran type \"REAL\" is equivalent to C type \"double\"\n");
    printf("FRL8=true\n");
    word_size=sizeof(double)*8;}
  else {
    printf("Could not determine which C type is equivalent to Fortran REAL type");
    printf("FRL8=null\n");}
     }  
     
printf("\nFortran REAL type is %d bits\n",word_size);
printf("FRL_SIZE=%d\n",word_size);

printf("\nFortran/C INTEGER types:\n");
if (word_size==64) {
  if (sizeof(long_integer)!=8) {
/*    printf("C integer type \"long\" is not 64 bits.\n");*/
    printf("Fortran type \"INTEGER\" is equivalent to C type \"long long\"\n");
    printf("INTLL=true\n"); }
  else { 
/*    printf("C integer type \"long\" is 64 bits.\n");*/
    printf("Fortran type \"INTEGER\" is equivalent to C type \"long\"\n");
    printf("INTLL=false\n"); } }
else if (word_size==32) {
/*  printf("C integer type \"int\" is 32 bits.\n"); */
  printf("Fortran type \"INTEGER\" is equivalent to C type \"int\"\n");
  printf("INTLL=false\n"); }
else {
  printf("Could not determine which C type is equivalent to Fortran INTEGER type");
/*  printf("C integer type could not be determined\n"); */
  printf("INTLL=null\n"); }

    
printf("\nThe byte ordering on this machine ");

one_byte=(unsigned char *) &full_int;
  
first_val=*one_byte;
one_byte+=sizeof(int)-1;
last_val=*one_byte; 

  if (first_val == 1) {
    printf("is Little Endian\n");
    printf("BIGEND=false\n"); }
  else
    { if (last_val == 1) {
    	printf("is Big Endian\n");
        printf("BIGEND=true\n"); }
   	  else {
            printf("could not be determined\n"); 
            printf("BIGEND=null\n"); }
   	}
}

