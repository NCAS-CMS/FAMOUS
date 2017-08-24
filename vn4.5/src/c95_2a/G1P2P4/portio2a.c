/* ******************************COPYRIGHT******************************
 * (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
 *
 * Use, duplication or disclosure of this code is subject to the
 * restrictions as set forth in the contract.
 * 
 *                Meteorological Office
 *                London Road
 *                BRACKNELL
 *                Berkshire UK
 *                RG12 2SZ
 * 
 * If no contract has been raised with this copy of the code, the use,
 * duplication or disclosure of it is strictly prohibited. Permission
 * to do so must first be obtained in writing from the Head of 
 * Numerical Modelling at the above address.
 * ******************************COPYRIGHT******************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>


#include <time.h>

     /* C language routines for portable version of UM */
     /*       Written by A.Dickinson 1/11/91           */
     /* Model    Modification history from version 3.0  */
     /* version  Date                                   */
     /* 3.1      15/01/93 Increase size of available    */
     /*          unit nos. from 1-100 to 1-199.         */
     /*          R.Rawlins                              */
     /*                                                 */
     /*          A range of C - Fortran interfaces      */
     /*          provided via *DEFS C_LOW & C_LOW_U.    */
     /*          See UM Doc Paper S5.                   */
     /*                                                 */
     /*          A.Dickinson 24/06/93                   */
     /*                                                 */
     /* 3.3      Error in dimension of array fno in     */
     /*          routine get_file corrected. Change type*/
     /*          in date_time to remove compilation     */
     /*          warnings. A.Dickinson 06/11/93         */
     /* 3.3       07/09/93 Set buffer in OPEN. R.Rawlins */
     /* 3.3      04/10/93 New routine FLUSH_BUFFER to   */
     /*          force buffer write. R. Rawlins         */
     /* 3.3       09/08/94 Include test in CLOSE to check*/
     /*           unit flag and by-pass close if unit    */
     /*           already closed or not yet opened.      */
     /*           Reset unit flag on normal close of     */
     /*           file. {close unit a then open unit b   */
     /*           then close unit a caused file pointer  */
     /*           for unit b to become unspecified.      */
     /*           R. Rawlins                             */
     /* 3.4      04/10/94 Improvement to code in getfile*/
     /*                   Tracey Smith                  */
  /* 3.5  23/03/95 Tidy up the interface between Fortran REALs */
  /*           and C float/double.                             */
  /*           Rename open and close to file_open + file_close */
  /*           T3D specific code for CHARACTER arguments       */
  /*           Change name of some routines for mpp code       */
  /*           Add routine FORT_GET_ENV to get environment     */
  /*           variables from Fortran                          */
  /*                Author : Paul Burton                       */
 /* 4.0 Error code argument added to setpos, close M.TURP */

   /*  Added perror call to BUFFIN and BUFFOUT */
   /*  Gordon Henderson */
  /* 4.1  11/04/96  FILE_OPEN: Check the return code from getenv    */
  /*                           and trap NULL                        */
  /*                FILE_CLOSE: Ditto                               */
  /*                GET_FILE: Checks getenv return code and if NULL */ 
  /*                          sets filename argument to blank       */
  /*                FORT_GET_ENV: Added code for T3D functionality  */
  /*      Paul Burton  */
  /*                                                           */
  /*           Make the necessary changes to use FFIO and      */
  /*           character variables correctly on all Cray       */
  /*           Platforms.                                      */
  /*                                                           */
  /*                Bob Carruthers, Cray Research U.K.         */
  /*                                                            */
  /* 4.4  17/06/97  Modify the write/print statements to use    */
  /*                a standard mechanism, set up via a #define  */
  /*                statement                                   */
  /*                  Author: Bob Carruthers, Cray Research     */
  /*                                                            */
  /* 4.4  17/06/97  Add code to accept and use the current      */
  /*                Length for a dumpfile                       */
  /*                  Author: Bob Carruthers, Cray Research     */
  /*                                                            */


  /*           Changes to GET_FILE to allow unit numbers     */
  /*           greater > 100 to be used as enviroment        */
  /*           variables for filenames.     Ian Edmond       */
  /* 4.5  01/04/98  Assorted mods to the C code:                */
  /*                BUFFIN32: New routine added, and check that */
  /*                          unit is open added to BUFFO32     */
  /*                SETPOS32: New routines added for setting/   */
  /*                GETPOS32: reading file pointer for 32bit    */
  /*                          word files                        */
  /*                                                            */
  /*                Authors: Bob Carruthers & Paul Burton       */
  /*                                                            */

  /*                                                            */
  /* 4.5  31/03/98  Modify various definitions to make porting  */
  /*                to platforms that fully support 32 and 64   */
  /*                integers easier.                            */
  /*                  Author: Bob Carruthers, Cray Research     */
  /*                                                            */
  /* 4.5  16/07/98  Add properties flag word for each unit -    */
  /*                initial use is to broadcasts in buffin for  */
  /*                the AC scheme.                              */
  /*                  Author: Bob Carruthers, Cray Research     */
  /*                                                            */
  /*                                                            */
  /* 4.5  17/08/98  Code to ensure that we do not over-index    */
  /*                character arrays in Fortran or C            */
  /*                  Author: Bob Carruthers, Cray Research     */
  /* 5.0  26/07/99  Fix error in declaration of fname in        */
  /*                FILE_OPEN and FILE_CLOSE                    */
  /*                Author: Paul Selwood                        */

typedef double real;
typedef long integer;

/* Fortran REAL is equivalent to C double */

/* Define the function that outputs the text string for         */
/* messages.  Note that the trailing newline character is       */
/* to be supplied by the print routine, and is no longer        */
/* in the string.  Similarly, leading new lines are now handled */
/* by the print routine - typically a newline is inserted for   */
/* each change of unit.                                         */

#define CALL_MESSAGE_PRINT(text) fprintf(stdout, "%s\n", text);        \
 fflush(stdout)

static char message[256];

#define MAX_UNITS 200

integer *the_unit;


FILE *pf[MAX_UNITS]=
                     {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL
                     };

static integer io_position[MAX_UNITS];
int open_flag[MAX_UNITS] =
                     {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/* Unit properies table - one word per unit, one bit per
   property at present */

integer file_properties[MAX_UNITS] =
                     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

/* Define the bits used in the properies table */

#define BCAST 1

extern void clear_unit_bcast_flag_();


void
buffin_single_
(unit, array, maxlen, length, status)
integer *unit;      /* Fortran unit                         */
real array[];       /* Array into which data is read        */
integer *maxlen;    /* Number of real numbers to be read    */
integer *length;    /* Number of real numbers actually read */
real *status;       /* Return code                          */
{
int k;
  the_unit=unit;

  if(open_flag[*unit]== 0){
    *length = fread(array,sizeof(real),*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      sprintf(message,
       "BUFFIN: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      sprintf(message,
       "BUFFIN: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

    io_position[*unit]=io_position[*unit]+*length;
}

void
buffout_single_
(unit, array, maxlen, length, status)
integer *unit;      /* Fortran unit                            */
real array[];       /* Array from which data is written        */
integer *maxlen;    /* Number of real numbers to be written    */
integer *length;    /* Number of real numbers actually written */
real *status;       /* Return code                             */
{
int k;
  the_unit=unit;

  if(open_flag[*unit]== 0){
    *length = fwrite(array,sizeof(real),*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      sprintf(message,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      sprintf(message,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

    io_position[*unit]=io_position[*unit]+*length;
}

void
setpos_single_
(unit, word_address,err)
integer *unit;      /* Fortran unit                         */
integer *word_address; /* Number of words into file            */
integer *err;          /* Error checking, err = 0 no errors, 
                                          err = 1 errors */
{
int k;
integer byte_address;
  the_unit=unit;

    byte_address=(*word_address)*sizeof(real);
    k = fseek(pf[*unit],byte_address,SEEK_SET);
    *err = 0;
    if(k!=0){
      perror("\nSETPOS: Seek Failed");
      sprintf(message,
       "SETPOS: Unit %d to Word Address %d Failed with Error Code %d",
       (int) *unit, (int) *word_address, k);
      CALL_MESSAGE_PRINT(message);
         *err = 1;
         abort();
    }
   io_position[*unit]=*word_address;

}

void
open_single_
(unit,file_name, char_len,intent,environ_var_flag,err)
char file_name[]; /* File name or environment variable    */
integer *unit;       /* Fortran unit                         */
integer *char_len;   /* No of chars in file name             */
integer *intent    ; /* =0 read only,!=0 read and write      */
integer *environ_var_flag; /* =0 file name in environment var, */
                     /*!=0 explicit file name            */
integer *err;        /* =0 file opened,                  */
                     /*!=0 file not opened               */
{
   char *fname;
   char *gname;
   int i;

   int readonly = 0;

   enum   filestat { old, new };
   enum   filestat filestatus;


   fname = calloc(*char_len + 1, 1);
   the_unit=unit;
   pf[*unit] = NULL;
/* convert file name to C format */

   strncpy( fname, file_name, *char_len );
   fname[ *char_len ] = '\0';
   sscanf( fname, "%s", fname );


   if ( *environ_var_flag == 0 )  /* File name held in environ var */
   {  gname = getenv( fname );
      if ( gname == NULL ) {
        sprintf(message,
         "OPEN:  WARNING: Environment variable %s not set",
         fname);
        CALL_MESSAGE_PRINT(message);
        open_flag[*unit]=1;
        *err=1;
        free (fname);
        return;
      }
   }
   else                           /* get file name from argmt fname */
      gname = fname;


   /* Check if file exists */

   if ( access( gname, 0 ) == 0 )  {   /* file exists */

      sprintf(message,
       "OPEN:  File %s to be Opened on Unit %d Exists",
       gname, (int) *unit );
      CALL_MESSAGE_PRINT(message);
      filestatus = old;
   }
   else  {   /* non-existent file */

      sprintf(message,
       "OPEN:  File %s to be Opened on Unit %d does not Exist",
       gname, (int) *unit);
      CALL_MESSAGE_PRINT(message);
      filestatus = new;
   }


   if ( filestatus == old )  {

      if ( *intent == readonly )  {

         if ( ( pf[*unit] = fopen( gname, "rb" ) ) == NULL )  {

            perror("OPEN:  File Open Failed");
            sprintf(message,
              "OPEN:  Unable to Open File %s for Reading", gname );
            CALL_MESSAGE_PRINT(message);
         }
      }
      else  {   /*  *intent == read_and_write )  */

         if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == NULL )  {

            perror("OPEN:  File Open Failed");
            sprintf(message,
              "OPEN:  Unable to Open File %s for Read/Write", gname );
            CALL_MESSAGE_PRINT(message);
         }
      }
   }


/* New file - check for write */
   if ( filestatus == new )  {

/* Initialise the file control word to NULL */
      pf[*unit] = NULL;

      if ( *intent == readonly )  {
         sprintf(message, "OPEN:  **WARNING: FILE NOT FOUND" );
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "OPEN:  Ignored Request to Open File %s for Reading",
           gname );
         CALL_MESSAGE_PRINT(message);
      }
      else  {        /*  *intent == read_and_write   */

/* File size not given - just open the file normally */
          if ( ( pf[*unit] = fopen( gname, "w+b" ) ) == NULL )  {

            perror("OPEN:  File Creation Failed");
            sprintf(message,
             "OPEN:  Unable to Open File %s for Read/Write", gname );
            CALL_MESSAGE_PRINT(message);
          }
          else  {
            sprintf(message, "OPEN:  File %s Created on Unit %d",
                     gname, (int) *unit );
            CALL_MESSAGE_PRINT(message);
          }
      }
   }



   /* Set error code and open flag used by buffin and buffout */

   if( pf[*unit] == NULL )  {

      *err = 1;
      open_flag[*unit]=1;
   }
   else  {
      if ( setvbuf( pf[*unit], NULL, _IOFBF, BUFSIZ ) != 0 )  {
         perror("\n**Warning: setvbuf failed");
         *err=1;
         open_flag[*unit]=1;
       }
       else
       {
         *err = 0;
         open_flag[*unit]=0;

/*    set buffer to default size to force buffer alloc on heap */
  /*  setvbuf(pf[*unit],NULL,_IOFBF,BUFSIZ);  See above */
        }
   }
   io_position[*unit]=0;

   clear_unit_bcast_flag_(unit);

free (fname);
}



void
close_single_
(unit,file_name,char_len,environ_var_flag,delete,err)
char file_name[];    /* File name or environment variable    */
integer *unit;       /* Fortran unit                         */
integer *char_len;   /* No of chars in file name             */
integer *delete;     /* =0 do not delete file,!=0 delete file*/
integer *environ_var_flag; /* =0 file name in environment var, */
                           /*!=0 explicit file name            */
integer *err; /* ERROR CHECKING err = 0 no errors, err = 1 Errors */
{
char *fname;
char *gname;
int i;
int k;

the_unit=unit;
fname = calloc(*char_len + 1, 1);
/* first check to see if unit has been closed already (or not opened)*/
if(open_flag[*unit]== 0){    /* unit currently open  */

/* close file */
      k=fclose(pf[*unit]);

/* convert file name to C format */
        strncpy(fname,file_name,*char_len);
        fname[*char_len] = '\0';
        for (i=0; i<*char_len; i++){

            if (fname[i] == ' '){
               fname[i] = '\0';
               break;
            }
         }

        if(*environ_var_flag == 0)
          { gname = getenv( fname );
            if ( gname == NULL ) {
            open_flag[*unit]=1;
            *err=1;
            free( fname );
            return;
            }
          }
        else
          gname=fname;

      if(k==0){

/* delete file */
        if(*delete != 0){
          k=remove(gname);
          if( k != 0){
            *err = 1;
            abort();
          }
          else{  /*normal end to delete so:*/
            open_flag[*unit]=1;     /* set unit flag to closed */
            *err = 0;
          }

        }
        else{
/* file closed */
           open_flag[*unit]=1;     /* set unit flag to closed */
        }
      }
/* file not closed */
}   /* end of test for unit already closed */

else {
/* unit either closed already or not open yet */
  sprintf(message,
   "CLOSE: WARNING: Unit %d Not Opened", (int) *unit);
  CALL_MESSAGE_PRINT(message);
}

free( fname );
}


void
date_time_
(year,month,day,hour,minute,second)
integer *year;       /* year                  */
integer *month;      /* month   Current date  */
integer *day;        /* day      and time     */
integer *hour;       /* hour                  */
integer *minute;     /* minute                */
integer *second;     /* second                */

{
char s[5];
time_t t,*r,a;
    r=&a;
    t=time(r);
    strftime(s,5,"%Y",localtime(r));
    *year=atoi(s);
    strftime(s,5,"%m",localtime(r));
    *month=atoi(s);
    strftime(s,5,"%d",localtime(r));
    *day=atoi(s);
    strftime(s,5,"%H",localtime(r));
    *hour=atoi(s);
    strftime(s,5,"%M",localtime(r));
    *minute=atoi(s);
    strftime(s,5,"%S",localtime(r));
    *second=atoi(s);
}
void
shell_
(command,command_len)
char command  []; /* Command to be executed               */
integer *command_len; /* No of chars in command               */
{
char *fname;
int i;

/* convert file name to C format */

fname = calloc( *command_len + 1, 1);
strncpy(fname,command,*command_len);
fname[*command_len]='\0';

/* execute command */
        i=system(fname);

free( fname );
}

void
buffo32_
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                            */
real array[];      /* Array from which data is written        */
integer *maxlen;   /* Number of real numbers to be written    */
integer *length;   /* Number of real numbers actually written */
real *status;      /* Return code                             */
{
int k;


    if (open_flag[*unit]==0){
    *length = fwrite(array,4,*maxlen,pf [*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFO32\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFO32\n");
      printf("Return code = %d\n",k);
      *status=1.0;
      }
    }
      else
        *status=3.0;

}

void
buffin32_
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                         */
real array[];      /* Array from which data is read        */
integer *maxlen;   /* Number of real numbers to be read    */
integer *length;   /* Number of real numbers actually read */
real *status;      /* Return code                          */
{
int k;

    if (open_flag[*unit]==0){

      *length = fread(array,4,*maxlen,pf [*unit]);

      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        printf("C I/O Error: failed in BUFFIN32\n");
        printf("Return code = %d\n",k);
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        printf("C I/O Error: failed in BUFFIN32\n");
        printf("Return code = %d\n",k);
        *status=1.0;
      }
    }
    else
      *status=3.0;
}

void
buffin8_
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                         */
char array[];      /* Array into which data is read        */
integer *maxlen;   /* Number of bytes to be read           */
integer *length;   /* Number of bytes actually read        */
real *status;      /* Return code                          */
{
int k;

  if(open_flag[*unit]== 0){
    *length = fread(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFIN8\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFIN8\n");
      printf("Return code = %d\n",k);
      *status=1.0;
    }
   }
   else
        *status=3.0;

}

void
buffou8_
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                            */
char   array[];    /* Array from which data is written        */
integer *maxlen;   /* Number of bytes to be written           */
integer *length;   /* Number of bytes actually written        */
real *status;      /* Return code                             */
{
int k;

  if(open_flag[*unit]== 0){
    *length = fwrite(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFOU8\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFOU8\n");
      printf("Return code = %d\n",k);
      *status=1.0;
    }
   }
   else
        *status=3.0;

}
void
setpos8_
(unit, byte_address)
integer *unit;     /* Fortran unit                         */
integer *byte_address; /* Number of bytes into file            */
{
int k;

    k = fseek(pf[*unit],*byte_address,SEEK_SET);
    if(k!=0){
         printf("ERROR detected in SETPOS\n");
         printf("word_address = %d\n", (int) *byte_address);
         printf("Return code from fseek = %d\n",k);
    }
}

void
getpos8_
(unit, byte_address)
integer *unit;     /* Fortran unit                         */
integer *byte_address; /* Number of bytes into file            */ 
{

    *byte_address = ftell(pf[*unit]);


}
void
setpos32_
(unit,word32_address,err)

integer *unit;            /* Fortran unit                         */
integer *word32_address;  /* Number of 32bit words into the file  */
integer *err;             /* 0: no error
                             1: error occured                     */

{
int k;
integer byte_address;

  *err=0;
  if (open_flag[*unit]==0) {
    byte_address=(*word32_address)*4;

    k=fseek(pf[*unit],byte_address,SEEK_SET);

    if (k != 0) {
      k=errno;
      perror("\nSETPOS32: Seek failed");
      sprintf(message,
        "SETPOS32: Unit %d to 32bit Word Address %d failed. Error: %d",
        (int) *unit, (int) *word32_address, k);
      the_unit=unit;
      CALL_MESSAGE_PRINT(message);
      *err=1;
      abort();
    }

  }
}

void
getpos32_
(unit,word32_address)

integer *unit;            /* Fortran unit                         */
integer *word32_address;  /* Number of 32bit words into the file  */

{
int byte_address;

  if (open_flag[*unit]==0) {
    byte_address=ftell(pf[*unit]);
    *word32_address = byte_address/4;
  }
}

void
getpos_
(unit, word_address)
integer *unit;     /* Fortran unit                         */
integer *word_address; /* Number of words into file            */
{
long byte_address;

   the_unit=unit;
     byte_address = ftell(pf[*unit]);
    *word_address = byte_address/sizeof(real);
      if(*word_address != io_position[*unit]) {
        sprintf(message,
         "GETPOS: IO_POSITION is %d, but FTELL gives %d",
         (int) io_position[*unit], (int) *word_address);
        the_unit=unit;
        CALL_MESSAGE_PRINT(message);
        abort();
      }

}

void
word_length_
(length)
integer *length;  /* Word length used by hardware         */
{

    *length=sizeof(real);

}
void
get_file_
(unit,filename,file_len,err)
char filename[]; /* File name                           */
integer *file_len; /* Dimension of filename               */
integer  *unit;    /* Fortran unit number                 */
integer  *err; /* Error checking err = 0 no errors, err = 1 errors */
{
char fname[16];
char fno[4];
char *gname;
int i;
int k;

the_unit=unit;
/* construct environment variable name         */
/* in form UNITnn, where nn is Fortran unit no */

       if ( *unit < 100){
       strcpy (fname, "UNIT");
       sprintf(fno,"%02i", (int) *unit);
       strcat(fname,fno);
       fname[6]='\0';
       }
       else{
       strcpy (fname, "UNIT");
       sprintf(fno,"%03i", (int) *unit);
       strcat(fname,fno);
       fname[7]='\0';
       }

/* get file name stored in environment variable UNITnn */
       gname=getenv(fname);
       if ( gname == NULL) {
         sprintf(message,
          "GET_FILE: Environment Variable %s not Set", fname);
         CALL_MESSAGE_PRINT(message);
         filename[0] = '\0';
         for (i=1; i<*file_len; i++){
                filename[i] = ' ';
         }
         return;               
       }
       k=strlen(gname);
       if(k >  *file_len){
         sprintf(message,
          "GET_FILE: File Name too long for Allocated Storage");
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "GET_FILE: Environment Variable %s", fname);
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "GET_FILE: File Name %s", gname);
         CALL_MESSAGE_PRINT(message);
         *err = 1;
         abort();
       }

/* convert file name to Fortran format */
          *err = 0;
          strncpy(filename,gname,k);
          for (i=k; i<*file_len; i++){
               filename[i] = ' ';
           }


}
void
abort_()
{
  abort();
}
     /*          Force i/o buffer to be written to file */
     /*          explicitly to prevent continuation run */
     /*          problems following 'hard' failures.    */
void
flush_buffer_
(unit, icode)
integer  *unit     ;  /* Fortran unit number             */
integer  *icode    ;  /* Integer return code             */
{
int  i         ;
  if(open_flag[*unit]== 0){
      i =   fflush(pf[*unit]);
      *icode = i;
  }
  else {
    if(pf[*unit] != NULL) {
      sprintf(message,
       "FLUSH_BUFFER: File Pointer for Unopened Unit %d is %16X",
       (int) *unit, (unsigned long) pf[*unit]);
      the_unit=unit;
      CALL_MESSAGE_PRINT(message);
      abort();
    }
  }
}
void
fort_get_env_
(env_var_name,ev_len,ev_contents,cont_len,ret_code)
char *env_var_name;  /* Name of environment variable   */
integer *ev_len;     /* length of name                 */
char *ev_contents;   /* contents of environment variable */
integer *cont_len;   /* length of contents              */
integer *ret_code;   /* return code: 0=OK  -1=problems  */

{
        integer minus_one=-1;

        char *c_env_var_name;

        char *value;
        int len,i;

        c_env_var_name = calloc(*ev_len + 1, 1);
        the_unit=&minus_one;
        strncpy(c_env_var_name,env_var_name,*ev_len);
        c_env_var_name[*ev_len]='\0';
        sscanf(c_env_var_name,"%s",c_env_var_name);

        value=getenv(c_env_var_name);
        if (value==NULL){
          *ret_code=-1;
          free( c_env_var_name );
          return;}
        else{
          *ret_code=0;}

        len=strlen(value);
        if (len > *cont_len){
         sprintf(message,
          "FORT_GET_ENV: Value too long for Allocated Storage");
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "FORT_GET_ENV: Environment Variable %s",
          c_env_var_name);
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "FORT_GET_ENV: Value %s", value);
         CALL_MESSAGE_PRINT(message);
                abort();
        }

        strncpy(ev_contents,value,len);
        for (i=len; i<*cont_len; i++){
                ev_contents[i]=' ';
        }
free( c_env_var_name );

}

void
CLOSE_ALL_FILES
()
{
}

void
FLUSH_ALL_FILES
()
{
}

/*                                                              */
/* Entry to accept the current File length prior                */
/* to an open request                                           */
/*                                                              */
void
set_dumpfile_length_
(unit, length)
integer *unit;
integer *length;
{
}
void
clear_unit_bcast_flag_
(unit)
integer *unit;
{
  integer minus_one=-1;

  file_properties[*unit]=(minus_one ^ BCAST) & file_properties[*unit];

}

void
set_unit_bcast_flag_
(unit)
integer *unit;
{

  file_properties[*unit]=BCAST | file_properties[*unit];

}


void
find_unit_bcast_flag_
(unit, flag)
integer *unit;
integer *flag; /* non-zero if set, otherwise 0 */
{

  *flag=BCAST & file_properties[*unit];

}


