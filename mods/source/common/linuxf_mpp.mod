*ID GCB3F404
*/ AV - one of Bob's mods for SGI. MPI must be initialised before
*/      opening any files. On the T3E, this doesn't matter.
*/
*/ AH - mod to delete lines relating to openwing + sending of 
*/      wakeup call to channel 8
*/      Also to delete closing of channel 8.  With these lines 
*/      intact the Linux mpi code crashes after the last timestep
*/      with a p4 error.
*/
*DC UMSHELL1
*D GPB0F305.5,GPB0F305.15
*D GRR2F305.293,GRR2F305.297

*I GPB0F305.163

!   Open file for UNIT 5 before initialisation of model. All runtime
!   control variables subsequently read in from UNIT 5 by namelist.
      CALL GET_FILE(5,FILENAME,80,ICODE)
      OPEN(5,FILE=FILENAME,IOSTAT=ISTATUS)
      IF(ISTATUS.NE.0) THEN
        ICODE=500
        WRITE(6,*) ' ERROR OPENING FILE ON UNIT 5'
        WRITE(6,*) ' FILENAME =',FILENAME
        WRITE(6,*) ' IOSTAT =',ISTATUS
        GOTO 999
      END IF

CL------------------------------------------------------------------
CL 0.1 Get submodel/internal model components of model run.
CL
      ICODE=0
      CALL UM_Submodel_Init(ICODE)
