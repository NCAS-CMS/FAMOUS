*ID ABORTFIX
*DECLARE UMSHELL1
*/
*/  Do not call the MPP timer if there has been an error.
*/  This is because the call can hang in situations where only one PE
*/  has returned to UMSHELL because of an error code, while all the other
*/  PEs are executing a different communications call.
*/
*/  Also, flush the output buffer before doing a potential abort;
*/  much improves the usefulness of the leave file.
*/
*/  Alan Iwi -- 4/2/02
*/
*D GSM1F401.25
*IF DEF,MPP,AND,DEF,C97_3A
      if (iCode .eq. 0) then
        CALL TIMER('UM_SHELL',2)
      else
        print *,'Model aborted; timings not available'
      endif
*ELSE
      CALL TIMER('UM_SHELL',2)
*ENDIF      
*/
*B UMSHELL1.101
      call flush(6)
