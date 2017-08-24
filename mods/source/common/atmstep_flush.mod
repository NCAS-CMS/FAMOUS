*ID FLUSH
*/
*/  Flush the buffers on unit 6 after printing the "ATMOS TIMESTEP" line.
*/  This means that the output file of a running job shows more accurately
*/  which timestep the job has reached.
*/
*DECLARE ATMSTEP1
*I ARB0F400.45
      Call flush(6)
