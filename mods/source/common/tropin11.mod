*ID TROPIN11
*/
*DECLARE QCTROP1A
*D QCTROP1A.373
c set minimum tropopause pressure to the match level structure
      TRP_MAX=6000.

*DECLARE TROPIN1A
*I TROPIN1A.132
c carried over from original low-res ozone fix - 
c not actually a good idea, overestimates things at high lats.?
c     DTI=MAX_TROP_LEVEL

*B TROPIN1A.362
c with a 0.002 LAPSE_TROP TROPIN seems to overestimate tp by 1 level or so
c c.f qctrop diagnostics (although they are capped too)
c     DO I=1,L2
c       if (IT(I).gt.MIN_TROP_LEVEL) IT(I)=IT(I)-1
c     end do

*DECLARE C_LAPSE
*/ tuned from 0.002 to suit 11 level setup.
*D C_LAPSE.5
      PARAMETER(LAPSE_TROP=0.003) !  TROPOPAUSE LAPSE RATE

*DECLARE RAD_CTL1
*/
*/  Slacken criterion for highest possible tropopause so it works with
*/        the old set of 11 levels.
*/
*D AWI1F402.24
          IF ( AKH(LEVEL)/PREF+BKH(LEVEL) .LT. .07 ) MAX_TROP = LEVEL
*/
*DECLARE CONVEC3C
*/
*/  Truly evil workspace defn & use needs fudge - 38 probably excessive !
*/
*D CONVEC3C.535,536
      REAL WORK(NPNTS,38),    !  WORK SPACE
     *     WORK2(NPNTS,38)
