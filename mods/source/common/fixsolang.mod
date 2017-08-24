*ID FIXSOLANG
*/
*/ Fix a problem which was leading to anomalous clear-sky SW fluxes
*/ at certain points on the equator. (vn4.5) 
*/     -- Alan Iwi, 5/12/02
*/      
*DECLARE SOLANG1A
*I SOLANG1A.70
      REAL EPS
      PARAMETER (EPS=1E-3) ! small number for floating-point comparisons      
*D SOLANG1A.134,137
          IF (DIFTIM.LT.0.) THEN
            IF (ABS(DIFTIM-(DTRAD-TWOPI)).LT.EPS) THEN
              !
              ! This is the genuine case where the sun sets and rises
              ! within the timestep.
              !
              DIFSIN = DIFSIN + 2. * SQRT(1.-COSHLD**2)
              DIFTIM = DIFTIM + 2. * HLD
            ELSEIF (DIFTIM.GT.-EPS) THEN
              !
              ! This is an artefact which can arise where the timestep
              ! ends exactly at sunrise (OMEGAE.eq.0).  This should give 
              ! OMEGA1.eq.-HLD, OMEGA2.eq.-HLD, hence DIFTIM.eq.0, but
              ! under compiler optimisations can sometimes give a tiny
              ! negative value.  (Seen with Intel compiler on linux.)
              ! Set the marker to indicate no sunshine during timestep.
              !
              DIFTIM=0.
            ELSE
              !
              ! If DIFTIM is negative for any other reason, 
              ! something's awry.
              !
              CALL ABORT
     $          ('Unexpected value for DIFTIM in SOLANG; please report')
            ENDIF
          ENDIF
