*IDENT CO2FIX
*/
*/ Ensure ocean biogechemistry sees the same OC2 levels as 
*/ the atmosphere for runs with a fully coupled carbon cycle. 
*/ 
*/ 
*DECLARE SWAPA2O2
*B SWAPA2O2.186
*CALL CRUNTIMC
*CALL UMSCALAR
*B SWAPA2O2.608
c**********************************************************************
c****PASS THE SPATIALLY CONSTANT ATMOSPHERIC pCO2 INTO THE OCEAN*******
c****(values copied across from TRANA2O)*******************************
      IF (.NOT. L_CO2_INTERACTIVE) THEN
        PCO2_ATM_0=CO2_MMR*1e6/1.5194
      END IF
c**********************************************************************

