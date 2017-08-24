*/ This modset is for version 4.4 of the Unified Model
*/ This is a FAMOUS mod. 
*/ **********************************************************************
*/ Modification set - MEDOUT44
*/ Author - Jonathan Palmer, 9th February 1998
*/ **********************************************************************
*/
*/ UM44 checks for the 3.75*2.5 degree ocean grid when setting
*/ up the parameters for the Mediterranean outflow parametrisation.
*/ The parameters set, however, a re harwired for HadCM2; they are
*/ timestep independent because HadCM2 mixes fully across the Straits
*/ of Gibraltar at every tstep. This modset makes the low res med
*/ outlow look like the hi-res verion ie HadCM3, by making tendfrc
*/ suitably grid and timestep dependent...
*/
*IDENT MEDOUT44
*/
*DECLARE OSETCON
*D OJG1F404.26
             tendfrc=9.6e-5 * DTTS/3600
*D OJG1F404.38
! Make tendfrc suitably dependent on timestep. There is one
! low res grid box to 6 hi-res boxes. But we are only mixing
! one pair here, instead of 2 pairs in the hi-res code above.
! Hence 9.6e-5*2/6 preserves mixing volume per second
             tendfrc=9.6e-5 * 2 / 6 * DTTS/3600
