*/NOT GREAT - cap pco2 and ocean co2 flux calcs that can (independently)
*/blow up easily and screw up ocean TCO2, and, in coupled carbon runs,
*/everything else besides. In practice doesn't appear to affect even some
*/pretty big perturbations to a modern climate, but far from ideal
*/
*DECLARE FLUX_CO2
*B FLUX_CO2.197
      DO I = IFROM_CYC, ITO_CYC
        if (co2_flux(i).ne.co2_flux(i)) then
!          write(6,*)"f_co2: NaN!",piston_vel(I),PCO2(I)
          co2_flux(i)=0.
        endif
        if (abs(co2_flux(i)).gt.30) then
!          write(6,*)"f_co2 cap",piston_vel(I),PCO2(I)
          co2_flux(i)=sign(30.,co2_flux(i))
        end if
      ENDDO
*DECLARE PPCO2
*B PPCO2.139
      if (PCO2(I).gt.1000) then
!        write(6,*)"ppco2 blown up:",PCO2(I)
        PCO2(I)=1000.
      endif
