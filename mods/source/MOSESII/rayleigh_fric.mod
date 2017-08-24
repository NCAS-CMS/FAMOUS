*DECLARE CRUNTIMC
*// declare array, put into common block
*B CRUNTIMC.52
     & ,RAY_TAU(MAX_P_LEVELS)
*B CRUNTIMC.63
     & ,RAY_TAU
*DECLARE READLSA1
*// make part of namelist, initialise. Will be read in with namelist
*I AJX3F405.151
     &  ,RAY_TAU
*I ADR1F305.178
        RAY_TAU(LEVEL) = 0.
*DECLARE ATMDYN1
*// common block declared here already, pass to difctl as arg
*I APB0F401.78
     &          RAY_TAU,
*DECLARE DIFCTL1A
*// pick up array, use where it has valid values
*I APB0F401.1446
     &                   RAY_TAU,
*I DIFCTL1A.74
     &,RAY_TAU(P_LEVELS)
*I DIFCTL1A.337
! Hack in simple Rayleigh friction to supplment GWD in FAMOUS
          if (abs(ray_tau(k)-0.) .gt. 1e-6) then 
            U(I,K)=U(I,K)-(U(I,K)*ADVECTION_TIMESTEP
     &                           /RAY_TAU(k)/86400.)
            V(I,K)=V(I,K)-(V(I,K)*ADVECTION_TIMESTEP
     &                           /RAY_TAU(K)/86400.)
          end if
