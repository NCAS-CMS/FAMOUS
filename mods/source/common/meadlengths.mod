*ID meadlengths
*/ work around problems in calculating input and output stash lengths
*/ for use with ocean MEAD diagnostics
*DC ADDRLN1
*D ORH5F403.296
      ELSE IF (IGPL.EQ.41 .OR. IGPL.EQ.43) THEN
*D ORH5F403.300
      ELSE IF (IGPL.EQ.42 .OR. IGPL.EQ.44) THEN
*DC OUTPTL1
*I GPB1F402.516
        ! for data on ocean zonal grid, subdomain returned by
        ! GLOBAL_TO_LOCAL_SUBDOMAIN does not include the halos, but needs to,
        ! so force the right answer
        IF (IGP.eq.43 .or. IGP.eq.44) THEN
          local_north=1
          local_south=g_lasize(2,mype)
        ENDIF
