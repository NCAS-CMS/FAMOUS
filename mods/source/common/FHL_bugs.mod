*DECLARE BL_CTL1
*D ACN1F405.42
            IF (ABS(D1(J_CO2_EMITS+I-1)) <=  ABS(0.999*RMDI)) THEN
*D ACN1F405.48
            IF (ABS(D1(J_CO2FLUX+I-1)) <=  ABS(0.999*RMDI)) THEN
*D ACN1F405.53
            IF (ABS(LAND_CO2(I)) <=  ABS(0.999*RMDI)) THEN
*DECLARE TRANO2A1
*D TRANO2A1.378 
        IF (ABS(UOUT(I,J)) <=  ABS(0.999*RMDI)) THEN
*D CCN1F405.337
            IF (ABS(CO2FLUXOUT(i,j))  <=  ABS(0.999*RMDI)) THEN
*DECLARE MICROB7A
*B MICROB7A.55
     &,STH_NEW                   ! WORK defines new regime for
!                                ! soil respiration dependance
!                                ! on moisture
*I MICROB7A.68
      STH_NEW = 0.0             ! FHL: NEC QUMP-like behaviour
*D MICROB7A.78,86
          IF (STH_SOIL(L) .LE. STH_WILT) THEN
            FSTH = 0.2
          ELSEIF (STH_SOIL(L) .GT. STH_WILT .AND.
     &            STH_SOIL(L) .LE. STH_NEW) THEN
            FSTH = 0.2 + 0.8 * ((STH_SOIL(L) - STH_WILT)
     &                        / (STH_NEW - STH_WILT))
          ELSEIF (STH_SOIL(L) .GT. STH_NEW .AND.
     &            STH_SOIL(L) .LE. STH_OPT) THEN
            FSTH = 1.0
          ELSEIF (STH_SOIL(L) .GT. STH_OPT) THEN
            FSTH = 1 - 0.8 * (STH_SOIL(L) - STH_OPT)
          ENDIF
