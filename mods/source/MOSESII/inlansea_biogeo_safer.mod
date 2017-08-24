*DECLARE ROWCALC
*I OOM1F405.90
     &,D1(joc_vflux_mask),RMDI
*DECLARE TRACER
*I OOM1F405.151
     &,vflux_mask,RMDI
*B TRACER.201
      REAL vflux_mask(IMT,JMT),RMDI
*B TRACER.1267
c do for actual inland seas, as specified by the salmask
      if (L_OCARBON) then
         DO I=1,IMT
           IF (vflux_mask(i,j).lt.0 .AND. FM(I,1).gt.0) THEN
               do n=3,nt
               do k=1,km
                 TA(I,k,n)=RMDI
               end do
               end do
               CO2_FLUX(I)=0.
               VTCO2_FLUX(I)=0.
               VALK_FLUX(I)=0.
               INVADE(I)=0.
               EVADE(I)=0.
           ENDIF
         END DO
      end if
