*DECLARE SETFIL1A
*// simple cap to stop the stratosphere of FAMOUS_XDBUA (no top friction)
*// from exploding in very cold/spiky/low GWD climates
*//
*I SETFIL1A.147   
      LOGICAL L_CAP
*I APB7F401.214
! CAP top-level tropical easterlies that can pull the model down
      L_CAP=.FALSE.
      DO  K=P_LEVELS-3,P_LEVELS
!      DO K=1,P_LEVELS
      DO ROW=1,P_ROWS
        global_row=ROW+datastart(2)-Offy-1  ! global row number
        if (global_row.ge.18 .AND. global_row.le.22) then
          DO I=(ROW-1)*ROW_LENGTH+1+Offx,ROW*ROW_LENGTH-Offx
             if (U(I,K).lt.-125.) then
               U(I,K)=-125.
               L_CAP=.TRUE.
             endif
          ENDDO
        endif
      ENDDO
      ENDDO
      if (L_CAP) write(6,*)"capping u@-125"


