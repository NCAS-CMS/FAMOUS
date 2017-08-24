*ID OFILTMPP
*/
*/ Replace shmem load balancing code for ocean filtering with 
*/ gcom version. Based on NEC benchmarking code.
*/
*/ Need to add code for free surface filtering.
*/
*DECLARE ARGOCFIL
*D ORH1F405.54
     &,SLAV_CNT_T,SLAV_CNT_F,MIN_K_U,MIN_K_T,MAX_K_U,MAX_K_T,MYSLAVE,
*/
*DECLARE TYPOCFIL
*I ORH1F403.2
       INTEGER MIN_K_U(JMT,0:O_NPROC-1)
     &,        MIN_K_T(JMT,0:O_NPROC-1)
     &,        MAX_K_U(JMT,0:O_NPROC-1)
     &,        MAX_K_T(JMT,0:O_NPROC-1)
       LOGICAL MYSLAVE(JMT,0:O_NPROC-1)
*/
*DECLARE TYPOCDPT
*I ORH1F405.35
     &, JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL
*/
*DECLARE COMOCDPT
*I ORH1F405.37
     &, JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL
*/
*DECLARE OCNARYPT
*I ORH1F405.30
      JP_MIKU=icount
      icount=icount+JMT*O_NPROC
      JP_MIKT=icount
      icount=icount+JMT*O_NPROC
      JP_MAKU=icount
      icount=icount+JMT*O_NPROC
      JP_MAKT=icount
      icount=icount+JMT*O_NPROC
      JP_MYSL=icount
      icount=icount+JMT*O_NPROC
*/
*DECLARE ARGOCTOP
*I ORH1F405.43
     &,O_SPCON(JP_MIKU),O_SPCON(JP_MIKT),O_SPCON(JP_MAKU)
     &,O_SPCON(JP_MAKT),O_SPCON(JP_MYSL)
*/
*DECLARE CALCFVN
*I ORH6F402.86
     &,O_NPROC
*I PXORDER.7
     &,O_NPROC           ! IN O_NPROC for ocean
*/
*DECLARE BLOKINIT
*I ORH6F402.85
     &,O_NPROC
*/
*DECLARE CTODUMP
*I ORH6F402.83
     & O_NPROC,
*I ORH6F402.84
     &,O_NPROC       ! IN O_NPROC for ocean
*/
*DECLARE OCNFRST1
*I ORH6F402.82
     & O_NPROC,
*/
*DECLARE OFILTR
*I ORH6F401.21
*CALL COCNINDX
*D ORH1F403.147
*/
*DECLARE OSETCON
*D ORH1F405.77
*IF DEF,MPP
*/
*DECLARE DECMFLTR
*I ORH1F405.97
!     03.07.03     J.Cole      Rewrite load balancing to use two-sided
!                              GCOM communications. Based on NEC 
!                              benchmarking code.
*I DECMFLTR.37
*CALL COCNINDX
*D DECMFLTR.41
*I ORH1F405.133

      INTEGER KLAST
      INTEGER IP0
      INTEGER MASTERCOUNT, TOTALWORK
      INTEGER NPROC_AV, ASSIGNEDCOUNT, IPX, ICOUNT

      INTEGER NWORKLOAD(0:O_NPROC-1)
      INTEGER SLAVELIST(0:O_NPROC-1,JMT_GLOBAL)   ! Master -> Slave1 
                                                  !        -> Slave2..
      INTEGER MYMASTER(0:O_NPROC-1,JMT_GLOBAL)    
      INTEGER ISLAVECOUNT(0:O_NPROC-1)

      LOGICAL MASTERLOAD(0:O_NPROC-1,JMT_GLOBAL)  ! TRUE: I am a master
      LOGICAL ASSIGNED(0:O_NPROC-1)

      REAL XMAX, WORKAVERAGE
      REAL SLAVECOUNT(0:O_NPROC-1)
*I ORH1F405.170


      DO J=1,JMT_GLOBAL
         DO PE=0,O_NPROC-1
            MASTERLOAD(PE,J) = .FALSE.
            SLAVELIST(PE,J) = -1
            MYMASTER(PE,J) = -1
         ENDDO
      ENDDO
      DO PE=0,O_NPROC-1
         DO J=1,JMT
            MIN_K_U(J,PE) = -1
            MAX_K_U(J,PE) = -2
            MIN_K_T(J,PE) = -1
            MAX_K_T(J,PE) = -2
            MYSLAVE(J,PE) = .FALSE.
         ENDDO
      ENDDO
*I ORH1F405.183
            NWORKLOAD(IPROC) = 0
*I ORH1F405.219
                           NWORKLOAD(IPROC) = NWORKLOAD(IPROC) + SEG_LEN
*I ORH1F405.257
                           NWORKLOAD(IPROC) = NWORKLOAD(IPROC) + SEG_LEN
*D ORH1F405.277,ORH1F405.459
         MASTERCOUNT=0
         TOTALWORK=0
         DO IPROC=0, O_NPROC - 1   ! For all PEs
            MASTERLOAD(IPROC,JIND) = NWORKLOAD(IPROC) .GT. 0
            ASSIGNED(IPROC)=.FALSE.
            IF (NWORKLOAD(IPROC) .GT. 0)  THEN
               WRITE(6,*)'JIND,IPROC,NWORKLOAD=',
     &                    JIND,IPROC,NWORKLOAD(IPROC)
               MASTERCOUNT=MASTERCOUNT+1
               TOTALWORK=TOTALWORK+NWORKLOAD(IPROC)
            ENDIF
         ENDDO
         
         IF (TOTALWORK .GT. 0) THEN
            ASSIGNEDCOUNT=0
            NPROC_AV=O_NPROC
            WORKAVERAGE=TOTALWORK/REAL(O_NPROC)
            WRITE(6,*) 'TOTALWORK,Workaverage=',TOTALWORK,Workaverage
            DO IPROC=0, O_NPROC - 1 
               SLAVECOUNT(IPROC)=0.0
               ISLAVECOUNT(IPROC)=0
               IF (MASTERLOAD(IPROC,JIND)) THEN
                  IF (NWORKLOAD(IPROC)/WORKAVERAGE .LT. 1.0) THEN
                     SLAVECOUNT(IPROC)=1.0
                     ISLAVECOUNT(IPROC)=1
                     NPROC_AV = NPROC_AV - 1
                     TOTALWORK=TOTALWORK-NWORKLOAD(IPROC)
                     WORKAVERAGE=TOTALWORK/NPROC_AV 
                     WRITE(6,*) 'TOTALWORK,Workaverage=',
     &                           TOTALWORK,Workaverage
                     ASSIGNEDCOUNT=ASSIGNEDCOUNT+1
                  ENDIF
               ENDIF
            ENDDO
            DO IPROC=0, O_NPROC - 1 
               IF (MASTERLOAD(IPROC,JIND) .AND. 
     &             ISLAVECOUNT(IPROC) .EQ. 0) THEN
                  SLAVECOUNT(IPROC)=NWORKLOAD(IPROC)/WORKAVERAGE
                  ISLAVECOUNT(IPROC)=SLAVECOUNT(IPROC)
                  ASSIGNEDCOUNT=ASSIGNEDCOUNT+ISLAVECOUNT(IPROC)
               ENDIF
            ENDDO

            IF (ASSIGNEDCOUNT .GT. O_NPROC) STOP 'ASSIGNEDCOUNT'
            DO ICOUNT=ASSIGNEDCOUNT,O_NPROC-1
               XMAX=-1.
               DO IPROC=0, O_NPROC - 1
                  IF (SLAVECOUNT(IPROC) .GT. 0) THEN
                     IF (SLAVECOUNT(IPROC)-ISLAVECOUNT(IPROC) .GT. 
     &                   XMAX) THEN
                        XMAX=SLAVECOUNT(IPROC)-ISLAVECOUNT(IPROC)
                        IPX=IPROC
                     ENDIF
                  ENDIF
               ENDDO
               IF (XMAX .LT. 0) STOP 'XMAX'
               ISLAVECOUNT(IPX)=ISLAVECOUNT(IPX)+1 
               ASSIGNEDCOUNT=ASSIGNEDCOUNT+1
               SLAVECOUNT(IPX)=0
            ENDDO

            IF (ASSIGNEDCOUNT .NE. O_NPROC) STOP 'ASSIGNEDCOUNT'

            DO IPROC=0, O_NPROC - 1
               IF (ISLAVECOUNT(IPROC) .NE. 0) THEN
                  WRITE(6,*) 'IPROC,ISLAVECOUNT,NWORKLOAD/ISLAVECOUNT=',
     &                        IPROC,ISLAVECOUNT(IPROC),
     &                        NWORKLOAD(IPROC)/ISLAVECOUNT(IPROC)
               ELSE
                  WRITE(6,*) 'IPROC,ISLAVECOUNT=',
     &                        IPROC,ISLAVECOUNT(IPROC)
               ENDIF
            ENDDO

            ICOUNT=0
            DO IPROC=0, O_NPROC - 1
               IF (ISLAVECOUNT(IPROC) .GT. 0) THEN
                  IPX = IPROC
                  ICOUNT = ISLAVECOUNT(IPROC)-1
                  DO IP0=0, O_NPROC - 1
                     IF (ICOUNT .GT. 0) THEN
                        IF ( .NOT. ( MASTERLOAD(IP0,JIND) .OR. 
     &                               ASSIGNED(IP0) )) THEN
                           ASSIGNED(IP0)=.TRUE.
                           SLAVELIST(IPX,JIND)=IP0
                           MYMASTER(IP0,JIND)=IPROC
                           IPX=IP0
                           ICOUNT=ICOUNT-1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

            DO IPROC=0, O_NPROC - 1
               IF (MASTERLOAD(IPROC,JIND)) THEN
                  WRITE(6,*) 'Master=',IPROC
                  IPX=IPROC
                  DO IP0=2,ISLAVECOUNT(IPROC)
                     WRITE(6,*) 'Slave=',SLAVELIST(IPX,JIND)
                     IPX=SLAVELIST(IPX,JIND)
                     IF (IPROC .EQ. O_MYPE) MYSLAVE(JIND,IPX) = .TRUE.
                  ENDDO
               ENDIF
            ENDDO

!           We now have decided who works for whom.
!           Now we go over all masterprocessors and assign adjacent 
!           work to the slaves which is about size of the average.

            WRITE(6,*)'JIND,SEG_CNT_U = ',JIND,SEG_CNT_U
            WRITE(6,*)'JIND,SEG_CNT_T = ',JIND,SEG_CNT_T

            ! For each PE
            DO PE = 0, O_NPROC-1

               IF (MASTERLOAD(PE,JIND)) THEN          

                  MAX_WORK = 0
                  TOT_WORK_LEFT = NWORKLOAD(PE)
                  WORKAVERAGE = NWORKLOAD(PE)/ISLAVECOUNT(PE)
                  WRITE(6,*)'JIND,PE,NWORKLOAD,WORKAVERAGE = ',
     &                       JIND,PE,NWORKLOAD(PE),WORKAVERAGE
                  IPX = SLAVELIST(PE,JIND)
                  ICOUNT = ISLAVECOUNT(PE)
                  IF (ICOUNT .EQ. 1) IPX=PE
                  KLAST = -1

                  ! For each U/V segment
                  DO L = 1, SEG_CNT_U

                     IF (WORK_PE_U(L) .EQ. PE) THEN
                        ! If this is my segment

                        IF ( WORK_K_U(L) .NE. KLAST ) THEN
!                          Change slave only when K changes
                           KLAST=WORK_K_U(L)
                           IF ( MAX_WORK .GE. WORKAVERAGE ) THEN
                              WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                                   IPX,WORKAVERAGE,MAX_WORK
                              IPX = SLAVELIST(IPX,JIND)
                              ICOUNT = ICOUNT-1
                              IF ( ICOUNT .LE. 1) IPX=PE
                              TOT_WORK_LEFT = TOT_WORK_LEFT - MAX_WORK
                              WORKAVERAGE = TOT_WORK_LEFT / 
     &                                      MAX(1,ICOUNT)
                              MAX_WORK=0
                           ENDIF 
                        ENDIF 

                        MAX_WORK = MAX_WORK + WORK_SIZE_U(L)

                        ! Assign to current slave
                        IF (IPX .EQ. O_MYPE) THEN
                           MAST_CNT_U(JIND) = MAST_CNT_U(JIND) + 1
                           MAST_PE_U(MAST_CNT_U(JIND),JIND)  = PE
                           MAST_K_U(MAST_CNT_U(JIND),JIND) = WORK_K_U(L)
                           MAST_SEG_U(MAST_CNT_U(JIND),JIND) = 
     &                        WORK_SEG_U(L)
                        ENDIF
                        IF (MIN_K_U(JIND,IPX) .LT. 0) 
     &                     MIN_K_U(JIND,IPX)=WORK_K_U(L)
                        MAX_K_U(JIND,IPX)=WORK_K_U(L)

                        IF (PE .EQ. O_MYPE) THEN
                           SLAV_CNT_U(JIND) = SLAV_CNT_U(JIND) + 1
                        ENDIF

                     ENDIF

                  ENDDO  ! Over L

                  KLAST=-1
                  DO L = 1, SEG_CNT_T

                     IF (WORK_PE_T(L) .EQ. PE) THEN
                        ! If this is my segment

                        IF ( WORK_K_T(L) .NE. KLAST ) THEN
!                          Change slave only when K changes
                           KLAST=WORK_K_T(L)
                           IF ( MAX_WORK .GE. WORKAVERAGE ) THEN
                              WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                                   IPX,WORKAVERAGE,MAX_WORK
                              IPX = SLAVELIST(IPX,JIND)
                              ICOUNT = ICOUNT -1
                              IF ( ICOUNT .LE. 1) IPX=PE
                              TOT_WORK_LEFT = TOT_WORK_LEFT - MAX_WORK
                              WORKAVERAGE = TOT_WORK_LEFT / 
     &                                      MAX(1,ICOUNT)
                              MAX_WORK=0
                           ENDIF
                        ENDIF

                        ! Assign to current slave
                        MAX_WORK = MAX_WORK + WORK_SIZE_T(L)

                        IF (IPX .EQ. O_MYPE) THEN
                           MAST_CNT_T(JIND) = MAST_CNT_T(JIND) + 1
                           MAST_PE_T(MAST_CNT_T(JIND),JIND)  = PE
                           MAST_K_T(MAST_CNT_T(JIND),JIND) = WORK_K_T(L)
                           MAST_SEG_T(MAST_CNT_T(JIND),JIND) = 
     &                        WORK_SEG_T(L)
                        ENDIF
                        IF (MIN_K_T(JIND,IPX) .LT. 0) 
     &                     MIN_K_T(JIND,IPX)=WORK_K_T(L)
                        MAX_K_T(JIND,IPX)=WORK_K_T(L)

                        IF (PE .EQ. O_MYPE) THEN
                           SLAV_CNT_T(JIND) = SLAV_CNT_T(JIND) + 1
                        ENDIF 

                     ENDIF

                  ENDDO  ! Over L
                  WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                       IPX,WORKAVERAGE,MAX_WORK

               ENDIF
            ENDDO  ! Over PE

            WRITE(6,*) 'JIND,SLAV_CNT_U,MAST_CNT_U=',
     &                  JIND,SLAV_CNT_U(JIND),MAST_CNT_U(JIND)
            WRITE(6,*) 'JIND,SLAV_CNT_T,MAST_CNT_T=',
     &                  JIND,SLAV_CNT_T(JIND),MAST_CNT_T(JIND)
            WRITE(6,*)
         ENDIF
*/
*DECLARE OFLTCN2A
*D OFLTCN2A.3
*IF -DEF,MPP
*/
*DECLARE OFLTCNTL
*D ORH1F405.522
*IF DEF,MPP
*I OFLTCNTL.45
!     03.07.03   J.Cole      Rewrite load balancing to use two-sided
!                            GCOM communications. Based on NEC 
!                            benchmarking code.
*D ORH1F405.527
*D OFLTCNTL.82,OFLTCNTL.85
       REAL FTARR(IMTIMT_FLT)
*I ORH1F405.535
     &,        ISTAT
*I ORH1F405.539
     &,     FIELD_TEMP(IMT*2,KM*2)
     &,     FIELD_FILT(IMT*2)
*D OFLTCNTL.87,OFLTCNTL.106
*IF DEF,DEBUG
       WRITE(6,*)'MAST_CNT_U(',J,') = ',MAST_CNT_U(J)
       WRITE(6,*)'SLAV_CNT_U(',J,') = ',SLAV_CNT_U(J)
       WRITE(6,*)'MAST_CNT_T(',J,') = ',MAST_CNT_T(J)
       WRITE(6,*)'SLAV_CNT_T(',J,') = ',SLAV_CNT_T(J)
*ENDIF
*D ORH1F405.562,ORH1F405.563
            ! I must set up my UA and VA fields in a temporary array
            ! ready to be send to the slaves
*D ORH1F405.570,ORH1F405.571
            ! Send U,V data to be filtered to the slave processors
            DO IPROC=0,NPROC-1
               IF (MYSLAVE(J,IPROC)) THEN
                  K=MIN_K_U(J,IPROC)
                  SIZEB=(MAX_K_U(J,IPROC)-K+1)*IMT*2
                  IF (SIZEB .GT. 0) THEN 
*IF DEF,DEBUG
                     WRITE(6,*)'SENDING U TO ',IPROC,SIZEB,
     &                                         K,MAX_K_U(J,IPROC)
*ENDIF
                     CALL GC_RSEND(100,SIZEB,IPROC,ISTAT,
     &                             FIELD_TEMP(1,K),FIELD_TEMP(1,K))
                     CALL GC_RSEND(101,1,IPROC,ISTAT,CS_FILT,CS(J))
                     CALL GC_RSEND(102,1,IPROC,ISTAT,PHI_FILT,PHI(J))
                  ENDIF
               ENDIF
            ENDDO
*I OFLTCNTL.190
         ! If I'm a slave for any processor other than myself then
         ! receive U,V data
         IF (MAST_CNT_U(J).GT.0 .AND. SLAV_CNT_U(J).EQ.0) THEN
            IPROC = MAST_PE_U(1,J)
            K=MIN_K_U(J,MYPE)
            SIZEB=(MAX_K_U(J,MYPE)-K+1)*IMT*2
*IF DEF,DEBUG
            WRITE(6,*)'RECEIVING U FROM ',IPROC,SIZEB,K,MAX_K_U(J,MYPE)
*ENDIF
            CALL GC_RRECV(100,SIZEB,IPROC,ISTAT,
     &                    FIELD_TEMP(1,K),FIELD_TEMP(1,K))
            CALL GC_RRECV(101,1,IPROC,ISTAT,CS_FILT,CS(J))
            CALL GC_RRECV(102,1,IPROC,ISTAT,PHI_FILT,PHI(J))
         ENDIF

         ! If I'm a slave for myself then I already have U,V data
         IF (MAST_CNT_U(J).GT.0 .AND. SLAV_CNT_U(J).GT.0) THEN
            CS_FILT = CS(J)
            PHI_FILT = PHI(J)
         ENDIF
*D ORH1F405.577
            ! temporary array ready to be send to the slaves.
*D ORH1F405.587,ORH1F405.588
            ! Send T,S data to be filtered to the slave processors
            DO IPROC=0,NPROC-1
               IF (MYSLAVE(J,IPROC)) THEN
                  K=MIN_K_T(J,IPROC)
                  SIZEB=(MAX_K_T(J,IPROC)-K+1)*IMT*2
                  IF (SIZEB .GT. 0) THEN 
*IF DEF,DEBUG
                     WRITE(6,*)'SENDING T TO ',IPROC,SIZEB,
     &                                         K,MAX_K_T(J,IPROC)
*ENDIF
                     CALL GC_RSEND(103,SIZEB,IPROC,ISTAT,
     &                            FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
                     CALL GC_RSEND(104,1,IPROC,ISTAT,CST_FILT,CST(J))
                     CALL GC_ISEND(105,1,IPROC,ISTAT,KMT1_FILT,KMT(1))
                  ENDIF
               ENDIF
            ENDDO
*I OFLTCNTL.209
         ! If I'm a slave for any processor other than myself then
         ! receive T,S data
         IF (MAST_CNT_T(J).GT.0 .AND. SLAV_CNT_T(J).EQ.0) THEN
            IPROC = MAST_PE_T(1,J)
            K=MIN_K_T(J,MYPE)
            SIZEB=(MAX_K_T(J,MYPE)-K+1)*IMT*2
*IF DEF,DEBUG
            WRITE(6,*)'RECEIVING T FROM ',IPROC,SIZEB,K,MAX_K_T(J,MYPE)
*ENDIF
            CALL GC_RRECV(103,SIZEB,IPROC,ISTAT,
     &                    FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
            CALL GC_RRECV(104,1,IPROC,ISTAT,CST_FILT,CST(J))
            CALL GC_IRECV(105,1,IPROC,ISTAT,KMT1_FILT,KMT(1))
         ENDIF

         ! If I'm a slave for myself then I already have T,S data
         IF (MAST_CNT_T(J).GT.0 .AND. SLAV_CNT_T(J).GT.0) THEN
            CST_FILT = CST(J)
            KMT1_FILT = KMT(J)
         ENDIF

*D ORH1F405.590
         CALL GC_GSYNC(O_NPROC,ISTAT)
*D ORH1F405.609,ORH1F405.613
*D ORH1F405.649,ORH1F405.653
*D ORH1F405.656,ORH1F405.660
*D ORH1F405.662,ORH1F405.676
           DO I=IS,IEA
               UDIF(I-ISM1,K)=-((FX*FIELD_TEMP(I,K))*SPSIN(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPCOS(I)
               VDIF(I-ISM1,K)= ((FX*FIELD_TEMP(I,K))*SPCOS(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPSIN(I)
           ENDDO

           IF (IE.GE.IMU)THEN
               DO I=2,IEB
                 UDIF(I+II,K)=((-FX*FIELD_TEMP(I,K))*SPSIN(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPCOS(I)
                 VDIF(I+II,K)=(( FX*FIELD_TEMP(I,K))*SPCOS(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPSIN(I)
               ENDDO
           ENDIF
*D ORH1F405.681,ORH1F405.688
           DO 720 I=IS,IEA
                    FIELD_TEMP(I,K)=FX*(-UDIF(I-ISM1 ,K)*SPSIN(I)
     &                  +VDIF(I-ISM1,K)*SPCOS(I))
                    FIELD_TEMP(I+IMT,K)=-UDIF(I-ISM1 ,K)*SPCOS(I)
     &                  -VDIF(I-ISM1,K)*SPSIN(I)
  720      CONTINUE

           IF(IE.GE.IMU) THEN
                    DO 722 I=2,IEB
                       FIELD_TEMP(I,K)=FX*(-UDIF(I+II,K)*SPSIN(I)
     &                    +VDIF(I+II,K)*SPCOS(I))
                       FIELD_TEMP(I+IMT,K)=-UDIF(I+II,K)*SPCOS(I)
     &                    -VDIF(I+II,K)*SPSIN(I)
  722               CONTINUE
           ENDIF
*D OFLTCNTL.324,ORH1F405.697
*D ORH1F405.716,OFLTCNTL.379
*I ORH1F405.734
           ISM1=IS-1
*D ORH1F405.742,ORH1F405.747
*D ORH1F405.750
*D ORH1F405.752,ORH1F405.753
              DO I=IS,IEA
                 FIELD_FILT(((MM-1-TLAST)*IMT)+I-ISM1) = 
     &           FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM)
              ENDDO
              IF (IE.GE.IMT) THEN
                 DO I=2,IEB
                    FIELD_FILT(((MM-1-TLAST)*IMT)+I+II) = 
     &              FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM)
                 ENDDO
              ENDIF
*D ORH1F405.758
     &           FTARR,FIELD_FILT(((MM-1-TLAST)*IMT)+1),IM,M,N,IDX)
*D ORH1F405.760,ORH1F405.768
              DO I=IS,IEA
                 FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM) = 
     &           FIELD_FILT(((MM-1-TLAST)*IMT)+I-ISM1)
              ENDDO
              IF (IE.GE.IMT) THEN
                 DO I=2,IEB
                    FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM) = 
     &              FIELD_FILT(((MM-1-TLAST)*IMT)+I+II)
                 ENDDO
              ENDIF
*D OFLTCNTL.541,OFLTCNTL.547
        ! Send filtered data from slaves to master processor
        IF (MAST_CNT_U(J) .GT. 0) THEN
           IPROC = MAST_PE_U(1,J)
           K=MAST_K_U(1,J)
           SIZEB=(MAST_K_U(MAST_CNT_U(J),J)-K+1)*IMT*2
           IF (IPROC .NE. MYPE) THEN
*IF DEF,DEBUG
              WRITE(6,*)'SENDING U TO ',IPROC,SIZEB,
     &                   K,MAST_K_U(MAST_CNT_U(J),J)
*ENDIF
              CALL GC_RSEND(106,SIZEB,IPROC,ISTAT,
     &                      FIELD_TEMP(1,K),FIELD_TEMP(1,K))
           ENDIF
        ENDIF
        IF (MAST_CNT_T(J) .GT. 0) THEN
           IPROC = MAST_PE_T(1,J)
           K=MAST_K_T(1,J)
           SIZEB=(MAST_K_T(MAST_CNT_T(J),J)-K+1)*IMT*2
           IF ( IPROC .NE. MYPE) THEN
*IF DEF,DEBUG
              WRITE(6,*)'SENDING T TO ',IPROC,SIZEB,
     &                   K,MAST_K_T(MAST_CNT_T(J),J)
*ENDIF
              CALL GC_RSEND(107,SIZEB,IPROC,ISTAT,
     &                      FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
           ENDIF
        ENDIF

        ! If master processor get filtered data from slaves
        DO IPROC=0,NPROC-1
           IF (MYSLAVE(J,IPROC) ) THEN
              K=MIN_K_U(J,IPROC)
              SIZEB=(MAX_K_U(J,IPROC)-K+1)*IMT*2
              IF (SIZEB .GT. 0) THEN
*IF DEF,DEBUG
                 WRITE(6,*)'RECEIVING U FROM ',IPROC,SIZEB,
     &                      K,MAX_K_U(J,IPROC)
*ENDIF
                 CALL GC_RRECV(106,SIZEB,IPROC,ISTAT,
     &                         FIELD_TEMP(1,K),FIELD_TEMP(1,K))
              ENDIF
           ENDIF
        ENDDO
        DO IPROC=0,NPROC-1
           IF ( MYSLAVE(J,IPROC) ) THEN
              K=MIN_K_T(J,IPROC)
              SIZEB=(MAX_K_T(J,IPROC)-K+1)*IMT*2
              IF ( SIZEB .GT. 0) THEN
*IF DEF,DEBUG
                 WRITE(6,*)'RECEIVING T FROM ',IPROC,SIZEB,
     &                      K,MAX_K_T(J,IPROC)
*ENDIF
                 CALL GC_RRECV(107,SIZEB,IPROC,ISTAT,
     &                         FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
              ENDIF
           ENDIF
        ENDDO

        ! Synchronize before further processing.
        CALL GC_GSYNC(O_NPROC,ISTAT)
*D ORH1F405.790
         CALL GC_GSYNC(O_NPROC,ISTAT)
