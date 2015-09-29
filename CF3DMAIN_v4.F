C * * Copyright (C) 2015 Wang & Stabenau
C     This program is free software: you can redistribute it
C     and/or modify it under the terms of the GNU General
C     Public License as published by the Free Software Foundation
C     version 3.
C     
C     This program is distributed in the hope that it will
C     be usefull, but WITHOUT ANY WARRANTY; without even the implied
C     warranty of MERCHANTABILITY, or FITNESS FOR A PARTICULAR
C     PURPOSE. See the GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public
C     License with this program. If not, see <http://www.gnu.org/licenses/>
C
C *CHECK* Items that should be checked before compiling
      INCLUDE 'ALL.CMM'
      CHARACTER*4 TITLE(40)
      DIMENSION NBN(12)
      COMMON/SORTNO/K1,K2,K3
      DATA NUMREC/1/,PI,TUPI/3.141592654,6.28318531/,
     1  ICHCK/0/
C     TRYING explicit declaration of several data types 
      DATA AVGCONST1/0.0/     
      IUOUT=2
C * * For 2-D hydro select NLEV = 1. 2-d advect-diff use NLEVAD=1
      READ(1,1002)
 1002 FORMAT()
      READ(1,*)IPRINT,NMEL,NMNP,NLEV,NLEVAD,IFRIC,IDEPTH,IEDVIS,
     1  IWIND,IRDHOT,ICNVEC,KSTART,KDINCR,IDNSTY,KDIHOT,NCONSTIT,
     2  ICSALT
C * ICSALT must be read here and not in STUPAD because inlet routines
C   use it whether salt is solved for or not!  The value of ICSALT is
C   the constituent number which is used for salt.  If ICSALT = 0
C   the dynamic coupling through density is turned off.  Salt could
C   still be one of the constituents in this case.
      NLEVM1=NLEV-1
      NMNP2=NMNP*2
      CALL DFSIGM
      NMHBN=0
      NMVBN=0
      NMQBN=0
      NMQ0BN=0
      IFLUX=0
      NBANDH=0
      CALL AMATZR(PSPLUS,NMNP)
      CALL AMATZR(ETA,MXNDSL)
      CALL AMATZR(SYSFH,MXNDSL)
      CALL AMATZR(QT,MXMQ)
      CALL AMATZR(QTPREV,MXMQ)
      CALL AMATZR(AREAND,MXNDSL)
      CALL AMATZR(AVTAUX,MXNDSL)
      CALL AMATZR(AVTAUY,MXNDSL)
      NZ=2*MXNDSL
      CALL AMATZR(TAUSTX,NZ)
      CALL AMATZR(TAUSTY,NZ)
      CALL AMATZR(ETAAV,NMNP)
      CALL AMATZR(TAUBLX,MXNDSL)
      CALL AMATZR(TAUBLY,MXNDSL)
      CALL AMATZR(TAUBGX,MXNDSL)
      CALL AMATZR(TAUBGY,MXNDSL)
      CALL AMATZR(DEN1DSH,MXNDSL)
      CALL AMATZR(DEN2DSH,MXNDSL)
      CALL AMATZR(SRC3DH,MXNDSL)
      CALL AMATZR(GXX,MXNDSL)
      CALL AMATZR(GXY,MXNDSL)
      CALL AMATZR(GYY,MXNDSL)
      CALL AMATZR(TAUX,NLEV)
      CALL AMATZR(TAUY,NLEV)
      DO 10 N=1,NLEV
      CALL AMATZR(QPREV(1,N),MXMQ)
      CALL AMATZR(Q(1,N),NMNP2)
      CALL AMATZR(U(1,N),NMNP)
      CALL AMATZR(V(1,N),NMNP)
      CALL AMATZR(OMEGA(1,N),NMNP)
C * All constituents must be initialized, so no need to zero
      CALL AMATZR(DEN1(1,N),NMNP)
      CALL AMATZR(DEN2(1,N),NMNP)
      CALL AMATZR(DRHODZ(1,N),NMNP)
 10   CONTINUE
      READ(1,1004)TITLE
 1004 FORMAT(20A4)
      WRITE(IUOUT,1006)TITLE
 1006 FORMAT(20A4/20A4)
      CALL SLINE(50)
      CALL WRITE1
      CALL BCOND
      CALL SLINE(15)
      WRITE(IUOUT,1008)NMHBN,NMQBN,NMQ0BN,NMVBN,IFLUX
 1008 FORMAT(5X,'NUMBER OF PRESCRIBED BOUNDARY AND INTERNAL FLUX NOD
     1ES'/10X,'PRESCRIBED HEIGHTS, NMHBN =',I6/,10X,
     2  'PRESCRIBED LOCAL X FLUX, NMQBN =',I6/,10X,
     3  'PRESCRIBED ZERO X FLUX, NMQ0BN =',I6/,10X,
     4  'PRESCRIBED X AND Y FLUX, NMVBN =',I6/,10X,
     5  'INTERNAL FLUX NODES, IFLUX =',I6)
      CALL SLINE(36)
      CALL ELINFO
      CALL GEOM
      CALL INTDRY
      CALL INLETS
      CALL BNDSEQ
      CALL SLINE(36)
      READ(1,*)ALATT,OMEG,GRAVT,PERIOD,REFDEN,VKRMN
c 1010 FORMAT(6F10.0)
      GD6=GRAVT/6.
      GD6DEN=GRAVT/6./REFDEN
C * * DEFINE DENSITY FIELD
      CALL DFDNST
      CORIO=2.*OMEG*SIN(ALATT*PI/180.)
      PHASE=TUPI/PERIOD
      WRITE(IUOUT,1012)ALATT,CORIO,GRAVT,OMEG,PERIOD,REFDEN,VKRMN
 1012 FORMAT(5X,'SYSTEM PROPERTIES.'/10X,
     1  'AVERAGE LATITUDE, ALATT = ',F7.2,2X,'(DEGREES N)'/,
     2  10X,'CORIOLIS PARAMETER, CORIO = 2*OMEG*SIN(ALATT) = ',E10.3,
     3  2X,'(SEC-1)'/,10X,'GRAVITATIONAL ACCELERATION, GRAVT = ',F6.3
     5  , 2X,'(M/SEC2)'/,10X,'ANGULAR VELOCITY OF EARTH ROTATION, OME
     6G = ',E10.3,2X,'(SEC-1)'/,10X,
     7  'PERIOD OF HARMONIC TIDAL EXCITATION, PERIOD = ',F6.0,2X,
     8  '(SEC)'/,10X,'REFERENCE DENSITY= ',F9.2,' (KG/M3)'/,
     9  10X,'VON KARMAN, VKRMN = ',F5.2)
      CALL SLINE(36)
C * * IADVDF DETERMINES WHETHER ADVECTION-DIFFUSION
C     IS COMPUTED.  0=NO, 1=YES, -1=Initialize but no comput.
C *    ISWAD DETERMINES ADVECTION/DIFFUSION. 1=adv, 2=adv+diff., 3=diff
C * * MULTAD is the # of external time steps per adv-dif time step;
C     it must be even
C     IMULT is the time step multiplier for the internal mode.
      READ(1,*)STRTIM,ENDTIM,TINC,N0,BOUND,ISMTH,ASMTH,
     1  IADVDF,MULTAD,IMULT,ISTRT,NOUT
c 1014 FORMAT(3F10.3,I5,F5.0,I5,F5.0,2I5,3I10)
      DO 30 I=1,NMNP
      AREAND(I)=3.*DIAGH(I)
 30   CONTINUE
      CALL QMAT(NMNP,DIAGH,DIAGQ,MXNDSL,MXMQ)
      TINC2=TINC/2.
      TINCEX=TINC
      TINCIN=TINC*IMULT
      TINCAD=TINC*MULTAD
      TIME=STRTIM
      WRITE(IUOUT,1016)STRTIM,ENDTIM,TINC,N0,BOUND,ISMTH,NOUT,IMULT
 1016 FORMAT(5X,'INTEGRATION PARAMETERS,'/10X,
     1  'START TIME OF INTEGRATION, STRTIM = ',F13.1,2X,'SEC'/
     2  ,10X,'END TIME OF INTEGRATION, ENDTIM = ',F13.1,2X,'SEC'/
     3  ,10X,'CONSTANT TIME INCREMENT, TINC = ',F7.1,2X,'SEC'/,10X
     4  ,'EXTERNAL NODE AT WHICH VARIATION IS BOUNDED BY BOUND, N0 = ',
     5  I4/,10X,'CRUDE STABILITY CONTROL,BOUND = ',F6.2/,10X,
     6  'SMOOTHING AT EVERY ',I12,' TIME STEPS'/,10X,
     7  'OUTPUT WILL BE PRINTED FOR EVERY ',I12,' TIMESTEPS'/,10X,
     8  'MULTIPLIER FOR INTERNAL TIME STEPPING, IMULT =',I12)
      WRITE(IUOUT,1018)IADVDF,MULTAD
 1018 FORMAT(5X,'DO ADVECTION-DIFFUSION, IADVDF = ',I3/,5X
     1    ,'A-D MULTIPLIER, MULTAD = ',I12)
C * * Note, you can make a diagnostic calculation using prescribed
C     densities, by using IADVDF=1 and a very large MULTAD
      ENDTIM=ENDTIM+0.001
      N0=NINT(N0)
      CALL SLINE(36)
      CALL OCNSGM
      CALL SLINE(36)
      READ(1,1020)NMBEL
 1020 FORMAT(16I5)
      NMBEL3=3*NMBEL
      READ(1,1022)(NBN(I),I=1,NMBEL3)
 1022 FORMAT(16I4)
      WRITE(IUOUT,1024)NMBEL
 1024 FORMAT(10X,'NUMBER OF BOUNDARY ELEMENTS, NMBEL=',I8)
      WRITE(IUOUT,1026)(NBN(I),I=1,NMBEL3)
 1026 FORMAT(10(I5,':',I4,'-',I4))
      DO 40 I=1,NMBEL
      I3=3*I
      NE=NBN(I3-2)
      NBN(I3-2)=NILM(NE)
      NE=NBN(I3-1)
      NBN(I3-1)=NINT(NE)
      NE=NBN(I3)
      NBN(I3)=NINT(NE)
 40   CONTINUE
      WRITE(IUOUT,1026)(NBN(I), I=1,NMBEL3)
      READ(1,1028)WINDSP,WDIRC,WTRAMP
 1028 FORMAT(3F10.0)
      READ(1,1030)DP,VF,R,ALPHA,XS,YS
 1030 FORMAT(6F10.0)
      WRITE(IUOUT,1032)WINDSP,WDIRC,WTRAMP
 1032 FORMAT(5X,'WIND SPEED (M/SEC) =',F10.3/,
     1  5X,'WIND DIRECTION (COMING FROM REL TO X-AXIS) =',F10.3,
     2 ' WTRAMP',F10.3)
      WDIRC=WDIRC/180.*PI+PI
      WRITE(IUOUT,1034)DP,VF,R,ALPHA,XS,YS
 1034 FORMAT('   STORM PARAMETERS    '/,
     1  '       DP = ',F10.3,' IN HG.'/,
     2  '       VF = ',F10.2,' KNOTS'/,
     3  '       R  = ',F10.1,' NAUTICAL MILES'/,
     4  '       ALPHA = ',F10.2,' DEGREES'/,
     5  '       XS = ',F10.2,' M'/,
     6  '       YS = ',F10.2,' M')
C     IF(R .GT. 0.)CALL STORM(WSPEED,DP,VF,R)
      READ(1,1020)IUNITQ,IUNITH,IURDHT,IUWTHT,IU3DQ
      WRITE(IUOUT,1036)IUNITQ,IUNITH,IURDHT,IUWTHT,IU3DQ
 1036 FORMAT(5X,'MASS STORAGE FILES FOR Q AND H: IUNITQ = ',I3,
     1  '  IUNITH = ',I3/,
     2  5X,'MASS STORAGE FILES FOR HOT START: IURDHT = ',I3,
     3  ' IUWTHT = ',I3/'  IU3DQ = ',I3)
       IF(IUNITQ*IUNITH*IURDHT*IUWTHT*IU3DQ .EQ. 0)THEN
         WRITE(*,*)'THERE IS A ZERO UNIT NUMBER',
     1     IUNITQ,IUNITH,IURDHT,IUWTHT,IU3DQ
         STOP
       ENDIF
       READ(1,*)NBTFRC,NODETA,NWNDST,NODCONSTITOUT,ICONSTITOUT
C      NMNP24=4*(NMNP2+NNLTS)
C      NMNP4=4*NMNP
C      OPEN(IUNITQ,FILE=' ',ACCESS='DIRECT',STATUS='UNKNOWN',
C     1  RECL=NMNP24)
C      OPEN(IUNITH,FILE=' ',ACCESS='DIRECT',STATUS='UNKNOWN',
C     1  RECL=NMNP4)
C * * Open IUWTHT for writing hot start variables, blocksize is in bytes.
      JBLCKSZ=4*NMNP
C     Skipping the use of a blocksize argument, let compiler optimize
C     Actually, changed from BLOCKSIZE to RECL and make variable start with J for integer
      OPEN(IUWTHT,STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=JBLCKSZ)
C      OPEN(IUWTHT,STATUS='UNKNOWN',FORM='UNFORMATTED')
C      Added *2 then tried +1 to NRECQ to fix "Write exceeds length of DIRECT access record" 
C      On John Wang's advice using 4*NMNP2 for NRECQ but this assumes
C      that the RECL is being used in words
C      NRECQ=(NMNP2+NNLTS)*2
C      NRECQ=(NMNP2+NNLTS)+1
C      NRECQ=NMNP2+NNLTS
      NRECQ=NMNP2*4
C     Added NRECH to allow RECL as an argument instead of BLOCKSIZE
C     in the following OPEN statements, using 4* as above...
      NRECH=NMNP*4
      OPEN(IUNITQ,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=NRECQ)
      OPEN(IUNITH,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=NRECH)
      IF(NLEV .GT. 1)THEN
        OPEN(IU3DQ,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=NMNP2)
        NMRC3D=1
      ENDIF
C * * SET INITIAL CONDITIONS
      IF(IRDHOT .EQ. 1)CALL INCOND
C * Allow using density for diagnostic comp. Use IADVDF<0 to set
C   initial density field without further solution of advect-diff.
      IF(IADVDF .NE. 0)THEN
        CALL STUPAD
        IF(MULTAD/2*2 .NE. MULTAD)THEN
          WRITE(*,*)'MULTAD NOT EVEN'
          STOP
        ENDIF
        ISAVEU=MULTAD/2
C * * ICNTAD counts number of external steps since last adv-dif step.
        ICNTAD=0
        DO 50 LEV=1,NLEVAD
          DO 60 I=1,NMNP
            IF(NLEVAD .EQ. 1)THEN
              UTMDT(I,LEV)=UT(I)
              VTMDT(I,LEV)=VT(I)
            ELSE
              UTMDT(I,LEV)=U(I,LEV)
              VTMDT(I,LEV)=V(I,LEV)
              WTMDT(I,LEV)=W(I,LEV)
            ENDIF
 60       CONTINUE
 50     CONTINUE
      ENDIF
      TIME=TIME+TINC2
      ITIME=0
      IF(IRDHOT .EQ. 0)THEN
        IF(NMHBN .GT. 0)CALL STETAB
        CALL CETA2D(ETA,DEPTH,HT,NMNP)
      ENDIF
C
C * * BEGIN EXTERNAL MODE.  ISWAD = 3 DO ONLY DIFFUSION.
C
 70   ITIME=ITIME+1
C * * Set evaporation so that it is available for advection and or diffusion
      CALL STEVAP
      IF(ISWAD .EQ. 3)THEN
C * Time is off by TINC2, see about 10 lines above here.
        TIME=TIME+TINCAD
        IF(TIME .GT. ENDTIM)STOP
        GOTO 670
      ENDIF
      IF(IFRIC .EQ. 0)CALL NTPTAU
 1038 FORMAT(I8,2e11.3,F7.3,2e11.3,2F7.3)
      TIME=TIME+TINC2
C call moved to another place      CALL CHKDRY
      CALL INTIME
      IF(R .LE. 0.)GO TO 90
C     CALL MVHUR(XS,YS,ALPHA,VF,TINC)
      WRITE(IUOUT,1040)XS,YS,ALPHA,VF,TIME
 1040 FORMAT('  THE HURRICANE CENTER IS AT (',
     1  F8.0,',',F8.0,')'/,
     2  '  DIRECTION OF TRAVEL, ALPHA = ',F6.1,
     3  '  FORWARD SPEED, VF =',F8.2,' KNOTS, @ TIME '
     4  ,F8.0,' SEC')
C *  DEFINE THE WIND STRESS/DENST
 90   IF(IWIND .EQ. 3)THEN
C        CALL WNDFLD(SPEED)
        CONTINUE
      ELSEIF(IWIND .EQ. 1)THEN
        CALL WNDGLO
      ELSEIF(IWIND .EQ. 2)THEN
        CALL VARWND
      ENDIF
C *  DEFINE BOTTOM STRESS/DENST. AT MOMENT MANNING IS USED IN FRC23Q
      CALL FRC23Q
C     CALL BOUNDE(QT,EDXX,EDYY,EDXY,XORD,YORD,SYSFQ,
C     1  NMBEL,NBN,NMBEL3)
      IF(NMQBN .GT. 0)CALL LOCGLO(QTBANG,NQN,SYSFQ,1.,NMQBN)
      IF(NMQ0BN .GT. 0)CALL LOCGLO(QT0ANG,NQ0N,SYSFQ,1.,NMQ0BN)
      CALL SOLVX(QT,SYSFQ,QTPREV,DIAGQ,NMNP2)
      IF(NMRAD .GT. 0)CALL RADIAT
      CALL STQB(QT,QTB,NQN,NVN)
      IF(NMQ0BN .GT. 0)CALL STQ0B(QT,NQ0N,NMQ0BN)
      IF(NMBDRY .GT. 0)CALL STQDRY
      DO 100 I=1,NMNP2
      QTPREV(I)=QT(I)
 100  CONTINUE
C * * Rotate to global coordinate system and determine UT, VT
      IF(NMQBN .GT. 0)CALL LOCGLO(QTBANG,NQN,QT,-1.,NMQBN)
      IF(NMQ0BN .GT. 0)CALL LOCGLO(QT0ANG,NQ0N,QT,-1.,NMQ0BN)
      CALL VELT
C * * Calculate ocean boundary discharges
C      CALL OCDSCH
C * * Advection-diffusion counter for assigning current
      ICNTAD=ICNTAD+1
C *   Need imult to be .5* multad, or need to
C     use half interpol of velocity. Latter not implemented yet.
C     There is a similar block for NLEVAD > 1 at end of internal sol.
      IF(ICNTAD .EQ. ISAVEU .AND. NLEVAD .EQ. 1)THEN
        DO 110 I=1,NMNP
          UTMHDT(I,NLEVAD)=UT(I)
          VTMHDT(I,NLEVAD)=VT(I)
 110    CONTINUE
      ENDIF
      IF(ITIME .LT. KSTART)GO TO 120
      IF(ITIME/KDINCR*KDINCR .NE. ITIME)GO TO 120
      WRITE(IUNITQ,REC=NUMREC)(QT(I), I=1,NMNP2),(QNLT12(I), I=1,NNLTS)
      WRITE(IUNITH,REC=NUMREC)(HT(I), I=1,NMNP)
      NUMREC=NUMREC+1
 120  CONTINUE
C * *  Allow print out to start after istrt time steps
      IF(ITIME. GE. ISTRT)THEN
        CALL OUTP2D(IOTP2D)
        IF(IOTP2D .GT. 0)THEN
          CALL WRTEXT
        ENDIF
      ENDIF
      CALL CHCK2D(HT(N0),DEPTH(N0),BOUND,ICHCK)
      IF(ICHCK .EQ. 1)THEN
        CALL WRTEXT
        STOP
      ENDIF
      IF(TIME+TINC .GT. ENDTIM)THEN
        IF(ILOOK*IADVDF .GT. 0)THEN
          REWIND(8)
          DO 340 N=1,NMNP
            WRITE(8,1054)N,LOOKUP(N,1),(LOOKUP(N,J), J=2,LOOKUP(N,1)+1)
 1054       FORMAT(I6,15I5)
 340      CONTINUE
        ENDIF
        STOP
      ENDIF
      TIME=TIME+TINC2
      CALL AMATZR(SYSFH,NMNP)
C * * change from qt to ubar/vbar, need to have u and v vary
C       linearly for compatibility with 3-D model.
      CALL FORCEH
      CALL CNLET
      IF(IFLUX .GT. 0)THEN
C * Source and sinks specified in input
        DO 150 I=1,IFLUX
          J=NFLUX(I)
          SYSFH(J)=SYSFH(J)+FLUX(I)
 150    CONTINUE
      ENDIF
      IF(IOPDRY .EQ. 1)THEN
C * * For hot start need to compute volsrc
        DO 20 I=1,NMNP
          SYSFH(I)=SYSFH(I)+VOLSRC(I)
 20     CONTINUE
      ENDIF
C * * Include canal structure discharges
      CALL STRCFL
      CALL GNSWAT
      CALL SOLVX(ETA,SYSFH,ETAPRV,DIAGH,NMNP)
C * * Smoothe
      IF(ITIME/ISMTH*ISMTH .EQ. ITIME)THEN
        IF(NMHBN .GT. 0)CALL STETAB
        CALL SMTH
      ENDIF
      IF(NMHBN .GT. 0)CALL STETAB
      CALL CETA2D(ETA,DEPTH,HT,NMNP)
      CALL CHKDRY
C * * Completion of one step of external mode.  If IMULT is set
C     do internal mode else go to advection-diffusion section
      IF(ITIME/IMULT*IMULT .NE. ITIME)GO TO 160
C
C * * Begin internal mode solution
C
      DO 170 I=1,NMNP
        ETAAV(I)=ETAAV(I)/IMULT
 170  CONTINUE
      CALL CETA2D(ETAAV,DEPTH,HTAV,NMNP)
C * * D ETA/ DT is found wrt external time step.
      CALL DETADT
      TINC=TINCIN
C * * Horizontal momentum solver
      CALL UV3DH
C * * qprev at prev int step. q at new step (n+*).
C       (q is really u). Both local.
C    Need to pass u/v in local cords&prev time step to vedstr. qprev has them.
C    Need u/v in local cords @ new time step (n+*) to mkrhs, q has them.
C    wndloc along with next two statments rotate global wind to local cord.
      CALL WNDLOC
      CALL LGTAUB(QBANG,NQN,TAUWLX,TAUWLY,1.,NMQBN)
      CALL LGTAUB(Q0BANG,NQ0N,TAUWLX,TAUWLY,1.,NMQ0BN)
C***  if btm strs is from ext mode, use the average after rotating it to local coords.
      IF(IFRIC .GT. 0)THEN
        DO 190 I=1,NMNP
         TAUBLX(I)=AVTAUX(I)/IMULT
         TAUBLY(I)=AVTAUY(I)/IMULT
 190    CONTINUE
        CALL LGTAUB(QBANG,NQN,TAUBLX,TAUBLY,1.,NMQBN)
        CALL LGTAUB(Q0BANG,NQ0N,TAUBLX,TAUBLY,1.,NMQ0BN)
        CALL AMATZR(AVTAUX,MXNDSL)
        CALL AMATZR(AVTAUY,MXNDSL)
      ENDIF
C * * The profile is calculated in local coordinates, since q and qprev
C       xfrc and yfrc are in local.
      CALL VEDSTR
c      IF(NMQBN. GT. 0)THEN
        DO 240 LEV=2,NLEV
          CALL STQB(Q(1,LEV),QB(1,LEV),NQN,NVN)
 240     CONTINUE
c      ENDIF
      IF(NMQ0BN .GT. 0)THEN
        DO 260 LEV=2,NLEV
          CALL STQ0B(Q(1,LEV),NQ0N,NMQ0BN)
 260     CONTINUE
      ENDIF
      IF(NMBDRY .GT. 0)CALL STQ3DR
C * * Adjust internal q's to external q in local system. Q is really u
      CALL QTALLY
C * * CALCULATE NEW BOTTOM STRESSES FOR EXTERNAL MODE
      IF(IFRIC .EQ. 0)THEN
        CALL BSTRSS
C***  if btm strs is from int mode, rotate it to global coords for
C***  frc23q which is in global. Use TAUSTX/Y. Keep TAUBLX/Y in local.
        CALL LGTAUB(QBANG,NQN,TAUSTX(1,2),TAUSTY(1,2),-1.,NMQBN)
        CALL LGTAUB(Q0BANG,NQ0N,TAUSTX(1,2),TAUSTY(1,2),-1.,NMQ0BN)
C        IF(NDTAU .GT. 0)THEN
C      WRITE(IUOUT,1052)ITIME,TAUSTX(NDTAU,2),TAUSTY(NDTAU,2)
C 1052 FORMAT(I8,7E10.3)
C        ENDIF
      ENDIF
      DO 270 LEV=2,NLEV
      DO 280 I=1,NMNP2
        QPREV(I,LEV)=Q(I,LEV)
 280  CONTINUE
 270  CONTINUE
C * * Rotate Q to global. QPREV remains in local.
      DO 290 LEV=2,NLEV
      IF(NMQBN .GT. 0)CALL LOCGLO(QBANG(1,LEV),NQN,Q(1,LEV),-1.,NMQBN)
      IF(NMQ0BN .GT. 0)CALL LOCGLO(Q0BANG(1,LEV),NQ0N,Q(1,LEV),-1.,
     1  NMQ0BN)
      CALL VEL(Q(1,LEV),U(1,LEV),V(1,LEV))
C * Save u(t-.5*dt) for advection solution
      IF(IADVDF .GT. 0 .AND. NLEVAD .GT. 1 .AND. ICNTAD .EQ. ISAVEU)THEN
        DO 295 I=1,NMNP
          UTMHDT(I,LEV)=U(I,LEV)
          VTMHDT(I,LEV)=V(I,LEV)
 295    CONTINUE
      ENDIF
 290  CONTINUE
C *041701* Compute deviation contributions to advective accl.
C          Global coord and at just calc internal mode step.
      IF(ICNVEC .EQ. 2)CALL VFEDBK(Q,QT,HT,DSIGMA,
     1    GXX,GXY,GYY,NMNP,NLEV)
C * * icnvec=2 means full 3-D of convective acc
C * * Find the average source flux over the IMULT external steps
      DO 350 I=1,NMNP
        SRC3DH(I)=SRC3DH(I)/IMULT
 350  CONTINUE
      NLEVM2=NLEV-2
      DO 300 LEV=1,NLEVM2
        CALL UVAVG
        CALL SLVOMG
C *  Omega is the vertical velocity relative to the sigma level
C    So define W as omega/H, ie. vertical velocity in sigma space relative
C      to sigma level.
        DO 310 I=1,NMNP
          W(I,LEV+1)=OMEGA(I,LEV+1)/HTAV(I)
C * *   W( ,1) and W( ,NLEV) are always zero.
 310    CONTINUE
        IF(IADVDF.GT.0 .AND. NLEVAD.GT.1 .AND. ICNTAD.EQ.ISAVEU)THEN
          DO 360 I=1,NMNP
            WTMHDT(I,LEV+1)=W(I,LEV+1)
 360      CONTINUE
        ENDIF
 300  CONTINUE
      CALL AMATZR(SRC3DH,MXNDSL)
C * * Writing q to disk
      IF(ITIME .GE. KSTART)THEN
        IF(ITIME/KDINCR*KDINCR .EQ. ITIME)THEN
          DO 330 LEV=2,NLEV
            WRITE(IU3DQ,REC=NMRC3D)(Q(I,LEV), I=1,NMNP2)
            NMRC3D=NMRC3D+1
 330      CONTINUE
        ENDIF
      ENDIF
C * *  zero etaav before starting averaging etas for next internal mode
      CALL AMATZR(ETAAV,NMNP)
 160  TINC=TINCEX
C * * End of internal mode hydro solution
C * * Begin advection-diffusion solution, so this is only done
C * * for every MULTAD*TINCEX
C * If ISWAD=3 then advd is done every time tincad (tinc*multad) step
C      (ICNTAD is not used)
 670  IF((IADVDF .GT. 0) .AND. (ICNTAD .EQ. MULTAD .OR. ISWAD .EQ. 3))
     1  THEN
        CALL DOADVD
        IF(NLEV .GT. 1)CALL TURBCL
      ELSEIF(ITIME/IMULT*IMULT .EQ. ITIME)THEN
        IF(NLEV .GT. 1)CALL TURBCL
      ENDIF
      CALL OUTP3D
C * Get turbulence parameters for next time step
C * * Save variables for hot start
      IF(ITIME/KDIHOT*KDIHOT .EQ. ITIME)THEN
        REWIND(IUWTHT)
C        IF(NMHBN .GT. 0)CALL STETAB
        CALL HOTCND
        WRITE(IUOUT,1050)NMNP,ITIME,ETA(1),ETA(NMNP)
 1050 FORMAT('  VARIABLES STORED FOR HOT START, NMNP,ITIME,ETA1,',
     1  'ETANMNP  ',I4,I8,2F8.5)
      ENDIF
      GO TO 70
      END
      SUBROUTINE AMATZR(AMAT,N)
      DIMENSION AMAT(N)
      DO 10 I=1,N
      AMAT(I)=0.
 10   CONTINUE
      RETURN
      END
      SUBROUTINE BCOND
      INCLUDE 'ALL.CMM'
      DIMENSION DUM(32)
      DO 10 I=1,NMNP
      READ(1,1002)NEXT(I),NBC(I),XORD(I),YORD(I),DEPTH(I),
     1  (DUM(J), J=1,4)
 1002 FORMAT(2I5,3F10.0,2F10.1,2F10.0)
C * * Since internal and external veloc. must have same rotation
C     because of vertical stress treatment, just set ang equal
C      IF(NBC(I) .GT. 0 .AND. NLEV .GT. 1)READ(1,1002)NEXT(I),NBC(I),
C     1  (DUM(J), J=5,11)
C      IF(NBC(I) .GT. 0 .AND. NLEV .GT. 3)READ(1,1002)NEXT(I),
C     1  NBC(I),(DUM(J), J=12,18)
C      IF(NBC(I) .GT. 0 .AND. NLEV .GT.7)READ(1,1002)NEXT(I),
C     1   NBC(I),(DUM(J), J=19,25)
C      IF(NBC(I) .GT. 0 .AND. NLEV .GT.10)READ(1,1002)NEXT(I),
C     1   NBC(I),(DUM(J), J=26,32)
      DEPTH(I)=-DEPTH(I)
      HT(I)=-DEPTH(I)
      IUNIV=NEXT(I)
      NINT(IUNIV)=I
      I1=NBC(I)
      IGOTO=I1+1
      NX=NEXT(I)
      IF(IDEPTH .EQ. 2) GO TO 20
      DEPTH(I)=DEPTH(1)
      HT(I)=-DEPTH(1)
 20   GO TO (30,40,50,60,70,80,90,30,30,100,110,120,130,
     1  30,140,30,30,150,160,170,180),IGOTO
 80   NMHBN=NMHBN+1
      NMQBN=NMQBN+1
      NMVBN=NMVBN+1
      NHN(NMHBN)=I
      NQN(NMQBN)=I
      NVN(NMVBN)=I
      HB(NMHBN)=DUM(1)
      ALAG(NMHBN)=DUM(2)
      QTB(NMQBN)=DUM(3)
      QTBANG(NMQBN)=DUM(4)
C      JD=4
      DO 190 LEV=1,NLEV
C      JDD=JD+2*LEV-1
C      QB(NMQBN,LEV)=DUM(JDD)
      QB(NMQBN,LEV)=DUM(3)
      QBANG(NMQBN,LEV)=DUM(4)
 190  CONTINUE
      GO TO 30
 70   NMQBN=NMQBN+1
      NMVBN=NMVBN+1
      NQN(NMQBN)=I
      NVN(NMVBN)=I
      QTB(NMQBN)=DUM(1)
      QTBANG(NMQBN)=DUM(2)
      DO 200 LEV=1,NLEV
      QB(NMQBN,LEV)=DUM(1)
      QBANG(NMQBN,LEV)=DUM(2)
 200  CONTINUE
      GO TO 30
 60   NMHBN=NMHBN+1
      NMQBN=NMQBN+1
      NHN(NMHBN)=I
      NQN(NMQBN)=I
      HB(NMHBN)=DUM(1)
      ALAG(NMHBN)=DUM(2)
      QTB(NMQBN)=DUM(3)
      QTBANG(NMQBN)=DUM(4)
      DO 210 LEV=1,NLEV
      QB(NMQBN,LEV)=DUM(3)
      QBANG(NMQBN,LEV)=DUM(4)
 210  CONTINUE
      GO TO 30
 50   NMHBN=NMHBN+1
      NHN(NMHBN)=I
      HB(NMHBN)=DUM(1)
      ALAG(NMHBN)=DUM(2)
      GO TO 30
 40   NMQBN=NMQBN+1
      NQN(NMQBN)=I
      QTB(NMQBN)=DUM(1)
      QTBANG(NMQBN)=DUM(2)
      DO 220 LEV=1,NLEV
      QB(NMQBN,LEV)=DUM(1)
      QBANG(NMQBN,LEV)=DUM(2)
 220  CONTINUE
      GO TO 30
 90   IFLUX=IFLUX+1
      NFLUX(IFLUX)=I
      FLUX(IFLUX)=DUM(1)
      GO TO 30
 100  NMQ0BN=NMQ0BN+1
      NQ0N(NMQ0BN)=I
      QT0ANG(NMQ0BN)=DUM(2)
      DO 230 LEV=1,NLEV
      Q0BANG(NMQ0BN,LEV)=DUM(2)
 230  CONTINUE
      GO TO 30
 110  NMQ0BN=NMQ0BN+1
      NMVBN=NMVBN+1
      NQ0N(NMQ0BN)=I
      NVN(NMVBN)=I
      QT0ANG(NMQ0BN)=DUM(2)
      DO 450 LEV=1,NLEV
        Q0BANG(NMQ0BN,LEV)=DUM(2)
 450  CONTINUE
      GO TO 30
 120  NMQBN=NMQBN+1
      NMRAD=NMRAD+1
      NQN(NMQBN)=I
      NRAD(NMRAD)=NMQBN
      QTB(NMQBN)=-999.
      QTBANG(NMQBN)=DUM(2)
      DO 240 LEV=1,NLEV
        QB(NMQBN,LEV)=-999.
        QBANG(NMQBN,LEV)=DUM(2)
 240  CONTINUE
      GO TO 30
 130  NMHBN=NMHBN+1
      NMQ0BN=NMQ0BN+1
      NMVBN=NMVBN+1
      NHN(NMHBN)=I
      NQ0N(NMQ0BN)=I
      NVN(NMVBN)=I
      HB(NMHBN)=DUM(1)
      ALAG(NMHBN)=DUM(2)
      QT0ANG(NMQ0BN)=DUM(4)
      DO 460 LEV=1,NLEV
        Q0BANG(NMQ0BN,LEV)=DUM(4)
 460  CONTINUE
      GO TO 30
 140  NMQBN=NMQBN+1
      NMVBN=NMVBN+1
      NMRAD=NMRAD+1
      NQN(NMQBN)=I
      NRAD(NMRAD)=NMQBN
      QTB(NMQBN)=-998.
      NVN(NMVBN)=I
      QTBANG(NMQBN)=DUM(2)
      DO 470 LEV=1,NLEV
        QB(NMQBN,LEV)=-998.
        QBANG(NMQBN,LEV)=DUM(2)
 470  CONTINUE
      GO TO 30
 150  IFLUX=IFLUX+1
      NMQBN=NMQBN+1
      NMRAD=NMRAD+1
      NRAD(NMRAD)=NMQBN
      QTB(NMQBN)=-999.
      NFLUX(IFLUX)=I
      NQN(NMQBN)=I
      FLUX(IFLUX)=DUM(3)
      QTBANG(NMQBN)=DUM(2)
      DO 510 LEV=1,NLEV
        QB(NMQBN,LEV)=-999.
        QBANG(NMQBN,LEV)=DUM(2)
 510  CONTINUE
      GO TO 30
 160  IFLUX=IFLUX+1
      NMQBN=NMQBN+1
      NMVBN=NMVBN+1
      NMRAD=NMRAD+1
      NRAD(NMRAD)=NMQBN
      QTB(NMQBN)=-999.
      NFLUX(IFLUX)=I
      NQN(NMQBN)=I
      NVN(NMVBN)=I
      FLUX(IFLUX)=DUM(3)
      QTBANG(NMQBN)=DUM(2)
      DO 480 LEV=1,NLEV
        QB(NMQBN,LEV)=-999.
        QBANG(NMQBN,LEV)=DUM(2)
 480  CONTINUE
      GO TO 30
 170  NMQ0BN=NMQ0BN+1
      NQ0N(NMQ0BN)=I
      IFLUX=IFLUX+1
      NFLUX(IFLUX)=I
      QT0ANG(NMQ0BN)=DUM(2)
      FLUX(IFLUX)=DUM(3)
      DO 490 LEV=1,NLEV
        Q0BANG(NMQ0BN,LEV)=DUM(2)
 490  CONTINUE
      GO TO 30
 180  NMHBN=NMHBN+1
      NMQ0BN=NMQ0BN+1
      IFLUX=IFLUX+1
      NFLUX(IFLUX)=I
      NHN(NMHBN)=I
      NQ0N(NMQ0BN)=I
      HB(NMHBN)=DUM(1)
      ALAG(NMHBN)=DUM(2)
      QT0ANG(NMQ0BN)=DUM(4)
      FLUX(IFLUX)=DUM(3)
      DO 500 LEV=1,NLEV
        Q0BANG(NMQ0BN,LEV)=DUM(4)
 500  CONTINUE
 30   IF(IPRINT.EQ.0)GOTO(10,250,10,250,250,250,10,10,10,260,
     1  10,250,260,10,250,10,10,250,250,260,260),IGOTO
      GOTO(270,280,290,300,280,300,310,270,270,320,330,280,
     1  300,270,280,270,270,340,340,350,360),IGOTO
 270  WRITE(IUOUT,1004)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I)
 1004 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3)
      GO TO 10
 280  WRITE(IUOUT,1006)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  QTB(NMQBN),QTBANG(NMQBN)
C     2  ,(QB(NMQBN,J),QBANG(NMQBN,LEV), LEV=1,NLEV)
 1006 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T56,F6.1,T64,F7.2)
      GO TO 250
 290  WRITE(IUOUT,1008)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  HB(NMHBN),ALAG(NMHBN)
 1008 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T44,F5.2,T51,F5.0)
      GO TO 10
 300  WRITE(IUOUT,1010)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),HB(NMHBN),
     1  ALAG(NMHBN),QTB(NMQBN),QTBANG(NMQBN)
C     2  ,(QB(NMQBN,LEV),QBANG(NMQBN,LEV), LEV=1,NLEV)
 1010 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T44,F5.2,T51,F5.0,T56,
     1  F6.1,T64,F7.2)
      GO TO 250
 310  WRITE(IUOUT,1012)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),FLUX(IFLUX)
 1012 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T80,F8.2)
      GO TO 10
 330  WRITE(IUOUT,1006)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  DUM(1),QT0ANG(NMQ0BN)
      GOTO 260
 340  WRITE(IUOUT,1014)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  QTB(NMQBN),QTBANG(NMQBN),FLUX(IFLUX)
 1014 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T56,F6.1,T64,F7.2,T80,F8.2)
      GOTO 250
 350  WRITE(IUOUT,1016)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  DUM(1),QT0ANG(NMQ0BN),FLUX(IFLUX)
 1016 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T56,F6.1,T64,F7.2,T80,F8.2)
      GOTO 260
 360  WRITE(IUOUT,1018)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),HB(NMHBN),
     1  ALAG(NMHBN),DUM(1),QT0ANG(NMQ0BN),FLUX(IFLUX)
 1018 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T44,F5.2,T51,F5.0,T56,
     1  F6.1,T64,F7.2,T80,F8.2)
      GO TO 260
 320  WRITE(IUOUT,1020)I,NX,XORD(I),YORD(I),DEPTH(I),NBC(I),
     1  QT0ANG(NMQ0BN)
C     2  ,(Q0BANG(NMQ0BN,LEV), LEV=1,NLEV)
 1020 FORMAT(2I6,T13,F9.1,F10.1,F8.2,T39,I3,T64,F7.2)
 260  QT0ANG(NMQ0BN)=QT0ANG(NMQ0BN)*0.01745329
      DO 370 LEV=1,NLEV
      Q0BANG(NMQ0BN,LEV)=Q0BANG(NMQ0BN,LEV)*0.01745329
 370  CONTINUE
      GO TO 10
 250  DO 380 LEV=1,NLEV
      QBANG(NMQBN,LEV)=QBANG(NMQBN,LEV)*0.01745329
 380  CONTINUE
      QTBANG(NMQBN)=QTBANG(NMQBN)*0.01745329
 10   CONTINUE
      RETURN
      END
      SUBROUTINE BNDSEQ
      INCLUDE 'ALL.CMM'
      NELRNG=-1
      IIBEG=1
      IIEND=NMEL
C * * LOOP OVER THE ELEMENTS
      NBOUND=0
      DO 50 I=1,NMEL
        DO 60 J=1,3
          J1=J+1
          IF(J .EQ. 3)J1=1
          I1=ICON(I,J)
          I2=ICON(I,J1)
C * * CHECK TO SEE IF THERE IS AN ELEMENT ON THE OPPOSITE SIDE
          IF(NELRNG .GT. 0)THEN
            IIBEG=I-NELRNG
            IIEND=I+NELRNG
            IF(IIBEG .LT. 1)IIBEG=1
            IF(IIEND .GT. NMEL)IIEND=NMEL
          ENDIF
          DO 70 II=IIBEG,IIEND
            DO 80 JJ=1,3
              IF(I1 .NE. ICON(II,JJ))GOTO 80
              JJ1=JJ+1
              IF(JJ .EQ. 3)JJ1=1
              IF(I2 .EQ. ICON(II,JJ1) .AND. II .NE. I)THEN
                WRITE(*,*)'DUPLICATE ELEMENT FOUND',I,II,I1,I2
                STOP
              ENDIF
              JJ1=JJ-1
              IF(JJ .EQ. 1)JJ1=3
              IF(I2 .EQ. ICON(II,JJ1))GOTO 60
 80         CONTINUE
 70       CONTINUE
C * * IF FALL THRU, NO ELEMENT ON OPPOSITE SIDE
          NBOUND=NBOUND+1
          IF(NBOUND .GT. MXBNDN)THEN
            WRITE(*,*)'TOO MANY BOUNDARY NODES',NBOUND,MXBNDN
            STOP
          ENDIF
          IBOUND(NBOUND,1)=I1
          IBOUND(NBOUND,2)=I2
C          IDONE(NBOUND)=0
C * * BPARMS, 1=DX1, 2=DY1, 3=S1, 4=DX2, 5=DY2, 6=S2
C * * CClockwise from node
        BPARMS(I1,4)=XORD(I2)-XORD(I1)
        BPARMS(I1,5)=YORD(I2)-YORD(I1)
        BPARMS(I1,6)=(BPARMS(I1,4)**2+BPARMS(I1,5)**2)
C * * Element number of side going cclockwise
        IELMSB(I1,2)=I
C * * Clockwise from node
        BPARMS(I2,1)=XORD(I1)-XORD(I2)
        BPARMS(I2,2)=YORD(I1)-YORD(I2)
        BPARMS(I2,3)=(BPARMS(I2,1)**2+BPARMS(I2,2)**2)
C * * Element number of side going clockwise
        IELMSB(I2,1)=I
 60     CONTINUE
 50   CONTINUE
      DO 110 I=1,NBOUND
        I1=IBOUND(I,1)
        I2=IBOUND(I,2)
 110  CONTINUE
      WRITE(2,1002)NBOUND
 1002 FORMAT('NUMBER OF BOUNDARY NODES, NBOUND = ',I7)
      RETURN
      END
      SUBROUTINE BOUNDE(Q,EDXX,EDYY,EDXY,XORD,YORD,
     1  SYSFQ,NMBEL,NBN,NMBEL3)
      INCLUDE 'PARM2.CMM'
      INCLUDE 'CGRID.CMM'
      INCLUDE 'CELEMN.CMM'
      DIMENSION Q(NMNP2),EDXX(NMEL),EDYY(NMEL),EDXY(NMEL),XORD(NMNP),
     2  YORD(NMNP),SYSFQ(NMNP2),NBN(NMBEL3)
      DO 10 I=1,NMBEL
      K=NBN(3*I-2)
      QXX=0.
      QXY=0.
      QYX=0.
      QYY=0.
      D=.25/AREA(K)
      DO 20 J=1,3
      KJ=ICON(K,J)
      QXX=B(K,J)*Q(2*KJ-1)+QXX
      QXY=A(K,J)*Q(2*KJ-1)+QXY
      QYX=B(K,J)*Q(2*KJ)+QYX
      QYY=A(K,J)*Q(2*KJ)+QYY
 20   CONTINUE
      VAR1=QXX*EDXX(K)*D
      VAR2=QXY*EDXY(K)*D
      VAR3=QYX*EDXY(K)*D
      VAR4=QYY*EDYY(K)*D
      K1=NBN(3*I-1)
      K2=NBN(3*I)
      ALY1=YORD(K2)-YORD(K1)
      ALX1=XORD(K1)-XORD(K2)
      VAR5=VAR1*ALY1+(VAR2+VAR3)*ALX1
      SYSFQ(2*K1-1)=SYSFQ(2*K1-1)+VAR5
      SYSFQ(2*K2-1)=SYSFQ(2*K2-1)+VAR5
      VAR5=(VAR2+VAR3)*ALY1+VAR4*ALX1
      SYSFQ(2*K1)=SYSFQ(2*K1)+VAR5
      SYSFQ(2*K2)=SYSFQ(2*K2)+VAR5
 10   CONTINUE
      RETURN
      END
      SUBROUTINE BSTRSS
C * * NOTE : Q IS ACTUALLY PREV and it is in local coords.
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
C * * TAUBLX&Y ARE ALWAYS IN LOCAL.  OUTSIDE OF THIS ROUTINE used for vert prof.
C     TAUBGX&Y used in external mode sol is interpolated from TAUSTX&Y in
C       NTPTAU.
C     Save the previous bottom stress
          TAUSTX(I,1)=TAUSTX(I,2)
          TAUSTY(I,1)=TAUSTY(I,2)
C * * EVALUATE NEW BTM STRESSES, IN LOCAL COORDS
C        Z=HTAV(I)*DSIGMA(1)+Z0
C  041601 Chg so that first sigma level is at Z0, also use HT rather than
C         HTAV, since this is stress at end of an internal mode step.
c        Z=HT(I)*(DSIGMA(1)+dsigma(2)+dsigma(3))
        Z=HT(I)*DSIGMA(1)
C *042001* Make z0 dependent on depth for very shallow water.
C          Limit the u/v star to 0.1* u/v
c        IF(Z/Z0 .LT. 35.)Z=Z0*35.
        ALOGZZ0=ALOG(Z/Z0)
        IF(ALOGZZ0 .LT. 3.)ALOGZZ0=3.
        USTAR=VKRMN*Q(2*I-1,2)/ALOGZZ0
        VSTAR=VKRMN*Q(2*I,2)/ALOGZZ0
        ABSVEL=SQRT(USTAR**2+VSTAR**2)
        TAUBLX(I)=USTAR*ABSVEL
        TAUBLY(I)=VSTAR*ABSVEL
        IF(IACTDR(I) .GT. 0)THEN
           TAUBLX(I)=0.
           TAUBLY(I)=0.
         ENDIF
C * * Taust( ,1) and taust( ,2) are used to calculate the global
C       bottom stress used in the external mode solution.  So they
C       are rotated to global upon exit of this routine and kept as global.
        TAUSTX(I,2)=TAUBLX(I)
        TAUSTY(I,2)=TAUBLY(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CETA2D(ETA,DEPTH,HT,NMNP)
      DIMENSION ETA(1),DEPTH(1),HT(1)
      DO 10 I=1,NMNP
      HT(I)=ETA(I)-DEPTH(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CHCK2D(H,D,BOUND,ICHCK)
      INCLUDE 'CINTEG.CMM'
      ICHCK=0
      IF(ABS(H+D) .LT. BOUND) RETURN
      WRITE(2,1002) H,D,BOUND,TIME
 1002 FORMAT(15X,'STABILITY CHECK: BOUND EXCEEDED AT NODE ',/
     1  ,5X,'HEIGHT WAS ',E11.4,3X,'DEPTH ',F8.2,5X,'BOUND = ',F6.2,
     2  5X,'TIME IS=',E13.5)
      ICHCK=1
      RETURN
      END
      SUBROUTINE DETADT
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
      ETADT(I)=(ETA(I)-ETAPRV(I))/TINCEX
  10  CONTINUE
      RETURN
      END
      SUBROUTINE DFDNST
      INCLUDE 'ALL.CMM'
      IF(IDNSTY.EQ.0)THEN
        READ(1,*)DENSTY(1,1)
        DO 10 LEV=1,NLEV
          DO 20 I=1,NMNP
            DENSTY(I,LEV)=DENSTY(1,1)
 20       CONTINUE
 10     CONTINUE
      ELSE
        IF(IDNSTY.EQ.1) THEN
          READ(1,*)(DENSTY(1,LEV), LEV=1,NLEV)
          DO 30 LEV=1,NLEV
            DO 40 I=2,NMNP
              DENSTY(I,LEV)=DENSTY(1,LEV)
 40         CONTINUE
 30       CONTINUE
        ELSE
          READ(1,*)((DENSTY(I,LEV),I=1,NMNP), LEV=1,NLEV)
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE DFSIGM
      INCLUDE 'ALL.CMM'
      READ(1,*)(SIGMA(LEV), LEV=1,NLEV)
C * * Bottom: Lev=1;  Surface: Lev=NLEV
      READ(1,*)EDYVIS,Z0,NDTAU
      WRITE(IUOUT,*)'SIGMAS ',(SIGMA(LEV), LEV=1,NLEV)
      WRITE(IUOUT,*)' VERT EDVIS =',EDYVIS,'Z0 =',Z0
      DO 10 LEV=1,NLEVM1
        DSIGMA(LEV)=SIGMA(LEV+1)-SIGMA(LEV)
 10   CONTINUE
C * * dummy dsigma at top most level. eg. needed in SOURCE for 1-layer
      DSIGMA(NLEV)=1.0
Ccc   get the coeff that appear in the d/dsigma (first deriv),due to unequal
Ccc   spacing of sigmas) and put in common fdash.
      DO 20 LEV=2,NLEVM1
        COF11=DSIGMA(LEV)+DSIGMA(LEV-1)
        COF12=DSIGMA(LEV-1)/DSIGMA(LEV)
        COF13=DSIGMA(LEV)*DSIGMA(LEV-1)
        COF14=DSIGMA(LEV)-DSIGMA(LEV-1)
        COF1(LEV)=COF14/COF13
        COF2(LEV)=COF12/COF11
        COF3(LEV)=1./COF12/COF11
 20   CONTINUE
      WRITE(IUOUT,*)'DSIGMAS ',(DSIGMA(LEV), LEV=1,NLEV)
      RETURN
      END
      SUBROUTINE DISCH(DIS,XORDN1,XORDN2,YORDN1,YORDN2,QXN1,QXN2,QYN1,
     1  QYN2)
      DEX=XORDN1-XORDN2
      DEY=YORDN2-YORDN1
      DISCX=DEY*(QXN1+QXN2)/2.
      DISCY=DEX*(QYN1+QYN2)/2.
      DIS=DISCX+DISCY
      RETURN
      END
      SUBROUTINE DNSTRM
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      NLEVM1=NLEV-1
C**   FOR ALL THE NODES GET THE VERTICALLY INTEGRATED AVERAGE OF
C**   SIGMA * D(DENSITY)/DSIGMA (=IN DEN1) AND DENSITY (=DEN2) ,
C**   INTEGRATED FROM ANY LEVEL SIGMA(NN) TO SIGMA(NLEV).
      IF(NLEV .GT. 1)THEN
        DO 20 I=1,NMNP
C**   GET SIGMA X D(DENSITY)/DSIGMA FOR EACH LEVEL FIRST.
          DO 40 LEV=2,NLEVM1
            DRHODSIG=DENSTY(I,LEV)*COF1(LEV)+
     1        DENSTY(I,LEV+1)*COF2(LEV)-DENSTY(I,LEV-1)*COF3(LEV)
            DRHODZ(I,LEV)=DRHODSIG/HT(I)
            SIGDEN(LEV)=SIGMA(LEV)*DRHODSIG
 40       CONTINUE
C * drho/dz is not needed at bottom because edzz and vdif are zero there.
C   Surface terms:
          SIGDEN(NLEV)=0.5*DSIGMA(NLEV-1)*(DENSTY(I,NLEV)-
     1             DENSTY(I,NLEVM1))/DSIGMA(NLEVM1)
          DRHODZ(I,NLEV)=(DENSTY(I,NLEV)-DENSTY(I,NLEVM1))/
     1      DSIGMA(NLEV-1)/HT(I)
C**   INTEGRATE BOTH THE TERMS, sigma drho/dsigma, and rho
          SUM1=0.0
          SUM2=0.0
          DO 50 LEV=1,NLEVM1
            N=NLEV-LEV
            SUM1=SUM1+(SIGDEN(N)+SIGDEN(N+1))*.5*DSIGMA(N)
            DEN1(I,N)=SUM1
            SUM2=SUM2+(DENSTY(I,N)+DENSTY(I,N+1))*.5*DSIGMA(N)
            DEN2(I,N)=SUM2
 50       CONTINUE
          SUM3=0.0
          SUM4=0.0
          DO 60 LEV=1,NLEVM1
            SUM3=SUM3+(DEN1(I,LEV)+DEN1(I,LEV+1))*.5*DSIGMA(LEV)
            SUM4=SUM4+(DEN2(I,LEV)+DEN2(I,LEV+1))*.5*DSIGMA(LEV)
 60       CONTINUE
          DEN1DSH(I)=SUM3
          DEN2DSH(I)=SUM4
 20     CONTINUE
      ELSE
C * * 2-D calculation, vertically const density.
        DO 70 I=1,NMNP
          DEN2DSH(I)=0.5*DENSTY(I,1)
          DEN1DSH(I)=0.
 70     CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE DRYND
C WET AND DRY
C * * SELECT NODES THAT POTENTIALLY CAN GO DRY
      INCLUDE 'ALL.CMM'
      NMBDRY=0
      DO 10 I=1,NMNP
C * * NOTE DEPTH AND HLIMIT ARE NEGATIVE.
        IACTDR(I)=0
        IF(DEPTH(I) .GT. HLIMIT)THEN
          NMBDRY=NMBDRY+1
          NDRY(NMBDRY)=I
          CALL NDELMC(I)
        ENDIF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE FORCEH
      INCLUDE 'ALL.CMM'
      COMMON/SORTNO/K(3)
      DO 10 I=1,NMEL
        VAR=0.
        DO 20 J=1,3
          K(J)=ICON(I,J)
          KJ=K(J)
          VAR=VAR+B(I,J)*QT(2*KJ-1)+A(I,J)*QT(2*KJ)
 20     CONTINUE
        VAR=-VAR/6.
C * * Evaporation (M/S), conservation of water volume
        EVAR=AREA(I)/3.
        DO 30 J=1,3
          KJ=K(J)
          SYSFH(KJ)=VAR+SYSFH(KJ)
          IF(IACTDR(KJ) .EQ. 0)SYSFH(KJ)=SYSFH(KJ)+EVAR*EVAP
 30   CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE FRC23Q
CTODO This routine could be changed to use vertically averaged flows, i.e. u
CTODO not q.  This would enhance compatibility with 3-D mode, I think.
      INCLUDE 'ALL.CMM'
C*** IFRIC=0 => BOTTOM STRSS FROM INTERNAL MODE, ALREADY IN TAUBGX&Y
C*** CALCULATED IN ROUTINE BSTRSS.  STRESS IS POSIT IN DIRECTION OF FLOW.
C    USE TAUBGX&TAUBGY AND SET CF=1.
      DIMENSION K(3),VB(3)
      IF(IFRIC .GT. 0)THEN
        DO 10 I=1,NMNP
          VABS=SQRT(UT(I)**2+VT(I)**2)
C * 2D
C          TAUBGX(I)=VABS*UT(I)
C          TAUBGY(I)=VABS*VT(I)
C * 2D end
C * * For 3D we need to have nodal bottom stresses. So CF becomes nodal
C * 3D
          TAUBGX(I)=VABS*UT(I)*CF(I)
          TAUBGY(I)=VABS*VT(I)*CF(I)
C * 3D AVTAU is zeroed in internal step
          AVTAUX(I)=AVTAUX(I)+TAUBGX(I)
          AVTAUY(I)=AVTAUY(I)+TAUBGY(I)
 10     CONTINUE
      ENDIF
      IF(ICNVEC .EQ. 0) GO TO 20
C * * icnvec=0 ignore convective acc
      DO 30 I=1,NMNP
      CXX(I)=UT(I)*QT(2*I-1)
      CYY(I)=VT(I)*QT(2*I)
      CXY(I)=UT(I)*QT(2*I)
 30   CONTINUE
 20   DO 40 I=1,NMEL
      PRESS=0.
      FXX=0.
      FXY=0.
      FYY=0.
      DRDY=0.
      DRDX=0.
      DHDXN=0.
      DHDYN=0.
      DO 50 J=1,3
        K(J)=ICON(I,J)
        KJ=K(J)
        PRESS=(ETA(KJ)*GRAVT+PSPLUS(KJ))*B(I,J)+PRESS
        FXX=FXX+UT(KJ)*B(I,J)
        FXY=FXY+(UT(KJ)*A(I,J)+VT(KJ)*B(I,J))
        FYY=FYY+VT(KJ)*A(I,J)
C2D
C        VB(J)=TAUWGX(KJ)+CORIO*QT(2*KJ)-TAUBGX(KJ)*CF(I)
C * 2D end
C3D
        VB(J)=TAUWGX(KJ)+CORIO*QT(2*KJ)-TAUBGX(KJ)
C * 3D end
C *043001*  Added the sigma d(dens)/dsigma (den1dsh) term
        DRDY=A(I,J)*DEN2DSH(KJ)+DRDY
        DRDX=B(I,J)*DEN2DSH(KJ)+DRDX
        DHDXN=DHDXN+B(I,J)*HT(KJ)
        DHDYN=DHDYN+A(I,J)*HT(KJ)
 50   CONTINUE
      PRESS=-PRESS/6.
      AD3=AREA(I)/3.
      A12=AREA(I)*12.
C * Smagorinski horizontal viscosity. EDXX etc. now should be
C   only a coefficient of order 0.2.  Deform is ds^2 * deformation
C   field. Here ds^2 is taken as 2*area.
      DEFORM=SQRT((FXX-FYY)*(FXX-FYY)+FXY*FXY)
      DO 120 J=1,3
        KJ=K(J)
c        VAR1=(EDXX(I)*FXX*B(I,J)+EDXY(I)*FXY*A(I,J))/A12
        VAR1=(EDXX(I)*FXX*B(I,J)+EDXY(I)*FXY*A(I,J))*DEFORM/A12
        SYSFQ(2*KJ-1)=SYSFQ(2*KJ-1)+AD3*VB(J)+(PRESS-VAR1)*HT(KJ)
     1   +(DHDXN*DEN1DSH(KJ)-DRDX*HT(KJ))*HT(KJ)*GD6DEN
 120  CONTINUE
      PRESS=0.
      DO 90 J=1,3
      KJ=K(J)
      PRESS=ETA(KJ)*GRAVT*A(I,J)+PRESS
C * 2D
C      VB(J)=TAUWGY(KJ)-CORIO*QT(2*KJ-1)-TAUBGY(KJ)*CF(I)
C * 2D end
C * 3D
      VB(J)=TAUWGY(KJ)-CORIO*QT(2*KJ-1)-TAUBGY(KJ)
C * 3D end
 90   CONTINUE
      PRESS=-PRESS/6.
      DO 130 J=1,3
        KJ=K(J)
c        VAR4=(EDXY(I)*FXY*B(I,J)+EDYY(I)*FYY*A(I,J))/A12
        VAR4=(EDXY(I)*FXY*B(I,J)+EDYY(I)*FYY*A(I,J))*DEFORM/A12
        SYSFQ(2*KJ)=SYSFQ(2*KJ)+(PRESS-VAR4)*HT(KJ)+AD3*VB(J)
     1   +(DHDYN*DEN1DSH(KJ)-DRDY*HT(KJ))*HT(KJ)*GD6DEN
 130  CONTINUE
C * * ICNVEC = 0 means don't do convective acc.
      IF(ICNVEC .EQ. 0) GO TO 40
      CNVCX=0.
      CNVCY=0.
      DO 80 J=1,3
        KJ=K(J)
        CNVCX=CNVCX+B(I,J)*(CXX(KJ)+GXX(KJ))+A(I,J)*(CXY(KJ)+GXY(KJ))
        CNVCY=CNVCY+A(I,J)*(CYY(KJ)+GYY(KJ))+B(I,J)*(CXY(KJ)+GXY(KJ))
 80   CONTINUE
      CNVCX=-CNVCX/6.
      CNVCY=-CNVCY/6.
      DO 110 J=1,3
        KJ=K(J)
        SYSFQ(2*KJ-1)=SYSFQ(2*KJ-1)+CNVCX
        SYSFQ(2*KJ)=SYSFQ(2*KJ)+CNVCY
 110  CONTINUE
 40   CONTINUE
      RETURN
      END
      SUBROUTINE FRC3DQ
C * * EXPLICIT VERSION
      INCLUDE 'ALL.CMM'
      DIMENSION USGMA(MXNDSL),VSGMA(MXNDSL)
C * *   FIRST DO ALL THE TERMS EXCEPT THE CONVECTIVE ACCELERATION TERMS
      DO 100 I=1,NMEL
      AD3=AREA(I)/3.
      DD2DX=0.0
      DD2DY=0.0
      FXX=0.
      FXY=0.
      FYY=0.
      DO 130 J=1,3
      KJ=ICON(I,J)
C * * DENSITY TERMS.  DEN1=INTEG(SIGMA*(PARTIAL RHO/PARTIAL SIGMA))
C     DEN2=INTEG(RHO)
      DD2DX=DD2DX+DEN2(KJ,LEV)*B(I,J)
      DD2DY=DD2DY+DEN2(KJ,LEV)*A(I,J)
      FXX=FXX+U(KJ,LEV)*B(I,J)
      FXY=FXY+(U(KJ,LEV)*A(I,J)+V(KJ,LEV)*B(I,J))
      FYY=FYY+V(KJ,LEV)*A(I,J)
 130  CONTINUE
      DD2DX=-DD2DX*GD6DEN
      DD2DY=-DD2DY*GD6DEN
      A12=AREA(I)*12.
C *  Smagorinski scheme used for horizontal viscosity. See FRC23Q
      DEFORM=SQRT((FXX-FYY)*(FXX-FYY)+FXY*FXY)
      if(DEFMAX .LT. DEFORM)DEFMAX=DEFORM
      DO 140 J=1,3
        KJ=ICON(I,J)
c        VAR1=(EDXX(I)*FXX*B(I,J)+EDXY(I)*FXY*A(I,J))/A12
        VAR1=(EDXX(I)*FXX*B(I,J)+EDXY(I)*FXY*A(I,J))*DEFORM/A12
C * * detadx really contains -deta/dx/rho*(g/6), get dhdx from uv3dh
C     HTAV is used.
      VARX=DETADX(I)+DHDX(I)*DEN1(KJ,LEV)+
     1  DD2DX*HTAV(KJ)+AD3*CORIO*U(KJ,LEV)
c      VAR4=(EDXY(I)*FXY*B(I,J)+EDYY(I)*FYY*A(I,J))/A12
      VAR4=(EDXY(I)*FXY*B(I,J)+EDYY(I)*FYY*A(I,J))*DEFORM/A12
      VARY=DETADY(I)+DHDY(I)*DEN1(KJ,LEV)+
     1  DD2DY*HTAV(KJ)-AD3*CORIO*V(KJ,LEV)
      SYSFQ(2*KJ-1)=SYSFQ(2*KJ-1)+VARX-VAR1
      SYSFQ(2*KJ)=SYSFQ(2*KJ)+VARY-VAR4
 140  CONTINUE
 100  CONTINUE
c      write(3,1002)itime,defmax
c 1002 format(i8,e12.4)
      defmax=0.
C * * IF CONVECTIVE ACCL. TERMS ARE TO BE IGNORED, RETURN.
      IF(ICNVEC.EQ.0)RETURN
C * * GET D/DSIGMA NODAL VALUES FOR U & V APPEARING IN THE CONVEC.
C * * TERMS. ASSUMED THAT MOM. EQN FOR LEV=1 IS NEVER DONE. SO FOR
C * * TOP MOST LEVEL ALONE, DO A ONE SIDED DIFFERENCE.
       IF(LEV .EQ. NLEV) THEN
         DO 300 I=1,NMNP
C           USGMA(I)=(U(I,NLEV)-U(I,NLEVM1))/DSIGMA(NLEVM1)
C           VSGMA(I)=(V(I,NLEV)-V(I,NLEVM1))/DSIGMA(NLEVM1)
C * Omega at NLEV must be zero
C           USGMA(I)=USGMA(I)*OMEGA(I,LEV)/HTAV(I)
C           VSGMA(I)=VSGMA(I)*OMEGA(I,LEV)/HTAV(I)
           USGMA(I)=0.
           VSGMA(I)=0.
 300     CONTINUE
       ELSE
         DO 310 I=1,NMNP
           USGMA(I)=U(I,LEV)*COF1(LEV)+
     1          U(I,LEV+1)*COF2(LEV)-U(I,LEV-1)*COF3(LEV)
           VSGMA(I)=V(I,LEV)*COF1(LEV)+
     1          V(I,LEV+1)*COF2(LEV)-V(I,LEV-1)*COF3(LEV)
           USGMA(I)=USGMA(I)*OMEGA(I,LEV)/HTAV(I)
           VSGMA(I)=VSGMA(I)*OMEGA(I,LEV)/HTAV(I)
 310     CONTINUE
       ENDIF
      DO 200 I=1,NMEL
      AD3=AREA(I)/3.
      DUDX=0.0
      DUDY=0.0
      DVDX=0.0
      DVDY=0.0
      DO 230 J=1,3
      KJ=ICON(I,J)
      DUDX=DUDX+U(KJ,LEV)*B(I,J)
      DVDX=DVDX+V(KJ,LEV)*B(I,J)
      DUDY=DUDY+U(KJ,LEV)*A(I,J)
      DVDY=DVDY+V(KJ,LEV)*A(I,J)
 230  CONTINUE
      DO 240 J=1,3
      KJ=ICON(I,J)
      VARX=(DUDX*U(KJ,LEV)+DUDY*V(KJ,LEV))/6.+AD3*USGMA(KJ)
      VARY=(DVDX*U(KJ,LEV)+DVDY*V(KJ,LEV))/6.+AD3*VSGMA(KJ)
      SYSFQ(2*KJ-1)=SYSFQ(2*KJ-1)-VARX
      SYSFQ(2*KJ)=SYSFQ(2*KJ)-VARY
 240  CONTINUE
 200  CONTINUE
      RETURN
      END
      SUBROUTINE GEOM
      INCLUDE 'ALL.CMM'
      COMMON/SORTNO/K(3)
      DIMENSION IPERM(3,2)
      DATA IPERM/3,1,2,2,3,1/
      IF(IPRINT .EQ. 2)WRITE(2,1002)
 1002 FORMAT(5X,'GEOMETRICAL RELATIONS'/,10X,
     2  'ELEMENT'/,10X,'NUMBER',4X,'A1',11X,'B1',
     3  11X,'A2',11X,'B2',11X,'A3',11X,'B3',10X,'AREA'/)
      DO 10 I=1,NMEL
      DO 20 J=1,3
      IUNIV=ICON(I,J)
      K(J)=NINT(IUNIV)
      ICON(I,J)=K(J)
 20   CONTINUE
      SRCHRDX(I)=0.
      SRCHRDY(I)=0.
      DO 40 J=1,3
      IPERM1=IPERM(J,1)
      IPERM2=IPERM(J,2)
      KI1=K(IPERM1)
      KI2=K(IPERM2)
      A(I,J)=XORD(KI1)-XORD(KI2)
      B(I,J)=YORD(KI2)-YORD(KI1)
C      IF(SRCHRDX(I) .LT. ABS(A(I,J)))SRCHRDX(I)=ABS(A(I,J))
C      IF(SRCHRDY(I) .LT. ABS(B(I,J)))SRCHRDY(I)=ABS(B(I,J))
      IF(SRCHRDX(I) .LT. ABS(A(I,J)))SRCHRDX(I)=5*ABS(A(I,J))
      IF(SRCHRDY(I) .LT. ABS(B(I,J)))SRCHRDY(I)=5*ABS(B(I,J))
 40   CONTINUE
      AREA(I)=0.5*(B(I,1)*A(I,2)-B(I,2)*A(I,1))
      IF(AREA(I) .GT. 0.) GO TO 30
      WRITE(2,1004) I,K(1),K(2),K(3)
 1004 FORMAT(5X,'NEGATIVE AREA IN ELEMENT:',I4,3I5)
      STOP
 30   VAR=AREA(I)/3.
      DO 50 J=1,3
        DIAGH(K(J))=DIAGH(K(J))+VAR
 50   CONTINUE
      IF(IPRINT .EQ. 2)WRITE(2,1006)I,A(I,1),B(I,1),A(I,2),B(I,2),
     1  A(I,3),B(I,3),AREA(I)
 1006 FORMAT(5X,I6,2X,7E13.6)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE HOTCND
C * * STORE CONDITIONS FOR HOT START. All writes must be in blocks of
C       NMNP because of open statement
      INCLUDE 'ALL.CMM'
      WRITE(IUWTHT)(ETAPRV(I), I=1,NMNP)
      WRITE(IUWTHT)(ETA(I), I=1,NMNP)
      WRITE(IUWTHT)(HT(I), I=1,NMNP)
      WRITE(IUWTHT)(QTPREV(I), I=1,NMNP)
      WRITE(IUWTHT)(QTPREV(I), I=NMNP+1,NMNP2)
      WRITE(IUWTHT)(QT(I), I=1,NMNP)
      WRITE(IUWTHT)(QT(I), I=NMNP+1,NMNP2)
      WRITE(IUWTHT)(UT(I), I=1,NMNP)
      WRITE(IUWTHT)(VT(I), I=1,NMNP)
      WRITE(IUWTHT)(TAUBLX(I), I=1,NMNP)
      WRITE(IUWTHT)(TAUBLY(I), I=1,NMNP)
      WRITE(IUWTHT)(TAUBGX(I), I=1,NMNP)
      WRITE(IUWTHT)(TAUBGY(I), I=1,NMNP)
      WRITE(IUWTHT)(IACTDR(I), I=1,NMNP)
      WRITE(IUWTHT)(VOLSRC(I), I=1,NMNP)
      DO 10 LEV=1,NLEV
        WRITE(IUWTHT)(Q(I,LEV), I=1,NMNP)
        WRITE(IUWTHT)(Q(I,LEV), I=NMNP+1,NMNP2)
        WRITE(IUWTHT)(QPREV(I,LEV), I=1,NMNP)
        WRITE(IUWTHT)(QPREV(I,LEV), I=NMNP+1,NMNP2)
        WRITE(IUWTHT)(OMEGA(I,LEV), I=1,NMNP)
        WRITE(IUWTHT)(W(I,LEV), I=1,NMNP)
C TODO Should save nconstit constituents, but makes it
C      difficult to change the number of constits. For now just 1
        WRITE(IUWTHT)(CONSTITP(I,LEV,1), I=1,NMNP)
        WRITE(IUWTHT)(EDZZ(I,LEV), I=1,NMNP)
        WRITE(IUWTHT)(VDIF(LEV,I), I=1,NMNP)
        WRITE(IUWTHT)(TQ(I,LEV), I=1,NMNP)
C * * For now at least don't need DENSTY because its calc based on sal
 10   CONTINUE
      IF(NLEV .GT. 1)THEN
        WRITE(IUWTHT)(TAUSTX(I,1), I=1,NMNP)
        WRITE(IUWTHT)(TAUSTX(I,2), I=1,NMNP)
        WRITE(IUWTHT)(TAUSTY(I,1), I=1,NMNP)
        WRITE(IUWTHT)(TAUSTY(I,2), I=1,NMNP)
        WRITE(IUWTHT)(DEN1DSH(I), I=1,NMNP)
        WRITE(IUWTHT)(DEN2DSH(I), I=1,NMNP)
        WRITE(IUWTHT)(GXX(I), I=1,NMNP)
        WRITE(IUWTHT)(GXY(I), I=1,NMNP)
        WRITE(IUWTHT)(GYY(I), I=1,NMNP)
      ENDIF
      WRITE(IUWTHT)TIME,ITIME,MULTAD,IMULT,NMBDRY
      RETURN
      END
      SUBROUTINE INCOND
C * * SET INITIAL CONDITIONS
      INCLUDE 'ALL.CMM'
      REAL*8 STRTIM1
      JBLCKSZ=4*NMNP
C      OPEN(IURDHT,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(IURDHT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1  RECL=JBLCKSZ)
c     DO 20 IRPT=1,2
      READ(IURDHT)(ETAPRV(I), I=1,NMNP)
      READ(IURDHT)(ETA(I), I=1,NMNP)
      READ(IURDHT)(HT(I), I=1,NMNP)
      READ(IURDHT)(QTPREV(I), I=1,NMNP)
      READ(IURDHT)(QTPREV(I), I=NMNP+1,NMNP2)
      READ(IURDHT)(QT(I), I=1,NMNP)
      READ(IURDHT)(QT(I), I=NMNP+1,NMNP2)
      READ(IURDHT)(UT(I), I=1,NMNP)
      READ(IURDHT)(VT(I), I=1,NMNP)
      READ(IURDHT)(TAUBLX(I), I=1,NMNP)
      READ(IURDHT)(TAUBLY(I), I=1,NMNP)
C * * Need TAUBL for internal mode VEDSTR and TAUBG for external
C       mode when IFRIC=0, i.e. bot. frict from internal mode solut.
      READ(IURDHT)(TAUBGX(I), I=1,NMNP)
      READ(IURDHT)(TAUBGY(I), I=1,NMNP)
      READ(IURDHT)(IACTDR(I), I=1,NMNP)
      READ(IURDHT)(VOLSRC(I), I=1,NMNP)
      DO 10 LEV=1,NLEV
        READ(IURDHT)(Q(I,LEV), I=1,NMNP)
        READ(IURDHT)(Q(I,LEV), I=NMNP+1,NMNP2)
        READ(IURDHT)(QPREV(I,LEV), I=1,NMNP)
        READ(IURDHT)(QPREV(I,LEV), I=NMNP+1,NMNP2)
        CALL VEL(Q(1,LEV),U(1,LEV),V(1,LEV))
        READ(IURDHT)(OMEGA(I,LEV), I=1,NMNP)
        READ(IURDHT)(W(I,LEV), I=1,NMNP)
        READ(IURDHT)(CONSTITP(I,LEV,1), I=1,NMNP)
        READ(IURDHT)(EDZZ(I,LEV), I=1,NMNP)
        READ(IURDHT)(VDIF(LEV,I), I=1,NMNP)
        READ(IURDHT)(TQ(I,LEV), I=1,NMNP)
C * * For now at least don't need DENSTY because its calc based on sal
 10   CONTINUE
      IF(NLEV .GT. 1)THEN
        READ(IURDHT)(TAUSTX(I,1), I=1,NMNP)
        READ(IURDHT)(TAUSTX(I,2), I=1,NMNP)
        READ(IURDHT)(TAUSTY(I,1), I=1,NMNP)
        READ(IURDHT)(TAUSTY(I,2), I=1,NMNP)
        READ(IURDHT)(DEN1DSH(I), I=1,NMNP)
        READ(IURDHT)(DEN2DSH(I), I=1,NMNP)
        READ(IURDHT)(GXX(I), I=1,NMNP)
        READ(IURDHT)(GXY(I), I=1,NMNP)
        READ(IURDHT)(GYY(I), I=1,NMNP)
      ENDIF
      READ(IURDHT)STRTIM1,ITIME1,MULTAD1,IMULT1,NMBDRY1
      WRITE(*,*)'HOT START, STORED PARAMS'
      WRITE(*,*)STRTIM1,ITIME1,MULTAD1,IMULT1,NMBDRY1
      IF(ABS(STRTIM1-STRTIM-.5*TINC) .GT. .1*TINC)THEN
        WRITE(*,*)'START TIMES DO NOT AGREE',STRTIM1,STRTIM
        WRITE(*,*)'CONTINUE ?  ANY KEY, CTRL C TO STOP'
        READ(*,*)ACHAR
      ENDIF
      IF(MULTAD1 .NE. MULTAD)THEN
        WRITE(*,*)'MULTAD CHANGED ',MULTAD1,MULTAD
        WRITE(*,*)'CONTINUE ?  ANY KEY, CTRL C TO STOP'
        READ(*,*)ACHAR
      ENDIF
      IF(IMULT1 .NE. IMULT)THEN
        WRITE(*,*)'IMULT CHANGED ',IMULT1,IMULT
        WRITE(*,*)'CONTINUE ?  ANY KEY, CTRL C TO STOP'
        READ(*,*)ACHAR
      ENDIF
      IF(NMBDRY1 .NE. NMBDRY)THEN
        WRITE(*,*)'NMBDRY CHANGED ',NMBDRY1,NMBDRY
        WRITE(*,*)'CONTINUE ?  ANY KEY, CTRL C TO STOP'
        READ(*,*)ACHAR
      ENDIF
C * * Need HT to calculate UT and VT, but UT and VT are calculated using
C     HT at previous time step.
C      CALL CETA2D(ETAPRV,DEPTH,HT,NMNP)
C      CALL VELT
C 20   CONTINUE
      RETURN
      END
      SUBROUTINE CNLET
C * * Calculate inlet flows. Valid for 1-layer and multilayer
C  Note CONSTITUENTS are calculated every external time step irrespective
C  of MULTAD and IADVDF to save boolean and averaging overhead. Only
C  a few points so not very costly, and no cost if no inlets.
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NNLTS
        NEL1=5
        ETANL1=0.
        DO 20 J=1,5
          IEL=ILNLT1(I,J)
          IF(IEL .EQ. 0)THEN
            NEL1=J-1
            GOTO 30
          ENDIF
          DO 40 K=1,3
            IN=ICON(IEL,K)
            ETANL1=ETANL1+ETAPRV(IN)
 40       CONTINUE
 20     CONTINUE
 30     ETA1AV=ETANL1/NEL1/3.
        NEL2=5
        ETANL2=0.
        DO 50 J=1,5
          IEL=ILNLT2(I,J)
          IF(IEL .EQ. 0)THEN
            NEL2=J-1
            GOTO 60
          ENDIF
          DO 70 K=1,3
            IN=ICON(IEL,K)
            ETANL2=ETANL2+ETAPRV(IN)
 70       CONTINUE
 50     CONTINUE
 60     ETA2AV=ETANL2/NEL2/3.
        DELTAE=ETA2AV-ETA1AV
        AVGE=.5*(ETA2AV+ETA1AV)
        R=AVGDPT(I)+AVGE
C * * Flow from 1 to 2.  Inhibit noise level flows.
        IF(ABS(DELTAE) .LT. 0.0005)DELTAE=0.
        QNLT12(I)=-SIGN(1.,DELTAE)*CRSRR(I)*SQRT(ABS(DELTAE)
     1    /ALNGTH(I))/R**.667/ANNLT(I)
C * * Subtract inlet flow from side 1
C       AVGCONST1=0.
       AVGCONST1=0.
        DO 80 J=1,NEL1
          IEL=ILNLT1(I,J)
          DO 90 K=1,3
            IN=ICON(IEL,K)
            SYSFH(IN)=SYSFH(IN)-QNLT12(I)/3./NEL1
C            SRC3DH(IN)=SRC3DH(IN)-QNLT12(I)/3./NEL1
            IF(NLEVAD .EQ. 1)THEN
C * Find average salinity over all elements on side 1.
              DO ICST=1,NCONSTIT
                AVGCONST1(ICST)=AVGCONST1(ICST)
     1           +CONSTITP(IN,NLEVAD,ICST)/3./NEL1
              ENDDO
            ELSE
C * Find vertical average constituent (ACONSTIT)
              ACONSTIT=0.
              DO 95 LEV=1,NLEVAD-1
                DO ICST=1,NCONSTIT
                  ACONSTIT(ICST)=ACONSTIT(ICST)+(CONSTITP(IN,LEV,ICST)
     1           +CONSTITP(IN,LEV+1,ICST))*0.5*DSIGMA(LEV)
                ENDDO
 95           CONTINUE
              DO ICST=1,NCONSTIT
                AVGCONST1(ICST)=AVGCONST1(ICST)+ACONSTIT(ICST)/3./NEL1
              ENDDO
            ENDIF
 90       CONTINUE
 80     CONTINUE
C * * Add inlet flow to side 2
        AVGCONST2=0.
        DO 100 J=1,NEL2
          IEL=ILNLT2(I,J)
          DO 110 K=1,3
            IN=ICON(IEL,K)
            SYSFH(IN)=SYSFH(IN)+QNLT12(I)/3./NEL2
C            SRC3DH(IN)=SRC3DH(IN)+QNLT12(I)/3./NEL2
            IF(NLEVAD .EQ. 1)THEN
              DO ICST=1,NCONSTIT
                AVGCONST2(ICST)=AVGCONST2(ICST)
     1            +CONSTITP(IN,NLEVAD,ICST)/3./NEL2
              ENDDO
            ELSE
C * Find vertical average sal (ASAL)
              ACONSTIT=0.
              DO 115 LEV=1,NLEVAD-1
                DO ICST=1,NCONSTIT
                  ACONSTIT(ICST)=ACONSTIT(ICST)+(CONSTITP(IN,LEV,ICST)
     1           +CONSTITP(IN,LEV+1,ICST))*0.5*DSIGMA(LEV)
                ENDDO
 115          CONTINUE
              DO ICST=1,NCONSTIT
                AVGCONST2(ICST)=AVGCONST2(ICST)+ACONSTIT(ICST)/3./NEL2
              ENDDO
            ENDIF
 110      CONTINUE
 100    CONTINUE
C * *   INSERT CORRESPONDING MASS FLUX TRANSFERS HERE
        IF(QNLT12(I) .GT. 0.)THEN
C * *   SIDE 1 IS DONOR
          DO 120 J=1,NEL2
            IEL=ILNLT2(I,J)
C * Total inflow to this element, simply assumed to be fraction
C   based on the number of elements
            FLOWIN=QNLT12(I)/NEL2
            DO 130 K=1,3
              IN=ICON(IEL,K)
              DO 140 LEV=1,NLEVAD
C * This works for both 1-layer and 3-d
                DO ICST=1,NCONSTIT
                  CONSTITP(IN,LEV,ICST)=(CONSTITP(IN,LEV,ICST)
     1             -AVGCONST1(ICST))*(1.-TINCEX*FLOWIN/HT(IN)
     2             /AREAND(IN))+AVGCONST1(ICST)
C * Flowin should be divided by 3 (3 nodes) but cancels with
C    1/3 AREAND.
                ENDDO
 140          CONTINUE
 130        CONTINUE
 120      CONTINUE
        ELSE
C * *  SIDE 2 IS DONOR, QNLT12 is negative
          DO 150 J=1,NEL1
            IEL=ILNLT1(I,J)
C * FLOWIN is flow to each element of side 1
            FLOWIN=QNLT12(I)/NEL1
            DO 170 K=1,3
              IN=ICON(IEL,K)
              DO 160 LEV=1,NLEVAD
                DO ICST=1,NCONSTIT
                  CONSTITP(IN,LEV,ICST)=(CONSTITP(IN,LEV,ICST)
     1             -AVGCONST2(ICST))*(1.+TINCEX*FLOWIN/HT(IN)
     2             /AREAND(IN))+AVGCONST2(ICST)
                ENDDO
 160          CONTINUE
 170        CONTINUE
 150      CONTINUE
        ENDIF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INLETS
      INCLUDE 'ALL.CMM'
C * * Read input for 1-D inlets.  NNLTS=# inlets
C     ALNGTH=length, CRSRR=x-sectional area, AVGDPT=avg depth,
C     ANNLT=Manning n for inlet,ILNLT1=Element numbers for inlet side 1
      READ(1,1001)NNLTS
 1001 FORMAT(I5)
      DO 10 I=1,NNLTS
        READ(1,1003)ALNGTH(I),CRSRR(I),AVGDPT(I),ANNLT(I),
     1  (ILNLT1(I,J), J=1,5),(ILNLT2(I,J), J=1,5)
 1003   FORMAT(4F10.0,10I5)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INTDRY
C WET AND DRY
      INCLUDE 'ALL.CMM'
C * * HLIMIT IS NEGATIVE
      READ(1,1001)IOPDRY,HLIMIT,HMIN,VOLFCT,HREDRY
 1001 FORMAT(I10,4F10.0)
      IF(IOPDRY .EQ. 1)THEN
        WRITE(2,1014)IOPDRY,HLIMIT,HMIN,VOLFCT,HREDRY
 1014 FORMAT('** DRYING AND WETTING OPTION CHOSEN **',I5/,
     1  'NODES WITH DEPTH > HLIMIT ARE CHOSEN, HLIMIT = ',F8.3/,
     2  'MIN DEPTH ALLOWED BEFORE DRYING, HMIN = ',F8.3/,
     3  'FACTOR DETERMINING FILLING OF NODE, VOLFCT = ',F10.4/,
     4  'INCREASE IN DEPTH TO BECOME WET, HREDRY = ',F7.4)
        CALL DRYND
        IF(NMBDRY .GT. MAXDRY)THEN
           WRITE(*,'(''NMBDRY TOO LARGE'', 2I7)')NMBDRY,MAXDRY
           STOP
        ENDIF
        IF(NMBDRY .GT. 0)WRITE(2,1012)NMBDRY,NDRY(1),NDRY(NMBDRY)
 1012 FORMAT('# DRYING NODES, NMBDRY = ',I5,' FIRST AND LAST ',2I5)
      ENDIF
      RETURN
      END
      SUBROUTINE INTIME
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
cc    taking the sum of all external mode etas that occur bet two internal
cc    modes for finding the average.
      ETAAV(I)=ETAAV(I)+ETA(I)
      ETAPRV(I)=ETA(I)
      HTPRV(I)=HT(I)
      SYSFQ(2*I-1)=0.0
      SYSFQ(2*I)=0.0
 10   CONTINUE
C*** IFRIC=0 => bottom stress from internal mode, set CF=1
C***      =2 => from external mode using cf(1)
C***      =1 => from external mode using variable cf
C    In order to preserve compatibility with 3D the CF must become
C    a nodal quantity.
C * * Update friction factor if manning n
      IF(IFRIC .EQ. 1 .AND. CFREF .GT. 5.)THEN
C * * There are 3 situations: 3D, CF is nodal (works for 2D too);
C     2D compatibility, CF is element quantity so must be
C     converted to equivalent nodal; and old CF, elemental and
C     remain elemental, this requires indicated changes in FRC23Q.
        AN=CFREF-5.
C3D        DO 20 I=1,NMNP
C3D           IF(HT(I) .LT. 0.)WRITE(*,*)'NEGATIVE DEPTH AT NODE',I,HT(I)
C3D          CF(I)=AN**2*9.81/(HT(I)**.333333)
C3D 20     CONTINUE
C * begin old 2d section, CF elemental, make chgs in FRC23Q
C        DO 125 I=1,NMEL
C          DE=0.
C          DO 128 J=1,3
C            J1=ICON(I,J)
C            DE=DE+HT(J1)
C 128      CONTINUE
C          IF(DE .LT. 0.)WRITE(*,*)'NEG DEPTH',I,DE
C          CF(I)=AN**2*9.81/(DE/3.)**.333333
C 125    CONTINUE
C * end old 2d section
C * Begin 2D compatibility, make element CF to nodal CF
        DO 140 I=1,NMNP
          CF(I)=0.
 140    CONTINUE
        DO 125 I=1,NMEL
          DE=0.
          DO 128 J=1,3
            DE=DE+HT(ICON(I,J))
 128      CONTINUE
C          IF(DE .LT. 0.)WRITE(*,*)'NEG DEPTH',I,DE,ICON(I,1),
C     1      ICON(I,2),ICON(I,3)
          IF(DE .LT. 0.)WRITE(*,*)'NEG DEPTH',I,DE,ICON(I,1),
     1      ICON(I,2),ICON(I,3),'WindSpd',WINDSP
C	Altered depth (DE/2.2) to simulate increased friction as
C       if shallower
          CFEL=AN*AN*GRAVT/(DE/2.2)**0.5
C         CFEL=AN*AN*GRAVT/(DE/3.)**0.333333
          DO 150 J=1,3
            CF(ICON(I,J))=CF(ICON(I,J))+CFEL*AREA(I)
 150      CONTINUE
 125    CONTINUE
        DO 130 I=1,NMNP
          CF(I)=CF(I)/AREAND(I)
 130    CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE LGTAUB(QBANG,NQN,X,Y,GLTOLO,NMQB)
      COMMON/ANGLE/S,C
      DIMENSION QBANG(NMQB),X(1),Y(1),NQN(NMQB)
      DO 10 I=1,NMQB
      I1=NQN(I)
      ANG=GLTOLO*QBANG(I)
      CALL TRIGO(ANG)
      CALL ROTV(X(I1),Y(I1))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE TURBCL
C * Compute turbulence closure using Mellor-Yamada level 2 scheme
C   The final product are the vertical viscosity and diffusion coefs.
      INCLUDE 'ALL.CMM'
      DATA A1MY,A2MY,B1MY,B2MY,C1MY/0.92,0.74,16.6,10.1,0.08/
      DATA GAMMA1,GAMMA2/0.2225,0.94096/
      DATA ALPHAL/0.1/
c * Assume global vel in Q
      DO 10 I=1,NMNP
C * Avoid qsum=0 initially, so set some small value
        QZSUM=1.E-9
        QSUM=1.E-9
        DO 30 LEV=1,NLEVM1
          QSUM=0.5*(TQ(I,LEV)+TQ(I,LEV+1))*DSIGMA(LEV)+QSUM
          QZSUM=0.5*(TQ(I,LEV)*(1.+SIGMA(LEV))+TQ(I,LEV+1)*
     1      (1.+SIGMA(LEV+1)))*DSIGMA(LEV)+QZSUM
 30     CONTINUE
        AL0MY=QZSUM/QSUM*ALPHAL*HT(I)
        DO 20 LEV=2,NLEVM1
          DUDZ=COF1(LEV)*Q(2*I-1,LEV)+COF2(LEV)*Q(2*I-1,LEV+1)
     1     -COF3(LEV)*Q(2*I-1,LEV-1)
          DUDZ=DUDZ/HT(I)
          DVDZ=COF1(LEV)*Q(2*I,LEV)+COF2(LEV)*Q(2*I,LEV+1)
     1     -COF3(LEV)*Q(2*I,LEV-1)
          DVDZ=DVDZ/HT(I)
C * The drho/dz term is computed in DNSTRM. Richardson gradient no
c Eliminate problems related to roundoff errors in drhodz and
c unstable situations e.g. caused by evap.
          IF(DRHODZ(I,LEV) .GT. -0.01)DRHODZ(I,LEV)=0.
          RICHG=-REFDEN*DRHODZ(I,LEV)/(DUDZ*DUDZ+DVDZ*DVDZ+1.E-10)
C *CHECK* For now assume unstable strat is like neutral strat
          RICHG=AMAX1(0.,RICHG)
C * Need the length scale (l Mellor and Yamada) to compute tq, -1<SIGMA<0
          ALMY=AL0MY*VKRMN*(1.+SIGMA(LEV))*HT(I)/
     1         (VKRMN*(1.+SIGMA(LEV))*HT(I)+AL0MY)
C * Arg is q^3
          ARG=ALMY*B1MY*(EDZZ(I,LEV)*(DUDZ*DUDZ+DVDZ*DVDZ)
     1      +GRAVT/REFDEN*DRHODZ(I,LEV)*VDIF(LEV,I))
          IF(ARG .LT. 0.)THEN
C            WRITE(*,*)'Total destruction of turb',i,itime
            TQ(I,LEV)=0.
          ELSE
            TQ(I,LEV)=ARG**0.33333
          ENDIF
C * Compute Richardson flux number w/out saving.
          RICHF=RICHG/EDZZ(I,LEV)*VDIF(LEV,I)
C * Compute SH and SM locally w/out saving.
          SH=3.*A2MY*(GAMMA1-(GAMMA1+GAMMA2)*RICHF)
     1       /(1.-RICHF)
C * Richf > gamma1/(gamma1+gamma2) means turbulence is destroyed
          SH=AMAX1(0.,SH)
          SM=A1MY/A2MY*(B1MY*(GAMMA1-C1MY)-(B1MY*(GAMMA1-C1MY)
     1      +6.*(A1MY+3.*A2MY))*RICHF)/(B1MY*GAMMA1-(B1MY*(GAMMA1
     2      +GAMMA2)-3.*A1MY)*RICHF)*SH
C * The fixed length AL0MY should be a function of wind and waves
C   Eddy viscosity should be at least kinematic viscosity
          EDZZ(I,LEV)=TQ(I,LEV)*ALMY*SM+EDYVIS
          VDIF(LEV,I)=TQ(I,LEV)*ALMY*SH+EDYVIS
 20     CONTINUE
C *  Determine VDIF and EDZZ at surface
        DUDZ=(Q(2*I-1,NLEV)-Q(2*I-1,NLEVM1))/DSIGMA(NLEVM1)/HT(I)
        DVDZ=(Q(2*I,NLEV)-Q(2*I,NLEVM1))/DSIGMA(NLEVM1)/HT(I)
        IF(DRHODZ(I,LEV) .GT. -0.01)DRHODZ(I,LEV)=0.
        RICHG=-REFDEN*DRHODZ(I,NLEV)/(DUDZ*DUDZ+DVDZ*DVDZ+1.E-20)
        RICHG=AMAX1(0.,RICHG)
        ALMY=AL0MY*VKRMN*HT(I)/(VKRMN*HT(I)+AL0MY)
        ARG=ALMY*B1MY*(EDZZ(I,NLEV)*(DUDZ*DUDZ+DVDZ*DVDZ)
     1      +GRAVT/REFDEN*DRHODZ(I,NLEV)*VDIF(NLEV,I))
        IF(ARG .LT. 0.)THEN
          TQ(I,NLEV)=0.
        ELSE
          TQ(I,NLEV)=ARG**0.33333
        ENDIF
        RICHF=RICHG/EDZZ(I,NLEV)*VDIF(NLEV,I)
        SH=3.*A2MY*(GAMMA1-(GAMMA1+GAMMA2)*RICHF)
     1    /(1.-RICHF)
        SH=AMAX1(0.,SH)
        SM=A1MY/A2MY*(B1MY*(GAMMA1-C1MY)-(B1MY*(GAMMA1-C1MY)
     1    +6.*(A1MY+3.*A2MY))*RICHF)/(B1MY*GAMMA1-(B1MY*(GAMMA1
     2    +GAMMA2)-3.*A1MY)*RICHF)*SH
        EDZZ(I,NLEV)=TQ(I,NLEV)*ALMY*SM+EDYVIS
        VDIF(NLEV,I)=TQ(I,NLEV)*ALMY*SH+EDYVIS
 10   CONTINUE
      RETURN
      END
      SUBROUTINE LOCGLO(QBANG,NQN,Q,GLTOLO,NMQB)
      COMMON/ANGLE/S,C
      DIMENSION QBANG(NMQB),Q(1),NQN(NMQB)
C * * GLTOLO = 1 From glob to loc.  = -1 From loc. to glob.
      DO 10 I=1,NMQB
      I1=NQN(I)
      ANG=GLTOLO*QBANG(I)
      CALL TRIGO(ANG)
      CALL ROTV(Q(2*I1-1),Q(2*I1))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE MKRHS
      INCLUDE 'ALL.CMM'
C * * DTAUX-YDZ CONTAIN D(TAU)/DSIGMA*H, Q IS VELOCITY FROM
C * * HORIZONTAL SOLUTION AT NEW TIME STEP (n+*), In local coords!
      DTAUXDZ(2)=DTAUXDZ(2)*TINC*0.5+Q(2*NK-1,2)
      DTAUYDZ(2)=DTAUYDZ(2)*TINC*0.5+Q(2*NK,2)
      DO 20 LEV=3,NLEV-1
        DTAUXDZ(LEV)=DTAUXDZ(LEV)*TINC*0.5+Q(2*NK-1,LEV)
        DTAUYDZ(LEV)=DTAUYDZ(LEV)*TINC*0.5+Q(2*NK,LEV)
 20   CONTINUE
      DTAUXDZ(NLEV)=DTAUXDZ(NLEV)*TINC*0.5+Q(2*NK-1,NLEV)+
     1               TAUWLX(NK)*TINC/HTAV(NK)/DSIGMA(NLEV-1)
      DTAUYDZ(NLEV)=DTAUYDZ(NLEV)*TINC*0.5+Q(2*NK,NLEV)+
     1               TAUWLY(NK)*TINC/HTAV(NK)/DSIGMA(NLEV-1)
      RETURN
      END
      SUBROUTINE NDELMC(I)
C WET AND DRY
      INCLUDE 'ALL.CMM'
C * * NODE ELEMENT CONNECTIVITY DETERMINES WHICH ELEMENTS SURROUND
C * * A NODE.  NDELM(N,1) = # CONNECTING ELEMENTS FOR NODE N
C * * NDELM(N,J) = CONNECTING ELEMENT NUMBERS, J=2,3,....
      ICOUNT=0
      DO 20 J=1,NMEL
        DO 30 K=1,3
          IF(ICON(J,K) .EQ. I)THEN
            ICOUNT=ICOUNT+1
            NDELM(NMBDRY,ICOUNT+1)=J
          ENDIF
 30     CONTINUE
 20   CONTINUE
      NDELM(NMBDRY,1)=ICOUNT
      IF(ICOUNT .GE. MAXNBR)THEN
        WRITE(*,*)'NDELMC: ICOUNT, MAXNBR',I,ICOUNT,MAXNBR
        STOP
      ENDIF
      RETURN
      END
      SUBROUTINE CHKDRY
C WET AND DRY
      INCLUDE 'ALL.CMM'
C * *  KEEP DRY UNTIL H>HMIN+HREDRY
      DO 20 I=1,NMNP
        VOLSRC(I)=0.
 20   CONTINUE
      DO 10 I=1,NMBDRY
        NODE=NDRY(I)
        HT(NODE)=ETA(NODE)-DEPTH(NODE)
        IF(HT(NODE) .LT. HMIN)THEN
C * * NEGATIVE DEPTH OR DEPTH TOO SMALL.  THIS CAN ONLY OCCUR
C * * THE FIRST TIME NODE GOES DRY BECAUSE THEN DIV(Q)=0
C          DH=HMIN-HT(NODE)
C * *  SET DEPTH TO SMALL POSITIVE VALUE (HMIN) TO AVOID DIVIDE BY ZERO
C * *  0.000001 IS ADDED FOR THE ABOVE LESS THAN TEST.
          HT(NODE)=HMIN+0.000001
          ETA(NODE)=HT(NODE)+DEPTH(NODE)
          IACTDR(NODE)=1
C * *  CORRECT FOR VOLUME ADDED
c          DVOL=DH*AREAND(NODE)*.333333/TINC
c          CALL SUBVOL(I,NODE,DVOL)
        ELSE
          IF(IACTDR(NODE) .EQ. 0)GOTO 10
C * *  KEEP DRY FOR AT LEAST TWO TIME STEPS
          IF(HT(NODE) .LT. HMIN+HREDRY .OR. IACTDR(NODE) .LT. 2)THEN
            IACTDR(NODE)=IACTDR(NODE)+1
            CALL REFLOD(I)
          ELSE
C * *  REVERT TO NORMAL NODE IF H>HMIN+REDRY
              IACTDR(NODE)=0
C            ENDIF
          ENDIF
        ENDIF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE OCDSCH
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMOSGM
      JEND=NMOBN(I)-1
      BDIS=0.
      DO 20 J=1,JEND
      N1=ICONO(I,J)
      N2=ICONO(I,J+1)
      XORDN1=XORD(N1)
      XORDN2=XORD(N2)
      YORDN1=YORD(N1)
      YORDN2=YORD(N2)
      QXN1=QT(2*N1-1)
      QXN2=QT(2*N2-1)
      QYN1=QT(2*N1)
      QYN2=QT(2*N2)
      IF(J .EQ. 1)QXN1=0.
      IF(J .EQ. 1)QYN1=0.
      IF(J .EQ. JEND)QXN2=0.
      IF(J .EQ. JEND)QYN2=0.
      DEX=XORDN1-XORDN2
      DEY=YORDN2-YORDN1
      DISCX=DEY*(QXN1+QXN2)/2.
      DISCY=DEX*(QYN1+QYN2)/2.
      DIS=DISCX+DISCY
C      DISCHS(J)=DIS
      BDIS=BDIS+DIS
 20   CONTINUE
C      CBDIS(I)=CBDIS(I)+BDIS*TINC
      CBDIS(I)=BDIS
C      IF(ITIME/NOUT*NOUT .NE. ITIME)GO TO 10
C      WRITE(2,1092)I,CBDIS(I),BDIS
C 1092 FORMAT(5X,'* * BOUNDARY',I4,5X,'CUMULATED DISCHARGE:',
C     1  E15.6,5X,'INSTANTANEOUS DISCHARGE:',E15.6)
C      WRITE(2,1094)(DISCHS(J), J=1,JEND)
C 1094 FORMAT(10E13.6)
 10   CONTINUE
      WRITE(31,1002)ITIME,(CBDIS(I), I=1,NMOSGM)
 1002 FORMAT(I8,12F8.1)
      RETURN
      END
      SUBROUTINE OCNSGM
      INCLUDE 'ALL.CMM'
      READ(1,1002)NMOSGM,(NMOBN(I), I=1,NMOSGM)
 1002 FORMAT(16I5)
      WRITE(IUOUT,1004) NMOSGM
 1004 FORMAT(10X,'OPEN BOUNDARY INFORMATION IS:'/15X,
     2  'NUMBER OF OCEAN SEGMENTS, NMOSGM = ',I5//)
      IF(NMOSGM .EQ. 0)RETURN
      DO 20 I=1,NMOSGM
      J1=NMOBN(I)
      READ(1,1006)(ICONO(I,J), J=1,J1)
 1006 FORMAT(20I4)
      WRITE(IUOUT,1008)I,NMOBN(I),(ICONO(I,J), J=1,J1)
 1008 FORMAT(5X,'SEGMENT ',I3,',  NUMBER OF NODES, NMOBN = ',
     1  I3,',  EXTERNAL NODE NUMBERS: ',15('-',I4))
      DO 30 J=1,J1
      IUNIV=ICONO(I,J)
      ICONO(I,J)=NINT(IUNIV)
 30   CONTINUE
      WRITE(IUOUT,1010) (ICONO(I,J), J=1,J1)
 1010 FORMAT(50X,'INTERNAL NODE NUMBERS: ',15('-',I4)//)
 20   CONTINUE
      RETURN
      END
      SUBROUTINE OUTP2D(IOTPT)
      INCLUDE 'ALL.CMM'
      IF(ITIME/NOUT*NOUT .EQ. ITIME)GO TO 10
      IOTPT=0
      RETURN
 10   IOTPT=1
      RETURN
      END
      SUBROUTINE OUTP3D
      INCLUDE 'ALL.CMM'
C * * NOTE IMULT ONLY CONTROLS 3-D.  SET VERY LARGE TO NOT CALCUL 3D.
      IF(ITIME .LT. ISTRT .OR . ITIME .LT. IMULT)RETURN
      IF(ITIME/NOUT*NOUT .EQ. ITIME)THEN
      WRITE(IUOUT,1034)TIME,TINC,ITIME
 1034 FORMAT(2X,'TIME = ',F15.1,' SEC',2X,'DELTA T, TINC = ',
     1  F8.2,1X,'SEC.',2X,'ITIME = ',I8)
      CALL VOLUME
      WRITE(IUOUT,1094)
 1094 FORMAT('NODE LEVEL#    U       V   OMEGA')
C *CHECK* limits of this do loop.
      DO 315 I=2,NMNP,7
C       USUM=0.
C       VSUM=0.
C       DO 328 LEV=1,NLEV-1
C        USUM=USUM+(U(I,LEV)+U(I,LEV+1))*.5*DSIGMA(LEV)
C        VSUM=VSUM+(V(I,LEV)+V(I,LEV+1))*.5*DSIGMA(LEV)
C 328  CONTINUE
C      WRITE(IUOUT,1032)I,NEXT(I),HT(I),ETAAV(I),UT(I),USUM,VT(I),VSUM,
C 1032 FORMAT (2I4,2F10.5,4F8.4/(I3,3E12.4))
       WRITE(IUOUT,1002)(I,
     1  LEV,U(I,LEV),V(I,LEV),OMEGA(I,LEV), LEV=2,NLEV)
 1002  FORMAT(I4,I6,2F8.4,E11.3)
 315  CONTINUE
      CALL SLINE(42)
      ENDIF
      RETURN
      END
      SUBROUTINE QMAT(NMNP,DIAGH,DIAGQ,MAXNOD,MAXMQ)
      DIMENSION DIAGH(MAXNOD),DIAGQ(MAXMQ)
      DO 10 I=1,NMNP
        DIAGQ(2*I-1)=DIAGH(I)
        DIAGQ(2*I)=DIAGH(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE QTALLY
      INCLUDE 'ALL.CMM'
C * * Q is u
      DO 10 I=1,NMNP
      QSUMX=0.
      QSUMY=0.
      DO 20 LEV=1,NLEV-1
      QSUMX=QSUMX+(Q(2*I-1,LEV+1)+Q(2*I-1,LEV))*0.5*DSIGMA(LEV)
      QSUMY=QSUMY+(Q(2*I,LEV+1)+Q(2*I,LEV))*0.5*DSIGMA(LEV)
 20   CONTINUE
C * * Rather than scaling the whole profile, a velocity is added.
C * * Because Q(  ,1) = 0. Correction is (1-.5*dsigma(1))
C * * Use QTPREV/HT instead of UT and VT because need them in local coord.
      DO 30 LEV=2,NLEV
        Q(2*I-1,LEV)=Q(2*I-1,LEV)+
     1  (QTPREV(2*I-1)/HT(I)-QSUMX)/(1.-.5*DSIGMA(1))
        Q(2*I,LEV)=Q(2*I,LEV)+
     1  (QTPREV(2*I)/HT(I)-QSUMY)/(1.-.5*DSIGMA(1))
 30   CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE RADIAT
      INCLUDE 'ALL.CMM'
C * * QB NODES
      DO 10 I=1,NMRAD
      NQ=NRAD(I)
      NODE=NQN(NQ)
      C=SQRT(GRAVT*HT(NODE))
C      IF(YORD(NODE) .GT. 5000.)THEN
C * * The following condition works well at the north boundary
C        etainn=.15*sin(phase*time)+.1*cos(phase*(time-6000.))
C         etainn=-0.03*cos(phase*(time+8000.))
C      ELSE
C         ETAINN=-.065*sin(phase*(time-1600.))
C     1  +.1*cos(phase*(time-6000.)))
C       ETAINN=.2*SIN(PHASE*TIME)
C      ENDIF
      ETAINN=0.
      QTB(NQ)=C*(ETA(NODE)-2.*etainn)
      DO 20 LEV=2,NLEV
        QB(NQ,LEV)=QTB(NQ)/HT(NODE)
 20   CONTINUE
C * * SIGN IS ACCOUNTED FOR IF QBANG IS OUTWARD NORMAL.  HAS TO BE!
 10   CONTINUE
      RETURN
      END
      SUBROUTINE READX(X,NMN,IUNIT)
      DIMENSION X(NMN)
      READ(IUNIT)(X(I), I=1,NMN)
      RETURN
      END
      SUBROUTINE ROTV(A,B)
      COMMON/ANGLE/ S,C
      X=A*C+B*S
      B=-A*S+B*C
      A=X
      RETURN
      END
      SUBROUTINE REFLOD(IDRY)
C WET AND DRY
      INCLUDE 'ALL.CMM'
C * * FIND SURROUNDING WET NODES AND TRANSFER VOLUME USING VOLSRC
C * * SHOULD ONLY BE CALLED FOR IACTDR>1.
      N=NDELM(IDRY,1)
      NODE=NDRY(IDRY)
      DO 10 I=1,N
        NEL=NDELM(IDRY,I+1)
        DO 20 J=1,3
          IJ=ICON(NEL,J)
C * * IF ADJACENT NODE HAS HIGHER ELEVATION AND IS NOT DRY, CREATE VOL SOURCE
          IF(ETA(IJ) .GT. ETA(NODE) .AND. IACTDR(IJ) .EQ. 0)THEN
            DETA=ETA(IJ)-ETA(NODE)
            DVOL=VOLFCT*DETA*AREA(NEL)/TINC
            VOLSRC(IJ)=VOLSRC(IJ)-DVOL
            SRC3DH(IJ)=SRC3DH(IJ)-DVOL
            VOLSRC(NODE)=VOLSRC(NODE)+DVOL
            SRC3DH(NODE)=SRC3DH(NODE)+DVOL
          ENDIF
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE SLINE(N)
      CHARACTER*3 STAR
      DATA STAR/'* *'/
      WRITE(2,1002)(STAR,I=1,N)
 1002 FORMAT(5X,42A3)
      RETURN
      END
      SUBROUTINE SLVOMG
      INCLUDE 'ALL.CMM'
      DO 50 I=1,NMNP
        OMEGA(I,LEV+1)=0.
 50   CONTINUE
      DO 10 I=1,NMEL
      VAR=0.
      DHAVDX=0.
      DUAVDX=0.
      DHAVDY=0.
      DVAVDY=0.
      DO 20 J=1,3
      KJ=ICON(I,J)
      DHAVDX=DHAVDX+B(I,J)*HTAV(KJ)
      DUAVDX=DUAVDX+B(I,J)*UAVG(KJ)
      DHAVDY=DHAVDY+A(I,J)*HTAV(KJ)
      DVAVDY=DVAVDY+A(I,J)*VAVG(KJ)
 20   CONTINUE
      DO 30 K=1,3
      KK=ICON(I,K)
      VAR=-(ETADT(KK)+((DUAVDX+DVAVDY)*HTAV(KK)+DHAVDX*UAVG(KK)+
     1     DHAVDY*VAVG(KK))*0.5/AREA(I))*AREA(I)/3.
      OMEGA(KK,LEV+1)=OMEGA(KK,LEV+1)+VAR
 30   CONTINUE
 10   CONTINUE
C ****  Thus RHS of eqn. is temporarily in omega. now solve for omega.
      DO 40 I=1,NMNP
C * * SRC3DH is total source flow over depth averaged over imult
C       external steps. It should be converted to source flow
C       for layer depth, so need to divide by HT and multiply
C       by dsigma*HT, or simply multiply by dsigma.
C       See \notes\source.tex
C      write(*,*)omega(i,lev+1),src3dh(i)
C      OMEGA(I,LEV+1)=OMEGA(I,LEV+1)+SRC3DH(I)*DSIGMA(LEV)
      OMEGA(I,LEV+1)=DSIGMA(LEV)*OMEGA(I,LEV+1)/DIAGH(I)+OMEGA(I,LEV)
 40   CONTINUE
      RETURN
      END
      SUBROUTINE SMTH
      INCLUDE 'ALL.CMM'
      DIMENSION X0(MXNDSL)
      DO 10 I=1,NMEL
      XCENT(I)=0.
      DO 20 J=1,3
      IJ=ICON(I,J)
      XCENT(I)=XCENT(I)+.3333333*ETA(IJ)
 20   CONTINUE
 10   CONTINUE
      DO 70 I=1,NMNP
      X0(I)=0.
 70   CONTINUE
      DO 40 I=1,NMEL
      DO 50 J=1,3
      IJ=ICON(I,J)
      X0(IJ)=X0(IJ)+XCENT(I)*AREA(I)
 50   CONTINUE
 40   CONTINUE
      DO 60 IJ=1,NMNP
      X0(IJ)=X0(IJ)/AREAND(IJ)
      ETA(IJ)=ASMTH*X0(IJ)+(1.-ASMTH)*ETA(IJ)
 60   CONTINUE
      RETURN
      END
      SUBROUTINE SOLVS
C * Horizontal diffusion sol.
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      DO 10 I=1,NMNP
      CONSTIT(I,LEV,ICONSTIT)=SYSF(I)/SYSM(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE SOLVX(X,SYSFX,XPREV,SYSMX,NMN)
      INCLUDE 'CINTEG.CMM'
      DIMENSION X(NMN),SYSFX(NMN),XPREV(NMN),SYSMX(NMN)
      DO 10 I=1,NMN
      X(I)=TINC*SYSFX(I)/SYSMX(I)+XPREV(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE STQ0B(Q,NQ0N,NMQ0BN)
      DIMENSION Q(1),NQ0N(NMQ0BN)
      DO 10 I=1,NMQ0BN
      J=NQ0N(I)
      Q(2*J-1)=0.
 10   CONTINUE
      RETURN
      END
      SUBROUTINE STQB(Q,QB,NQN,NVN)
      INCLUDE 'PARM2.CMM'
      INCLUDE 'CGRID.CMM'
      INCLUDE 'CINTEG.CMM'
      DIMENSION Q(1),QB(1),NQN(1),NVN(1)
      DO 10 I=1,NMQBN
      I1=2*NQN(I)-1
C     Q(I1)=QB(I)*(-COS(PHASE*TIME))
      Q(I1)=QB(I)
 10   CONTINUE
      IF(NMVBN .EQ. 0) GO TO 30
      DO 20 I=1,NMVBN
      I1=2*NVN(I)
      Q(I1)=0.
 20   CONTINUE
 30   RETURN
      END
      SUBROUTINE STQDRY
C WET AND DRY, EXTERNAL MODE
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMBDRY
        NODE=NDRY(I)
        IF(IACTDR(NODE) .EQ. 0)GOTO 10
        N=NDELM(I,1)
        DO 20 J=1,N
          NEL=NDELM(I,J+1)
          DO 30 K=1,3
            IJ=ICON(NEL,K)
            QT(2*IJ-1)=0.
            QT(2*IJ)=0.
C * * Also zero internal velocities so that internal mode uses
C       zeroed velocities for gradients.
            DO 40 IL=2,NLEV
              Q(2*IJ-1,IL)=0.
              Q(2*IJ,IL)=0.
 40         CONTINUE
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE STQ3DR
C WET AND DRY, INTERNAL MODE
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMBDRY
        NODE=NDRY(I)
        IF(IACTDR(NODE) .EQ. 0)GOTO 10
        N=NDELM(I,1)
        DO 20 J=1,N
          NEL=NDELM(I,J+1)
          DO 30 K=1,3
            IJ=ICON(NEL,K)
            DO 40 LEV=2,NLEV
              Q(2*IJ-1,LEV)=0.
              Q(2*IJ,LEV)=0.
 40         CONTINUE
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE SUBVOL(IDRY,NODE,DVOL)
C WET AND DRY
C     No longer used.  Since qt's are set to zero, the eta
C     initiating drying is not meaningful. So computing a volume
C     based on it is also meaningless.
C * * THIS ROUTINE PRESERVES VOLUME THE FIRST TIME STEP WHEN A
C * * NODE GOES DRY.  THE VOLUME NEEDED IS TO BRING H TO HMIN.
C * * FIND SURROUNDING NODES WITH POSITIVE DEPTH AND SAVE IN JTMPST
C * * THEN CREATE VOLSRC BY SUBTRACTING DVOL TO CONSERVE MASS
      INCLUDE 'ALL.CMM'
      DIMENSION JTMPST(20)
      N=NDELM(IDRY,1)
      NUMB=0
      DO 10 I=1,N
        NEL=NDELM(IDRY,I+1)
        DO 20 J=1,3
          IJ=ICON(NEL,J)
          IF(HT(IJ) .GT. HT(NODE))THEN
            NUMB=NUMB+1
            JTMPST(NUMB)=IJ
          ENDIF
 20     CONTINUE
 10   CONTINUE
      DO 30 I=1,NUMB
        VOLSRC(JTMPST(I))=-DVOL/NUMB
        SRC3DH(JTMPST(I))=VOLSRC(JTMPST(I))+SRC3DH(JTMPST(I))
 30   CONTINUE
      RETURN
      END
      SUBROUTINE TRIGO(A)
      COMMON/ANGLE/S,C
      S=SIN(A)
      C=COS(A)
      RETURN
      END
      SUBROUTINE VEL(Q,U,V)
      INCLUDE 'PARM2.CMM'
      INCLUDE 'CGRID.CMM'
      DIMENSION Q(1),U(1),V(1)
      DO 10 I=1,NMNP
      U(I)=Q(2*I-1)
      V(I)=Q(2*I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE VELT
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
CTODO   Change from qt to ubar/vbar for compatibility with 3-D
CTODO      UT(I)=QT(2*I-1)
CTODO      VT(I)=QT(2*I)
      UT(I)=QT(2*I-1)/HT(I)
      VT(I)=QT(2*I)/HT(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE UV3DH
      INCLUDE 'ALL.CMM'
C * * Levels start with 1 at bottom, NLEV is surface.
      DO 10 I=1,NMEL
        DETADX(I)=0.0
        DETADY(I)=0.0
        DHDX(I)=0.
        DHDY(I)=0.
        DO 20 J=1,3
          KJ=ICON(I,J)
          DETADX(I)=DETADX(I)+ETAAV(KJ)*B(I,J)
          DETADY(I)=DETADY(I)+ETAAV(KJ)*A(I,J)
          DHDX(I)=DHDX(I)+HTAV(KJ)*B(I,J)
          DHDY(I)=DHDY(I)+HTAV(KJ)*A(I,J)
 20     CONTINUE
        DETADX(I)=-GD6*DETADX(I)
        DETADY(I)=-GD6*DETADY(I)
        DHDX(I)=DHDX(I)*GD6DEN
        DHDY(I)=DHDY(I)*GD6DEN
 10   CONTINUE
      DO 30 LEV=2,NLEV
      CALL AMATZR(SYSFQ,NMNP2)
      CALL FRC3DQ
C     CALL BOUNDE(Q(1,LEV),EDXX,EDYY,EDXY,XORD,YORD,SYSFQ,
C    1  NMBEL,NBN,NMBEL3)
      CALL LOCGLO(QBANG(1,LEV),NQN,SYSFQ,1.,NMQBN)
      CALL LOCGLO(Q0BANG(1,LEV),NQ0N,SYSFQ,1.,NMQ0BN)
C * * Save forces for getting stress in vertical mode.
      IF(LEV .EQ. 2)THEN
C * * XFRC/yfrc are in local cord. so get sysfq already in local.
      DO 40 J=1,NMNP
        XFRC(J)=SYSFQ(2*J-1)/DIAGQ(2*J-1)
        YFRC(J)=SYSFQ(2*J)/DIAGQ(2*J)
 40   CONTINUE
      ENDIF
C * * The variables Q and QPREV really contain u and uprev
      CALL SOLVX(Q(1,LEV),SYSFQ,QPREV(1,LEV),DIAGQ,NMNP2)
C      IF(NMQBN. GT. 0)CALL STQB(Q(1,LEV),QB(1,LEV),NQN,NVN)
C      IF(NMQ0BN .GT. 0)CALL STQ0B(Q(1,LEV),NQ0N,NMQ0BN)
C      IF(NMBDRY .GT. 0)CALL STQ3DR
 30   CONTINUE
      RETURN
      END
      SUBROUTINE UVAVG
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
      UAVG(I)=0.5*(U(I,LEV)+U(I,LEV+1))
      VAVG(I)=0.5*(V(I,LEV)+V(I,LEV+1))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE VEDSTR
C * 073001 Make this the controlling routine for vertical diffusion
C   of momentum. Also change to one external loop over the nodes.
C * * Computes (dtau/dsigma)/(rhoH) at sigma levels at previous
C     internal time step
      INCLUDE 'ALL.CMM'
      NLEVM1=NLEV-1
      DTAUXDZ(1)=0.
      DTAUYDZ(1)=0.
      DO 10 NK=1,NMNP
C * * XFRC,YFRC CONTAINS 3D FORCE at sigma level 2, eg. -G(DETA/DX)
C       at n+* time level (just calculated in UV3DH.
C     TAUX=taux/rho is at the middle between sigma levels
C     If ifric=0, TAUBL is the bottom stress from the previous
C       internal time step, local coord.
C * * Ustar is taken as the (taub/rho+tauw/rho)^.5
        IF(IACTDR(NK) .GT. 0)GOTO 10
C * Change to Mellor-Yamada level 2 scheme. Make absust=1.
C        ABSUST(NK)=SQRT(TAUBLX(NK)*TAUBLX(NK)+TAUBLY(NK)*TAUBLY(NK))
C        ABSUST(NK)=ABSUST(NK)+SQRT(TAUWLX(NK)*TAUWLX(NK)+
C     1     TAUWLY(NK)*TAUWLY(NK))
C        ABSUST(NK)=SQRT(ABSUST(NK))
        ABSUST(NK)=1.
C * * Here taux and tauy are stress/rho. Use HTAV since HTPRV
C      at previous internal time step does not exist
          TAUX(1)=TAUBLX(NK)-HTAV(NK)*XFRC(NK)*DSIGMA(1)*.5
          TAUY(1)=TAUBLY(NK)-HTAV(NK)*YFRC(NK)*DSIGMA(1)*.5
          DO 30 LEV=2,NLEVM1
C * * QPREV is previous internal time step. EDZZ is mixing length, ie z.
C     in old closure scheme. In MY level 2, it is real EDZZ, 02/04/02
            TAUX(LEV)=(QPREV(2*NK-1,LEV+1)-QPREV(2*NK-1,LEV))/
C For MY level 2, 02/04/02, set absust=1. Done in vedstr
     1                DSIGMA(LEV)*.5*ABSUST(NK)*
     2                 (EDZZ(NK,LEV+1)+EDZZ(NK,LEV))
            TAUY(LEV)=(QPREV(2*NK,LEV+1)-QPREV(2*NK,LEV))/
     1                DSIGMA(LEV)*.5*ABSUST(NK)*
     2                 (EDZZ(NK,LEV+1)+EDZZ(NK,LEV))
C * * Should use HTPREV, but does not exist. So HTAV will have to do.
            TAUX(LEV)=TAUX(LEV)/HTAV(NK)
            TAUY(LEV)=TAUY(LEV)/HTAV(NK)
 30       CONTINUE
C  DTAU/DSIGMA is not reqd for lev=1, because u,v are zero.
C       DTAUXDZ(1)=2.*(TAUX(NK,1)-TAUBLX(NK))/(DSIGMA(1)*HTAV(NK))
C       DTAUYDZ(1)=2.*(TAUY(NK,1)-TAUBLY(NK))/(DSIGMA(1)*HTAV(NK))
C *  dtauxdz=(d(taux/rho)/dsigma)/H at sigma levels.  dtauydz similar for y
         DO 40 LEV=2,NLEVM1
           DTAUXDZ(LEV)=2.*(TAUX(LEV)-TAUX(LEV-1))/(DSIGMA(LEV)+
     1      DSIGMA(LEV-1))
           DTAUYDZ(LEV)=2.*(TAUY(LEV)-TAUY(LEV-1))/(DSIGMA(LEV)+
     1      DSIGMA(LEV-1))
           DTAUXDZ(LEV)=DTAUXDZ(LEV)/HTAV(NK)
           DTAUYDZ(LEV)=DTAUYDZ(LEV)/HTAV(NK)
 40      CONTINUE
C * Half sigma below surface, use one-sided difference
         DTAUXDZ(NLEV)=2.*(TAUWLX(NK)-TAUX(NLEVM1))
     1    /(DSIGMA(NLEVM1)*HTAV(NK))
         DTAUYDZ(NLEV)=2.*(TAUWLY(NK)-TAUY(NLEVM1))
     1    /(DSIGMA(NLEVM1)*HTAV(NK))
         CALL MKRHS
         CALL TRIDS
C       ENDIF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE VFEDBK(Q,QT,HT,DSIGMA,GXX,GXY,GYY,NMNP,NLEV)
      INCLUDE 'PARM2.CMM'
C *041701* These terms should be calc in global coords. Fixed.
      REAL*4 UPRIME(MXNDSL),VPRIME(MXNDSL),UVPRIM(MXNDSL),
     1  GXX(MXNDSL),GXY(MXNDSL),GYY(MXNDSL),HT(MXNDSL),
     2  Q(MXMQ,MXLVLS),QT(MXMQ),DSIGMA(MXLVLS)
      DO 10 I=1,NMNP
        SUMXX=0.0
        SUMXY=0.0
        SUMYY=0.0
        DO 20 LEV=1,NLEV
          UDASH=Q(2*I-1,LEV)-QT(2*I-1)/HT(I)
          VDASH=Q(2*I,LEV)-QT(2*I)/HT(I)
          UPRIME(LEV)=UDASH*UDASH
          VPRIME(LEV)=VDASH*VDASH
          UVPRIM(LEV)=UDASH*VDASH
 20     CONTINUE
        DO 30 LEV=1,NLEV-1
          SUMXX=SUMXX+(UPRIME(LEV)+UPRIME(LEV+1))*.5*DSIGMA(LEV)
          SUMYY=SUMYY+(VPRIME(LEV)+VPRIME(LEV+1))*.5*DSIGMA(LEV)
          SUMXY=SUMXY+(UVPRIM(LEV)+UVPRIM(LEV+1))*.5*DSIGMA(LEV)
 30     CONTINUE
      GXX(I)=SUMXX*HT(I)
      GXY(I)=SUMXY*HT(I)
      GYY(I)=SUMYY*HT(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE VOLUME
      INCLUDE 'ALL.CMM'
      VOL=0.
      DO 10 I=1,NMEL
      SUM=0.
      DO 20 J=1,3
      IJ=ICON(I,J)
      SUM=SUM+ETA(IJ)
 20   CONTINUE
      VOL=AREA(I)*SUM/3.+VOL
 10   CONTINUE
      WRITE(IUOUT,1035)VOL
 1035 FORMAT('NET EXCESS VOLUME',' = ',E14.6)
      CALL SLINE(25)
      RETURN
      END
      SUBROUTINE WNDLOC
      INCLUDE 'ALL.CMM'
C * *  TAUWLX/Y IS IN LOCAL COORDS FOR USE IN INTERNAL MODE AND IS
C * *  ROTATED JUST AFTER THIS ROUTINE IN THE MAIN PROG.
      DO 10 I=1,NMNP
      TAUWLX(I)=TAUWGX(I)
      TAUWLY(I)=TAUWGY(I)
  10  CONTINUE
      RETURN
      END
      SUBROUTINE WRITE1
      INCLUDE 'ALL.CMM'
      WRITE(IUOUT,1004)NMEL,NMNP,NLEV,IDEPTH,IEDVIS,
     1  IWIND,IRDHOT,ICNVEC,NCONSTIT,ICSALT
 1004 FORMAT(5X,'THIS PROBLEM HAS THE FOLLOWING CHARACTERISTICS:'//
     1  10X,'NUMBER OF ELEMENTS, NMEL = ',I6/,
     2  10X,'NUMBER OF NODES, NMNP = ',I6/,
     3  10X,'THE MODEL APPLIED HAS',I3,' LAYERS'/,
     4  10X,'IT IS ASSUMED THAT SPATIALLY,'/,
     6  40X,'MEAN LOW WATER DEPTH: IDEPTH=',I2,'; 1=CONST, 2=VAR'/,
     7  40X,'EDDY VISCOSITY: IEDVIS =  ',I2,'; 1=CONST, 2=VAR'/,
     8  40X,'WIND STRESS: IWIND = ',I2,'; 1=CONST,2=VAR'/,
     9  10X,'HOT START, IRDHOT = ',I2,'; 0=NO; 1=YES'/,
     A  10X,'CONVECTIVE ACCELERATIONS: ICNVEC =  ',I3,'; 0=NO, 1=EXT',
     B  ', 2=INT'/,10X,'NUMBER OF CONSTITUENTS, NCONSTIT = ',I5/,
     1  10X,'SALT IS ICSALT = ',I3)
      CALL SLINE(50)
      WRITE(IUOUT,1080)IFRIC
 1080 FORMAT('FRICTION SWITCH IFRIC =',I3/,
     1  'IFRIC=0 => bottom stress using law of the wall'/,
     2  'IFRIC=1 => bottom stress from external mode,',
     3  ' variable for each node'/,
     4  'IFRIC=2 => like 1, but constant for all nodes')
      WRITE(IUOUT,1020)KSTART,KDINCR,KDIHOT
 1020 FORMAT(15X,'START ON DISC, KSTART = ',I12/,15X,
     1  'INCREMENT BETWEEN OUTPUT TO DISC, KDINCR = ',I12/,15X,
     2  'TIME STEP STORED ON DISC FOR HOT START, KDIHOT = ',I12)
      IF (IPRINT .NE.0) WRITE(IUOUT,1010)
 1010 FORMAT(T5,'NODE',T14,'X-',T25,'Y-',T32,'DEPTH',T39,'NODE',T44,
     1  'TIDE',T51,'TIME',T57,'LOC X',T65,'FLUX',T73,'LOC Y',T80,
     3  'SOURCE',T88,'LEVEL'/,T4,'NUMBER',T13,'COORD.',T23,
     4  'COORD.',T32,'(MLW)',T39,'CODE',T44,'AMPL.',T51,'LAG',
     5  T57,'FLUX',T65,'ANG',T73,'FLUX',T80,'FLUX'/,T14,'(M)',
     6  T25,'(M)',T32,'(M)',T44,'(M)',T51,'(SEC)',T57,'(M2/S)',T65,
     7  '(DEGR.)',T73,'(M2/S)',T80,'(M/S)'/)
      RETURN
      END
      SUBROUTINE DOADVD
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
C *   For 1-layer solution, nlevad=1, else nlevad=nlev
      IADTM=IADTM+1
C * * This is previous time step's IT (number of iter. in diffusion)
      WRITE(*,1002)IADTM,NLEVAD,IT,TINCAD,TIME
      DO 20 LEV=1,NLEVAD
        DO 10 I=1,NMNP
          NDSTSL(I)=0
 10     CONTINUE
C * Advection is done for ISWAD < 3
        IF(ISWAD .LT. 3)THEN
C       there is CONSTIT at bottom, but there is no advection. Need to set
C       CONSTIT at bottom for 3D solution. Won't hurt 2D solution.
          IF(LEV .EQ. 1)THEN
            DO 40 I=1,NMNP
              DO 42 ICONSTIT=1,NCONSTIT
                CONSTIT(I,LEV,ICONSTIT)=CONSTITP(I,LEV,ICONSTIT)
 42           CONTINUE
 40         CONTINUE
          ENDIF
          IF(NLEVAD .EQ. 1 .OR. LEV .GT. 1)THEN
            CALL TRCBK
            DO 44 ICONSTIT=1,NCONSTIT
              CALL GTCONC
 44         CONTINUE
          ENDIF
C *CHECK*
C * With source flow distributed over depth, 2D and 3D is
C     treated the same, see \notes\source.tex
          DO 46 ICONSTIT=1,NCONSTIT
C *CHECK* Source only set up for salinity so far
            CALL SOURCE
            CALL GNSSAL
 46       CONTINUE
C * If source flow is forced to surface in 3D then use SRC3D instead
C     of SOURCE
C          CALL SRC3D
C * Evaporation is done. Groundwater could be done similarly
C      CALL PRCS1(NMNPL,DT)
C      DO 510 I=1,NMBN
C        NK=NB(I)
C        IF(NDSTSQ(NK) .LT. 0)CALL STCBC
C 510  CONTINUE
C *  Advance utmdt=u(t-dt)
        DO 60 I=1,NMNP
          IF(NLEVAD .GT. 1)THEN
            UTMDT(I,LEV)=U(I,LEV)
            VTMDT(I,LEV)=V(I,LEV)
            WTMDT(I,LEV)=W(I,LEV)
C * * Here W is omega/ht
          ELSE
            UTMDT(I,LEV)=UT(I)
            VTMDT(I,LEV)=VT(I)
          ENDIF
 60     CONTINUE
        ENDIF
C     Store newly calculated conc from advection into previous conc.
C     and apply first order decay. Always done so decay is applied
C     even if advection is not solved for.
C *CHECK* DECAY is hardwired here
        DECAY=0.
C        DECAY=0.6/3600.
        DO 30 I=1,NMNP
          DO 32 ICONSTIT=1,NCONSTIT
            CONSTITP(I,LEV,ICONSTIT)=CONSTIT(I,LEV,ICONSTIT)*(1.
     1         -DECAY*TINCAD)
            CONSTIT(I,LEV,ICONSTIT)=CONSTITP(I,LEV,ICONSTIT)
 32       CONTINUE
 30     CONTINUE
C * Diffusion for ISWAD > 1
        IF(ISWAD .GT. 1)THEN
          DO 65 ICONSTIT=1,NCONSTIT
            CALL DIFFUS
 65       CONTINUE
C * New solution is in CONSTIT
        ENDIF
 20   CONTINUE
      IF(NLEVAD .GT. 1 .AND. ISWAD .GT. 1)THEN
        DO 75 ICONSTIT=1,NCONSTIT
          CALL VDIFUS
 75     CONTINUE
      ENDIF
C TODO You could apply half the decay after advection and half here
C      after diffusion.
      DO 90 LEV=1,NLEVAD
        DO 70 I=1,NMNP
          DO 77 ICONSTIT=1,NCONSTIT
            CONSTITP(I,LEV,ICONSTIT)=CONSTIT(I,LEV,ICONSTIT)
 77       CONTINUE
 70     CONTINUE
C Write CONSTIT to disk
       IF(IADTM/NDISK*NDISK .EQ. IADTM)THEN
          DO ICONSTIT=1,NCONSTIT
            NRCDSK=NRCDSK+1
        WRITE(16,REC=NRCDSK)(CONSTITP(I,LEV,ICONSTIT),I=1,NMNP)
          ENDDO
        ENDIF
        IF(ICSALT .GT. 0)THEN
          DO I=1,NMNP
C          DENSTY(I,LEV)=999.84+.8245*CONSTITP(I,LEV,ICSALT)-REFDEN
            DENSTY(I,LEV)=CONSTITP(I,LEV,ICSALT)*0.7
          ENDDO
        ELSE
          DENSTY=1000.
        ENDIF
        IF(IADTM/NADOUT*NADOUT .EQ. IADTM .AND. ITIME .GE. ISTRT)THEN
          WRITE(IUOUT,1002)IADTM,LEV,IT,TINCAD,TIME
 1002     FORMAT('ADV-DIF, IADTM=',I8,' LEVEL=',I2,
     1     ', ITS.=',I3,' TINCAD=',F6.1,' TIME=',F15.1)
           WRITE(IUOUT,1004)(CONSTITP(I,LEV,1), I=1,NMNP,1)
 1004     FORMAT(F8.2,9F7.2)
c          write(17,1088)(CONSTITp(i,lev,ICSALT), i=2,nmnp,3)
c          write(17,1088)(q(2*i,lev), i=2,nmnp,3)
 1088     format(50f6.2)
c 1088     format(50f7.3)
        ENDIF
 90   CONTINUE
C**** now calc. the Density terms for next internal and external mode.
      CALL DNSTRM
      ICNTAD=0
      RETURN
      END
      SUBROUTINE FNDNBN
C * * Find nearest boundary node.
      INCLUDE 'ALL.CMM'
      DSTMIN=1.E+30
      NDNEAR=0
      DO I=1,NMEL
        DO J=1,3
          IJ=ICON(I,J)
          IF(NBC(IJ) .EQ. 0)CYCLE
          XDIF=XN-XORD(IJ)
          IF(ABS(XDIF) .GT. SRCHRDX(I))EXIT
          YDIF=YN-YORD(IJ)
          IF(ABS(YDIF) .GT. SRCHRDY(I))EXIT
          DIST=XDIF*XDIF+YDIF*YDIF
          IF(DIST .LT. DSTMIN)THEN
            DSTMIN=DIST
            NDNEAR=IJ
          ENDIF
        ENDDO
      ENDDO
      IF(NDNEAR .EQ. 0)THEN
C * If xn, yn are far outside of any element, it may not be found
C   Then try to increase the srchrd in GEOM.
        WRITE(*,*)'NDNEAR NOT FOUND IN FNDNBN',NK,XN,YN,LEV
        STOP
      ELSE
        NDSTSL(NK)=-NDNEAR
      ENDIF
      RETURN
      END
      SUBROUTINE GTCONC
C * Sets boundary concentrations for advection step
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
C * 050901  3-D ready. Called for each ICONSTIT,but check SALBC and LINTPL,
      CALL CONSTITBC
      DO 10 NK=1,NMNP
        IF(NDSTSL(NK) .LT. 0)THEN
          CONSTIT(NK,LEV,ICONSTIT)=CONSTITP(-NDSTSL(NK),LEV,ICONSTIT)
C * *  If NDSTSL is neg. the node traced out of domain and
C        NDSTSL contains the nearest boundary node.
C        If the nearest boundary node is the node itself, it
C        needs an external boundary condition.
          IF(NK .EQ. -NDSTSL(NK))CONSTIT(NK,LEV,ICONSTIT)=CONSTITBND(NK)
        ELSE
          IEL=NDSTSL(NK)
          CALL LINTPL
        ENDIF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE LINTPL
C * 050901 Only interpolates in fixed sigma
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      IJ1=ICON(IEL,1)
      IJ2=ICON(IEL,2)
      IJ3=ICON(IEL,3)
      A1=((XORD(IJ3)-XORD(IJ2))*(YOL(NK)-YORD(IJ2)))
     1   -((XOL(NK)-XORD(IJ2))*(YORD(IJ3)-YORD(IJ2)))
      A2=((XORD(IJ1)-XORD(IJ3))*(YOL(NK)-YORD(IJ3)))
     1   -((XOL(NK)-XORD(IJ3))*(YORD(IJ1)-YORD(IJ3)))
      A3=((XORD(IJ2)-XORD(IJ1))*(YOL(NK)-YORD(IJ1)))
     1   -((XOL(NK)-XORD(IJ1))*(YORD(IJ2)-YORD(IJ1)))
      EA2=AREA(IEL)*2.
      XI1=A1/EA2
      XI2=A2/EA2
      XI3=A3/EA2
      IF(NLEVAD .EQ. 1)THEN
        CONSTIT(NK,LEV,ICONSTIT)=XI1*CONSTITP(IJ1,LEV,ICONSTIT)+
     1          XI2*CONSTITP(IJ2,LEV,ICONSTIT)+
     2          XI3*CONSTITP(IJ3,LEV,ICONSTIT)
C * * Evaporation. Should use time centered values on RHS.
C     At the moment EVAP is the last value from the external
C     solution.  It would be better to use an average over the
C     last TINCAD sec.
C     For this to work properly, EVAP should be nodal
C     evaporation in m/sec. It is usually negative
C *   The following statement works for 1-layer situation
        CONSTIT(NK,NLEVAD,ICONSTIT)=CONSTIT(NK,NLEVAD,ICONSTIT)*
     1    (HT(NK)-EVAP*TINCAD)/HT(NK)
      ELSE
C * 3D. Find levzn
C * Check to see if above surf or below bott.
C   LEVZN is the sigma level immediately below zn.
        DO 20 ILEV=1,NLEVAD
          IF(ZOL(NK) .GT. SIGMA(ILEV))GOTO 20
          LEVZN=ILEV-1
          GOTO 30
 20     CONTINUE
 30     SIGINT=(ZOL(NK)-SIGMA(LEVZN))/DSIGMA(LEVZN)
        CONSTITL=XI1*CONSTITP(IJ1,LEVZN,ICONSTIT)+
     1          XI2*CONSTITP(IJ2,LEVZN,ICONSTIT)+
     2          XI3*CONSTITP(IJ3,LEVZN,ICONSTIT)
        CONSTITU=XI1*CONSTITP(IJ1,LEVZN+1,ICONSTIT)+
     1          XI2*CONSTITP(IJ2,LEVZN+1,ICONSTIT)+
     2          XI3*CONSTITP(IJ3,LEVZN+1,ICONSTIT)
        CONSTIT(NK,LEV,ICONSTIT)=CONSTITL+(CONSTITU-CONSTITL)*SIGINT
C *  Evaporation, only top salinity is changed, see \notes\source.tex
        IF(LEV .EQ. NLEVAD)THEN
          SUM1=0.
          SUM2=0.
          DO 40 ILEV=1,NLEVAD-2
            SUM1=SUM1+DSIGMA(ILEV)*CONSTIT(NK,ILEV,ICONSTIT)
            SUM2=SUM2+DSIGMA(ILEV)*CONSTIT(NK,ILEV+1,ICONSTIT)
 40       CONTINUE
          SUM1=SUM1+DSIGMA(NLEVAD-1)*CONSTIT(NK,NLEVAD-1,ICONSTIT)
          CONSTIT(NK,NLEVAD,ICONSTIT)=CONSTIT(NK,NLEVAD,ICONSTIT)*
     1      (HT(NK)-EVAP*TINCAD)/HT(NK)
     2      -EVAP*TINCAD*(SUM1+SUM2)/HT(NK)/DSIGMA(NLEVAD-1)
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE MAKEM(I,HTOTL)
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      DO 10 LEV=1,NLEVM1
        DEL6=DSIGMA(LEV)*HTOTL/6.
        ESUM=.25*TINCIN*(VDIF(LEV,I)+VDIF(LEV+1,I))/DSIGMA(LEV)/HTOTL
        SYSMVT(LEV,2)=SYSMVT(LEV,2)+2.*DEL6+ESUM
        SYSMVT(LEV,3)=SYSMVT(LEV,3)+DEL6-ESUM
        SYSMVT(LEV+1,1)=SYSMVT(LEV+1,1)+DEL6-ESUM
        SYSMVT(LEV+1,2)=SYSMVT(LEV+1,2)+2.*DEL6+ESUM
 10   CONTINUE
      RETURN
      END
      SUBROUTINE MAKERH(I,HTOTL)
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      DO 10 LEV=1,NLEVM1
        ESUM=-.25*TINCIN*(VDIF(LEV,I)+VDIF(LEV+1,I))/DSIGMA(LEV)/HTOTL
        TVAR=(CONSTIT(I,LEV,ICONSTIT)-CONSTIT(I,LEV+1,ICONSTIT))*ESUM
        RHS(LEV)=RHS(LEV)+TVAR
        RHS(LEV+1)=RHS(LEV+1)-TVAR
        DEL6=DSIGMA(LEV)*HTOTL/6.
        RHS(LEV)=RHS(LEV)+DEL6*(2.*CONSTIT(I,LEV,ICONSTIT)
     1    +CONSTIT(I,LEV+1,ICONSTIT))
        RHS(LEV+1)=RHS(LEV+1)+DEL6*(CONSTIT(I,LEV,ICONSTIT)
     1    +2.*CONSTIT(I,LEV+1,ICONSTIT))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE NORM
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      ANORM=0.
      AVG=0.
      DO 10 I=1,NMNP
        ANORM=ANORM+(CONSTITIT(I)-CONSTIT(I,LEV,ICONSTIT))*
     1    (CONSTITIT(I)-CONSTIT(I,LEV,ICONSTIT))
        AVG=AVG+CONSTIT(I,LEV,ICONSTIT)*CONSTIT(I,LEV,ICONSTIT)
 10   CONTINUE
      ANORM=SQRT(ANORM/AVG)
      RETURN
      END
      SUBROUTINE ONSIDE
C * * Force on adjacent side
      INCLUDE 'ALL.CMM'
      SCALR=(UT(NK)*BPARMS(NK,1)+VT(NK)*BPARMS(NK,2))*TINCAD
      IF(SCALR .GE. 0.)THEN
C * * Clockwise from node
        FRACT=SQRT(SCALR/BPARMS(NK,3))
        IF(FRACT .GT. 1.)FRACT=1.
        XOL(NK)=XORD(NK)+FRACT*BPARMS(NK,1)
        YOL(NK)=YORD(NK)+FRACT*BPARMS(NK,2)
        NDSTSL(NK)=IELMSB(NK,1)
      ELSE
        SCALR=(UT(NK)*BPARMS(NK,4)+VT(NK)*BPARMS(NK,5))*TINCAD
        IF(SCALR .GE. 0)THEN
C * * CClockwise from node
          FRACT=SQRT(SCALR/BPARMS(NK,6))
          IF(FRACT .GT. 1.)FRACT=1.
          XOL(NK)=XORD(NK)+FRACT*BPARMS(NK,4)
          YOL(NK)=YORD(NK)+FRACT*BPARMS(NK,5)
          NDSTSL(NK)=IELMSB(NK,2)
        ELSE
          WRITE(*,*)'BOTH SCALARS ARE NEG.',NK
          XOL(NK)=XORD(NK)
          YOL(NK)=YORD(NK)
          NDSTSL(NK)=IELMSB(NK,1)
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE PUSH
      INCLUDE 'ALL.CMM'
C * * LOOKUP(NK,1) HAS THE NUMBER OF ACTIVE ELEMENTS
      N=LOOKUP(NK,1)+1
      DO 10 I=N,3,-1
        LOOKUP(NK,I)=LOOKUP(NK,I-1)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE PPUSH(LK)
      INCLUDE 'ALL.CMM'
C * * PARTIAL PUSH TO MOVE ELEMENT TO FRONT KEEPING THE SEQUENCE
      SAVE NRCDSK      
      DIMENSION LTMP(MXLKUP)
      DO 20 I=1,LK-1
        LTMP(I)=LOOKUP(NK,I+1)
 20   CONTINUE
      N=LOOKUP(NK,1)
      DO 10 I=LK+1,N+1
        LOOKUP(NK,I-LK+1)=LOOKUP(NK,I)
 10   CONTINUE
      NMOVED=N-LK+1
      DO 30 I=1,LK-1
        LOOKUP(NK,NMOVED+I+1)=LTMP(I)
 30   CONTINUE
      RETURN
      END
C===============================================================C
C     SUB PROGRAM FOR SOLVING DIFFUSION IN FE GRID    C
C===============================================================C
      SUBROUTINE DIFFUS
C * Compute "horizontal" diffusion at one level, LEV, and for ICONSTIT
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
C       ASSEMBLING DIAGONAL GLOBAL MATRIX
      DO 300 I=1,NMNP
        SYSM(I)=0.
        SYSF(I)=0.
 300  CONTINUE
C * * Predictor step
      DO 310 M=1,NMEL
        AREA3=AREA(M)/3.
        DO 160 I=1,3
          KI=ICON(M,I)
          SYSM(KI)=SYSM(KI)+AREA3
          SYSF(KI)=SYSF(KI)+AREA3*CONSTITP(KI,LEV,ICONSTIT)
C * *  DIFXX contains diff dt(adv/dif)DXX/(8area)
          DO 170 J=1,3
            KJ=ICON(M,J)
            SYSM(KI)=SYSM(KI)+(DIFXX(M)*B(M,I)*B(M,J)
     1               +DIFYY(M)*A(M,I)*A(M,J))
            SYSF(KI)=SYSF(KI)-(DIFXX(M)*B(M,I)*B(M,J)
     1               +DIFYY(M)*A(M,I)*A(M,J))*CONSTITP(KJ,LEV,ICONSTIT)
 170      CONTINUE
 160    CONTINUE
 310  CONTINUE
C        WRITE(2,1004)(SYSM(I), I=1,NMNP)
 1004   FORMAT(10E11.3)
C        WRITE(2,1004)(SYSF(I), I=1,NMNP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SOURCE ADDITION TAKES PLACE HERE. ***THIS IS TO BE CUSTOMISED***  C
C     A SAMPLE COSINE LOADING IS SHOWN BELOW                            C
C     GIVE YOUR ELEMNT # WHERE SOURCE IS TO BE PUT = ISRCEL             C
C     ALSO A SOURCE RATE IN MASS/AREA/TIME = SRATE                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TIMER=TIME/TINCAD
C****FOLL ALL TEMPRY
C     SRATE=0.5*(1.-COS(PI*TIMER/18.))
C     ISRCE1=10
C     ISRCE2=11
C     IF(LEV.EQ.6) THEN ! LAYER 6
C     DO 220 J=4,6
C     SYSF(ICONQ(ISRCE1,J))=SYSF(ICONQ(ISRCE1,J)) +
C    1                      1./3.*TINCAD*SRATE*AREA(ISRCE1)
C220  SYSF(ICONQ(ISRCE2,J))=SYSF(ICONQ(ISRCE2,J)) +
C    1                      1./3.*TINCAD*SRATE*AREA(ISRCE2)
C      ENDIF
      CALL SOLVS
      CALL STPRCC
C * * Corrector steps
      IT=0
 450  IT=IT+1
      DO 390 I=1,NMNP
        SYSM(I)=0.
        SYSF(I)=0.
        CONSTITIT(I)=CONSTIT(I,LEV,ICONSTIT)
 390  CONTINUE
      DO 410 M=1,NMEL
        AREA12=AREA(M)/12.
        DO 420 I=1,3
          KI=ICON(M,I)
          SYSM(KI)=SYSM(KI)+AREA12*2.+DIFXX(M)*B(M,I)*B(M,I)+
     1               DIFYY(M)*A(M,I)*A(M,I)
          SYSF(KI)=SYSF(KI)+AREA12*2.*CONSTITP(KI,LEV,ICONSTIT)
          DO 430 J=1,3
            KJ=ICON(M,J)
            IF(J .EQ. I)THEN
              SYSF(KI)=SYSF(KI)-(DIFXX(M)*B(M,I)*B(M,I)
     1             +DIFYY(M)*A(M,I)*A(M,I))*CONSTITP(KI,LEV,ICONSTIT)
            ELSE
              SYSF(KI)=SYSF(KI)-(DIFXX(M)*B(M,I)*B(M,J)
     1             +DIFYY(M)*A(M,I)*A(M,J))*(CONSTITP(KJ,LEV,ICONSTIT)
     2             +CONSTITIT(KJ))
     3             +AREA12*(CONSTITP(KJ,LEV,ICONSTIT)-CONSTITIT(KJ))
            ENDIF
 430      CONTINUE
 420    CONTINUE
 410  CONTINUE
      CALL SOLVS
      CALL STPRCC
      CALL NORM
      IF(ANORM .GT. EPS .AND. IT .LT. MAXIT)GOTO 450
      RETURN
      END
      SUBROUTINE QVOLUM(VOL,X)
      INCLUDE 'ALL.CMM'
      DIMENSION X(1)
      VOL=0.
      DO 10 M=1,NMEL
      TEMPV=0.
      DO 20 I=1,3
      TEMPV=TEMPV+X(ICON(M,I))
 20   CONTINUE
      VOL=VOL+TEMPV*AREA(M)*.3333333
 10   CONTINUE
      RETURN
      END
C===========================================C
C      ROTATION OF DIFFUSION COEF. TENSOR   C
C===========================================C
      SUBROUTINE RTDIFT(NMNPQD,N)
      INCLUDE 'ALL.CMM'
      COMMON/DIFCOF/DISPCF(MXNDSL,2,2),DIF1(MXNDSL,MXLVLS),
     1  DIF2(MXNDSL,MXLVLS),DIFANG(MXNDSL,MXLVLS)
      DO 10 I=1,NMNPQD
      ANG=DIFANG(I,N)
      CANG=COS(ANG)
      SANG=SIN(ANG)
C * * DIF1 IS LONGITUDINAL DIF. COEFF, DIF2 IS LATERAL
C * * DIFANG = ANGLE FROM X AXIS (COUNTERCLKW) TO POSITIVE LONG AX
      DISPCF(I,1,1)=DIF1(I,N)*CANG*CANG+DIF2(I,N)*SANG*SANG
      DISPCF(I,1,2)=(DIF2(I,N)-DIF1(I,N))*CANG*SANG
      DISPCF(I,2,1)=DISPCF(I,1,2)
      DISPCF(I,2,2)=DIF1(I,N)*SANG*SANG+DIF2(I,N)*CANG*CANG
C      WRITE(*,*) DISPCF(I,1,1),DISPCF(I,1,2),DISPCF(I,2,1),DISPCF(I,2,2)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE SFLUX(I)
C * Boundary fluxes for vertical diffusion sol.
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      BTFLX(I)=0.
      TPFLX(I)=0.
      RHS(1)=RHS(1)+BTFLX(I)
      RHS(NLEVAD)=RHS(NLEVAD)+TPFLX(I)
      RETURN
      END
      SUBROUTINE SOLVC(I)
C *  Vertical diffusion sol.  See MAKERH for comment about generalization
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      CALL TRIDSOL
      DO 10 LEV=1,NLEV
        CONSTIT(I,LEV,ICONSTIT)=RHS(LEV)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE STCBC
      INCLUDE 'ALL.CMM'
      CONSTIT(NK,LEV,ICSALT)=36.0
      RETURN
      END
      SUBROUTINE STPRCC
C * Sets prescribed concentrations for diffusion step
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
C * * Only specify conc if coming in.  Going out, assume diff flux=0
C     If advection is done, CONSTITBND are already set
C     else call CONSTITBC here
      DO 10 NK=1,NMNP
C * * If NDSTSL is neg. it contains the nearest boundary node
C        (which could be itself, in which case an external BC is required)
        IF(NDSTSL(NK) .LT. 0)THEN
          CONSTIT(NK,LEV,ICONSTIT)=CONSTITP(-NDSTSL(NK),LEV,ICONSTIT)
          IF(NK .EQ. -NDSTSL(NK))CONSTIT(NK,LEV,ICONSTIT)=CONSTITBND(NK)
        ENDIF
 10   CONTINUE
C * * Set prescribed conc
C *CHECK* Need to handle different constituents
      DO 20 I=1,NMPBN(LEV)
        NK=NB(I,LEV)
        CONSTIT(NK,LEV,ICONSTIT)=PRCONC(I,LEV)
 20   CONTINUE
      RETURN
      END
      SUBROUTINE STUPAD
      INCLUDE 'ALL.CMM'
      CHARACTER*1 ACHAR
      CHARACTER*80 TITLE
      SAVE NRCDSK 
      WRITE(2,1002)
 1002 FORMAT('ADVECTION-DIFFUSION ACTIVATED'/)
C * * For now make prescribed nodes at all levels the same.
      READ(9,*)NMPBN(1),IDFVAR,ISWAD,NDISK,ICONC1,ILOOK,
     1  IHOT,IHALF
      DO 70 LEV=2,NLEVAD
        NMPBN(LEV)=NMPBN(1)
 70   CONTINUE
C================================================================
C      NOTE THAT INTERNAL & EXTERNAL NODE NUMBERS MUST BE IDENTICAL
C
C================================================================
C
C     trying 4 byte word size for recl by mult NMNP*4
      NMNPX4=NMNP*4
      OPEN(16,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 FORM='UNFORMATTED',RECL=NMNPX4)
      NRCDSK=0
      READ(9,1004)TITLE
 1004 FORMAT(/A80)
      WRITE(2,1006)TITLE
 1006 FORMAT(A80)
      WRITE(2,1008)NMPBN(1),IDFVAR,ISWAD,NDISK,ICONC1,
     1  ILOOK,IHOT,IHALF
 1008 FORMAT(10X,
     1  'NUMBER OF PRESCRIBED CONCENTRATION NODES, NMPBN = ',I4/,10X,
     2  'VARIABLE DIFFUSION COEFFICIENTS, IDFVAR = ',I3/,10X,
     3  'SWITCH, ISWAD = ',I2,' (=1 => ADV. ONLY; =2 => ADV.+DIF;',
     3  ' =3 => DIF ONLY)'/,
     4  10X,'STORE CONCENTRATIONS IN DISK FILE AT',I6,'TIME STEPS'/,
     5  10X,'ICONC1 = 0 => CONST CONC.;  ICONC1 = 1, INPUT VARIABLE',
     6  ' CONC.; ICONC1 = ',I3/,10X,'ILOOK =',I3,' =1:W, =2:R/W'
     7  /,10X,'IHOT = 0(COLD); = 1(HOT)',I4/
     8  ,10X,'IHALF (CREATE U&V AT HALF STEP BY INTERP)',I3)
C * * Define the horizontal diffusion coefficients
      CALL DEFDIF
      IADTM=0
      READ(9,*)NADOUT,MAXIT,EPS
C * * Use 0 for NADOUT to let NOUT determine.  Else NADOUT is used
C     with IADTM to determine when to write sal to IUOUT.
C * * The TINCAD used here is the time step between two velocity fields
C     Thus TINCAD would typically be TWICE the time difference between
C     velocity records.  In the coupled model, velocities are kept
C     for previous time step, half time step and present time step
C     The half time step is a .5*multad steps.  Don't use ihalf=1
      IF(NADOUT .EQ. 0)THEN
        IF(NOUT/MULTAD*MULTAD .NE. NOUT)THEN
          WRITE(*,*)'WARNING NOUT NOT MULTAD MULTIPLE. Hydro and sal'
          WRITE(*,*)'  will not be written to IUOUT at the same times.'
          WRITE(*,*)'ENTER ANY CHARACTER TO CONTINUE OR CTRL C TO STOP'
          READ(*,'(A1)')ACHAR
        ENDIF
        NADOUT=NOUT/MULTAD
        IF(NADOUT .LT. 1)NADOUT=1
      ENDIF
      DTD2AD=.5*TINCAD
      WRITE(2,1010)TINCAD,NADOUT,MAXIT
 1010 FORMAT(10X,
     2  'TIME STEP OF FULL 4TH ORD RK INT, TINCAD = ',F12.4/,10X,
     3  'OUTPUT CONTROL, NADOUT = ',I8
     4  /,'ITERATION CONTROL, MAXIT =',I4)
C============================================C
C     INITIAL CONDITION FOR TEMP & SALINITY  C
C============================================C
      IF(IRDHOT .EQ. 0)THEN
        READ(9,*)SAL1
        CALL INITTS(SAL1)
      ENDIF
C * If salt is solved for and affects density, ICSALT > 0, update density
      IF(ICSALT .GT. 0)THEN
        DO LEV=1,NLEVAD
          DO I=1,NMNP
            DENSTY(I,LEV)=CONSTITP(I,LEV,ICSALT)*0.7
          ENDDO
        ENDDO
      ELSE
        DENSTY=1000.
      ENDIF
      CALL DNSTRM
C**** READ IN PRESCRIBED TEMP/SAL , IF ANY
      IF(NMPBN(LEV) .GT. 0)THEN
        READ(9,*)(NB(I,LEV),PRCONC(I,LEV), I=1,NMPBN(LEV))
      ENDIF
C****  IF HEAT/SAL FLUX SWITCH IHSFLX=1 , THEN
C**** READ IN PRESCRIBED HEAT/SAL FLUXES AT THE SURFACE & BOTTOM
C****   IF (IHSFLX.NE.1) THEY ARE ASSUMED ZERO.
C      IF(IHSFLX.EQ.1) THEN
C****  1ST HEAT FLX AT BOTTOM
C         READ(9,*) (HFLXBT(I),I=1,NMNP)
C****  HEAT FLX AT SURFACE
C         READ(9,*) (HFLXTP(I),I=1,NMNP)
C****  SALINITY FLX AT BOTTOM
C         READ(9,*) (SFLXBT(I),I=1,NMNP)
C****  SALINITY FLX AT SURFACE
C         READ(9,*) (SFLXTP(I),I=1,NMNP)
C      ELSE
C        DO 720  I=1,NMNP
C           HFLXBT(I)=0.0
C           HFLXTP(I)=0.0
C           SFLXBT(I)=0.0
C           SFLXTP(I)=0.0
C 720    CONTINUE
C      ENDIF
      WRITE(2,1016)CONSTITP(1,1,1),CONSTITP(NMNP,1,1)
 1016 FORMAT('INITIAL CONC ',4E13.5)
      IF(ILOOK .GT. 1)THEN
      DO 40 N=1,NMNP
        READ(8,1018)NK,LOOKUP(NK,1),(LOOKUP(NK,J), J=2,LOOKUP(NK,1)+1)
 1018   FORMAT(I6,15I5)
 40   CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE TRCBK
C * 050900 2-D and 3-D ignoring vertical motion
C * 060401 3-D version
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
C *  In RK first index is variable (u,v,w) and second is RK step -1
      DIMENSION RK(3,3)
C The level is transferred in LEV of common CGRID
      DO 10 NK=1,NMNP
        IF(NLEVAD .EQ. 1)THEN
          UN=UT(NK)
          VN=VT(NK)
        ELSE
          UN=U(NK,LEV)
          VN=V(NK,LEV)
C * W is vertical vel relative to local sigma level and in sigma units.
          WN=W(NK,LEV)
c           wn=0.
        ENDIF
C * * First step of RK complete.
        DO 60 NRK=2,4
C *  NRK COUNTS THE R-K STEP. Each step consists of specifying
C    the location and then determining the velocity at that loc.
C    Determine location based on previous step
          XN=XORD(NK)-UN*DTD2AD
          YN=YORD(NK)-VN*DTD2AD
          ZN=SIGMA(LEV)-WN*DTD2AD
          IF(NLEVAD .GT. 1)THEN
C * Check to see if above surf or below bott.
            IF(ZN .LE. -1.)THEN
              ZN=-0.99999
            ELSEIF(ZN .GE. 0.)THEN
              ZN=-0.00001
            ENDIF
C   LEVZN is the sigma level immediately below zn. For bott levzn=1
C     For surf levzn=nlev-1
            DO 20 ILEV=1,NLEVAD
              IF(ZN .GT. SIGMA(ILEV))GOTO 20
              LEVZN=ILEV-1
              GOTO 30
 20         CONTINUE
          ELSE
            LEVZN=1
          ENDIF
C * Nlocat only locates the horizontal position.
 30       CALL NLOCAT(IONSIDE)
C *  If node left domain, cannot track anymore, NDSTSL=-nearest node
C    and xol,yol are set. Jump to next node, so take care of ZOL for
C    these nodes (ndstsl < 0, and ionside =1)
          IF(NDSTSL(NK) .LT. 0.)THEN
            ZOL(NK)=SIGMA(LEV)
            GOTO 10
          ENDIF
          IF(IONSIDE .EQ. 1)THEN
            ZOL(NK)=ZN
            GOTO 10
          ENDIF
          IE=NDSTSL(NK)
          A1=((XORD(ICON(IE,3))-XORD(ICON(IE,2)))*(YN-YORD(ICON(IE,2))))
     1      -((XN-XORD(ICON(IE,2)))*(YORD(ICON(IE,3))-YORD(ICON(IE,2))))
          A2=((XORD(ICON(IE,1))-XORD(ICON(IE,3)))*(YN-YORD(ICON(IE,3))))
     1      -((XN-XORD(ICON(IE,3)))*(YORD(ICON(IE,1))-YORD(ICON(IE,3))))
          A3=((XORD(ICON(IE,2))-XORD(ICON(IE,1)))*(YN-YORD(ICON(IE,1))))
     1      -((XN-XORD(ICON(IE,1)))*(YORD(ICON(IE,2))-YORD(ICON(IE,1))))
          EA2=AREA(IE)*2.
          XI1=A1/EA2
          XI2=A2/EA2
          XI3=A3/EA2
          IF(NRK .LT. 4)THEN
C * L and U refer to lower and upper sigma levels
            UNL=XI1*UTMHDT(ICON(IE,1),LEVZN)+
     1         XI2*UTMHDT(ICON(IE,2),LEVZN)+
     2         XI3*UTMHDT(ICON(IE,3),LEVZN)
            VNL=XI1*VTMHDT(ICON(IE,1),LEVZN)+
     1         XI2*VTMHDT(ICON(IE,2),LEVZN)+
     2         XI3*VTMHDT(ICON(IE,3),LEVZN)
            WNL=XI1*WTMHDT(ICON(IE,1),LEVZN)+
     1         XI2*WTMHDT(ICON(IE,2),LEVZN)+
     2         XI3*WTMHDT(ICON(IE,3),LEVZN)
            UNU=XI1*UTMHDT(ICON(IE,1),LEVZN+1)+
     1         XI2*UTMHDT(ICON(IE,2),LEVZN+1)+
     2         XI3*UTMHDT(ICON(IE,3),LEVZN+1)
            VNU=XI1*VTMHDT(ICON(IE,1),LEVZN+1)+
     1         XI2*VTMHDT(ICON(IE,2),LEVZN+1)+
     2         XI3*VTMHDT(ICON(IE,3),LEVZN+1)
            WNU=XI1*WTMHDT(ICON(IE,1),LEVZN+1)+
     1         XI2*WTMHDT(ICON(IE,2),LEVZN+1)+
     2         XI3*WTMHDT(ICON(IE,3),LEVZN+1)
            ZIGINT=(ZN-SIGMA(LEVZN))/DSIGMA(LEVZN)
            RK(1,NRK-1)=UNL+(UNU-UNL)*ZIGINT
            RK(2,NRK-1)=VNL+(VNU-VNL)*ZIGINT
            RK(3,NRK-1)=WNL+(WNU-WNL)*ZIGINT
c             rk(3,nrk-1)=0.
            IF(NRK .EQ. 3)THEN
C *  Apply 2 factor to compensate for using 1/2 dt in formula
              UN=2.*RK(1,NRK-1)
              VN=2.*RK(2,NRK-1)
              WN=2.*RK(3,NRK-1)
            ELSE
              UN=RK(1,NRK-1)
              VN=RK(2,NRK-1)
              WN=RK(3,NRK-1)
            ENDIF
          ELSE
C * Final RK step
            UNL=XI1*UTMDT(ICON(IE,1),LEVZN)+
     1        XI2*UTMDT(ICON(IE,2),LEVZN)+
     2        XI3*UTMDT(ICON(IE,3),LEVZN)
            VNL=XI1*VTMDT(ICON(IE,1),LEVZN)+
     1        XI2*VTMDT(ICON(IE,2),LEVZN)+
     2        XI3*VTMDT(ICON(IE,3),LEVZN)
            WNL=XI1*WTMDT(ICON(IE,1),LEVZN)+
     1        XI2*WTMDT(ICON(IE,2),LEVZN)+
     2        XI3*WTMDT(ICON(IE,3),LEVZN)
            UNU=XI1*UTMDT(ICON(IE,1),LEVZN+1)+
     1        XI2*UTMDT(ICON(IE,2),LEVZN+1)+
     2        XI3*UTMDT(ICON(IE,3),LEVZN+1)
            VNU=XI1*VTMDT(ICON(IE,1),LEVZN+1)+
     1        XI2*VTMDT(ICON(IE,2),LEVZN+1)+
     2        XI3*VTMDT(ICON(IE,3),LEVZN+1)
            WNU=XI1*WTMDT(ICON(IE,1),LEVZN+1)+
     1        XI2*WTMDT(ICON(IE,2),LEVZN+1)+
     2        XI3*WTMDT(ICON(IE,3),LEVZN+1)
            ZIGINT=(ZN-SIGMA(LEVZN))/DSIGMA(LEVZN)
            RK(1,3)=UNL+(UNU-UNL)*ZIGINT
            RK(2,3)=VNL+(VNU-VNL)*ZIGINT
            RK(3,3)=WNL+(WNU-WNL)*ZIGINT
c             rk(3,3)=0.
            IF(NLEVAD .EQ. 1)THEN
              XN=XORD(NK)-TINCAD*(UT(NK)+2.*(RK(1,1)+RK(1,2))
     1          +RK(1,3))/6.
              YN=YORD(NK)-TINCAD*(VT(NK)+2.*(RK(2,1)+RK(2,2))
     1          +RK(2,3))/6.
            ELSE
              XN=XORD(NK)-TINCAD*(U(NK,LEV)+2.*(RK(1,1)+RK(1,2))
     1          +RK(1,3))/6.
              YN=YORD(NK)-TINCAD*(V(NK,LEV)+2.*(RK(2,1)+RK(2,2))
     1          +RK(2,3))/6.
C * Should not be necessary to do for lev=nlev since W would be zero
C   But do it anyway just to simplify things
              ZN=SIGMA(LEV)-TINCAD*(W(NK,LEV)+2.*(RK(3,1)+RK(3,2))
     1          +RK(3,3))/6.
            ENDIF
            CALL NLOCAT(IONSIDE)
            IF(NDSTSL(NK) .GT. 0 .AND. IONSIDE .EQ. 0)THEN
C * Check to see if above surf or below bott.
              IF(ZN .LE. -1.)THEN
                ZN=-0.99999
              ELSEIF(ZN .GE. 0.)THEN
                ZN=-0.00001
              ENDIF
              ZOL(NK)=ZN
              XOL(NK)=XN
              YOL(NK)=YN
            ELSE
              ZOL(NK)=SIGMA(LEV)
            ENDIF
          ENDIF
 60     CONTINUE
 10   CONTINUE
C * * UPON RETURN THE TRACEBACK LOCATION XOL,YOL,ZOL IS KNOWN. NDSTSL CONTAINS
C     ELEMENT WHERE THIS IS OR A NEG. NODE NUMBER (NEAREST BOUNDARY NODE).
      RETURN
      END
      SUBROUTINE TRIDS
      INCLUDE 'ALL.CMM'
      DIMENSION COFMAT(MXLVLS+1,3),ESUM1(MXLVLS),
     1  ESUM2(MXLVLS)
      DSUM=TINC/(DSIGMA(1)+DSIGMA(2))
      COFMAT(1,1)=.0
      COFMAT(1,2)=1.
      COFMAT(1,3)=.0
      COFMAT(2,1)=0.
      COFMAT(NLEV,3)=0.
C * * EDZZ has the mixing length (ie z) in it.
C * 021502 Using the Mellor Yamada turb cl, EDZZ is just EDZZ, so
C   ABSUST should be one.
      ESUM=0.5*ABSUST(NK)*(EDZZ(NK,2)+EDZZ(NK,3))/DSIGMA(2)
      ESUM3=0.5*TINC*ABSUST(NK)*(EDZZ(NK,NLEV)+EDZZ(NK,NLEV-1))/
     1  DSIGMA(NLEV-1)**2.
        DO 10 LEV=3,NLEV-1
          ESUM1(LEV)=0.5*ABSUST(NK)*(EDZZ(NK,LEV)+EDZZ(NK,LEV-1))/
     1      DSIGMA(LEV-1)/(DSIGMA(LEV-1)+DSIGMA(LEV))
          ESUM2(LEV)=0.5*ABSUST(NK)*(EDZZ(NK,LEV)+EDZZ(NK,LEV+1))/
     1      DSIGMA(LEV)/(DSIGMA(LEV-1)+DSIGMA(LEV))
 10     CONTINUE
C  041601 Chg so that first sigma level is at z0
C      Z=HTAV(NK)*DSIGMA(1)+Z0
      Z=HTAV(NK)*DSIGMA(1)
C *042001* Keep Ks/H (relative roughness) less than 0.05. z0=Ks/30
      ALOGZZ0=ALOG(Z/Z0)
      IF(ALOGZZ0 .LT. 3.)ALOGZZ0=3.
      SPD2=SQRT(QPREV(2*NK-1,2)*QPREV(2*NK-1,2)
     1  +QPREV(2*NK,2)*QPREV(2*NK,2))
      COFMAT(2,2)=1.+DSUM*ESUM/HTAV(NK)**2.
     1          + DSUM/HTAV(NK)*(VKRMN/ALOGZZ0)**2.*SPD2
      COFMAT(2,3)=-DSUM*ESUM/HTAV(NK)**2.
      DO 20 LEV=3,NLEV-1
        COFMAT(LEV,1)=-TINC*ESUM1(LEV)/HTAV(NK)**2.
        COFMAT(LEV,2)=1.+TINC*(ESUM1(LEV)+ESUM2(LEV))/HTAV(NK)**2.
        COFMAT(LEV,3)=-TINC*ESUM2(LEV)/HTAV(NK)**2.
 20   CONTINUE
      COFMAT(NLEV,1)=-ESUM3/HTAV(NK)**2.
      COFMAT(NLEV,2)=1.-COFMAT(NLEV,1)
      DO 30 LEV=1,NLEV-1
        FA=COFMAT(LEV+1,1)/COFMAT(LEV,2)
        COFMAT(LEV+1,2)=COFMAT(LEV+1,2)-FA*COFMAT(LEV,3)
        DTAUXDZ(LEV+1)=DTAUXDZ(LEV+1)-FA*DTAUXDZ(LEV)
        DTAUYDZ(LEV+1)=DTAUYDZ(LEV+1)-FA*DTAUYDZ(LEV)
 30   CONTINUE
      DTAUXDZ(NLEV)=DTAUXDZ(NLEV)/COFMAT(NLEV,2)
      DTAUYDZ(NLEV)=DTAUYDZ(NLEV)/COFMAT(NLEV,2)
      DO 40 LEV=1,NLEV-1
        DTAUXDZ(NLEV-LEV)=(DTAUXDZ(NLEV-LEV)-COFMAT(NLEV-LEV,3)*
     1                         DTAUXDZ(NLEV-LEV+1))/COFMAT(NLEV-LEV,2)
        DTAUYDZ(NLEV-LEV)=(DTAUYDZ(NLEV-LEV)-COFMAT(NLEV-LEV,3)*
     1                         DTAUYDZ(NLEV-LEV+1))/COFMAT(NLEV-LEV,2)
 40   CONTINUE
      DO 210 LEV=2,NLEV
        Q(2*NK-1,LEV)=DTAUXDZ(LEV)
        Q(2*NK,LEV)=DTAUYDZ(LEV)
 210  CONTINUE
      RETURN
      END
      SUBROUTINE TRIDSOL
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      RHS(1)=RHS(1)/SYSMVT(1,2)
      DO 10 LEV=1,NLEVM1
C * * DIAGONALS ARE MADE INTO ONES
      SYSMVT(LEV,3)=SYSMVT(LEV,3)/SYSMVT(LEV,2)
      RHS(LEV+1)=RHS(LEV+1)-SYSMVT(LEV+1,1)*RHS(LEV)
      SYSMVT(LEV+1,2)=SYSMVT(LEV+1,2)-SYSMVT(LEV+1,1)*SYSMVT(LEV,3)
      RHS(LEV+1)=RHS(LEV+1)/SYSMVT(LEV+1,2)
 10   CONTINUE
      DO 20 LEV=NLEVM1,1,-1
      RHS(LEV)=RHS(LEV)-SYSMVT(LEV,3)*RHS(LEV+1)
 20   CONTINUE
      RETURN
      END
C**** 1-D MASS DIFFUSION IN THE VERTICAL (ALONG SIGMA LEVEL AT ANY NODE)
      SUBROUTINE VDIFUS
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      DO 10 I=1,NMNP
       DO 20 LEV=1,NLEVAD
        RHS(LEV)=0.0
       DO 20 K=1,3
        SYSMVT(LEV,K)=0.0
 20    CONTINUE
        CALL MAKEM(I,HT(I))
        CALL MAKERH(I,HT(I))
        CALL SFLUX(I)
        CALL SOLVC(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE WRTEXT
      INCLUDE 'ALL.CMM'
      IF(NMHBN .GT. 0)CALL STETAB
      WRITE(IUOUT,1002) TIME,TINC,ITIME
 1002 FORMAT(10X,'TIME = ',F15.1,' SEC',2X,'DELTA T WAS, TINC = ',
     1  F8.2,1X,'SEC.',2X,'TIME STEP, ITIME = ',I8/10X,
     2  'NET VOLUME ABOVE MLW, VOL = ',E13.6)
      CALL VOLUME
      CALL SLINE(15)
      WRITE(IUOUT,1004)
 1004 FORMAT(1X,'INT',3X,'EXT',5X,'H',11X,'QX',10X,'QY',10X,
     1  'ETA',10X,'U',11X,'V',9X,'SYSFH',8X,'SYSFQ',8X,'SYSFQ'/,2X,
     2  'NODE',6X,'NODE'/,1X,'NUMBER',4X,'NUMBER'/)
      DO 10 I=1,NMNP
        WRITE(IUOUT,1006)I,NEXT(I),HT(I),QT(2*I-1),QT(2*I),ETA(I),
     1    UT(I),VT(I),SYSFH(I),SYSFQ(2*I-1),SYSFQ(2*I)
 1006   FORMAT(2X,I4,2X,I4,2X,6(E12.6,1X),2(E12.5,1X),E12.5)
 10   CONTINUE
      CALL SLINE(36)
      RETURN
      END
      SUBROUTINE NTPTAU
C * * Interpolate the bottom stress from the internal mode for use
C       in external mode solution.  The stress lags by an internal time step.
      INCLUDE 'ALL.CMM'
      NT=ITIME-ITIME/IMULT*IMULT
      IF(NT .EQ. 0)NT=IMULT
      ANT=NT
      ANT=ANT/IMULT
      DO 10 I=1,NMNP
        TAUBGX(I)=TAUSTX(I,1)*(1.-ANT)+ANT*TAUSTX(I,2)
        TAUBGY(I)=TAUSTY(I,1)*(1.-ANT)+ANT*TAUSTY(I,2)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE NLOCAT(IONSIDE)
c *  Given point coordinates xn,yn find which element it is in.
C    xn,yn and nk are in common
      INCLUDE 'ALL.CMM'
      SAVE NRCDSK
      REAL*8 DDDIST,DDXN,DDYN,DDX1,DDX2,DDY1,DDY2
C *  First check lookup table for this node, nk
      IONSIDE=0
      IF(LOOKUP(NK,1) .GT. 0)THEN
        DO 210 LK=1,LOOKUP(NK,1)
          IEL=LOOKUP(NK,LK+1)
          DO 220 J=1,3
            K1=ICON(IEL,J)
            IF(J .EQ. 3)THEN
              K2=ICON(IEL,1)
            ELSE
              K2=ICON(IEL,J+1)
            ENDIF
C * *  Single precision can cause problems when xn,yn is close to a node.
            DDXN=XN
            DDX1=XORD(K1)
            DDY1=YORD(K1)
            DDY2=YORD(K2)
            DDYN=YN
            DDX2=XORD(K2)
            DDDIST=(DDXN-DDX1)*(DDY1-DDY2)+(DDYN-DDY1)*(DDX2-DDX1)
            IF(DDDIST .LT. -0.01) GO TO 210
 220      CONTINUE
C * * FALL THRU MEANS PARTICLE FOUND IN THIS ELEMENT
          NDSTSL(NK)=IEL
C * * UPDATE LOOKUP TABLE BY CYCLING ELEMENTS.
          IF(LK .NE. 1)THEN
            NLKTMP=LOOKUP(NK,LK+1)
            CALL PPUSH(LK)
            LOOKUP(NK,2)=NLKTMP
          ENDIF
c *  Point found in lookup table element, status set, so return
          RETURN
 210    CONTINUE
      ENDIF
C * * If control passes to here, point is not in lookup table.
C     Start an element search.  Simple check on nodal coordinates
      DO 80 IEL=1,NMEL
        DO 90 J=1,3
          XJ=XORD(ICON(IEL,J))
          IF(ABS(XJ-XN) .GT. 1.1*SRCHRDX(IEL))GOTO 80
          YJ=YORD(ICON(IEL,J))
          IF(ABS(YJ-YN) .GT. 1.1*SRCHRDY(IEL))GOTO 80
 90     CONTINUE
C * * FALL THRU MEANS ALL DIFFS. ARE .LE. SRCHRD. CHECK THIS ELEMENT
        DO 320 J=1,3
          K1=ICON(IEL,J)
          IF(J .EQ. 3)THEN
            K2=ICON(IEL,1)
          ELSE
            K2=ICON(IEL,J+1)
          ENDIF
C * *  Single precision can cause problems when xn,yn is close to a node.
          DDXN=XN
          DDX1=XORD(K1)
          DDY1=YORD(K1)
          DDY2=YORD(K2)
          DDYN=YN
          DDX2=XORD(K2)
          DDDIST=(DDXN-DDX1)*(DDY1-DDY2)+(DDYN-DDY1)*(DDX2-DDX1)
C          DIST=((XN-XORD(K1))*(YORD(K1)-YORD(K2)))+((YN-YORD(K1))*
C     1      (XORD(K2)-XORD(K1)))
          IF(DDDIST .LT. -0.01)GO TO 80
 320    CONTINUE
C * * FALL THRU MEANS POSITION IN THIS ELEMENT
        NDSTSL(NK)=IEL
C * * ADD ELEMENT TO SEARCH ARAY
        LOOKUP(NK,1)=LOOKUP(NK,1)+1
        IF(LOOKUP(NK,1) .GT. MXLKUP-1)THEN
          WRITE(*,*)'TOO MANY LOOKUP ELEMENTS ',NK,IEL,XN,YN
          STOP
        ENDIF
        CALL PUSH
        LOOKUP(NK,2)=IEL
C        WRITE(*,*)NK,NDSTSL(NK),LOOKUP(NK,1),LOOKUP(NK,2)
C *  Point located in element, status saved, lookup table extended.
C    So return
        RETURN
 80   CONTINUE
C * * FALL THRU MEANS PARTICLE NOT IN MODEL AREA.  IF CODE IS 0,2, OR 3
C     FIND NEAREST BOUNDARY NODE AND USE ITS CONC. SET NDSTSL TO MINUS NEAREST.
      IF(NBC(NK) .EQ. 2 .OR. NBC(NK) .EQ. 3 .OR. NBC(NK) .EQ. 0
     1  .OR. NBC(NK) .EQ. 14 .OR. NBC(NK) .EQ. 11)THEN
        CALL FNDNBN
        XOL(NK)=XORD(NDNEAR)
        YOL(NK)=YORD(NDNEAR)
C * * If land boundary (1 or 9) then force on side.
      ELSEIF(NBC(NK) .EQ. 1 .OR. NBC(NK) .EQ. 9 .OR.
     1  NBC(NK) .EQ. 19)THEN
        CALL ONSIDE
C * Shortcut, if traced to side don't do Runge-Kutta, just get location
C    and go to next node
        IONSIDE=1
C * * If 4 or 10 node force to node itself.
      ELSEIF(NBC(NK) .EQ. 4 .OR. NBC(NK) .EQ. 10)THEN
        XOL(NK)=XORD(NK)
        YOL(NK)=YORD(NK)
        NDSTSL(NK)=-NK
      ENDIF
      RETURN
      END

