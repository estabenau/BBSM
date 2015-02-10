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
C * * These subroutines are set up for 00:00UTC 1996 as time 0.
      SUBROUTINE DEFDIF
      INCLUDE 'ALL.CMM'
      INTEGER*4 ICNMTH(8)
      DATA ICNMTH/552,1860,2083,2303,2425,2801,2992,3189/
C * * Here we just define constant diff coefficients *DT/(8*AREA)
C * * To use variable, define DIFXX and DIFYY here.
      TINC8=TINCAD/8.
      READ(9,*)DXX,DYY,VDIF(1,1),HDIFFUS,HUPPERLIM
      WRITE(IUOUT,1002)DXX,DYY,VDIF(1,1),HDIFFUS,HUPPERLIM
 1002 FORMAT('Spatially constant diffusion coefficients are chosen'
     1  /,'     DIFXX, DIFYY = ',2F10.2/,
     2    '     Vertical diff. coeff. = ',F12.5/,
     3  '     Depth at which dxx and dyy apply, HDIFFUS = ',F8.2/,
     4  '     Depth at which diff no longer depend on depth, ',
     5  ' HUPPERLIM = ',F8.2)
C * Scale diffusion to hdiffus, square function
      DXX=DXX/(HDIFFUS*HDIFFUS)
      DYY=DYY/(HDIFFUS*HDIFFUS)
      DO 10 I=1,NMEL
        HSUM=0.
        DO J=1,3
          IJ=ICON(I,J)
          HSUM=HSUM+HT(IJ)
        ENDDO
        HSUM=HSUM/3.
        IF(HSUM .GT. HUPPERLIM)HSUM=HUPPERLIM
C * Make dif coef depth dependent. square function
        DIFXX(I)=DXX*HSUM*HSUM/AREA(I)*TINC8
        DIFYY(I)=DYY*HSUM*HSUM/AREA(I)*TINC8
 10   CONTINUE
C**** READ IN THE VERTICAL DIFFUSIVITIES
C      IF(IVDFVAR.EQ.0) THEN
C**** IVDFVAR=0 => THE DIF COEF ARE ALL EQUAL
        DO 210 I=1,NMNP
          DO 210 LEV=1,NLEV
             VDIF(LEV,I)=VDIF(1,1)
 210       CONTINUE
C      ELSE
C       IF(IVDFVAR.EQ.1) THEN
C**** IVDFVAR=1 => ALL NODAL DIF COEF ARE EQUAL AT ANY LEVEL BUT
C**** CHANGE WITH SIGMA LEVEL
C        READ(9,1052) (VDIF(1,LEV), LEV=1,NLEV)
C 1052  FORMAT(8F10.2)
C        DO 220 LEV=1,NLEV
C           DO 230 I=2,NMNP
C              VDIF(LEV,I)=VDIF(LEV,1)
C 230        CONTINUE
C 220     CONTINUE
C       ELSE
C**** IVDFVAR=2 => THE DIF COEF ARE ALL UNEQUAL
C        DO 240 LEV=1,NLEV
C         READ(9,*) (VDIF(LEV,I), I=1,NMNP)
C 240    CONTINUE
C       ENDIF
C      ENDIF
C      DO 20 I=1,NMEL
C        XCNTRD=0.
C        YCNTRD=0.
C        DO 30 J=1,3
C          NJ=ICON(I,J)
C          XCNTRD=XCNTRD+XORD(NJ)
C          YCNTRD=YCNTRD+YORD(NJ)
C 30     CONTINUE
C        XCNTRD=0.3333*XCNTRD
C        YCNTRD=0.3333*YCNTRD
C        DO 40 IC=1,8
C          XREF=XORD(ICNMTH(IC))
C          YREF=YORD(ICNMTH(IC))
C          RADIUS=SQRT((XCNTRD-XREF)*(XCNTRD-XREF)+
C     1      (YCNTRD-YREF)*(YCNTRD-YREF))
C          IF(RADIUS .LT. 3500.)THEN
C            DFCOF=14.*(1.-RADIUS/3500.)+1.
C            IF(DIFXX(I) .LT. DFCOF)THEN
C              DIFXX(I)=DFCOF/AREA(I)*TINC8
C              DIFYY(I)=DFCOF/AREA(I)*TINC8
C            ENDIF
C          ENDIF
C 40     CONTINUE
C 20   CONTINUE
      RETURN
      END
      SUBROUTINE ELINFO
      INCLUDE 'ALL.CMM'
      DO 15 L=1,NMEL
      READ(1,1007)NELM(L),(ICON(L,J),J=1,3),CF(L),EDXX(L),EDYY(L),
     1  EDXY(L)
      NILM(NELM(L))=L
 1007 FORMAT(4I5,4F10.0)
 15   CONTINUE
C * * SET VERTICAL EDDY VISCOSITY, edzz=ustar*edyvis*z z is dist from bot
C     ustar is calculated and multipl later.
      DO 47 LEV=1,NLEV
      DO 44 I=1,NMNP
C      EDZZ(I,LEV)=EDYVIS*(1.+SIGMA(LEV))*(-DEPTH(I))
C * 021502 For Mellor Yamada, set kinematic visc.
       EDZZ(I,LEV)=EDYVIS
  44  CONTINUE
  47  CONTINUE
C*** IFRIC=0 => bottom stress from internal mode, set CF=1
C***      =2 => from external mode using cf(1)
C***      =1 => from external mode using variable cf
C    In order to preserve compatibility with 3D the CF must become
C    a nodal quantity.
      CFREF=CF(1)
      IF(IFRIC .EQ. 2)THEN
C * * End loop at NMEL just in case CF is element quantity
        DO 120 I=2,NMEL
          CF(I)=CF(1)
 120    CONTINUE
      ELSEIF(IFRIC .EQ. 0)THEN
C * * 3D bottom stress determined from internal mode, so set CF=1
        DO 500 I=1,NMEL
          CF(I)=1.
 500    CONTINUE
      ENDIF
 135  IF (IEDVIS .EQ. 2) GO TO 200
      DO 210 I=2,NMEL
      EDXX(I)=EDXX(1)
      EDYY(I)=EDYY(1)
      EDXY(I)=EDXY(1)
 210  CONTINUE
 200  IF (IPRINT.EQ.0) GO TO 55
      WRITE(IUOUT,1022)
 1022 FORMAT(5X,'ELEMENT CONNECTIVITIES:'/'ELEMENT #',
     1  2X,'NODE 1',2X,'NODE 2',2X,'NODE 3',3X,'EDDY XX',5X,'EDDY YY',
     3  5X,'EDDY XY',5X,'DRAG COEF'//)
      DO 220 I=1,NMEL
      WRITE(IUOUT,1024)I,(ICON(I,J),J=1,3),EDXX(I),EDYY(I),EDXY(I),CF(I)
 1024 FORMAT(I5,4X,I5,3X,I5,3X,I5,4X,F7.2,5X,F7.2,5X,
     1  F7.2,5X,F7.3)
 220  CONTINUE
 55   RETURN
      END
      SUBROUTINE INITTS(SAL1)
      INCLUDE 'ALL.CMM'
      DO 10 LEV=1,NLEV
      DO 20 I=1,NMNP
        CONSTITP(I,LEV,1)=SAL1
C        TEMPRV(I,LEV)=TEMP1
c        IF(YORD(I) .LT. 16903.3)THEN
c          IF(YORD(I) .GT. 10258.)THEN
c            DIST=YORD(I)-10258.
c            DISTN=DIST/(20000.-10258.)
c            CONSTITP(I,LEV,1)=DISTN*SAL1+(1.-DISTN)*19.
c          ELSE
c            CONSTITP(I,LEV,1)=19.
c          ENDIF
c        ENDIF
 20   CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CONSTITBC
C * * Set prescribed salinities at boundaries.
C     Input file 19 contains, time [day] and salinity [ppt]
      INCLUDE 'ALL.CMM'
      COMMON/CSALBC/TIMESF(194),SFSAL(194),TIMECC(4283),CCSAL(4283),
     1  TIMERC(194),RCSAL(194),TIMEBC(194),BCSAL(194)
      INTEGER ICOUNTRC,ICOUNTBC,ICOUNTSF,ICOUNTCC
      SAVE ICOUNTRC,ICOUNTBC,ICOUNTSF,ICOUNTCC
      IF(IADTM .EQ. 1)THEN
       OPEN(19,FILE='./Input/Salinity/SFSalinity.dat',
     1     STATUS='OLD')
        DO I=1,194
          READ(19,*,IOSTAT=IOS)TIMESF(I),SFSAL(I)
          IF(IOS .NE. 0)EXIT
        ENDDO
        CLOSE(19)
       OPEN(19,FILE='./Input/Salinity/CCSalinity.dat',
     1     STATUS='OLD')
        DO I=1,4283
          READ(19,*,IOSTAT=IOS)TIMECC(I),CCSAL(I)
          IF(IOS .NE. 0)EXIT
        ENDDO
        CLOSE(19)
       OPEN(19,FILE='./Input/Salinity/BB26Salinity.dat',
     1    STATUS='OLD')
        DO I=1,194
          READ(19,*,IOSTAT=IOS)TIMERC(I),RCSAL(I)
          IF(IOS .NE. 0)EXIT
        ENDDO
        CLOSE(19)
       OPEN(19,FILE='./Input/Salinity/BB28Salinity.dat',
     1    STATUS='OLD')
        DO I=1,194
          READ(19,*,IOSTAT=IOS)TIMEBC(I),BCSAL(I)
          IF(IOS .NE. 0)EXIT
        ENDDO
        CLOSE(19)
        ICOUNTSF=1
        ICOUNTBC=1
        ICOUNTRC=1
        ICOUNTCC=1
      ENDIF
C * * Convert time in sec to days.
      DAY=TIME/86400.+1.
 10   IF(DAY .GT. TIMERC(ICOUNTRC))THEN
        ICOUNTRC=ICOUNTRC+1
        GOTO 10
      ENDIF
 20   IF(DAY .GT. TIMEBC(ICOUNTBC))THEN
        ICOUNTBC=ICOUNTBC+1
        GOTO 20
      ENDIF
 30   IF(DAY .GT. TIMESF(ICOUNTSF))THEN
        ICOUNTSF=ICOUNTSF+1
        GOTO 30
      ENDIF
 40   IF(DAY .GT. TIMECC(ICOUNTCC))THEN
        ICOUNTCC=ICOUNTCC+1
        GOTO 40
      ENDIF
      SALBRC=RCSAL(ICOUNTRC-1)+(DAY-TIMERC(ICOUNTRC-1))/
     1  (TIMERC(ICOUNTRC)-TIMERC(ICOUNTRC-1))*
     2  (RCSAL(ICOUNTRC)-RCSAL(ICOUNTRC-1))
      SALBBC=BCSAL(ICOUNTBC-1)+(DAY-TIMEBC(ICOUNTBC-1))/
     1  (TIMEBC(ICOUNTBC)-TIMEBC(ICOUNTBC-1))*
     2  (BCSAL(ICOUNTBC)-BCSAL(ICOUNTBC-1))
      SALBSF=SFSAL(ICOUNTSF-1)+(DAY-TIMESF(ICOUNTSF-1))/
     1  (TIMESF(ICOUNTSF)-TIMESF(ICOUNTSF-1))*
     2  (SFSAL(ICOUNTSF)-SFSAL(ICOUNTSF-1))
      SALBCC=CCSAL(ICOUNTCC-1)+(DAY-TIMECC(ICOUNTCC-1))/
     1  (TIMECC(ICOUNTCC)-TIMECC(ICOUNTCC-1))*
     2  (CCSAL(ICOUNTCC)-CCSAL(ICOUNTCC-1))
C * * Rickenbacker, nodes 3595,3600,1148,1147,1146
      CONSTITBND(3595)=SALBRC
      CONSTITBND(3600)=SALBRC
      CONSTITBND(1148)=SALBRC
      CONSTITBND(1147)=SALBRC
      CONSTITBND(1146)=SALBRC
C * * Bear Cut, nodes 1126,1125,1124
C        SALBN3=33.80+0.83*COS(1.99107579D-07*TIME-1.86401164D+00)
      CONSTITBND(1126)=SALBBC
      CONSTITBND(1125)=SALBBC
      CONSTITBND(1124)=SALBBC
C * * Safety Valve 1073-1099 , Boca Chita 3679, Sands Cut 1070
      DO I=1073,1099
        CONSTITBND(I)=SALBSF
      ENDDO
      CONSTITBND(3679)=SALBSF
      CONSTITBND(1070)=SALBSF
C * * ABCs 1013 1016-1018, 1032-1034,
      CONSTITBND(1013)=SALBCC
      CONSTITBND(1016)=SALBCC
      CONSTITBND(1017)=SALBCC
      CONSTITBND(1018)=SALBCC
      CONSTITBND(1032)=SALBCC
      CONSTITBND(1033)=SALBCC
      CONSTITBND(1034)=SALBCC
      RETURN
      END
      SUBROUTINE SOURCE
C * * This routine should only be called for 1-layer solutions
      INCLUDE 'ALL.CMM'
C************************************************
C * * The arrays QSTR contain date field and 10 discharges or withdrawals
C     This is an update from v2 to v3 of the model
C     QSTR is current value in m/sec calculated in STRCFL
      INTEGER*4 ISTR2Q(10)
      DATA ISTR2Q/3,4,5,6,9,8,10,2,7,11/
      DO I=1,NMSTRC
        DO K=1,NUMSTRCELEM(I)
          IEL=ISRCEL(K,I)
          HEL=0.
          DO J=1,3
            HEL=HT(ICON(IEL,J))+HEL
          ENDDO
          DO  J=1,3
            KJ=ICON(IEL,J)
C * The structure order and flow data input order are messed up
C   Use ISTR2Q to fix.  QSTR here is really Q/AREAC calc in STRCFL
            FLOWIN=QSTR(ISTR2Q(I))*AREA(IEL)/AREAND(KJ)
C   v3 update If FLOWIN is negative then we are looking at loss of water from system
C       not a gain from a canal flow.  In that case, salinty is not recalculated
            IF(FLOWIN .GT. 0)THEN
              CONSTIT(KJ,LEV,ICONSTIT)=CONSTIT(KJ,LEV,ICONSTIT)/
     1         (1.+FLOWIN*TINCAD*3./HEL)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SRC3D
C * * This routine is used to specify mass source input when
C     using characteristic advection model.
C      INCLUDE 'ALL.CMM'
      RETURN
      END
      SUBROUTINE STETAB
      INCLUDE 'ALL.CMM'
C * * ROUTINE TO CALCULATE REAL TIDE IN BISCAYNE BAY
C * * BASED ON HARMONIC ANALYSIS OF VIRGINIA KEY DATA
C * * Harmonic constants for t=0 at 00:00GMT Jan 1, 1996
C  See ..\input\boundaryTide.xls
      DATA ACOSM2,PHASM2,OMEGM2/0.298,4.9,1.40518903E-04/,
     1  ACOSN2,PHASN2,OMEGN2/0.062,0.46,1.37879699E-04/,
     2  ACOSK1,PHASK1,OMEGK1/0.027,4.46,7.29211584E-05/,
     3  ACOSO1,PHASO1,OMEGO1/0.025,3.16,6.75977442E-05/,
     4  ACOSS2,PHASS2,OMEGS2/0.049,1.23,1.45444104E-04/,
     5  ACOSK2,PHASK2,OMEGK2/0.011,3.83,1.45842317E-04/
C  START OF 1-COS SECTION.
      IF(IRDHOT .EQ. 0 .AND. (TIME-STRTIM) .LT. WTRAMP)THEN
        RAMP=0.5*(1.-COS(3.1415926/WTRAMP*(TIME-STRTIM)))
      ELSE
        RAMP=1.
      ENDIF
C      DO 20 I=1,NMHBN
C      I1=NHN(I)
C      ETA(I1)=HB(I)*COS(PHASE*(TIME-ALAG(I)))
      VK=ACOSM2*COS(OMEGM2*TIME-PHASM2)
     1     +ACOSN2*COS(OMEGN2*TIME-PHASN2)
     2     +ACOSO1*COS(OMEGO1*TIME-PHASO1)
     3     +ACOSK1*COS(OMEGK1*TIME-PHASK1)
     4     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     5     +ACOSK2*COS(OMEGK2*TIME-PHASK2)
      ETA(1126)=VK*RAMP
      ETA(1125)=ETA(1126)
      ETA(1124)=ETA(1126)
      RC=1.02*ACOSM2*COS(OMEGM2*(TIME-.69*3600.)-PHASM2)
     1     +1.06*ACOSN2*COS(OMEGN2*(TIME-1.3*3600.)-PHASN2)
     2     +1.06*ACOSO1*COS(OMEGO1*(TIME-1.94*3600.)-PHASO1)
     3     +0.83*ACOSK1*COS(OMEGK1*TIME-PHASK1)
     4     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     5     +ACOSK2*COS(OMEGK2*TIME-PHASK2)
      ETA(3595)=RC*RAMP
      ETA(3600)=ETA(3595)
      ETA(1148)=ETA(3600)
      ETA(1147)=ETA(1148)
      ETA(1146)=ETA(1147)
      CAPEFL=0.98*ACOSM2*COS(OMEGM2*(TIME+.345*3600.)-PHASM2)
     1     +1.03*ACOSN2*COS(OMEGN2*(TIME+.387*3600.)-PHASN2)
     2     +1.08*ACOSO1*COS(OMEGO1*(TIME-.27*3600.)-PHASO1)
     3     +1.01*ACOSK1*COS(OMEGK1*(TIME-1.38*3600.)-PHASK1)
     4     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     5     +ACOSK2*COS(OMEGK2*TIME-PHASK2)
      SK=1.02*ACOSM2*COS(OMEGM2*(TIME+.04*3600.)-PHASM2)
     1     +1.07*ACOSN2*COS(OMEGN2*(TIME+.09*3600.)-PHASN2)
     2     +1.14*ACOSO1*COS(OMEGO1*(TIME+.18*3600.)-PHASO1)
     3     +1.03*ACOSK1*COS(OMEGK1*(TIME+.88*3600.)-PHASK1)
     4     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     5     +ACOSK2*COS(OMEGK2*TIME-PHASK2)
       ETA(1099)=CAPEFL*RAMP
       ETA(1083)=SK*RAMP
       ETA(1084)=ETA(1083)
       ETA(1098)=14./15.*ETA(1099)+1./15.*ETA(1084)
       ETA(1097)=13./15.*ETA(1099)+2./15.*ETA(1084)
       ETA(1096)=12./15.*ETA(1099)+3./15.*ETA(1084)
       ETA(1095)=11./15.*ETA(1099)+4./15.*ETA(1084)
       ETA(1094)=10./15.*ETA(1099)+5./15.*ETA(1084)
       ETA(1093)=9./15.*ETA(1099)+6./15.*ETA(1084)
       ETA(1092)=8./15.*ETA(1099)+7./15.*ETA(1084)
       ETA(1091)=7./15.*ETA(1099)+8/15.*ETA(1084)
       ETA(1090)=6./15.*ETA(1099)+9./15.*ETA(1084)
       ETA(1089)=5./15.*ETA(1099)+10/15.*ETA(1084)
       ETA(1088)=4./15.*ETA(1099)+11./15.*ETA(1084)
       ETA(1087)=3./15.*ETA(1099)+12./15.*ETA(1084)
       ETA(1086)=2./15.*ETA(1099)+13./15.*ETA(1084)
       ETA(1085)=1./15.*ETA(1099)+14./15.*ETA(1084)
      BC=0.76*ACOSM2*COS(OMEGM2*(TIME-.75*3600.)-PHASM2)
     1     +0.75*ACOSN2*COS(OMEGN2*(TIME-.71*3600.)-PHASN2)
     2     +0.9*ACOSO1*COS(OMEGO1*(TIME-1.17*3600.)-PHASO1)
     3     +0.7*ACOSK1*COS(OMEGK1*(TIME-.01*3600.)-PHASK1)
     4     +0.75*ACOSS2*COS(OMEGS2*(TIME-.66*3600.)-PHASS2)
     5     +0.75*ACOSK2*COS(OMEGK2*(TIME-.92*3600.)-PHASK2)
      ETA(1070)=BC*RAMP
      ETA(1073)=10./13.*ETA(1070)+3./13.*ETA(1083)
      ETA(1074)= 9./13.*ETA(1070)+4./13.*ETA(1083)
      ETA(1075)= 8./13.*ETA(1070)+5./13.*ETA(1083)
      ETA(1076)= 7./13.*ETA(1070)+6./13.*ETA(1083)
      ETA(1077)= 6./13.*ETA(1070)+7./13.*ETA(1083)
      ETA(1078)= 5./13.*ETA(1070)+8./13.*ETA(1083)
      ETA(1079)= 4./13.*ETA(1070)+9./13.*ETA(1083)
      ETA(1080)= 3./13.*ETA(1070)+10./13.*ETA(1083)
      ETA(1081)= 2./13.*ETA(1070)+11./13.*ETA(1083)
      ETA(1082)= 1./13.*ETA(1070)+12./13.*ETA(1083)
C      ETA(3679)=ETA(1070)-0.03*(YORD(1060)-YORD(1125))/50000.
      ETA(3679)=ETA(1070)
C * * Caesar
      ETA(1033)=0.5*(ACOSM2*COS(OMEGM2*(TIME-3600.)-PHASM2)
     1     +ACOSN2*COS(OMEGN2*(TIME-3600.)-PHASN2)
     2     +ACOSO1*COS(OMEGO1*(TIME-3600.)-PHASO1)
     3     +ACOSK1*COS(OMEGK1*TIME-PHASK1)
     4     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     5     +ACOSK2*COS(OMEGK2*TIME-PHASK2))*RAMP
      ETA(1034)=ETA(1033)
      ETA(1032)=ETA(1033)
C * * Broad and Angelfish both use 1018
      ETA(1018)=0.45*(ACOSM2*COS(OMEGM2*(TIME-3800.)-PHASM2)
     2     +ACOSN2*COS(OMEGN2*(TIME-3800.)-PHASN2)
     4     +ACOSO1*COS(OMEGO1*(TIME-3800.)-PHASO1)
     6     +ACOSK1*COS(OMEGK1*TIME-PHASK1)
     8     +ACOSS2*COS(OMEGS2*TIME-PHASS2)
     A     +ACOSK2*COS(OMEGK2*TIME-PHASK2))*RAMP
      ETA(1017)=ETA(1018)
      ETA(1016)=ETA(1018)
      ETA(1013)=ETA(1018)
      RETURN
      END
      SUBROUTINE STEVAP
C * * Set evaporation in m/sec, could be different for each element
      INCLUDE 'ALL.CMM'
      SAVE DAY1,DAY2,EV1,EV2
      IF(ITIME .EQ. 1)THEN
        OPEN(71,FILE='./Input/RAINMEVAP_run_file.dat',STATUS='OLD')
        READ(71,*)DAY1,EV1
C * * P-E values are in m/day,  so convert to m/sec. Evap is normally neg.
        EV1=EV1/86400.
        READ(71,*)DAY2,EV2
        EV2=EV2/86400.
      ENDIF
C * * Convert time in sec to Julian days.
      DAY=TIME/86400.+1.
 10   IF(DAY .GT. DAY2)THEN
        EV1=EV2
        DAY1=DAY2
        READ(71,*)DAY2,EV2
        EV2=EV2/86400.
        GOTO 10
      ELSE
        EVAP=EV1+(DAY-DAY1)/(DAY2-DAY1)*(EV2-EV1)
      ENDIF
      RETURN
      END
      SUBROUTINE STRCFL
      INCLUDE 'ALL.CMM'
C * * The arrays QSTR contain day and NMSTRC structure discharges
C     QSTR is current value, QSTR1 is previous daily value,
C     QSTR is changed to m/sec in this routine by dividing by element areas.
C     QSTR2 is next daily value. m3/sec.
C     ISRCEL(j,i) contains the elements for structure i in which
C     the flow is distributed. AREAC is the area of those elements
C     NOTE order of the columns
      INTEGER*4 ISTR2Q(10)
      DATA ISTR2Q/3,4,5,6,9,8,10,2,7,11/
      DATA ISWITCH/1/
      SAVE NMSTRCP1,ISWITCH
      IF(ISWITCH .EQ. 1)THEN
        ISRCEL=0
       OPEN(921,FILE='./Input/GRIDV3/CANAL_ELEMENTSv3.DAT',
     1    FORM='FORMATTED')
        READ(921,*)NMSTRC
C * Number of structures. Add one to allow for time stamp
        NMSTRCP1=NMSTRC+1
        DO I=1,NMSTRC
          READ(921,*)NUMSTRCELEM(I)
          DO J=1,NUMSTRCELEM(I)
            READ(921,*)ISRCEL(J,I)
          ENDDO
        ENDDO
        CLOSE(921)
       OPEN(7,FILE='./Input/Canal_flow_run_file.dat',
     1    STATUS='OLD')
C * Read 1 header line.
        READ(7,*)
        READ(7,1001)(QSTR1(J),J=1,NMSTRCP1)
 1001   FORMAT(11F8.2)
        READ(7,1001)(QSTR2(J),J=1,NMSTRCP1)
C * Find the area for each structure discharge, careful, discharges
C     and isrcl are in messed up order. Fixed with istr2q which gives
C     q input order as fct of strcture numbers 1 thru 9, eg.
C     structure 7 is the 10th input q.
        AREAC=0.
        DO I=1,NMSTRC
          DO J=1,NUMSTRCELEM(I)
            AREAC(ISTR2Q(I))=AREAC(ISTR2Q(I))+AREA(ISRCEL(J,I))
          ENDDO
        ENDDO
        ISWITCH=0
      ENDIF
 10   DAY=TIME/86400.+1.
      IF(DAY .GT. QSTR2(1))THEN
        DO 20 J=1,NMSTRCP1
          QSTR1(J)=QSTR2(J)
 20     CONTINUE
        READ(7,1001)(QSTR2(J),J=1,NMSTRCP1)
        GOTO 10
      ELSE
        DO 30 J=2,NMSTRCP1
        QSTR(J)=QSTR1(J)+(DAY-QSTR1(1))/(QSTR2(1)-QSTR1(1))
     1    *(QSTR2(J)-QSTR1(J))
C * Make qstr m/sec spread out over areac. Areac has been put in
C     q input order
        QSTR(J)=QSTR(J)/AREAC(J)
 30     CONTINUE
      ENDIF
      DO I=1,NMSTRC
        DO J=1,NUMSTRCELEM(I)
          IEL=ISRCEL(J,I)
          EVAR=AREA(IEL)/3.
          DO K=1,3
            KJ=ICON(IEL,K)
C * Convert qstr from q input order to structure order
            SYSFH(KJ)=SYSFH(KJ)+QSTR(ISTR2Q(I))*EVAR
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE VARWND
      INCLUDE 'ALL.CMM'
      SAVE DAY1,WU1,WV1,DAY2,WU2,WV2
      IF(ITIME .EQ. 1)THEN
C * * x and y components of wind vector (blowing to). m/s, utc
        OPEN(4,FILE='./Input/VKWindUV_run_file.dat',STATUS='OLD')
        READ(4,1003)DAY1,WU1,WV1
 1003   FORMAT(F10.4,F8.2,F8.2)
        WRITE(IUOUT,*)'DAY1',DAY1,'WU1',WU1,'WV1',WV1
        READ(4,1003)DAY2,WU2,WV2
      ENDIF
C * * Convert time in sec to Julian days.
 10   DAY=TIME/86400.+1.
      IF(DAY .GT. DAY2)THEN
        WU1=WU2
        WV1=WV2
        DAY1=DAY2
        READ(4,1003)DAY2,WU2,WV2
        GOTO 10
      ELSE
        WU=WU1+(DAY-DAY1)/(DAY2-DAY1)*(WU2-WU1)
        WV=WV1+(DAY-DAY1)/(DAY2-DAY1)*(WV2-WV1)
        WINDSP=SQRT(WU*WU+WV*WV)
        IF(WU .EQ. 0.)WU=.01
        WDIRC=ATAN2(WV,WU)
      ENDIF
      CALL WNDGLO
      RETURN
      END
      SUBROUTINE WNDINT
      INCLUDE 'ALL.CMM'
      DO 10 I=1,NMNP
        WNDA(I)=1.
        WNDB(I)=0.
 10   CONTINUE
      RETURN
      END
      SUBROUTINE WNDGLO
      INCLUDE 'ALL.CMM'
C * *   THIS GETS WIND STRESS (TAUWG X/Y) IN GLOBAL COORDS FOR USE IN FRC2DQ
C * *   REFDEN IS THE REFERENCE DENSITY
      RR=1.2/REFDEN
C      CD=(1.1+0.0536*WINDSP)*0.001
C LARGE & POND
C Section to modify drag if windspeed is 'high'
      IF(WINDSP .GT. 11)THEN
        CD=(0.49+0.065*WINDSP)*0.001
      ELSE
        CD=1.2*0.001
      ENDIF
C REDUCED ROUGHNESS IN BISCAYNE BAY
      CD=CD*0.5
C      CD=CD*0.4
C Further reduced during several storm events here - work-around to a problem
      if(time .gt. 273844800. .and. time .lt. 274060800.)cd=cd*0.8
      if(time .gt. 309571200. .and. time .lt. 309787200.)cd=cd*.3
      TAUWGX(1)=RR*CD*WINDSP*WINDSP*COS(WDIRC)
      TAUWGY(1)=RR*CD*WINDSP*WINDSP*SIN(WDIRC)
      IF(IRDHOT .EQ. 0 .AND. (TIME-STRTIM) .LT. WTRAMP)THEN
       TAUWGX(1)=TAUWGX(1)*(TIME-STRTIM)/WTRAMP
       TAUWGY(1)=TAUWGY(1)*(TIME-STRTIM)/WTRAMP
      ENDIF
      DO 10 I=2,NMNP
      TAUWGX(I)=TAUWGX(1)
      TAUWGY(I)=TAUWGY(1)
  10  CONTINUE
      RETURN
      END
      SUBROUTINE GNSWAT
      INCLUDE 'ALL.CMM'
      DATA ISWITCH/1/
      SAVE ISWITCH,NMGNSZONESP1
      IF(ISWITCH .EQ. 1)THEN
       OPEN(921,FILE='./Input/GRIDV3/GNSWCERPELEM.DAT',STATUS='OLD')
        IGSWEL=0
        READ(921,*)NMGNSZONES
C * Number of gns zones. Add one to allow for time stamp
        NMGNSZONESP1=NMGNSZONES+1
        DO I=1,NMGNSZONES
          READ(921,*)NGNSELEM(I)
          DO J=1,NGNSELEM(I)
            READ(921,*)IGSWEL(J,I)
          ENDDO
        ENDDO
        CLOSE(921)
       OPEN(17,FILE='./Input/scenarios/CERPGnSW_run_file.dat',
     1 STATUS='OLD')
C * Skip one header line
        READ(17,*)
        READ(17,1001)(GNSW1(J), J=1,NMGNSZONESP1)
 1001   FORMAT(F7.2,22F7.2)
        READ(17,1001)(GNSW2(J), J=1,NMGNSZONESP1)
C       FIXX
	print*,'last surface flow=',GNSW1(23),GNSW2(23)
C * * Compute areas.
        ACELL=0
        DO I=1,NMGNSZONES
          DO J=1,NGNSELEM(I)
C * First element is for time
            ACELL(I+1)=ACELL(I+1)+AREA(IGSWEL(J,I))
          ENDDO
        ENDDO
        ISWITCH=0
      ENDIF
C * * The time in the GnSW.DAT file is in YEAR day UTC re 1/1/96.
C     Use step function (not linear interp) to preserve volumes.
 10   DAY=TIME/86400.+1.
      IF(DAY .GT. GNSW2(1))THEN
        DO J=1,NMGNSZONESP1
          GNSW1(J)=GNSW2(J)
        ENDDO
        READ(17,1001)(GNSW2(J),J=1,NMGNSZONESP1)
        GOTO 10
      ENDIF
C * * 17 zones corresponding to the WMM cells, 5 zone to CERP,YD is in first elem.
C * * Each zone consists of NEL finite elements, given in IGSWEL
      DO I=1,NMGNSZONES
        GNSW(I+1)=GNSW1(I+1)/ACELL(I+1)
        DO J=1,NGNSELEM(I)
          IEL=IGSWEL(J,I)
          EVAR=AREA(IEL)/3.
          DO K=1,3
            KJ=ICON(IEL,K)
C * * SOLVX multiplies SYSFH by TINC
            SYSFH(KJ)=SYSFH(KJ)+GNSW(I+1)*EVAR
          ENDDO
        ENDDO
      ENDDO
      END
      SUBROUTINE GNSSAL
      INCLUDE 'ALL.CMM'
      DO I=1,NMGNSZONES
        DO K=1,NGNSELEM(I)
          IEL=IGSWEL(K,I)
          HEL=0.
          DO J=1,3
            HEL=HT(ICON(IEL,J))+HEL
          ENDDO
          DO J=1,3
            KJ=ICON(IEL,J)
            FLOWIN=GNSW(I+1)*AREA(IEL)/AREAND(KJ)
            CONSTIT(KJ,LEV,ICONSTIT)=CONSTIT(KJ,LEV,ICONSTIT)/
     1        (1.+FLOWIN*TINCAD*3./HEL)
C	  if (i .eq. 22) print*,'check last surface flow',i,flowin
C	  if (i .eq. 22) print*,CONSTIT(KJ,LEV,ICONSTIT),hel,kj,area(iel)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
