      SUBROUTINE IMICRO
c
c  Changes by Martin Koehler for 
c  I41: Ice Cloud Decay Sensitivity Experiment (ICDS-41)
c  * all thresholds neglecting small mixing ratio cases: * 10^-2
c  - threshold in any condensate to run microphysics: 10^-8 (normally 10^-6)
c  - threshold in each condensate to define its existance: 5*10^-8 (5*10^-6)
c  - other thresholds (CRIT, QVK/QS...-1) * 10^-2  (increased accuracy)
c
c *************************************************************************
c MODIFICATIONS BELOW ARE MADE TO LORDS ORIGINAL CODE (Q. Fu & S. Krueger)
c**************************************************************************
c	Updates (Jun. 8, 1993, cfu comments):
c	  1. RMINUC = QI * RHO / NC;  this  will give the depositional
c            growth of cloud ice at expense of cloud water  properly (
c            PIDW).
c         2. Replacing R50 by R100; then DT1 is the growth time needed
c            for ice crystal to grow from 40 to 100um (PSFI and PSFW).
c         3. Setting Psfi = 0.0 when QL is less than 5e-6 (PSFI).
c         4. A1 and A2 are obtained by using interpolation; A1T(14) is
c            0.1725e-9 for function IDW ( PSFI, PSFW and PIDW).
c         5. Setting QI0 = .6e-3 (threshold for cloud ice aggregation)
c            and QS0 = 1.0e-3 (threshold for snow aggregation ) (PSAUT
c            and PGAUT) (private communication with S. Chin).
c         6. GAM625 = 184.860962 (PWACS)
c         7. Esi and Egs ( T < T0 ) are set to be  0.1 following  RH84 
c            ( PSACI and PGACS ).
c**********************************************************************
c       Updates (Jun. 9, 1993, crh comments):
c         1. Modifying intercept parameters and densities.
c         2. Replacing Lin's hail fall speed relation with an average
c            graupel relation (Ug=40.74*Dg**0.5) based on Locatelli &
c            Hobbs (1974).
c**********************************************************************
c       Updates (Jun. 10, 1993, clin comments):
c         1. Correcting C2BRG.
c         2. Egs = 1.0 for T > T0 and wet growth. Note:  modifications
c            here should be consistent to (7) under the comment 'cfu'.
c**********************************************************************
c       Updates (Jun. 10, 1993, csk comments):
c         1. Lin et al. (1983) version is  installed for the following
c            processes: PSMLT, PGWET, PGMLT (constants only).
c            Change DIMENSION CSMLT(3),CGWET(2) => CSMLT(5), CGWET(4).
c         2. PWACS = 0.0 following Lin et al. (1983).
c         3. And more.
c ************************************************************************
C
      COMMON/VTERM/VTS,VTG,VTR
      COMMON/MVAR/T,PRESR,QMIX(6),TSS(7),PSS(26),RHOK
      COMMON /STK/DT,DTL,DZ,TDT,TDZ,DZ2,DX,TDX,X0,ALFA,
     1  NND,IND,IBD,JND,JBD,JM1,Z1,ZB,RHOS,RMAX
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
C
      COMMON/VCON/VCONR,VCONS,VCONG
C
      COMMON/RHOCON/RHOFAC,RHO
C
      DIMENSION C3RACS(3)
     *,C3SACR(3),C3GACR(3),C3GACS(3)
C
      DIMENSION QRHO(6),QLOGIC(6),THRSH(6),SMAX(6)
C
      REAL IDW
C
      LOGICAL QLOGIC,VAPR,LIQ,ICE,SNO,GRAUP,RAIN
C
      EQUIVALENCE (QMIX(1),QVK),(QMIX(2),QLK),(QMIX(3),QIK),
     *  (QMIX(4),QSK),(QMIX(5),QGK),(QMIX(6),QRK),(QRHO(1),QVR),
     *  (QRHO(2),QLR),(QRHO(3),QIR),(QRHO(4),QSR),(QRHO(5),QGR),
     *  (QRHO(6),QRR),(QLOGIC(1),VAPR),(QLOGIC(2),LIQ),
     *  (QLOGIC(3),ICE),(QLOGIC(4),SNO),(QLOGIC(5),GRAUP),
     *  (QLOGIC(6),RAIN),(TSSV,TSS(1)),(TSSL,TSS(2)),(TSSI,TSS(3)),
     *  (TSSS,TSS(4)),(TSSG,TSS(5)),(TSSR,TSS(6)),(TTMP,TSS(7)),
     *  (SVMAX,SMAX(1)),(SLMAX,SMAX(2)),(SIMAX,SMAX(3)),
     *  (SSMAX,SMAX(4)),(SGMAX,SMAX(5)),(SRMAX,SMAX(6))
C
      EQUIVALENCE (PRAUT,PSS(1)),(PRACW,PSS(2)),(PRACS,PSS(3))
     *,(PSACW,PSS(4)),(PWACS,PSS(5)),(PGACW,PSS(6)),(PGMLT,PSS(7))
     *,(PSMLT,PSS(8)),(PREVP,PSS(9)),(PIACR,PSS(10)),(PSACR,PSS(11))
     *,(PGACR,PSS(12)),(PGFR,PSS(13)),(PGACS,PSS(14)),(PGAUT,PSS(15))
     *,(PGACI,PSS(16)),(PGWORD,PSS(17)),(PRACI,PSS(18)),(PSAUT,PSS(19))
     *,(PSACI,PSS(20)),(PSSUB,PSS(21)),(PGSUB,PSS(22)),(PSFW,PSS(23))
     *,(PSFI,PSS(24)),(PIDW,PSS(25)),(PGWET,PSS(26))
C
      EQUIVALENCE (C3RACS(1),ACCO(1,1)),(C3SACR(1),ACCO(1,2))
     *,(C3GACR(1),ACCO(1,3)),(C3GACS(1),ACCO(1,4))
C
cxmk  DATA THRSH/0.0,5*1.E-6/
      DATA THRSH/0.0,5*1.E-8/          !threshold / 100
C
C     ***************************
C     ***************************
C     ****                   ****
C     ****   PRELIMINARIES   ****
C     ****                   ****
C     ***************************
C     ***************************
C
C
C     LOCAL VARIABLES -- TEMP, DENSITY, TEMP. DEPENDENT EFFICIENCIES
C     AND SATURATION VARIABLES
C
      RHO = RHOK
      RHOFAC = SQRT(RHOS/RHO)
      QL0RHO = QL0 * RHO
      TC = T-TICE
      TSQ = T**2
C
Cbloss
      EXPT1 = AMIN1(EXP(0.09*amin1(TC,0.0)),1.0)
      EXPT2 = AMIN1(EXP(0.025*amin1(TC,0.0)),1.0)
C
cxpb  ESW   = ESATW(T)
cxpb  QSW   = EPS*ESW/(PRESR-ESW)
      QSW   = qsatw(T,PRESR)
      DQS   = QVK-QSW
      DQS0  = CES0/(PRESR-ES0)-QVK
C
C     ZERO SOURCE AND SINK TERMS
C
      DO 1 K=1,26
    1 PSS(K) = 0.0
C
C     DEFINE MIXING RATIOS GREATER THAN THRESHOLD FOR
C     ACCRETION PROCESSES AND MASS DENSITIES
C
C           1:  WATER VAPOR
C           2:  CLOUD (SUSPENDED LIQUID) WATER
C           3:  CLOUD ICE CRYSTALS
C           4:  SNOW
C           5:  GRAUPEL
C           6:  RAIN
C
      DO 10 K=1,6
      TSS(K) = 0.0
      QLOGIC(K) = QMIX(K) .GT. THRSH(K)
      SMAX(K) = QMIX(K)/TDT
   10 QRHO(K) = QMIX(K) * RHO
      TSS(7)  = 0.0
C
C     TERMINAL VELOCITIES
C
      IF(SNO)    VTS = VTRS(QSR)
      IF(GRAUP)  VTG = VTRG(QGR)
      IF(RAIN)   VTR = VTRR(QRR)
C
C
C     *********************************************
C     ****                                     ****
C     ****   PROCESSES CALLED INDEPENDENT OF   ****
C     ****   TEMPERATURE ARE:                  ****
C     ****                                     ****
C     ****        GACW    SACW    RAUT         ****
C     ****        RACW    RACS    GACS         ****
C     ****                REVP                 ****
C     ****                                     ****
C     *********************************************
C     *********************************************
C
      IF(.NOT. LIQ)  GO TO 150
      IF(.NOT. GRAUP)  GO TO 110
C
C     PGACW
C
      PGACW = AMIN1(ACR2(QLK,QGR,CGACW,0.875),SLMAX)
C
  110 IF(.NOT. SNO)  GO TO 120
C
C     PSACW
C
      PSACW = AMIN1(ACR1(QLK,QSR,CSACW,0.8125),SLMAX)
C
  120 IF(QLK .LE. QL0)  GO TO 130
C
C     PRAUT
C
      QEXRHO = QLR - QL0RHO
      PRAUT = AMIN1(RAUT(CRAUT,QEXRHO),SLMAX)
C
  130 IF(.NOT. RAIN)  GO TO 200
C
C     PRACW
C
      PRACW = AMIN1(ACR1(QLK,QRR,CRACW,0.95),SLMAX)
C
  150 IF(.NOT. RAIN)  GO TO 200
C
  170 IF(.NOT. SNO)  GO TO 290
C
C     PRACS
C
      PRACS = AMIN1(ACR3(VTR,VTS,QSR,QRR,CRACS,C3RACS),SSMAX)
C
      GO TO 210
C
  200 IF(.NOT. SNO)  GO TO 290
C
  210 IF(.NOT. GRAUP)  GO TO 290
C
C     PGACS
C
cfu ********************************************************************
cfu   CTGACS = CGACS*EXPT1
      CTGACS = CGACS*0.1
cfu ********************************************************************
      PGACS = AMIN1(ACR3(VTG,VTS,QSR,QGR,CTGACS,C3GACS),SSMAX)
C
  290 IF(QRK .EQ. 0.0 .OR. ICE .OR. DQS .GE. 0.0)  GO TO 300
C
C     PREVP
C
      SW = QVK/QSW
      PREVP = AMIN1(REVP(SW,TSQ,QSW,QRR,CREVP),SRMAX)
C
  300 IF(TC .LT. 0.0)  GO TO 400
C
C     ***********************************
C     ***********************************
C     ****                           ****
C     ****     TC >= 0 PROCESSES     ****
C     ****                           ****
C     ****    SMLT   WACS   GMLT     ****
C     ****                           ****
C     ***********************************
C     ***********************************
C
csk ********************************************************************
      IF(QSK .EQ. 0.0) GO TO 305
C
C     PSMLT
C
csk   PSMLT = AMIN1(SMLT(TC,QSR,CSMLT),SSMAX)
      IF ( SNO .AND. RAIN ) THEN
csk
csk   PSACR CALLED FOR SMLT HERE    (T < 0 PROCESS)
csk
      PSACR = AMIN1(ACR3(VTS,VTR,QRR,QSR,CSACR,C3SACR),SRMAX)
      END IF
csk
csk   PSMLT (Follow Lin et al. 1983)
csk
      PSMLT = AMIN1(SMLT(TC,DQS0,QSR,PSACW,PSACR,CSMLT),SSMAX)
      PSMLT = AMAX1(PSMLT,0.0)
      PSACR = 0. ! Not used except in PSMLT
csk ********************************************************************
C
  305 IF(.NOT. LIQ .OR. .NOT. SNO)  GO TO 310
C
C     PWACS
C
csk ******************************************************************
      PWACS = 0.0
csk   PWACS = AMIN1(ACR1(QLK,QSR,CWACS,1.5625),SSMAX)
csk ******************************************************************
C
C
  310 IF(.NOT. GRAUP .OR. .NOT. RAIN)  GO TO 320
C
C     GACR CALLED FOR GMLT HERE
C
      PGACR = AMIN1(ACR3(VTG,VTR,QRR,QGR,CGACR,C3GACR),SRMAX)
C
C     GMLT:  GACW HAS ALREADY BEEN CALLED IF APPROPRIATE
C            GUARD AGAINST NEGATIVE VALUES AT TEMP CLOSE TO 0 DEG C
C
  320 IF(QGK .EQ. 0.0)  GO TO 330
C
      PGMLT = AMIN1(GMLT(TC,DQS0,QGR,PGACW,PGACR,CGMLT),SGMAX)
      PGMLT = AMAX1(PGMLT,0.0)
csk ********************************************************************
      PGACR = 0. ! Not used except in PGMLT
csk ********************************************************************
C
clin *******************************************************************
      PGACS = AMIN1(10.*PGACS,SSMAX)
clin *******************************************************************
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****             T>=0             ****
C     ****                              ****
C     **************************************
C
  330 CONTINUE
C
      TSSV = PREVP
      TSSL = -(PRAUT+PRACW+PSACW+PGACW)
      TSSS = -(PGACS+PRACS+PWACS+PSMLT)
      TSSG = PGACS-PGMLT
      TSSR = PRAUT+PRACW+PRACS+PSACW+PWACS+PGACW+PSMLT+PGMLT-PREVP
      TTMP = (-HLTF*(PRACS+PWACS+PGMLT+PSMLT)-HLTC*PREVP)/CP
C
C     WRITE(6,111) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  111 FORMAT('0',30X,'T>=0  '//10X,'TSSV=',E16.8/
     *,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8,
     */10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
      TEM = 0.0
      DO 1112 K=1,6
 1112 TEM = TEM+TSS(K)
C     WRITE(6,1113)  TEM
 1113 FORMAT('0',10X,'TOTAL SOURCE AND SINK=',E16.8)
C
      GO TO 1000
C
C     ************************************
C     ************************************
C     ****                            ****
C     ****      TC < 0 PROCESSES      ****
C     ****                            ****
C     ****   IDW    BERGRN    RACI    ****
C     ****   GACI      SACI    SAUT   ****
C     ****   GFR       SACR    GACR   ****
C     ****   SSUB      GAUT    GSUB   ****
C     ****   WORD                     ****
C     ****   NOTE: IACR IS NOT USED   ****
C     ****                            ****
C     ************************************
C     ************************************
C
  400 IF(.NOT. LIQ)  GO TO 410
C
C     PIDW
C
cfu *******************************************************************
      PIDW = AMIN1(IDW(TC,qik),SLMAX)
cfu   PIDW = AMIN1(IDW(TC),SLMAX)
cfu *******************************************************************
C
  410 IF(.NOT. ICE)  GO TO 450
C
C     BERGERON PROCESSES -- PSFW AND PSFI
C
cfu *******************************************************************
      CALL BERGRN(TC,QLK,QIK,QLR,PSFW,PSFI)
cfu   CALL BERGRN(TC,QIK,QLR,PSFW,PSFI)
cfu *******************************************************************
      PSFW = AMIN1(PSFW,SLMAX)
      PSFI = AMIN1(PSFI,SIMAX)
C
      IF(.NOT. RAIN)  GO TO 420
C
C     PRACI
C
      PRACI = AMIN1(ACR1(QIK,QRR,CRACI,0.95),SIMAX)
C
  420 IF(.NOT. GRAUP)  GO TO 430
C
C     PGACI
C
      PGACI = AMIN1(ACR2(QIK,QGR,CGACI,0.875),SIMAX)
C
  430 IF(.NOT. SNO)  GO TO 440
C
C     PSACI
C
cfu ******************************************************************
cfu   QIKT = QIK*EXPT2
      QIKT = QIK*0.1
cfu ******************************************************************
      PSACI = AMIN1(ACR1(QIKT,QSR,CSACI,0.8125),SIMAX)
C
  440 IF(QIK .LE. QI0)  GO TO 450
C
C     PSAUT
C
      C = 0.1*EXPT2
      PSAUT = AMIN1(AUT(C,QIK,QI0),SIMAX)
C
  450 IF(.NOT. RAIN)  GO TO 470
C
C     PGFR
C
      PGFR = AMIN1(GFR(TC,QRR,CGFR),SRMAX)
C
      IF(.NOT. SNO)  GO TO 460
C
C     PSACR
C
      PSACR = AMIN1(ACR3(VTS,VTR,QRR,QSR,CSACR,C3SACR),SRMAX)
C
  460 IF(.NOT. GRAUP)  GO TO 470
C
C     PGACR
C
      PGACR = AMIN1(ACR3(VTG,VTR,QRR,QGR,CGACR,C3GACR),SRMAX)
C
cxpb  470 ESI = ESATI(T)
cxpb      QSI = EPS*ESI/(PRESR-ESI)
  470 QSI = qsati(T,PRESR)
      SI = QVK/QSI
C
      IF(QSK .EQ. 0.0)  GO TO 480
C
C     PSSUB: CAN BE EITHER POSITIVE OR NEGATIVE
C
      PSSUB = AMIN1(SSUB(SI,TSQ,QSI,QSR,CSSUB),SVMAX)
      PSSUB = AMAX1(PSSUB,-SSMAX)
C
      IF(QSK .LE. QS0)  GO TO 480
C
C     PGAUT
C
      C = 1.E-3*EXPT1
      PGAUT = AUT(C,QSK,QS0)
C
  480 IF(QGK .EQ. 0.0)  GO TO 520
      IF(QIK .NE. 0.0 .OR. QLK .NE. 0.0 .OR. QVK .GE. QSI)  GO TO 500
C
C     PGSUB: NEGATIVE VALUES ONLY
C
      PGSUB = AMAX1(GSUB(SI,TSQ,QSI,QGR,CGSUB),-SGMAX)
C
C     PGDRY OR PGWET
C
  500 PGDRY = PGACW + PGACI + PGACR + PGACS
      IF(.NOT.LIQ .AND. .NOT.RAIN)  GO TO 510
C
clin *******************************************************************
      PGWET = GWET(DQS0,TC,QGR,PGACI,PGACS,CGWET,SIMAX,SSMAX)
clin  PGWET = GWET(DQS0,TC,QGR,PGACI,PGACS,CGWET)
clin *******************************************************************
C
C
      IF(PGDRY .LE. PGWET)  GO TO 510
C
clin  PGWET INVOLVES REDEFINITION OF PGACI, PGACS (INCREASE OF
C     COLLECTION EFFICIENCY BY FACTOR OF 10), AND RECALCULATION
C     OF PGACR (WHICH CAN NOW BECOME NEGATIVE DUE TO SHEDDING
C     OF CLOUD LIQUID WATER WHICH IS CONVERTED INTO RAINWATER)
C
clin *******************************************************************
clin  PGACI = PGACI * 10.0
      PGACI = AMIN1(PGACI * 10.0,SIMAX)
      PGACS = AMIN1(10.*PGACS,SSMAX)
clin  PGACR = PGWET - PGACW - PGACI - PGACS
      PGACR = AMIN1(PGWET - PGACW - PGACI - PGACS,SRMAX)
clin *******************************************************************
C
      PGWORD = PGWET
      GO TO 520
C
  510 PGWORD = PGDRY
C
C
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****              T<0             ****
C     ****                              ****
C     **************************************
C
  520 CONTINUE
      TSSV = -(PSSUB+PGSUB)+PREVP
      TSSL = -(PSACW+PSFW+PRAUT+PRACW+PGACW+PIDW)
      TSSI = -(PSAUT+PSACI+PSFI+PRACI+PGACI)+PIDW
      TSSS = PSAUT+PSACI+PSACW+PSFW+PSFI+PSSUB-(PGACS+PGAUT)
C
      IF(QRK .LT. 1.E-4 .AND. QSK .LT. 1.E-4)  GO TO 530
      TSSS = TSSS-PRACS
      TSSG = PSACR+PRACS
      GO TO 540
  530 TSSS = TSSS+PSACR
  540 IF(QRK .GE. 1.E-4)  GO TO 550
      TSSS = TSSS+PRACI
      GO TO 560
  550 TSSG = PRACI+TSSG
  560 TSSG = TSSG+PGSUB+PGAUT+PGFR+PGWORD
      TSSR = PRAUT+PRACW-(PSACR+PGACR+PREVP+PGFR)
      TTMP = (HLTF*(PSACW+PGACW+PSACR+PGACR+PGFR+PSFW+PIDW)
     * -HLTC*PREVP+HLTS*(PSSUB+PGSUB))/CP
C
C     WRITE(6,666) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  666 FORMAT('1',30X,'T<0  '//10X,'TSSV=',E16.8/
     *,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8,
     */10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
      TEM = 0.0
      DO 6661 K=1,6
 6661 TEM = TEM+TSS(K)
C     WRITE(6,1113)  TEM
C
C
 1000 CONTINUE
C
C     PRINT VALUES FOR EACH PROCESS
C
C     WRITE(6,2000)'ALL TEMP','PGACW',PGACW,'PSACW',PSACW,'PRAUT',PRAUT,
C    1   'PRACW',PRACW,'PRACS',PRACS,'PGACS',PGACS,'PREVP',PREVP
C     WRITE(6,2000)'TC > 0','PSMLT',PSMLT,'PWACS',PWACS,'PGMLT',PGMLT
C     WRITE(6,2000) 'TC < 0','PSFW',PSFW,'PSFI',PSFI,'PRACI',PRACI,
C    1   'PIACR',PIACR,'PGACI',PGACI,'PSACI',PSACI,'PSAUT',PSAUT,
C    2   'PGFR',PGFR,'PSACR',PSACR,'PGACR',PGACR,'PSSUB',PSSUB,
C    3   'PGAUT',PGAUT,'PGSUB',PGSUB,'PIDW',PIDW
C     WRITE(6,2000)'WET-DRY','PGWET',PGWET,'PGDRY',PGDRY,'PGWORD',PGWORD
      RETURN
      ENTRY ERMCRO
      WRITE(6,*) '********* ERROR REPORT FROM MICRO *********'
      WRITE(6,666) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
C
C     PRINT VALUES FOR EACH PROCESS
C
      WRITE(6,2000)'ALL TEMP','PGACW',PGACW,'PSACW',PSACW,'PRAUT',PRAUT,
     1   'PRACW',PRACW,'PRACS',PRACS,'PGACS',PGACS,'PREVP',PREVP
      WRITE(6,2000)'TC > 0','PSMLT',PSMLT,'PWACS',PWACS,'PGMLT',PGMLT
      WRITE(6,2000) 'TC < 0','PSFW',PSFW,'PSFI',PSFI,'PRACI',PRACI,
     1   'PIACR',PIACR,'PGACI',PGACI,'PSACI',PSACI,'PSAUT',PSAUT,
     2   'PGFR',PGFR,'PSACR',PSACR,'PGACR',PGACR,'PSSUB',PSSUB,
     3   'PGAUT',PGAUT,'PGSUB',PGSUB,'PIDW',PIDW
      WRITE(6,2000)'WET-DRY','PGWET',PGWET,'PGDRY',PGDRY,'PGWORD',PGWORD
 2000 FORMAT(1H0,'PROCESSES:  ',A10/
     1   3(1H0,5(A6,'=',E12.5,4X)/))
      RETURN
      END
      SUBROUTINE SAT(JPOINT,KPOINT,T,P,QVK,QLK,QIK)
C
C     QS REDEFINED FOR QI = QL = 0 -- 12/2/86
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
      COMMON/PATHO/NOCONV,NCTOT
C **DEBUG
      COMMON/ERRDAT/NERR,IERR,JERR,PRESS,TO,QVO,QLO,QIO
C
C
C     *************************************************************
C     *************************************************************
C     ****                                                     ****
C     ****    PRELIMINARY ADJUSTMENT CONSISTS OF FREEZING      ****
C     ****    SUSPENDED LIQUID WATER (SLW) OR MELTING          ****
C     ****    SUSPENDED ICE CRYSTALS (SIW) AND DEFINING        ****
C     ****    SATURATION MIXING RATIO (QS).                    ****
C     ****      IF    TC > 0    NO SIW --> QS = QSW            ****
C     ****            TC < -40  NO SLW --> QS = QSI            ****
C     ****      -40 < TC < 0    SLW AND SIW EXIST --> QS IS    ****
C     ****                      WEIGHTED AVG. OF QSW AND QSI.  ****
C     ****      -40 < TC < 0    NO SLW AND NO SIW --> QS IS A  ****
C     ****                      WEIGHTED AVG. OF QSW AND QSI   ****
C     ****                      THAT DEPENDS ON TC ONLY.       ****
C     ****                                                     ****
C     *************************************************************
C     *************************************************************
C
      IERR = JPOINT
      JERR = KPOINT
      PRESS = P
C
      TO = T
      QVO = QVK
      QLO = QLK
      QIO = QIK
      ISTOP = 0
      NOCONV = 0
      NLOOP  = 0
C
      IF (QVK .LT. 1.E-15) QVK=0.
cxpb  IF (QLK .LT. 1.E-15) QLK=0.
cxpb  IF (QIK .LT. 1.E-15) QIK=0.
      IF (QLK .LT. 1.E-8) THEN
         QVK = QVK + QLK
         QLK=0.
      END IF
      IF (QIK .LT. 1.E-8) THEN
         QVK = QVK + QIK
         QIK=0.
      END IF
C
      IF(T .GT. TICE)  GO TO 100
C
      IF(QLK .EQ. 0.0)  GO TO 30
      IF(T .LT. TTFRZ)  GO TO 10
      IF(QIK .EQ. 0.0)  GO TO 110
      GO TO 130
C
C          T < -40 DEG. C
C     FREEZE ALL SLW UNLESS HEAT RELEASED WOULD
C     RAISE TEMP. ABOVE TTFRZ.
C
C
   10 DTMP = CPLF*QLK
      DTFR = TTFRZ-T
C
      IF(DTFR .GT. DTMP)  GO TO 20
      QIK = QIK+DTFR/CPLF
      QLK = QLK-DTFR/CPLF
      T = T+DTFR
C     WRITE(6,*) 'FREEZE LIQ TIL -40   DTMP,DTFR',DTMP,DTFR
      GO TO 130
C
   20 QIK = QIK+QLK
      QLK = 0.0
      T = T+DTMP
C     WRITE(6,*) 'FREEZE ALL LIQ   DTMP,DTFR ',DTMP,DTFR
C
   30 IF ( QIK .EQ. 0.0 ) GO TO 40
C
C     NO SLW --> QS = QSI
C
cxpb      ESI = ESATI(T)
cxpb      QSI = EPS*ESI/(P-ESI)
      QSI = qsati(T,P)
      QS  = QSI
C     WRITE(6,*) ' ICE BUT NO LIQ WATER '
      GO TO 200
C
C     NO SLW AND NO SIW -- QS DEPENDS ONLY ON T
C
cxpb   40 ESW = ESATW(T)
cxpb  ESI = ESATI(T)
cxpb  QSW = EPS*ESW/(P-ESW)
cxpb  QSI = EPS*ESI/(P-ESI)
   40 QSW = qsatw(T,P)
      QSI = qsati(T,P)
      QS  = ( QSW * AMAX1( T-TTFRZ ,0.0 )
     $      + QSI * AMIN1( TICE-T ,TDIF ) ) / TDIF
C     WRITE(6,*) ' NO LIQ WATER AND NO ICE'
      GO TO 200
C
C                  T > 0
C     MELT ALL SIW UNLESS HEAT ABSORBED WOULD
C     LOWER TEMP. BELOW TICE.
C
C
  100 IF(QIK .EQ. 0.0)  GO TO 110
      DTMP = -CPLF*QIK
      DTMLT = TICE-T
C
      IF(DTMLT .GE. DTMP)  GO TO 120
      QLK = QLK+QIK
      QIK = 0.0
      T = T+DTMP
C     WRITE(6,*) ' MELT ALL ICE   DTMP,DTMLT ',DTMP,DTMLT
C
C     NO SIW --> QS = QSW
C
cxpb  110 ESW = ESATW(T)
cxpb      QSW = EPS*ESW/(P-ESW)
  110 QSW = qsatw(T,P)
      QS  = QSW
C     WRITE(6,*) ' LIQ WATER BUT NO ICE '
      GO TO 200
C
  120 QLK = QLK-DTMLT/CPLF
      QIK = QIK+DTMLT/CPLF
      T = T+DTMLT
C     WRITE(6,*) ' MELT ICE TILL T=0   DTMP,DTMLT ',DTMP,DTMLT
C
C     MIXED WATER AND ICE CASE
C
cxpb  130 ESW = ESATW(T)
cxpb      ESI = ESATI(T)
cxpb      QSW = EPS*ESW/(P-ESW)
cxpb      QSI = EPS*ESI/(P-ESI)
  130 QSW = qsatw(T,P)
      QSI = qsati(T,P)
      QS  = (QLK*QSW+QIK*QSI)/(QLK+QIK)
C     WRITE(6,*) ' ICE AND LIQ WATER '
C
C     ****************************************************
C     ****************************************************
C     ****                                            ****
C     ****         ADJUSTMENT TO SATURATION           ****
C     ****                                            ****
C     ****     INITIALLY UNSATURATED                  ****
C     ****          CONDENSATION OF VAPOR TO SLW      ****
C     ****          SUBLIMATION OF SIW                ****
C     ****                                            ****
C     ****     INITIALLY SATURATED                    ****
C     ****          CONDENSATION OF VAPOR TO SLW      ****
C     ****          SUBLIMATION  OF VAPOR TO SIW      ****
C     ****                                            ****
C     ****     T > TICE ......... CONDENSATON ONLY    ****
C     ****     T < TTFRZ ........ SUBLIMATION ONLY    ****
C     ****     TTFRZ < T < TICE . BOTH PROCESSES      ****
C     ****                        OCCUR. PRODUCTION   ****
C     ****                        OF SLW AND SIW      ****
C     ****                        DEPENDS LINEARLY    ****
C     ****                        ON TEMPERATURE      ****
C     ****                                            ****
C     ****************************************************
C     ****************************************************
C
  200 CONTINUE
C     WRITE(6,201)  QVK,QLK,QIK,QS,T
  201 FORMAT('0',20X,'200'/10X,'QVK',13X,
     1 'QLK',13X,'QIK',13X,'QS',13X,'T',
     2 /5X,4E16.8,F12.7)
cxmk  IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
cxpb  IF(ABS(QVK/QS-1.0) .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
      IF(QVK .GT. QS)  GO TO 500
C
C     ***********************************
C     ***********************************
C     ****                           ****
C     ****     UNSATURATED CASE      ****
C     ****                           ****
C     ***********************************
C     ***********************************
C
  210 CONTINUE
C     WRITE(6,*) 'UNSATURATED CASE'
      IF(QLK .EQ. 0.0)  GO TO 300
C
C     SLW PRESENT --> EVAPORATE ALL SLW
C
      T = T-CPLC*QLK
      QVK = QVK+QLK
      QLK = 0.0
C     WRITE(6,*) ' EVAPORATE ALL LIQ WATER '
C     WRITE(6,*) 'TEMP',T
C
C     TEST FOR SATURATION
C
      IF(QIK .GT. 0.0)  GO TO 250
C
C     WATER SATURATION
C     (NO ICE PRESENT)
C
cxpb      ESW = ESATW(T)
cxpb      QSW = EPS*ESW/(P-ESW)
      QSW = qsatw(T,P)
C     WRITE(6,*) ' QSW ',QSW,' QVK ',QVK
C
C
C     ADJUST TO SATURATION (CONDENSATION)
C
cxmk  IF(QVK/QSW-1.0 .LE. 1.E-6)  GO TO 1000
cxpb  IF(QVK/QSW-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(QVK/QSW-1.0 .LE. 1.E-6)  GO TO 1000
      CALL CONDNS(T,P,QVK,QSW,ESW,QLK,ISTOP)
      GO TO 280
C
C     MIXED WATER AND ICE SATURATION
C
cxpb  250 ESW = ESATW(T)
cxpb      ESI = ESATI(T)
cxpb      QSW = EPS*ESW/(P-ESW)
cxpb      QSI = EPS*ESI/(P-ESI)
250   QSW = qsatw(T,P)
      QSI = qsati(T,P)
      QSIW = QSI
C     WRITE(6,*) ' QSIW ',QSIW,' QVK ',QVK
      IF(QVK .LT. QSIW)  GO TO 310
C
C     ADJUST TO SATURATION (CONDENSATION)
C
cxmk  IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
cxpb  IF(QVK/QSI-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
      N = 0
cxmk  CRIT = 1.E-5
cxpb  CRIT = 1.E-7                            !thresholds / 100
      CRIT = 1.E-5
      GAMFAC = EPS*CPLC*P
C     WRITE(6,*) ' 260  ADJUST TO SAT---ICE AND WATER   CRIT=',CRIT
C
  260 N = N+1
C     WRITE(6,261)  QSW,QSI,QS,QSIW,QVK,QLK,QIK
  261 FORMAT('0','   260'/5X,'QSW',13X,'QSI',13X,'QS',13X
     1,'QSIW',13X,'QVK',13X,'QLK',13X,'QIK'/7E16.8)
      IF(N .GT. 15)  GO TO 290
      SS = QVK-QSIW
cxpb      GAMMAW = GAMFAC*DTESATW(T)/(P-ESW)**2
cxpb      GAMMAI = GAMFAC*dtesati(T)/(P-ESI)**2
      GAMMAW = CPLC*dtqsatw(T,P)
      GAMMAI = CPLC*dtqsati(T,P)
      TEM = 1.+(GAMMAW*QLK+GAMMAI*QIK)/(QLK+QIK)+QIK*(QSW-QSI)/
     1 (QLK+QIK)**2
cxpb      EX1 = SS/TEM
      EX1 = MIN(SS/TEM,0.5*QVK)
      IF(N .EQ.  1)  GO TO 270
      IF(ABS(EX1/QLK) .LT. CRIT)  GO TO 280
  270 T = T+CPLC*EX1
cxpb      ESW = ESATW(T)
cxpb      ESI = ESATI(T)
cxpb      QSW = EPS*ESW/(P-ESW)
cxpb      QSI = EPS*ESI/(P-ESI)
      QSW = qsatw(T,P)
      QSI = qsati(T,P)
      QVK = QVK-EX1
      QLK = QLK+EX1
      QSIW = (QLK*QSW+QIK*QSI)/(QLK+QIK)
C     WRITE(6,111) N,T,SS,EX1,QSIW,QVK,QLK,TEM
  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO 260
C
C     CHECK THAT TEMPERATURE DOES NOT FALL BELOW TTFRZ WITH
C     LIQUID WATER PRESENT.
C
  280 IF(T .GE. TTFRZ)  GO TO 1000
C     WRITE(6,*) 'T<-40 WITH LIQ  NLOOP ',NLOOP
      IF(NLOOP .GT. 0)  GO TO 295
      NLOOP = 1
      GO TO 10
C
  290 ISTOP = 290
      GO TO 280
C
  295 ISTOP = 295
      GO TO 2000
C
C     SUBLIMATION OF ALL SIW
C
  300 IF(QIK .EQ. 0.0)  GO TO 1000
C
C     SIW PRESENT --> SUBLIMATE ALL SIW
C
  310 T = T-QIK*CPLS
      QVK = QVK+QIK
      QIK = 0.0
C     WRITE(6,*) ' NO WATER---SUBLIME ALL ICE '
C     WRITE(6,*) 'TEMP',T
C
C     TEST FOR ICE SATURATION
C         (NO SLW PRESENT)
C
cxpb      ESI = ESATI(T)
cxpb      QSI = EPS*ESI/(P-ESI)
      QSI = qsati(T,P)
C     WRITE(6,*) ' QSI ',QSI,' QVK ',QVK
C
C     ADJUST TO SATURATION (SUBLIMATION)
C
cxmk  IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
cxpb  IF(QVK/QSI-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
      CALL SUBVAP(T,P,QVK,QSI,ESI,QIK,ISTOP)
      GO TO 1000
C
C
C     *********************************
C     *********************************
C     ****                         ****
C     ****     SATURATED CASE      ****
C     ****                         ****
C     *********************************
C     *********************************
C
C
  500 IF(T .GT. TICE)  GO TO 700
      IF(T .GE. TTFRZ)  GO TO 520
C     WRITE(6,*) 'SUPERSAT  T<-40  ICE ONLY'
C
C     SUBLIMATE VAPOR UNTIL ICE SATURATION
C
cxmk  IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
cxpb  IF(QVK/QSI-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(QVK/QSI-1.0 .LE. 1.E-6)  GO TO 1000
      TQIK = 0.0
      QVKSV = QVK
      TSV   = T
      CALL SUBVAP(T,P,QVK,QSI,ESI,TQIK,ISTOP)
C     WRITE(6,*) 'TEMP AFTER VAP-->ICE ',T
C
      IF(T .GT. TTFRZ)  GO TO 510
C
      QIK = QIK+TQIK
      GO TO 1000
C
C     LATENT HEAT RELEASE PRODUCES INCONSISTENCY IN INITIAL ASSUMPTION
C     OF T <= TTFRZ.  THEREFORE SUBLIMATE ENOUGH VAPOR TO BRING
C     T = TTFRZ AND THEN PROCEED WITH CONDENSATION AND SUBLIMATION
C     FOR MIXED WATER AND ICE SATURATION CASE.
C
  510 T = TTFRZ
C     WRITE(6,*) 'RE-DO VAP-->ICE UNTIL T=-40'
      DQI = (T-TSV)/CPLS
      QIK = QIK+DQI
      QVK = QVKSV-DQI
      QSW = EPS*ESW00/(P-ESW00)
      QSI = EPS*ESI00/(P-ESI00)
      QS  = (QLK*QSW+QIK*QSI)/(QLK+QIK)
C
cxmk  IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
cxpb  IF(ABS(QVK/QS-1.0) .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
      IF(QVK .LT. QS)  GO TO 600
C
C     CONDENSE AND SUBLIMATE VAPOR UNTIL MIXED WATER AND ICE SATURATION.
C     CONDENSATION VERSUS SUBLIMATION DEPENDS ON TEMPERATURE.
C
C     METHOD IS AS FOLLOWS --
C
C       INITIAL STEP (TAG = 1) IS SUPERSATURATED (SS1 > 0)
C       SECOND GUESS (TAG = 2) IS MADE BY CONDENSING AND
C            SUBLIMATING SOME ARBITRARY AMOUNT OF VAPOR
C            (SST).  THIS GUESS MAY GIVE SS2 < SS1 (I.E.
C            CONVERGENCE).  IF SO, THIRD GUESS IS CALCU-
C            LATED ITERATIVELY AS EXPLAINED BELOW.  IF
C            SS2 > SS1 (I.E. NON-CONVERGENCE) MORE VAPOR
C            IS CONDENSED UNTIL CONVERGENCE OCCURS.
C       THIRD GUESS (TAG = 3) IS BY LINEAR INTERPOLATION IN
C            MIXING RATIO (QV3).
C       IF SS3 > 0, TAG = 1 IS REPLACED BY TAG = 3.
C       IF SS3 < 0, TAG = 2 IS REPLACED BY TAG = 3.
C
C
cxmk
c 520 CRIT = 5.E-6
cxpb
c 520 CRIT = 5.E-8     !thresholds / 100
  520 CRIT = 5.E-6
C     WRITE(6,*) 'SUPERSAT:  COND VAPOR TO ICE/LIQ SAT '
      N   = 0
cxpb      ESI = ESATI(T)
cxpb      QSI = EPS*ESI/(P-ESI)
      QSI = qsati(T,P)
      SST = AMIN1(5.*(QVK-QS),QVK-QSI)
  530 QV1 = QVK
      QL1 = QLK
      QI1 = QIK
      SS1 = QVK-QS
      T1  = T
      QSIW1 = QS
C
      DQL = SST*AMIN1(AMAX1(T-TTFRZ,0.0),TDIF)/TDIF
      DQI = SST*AMIN1(AMAX1(TICE-T,0.0),TDIF)/TDIF
      QV2 = QV1-SST
      QL2 = QL1+DQL
      QI2 = QI1+DQI
      T2  = T1+DQL*CPLC+DQI*CPLS
cxpb      ESI = ESATI(T2)
cxpb      ESW = ESATW(T2)
cxpb      QSIT = EPS*ESI/(P-ESI)
cxpb      QSWT = EPS*ESW/(P-ESW)
      QSIT = qsati(T2,P)
      QSWT = qsatw(T2,P)
      QSIW2 = (QL2*QSWT+QI2*QSIT)/(QL2+QI2)
      SS2 = QV2-QSIW2
      IF(SS2 .LT. 0.0) GO TO 540
C
C     CONDENSATION AND SUBLIMATION PRODUCE LARGER SUPERSATURATION  --
C     TRY AGIN
C
  531 T   = T2
      QVK = QV2
      QLK = QL2
      QIK = QI2
      QS  = QSIW2
C     WRITE(6,11)  N,T1,T2,SS1,SS2,QV1,QV2,QSIW1,QSIW2,QL1,QL2,
C    1QI1,QI2,DQL,DQI
   11 FORMAT('0'///'  SS2 > SS1 DIVERGENT CASE'/I6,2F14.7/
     1 (2E16.8/))
      IF(SS2/QSIW2 .LE. CRIT)  GO TO 570
      N = N+1
      IF(N .GE. 5)  GO TO 532
      GO TO 530
C
C     IF (QVK-QSI) WAS CHOSEN, THERE APPEARS TO BE NO HOPE TO FIND
C     THE SOLUTION.  HOWEVER, IF 5.*(QVK-QS) IS CHOSEN, THERE IS
C     SOME SLIM HOPE THAT THE LARGER VALUE OF (QVK-QSI) WILL PRODUCE
C     A NEGATIVE SS2.
C
  532 IF(SST .EQ. QVK-QSI)  GO TO 533
      SST = QVK-QSI
      GO TO 530
C
C     ALAS, THE SITUATION APPEARS HOPELESS
C
  533 ISTOP = 533
      GO TO 570
C
C
C     ITERATION BEGINS
C
  540 QV3 = (QV2*SS1-QV1*SS2)/(SS1-SS2)
      DQ  = QV1-QV3
      DQL = DQ*AMIN1(AMAX1(T1-TTFRZ,0.0),TDIF)/TDIF
      DQI = DQ*AMIN1(AMAX1(TICE-T1,0.0),TDIF)/TDIF
      QL3 = QL1+DQL
      QI3 = QI1+DQI
      T3  = T1+DQL*CPLC+DQI*CPLS
cxpb      ESI = ESATI(T3)
cxpb      ESW = ESATW(T3)
cxpb      QSIT = EPS*ESI/(P-ESI)
cxpb      QSWT = EPS*ESW/(P-ESW)
      QSIT = qsati(T3,P)
      QSWT = qsatw(T3,P)
      QSIW3 = (QL3*QSWT+QI3*QSIT)/(QL3+QI3)
      SS3 = QV3-QSIW3
C     WRITE(6,112) N,T1,T2,T3,SS1,SS2,SS3,QV1,QV2,QV3,QSIW1,QSIW2
C    1,QSIW3,QL1,QL2,QL3,QI1,QI2,QI3,DQL,DQI
  112 FORMAT('0',I5,3F14.7/(3E16.8/))
C
      IF(ABS(SS3)/QSIW3 .LE. CRIT)  GO TO 560
C
      N = N+1
      IF(N .EQ. 30)  GO TO 565
      IF(SS3 .GT. 0.0)  GO TO 550
C
      QV2 = QV3
      SS2 = SS3
      GO TO 540
C
  550 T1  = T3
      QV1 = QV3
      QL1 = QL3
      QI1 = QI3
      SS1 = SS3
      GO TO 540
C
  560 T = T3
      QVK = QV3
      QLK = QL3
      QIK = QI3
      QSIW = QSIW3
      GO TO 570
C
C     NON-CONVERGENT CASE.  SITUATION IS AS FOLLOWS --
C
C     SOME SLW AND SIW HAVE BEEN PRODUCED BUT PRESENCE OF SLW
C     RAISES SATURATION VAPOR PRESSURE ABOVE INITIAL (TAG = 1)
C     VALUE.  THEREFORE, EVAPORATE ALL SLW AND ITERATE TO SAT-
C     URATION (I.E. GO TO #210).
C
  565 QVK = QV3
      QLK = QL3
      QIK = QI3
      T   = T3
C     WRITE(6,22)  T,QVK,QSIW3,QLK,QIK
   22 FORMAT('0'///10X,'CONDENSATION AND SUBLIMATION TO MIXED
     1WATER AND ICE SATURATION FAILS'/5X,'T',13X,'QVK',13X,
     2'QSIW3',13X,'QLK',13X,'QIK'/F12.7,4E16.8)
      IF(SS3 .LT. 0.0)  GO TO 210
      NOCONV = 1
      ISTOP = 565
C
  570 IF(T .LE. TICE)  GO TO 1000
C
C     LATENT HEAT RELEASE PRODUCED INCONSISTENCY IN INITIAL
C     ASSUMPTION OF T <= TICE.  THEREFORE MELT ENOUGH ICE
C     TO SET T = TICE.
C
      DQLK = (T-TICE)/CPLF
C     WRITE(6,*) 'T>0  ICE AND WATER  SATD'
C     WRITE(6,*) 'DQLK,QIK',DQLK,QIK
      IF(DQLK .LE. QIK)  GO TO 580
C
C     ALL ICE REMELTED BUT T > TICE
C
C     WRITE(6,*) ' ALL ICE REMELTED BUT T>0'
      DQLK = QIK
      QIK = 0.0
      T = T-DQLK*CPLF
cxpb      ESW = ESATW(T)
cxpb      QSW  = EPS*ESW/(P-ESW)
      QSW = qsatw(T,P)
      QS = QSW
      QLK = QLK+DQLK
C     WRITE(6,571) T,DQLK,QS,QLK
  571 FORMAT('0'///'   571',F14.7,3E16.8)
C
cxmk  IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
cxpb  IF(ABS(QVK/QS-1.0) .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
      IF(QVK .GT. QS)  GO TO 700
      IF(QVK .LT. QS)  GO TO 210
C
C     SOME ICE MELTS AND T = TICE.
C
  580 T = TICE
C     WRITE(6,*) 'MELT ENOUGH ICE TO MAKE T=0'
      QLK = QLK+DQLK
      QIK = QIK-DQLK
cxpb      ESW = ESATW(T)
cxpb      QSW = EPS*ESW/(P-ESW)
      QSW = qsatw(T,P)
      QS = QSW
C     WRITE(6,581) T,QLK,QIK,QS
  581 FORMAT('0'///'   581',F14.7,3E16.8)
C
cxmk  IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
cxpb  IF(ABS(QVK/QS-1.0) .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
      IF(QVK .LT. QS)  GO TO 210
C
C     STILL SUPERSATURATED AT T = TICE WITH ICE PRESENT.
C     CONDENSE SUPERSATURATE AND USE HEAT TO MELT ICE
C     (T IS HELD CONSTANT).  IF THERE IS STILL EXCESS HEAT,
C     ADJUST TO SATURATION.
C
      DQV = QVK-QS
      DQI = DQV*HLTC/HLTF
C     WRITE(6,*) 'SUPERSAT  T=0  ICE AND LIQ '
C     WRITE(6,*) 'DQI,QIK',DQI,QIK
      IF(DQI .GT. QIK)  GO TO 590
C
C     SOME ICE MELTS
C
C     WRITE(6,*) 'MELT ENOUGH ICE TO USE HEAT  SATD'
      QVK = QVK-DQV
      QIK = QIK-DQI
      QLK = QLK+DQV+DQI
C     WRITE(6,589) T,DQV,DQI,QVK,QIK,QLK
  589 FORMAT('0'///'   589',F14.7,5E16.8)
      GO TO 1000
C
C     ALL ICE MELTS
C
  590 DQV = QIK*HLTF/HLTC
      QVK = QVK-DQV
      QLK = QLK+DQV+QIK
      QIK = 0.0
C     WRITE(6,*) 'MELT ALL ICE TO USE HEAT  SUPERSAT'
C     WRITE(6,591) T,DQV,QVK,QIK,QLK
  591 FORMAT('0'///'   591',F14.7,4E16.8)
      GO TO 700
C
  600 WRITE(6,601)  T,QVK,QS,QSW,QSI,QLK,QIK
  601 FORMAT('1',10X,'WILD AND CRAZY CASE'/F12.7,6E16.8)
      ISTOP = 600
      GO TO 2000
C
cxmk
c 700 IF(QVK/QS-1.0 .LE. 1.E-6)  GO TO 1000
cxpb
c 700 IF(QVK/QS-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
  700 IF(QVK/QS-1.0 .LE. 1.E-6)  GO TO 1000
C     WRITE(6,*)' T>=0  LIQ ONLY  SUPERSAT'
      TQLK = 0.0
      CALL CONDNS(T,P,QVK,QSW,ESW,TQLK,ISTOP)
      QLK = QLK+TQLK
C
 1000 CONTINUE
      IF(ISTOP .NE. 0)  NOCONV = 1
C     IF (NOCONV.EQ.0) GO TO 1001
C     WRITE(6,*)'NOCONV 1000','I,J',IERR,JERR
C     WRITE(6,2001)  ISTOP,TO,P,QVO,QLO,QIO
 1001 CONTINUE

Cbloss -- ensure only non-negative qik, qlk
      if (qik.lt.0.) qvk = qvk + qik
      qik = max(qik,0.)
      if (qlk.lt.0.) qvk = qvk + qlk
      qlk = max(qlk,0.)

      RETURN
C
 2000 CONTINUE
      WRITE(6,2001)  ISTOP,TO,P,QVO,QLO,QIO
 2001 FORMAT('0','   PROBLEM',I8,' IN SAT'/15X,'T=',18X,'P=',
     117X,'QVO=',17X,'QLO=',17X,'QIO='/5E20.10)
      STOP
      END
c$$$      FUNCTION ESATW(T)
c$$$C
c$$$      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
c$$$     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
c$$$C
c$$$      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
c$$$     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
c$$$csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
c$$$     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
c$$$     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
c$$$     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
c$$$     *,TT,TTD
c$$$C
c$$$      TTT = T-TT
c$$$      I = INT(TTT)
c$$$      IF(I.LT.1 .OR. I.GT.151) CALL BOMB(I,T)
c$$$      ESATW = ESWT(I)+(TTT-AINT(TTT))*DESWT(I)
c$$$      RETURN
c$$$C
c$$$      ENTRY ESATI ( T )
c$$$      TTT = T-TT
c$$$      I = INT(TTT)
c$$$      IF(I.LT.1 .OR. I.GT.151) CALL BOMB(I,T)
c$$$      ESATW = ESIT(I)+(TTT-AINT(TTT))*DESIT(I)
c$$$      RETURN
c$$$C
c$$$      ENTRY DESWDT ( T )
c$$$      TTT = T-TTD
c$$$      I = INT(TTT)
c$$$      IF(I.LT.1 .OR. I.GT.150) CALL BOMB(I,T)
c$$$      ESATW = DESWT(I)+(TTT-AINT(TTT))*(DESWT(I+1)-DESWT(I))
c$$$      RETURN
c$$$C
c$$$      ENTRY DESIDT ( T )
c$$$      TTT = T-TTD
c$$$      I = INT(TTT)
c$$$      IF(I.LT.1 .OR. I.GT.150) CALL BOMB(I,T)
c$$$      ESATW = DESIT(I)+(TTT-AINT(TTT))*(DESIT(I+1)-DESIT(I))
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE BOMB(I,T)
c$$$      COMMON/ERRDAT/NERR,IERR,JERR,PRESS,TO,QVO,QLO,QIO
c$$$      WRITE(6,*) 'BOMB: I, T',I,T
c$$$      WRITE(6,*) 'BOMB: J,K,P,T,QV,QL,QI',IERR,JERR,PRESS,TO,QVO,QLO,QIO
c$$$c$$$      CALL OUTPUT ( 0 )
c$$$c$$$      CALL TSPUT
c$$$      STOP 911
c$$$      END
      SUBROUTINE SETUPM ( DT, SFCRHO )
C
C     CONSTANTS FOR MICROPHYSICS SOURCE AND SINK TERMS
C     AND TERMINAL VELOCITIES
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,
     1      CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),
csk  2      CREVP(5),CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),
     2      CREVP(5),CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),
     3      QI0,QS0,QL0,ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,
     4      CPLC,CPLS,CPLF,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),
     5      DESWT(151),DESIT(151),TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
      COMMON/VCON/VCONR,VCONS,VCONG
C
      COMMON /STK/DTS,DTL,DZ,TDT,TDZ,DZ2,DX,TDX,X0,ALFA,
     1  NND,IND,IBD,JND,JBD,JM1,Z1,ZB,RHOSFC,RMAX
C
      DIMENSION ACT(8),ACC(3)
C
C     PHYSICAL CONSTANTS (MKS)
C
C*    DATA PIE/3.14159265/GRAV/9.81/VDIFU/2.11E-5/
C*    DATA TCOND/2.36E-2/RVAPR/4.615E2/RDRYA/2.87E2/VISK/1.259E-5/
C*    DATA HLTS/2.8336E6/HLTC/2.5E6/HLTF/3.336E5/CH2O/4.1855E3/
C*    DATA CICE/2.093E3/TICE/273.16/EPS/0.62197/TTFRZ/233.1/
C
C     LIN'S CONSTANTS(MKS) EXCEPT RMI50,RMI40 (CGS)
C
crh ********************************************************************
crh   DATA RNZR/20.E6/RNZS/3.E6/RNZG/4.E4/
crh   DATA RHOR/1.E3/RHOS/1.E2/RHOG/0.3E3/
      DATA RNZR/8.E6/  ! Lin83
      DATA RNZS/3.E6/  ! Lin83
      DATA RNZG/4.E6/  ! RH84
      DATA RHOR/1.0E3/ ! Lin83
      DATA RHOS/0.1E3/ ! Lin83
      DATA RHOG/0.4E3/ ! RH84
crh ********************************************************************
      DATA ALIN/841.99667/CLIN/4.836071224/
C*    DATA RMI50/4.80E-7/RMI40/2.46E-7/RI50/5.E-5/
C*    DATA QI0/1.E-3/QS0/0.6E-3/QL0/0.5E-3/
C*    DATA CRAUT/1.2E-1,1.064E-2/

      DATA ACC/5.0,2.0,0.5/
      DATA GAM263/1.456943/GAM275/1.608355/GAM290/1.827363/
     1     GAM325/2.54925/GAM350/3.323363/GAM380/4.694155/
     2     GAM425/8.285063/GAM450/11.631769/GAM480/17.837789/
     3     GAM625/184.860962/GAM680/496.604067/
C
      PIE=3.14159265
      GRAV=9.81
      VDIFU=2.11E-5
      TCOND=2.36E-2
      RVAPR=4.615E2
      RDRYA=2.87E2
      VISK=1.259E-5
      HLTS=2.8336E6
      HLTC=2.5E6
      HLTF=3.336E5
      CH2O=4.1855E3
      CICE=2.093E3
      TICE=273.16
      EPS=0.62197
      TTFRZ=233.1
cfu ******************************************************************
      RMI50=3.84E-6
cfu   RMI50=4.80E-7
cfu ******************************************************************
      RMI40=2.46E-7
cfu ******************************************************************
      RI50=1.E-4
cfu   RI50=5.E-5
cfu ******************************************************************
cfu ******************************************************************
cfu   QI0=1.E-3
cfu   QS0=0.6E-3
      QS0=1.E-3
      QI0=0.6E-3
cfu ******************************************************************
      QL0=0.5E-3
      CRAUT(1)=1.2E-1
      CRAUT(2)=1.064E-2
C
C     PARAMETERS PASSED FROM MODEL WHICH USES ADAMS-BASHFORTH SCHEME
C
      TDT = DT
      DTL = DT
      RHOSFC = SFCRHO
C
      CP = 3.5*RDRYA
      CPLC = HLTC/CP
      CPLS = HLTS/CP
      CPLF = HLTF/CP
      PISQ = PIE*PIE
      SCM3 = (VISK/VDIFU)**(1./3.)
C
C     ACR3:  FOUR LEAD CONSTANTS REQUIRED, THREE FOR EACH SUMMATION
C            FOUR SEPARATE PROCESSES:  RACS,SACR,GACR,GACS
C
      CRACS = PISQ*RNZR*RNZS*RHOS
      CSACR = PISQ*RNZR*RNZS*RHOR
      CGACR = PISQ*RNZR*RNZG*RHOR
      CGACS = PISQ*RNZG*RNZS*RHOS
C
C     ACT:  1-2:RACS(S-R); 3-4:SACR(R-S);
C           5-6:GACR(R-G); 7-8:GACS(S-G)
C
      ACT(1) = PIE * RNZS * RHOS
      ACT(2) = PIE * RNZR * RHOR
      ACT(6) = PIE * RNZG * RHOG
      ACT(3) = ACT(2)
      ACT(4) = ACT(1)
      ACT(5) = ACT(2)
      ACT(7) = ACT(1)
      ACT(8) = ACT(6)
C
      DO 10 I=1,3
      DO 10 K=1,4
      ACCO(I,K) = ACC(I)/(ACT(2*K-1)**((7-I)*0.25)*ACT(2*K)**(I*0.25))
   10 CONTINUE
C
C     TERMINAL VELOCITY CONSTANTS
C
      VCONR = ALIN*GAM480/(6.*ACT(2)**0.20)
      VCONS = CLIN*GAM425/(6.*ACT(1)**.0625)
crh ********************************************************************
      aa22 = 40.74 * 40.74
      CD = 4. * GRAV * RHOG / 3.0 / RHOSFC / aa22
      GCON  = SQRT(4.*GRAV*RHOG/3.0/CD)
crh   GCON  = SQRT(4.*GRAV*RHOG/1.8)
crh ********************************************************************
      VCONG = GAM450*GCON/(6.*ACT(6)**0.125)
C
C     ACR1:  SINGLE CONSTANT REQUIRED
C     FIVE SEPARATE PROCESSES:  SACW,WACS,IACR,RACI,SACI
C
      CSACW = PIE*RNZS*CLIN*GAM325/(4.*ACT(1)**0.8125)
      CWACS = PISQ*RHOS*RNZS*CLIN*GAM625/(1.0056E-10*ACT(1)**1.5625)
      CIACR = PISQ*RHOR*RNZR*ALIN*GAM680/(1.0056E-11*ACT(2)**1.7)
      CRACI = PIE*RNZR*ALIN*GAM380/(4.*ACT(2)**0.95)
      CSACI = CSACW
C
C     ACR2:  SINGLE CONSTANT REQUIRED
C     CAGCI IS FOR DRY GROWTH
C
      CGACW = PIE*RNZG*GAM350*GCON/(4.*ACT(6)**0.875)
      CGACI = CGACW*0.1
C
C     RACW
C
      CRACW = CRACI
C
C     SUBL AND REVP:  FIVE CONSTANTS FOR THREE SEPARATE PROCESSES
C
      CSSUB(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZS
      CGSUB(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZG
      CREVP(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZR
      CSSUB(2) = 0.78/SQRT(ACT(1))
      CGSUB(2) = 0.78/SQRT(ACT(6))
      CREVP(2) = 0.78/SQRT(ACT(2))
      CSSUB(3) = 0.31*SCM3*GAM263*SQRT(CLIN/VISK)/ACT(1)**0.65625
      CGSUB(3) = 0.31*SCM3*GAM275*SQRT(GCON/VISK)/ACT(6)**0.6875
      CREVP(3) = 0.31*SCM3*GAM290*SQRT(ALIN/VISK)/ACT(2)**0.725
      CSSUB(4) = TCOND*RVAPR
      CSSUB(5) = HLTS**2*VDIFU
      CGSUB(4) = CSSUB(4)
      CREVP(4) = CSSUB(4)
      CGSUB(5) = CSSUB(5)
      CREVP(5) = HLTC**2*VDIFU
C
C     GFR:  TWO CONSTANTS
C
      CGFR(1) = 20.E2*PISQ*RNZR*RHOR/ACT(2)**1.75
      CGFR(2) = 0.66
C
csk   SMLT:  FIVE CONSTANTS ( Lin et al. 1983 )
C
csk ********************************************************************
csk   CSMLT(1) = 2.*PIE*TCOND*RNZS/HLTF
csk   CSMLT(2) = CSSUB(2)
csk   CSMLT(3) = CSSUB(3)
      CSMLT(1) = 2.*PIE*TCOND*RNZS/HLTF
      CSMLT(2) = 2.*PIE*VDIFU*RNZS*HLTC/HLTF
      CSMLT(3) = CSSUB(2)
      CSMLT(4) = CSSUB(3)
      CSMLT(5) = CH2O/HLTF
csk ********************************************************************
C
C     GMLT:  FIVE CONSTANTS
C
      CGMLT(1) = 2.*PIE*TCOND*RNZG/HLTF
      CGMLT(2) = 2.*PIE*VDIFU*RNZG*HLTC/HLTF
csk ********************************************************************
      CGMLT(3) = CGSUB(2)
      CGMLT(4) = CGSUB(3)
csk   CGMLT(3) = 1.6/SQRT(ACT(6))
csk   CGMLT(4) = CGSUB(3)/SCM3
csk ********************************************************************
      CGMLT(5) = CH2O/HLTF
C
C     GWET:  TWO CONSTANTS PLUS CALC. OF SAT. VAPOR PRESSURE AT TC=0.0
C     VENTILATION COEFFICIENT OF 10 INCLUDED
C
csk ********************************************************************
csk   CGWET(1) = 20.*PIE*RNZG/SQRT(ACT(6))*HLTC*VDIFU
csk   CGWET(2) = 20.*PIE*RNZG/SQRT(ACT(6))*TCOND
      CGWET(1) = CGMLT(1)
      CGWET(2) = CGMLT(2)
      CGWET(3) = CGMLT(3)
      CGWET(4) = CGMLT(4)
csk ********************************************************************
      ES0 = 6.107799961
      CES0 = EPS*ES0
C
C     BERGRN: TWO CONSTANTS
C     C2BRG HAS CONVERSION FACTOR OF 10**3
C
      C1BRG = DTL/RMI50
clin ********************************************************************
clin  C2BRG = RI50**2*1.E3 ! Error
      C2BRG = PIE*RI50**2*1.E3
clin ********************************************************************
C
C     SATURATION VAPOR PRESSURE VALUES ARE TABULATED EACH WHOLE
C     DEG. (C) FOR -100 < TC < 50 FOR WATER AND -100 < TC < 0 FOR ICE.
C     DERIVATIVES ARE CENTERED AT EACH HALF DEG. (C).
C
      TT  = TICE-101.
      TTD = TICE-100.5
      ESWT(1) = GGESW(TICE-100.)
      ESIT(1) = GGESI(TICE-100.)
      ITC = -99
      DO 20 I=2,151
      ESWT(I)    = GGESW(TICE+FLOAT(ITC))
      DESWT(I-1) = ESWT(I)-ESWT(I-1)
      ITC = ITC+1
   20 CONTINUE
      DESWT(151) = 6.25
C
      ITC = -99
      DO 30 I=2,101
      ESIT(I)    = GGESI(TICE+FLOAT(ITC))
      DESIT(I-1) = ESIT(I)-ESIT(I-1)
      ITC = ITC+1
   30 CONTINUE
      DESIT(101) = DESWT(101)
      DO 40 I=102,151
      ESIT(I) = ESWT(I)
      DESIT(I) = DESWT(I)
   40 CONTINUE
C
C     SATURATION ADJUSTMENT:  RANGE OF MIXED ICE-WATER SATURATION
C     AND SATURATION VAPOR PRESSURE OVER WATER AND ICE AT T = -40 C.
C
      TDIF  = TICE-TTFRZ
      ESW00 = GGESW(TICE-TDIF)
      ESI00 = GGESI(TICE-TDIF)
C
      RETURN
      END
      FUNCTION GGESI(T)
C
C     SATURATION VAPOR PRESSURE OVER ICE
C      (GOFF AND GRATCH)
C
C     ESI     SATURATION VAPOR PRESSURE  (MB)
C     T       TEMP  (KELVIN)
C
      DATA C1/-9.09718/C2/-3.56654/C3/0.876793/C4/0.78583503/
C
      A = 273.16/T
      B = C1*(A-1.0)+C2*ALOG10(A)+C3*(1.0-1.0/A)+C4
      GGESI = 10.0**B
      RETURN
      END
      FUNCTION GGESW (TA)
C
C     SATURATION VAPOR PRESSURE OVER WATER
C          (GOFF AND GRATCH)
C
C     TA IS TEMPERATURE IN DEG KELVIN
C     ES IS SATURATION VAPOR PRESSURE IN MB
C
      DATA TS/373.16/F5/3.0057149/
      E1 = 11.344*(1.0 - TA/TS)
      E2 = -3.49149*(TS/TA - 1.0)
      F1 = -7.90298*(TS/TA - 1.0)
      F2 = 5.02808*ALOG10(TS/TA)
      F3 = -1.3816E-7*(10.0**E1 -1.0)
      F4 = 8.1328E-3*(10.0**E2 - 1.0)
      F = F1 + F2 + F3 + F4 + F5
      GGESW = 10.0**F
      RETURN
      END
      SUBROUTINE CONDNS(T,P,QVK,QSW,ESW,QLK,ISTOP)
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
cxmk  DATA CRIT/1.E-6/
cxpb  DATA CRIT/1.E-8/    !thresholds / 100
      DATA CRIT/1.E-6/
C
      N = 0
      GAMFAC = EPS*CPLC*P
C     WRITE(6,*) '  10  ADJUST TO SAT---NO ICE   CRIT=',CRIT
   10 N = N+1
      IF(N .GT. 10)  GO TO 1000
      SS = QVK-QSW
cxpb      GAMMAW = GAMFAC*DTESATW(T)/(P-ESW)**2
      GAMMAW = CPLC*dtqsatw(T,P)
cxpb      EX = SS/(1.+GAMMAW)
      EX = MIN(SS/(1.+GAMMAW),0.5*QVK)
      IF(N .EQ.  1)  GO TO  20
      IF(ABS(EX/QLK) .LT. CRIT)  RETURN
   20 T = T+CPLC*EX
cxpb      ESW = ESATW(T)
cxpb      QSW = EPS*ESW/(P-ESW)
      QSW = qsatw(T,P)
      QVK = QVK-EX
      QLK = QLK+EX
C     WRITE(6,111) N,T,SS,EX,QSW,QVK,QLK,GAMMAW
  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO  10
C
 1000 ISTOP = 1111
      RETURN
      END
      SUBROUTINE SUBVAP(T,P,QVK,QSI,ESI,QIK,ISTOP)
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
cxmk  DATA CRIT/0.5E-7/
cxpb  DATA CRIT/0.5E-9/    !thresholds / 100
      DATA CRIT/0.5E-7/
C
      N = 0
      GAMFAC = EPS*CPLS*P
C     WRITE(6,*) '      ADJUST TO SAT---ICE ONLY   CRIT=',CRIT
   10 N = N+1
      IF(N .GT. 10)  GO TO 1000
      SS = QVK-QSI
cxpb      GAMMAI = GAMFAC*dtesati(T)/(P-ESI)**2
      GAMMAI = CPLS*dtqsati(T,P)
cxpb      EX = SS/(1.+GAMMAI)
      EX = MIN(SS/(1.+GAMMAI),0.5*QVK)
      IF(N .EQ.  1)  GO TO 20
      IF(ABS(EX/QIK) .LT. CRIT)  RETURN
   20 T = T+CPLS*EX
cxpb      ESI = ESATI(T)
cxpb      QSI = EPS*ESI/(P-ESI)
      QSI = qsati(T,P)
      QVK = QVK-EX
      QIK = QIK+EX
C     WRITE(6,111) N,T,SS,EX,QSI,QVK,QIK,GAMMAI
  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO 10
C
 1000 ISTOP = 2222
      RETURN
      END
cfu *******************************************************************
      SUBROUTINE BERGRN(TC,QL,QI,QLRHO,PSFW,PSFI)
cfu   SUBROUTINE BERGRN(TC,QI,QLRHO,PSFW,PSFI)
cfu ******************************************************************
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
cfu ******************************************************************
      DIMENSION A1T(0:31),A2T(0:31)
cfu   DIMENSION A1T(31),A2T(31)
cfu ******************************************************************
      DATA A1T/0.,0.7939E-7,0.7841E-6,0.3369E-5,0.4336E-5,0.5285E-5,
     1         0.3728E-5,0.1852E-5,0.2991E-6,0.4248E-6,0.7434E-6,
     2         0.1812E-5,0.4394E-5,0.9145E-5,0.1725E-4,0.3348E-4,
     3         0.1725E-4,0.9175E-5,0.4412E-5,0.2252E-5,0.9115E-6,
     4         0.4876E-6,0.3473E-6,0.4758E-6,0.6306E-6,0.8573E-6,
     5         0.7868E-6,0.7192E-6,0.6513E-6,0.5956E-6,0.5333E-6,
     6         0.4834E-6/
      DATA A2T/0.,0.4006,0.4831,0.5320,0.5307,0.5319,
     1         0.5249,0.4888,0.3894,0.4047,0.4318,
     2         0.4771,0.5183,0.5463,0.5651,0.5813,
     3         0.5655,0.5478,0.5203,0.4906,0.4447,
     4         0.4126,0.3960,0.4149,0.4320,0.4506,
     5         0.4483,0.4460,0.4433,0.4413,0.4382,
     6         0.4361/
C
      TC1 = AMAX1(TC,-30.0)
cfu ****************************************************************
      if ( tc1 .gt. -1.0 ) then
      a1 = a1t(1)
      a2 = a2t(1)
      else
      A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+
     *      A1T(-INT(TC1)+1)
      A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+
     *      A2T(-INT(TC1)+1)
      endif
cfu   A1  = A1T(-INT(TC1)+1)
cfu   A2  = A2T(-INT(TC1)+1)
cfu ****************************************************************
      A21 = 1.0 - A2
C
      DT1 = (RMI50**A21 - RMI40**A21) / (A1*A21)
C
C     NOTE:  MKS UNITS, UI50=1.0 M/SEC, EIW=1.0
C
      PSFW = C1BRG*QI/DT1*(A1*RMI50**A2+QLRHO*C2BRG)
cfu *******************************************************************
      if ( QL .gt. 5.0e-6 ) then
        PSFI = QI/DT1
      else
ctry    psfi = 0.0
        PSFI = QI/DT1
      endif
cfu ******************************************************************
      RETURN
      END
cfu *******************************************************************
      REAL FUNCTION IDW(TC,qik)
cfu   REAL FUNCTION IDW(TC,)
cfu *******************************************************************
      COMMON/RHOCON/RHOFAC,RHO
cfu *******************************************************************
      DIMENSION A1T(0:31),A2T(0:31)
cfu   DIMENSION A1T(31),A2T(31)
cfu *******************************************************************
C
C     NOTE -- A2T IS IDENTICAL TO A2T IN SUBROUTINE BERGRN, BUT A1T
C     IS MULTIPLIED BY A FACTOR OF 1.E-5.
C
cfu *******************************************************************
      DATA A1T/.0,.7939E-12,.7841E-11,.3369E-10,0.4336E-10,0.5285E-10,
cfu   DATA A1T/0.7939E-12,0.7841E-11,0.3369E-10,0.4336E-10,0.5285E-10,
cfu *******************************************************************
     1         0.3728E-10,0.1852E-10,0.2991E-11,0.4248E-11,0.7434E-11,
cfu *******************************************************************
     2         0.1812E-10,0.4394E-10,0.9145E-10,0.1725E-9 ,0.3348E-9 ,
cfu  2         0.1812E-10,0.4394E-10,0.9145E-10,0.1725E-11,0.3348E-9 ,
cfu *******************************************************************
     3         0.1725E-9 ,0.9175E-10,0.4412E-10,0.2252E-10,0.9115E-11,
     4         0.4876E-11,0.3473E-11,0.4758E-11,0.6306E-11,0.8573E-11,
     5         0.7868E-11,0.7192E-11,0.6513E-11,0.5956E-11,0.5333E-11,
     6         0.4834E-11/
cfu *******************************************************************
      DATA A2T/0.0,0.4006,0.4831,0.5320,0.5307,0.5319,
cfu   DATA A2T/0.4006,0.4831,0.5320,0.5307,0.5319,
cfu *******************************************************************
     1         0.5249,0.4888,0.3894,0.4047,0.4318,
     2         0.4771,0.5183,0.5463,0.5651,0.5813,
     3         0.5655,0.5478,0.5203,0.4906,0.4447,
     4         0.4126,0.3960,0.4149,0.4320,0.4506,
     5         0.4483,0.4460,0.4433,0.4413,0.4382,
     6         0.4361/
      DATA RMINUC/1.05E-15/
C
      TC1 = AMAX1(TC,-30.0)
cfu ********************************************************************
      A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+
     *      A1T(-INT(TC1)+1)
      A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+
     *      A2T(-INT(TC1)+1)
cfu   A1  = A1T(-INT(TC1)+1)
cfu   A2  = A2T(-INT(TC1)+1)
cfu ********************************************************************
C
cfu ********************************************************************
      fnc = 0.01 * exp ( - 0.6 * tc )
      fmi = rho * qik / fnc * 1000.0
      IDW = EXP(-0.6*TC)*A1*fmi**A2/RHO
cfu   IDW = EXP(-0.6*TC)*A1*RMINUC**A2/RHO
cfu ********************************************************************
      RETURN
      END
      FUNCTION ACR1(Q,QRHO,C,PWR)
      COMMON/RHOCON/RHOFAC,RHO
C
C     RHO FACTOR FOR SOURCE TERMS:  PSACW,PWACS,PRACI,PSACI
C
      FAC = RHOFAC
      GO TO 10
C
C     RHO FACTOR FOR SOURCE TERMS:  PGACW,PGACI
C
      ENTRY ACR2(Q,QRHO,C,PWR)
      FAC = 1./SQRT(RHO)
C
C     CALC ACCRETION
C
   10 CONTINUE
      ACR1 = C * Q * QRHO**PWR * FAC
      RETURN
      END
      FUNCTION ACR3(V1,V2,QRHO1,QRHO2,C,CAC)
C
C      SOURCE TERMS:  PRACS,PSACR,PGACR,PGACS
C
      DIMENSION CAC(3)
      COMMON/RHOCON/RHOFAC,RHO
      A=0.0
      DO 10 K=1,3
      A = A + CAC(K) * (QRHO1**((7-K)*0.25) * QRHO2**(K*0.25))
   10 CONTINUE
      ACR3 = C * ABS(V1-V2) * A / RHO
      RETURN
      END
      FUNCTION REVP(S,TSQ,QS,QRHO,C)
      DIMENSION C(5)
      COMMON/RHOCON/RHOFAC,RHO
C
C     EVAPORATION
C
      REVP = C(1)*TSQ*QS*(1.0-S) * (C(2)*SQRT(QRHO)+C(3)*QRHO**0.725)
     1   / (C(4)*TSQ+C(5)*QS*RHO)
      RETURN
      END
      FUNCTION AUT(C,Q,Q0)
C
C      SOURCE TERMS:  PGAUT,PSAUT
C
      AUT = AMAX1(C*(Q-Q0),0.0)
      RETURN
      END
      FUNCTION RAUT(C,QEXRHO)
      DIMENSION C(2)
      COMMON/RHOCON/RHOFAC,RHO
C
C     SEE ORVILLE AND KOPP (1977)
C
      RAUT = QEXRHO**3 / ((C(1)*QEXRHO+C(2))*RHO)
      RETURN
      END
      FUNCTION GFR(TC,QRR,C)
      DIMENSION C(2)
      COMMON/RHOCON/RHOFAC,RHO
      GFR = C(1)*(EXP(-C(2)*TC)-1.0)*QRR**1.75/RHO
      RETURN
      END
csk ********************************************************************
csk   FUNCTION SMLT(TC,QRHO,C)
      FUNCTION SMLT(TC,DQS,QSRHO,PSACW,PSACR,C)
csk   DIMENSION C(3)
      DIMENSION C(5)
      COMMON/RHOCON/RHOFAC,RHO
C
csk   SMLT = C(1)*TC*(C(2)*SQRT(QRHO)+C(3)*SQRT(RHOFAC)*QRHO**0.65625)
csk  1/RHO
      SMLT = (C(1)*TC/RHO-C(2)*DQS) * (C(3)*SQRT(QSRHO)+
     1   C(4)*QSRHO**0.65625*SQRT(RHOFAC)) + C(5)*TC*(PSACW+PSACR)
csk ********************************************************************
      RETURN
      END
      FUNCTION GMLT(TC,DQS,QGRHO,PGACW,PGACR,C)
      DIMENSION C(5)
      COMMON/RHOCON/RHOFAC,RHO
C
C     NOTE:  PGACW AND PGACR MUST BE CALC BEFORE GMLT IS CALLED
C
      GMLT = (C(1)*TC/RHO-C(2)*DQS) * (C(3)*SQRT(QGRHO)+
     1   C(4)*QGRHO**0.6875/RHO**0.25) + C(5)*TC*(PGACW+PGACR)
      RETURN
      END
clin *******************************************************************
      FUNCTION GWET(DQS,TC,QGRHO,PGACI,PGACS,C,SIMAX,SSMAX)
clin  FUNCTION GWET(DQS,TC,QGRHO,PGACI,PGACS,C)
clin *******************************************************************
      DIMENSION C(4) !csk
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
      COMMON/RHOCON/RHOFAC,RHO
C
C     NOTE:  PGACI AND PGACS MUST BE CALC BEFORE GWET IS CALLED
C      FACTOR OF 10 CONVERTS PGACI TO PGACI PRIME
clin   FACTOR OF 10 CONVERTS PGACS TO PGACS PRIME
C
clin *******************************************************************
      x = AMIN1(PGACI * 10.0,SIMAX)
      y = AMIN1(10.*PGACS,SSMAX)
      GWET = AMAX1((C(2)*DQS-C(1)*TC/RHO)*HLTF/(HLTF+CH2O*TC)* !csk
     *   (C(3)*SQRT(QGRHO)+C(4)*QGRHO**0.6875/RHO**0.25) +  ! csk
     1   (x+y)*(1.0-CICE*TC/(HLTF+CH2O*TC)),0.0)
clin 1   (10.*PGACI+PGACS)*(1.0-CICE*TC/(HLTF+CH2O*TC)),0.0)
clin *******************************************************************
      RETURN
      END
      FUNCTION GSUB(S,TSQ,QS,QRHO,C)
      DIMENSION C(5)
      COMMON/RHOCON/RHOFAC,RHO
C
C     GRAUPEL SUBLIMATION
C
      GSUB = C(1)*TSQ*QS*(S-1.0) * (C(2)*SQRT(QRHO) +
     1       C(3)*QRHO**0.6875/RHO**0.25) / (C(4)*TSQ+C(5)*QS*RHO)
      RETURN
C
C     SNOW SUBLIMATION
C
      ENTRY SSUB(S,TSQ,QS,QRHO,C)
      GSUB = C(1)*TSQ*QS*(S-1.0) * (C(2)*SQRT(QRHO) +
     1       C(3)*QRHO**0.65625*SQRT(RHOFAC)) / (C(4)*TSQ+C(5)*QS*RHO)
      RETURN
      END
      FUNCTION VTRR(QRHO)
C
C     CALC MASS-WEIGHTED MEAN TERMINAL VELOCITIES
C
      COMMON/VCON/VCONR,VCONS,VCONG
      COMMON/RHOCON/RHOFAC,RHO
C
C     RAIN
C
      VTRR = AMIN1(VCONR * QRHO**0.2 * RHOFAC,10.0)
      RETURN
C
C    C SNOW
C
      ENTRY VTRS ( QRHO )
      VTRR = VCONS * QRHO**0.0625 * RHOFAC
      RETURN
C
C     GRAUPEL
C
      ENTRY VTRG ( QRHO )
      VTRR = AMIN1(VCONG * QRHO**0.125 / SQRT(RHO),20.0)
      RETURN
      END
      SUBROUTINE RMICRO
C
      COMMON/VTERM/VTS,VTG,VTR
      COMMON/MVAR/T,PRESR,QMIX(6),TSS(7),PSS(26),RHOK
      COMMON /STK/DT,DTL,DZ,TDT,TDZ,DZ2,DX,TDX,X0,ALFA,
     1  NND,IND,IBD,JND,JBD,JM1,Z1,ZB,RHOS,RMAX
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
C
      COMMON/VCON/VCONR,VCONS,VCONG
C
      COMMON/RHOCON/RHOFAC,RHO
C
      DIMENSION C3RACS(3)
     *,C3SACR(3),C3GACR(3),C3GACS(3)
C
      DIMENSION QRHO(6),QLOGIC(6),THRSH(6),SMAX(6)
C
C
      LOGICAL QLOGIC,VAPR,LIQ,ICE,SNO,GRAUP,RAIN
C
      EQUIVALENCE (QMIX(1),QVK),(QMIX(2),QLK),(QMIX(3),QIK),
     *  (QMIX(4),QSK),(QMIX(5),QGK),(QMIX(6),QRK),(QRHO(1),QVR),
     *  (QRHO(2),QLR),(QRHO(3),QIR),(QRHO(4),QSR),(QRHO(5),QGR),
     *  (QRHO(6),QRR),(QLOGIC(1),VAPR),(QLOGIC(2),LIQ),
     *  (QLOGIC(3),ICE),(QLOGIC(4),SNO),(QLOGIC(5),GRAUP),
     *  (QLOGIC(6),RAIN),(TSSV,TSS(1)),(TSSL,TSS(2)),(TSSI,TSS(3)),
     *  (TSSS,TSS(4)),(TSSG,TSS(5)),(TSSR,TSS(6)),(TTMP,TSS(7)),
     *  (SVMAX,SMAX(1)),(SLMAX,SMAX(2)),(SIMAX,SMAX(3)),
     *  (SSMAX,SMAX(4)),(SGMAX,SMAX(5)),(SRMAX,SMAX(6))
C
      EQUIVALENCE (PRAUT,PSS(1)),(PRACW,PSS(2)),(PRACS,PSS(3))
     *,(PSACW,PSS(4)),(PWACS,PSS(5)),(PGACW,PSS(6)),(PGMLT,PSS(7))
     *,(PSMLT,PSS(8)),(PREVP,PSS(9)),(PIACR,PSS(10)),(PSACR,PSS(11))
     *,(PGACR,PSS(12)),(PGFR,PSS(13)),(PGACS,PSS(14)),(PGAUT,PSS(15))
     *,(PGACI,PSS(16)),(PGWORD,PSS(17)),(PRACI,PSS(18)),(PSAUT,PSS(19))
     *,(PSACI,PSS(20)),(PSSUB,PSS(21)),(PGSUB,PSS(22)),(PSFW,PSS(23))
     *,(PSFI,PSS(24)),(PIDW,PSS(25)),(PGWET,PSS(26))
C
      EQUIVALENCE (C3RACS(1),ACCO(1,1)),(C3SACR(1),ACCO(1,2))
     *,(C3GACR(1),ACCO(1,3)),(C3GACS(1),ACCO(1,4))
C
cxmk  DATA THRSH/0.0,5*1.E-6/
      DATA THRSH/0.0,5*1.E-8/     !threshold / 100
C
C     ***************************
C     ***************************
C     ****                   ****
C     ****   PRELIMINARIES   ****
C     ****                   ****
C     ***************************
C     ***************************
C
C
C     LOCAL VARIABLES -- TEMP, DENSITY, TEMP. DEPENDENT EFFICIENCIES
C     AND SATURATION VARIABLES
C
      RHO = RHOK
      RHOFAC = SQRT(RHOS/RHO)
      QL0RHO = QL0 * RHO
      TSQ = T**2
C
      ESW   = ESATW(T)
      QSW   = EPS*ESW/(PRESR-ESW)
      DQS   = QVK-QSW
C
C     ZERO SOURCE AND SINK TERMS
C
      DO 1 K=1,26
    1 PSS(K) = 0.0
C
C     DEFINE MIXING RATIOS GREATER THAN THRESHOLD FOR
C     ACCRETION PROCESSES AND MASS DENSITIES
C
C           1:  WATER VAPOR
C           2:  CLOUD (SUSPENDED LIQUID) WATER
C           3:  CLOUD ICE CRYSTALS
C           4:  SNOW
C           5:  GRAUPEL
C           6:  RAIN
C
      DO 10 K=1,6
      TSS(K) = 0.0
      QLOGIC(K) = QMIX(K) .GT. THRSH(K)
      SMAX(K) = QMIX(K)/TDT
   10 QRHO(K) = QMIX(K) * RHO
      TSS(7)  = 0.0
C
C     TERMINAL VELOCITIES
C
      IF(RAIN)   VTR = VTRR(QRR)
C
      IF(.NOT. LIQ)  GO TO 150
      IF(QLK .LE. QL0)  GO TO 130
C
C     PRAUT
C
      QEXRHO = QLR - QL0RHO
      PRAUT = AMIN1(RAUT(CRAUT,QEXRHO),SLMAX)
C
  130 CONTINUE
      IF(.NOT. RAIN)  GO TO 150
C
C     PRACW
C
      PRACW = AMIN1(ACR1(QLK,QRR,CRACW,0.95),SLMAX)
C
  150 CONTINUE
      IF(QRK .EQ. 0.0 .OR. DQS .GE. 0.0)  GO TO 330
C
C     PREVP
C
      SW = QVK/QSW
      PREVP = AMIN1(REVP(SW,TSQ,QSW,QRR,CREVP),SRMAX)
C
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****                              ****
C     **************************************
C
  330 CONTINUE
C
      TSSV = PREVP
      TSSL = -(PRAUT+PRACW)
      TSSR = PRAUT+PRACW-PREVP
      TTMP = -HLTC*PREVP/CP
C
C     WRITE(6,111) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  111 FORMAT('0',30X,'T>=0  '//10X,'TSSV=',E16.8/
     *,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8,
     */10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
      TEM = 0.0
      DO 1112 K=1,6
 1112 TEM = TEM+TSS(K)
C     WRITE(6,1113)  TEM
 1113 FORMAT('0',10X,'TOTAL SOURCE AND SINK=',E16.8)
C
C
C     PRINT VALUES FOR EACH PROCESS
C
C     WRITE(6,2000)'ALL TEMP','PGACW',PGACW,'PSACW',PSACW,'PRAUT',PRAUT,
C    1   'PRACW',PRACW,'PRACS',PRACS,'PGACS',PGACS,'PREVP',PREVP
C     WRITE(6,2000)'TC > 0','PSMLT',PSMLT,'PWACS',PWACS,'PGMLT',PGMLT
C     WRITE(6,2000) 'TC < 0','PSFW',PSFW,'PSFI',PSFI,'PRACI',PRACI,
C    1   'PIACR',PIACR,'PGACI',PGACI,'PSACI',PSACI,'PSAUT',PSAUT,
C    2   'PGFR',PGFR,'PSACR',PSACR,'PGACR',PGACR,'PSSUB',PSSUB,
C    3   'PGAUT',PGAUT,'PGSUB',PGSUB,'PIDW',PIDW
C     WRITE(6,2000)'WET-DRY','PGWET',PGWET,'PGDRY',PGDRY,'PGWORD',PGWORD
      RETURN
      ENTRY ERRMIC
      WRITE(6,*) '********* ERROR REPORT FROM RMICRO *********'
      WRITE(6,666) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  666 FORMAT('1',30X,'T<0  '//10X,'TSSV=',E16.8/
     *,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8,
     */10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
C
C     PRINT VALUES FOR EACH PROCESS
C
      WRITE(6,2000)'ALL TEMP','PGACW',PGACW,'PSACW',PSACW,'PRAUT',PRAUT,
     1   'PRACW',PRACW,'PRACS',PRACS,'PGACS',PGACS,'PREVP',PREVP
      WRITE(6,2000)'TC > 0','PSMLT',PSMLT,'PWACS',PWACS,'PGMLT',PGMLT
      WRITE(6,2000) 'TC < 0','PSFW',PSFW,'PSFI',PSFI,'PRACI',PRACI,
     1   'PIACR',PIACR,'PGACI',PGACI,'PSACI',PSACI,'PSAUT',PSAUT,
     2   'PGFR',PGFR,'PSACR',PSACR,'PGACR',PGACR,'PSSUB',PSSUB,
     3   'PGAUT',PGAUT,'PGSUB',PGSUB,'PIDW',PIDW
 2000 FORMAT(1H0,'PROCESSES:  ',A10/
     1   3(1H0,5(A6,'=',E12.5,4X)/))
      RETURN
      END
      SUBROUTINE WATSAT(JPOINT,KPOINT,T,P,QVK,QLK)
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
      COMMON/PATHO/NOCONV,NCTOT
C **DEBUG
      COMMON/ERRDAT/NERR,IERR,JERR,PRESS,TO,QVO,QLO,QIO
C
      IERR = JPOINT
      JERR = KPOINT
      PRESS = P
C
      TO = T
      QVO = QVK
      QLO = QLK
      QIO = 0.0
      ISTOP = 0
      NOCONV = 0
C
      IF (QVK .LT. 1.E-15) QVK=0.
      IF (QLK .LT. 1.E-15) QLK=0.
C
      ESW = ESATW(T)
      QSW = EPS*ESW/(P-ESW)
      QS  = QSW
C     IF(IERR.EQ.14 .AND. JERR.EQ.4)WRITE(6,201)  QVK,QLK,QS,T
  201 FORMAT('0',20X,'200'/10X,'QVK',13X,
     1 'QLK',13X,'QS',13X,'T',
     2 /5X,3E16.8,F12.7)
C
C SEPARATE SATURATED, UNSATURATED, SUPERSATURATED CASES
C
cxmk  IF(ABS(QVK/QS-1.0) .LE. 1.E-6)  GO TO 1000
      IF(ABS(QVK/QS-1.0) .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      IF(QVK .GT. QS)  GO TO 700
      IF (QLK .EQ. 0.0) GO TO 1000
C
C     UNSATURATED CASE
C     SLW PRESENT --> EVAPORATE ALL SLW
C
      T = T-CPLC*QLK
      QVK = QVK+QLK
      QLK = 0.0
C     WRITE(6,*) ' EVAPORATE ALL LIQ WATER '
C     WRITE(6,*) 'TEMP',T
C
      ESW = ESATW(T)
      QSW = EPS*ESW/(P-ESW)
C     WRITE(6,*) ' QSW ',QSW,' QVK ',QVK
C
C
C     ADJUST TO SATURATION (CONDENSATION)
C
cxmk  IF(QVK/QSW-1.0 .LE. 1.E-6)  GO TO 1000
      IF(QVK/QSW-1.0 .LE. 1.E-8)  GO TO 1000  !thresholds / 100
      CALL CONDNS(T,P,QVK,QSW,ESW,QLK,ISTOP)
C     WRITE(6,*) 'QV',QVK,'QSW',QSW,'QL',QLK
      GO TO 1000
  700 CONTINUE
C
C SUPERSATURATED CASE --- ADJUST TO SATURATION
C
C     WRITE(6,*)' T>=0  LIQ ONLY  SUPERSAT'
      TQLK = 0.0
      CALL CONDNS(T,P,QVK,QSW,ESW,TQLK,ISTOP)
      QLK = QLK+TQLK
C     IF(IERR.EQ.14 .AND. JERR.EQ.4)
C    1   WRITE(6,*) 'QV',QVK,'QSW',QSW,'QL',QLK
C
 1000 CONTINUE
      IF(ISTOP .NE. 0)  NOCONV = 1
      RETURN
      END
      SUBROUTINE FMICRO
C
      COMMON/VTERM/VTS,VTG,VTR
      COMMON/MVAR/T,PRESR,QMIX(6),TSS(7),PSS(26),RHOK
      COMMON /STK/DT,DTL,DZ,TDT,TDZ,DZ2,DX,TDX,X0,ALFA,
     1  NND,IND,IBD,JND,JBD,JM1,Z1,ZB,RHOS,RMAX
C
      COMMON/SSCON/CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR
     *,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)
csk  *,CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2),QI0,QS0,QL0,ES0
     *,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0
     *,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF
     *,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),DESWT(151),DESIT(151)
     *,TT,TTD
C
      COMMON/FIZCON/PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC
     *,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
C
C
      COMMON/VCON/VCONR,VCONS,VCONG
C
      COMMON/RHOCON/RHOFAC,RHO
C
      DIMENSION C3RACS(3)
     *,C3SACR(3),C3GACR(3),C3GACS(3)
C
      DIMENSION QRHO(6),QLOGIC(6),THRSH(6),SMAX(6)
C
C
      LOGICAL QLOGIC,VAPR,LIQ,ICE,SNO,GRAUP,RAIN
C
      EQUIVALENCE (QMIX(1),QVK),(QMIX(2),QLK),(QMIX(3),QIK),
     *  (QMIX(4),QSK),(QMIX(5),QGK),(QMIX(6),QRK),(QRHO(1),QVR),
     *  (QRHO(2),QLR),(QRHO(3),QIR),(QRHO(4),QSR),(QRHO(5),QGR),
     *  (QRHO(6),QRR),(QLOGIC(1),VAPR),(QLOGIC(2),LIQ),
     *  (QLOGIC(3),ICE),(QLOGIC(4),SNO),(QLOGIC(5),GRAUP),
     *  (QLOGIC(6),RAIN),(TSSV,TSS(1)),(TSSL,TSS(2)),(TSSI,TSS(3)),
     *  (TSSS,TSS(4)),(TSSG,TSS(5)),(TSSR,TSS(6)),(TTMP,TSS(7)),
     *  (SVMAX,SMAX(1)),(SLMAX,SMAX(2)),(SIMAX,SMAX(3)),
     *  (SSMAX,SMAX(4)),(SGMAX,SMAX(5)),(SRMAX,SMAX(6))
C
      EQUIVALENCE (PRAUT,PSS(1)),(PRACW,PSS(2)),(PRACS,PSS(3))
     *,(PSACW,PSS(4)),(PWACS,PSS(5)),(PGACW,PSS(6)),(PGMLT,PSS(7))
     *,(PSMLT,PSS(8)),(PREVP,PSS(9)),(PIACR,PSS(10)),(PSACR,PSS(11))
     *,(PGACR,PSS(12)),(PGFR,PSS(13)),(PGACS,PSS(14)),(PGAUT,PSS(15))
     *,(PGACI,PSS(16)),(PGWORD,PSS(17)),(PRACI,PSS(18)),(PSAUT,PSS(19))
     *,(PSACI,PSS(20)),(PSSUB,PSS(21)),(PGSUB,PSS(22)),(PSFW,PSS(23))
     *,(PSFI,PSS(24)),(PIDW,PSS(25)),(PGWET,PSS(26))
C
      EQUIVALENCE (C3RACS(1),ACCO(1,1)),(C3SACR(1),ACCO(1,2))
     *,(C3GACR(1),ACCO(1,3)),(C3GACS(1),ACCO(1,4))
C
cxmk  DATA THRSH/0.0,5*1.E-6/
      DATA THRSH/0.0,5*1.E-8/    !threshold / 100
C
C     ***************************
C     ***************************
C     ****                   ****
C     ****   PRELIMINARIES   ****
C     ****                   ****
C     ***************************
C     ***************************
C
C
C     LOCAL VARIABLES -- TEMP, DENSITY, TEMP. DEPENDENT EFFICIENCIES
C     AND SATURATION VARIABLES
C
      RHO = RHOK
      RHOFAC = SQRT(RHOS/RHO)
      QL0RHO = QL0 * RHO
      TC = T - TICE
      TSQ = T**2
C
      ESW   = ESATW(T)
      QSW   = EPS*ESW/(PRESR-ESW)
      DQS   = QVK-QSW
      DQS0 = CES0/(PRESR-ES0) - QVK
C
C     ZERO SOURCE AND SINK TERMS
C
      DO 1 K=1,26
    1 PSS(K) = 0.0
C
C     DEFINE MIXING RATIOS GREATER THAN THRESHOLD FOR
C     ACCRETION PROCESSES AND MASS DENSITIES
C
C           1:  WATER VAPOR
C           2:  CLOUD (SUSPENDED LIQUID) WATER
C           3:  CLOUD ICE CRYSTALS
C           4:  SNOW
C           5:  GRAUPEL
C           6:  RAIN
C
      DO 10 K=1,6
      TSS(K) = 0.0
      QLOGIC(K) = QMIX(K) .GT. THRSH(K)
      SMAX(K) = QMIX(K)/TDT
   10 QRHO(K) = QMIX(K) * RHO
      TSS(7)  = 0.0
C
C     TERMINAL VELOCITIES
C
      IF(RAIN)   VTR = VTRR(QRR)
      IF(GRAUP)  VTG = VTRG(QGR)
C
      IF(.NOT. LIQ)  GO TO 150
      IF(.NOT. GRAUP)  GO TO 110
C
C     PGACW
C
      PGACW = AMIN1(ACR2(QLK,QGR,CGACW,0.875),SLMAX)
C
  110 CONTINUE
      IF(QLK .LE. QL0)  GO TO 130
C
C     PRAUT
C
      QEXRHO = QLR - QL0RHO
      PRAUT = AMIN1(RAUT(CRAUT,QEXRHO),SLMAX)
C
  130 CONTINUE
      IF(.NOT. RAIN)  GO TO 150
C
C     PRACW
C
      PRACW = AMIN1(ACR1(QLK,QRR,CRACW,0.95),SLMAX)
C
  150 CONTINUE
      IF(QRK .EQ. 0.0 .OR. DQS .GE. 0.0)  GO TO 300
C
C     PREVP
C
      SW = QVK/QSW
      PREVP = AMIN1(REVP(SW,TSQ,QSW,QRR,CREVP),SRMAX)
  300 IF(TC .LT. 0.0)  GO TO 400
C
C     ***********************************
C     ***********************************
C     ****                           ****
C     ****     TC >= 0 PROCESSES     ****
C     ****                           ****
C     ****           GMLT            ****
C     ****                           ****
C     ***********************************
C     ***********************************
C
      IF(.NOT. GRAUP .OR. .NOT. RAIN)  GO TO 320
C
C     GACR CALLED FOR GMLT HERE
C
      PGACR = AMIN1(ACR3(VTG,VTR,QRR,QGR,CGACR,C3GACR),SRMAX)
C
C     GMLT:  GACW HAS ALREADY BEEN CALLED IF APPROPRIATE
C            GUARD AGAINST NEGATIVE VALUES AT TEMP CLOSE TO 0 DEG C
C
  320 IF(QGK .EQ. 0.0)  GO TO 330
C
      PGMLT = AMIN1(GMLT(TC,DQS0,QGR,PGACW,PGACR,CGMLT),SGMAX)
      PGMLT = AMAX1(PGMLT,0.0)
C
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****             T>=0             ****
C     ****                              ****
C     **************************************
C
  330 CONTINUE
      TSSG = -PGMLT
      TSSR = PGMLT + PRAUT + PRACW +PGACW - PREVP
      TTMP = -(HLTF*PGMLT+HLTC*PREVP)/CP
      GO TO 500
C
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****             T<0              ****
C     ****                              ****
C     **************************************
C
  400 CONTINUE
      TSSG = PGACW
      TSSR = PRAUT + PRACW - PREVP
      TTMP = (-HLTC*PREVP + HLTF*PGACW)/CP
C
C     **************************************
C     ****                              ****
C     ****     ADD SOURCES AND SINKS    ****
C     ****           ALL TEMP           ****
C     ****                              ****
C     **************************************
C
  500 CONTINUE
      TSSV = PREVP
      TSSL = -(PRAUT+PRACW+PGACW)
C
C     WRITE(6,111) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  111 FORMAT('0',30X,'T>=0  '//10X,'TSSV=',E16.8/
     *,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8,
     */10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
      TEM = 0.0
      DO 1112 K=1,6
 1112 TEM = TEM+TSS(K)
C     WRITE(6,1113)  TEM
 1113 FORMAT('0',10X,'TOTAL SOURCE AND SINK=',E16.8)
C
C
C     PRINT VALUES FOR EACH PROCESS
C
C     WRITE(6,2000)'ALL TEMP','PGACW',PGACW,'PSACW',PSACW,'PRAUT',PRAUT,
C    1   'PRACW',PRACW,'PRACS',PRACS,'PGACS',PGACS,'PREVP',PREVP
C     WRITE(6,2000)'TC > 0','PSMLT',PSMLT,'PWACS',PWACS,'PGMLT',PGMLT
C     WRITE(6,2000) 'TC < 0','PSFW',PSFW,'PSFI',PSFI,'PRACI',PRACI,
C    1   'PIACR',PIACR,'PGACI',PGACI,'PSACI',PSACI,'PSAUT',PSAUT,
C    2   'PGFR',PGFR,'PSACR',PSACR,'PGACR',PGACR,'PSSUB',PSSUB,
C    3   'PGAUT',PGAUT,'PGSUB',PGSUB,'PIDW',PIDW
C     WRITE(6,2000)'WET-DRY','PGWET',PGWET,'PGDRY',PGDRY,'PGWORD',PGWORD
      RETURN
      END
