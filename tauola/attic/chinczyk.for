      SUBROUTINE DPHNPI(DGAMT,PN,PR,PPI,JNPI)
C ----------------------------------------------------------------------
C IT SIMULATES MULTIPI DECAY IN TAU REST FRAME WITH
C Z-AXIS OPPOSITE TO NEUTRINO MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      PARAMETER (NMODE=4)
      COMMON / TAUNPI / CBRNPI       ,AMAS
     &                 ,KPI(6,NMODE) ,MULT(NMODE)
      REAL*4            CBRNPI(NMODE),AMAS(6,NMODE)
C
      REAL  PN(4),PR(4),PPI(4,6)
      REAL  PV(5,6),PT(4),UE(3),BE(3)
      REAL*4 RRR(6)
C
      DATA PI /3.141592653589793238462643/
C
C     PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.*A)
      PAWT(A,B,C)=SQRT(MAX(0.,(A**2-(B+C)**2)*(A**2-(B-C)**2)))/(2.*A)
C
C
C
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
C MASS OF VIRTUAL W
      ND=MULT(JNPI)
      PS=0.
      PHSPAC = 1./2.**5 /PI**2
      DO 4 I=1,ND
4     PS  =PS+AMAS(I,JNPI)
      CALL RANMAR(RR1,1)
      AMS1=PS**2
      AMS2=(AMTAU-AMNUTA)**2
C
CAM   FLAT PHASE SPACE  !!!
      AMX2=AMS1+   RR1*(AMS2-AMS1)
      AMX =SQRT(AMX2)
      AMW =AMX
      PHSPAC=PHSPAC * (AMS2-AMS1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX2)
      PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C W MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC * (4.*PI) * (2.*PR(3)/AMTAU)
 
 
C AMPLITUDE  (cf YS.Tsai Phys.Rev.D4,2821(1971)
C    or F.Gilman SH.Rhie Phys.Rev.D31,1066(1985)
C
        PXQ=AMTAU*PR(4)
        PXN=AMTAU*PN(4)
        QXN=PR(4)*PN(4)-PR(1)*PN(1)-PR(2)*PN(2)-PR(3)*PN(3)
 
        BRAK=2*(GV**2+GA**2)*(2*PXQ*QXN+AMX2*PXN)
     &      -6*(GV**2-GA**2)*AMTAU*AMNUTA*AMX2
 
CAM     Assume neutrino mass=0. and sum over final polarisation
C     BRAK= 2*(AMTAU**2-AMX2) * (AMTAU**2+2.*AMX2)
      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AMX2*SIGEE(AMX2,JNPI)
      DGAMT=1./(2.*AMTAU)*AMPLIT*PHSPAC
 
 
C   ISOTROPIC W DECAY IN W REST FRAME
      PHSMAX = 1.
      DO 200 I=1,4
  200 PV(I,1)=PR(I)
      PV(5,1)=AMW
      PV(5,ND)=AMAS(ND,JNPI)
C    COMPUTE MAX. PHASE SPACE FACTOR
      PMAX=AMW-PS+AMAS(ND,JNPI)
      PMIN=.0
      DO 220 IL=ND-1,1,-1
      PMAX=PMAX+AMAS(IL,JNPI)
      PMIN=PMIN+AMAS(IL+1,JNPI)
  220 PHSMAX=PHSMAX * PAWT(PMAX,PMIN,AMAS(IL,JNPI))/PMAX
 
CC----------------------------------------------------------------------C
C YESW  FACTOR RATE WAS MISSING (3.02.94)
CC----------------------------------------------------------------------C
      RATE = 1.
      AMX  = AMW
      DO 222 IL=1,ND-2
      AMS1=.0
      DO 223 JL=IL+1,ND
  223 AMS1=AMS1+AMAS(JL,JNPI)
      AMS1=AMS1**2
      AMX = AMX-AMAS(IL,JNPI)
      AMS2= AMX**2
      RATE = RATE * (AMS2-AMS1)
  222 CONTINUE
 
      PHSMAX=PHSMAX * RATE
 
  100 CONTINUE
CAM  GENERATE ND-2 EFFECTIVE MASSES
      PHS=1.
      RATE = 1.
      PHSPAC = 1./2.**(6*ND-7) /PI**(3*ND-4)
      AMX=AMW
      CALL RANMAR(RRR,ND-2)
      DO 230 IL=1,ND-2
      AMS1=.0
      DO 231 JL=IL+1,ND
  231 AMS1=AMS1+AMAS(JL,JNPI)
      AMS1=AMS1**2
      AMS2=(AMX-AMAS(IL,JNPI))**2
      RATE= RATE * (AMS2-AMS1)
      RR1=RRR(IL)
      AMX2=AMS1+  RR1*(AMS2-AMS1)
      AMX=SQRT(AMX2)
      PV(5,IL+1)=AMX
      PHSPAC=PHSPAC * (AMS2-AMS1)
      IF(PV(5,IL).LT.(PV(5,IL+1)+AMAS(IL,JNPI))) GOTO 100
      PA=PAWT(PV(5,IL),PV(5,IL+1),AMAS(IL,JNPI))
      PHS   =PHS    *PA/PV(5,IL)
  230 CONTINUE
C ZW 3.02.94 ¦¦¦¦
      PA=PAWT(PV(5,ND-1),AMAS(ND-1,JNPI),AMAS(ND,JNPI))
      PHS   =PHS    *PA/PV(5,ND-1)
C
      PHS = PHS*RATE
 
      CALL RANMAR(RN,1)
      IF(RN*PHSMAX.GT.PHS) GO TO 100
 
C...PERFORM SUCCESSIVE TWO-PARTICLE DECAYS IN RESPECTIVE CM FRAME
  280 DO 300 IL=1,ND-1
      PA=PAWT(PV(5,IL),PV(5,IL+1),AMAS(IL,JNPI))
      CALL RANMAR(RRR,2)
      UE(3)=2.*RRR(1)-1.
      PHI=2.*PI*RRR(2)
      UE(1)=SQRT(1.-UE(3)**2)*COS(PHI)
      UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI)
      DO 290 J=1,3
      PPI(J,IL)=PA*UE(J)
  290 PV(J,IL+1)=-PA*UE(J)
      PPI(4,IL)=SQRT(PA**2+AMAS(IL,JNPI)**2)
      PV(4,IL+1)=SQRT(PA**2+PV(5,IL+1)**2)
      PHSPAC=PHSPAC *(4.*PI)*(2.*PA/PV(5,IL))
  300 CONTINUE
C...LORENTZ TRANSFORM DECAY PRODUCTS TO TAU FRAME
      DO 310 J=1,4
  310 PPI(J,ND)=PV(J,ND)
      DO 340 IL=ND-1,1,-1
      DO 320 J=1,3
  320 BE(J)=PV(J,IL)/PV(4,IL)
      GAM=PV(4,IL)/PV(5,IL)
      DO 340 I=IL,ND
      BEP=BE(1)*PPI(1,I)+BE(2)*PPI(2,I)+BE(3)*PPI(3,I)
      DO 330 J=1,3
  330 PPI(J,I)=PPI(J,I)+GAM*(GAM*BEP/(1.+GAM)+PPI(4,I))*BE(J)
      PPI(4,I)=GAM*(PPI(4,I)+BEP)
  340 CONTINUE
C
      RETURN
      END
 
