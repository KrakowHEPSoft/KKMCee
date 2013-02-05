      SUBROUTINE KINE5IIZ0(
     # AMTAx,AMP3,AMP2,AMP1,AMNUTA,amnut2,PIM3,PIM2,PIM1,PN,pn2,WT)
C     *******************************************************************
C generator of 5 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA,amnut2  - masses of particles
C         PIM1,PIM2,PIM3,PN,PN2  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C         pn2-flat new additional particle. Not yet integrated.
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(10)
      COMMON /POMOC2/ Y,Z,DELZ,DEL1,DEL2,BETA,C,SPHI,C1,S1,CJAC

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
!      PHSPAC=1.D0/2**17/PI**8
      PHSPAC=1.D0/2**23/PI**11

C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAx
C
      CALL VARRAN(RRR,13)
! flat phase space
      ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
      ams2=(amtax-amnut2)**2
      amx4=ams1+rrr(11)*(ams2-ams1)
      amtau=sqrt(amx4)
      phspac=phspac*(AMS2-AMS1)
C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in (PIM1+PIM2)**2 
C GENERATING MASS of l+l- pair
        AMS1=(AMP2+AMP3)**2       
c...        AMS2=AMTAU**2-(AMP1+AMNUTA)**2 
c...  xk0 defines min. photon energy in photon-photon rest frame
        XK0=0.01D0 
        AMS2=AMTAU**2-(2D0*XK0*AMTAU)**2       
        Y1  = 1D0-AMS1/AMTAU**2
        Y0  = 1D0-AMS2/AMTAU**2
        IF(Y1.LT.Y0) Y1=Y0
        AZ0=AMAZ/AMTAU
        BZ0=AMAZ*GAMMZ/AMTAU**2
        XP2=1D0/BZ0*
     #   (ATAN((1D0-Y0-AZ0**2)/BZ0)-ATAN((1D0-Y1-AZ0**2)/BZ0))
        XP3=DLOG((1D0-Y0)/(1D0-Y1))
        IF(RRR(9).LE.(XP2)/(XP2+XP3)) THEN
          Y=1D0-AZ0**2-BZ0*TAN(
     #            +RRR(1) *ATAN((1D0-Y1-AZ0**2)/BZ0)
     #      +(1D0 -RRR(1))*ATAN((1D0-Y0-AZ0**2)/BZ0)
     #                       )
        ELSEIF(RRR(9).GT.(XP2)/(XP2+XP3)) THEN
           Y=1D0-(1D0-Y0)*((1D0-Y1)/(1D0-Y0))**RRR(1)
        ENDIF
        XJAC2=ABS(BZ0*( ATAN((1D0-Y1-AZ0**2)/BZ0)
     #           -ATAN((1D0-Y0-AZ0**2)/BZ0)   )
     #         /COS(ATAN((1D0-Y-AZ0**2)/BZ0))**2)
        XJAC3=ABS((1D0-Y)*DLOG((1D0-Y0)/(1D0-Y1)))
        XJAC=(XP2+XP3)/(XP2/XJAC2+XP3/XJAC3)
        AM3SQ=(1D0-Y)*AMTAU**2
        AM3=SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMTAU**2*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900
* AM3 RESTFRAME, DEFINE PIM2 AND PIM3
        ENQ1=(AM3SQ-AMP2**2+AMP3**2)/(2*AM3)
        ENQ2=(AM3SQ+AMP2**2-AMP3**2)/(2*AM3)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM3)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PIM3)
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
C MASS OF gam-gam pair 
        AMS1=(AMP1+AMNUTA)**2
        AMS2=(AMTAU-AM3)**2
        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
        PHSPAC=PHSPAC*(AMS2-AMS1)
* AM2 RESTFRAME, DEFINE PIM1 AND PN
        ENQ1=(AM2SQ-AMP1**2+AMNUTA**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP1**2-AMNUTA**2)/(2*AM2)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PN)
        PN(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 40 I=1,3
 40     PIM1(I)=-PN(I)
        PIM1(4)= ENQ2
* NOW TO THE TAU REST FRAME, DEFINE AM2 AND AM3 MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AM2**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
      PBB(1)=0.D0
      PBB(2)=0.D0
      PBB(4)=1.D0/(2*AMTAU)*(AMTAU**2+AM2**2-AM3**2)
      PBB(3)=-PAA(3)
* ALL PIONS BOOSTED  TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO AM2 MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      EXE=(PBB(4)+PBB(3))/AM2
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(5))
      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)



* now to the tau rest frame, define paa and neutrino momenta            
* paa  momentum                                                         
      paa(1)=0                                                          
      paa(2)=0                                                          
      paa(4)=1./(2*amtax)*(amtax**2-amnut2**2+amtau**2)                    
      paa(3)= sqrt(abs(paa(4)**2-amtau**2))                                            
      phspac=phspac*(4*pi)*(2*paa(3)/amtax)                             
* tau-neutrino momentum                                                 
      pn2(1)=0                                                           
      pn2(2)=0                                                           
      pn2(4)=1./(2*amtax)*(amtax**2+amnut2**2-amtau**2)                    
      pn2(3)=-paa(3)  

      exe=(paa(4)+paa(3))/amtau                                           
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)
                                                  
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END
