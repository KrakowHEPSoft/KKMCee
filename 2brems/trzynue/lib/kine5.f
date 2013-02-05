
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
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),
     $          pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(17)

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
      CALL VARRAN(RRR,17)
! flat phase space
      ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
      ams2=(amtax-amnut2)**2
!     amx4=ams1+rrr(11)*(ams2-ams1)
!      phspac=phspac*(AMS2-AMS1)
        if (amtax.gt.ams1) goto 5
        amro=100
        gamro=100
        prob1=.5
        prob2=.5
        PROB3=.0
        PROB4=.0
        rr2=rrr(11)
        ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
        ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
        IF (RRR(14).LT.PROB1) THEN
         AM2SQ=AMS1+   RR2*(AMS2-AMS1)
         AM2 =SQRT(AM2SQ)
        elseIF (RRR(14).LT.(PROB1+PROB2)) THEN  
         B=LOG(AMS1)
         A=LOG(AMS2)
         AM2SQ=AMS2*EXP((B-A)*RR2)
         am2sq=ams2+ams1-am2sq
         AM2 =SQRT(AM2SQ)     
        ELSEIF (RRR(14).LT.(PROB1+PROB2+PROB3)) THEN
         ALP=ALP1+RR2*(ALP2-ALP1)
         AM2SQ=AMRO**2+AMRO*GAMRO*TAN(ALP)
         AM2 =SQRT(AM2SQ)
        ELSE
         n=1
          if(n.eq.1) then
         AM2SQ=AMS1/(1D0-RR2*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQ=AMS1/sqrt(1D0-RR2*(1-(ams1/ams2)**n))
          else
         AM2SQ=AMS1*(1D0-RR2*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
         AM2 =SQRT(AM2SQ)
         if (am2sq.gt.ams2) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.gt.ams2) stop
         if (am2sq.lt.ams1) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.lt.ams1) stop

        ENDIF

        XJ1=(AMS2-AMS1)
         B=LOG(AMS1)
         A=LOG(AMS2)
         am2sqx=ams2+ams1-am2sq
        xj2=AM2SQx*(A-B)
        xj3=((AM2SQ-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
        xj3=xj3*(ALP2-ALP1)
        n=1
        xj4=am2SQ**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)
!        sum=Sum+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
!        sum2=Sum2+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)**2
!        enddo
!        sum=sum/nn
!        sum2=sum2/nn
!        err=sqrt((sum2-sum**2)/nn)
!        write(*,*) sum,'+-',err
!        write(*,*) '28761.2547837270613 +- 0'
!        stop
        PHSPAC=PHSPAC/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
 5     continue

C>>> 1 soft photon
C>>>          RRR(11)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in pn2(4)
C GENERATING MASS of l+l- gamam gamma 
c...  xk0 defines min. photon pn2 energy in lab 
        ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
        XK0=0.0051D0 
        AMS2=AMTAx**2

          RR1=RRR(11)
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0)
          XL0=LOG(2*XK0)
          XK=EXP(XL1*RR1+XL0)
          AM3SQ=(1-XK)*AMS2
          AM3 =SQRT(AM3SQ)
          PHSPAC=PHSPAC*AMS2*XL1*XK
        IF(PHSPAC.EQ.0D0) GOTO 900
       
      amtau=sqrt(am3sq)


C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in (PIM1+PIM2)**2 
C GENERATING MASS of l+l- pair
        AMS1=(AMP2+AMP3)**2       
        AMS2=(AMTAU-AMP1-AMNUTA)**2 
        ALP1=ATAN((AMS1-AMaz**2)/AMaz/GAMmz)
        ALP2=ATAN((AMS2-AMaz**2)/AMaz/GAMmz)
        prob1=0.5
        prob2=0.5
        IF(RRR(9).LE.prob1) THEN
         am2sq=ams1+(ams2-ams1)*rrr(1)
        ELSE
         ALP=ALP1+RRr(1)*(ALP2-ALP1)
         AM2sq=AMaz**2+AMaz*GAMmz*TAN(ALP) 
        ENDIF
        XJAC1=ams2-ams1
        XJAC2=((AM2sq-AMaz**2)**2+(AMaz*GAMmz)**2)
     $       /(AMaz*GAMmz)*(ALP2-ALP1)
        XJAC=1d0/(prob1/XJAC1+prob2/XJAC2)
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900

C MASS OF ll gam
        AMS1=(am2+amp1)**2
        AMS2=(AMTAU-AMnuta)**2
        AM3SQ=AMS1+   RRR(2)*(AMS2-AMS1)
        AM3 =SQRT(AM3SQ)
        PHSPAC=PHSPAC*(AMS2-AMS1)

* AM2 RESTFRAME, DEFINE PIM2 AND PIM3
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PIM3)
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
* NOW boost TO THE am3 REST FRAME
      paa4=1./(2*am3)*(am3**2-amp1**2+am2**2)                    
      paa3= sqrt(abs(paa4**2-am2**2))                                            
      exe=(paa4+paa3)/am2
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)
      eee=0
      do k=1,3
       pim1(k)=-pim2(k)-pim3(k)
       eee=eee+pim1(k)**2
      enddo
      pim1(4)=sqrt(eee)
        PHSPAC=PHSPAC*(4*PI)*(2*pim1(4)/AM3)
* ALL ROTATED IN THE am3 REST FRAME
      THET =ACOS(-1.D0+2*RRR(16))
      PHI = 2*PI*RRR(17)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)

* AM4/amtau RESTFRAME, DEFINE PN
        ENQ1=(Amtau**2-AM3**2+AMNUTA**2)/(2*AMtau)
        ENQ2=(AMtau**2+AM3**2-AMNUTA**2)/(2*AMtau)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMtau)
* PI MINUS MOMENTUM IN RHO REST FRAME
* NOW boost TO THE am4 REST FRAME
      paa4=1./(2*amtau)*(amtau**2-amnuta**2+am3**2)                    
      paa3= sqrt(abs(paa4**2-am3**2))                                            
      exe=(paa4+paa3)/am3
      CALL BOSTd3(EXE,PIm1,PIm1)
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)

        pn(1)=0d0
        pn(2)=0d0
        pn(3)=-pppi
        PN(4)=ENQ1
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

        prob1=1.
        prob2=.0
        prob3=.0
        prob4=0
        prob5=0
        EPS=(AMp3/AMtax)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
      IF    (RRR(15).lt.PROB1) then 
       THET =ACOS(-1.D0+2*RRR(12))
       CTHET=COS(THET)
      elseIF(RRR(15).lt.(PROB1+PROB2)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=-(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
!         CTHET=-cthet
        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)

        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3+PROB4)) then
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=AM2SQX-2D0+sqrt(1d0-eps)
        THET =ACOS(CTHET)
      else
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=-AM2SQX+2D0-sqrt(1d0-eps)
        THET =ACOS(CTHET)
      endif
      if (cthet**2.gt.1d0) then
       cthet=cthet/cthet**2
       write(*,*) 'cthet error -- arbi action'
       write(*,*) cthet,rrr(12),rrr(15)
       write(*,*) ams1,ams2,am2sq
        THET =ACOS(CTHET)
      endif
      eta1=1+eps+cthet
      eta2=1+eps-cthet
      xx=eps
      beta=sqrt(1d0-eps)
      xx=eps
      xlog=-log((1+beta)**2/xx)
      xlog1=-log(16D0/xx)
      ct=-cthet

      xccos1=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))
      ct=cthet
      xccos2=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))

      xccos1=1d0/(XL1/2*ETA1)
      xccos2=1d0/(XL1/2*ETA2)

        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
        n=1
        AM2SQX= CTHET+2D0-sqrt(1d0-eps)
        xj4=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
        n=1
        AM2SQX=-CTHET+2D0-sqrt(1d0-eps)
        xj5=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
       
      FF=PROB1/1d0+PROB2*xccos1+PROB3*xccos2+PROB4/XJ4+PROB5/XJ5

      PHSPAC=PHSPAC/FF                                                   
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
!      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)


C FINAL WEIGHT
      WT = PHSPAC
C THE STATISTICAL FACTOR FOR IDENTICAL gammas 
C is replaced with ordering in p_T
      pt3=pn2(1)**2+pn2(2)**2
      pt2=pn(1)**2+pn(2)**2
      pt1=pim1(1)**2+pim1(2)**2
      if (pt1.lt.pt2) wt=0
      if (pt2.lt.pt3) wt=0
      RETURN
 900  WT=0D0

      END








