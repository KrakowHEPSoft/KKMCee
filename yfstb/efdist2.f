

      SUBROUTINE yostes(mode,wtx2)   
C     ***********************    
C          ===================================      
C          =========== yostes  ===============      
C          ===================================      
***
* New test of amplitude with L3 cut, amplitude by Scott Yost 
*** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      PARAMETER(PI=3.1415926535897932D0)        
      COMMON / INOUT  / NINP,NOUT           
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT     
      COMMON / INPARM / CMSENE,AMFIN,KEYRAD      
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)      
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2        
      COMMON /trmkey/ DCSkey, LLkey, Zkey
      integer         DCSkey, LLkey, Zkey
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /,/ WGTD2  /,
     $     /trmkey/
      SAVE IEV,ievac,p1,p2
      DIMENSION ph1(4),ph2(4),sph(4),p1(4),p2(4)
      DATA iev,ievac /0,0/

      IF (mode.eq.-1) THEN
!     ====================
      DCSkey = 100
      LLkey  = 100
      Zkey   = 1

c      CALL AMPLIT(-1)           
      IEV  =0
      IEVac=0
      p1(1)=0
      p1(2)=0
      p1(3)= cmsene/2
      p1(4)= cmsene/2
      p2(1)=0
      p2(2)=0
      p2(3)=-cmsene/2
      p2(4)= cmsene/2

      ELSEIF (mode.eq.0) THEN
!     =======================

      if(nphot.ne.2) RETURN

      DO 30 k=1,4
      ph1(k) = sphot(1,k)
      ph2(k) = sphot(2,k)
      sph(k) = ph1(k)+ph2(k)
 30   CONTINUE
      sumeng  = sph(4)
      if(sumeng.gt.0.01) RETURN

c      CALL amplit( 0)
      wtx2 = wt2

      CALL twopho(dcs, dcsll, p1, qf1, p2, qf2, ph1, ph2)

C.. cross section compact formula from YFS version january 1990 
      CALL dfinap2(cmsene, p1, qf1, p2, qf2, ph1, ph2, compfi)  

      ievac=ievac+1
      IF(ievac.lt.10) THEN
      write(6,*) "--------------- iev,ievac",iev,ievac     
c      CALL dumps(6)
      write(6,"(a,3f15.10)") "ph1(4)+ph2(4)=",sumeng
      write(6,*) "COMPFI=",COMPFI
      write(6,*) "DCS   =",DCS,    DCS/COMPFI
      write(6,*) "DCSll =",DCSll,DCSll/COMPFI
      ENDIF

      ENDIF

      END 

      SUBROUTINE dfinap2(cmsene,p1,p2,q1,q2,pk1,pk2,xection)  
!     ******************************************************
! Compact form for double final bremstr cross section from yfs3 M.C.
! Note BHLUMI notation for  p1,q1,p2,q2.
C     ***************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      PARAMETER(pi=3.1415926535897932d0, alfinv=137.03604d0)        
      DIMENSION p1(4),q1(4),p2(4),q2(4),pk1(4),pk2(4),xx(4)    
        
      alfa=1d0/alfinv       
      s0=cmsene**2  
      s1=(p2(4)+q2(4))**2-(p2(1)+q2(1))**2-(p2(2)+q2(2))**2     
     $   -(p2(3)+q2(3))**2  
      DO 10 i=1,4   
      xx(i)=p1(i)+q1(i)     
  10  CONTINUE  
        
      CALL efdist2(xx,p1,q1,p2,q2,pk1,pk2,dist2)        
      xection=
     $   2**12 *pi**8 *alfa**2
     $   *S1/S0
     $   *4*(ALFA/4D0/PI**2)**2
     $   *4*(2*pi*alfa)**2
     $   *dist2    
      END       


      SUBROUTINE EFDIST2(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)  
C     ***********************************************   
C Provides double bremsstrahlung distribution - FINAL state brem.       
C INPUT:  P1,P2,Q1,Q2,PH1,PH2, four momenta 
C OUTPUT: DIST2     double bremsstrahlung distribution  
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)     
      DIMENSION PR1(4),PR2(4),PH1R(4),PH2R(4),QR1(4),QR2(4)     
        
      CALL EREDUZ2(QQ,Q1,Q2,PH1,PH2,QR1,QR2,PH1R,PH2R)  
      CALL EREDUZ0(QQ,P1,P2,PR1,PR2)        
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2       
C infrared factors from reduced momenta     
C double bremsstrahlung Xsect in next-to-leading log approx.    
      CALL EGSFIN2(Q1,Q2,PH1,PH2,GF1,GF2)   
      CALL EGTHET1(QR1,QR2,PR1,COSTH1,COSTH2)   
      ANDI11= BORNVV(SVAR1,COSTH1)  
      ANDI12= BORNVV(SVAR1,COSTH2)  
      DIST2 =   GF1*ANDI11+   GF2*ANDI12    
      END       


      SUBROUTINE EGTHET1(P1,P2,Q1,COSTH1,COSTH2)        
C     ***************************************** 
C Calculates CosTh1 and CosTh2 between BEAM amd FINAL   
C fermion momenta in final fermion rest frame Q1(4)+Q2(4)=0     
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION P1(*),P2(*),Q1(*)   
      COSTH1 = (P1(1)*Q1(1) +P1(2)*Q1(2) +P1(3)*Q1(3))  
     $    /SQRT((Q1(1)**2 +Q1(2)**2 +Q1(3)**2)  
     $ *(P1(1)**2 +P1(2)**2 +P1(3)**2)) 
      COSTH2 =-(P2(1)*Q1(1) +P2(2)*Q1(2) +P2(3)*Q1(3))  
     $    /SQRT((Q1(1)**2 +Q1(2)**2 +Q1(3)**2)  
     $ *(P2(1)**2 +P2(2)**2 +P2(3)**2)) 
      END       
      SUBROUTINE EREDUZ0(QQ,P1,P2,PR1,PR2)  
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION QQ(4),P1(4),P2(4),PR1(4),PR2(4) 
      DO 20 K=1,4   
      PR1(K)=P1(K)  
 20   PR2(K)=P2(K)  
      END       
      SUBROUTINE EREDUZ2(QQ,P1,P2,PH1,PH2,PR1,PR2,PH1R,PH2R)    
C     *****************************************************     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION QQ(*), P1(*),  P2(*),  PH1(*),  PH2(*)  
      DIMENSION        PR1(*), PR2(*), PH1R(*), PH2R(*) 
      DO 20 K=1,4   
      PH1R(K)=PH1(K)        
      PH2R(K)=PH2(K)        
      PR1(K)=P1(K)  
 20   PR2(K)=P2(K)  
      END       

      SUBROUTINE EGSFIN2(P1,P2,PH1,PH2,F1,F2)   
C     **************************************    
C CALCULATES INGREDIENTS FOR REAL DOUBLE PHOTON DIFF. XSECTION  
C     ***************************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION PH1(*),PH2(*),P1(*),P2(*)   
C       
      WM (A  )=     (1D0-A)**2      
      WMS(A,B)=     ((1D0-A)**2+(1D0-B)**2) 
      WWM(A,B)=     
     $   1D0-AM*2D0*(1D0-A)*(1D0-B)/((1D0-A)**2+(1D0-B)**2)*(A/B+B/A)   
C       
c       call dumpt(6,"===p1===",p1)
c       call dumpt(6,"   p2   ",p2)
c       call dumpt(6,"--ph1---",ph1)
c       call dumpt(6,"  ph2   ",ph2)
      PP = P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)      
      AM2= P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2  
      AM = AM2/(2D0*PP)     
      BB1=ABS(P2(4)*PH1(4)-P2(1)*PH1(1)-P2(2)*PH1(2)-P2(3)*PH1(3))/PP   
      AA1=ABS(P1(4)*PH1(4)-P1(1)*PH1(1)-P1(2)*PH1(2)-P1(3)*PH1(3))/PP   
      BB2=ABS(P2(4)*PH2(4)-P2(1)*PH2(1)-P2(2)*PH2(2)-P2(3)*PH2(3))/PP   
      AA2=ABS(P1(4)*PH2(4)-P1(1)*PH2(1)-P1(2)*PH2(2)-P1(3)*PH2(3))/PP   
      AA1P= AA1/(1D0+AA2)   
      BB1P= BB1/(1D0+BB2)   
      AA2P= AA2/(1D0+AA1)   
      BB2P= BB2/(1D0+BB1)   
      A1  = AA1/(1+AA1+BB1) 
      A2  = AA2/(1+AA2+BB2) 
      B1  = BB1/(1+AA1+BB1) 
      B2  = BB2/(1+AA2+BB2) 
      A1P = AA1P/(1+AA1P+BB1P)      
      A2P = AA2P/(1+AA2P+BB2P)      
      B1P = BB1P/(1+AA1P+BB1P)      
      B2P = BB2P/(1+AA2P+BB2P)      
      SFAC1  =  2D0/(PP*AA1*BB1)*WWM(A1,B1) 
      SFAC2  =  2D0/(PP*AA2*BB2)*WWM(A2,B2) 
      IF((A1+B1).GT.(A2+B2)) THEN   
        X1=WM (A1   )*WMS(A2P,B2P) +WM (A1P    )*WMS(A2,B2)     
        X2=WM (   B1)*WMS(A2P,B2P) +WM (    B1P)*WMS(A2,B2)     
      ELSE      
        X1=WM (A2   )*WMS(A1P,B1P) +WM (A2P    )*WMS(A1,B1)     
        X2=WM (   B2)*WMS(A1P,B1P) +WM (    B2P)*WMS(A1,B1)     
      ENDIF     
      F1 = X1*SFAC1*SFAC2/8D0       
      F2 = X2*SFAC1*SFAC2/8D0       
C.. correction E.W. november 1989................................       
C.. temporary solution      
C.. this correction reconstructs double collinear limit with 50%        
C.. and affects below photon-fermion angle  <0.1 ammi/ene       
C.. without correction error in this limit 1000%        
      SFAC1  =  2D0/(PP*AA1*BB1)    
      SFAC2  =  2D0/(PP*AA2*BB2)    
      DELT=(AM2/(2D0*PP))**2*(B2**2*A1**2+A2**2*B1**2)* 
     #  ( B1*B2/(A1*A2)/(A1+A2)**2  
     #   +A1*A2/(B1*B2)/(B1+B2)**2  )       
      WMINF=2D0*DELT/(X1+X2)        
      WMM=WWM(A1,B1)*WWM(A2,B2)+WMINF       
      F1 = X1*SFAC1*SFAC2/8D0*WMM   
      F2 = X2*SFAC1*SFAC2/8D0*WMM   
C...end of correction............................................       
      END       

      FUNCTION BORNVV(SVARI,COSTHE) 
C     ***********************************   
C THIS ROUTINE PROVIDES BORN DIFFERENTIAL CROSS SECTION 
C a version without COMPLEX*16      
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
C     COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF     
      COMMON / BHPAR1 / CMS,AMFIN   
      COMMON / BHPAR2 / CMSENE,AMEL 
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE    
      COMMON /COEFF/ VE,AE  
      SAVE / WEKING /,/ BHPAR1 /,/ BHPAR2 /,/ WEKINP /,/COEFF/
        
      VF=VE     
      AF=AE     
      QE= -1D0  
      QF= -1D0  
C     AA= 4D0*SQRT(SINW2*(1D0-SINW2))       
C     VE= (-1D0+4*SINW2)/AA 
C     AE= 1D0/AA    
C     VF= (-1D0+4*SINW2)/AA 
C     AF= 1D0/AA    
      S = SVARI     
      CHI2 = S**2/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)      
      RECHI=(S-AMAZ**2)*S/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)      
      XE= VE**2 +AE**2      
      XF= VF**2 +AF**2      
      YE= 2*VE*AE   
      YF= 2*VF*AF   
      FF0= QE**2*QF**2 +2*RECHI*QE*QF*VE*VF +CHI2*XE*XF 
      FF1=     +2*RECHI*QE*QF*AE*AF +CHI2*YE*YF 
C     BORN    = (1D0+ COSTHE**2)*FF0 +2D0*COSTHE*FF1    
      BORN    = (1D0+ COSTHE**2     
     1 +4D0*(AMEL**2+AMFIN**2)/S*(1D0-COSTHE**2)        
     1 +16*AMEL**2*AMFIN**2/S**2*COSTHE**2)*FF0 
     1 +2D0*COSTHE*FF1   
![[[[[ gamma only   
      BORN    = 1D0+ COSTHE**2     
C************   
C THIS IS A BIT CRUDE METHOD OF INTRODUCING THRESHOLD BEHAVIOUR 
      IF(    SVARI.LE. 4D0*AMFIN**2) THEN   
        THRESH=0D0  
      ELSEIF(SVARI.LE.16D0*AMFIN**2) THEN   
        AMX2=4D0*AMFIN**2/SVARI     
        THRESH=SQRT(1D0-AMX2)*(1D0+AMX2/2D0)    
      ELSE      
        THRESH=1D0  
      ENDIF     
      BORNVV= BORN*THRESH   
      END       
        
