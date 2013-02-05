*\input fortex \plain
!!!!!!!!!!! Z is off!!!!!    
 
c===========================================================
c================== St.J. Corrections JULY 92 ==============
c 1. New glibk version 1.05 is attached,
c 2. SAVE added all over through subprograms -- it was checked
c    that results do not change when compilation option -K removed.
c 3. Statistics 100k events on HP-720 requires 8 minutes

      PROGRAM  MAIN       
C     *************       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      common / cglib / b(20000)
      CHARACTER*6 DNAME   
      COMMON / INOUT  / NINP,NOUT     
        
      CALL GLIMIT(20000)  
      NINP= 5    
      NOUT=16    
      OPEN(NOUT,file="./tb.output")
      CALL GOUTPU(NOUT)     
C     ********** 
      CALL SPLOT 
C     ********** 
      CLOSE(NOUT)         
C ------------WRITING HISTOS ON THE DISK ------------------------      
      NOUTH=10   
      OPEN(NOUTH,file="./tb.hst")
      CALL GRFILE(NOUTH,DNAME,'N') 
      CALL GROUT( 0,ICY,' ')       
      CALL GREND(DNAME)   
C ------------THE END OF HISTO WRITING -------------------------       
      END        
        
      LOGICAL FUNCTION KONIEC(IEV,NEVLIM,LIMTIM)    
C     ******************************************    
      TLIMIT=LIMTIM    
C THIS IS SLAC TIMER   
CCC   CALL LEFT1A(TIMLFT)        
C THIS IS CERN TIMER   
cccccc      CALL TIMEL(TIMLFT)
      TIMLFT = 1E15
      KONIEC = (IEV.GE.NEVLIM).OR.(TIMLFT.LT.TLIMIT)         
      END     
    
      SUBROUTINE SPLOT 
C     **************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT      
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / INOUT  / NINP,NOUT   
      SAVE / MOMSET /,/ INPARM  /,/ INOUT  /
      DIMENSION  XPAR(100),NPAR(100),XDMEM(100)       
      CHARACTER*80 TESTIT
      LOGICAL KONIEC   
    
      LIMTIM= 20       
 3000 FORMAT(8I2)      
 3001 FORMAT(I10)      
 3002 FORMAT(F10.0)    
 3010 FORMAT(A80)
      READ( NINP,3001)  NEVT
      WRITE( 6,*)       NEVT,' REQUESTED EVENTS '
      READ( NINP,3010)  TESTIT
      WRITE(NOUT,3010)  TESTIT
      READ( NINP,3000)  KAT1,KAT2,KAT3,KAT4,KAT5,KAT6,KAT7,KAT8
      READ( NINP,3001)  KEYBRM     
      READ( NINP,3002)  CMSENE,SINW2,AMAZ,GAMMZ       
      READ( NINP,3002)  AMFIN,VVMIN,VVMAX   
      READ( NINP,3001)  KEYRED,KEYWGT,KEYZET          
C Control output
      WRITE(NOUT,'(8A6/8I6)')
     $ 'KAT1','KAT2','KAT3','KAT4','KAT5','KAT6','KAT7','KAT8',
     $  KAT1 , KAT2 , KAT3 , KAT4 , KAT5 , KAT6 , KAT7 , KAT8
      WRITE(NOUT,'(4A6/4I6)')
     $ 'KAT1','KAT2','KAT3','KAT4',
     $  KAT1 , KAT2 , KAT3 , KAT4 
      WRITE(NOUT,'(4A12/4I12)')
     $  'NEVT','KEYBRM','KEYWGT','KEYZET',
     $   NEVT,  KEYBRM , KEYWGT , KEYZET
      WRITE(NOUT,'(7A12/7F12.6)')
     $ 'CMSENE','SINW2','AMAZ','GAMMZ','AMFIN','VVMIN','VVMAX',
     $  CMSENE , SINW2 , AMAZ , GAMMZ , AMFIN , VVMIN , VVMAX
      NPAR(1)=KEYBRM   
      NPAR(2)=KEYRED   
      NPAR(3)=KEYWGT   
      NPAR(4)=KEYZET   
      XPAR(1)=CMSENE
      XPAR(2)=AMAZ     
      XPAR(3)=SINW2    
      XPAR(4)=GAMMZ    
      XPAR(5)=AMFIN    
      XPAR(6)=VVMIN    
      XPAR(7)=VVMAX    

      CALL EXPAND(-1,XPAR,NPAR)    
C     *************************    
      IF(KAT1.EQ.1) CALL ROBOL1(-1)       
      IF(KAT2.EQ.1) CALL ROBOL2(-1)       
      IF(KAT3.EQ.1) CALL ROBOL3(-1)       
      IF(KAT4.EQ.1) CALL ROBOL4(-1) 

C----------- TIME AND EVENT LIMITS FOR THE MAIN LOOP ----------------- 
      IEV=0   
C----------------------------------------------------------------------
      DO 100 IJ=1,10000000      
      DO 100 JK=1,10000000      
      IEV=IEV+1        
      IF(MOD(IEV,5000).EQ.1) WRITE(6,*) '  IEV= ',IEV        
      CALL EXPAND( 0,XPAR,NPAR)    
C     ********************************  
      IF(IEV.LE.10) CALL DUMPS(6)         
CCC   IF(IEV.LE.10) CALL DUMPI(NOUT)      
CCC   IF(IEV.LE.10) CALL DUMPF(NOUT)      
      IF(IEV.LE.25) CALL DUMPS(NOUT)      
      IF(KONIEC(IEV,NEVT,LIMTIM)) CALL DUMPS(NOUT)      
      IF(IEV.LE. 5) CALL DUMPBT(NOUT)     
      IF(KAT1.EQ.1) CALL ROBOL1( 0)       
      IF(KAT2.EQ.1) CALL ROBOL2( 0)       
      IF(KAT3.EQ.1) CALL ROBOL3( 0)       
      IF(KAT4.EQ.1) CALL ROBOL4( 0)       
C------------------- EXIT FROM THE MAIN LOOP --------------------------
      IF(KONIEC(IEV,NEVT,LIMTIM)) GOTO 103          
  100 CONTINUE         
  103 CONTINUE         
      write(6,*)  ' =======> generation finished'
C-----------------------------------------------------------------------
      CALL EXPAND( 2,XPAR,NPAR)    
C Store usefull information in first histo ID=1
      NEVT   = NPAR(10)   
      XSECNB = XPAR(10)   
      XSCRNB = XPAR(20)   
      ERREL  = XPAR(21)   
      XSECR  = XPAR(22)   
      DO 500 K=1,100
 500  XDMEM(K)=0D0
      XDMEM( 1) = NEVT
      XDMEM(10) = XSECNB
      XDMEM(20) = XSCRNB
      XDMEM(21) = ERREL
      XDMEM(30) = CMSENE
      CALL GBOOK1(1,' Memory  $',100 ,0D0,1D0)   
      CALL GPAK(  1, XDMEM)
C-----------------------------------------------------------------------
      IF(KAT1.EQ.1) CALL ROBOL1(1)        
      IF(KAT2.EQ.1) CALL ROBOL2(1)        
      IF(KAT3.EQ.1) CALL ROBOL3(1)        
      IF(KAT4.EQ.1) CALL ROBOL4(1)       
      END     

      SUBROUTINE ROBOL1(MODE)    
C     ***********************    
C          ===================================      
C          ============ ROBOL1 ===============      
C          ===================================      
***
* New test of amplitude with L3 cut, amplitude by Scott Yost 
*** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      PARAMETER(PI=3.1415926535897932D0)        
      COMMON / INOUT  / NINP,NOUT           
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT     
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100) 
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /,/ WGTD2  /
      SAVE iev
      DIMENSION ph1(4),ph2(4),sph(4)
C         
      IF(MODE.EQ.-1) THEN     
C     ===================     
C---------------------        
      WRITE(NOUT,*) '======ROBOL1 initialization==============='
      WRITE(   6,*) '======ROBOL1 initialization==============='

      CALL amptes(-1,dum1,dum2,dum3,dum4)

* booking
      call gbook1(1000,"2-nd phot. energy, wt      $", 30, 1d0,46d0)
      call gbook1(1001,"2-nd phot. energy, wt*R    $", 30, 1d0,46d0)
      call gbook1(1002,"2-phot. mass, cut, wt      $", 30, 5d0,95d0)
      call gbook1(1003,"2-phot. mass, cut, wt*R    $", 30, 5d0,95d0)
**********
      call gbook1(2000,"ph-ph mass   no cut,       $",100, 0d0,100d0)
      call gbook1(2001,"ph-ph mass with cut,       $",100, 0d0,100d0)
      call gbook1(2010,"log10(theta-gamma-mu)   no cut$",70,-6d0,1d0)
      call gbook1(2011,"log10(theta-gamma-mu) with cut$",70,-6d0,1d0)
      call gbook1(2020,"log10(ene-gamma)    no cut    $",70,-5d0,2d0)
      call gbook1(2021,"log10(ene-gamma)  with cut    $",70,-5d0,2d0)
      call gbook1(2030,"cos-theta(mu-beam)  no cut    $",50,-1d0,1d0)
      call gbook1(2031,"cos-theta(mu-beam) with cut   $",50,-1d0,1d0)
      call gbook1(2040,"cos-theta(ph-beam)  no cut   $",100,-1d0,1d0)
      call gbook1(2041,"cos-theta(ph-beam) with cut  $",100,-1d0,1d0)
**********

      ELSEIF(MODE.EQ.0) THEN    
C     ======================
      IEV=IEV+1

* events outside phase space
      if(WTCRU1*WTCRU2.eq.0d0)  goto 500
*------------------------------------------------
*---------------trigger--------------------------
*------------------------------------------------
      if(nphot.ne.2)    goto 500

* photon pair mass [GeV]
      amgmin  = 1
* muon-gamma angular separation [deg]
      thmuga = 10
* muon separation from beam [cos(theta)]
      cmuon  = 0.80
* photon separation from beam [cos(theta)]
      cphot  = 0.95

      do 30 k=1,4
      ph1(k) = sphot(1,k)
      ph2(k) = sphot(2,k)
      sph(k) = ph1(k)+ph2(k)
 30   continue
* photon-photon mass
      amph2  = sph(4)**2-sph(3)**2-sph(2)**2-sph(1)**2
      amph   = sqrt(abs(amph2))
* muon-gamma angular separation [deg]
      pg1 = sqrt(ph1(1)**2+ph1(2)**2+ph1(3)**2)
      pg2 = sqrt(ph2(1)**2+ph2(2)**2+ph2(3)**2)
      pm1 = sqrt(qf1(1)**2+qf1(2)**2+qf1(3)**2)
      pm2 = sqrt(qf2(1)**2+qf2(2)**2+qf2(3)**2)
      cgmu11 = (ph1(1)*qf1(1)+ph1(2)*qf1(2)+ph1(3)*qf1(3))/pg1/pm1
      cgmu12 = (ph1(1)*qf2(1)+ph1(2)*qf2(2)+ph1(3)*qf2(3))/pg1/pm2
      cgmu21 = (ph2(1)*qf1(1)+ph2(2)*qf1(2)+ph2(3)*qf1(3))/pg2/pm1
      cgmu22 = (ph2(1)*qf2(1)+ph2(2)*qf2(2)+ph2(3)*qf2(3))/pg2/pm2
      tgmu11 = acos(min(cgmu11,0.9999999999d0))*180/pi
      tgmu12 = acos(min(cgmu12,0.9999999999d0))*180/pi
      tgmu21 = acos(min(cgmu21,0.9999999999d0))*180/pi
      tgmu22 = acos(min(cgmu22,0.9999999999d0))*180/pi
      tgmin1 = min(tgmu11,tgmu12)
      tgmin2 = min(tgmu21,tgmu22)
      tgmin  = min(tgmin1,tgmin2)
* muon separation from beam [cos(theta)]
      cmu1 = qf1(3)/pm1
      cmu2 = qf2(3)/pm2
      cmux = max(abs(cmu1),abs(cmu2))
* photon separation from beam [cos(theta)]
      cgm1 = ph1(3)/pg1
      cgm2 = ph2(3)/pg2
      cgmx = max(abs(cgm1),abs(cgm2))
************
      ene = ph2(4)
      wt  = wtmod
      call gf1(2000, amph,wt)
      call gf1(2010,dlog10(tgmin1*pi/180),wt)
      call gf1(2010,dlog10(tgmin2*pi/180),wt)
      call gf1(2020,dlog10(pg1),wt)
      call gf1(2020,dlog10(pg2),wt)
      call gf1(2030,cmu1,wt)
      call gf1(2030,cmu2,wt)
      call gf1(2040,cgm1,wt)
      call gf1(2040,cgm2,wt)
************
      if(tgmin.lt.thmuga) goto 500
*     ----------------------------

      if(cmux.gt.cmuon)   goto 500
*     ----------------------------
      if(cgmx.gt.cphot)   goto 500
*     ----------------------------
      if(amph.lt.amgmin)   goto 500
*     ----------------------------


****************
      CALL amptes(0,tgmin,amph,wtsyo,wterw)
      call gf1(2001, amph,wt)
      call gf1(2011,dlog10(tgmin1*pi/180),wt)
      call gf1(2011,dlog10(tgmin2*pi/180),wt)
      call gf1(2021,dlog10(pg1),wt)
      call gf1(2021,dlog10(pg2),wt)
      call gf1(2031,cmu1,wt)
      call gf1(2031,cmu2,wt)
      call gf1(2041,cgm1,wt)
      call gf1(2041,cgm2,wt)
****************
      call gf1(1000, ene,wt)
      call gf1(1001, ene,wt*wterw)
      call gf1(1002, amph,wt)
      call gf1(1003, amph,wt*wterw)
 500  continue
*     -----

      ELSE            
C     ===========     
      WRITE(NOUT,*) '==================================='     
      WRITE(NOUT,*) '============ ROBOL1 ==============='     
      WRITE(NOUT,*) '==================================='     
      call gprint(1000)
      call gprint(1001)
      call gprint(1002)
      call gprint(1003)
************
      call gprint(2000)
      call gprint(2001)
      call gprint(2010)
      call gprint(2011)
      call gprint(2020)
      call gprint(2021)
      call gprint(2030)
      call gprint(2031)
      call gprint(2040)
      call gprint(2041)
************
      ENDIF          
C     =====          
      END     




      SUBROUTINE ROBOL2(MODE)    
C     ***********************    
C          ===================================      
C          ============ ROBOL2 ===============      
C          ===================================     
***
* This is test of dsigma/dM and of amplitued (febr 93)
*** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      PARAMETER(PI=3.1415926535897932D0)        
      COMMON / INOUT  / NINP,NOUT           
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT     
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)      
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2,COMPFI,AMPLFI   
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /,/ WGTD2  /
      SAVE IEV,ievac,amph1
      dimension ph1(4),phac(20,4),tgmin(20)

C         
      IF(MODE.EQ.-1) THEN     
C     ===================     
C---------------------        
      WRITE(NOUT,*) '======ROBOL2 initialization==============='
      WRITE(   6,*) '======ROBOL2 initialization==============='
* booking
      call gbook1(1000,"2-nd phot. energy, wt      $", 30, 1d0,46d0)
      amph1 = 5.
      call gbook1(1002,"2-phot. mass, cut, wt  $", 30, amph1,90+amph1)
      call gbook1(1003,"2-phot. mass, cut, wt*R$", 30, amph1,90+amph1)
      call gbook1(1004,"2-phot. mass, cut, wt*Y$", 30, amph1,90+amph1)
**********
      call gbook1(2000,"ph-ph mass   no cut,       $",100, 0d0,100d0)
      call gbook1(2001,"ph-ph mass with cut,       $",100, 0d0,100d0)
      call gbook1(2010,"log10(theta-gamma-mu)   no cut$",70,-6d0,1d0)
      call gbook1(2011,"log10(theta-gamma-mu) with cut$",70,-6d0,1d0)
      call gbook1(2020,"log10(ene-gamma)    no cut    $",70,-5d0,2d0)
      call gbook1(2021,"log10(ene-gamma)  with cut    $",70,-5d0,2d0)
      call gbook1(2030,"cos-theta(mu-beam)  no cut    $",50,-1d0,1d0)
      call gbook1(2031,"cos-theta(mu-beam) with cut   $",50,-1d0,1d0)
      call gbook1(2040,"cos-theta(ph-beam)  no cut   $",100,-1d0,1d0)
      call gbook1(2041,"cos-theta(ph-beam) with cut  $",100,-1d0,1d0)
**********
      call amptes(-1,dum1,dum2,dum3,dum4)
      IEV  =0
      IEVac=0

      ELSEIF(MODE.EQ.0) THEN    
C     ======================
      IEV=IEV+1

      wt  = wtmod
* events outside phase space
      if(WTCRU1*WTCRU2.eq.0d0)  goto 500
*------------------------------------------------
*---------------trigger--------------------------
*------------------------------------------------
* photon minimum energy [GeV]
      ephmin = 1
* photon pair minimum mass [GeV]
      amgmin  = 0
* muon-gamma angular separation [deg]
      thmuga = 10
* muon separation from beam [cos(theta)]
      cmuon  = 0.80
* photon separation from beam [cos(theta)]
      cphot  = 0.95

* multiplicity and list of photons separated from muons
* and with minimum energy
      nphac=0
      pm1 = sqrt(qf1(1)**2+qf1(2)**2+qf1(3)**2)
      pm2 = sqrt(qf2(1)**2+qf2(2)**2+qf2(3)**2)
*-------------
      do 100 IPH=1,nphot
        call gf1(2020,dlog10(sphot(iph,4)),wt)
        if(sphot(iph,4).lt.ephmin) goto 100
        do  30 k=1,4
        ph1(k) = sphot(iph,k)
 30     continue
        pg1 = sqrt(ph1(1)**2+ph1(2)**2+ph1(3)**2)
* muon-gamma angular separation [deg]
        cgmu11 = (ph1(1)*qf1(1)+ph1(2)*qf1(2)+ph1(3)*qf1(3))/pg1/pm1
        cgmu12 = (ph1(1)*qf2(1)+ph1(2)*qf2(2)+ph1(3)*qf2(3))/pg1/pm2
        tgmu11 = acos(min(cgmu11,0.9999999999d0))*180/pi
        tgmu12 = acos(min(cgmu12,0.9999999999d0))*180/pi
        tgmin1 = min(tgmu11,tgmu12)
        call gf1(2010,dlog10(tgmin1*pi/180),wt)
        cgm1 = ph1(3)/pg1
        call gf1(2040,cgm1,wt)
* photon with one of angles below thmuga not accepted
        if(tgmin1.lt.thmuga) goto 100
*       -----------------------------
* photon separation from beam [cos(theta)]
        if(cgm1.gt.cphot)   goto 100
*       ----------------------------
        nphac=nphac+1
        do 40 k=1,4
          phac(nphac,k)=sphot(iph,k)
          tgmin(nphac) =tgmin1 
 40     continue
        call gf1(2041,cgm1,wt)
 100  continue
*----------------
* muon separation from beam [cos(theta)]
      cmu1 = qf1(3)/pm1
      cmu2 = qf2(3)/pm2
      cmux = max(abs(cmu1),abs(cmu2))
      call gf1(2030,cmu1,wt)
      call gf1(2030,cmu2,wt)

*  only events with 2 photons in the trigger accepted!!!!
      if(nphac.ne.2)   goto 500
*     -------------------------
* photon-photon mass
      amph2  = (phac(1,4)+phac(2,4))**2-(phac(1,3)+phac(2,3))**2
     $        -(phac(1,2)+phac(2,2))**2-(phac(1,1)+phac(2,1))**2
      amph   = sqrt(abs(amph2))
      call gf1(2000,amph,wt)

      if(amph.lt.amgmin)   goto 500
*     ----------------------------

      ievac=ievac+1
      if(ievac.lt.10) then
      write(6,*) "--------------- iev,ievac",iev,ievac     
      CALL DUMPS(6)
      write(6,"(4f15.9)") tgmu11,tgmu12,tgmin1
      write(6,*) "   amph=",amph
      write(6,*) "phac(1,4),phac(2,4)=",phac(1,4),phac(2,4)
      endif
**************** mass and angle distributions
      call gf1(1000, phac(2,4),wt)
      call gf1(2001, amph,wt)
      call gf1(2011, dlog10(tgmin(1)*pi/180),wt)
      call gf1(2011, dlog10(tgmin(2)*pi/180),wt)
      call gf1(2021, dlog10(phac(1,4)),wt)
      call gf1(2021, dlog10(phac(2,4)),wt)
**************** amlitude tests for exactly two photons
      if(nphot.ne.2)    goto 500
      if(amph.lt.amph1) goto 500
      tgminx = min(tgmin(1),tgmin(2))
      call amptes(0,tgminx,amph,wtsyo,wterw)
      call gf1(1002, amph,wt)
      call gf1(1003, amph,wt*wterw)
      call gf1(1004, amph,wt*wtsyo)
 500  continue
*     -----

      ELSE            
C     ===========     
      WRITE(NOUT,*) '==================================='     
      WRITE(NOUT,*) '============ ROBOL3 ==============='     
      WRITE(NOUT,*) '==================================='     
      CALL AMPLIT( 1)
      call gprint(1000)
      call gprint(1002)
      call gprint(1003)
      call gprint(1004)
************
      call gprint(2000)
      call gprint(2001)
      call gprint(2010)
      call gprint(2011)
      call gprint(2020)
      call gprint(2021)
      call gprint(2030)
      call gprint(2031)
      call gprint(2040)
      call gprint(2041)
************
      ENDIF          
C     =====          
 
      END     

      SUBROUTINE ROBOL3(MODE)    
C     ***********************    
C          ===================================      
C          ============ ROBOL3 ===============      
C          ===================================   
***
* Test of amplitude with L3 cut, amplitude by Ela Richter-Was 
*** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      PARAMETER(PI=3.1415926535897932D0)        
      COMMON / INOUT  / NINP,NOUT           
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT     
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)      
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2,COMPFI,AMPLFI   
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /,/ WGTD2  /
      SAVE IEV,ievac
      dimension ph1(4),ph2(4),sph(4)
C         
      IF(MODE.EQ.-1) THEN     
C     ===================     
C---------------------        
      WRITE(NOUT,*) '======ROBOL3 initialization==============='
      WRITE(   6,*) '======ROBOL3 initialization==============='
* booking
      call gbook1(1000,"2-nd phot. energy, wt      $", 30, 1d0,46d0)
      call gbook1(1001,"2-nd phot. energy, wt*R    $", 30, 1d0,46d0)
      call gbook1(1002,"2-phot. mass, cut, wt      $", 30, 5d0,95d0)
      call gbook1(1003,"2-phot. mass, cut, wt*R    $", 30, 5d0,95d0)
**********
      call gbook1(2000,"ph-ph mass   no cut,       $",100, 0d0,100d0)
      call gbook1(2001,"ph-ph mass with cut,       $",100, 0d0,100d0)
      call gbook1(2010,"log10(theta-gamma-mu)   no cut$",70,-6d0,1d0)
      call gbook1(2011,"log10(theta-gamma-mu) with cut$",70,-6d0,1d0)
      call gbook1(2020,"log10(ene-gamma)    no cut    $",70,-5d0,2d0)
      call gbook1(2021,"log10(ene-gamma)  with cut    $",70,-5d0,2d0)
      call gbook1(2030,"cos-theta(mu-beam)  no cut    $",50,-1d0,1d0)
      call gbook1(2031,"cos-theta(mu-beam) with cut   $",50,-1d0,1d0)
      call gbook1(2040,"cos-theta(ph-beam)  no cut   $",100,-1d0,1d0)
      call gbook1(2041,"cos-theta(ph-beam) with cut  $",100,-1d0,1d0)
**********
      CALL AMPLIT(-1)           
      IEV  =0
      IEVac=0

      ELSEIF(MODE.EQ.0) THEN    
C     ======================
      IEV=IEV+1

* events outside phase space
      if(WTCRU1*WTCRU2.eq.0d0)  goto 500
*------------------------------------------------
*---------------trigger--------------------------
*------------------------------------------------
      if(nphot.ne.2)    goto 500
* photon pair mass [GeV]
      amgmin  = 1
* muon-gamma angular separation [deg]
      thmuga = 10
* muon separation from beam [cos(theta)]
      cmuon  = 0.80
* photon separation from beam [cos(theta)]
      cphot  = 0.95

      do 30 k=1,4
      ph1(k) = sphot(1,k)
      ph2(k) = sphot(2,k)
      sph(k) = ph1(k)+ph2(k)
 30   continue
* photon-photon mass
      amph2  = sph(4)**2-sph(3)**2-sph(2)**2-sph(1)**2
      amph   = sqrt(abs(amph2))
* muon-gamma angular separation [deg]
      pg1 = sqrt(ph1(1)**2+ph1(2)**2+ph1(3)**2)
      pg2 = sqrt(ph2(1)**2+ph2(2)**2+ph2(3)**2)
      pm1 = sqrt(qf1(1)**2+qf1(2)**2+qf1(3)**2)
      pm2 = sqrt(qf2(1)**2+qf2(2)**2+qf2(3)**2)
      cgmu11 = (ph1(1)*qf1(1)+ph1(2)*qf1(2)+ph1(3)*qf1(3))/pg1/pm1
      cgmu12 = (ph1(1)*qf2(1)+ph1(2)*qf2(2)+ph1(3)*qf2(3))/pg1/pm2
      cgmu21 = (ph2(1)*qf1(1)+ph2(2)*qf1(2)+ph2(3)*qf1(3))/pg2/pm1
      cgmu22 = (ph2(1)*qf2(1)+ph2(2)*qf2(2)+ph2(3)*qf2(3))/pg2/pm2
      tgmu11 = acos(min(cgmu11,0.9999999999d0))*180/pi
      tgmu12 = acos(min(cgmu12,0.9999999999d0))*180/pi
      tgmu21 = acos(min(cgmu21,0.9999999999d0))*180/pi
      tgmu22 = acos(min(cgmu22,0.9999999999d0))*180/pi
      tgmin1 = min(tgmu11,tgmu12)
      tgmin2 = min(tgmu21,tgmu22)
      tgmin  = min(tgmin1,tgmin2)
* muon separation from beam [cos(theta)]
      cmu1 = qf1(3)/pm1
      cmu2 = qf2(3)/pm2
      cmux = max(abs(cmu1),abs(cmu2))
* photon separation from beam [cos(theta)]
      cgm1 = ph1(3)/pg1
      cgm2 = ph2(3)/pg2
      cgmx = max(abs(cgm1),abs(cgm2))
************
      ene = ph2(4)
      wt  = wtmod
      call gf1(2000, amph,wt)
      call gf1(2010,dlog10(tgmin1*pi/180),wt)
      call gf1(2010,dlog10(tgmin2*pi/180),wt)
      call gf1(2020,dlog10(pg1),wt)
      call gf1(2020,dlog10(pg2),wt)
      call gf1(2030,cmu1,wt)
      call gf1(2030,cmu2,wt)
      call gf1(2040,cgm1,wt)
      call gf1(2040,cgm2,wt)
************
      if(amph.lt.amgmin)   goto 500
*     ----------------------------
      if(tgmin.lt.thmuga) goto 500
*     ----------------------------
      if(cmux.gt.cmuon)   goto 500
*     ----------------------------
      if(cgmx.gt.cphot)   goto 500
*     ----------------------------

      CALL AMPLIT( 0)
      ievac=ievac+1
      if(ievac.lt.10) then
      write(6,*) "--------------- iev,ievac",iev,ievac     
c      CALL DUMPS(6)
      write(6,"(4f15.9)") tgmu11,tgmu12,tgmu21,tgmu22
      write(6,*) "tgmin====>",tgmin, "   amph=",amph
      write(6,*) "wt2=",wt2
      write(6,*) "ph1(4),ph2(4)=",ph1(4),ph2(4)
      endif
****************
      call gf1(2001, amph,wt)
      call gf1(2011,dlog10(tgmin1*pi/180),wt)
      call gf1(2011,dlog10(tgmin2*pi/180),wt)
      call gf1(2021,dlog10(pg1),wt)
      call gf1(2021,dlog10(pg2),wt)
      call gf1(2031,cmu1,wt)
      call gf1(2031,cmu2,wt)
      call gf1(2041,cgm1,wt)
      call gf1(2041,cgm2,wt)
****************
      call gf1(1000, ene,wt)
      call gf1(1001, ene,wt*wt2)
      call gf1(1002, amph,wt)
      call gf1(1003, amph,wt*wt2)
 500  continue
*     -----

      ELSE            
C     ===========     
      WRITE(NOUT,*) '==================================='     
      WRITE(NOUT,*) '============ ROBOL3 ==============='     
      WRITE(NOUT,*) '==================================='     
      CALL AMPLIT( 1)
      call gprint(1000)
      call gprint(1001)
      call gprint(1002)
      call gprint(1003)
************
      call gprint(2000)
      call gprint(2001)
      call gprint(2010)
      call gprint(2011)
      call gprint(2020)
      call gprint(2021)
      call gprint(2030)
      call gprint(2031)
      call gprint(2040)
      call gprint(2041)
************
      ENDIF          
C     =====          
      END    


      SUBROUTINE ROBOL4(MODE)           
C     ***********************           
C          ===================================    
C          ============ ROBOL4 ===============    
C          ===================================    
C study on COMPACT/AMPLITUD  D2 matrix element,2 photons  old one  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON / INOUT  / NINP,NOUT           
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT     
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)      
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2,COMPFI,AMPLFI   
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /,/ WGTD2  /
      DIMENSION  PT(100),PC(100),IND(100)         
      SAVE NBX,XMIN,XMAX
C         
      IF(MODE.EQ.-1) THEN     
C     ===================     
C---------------------        
      WRITE(NOUT,*) '======ROBOL4 initialization==============='
      WRITE(   6,*) '======ROBOL4 initialization==============='
      NBX =  50       
      XMIN=  0D0      
      XMAX=  1D0      
      DO 13 J=1,2     
      CALL GBOOK1( 170+J,'compact/amplit VERS E/EEAM   170+J PHOTON $'
     $  ,NBX,XMIN,XMAX)         
      CALL GBOOK1(1170+J,'count on beam entry E/EEAM   170+J PHOTON $'
     $  ,NBX,XMIN,XMAX)         
      CALL GBOOK1( 180+J,'compact/amplit VERS P_T/EEAM 180+J PHOTON $'  
     $  ,NBX,XMIN,XMAX)         
      CALL GBOOK1(1180+J,'count on beam entry P_T/EEAM 180+J PHOTON $'
     $  ,NBX,XMIN,XMAX)         
      CALL GBOOK1( 190+J,'compact/amplit VERS P_C/E_ph 190+J PHOTON $'  
     $  ,NBX,XMIN,XMAX)         
  13  CALL GBOOK1(1190+J,'count on beam entry P_C/E_ph 190+J PHOTON $'
     $  ,NBX,XMIN,XMAX)         
      CALL AMPLIT(-1)           
      ELSEIF(MODE.EQ.0) THEN    
C     ======================    
      CALL AMPLIT( 0)           
C.....only 2 hard photons accepted        
      IF(WTHARD.EQ.0D0) GOTO 53           
      WT=WTMOD        
C-----------------------------------------------------------------------
            
C     ordering pC/E_ph in respect to final states   
      CALL ORDPCF(PC,IND)       
C     ordering pT in respect to final states        
      CALL ORDPTF(PT,IND)       
      DO 18 J=1,2     
      IF(NPHOT.GE.J) THEN       
C hardest photons     
       ZE=      2*SPHOT(J,4)/CMSENE       
       CALL GF1( 170+J,ZE,WT2)            
       CALL GF1(1170+J,ZE,1D0)            
C photons with highest pT       
       ZPT=      2*PT(J)/CMSENE           
       CALL GF1( 180+J,ZPT, WT2)          
       CALL GF1(1180+J,ZPT,1D0)           
C photons with highest pC colinearity     
       ZPC=        PC(J)        
       CALL GF1( 190+J,ZPC, WT2)          
       CALL GF1(1190+J,ZPC,1D0)           
      ENDIF           
   18 CONTINUE        
C---------------------------------------------------------------------- 
C calculating ratio of compact to amplit (from spin amplitudes)         
C version of function D2        
            
  53  CONTINUE        
      ELSE            
C     ===========     
      WRITE(NOUT,*) '==================================='     
      WRITE(NOUT,*) '============ ROBOL4 ==============='     
      WRITE(NOUT,*) '==================================='     
      CALL AMPLIT( 1)           
C     ********************************    
C           
      DO 312 J=1,2              
      DO 312 K=0,1              
        CALL GPRINT (170+1000*K+J)      
        CALL GPRINT (180+1000*K+J)      
        CALL GPRINT (190+1000*K+J)  
 312  CONTINUE

      ENDIF          
C     =====          
      END            

      SUBROUTINE ORDPTF(PT,IND)          
C     ************************           
C ordering photons according to pT, PT is list of ordered pT 
C and IND is the adress list of photons (indices in SPHOT)   
C pT is calculated in respect to final fermions    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
      DIMENSION PT(*),IND(*)             
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT      
      SAVE   / MOMSET /
      DIMENSION VERS(3)        
           
      DO 5 I=1,3     
   5  VERS(I)=QF1(I)-QF2(I)    
      DO 10 I=1,NPHOT          
      PL   =(SPHOT(I,1)*VERS(1)+SPHOT(I,2)*VERS(2)+SPHOT(I,3)*VERS(3)) 
     %     /DSQRT(VERS(1)**2+VERS(2)**2+VERS(3)**2)          
      PT(I)=SQRT(SPHOT(I,4)**2-PL**2)    
   10 IND(I)=I       
      IF(NPHOT.LE.1) RETURN    
      DO 30 I=2,NPHOT          
      DO 30 J=NPHOT,I,-1       
      IF(PT(J).GT.PT(J-1)) THEN          
        X=PT(J)      
        PT(J)=PT(J-1)          
        PT(J-1)=X              
        L=IND(J)     
        IND(J)=IND(J-1)        
        IND(J-1)=L             
      ENDIF          
   30 CONTINUE       
      END            

      SUBROUTINE ORDPCF(PC,IND)          
C     ************************           
C ordering photons according to pC, PC is list of ordered pC 
C and IND is the adress list of photons (indices in SPHOT)   
C pC collinearity is calculated in respect to final fermions 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
      DIMENSION PC(*),IND(*)             
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT      
      SAVE   / MOMSET /
      DIMENSION VERS(3)        
           
      DO 5 I=1,3     
   5  VERS(I)=QF1(I)-QF2(I)    
      DO 10 I=1,NPHOT          
      PL   =(SPHOT(I,1)*VERS(1)+SPHOT(I,2)*VERS(2)+SPHOT(I,3)*VERS(3)) 
     %     /DSQRT(VERS(1)**2+VERS(2)**2+VERS(3)**2)          
      PC(I)=PL/SPHOT(I,4)      
   10 IND(I)=I       
      IF(NPHOT.LE.1) RETURN    
      DO 30 I=2,NPHOT          
      DO 30 J=NPHOT,I,-1       
      IF(PC(J).GT.PC(J-1)) THEN          
        X=PC(J)      
        PC(J)=PC(J-1)          
        PC(J-1)=X              
        L=IND(J)     
        IND(J)=IND(J-1)        
        IND(J-1)=L             
      ENDIF          
   30 CONTINUE       
      END        !\end



      SUBROUTINE amptes(mode,tgmin,amph,wtsyo,wterw)
C     ***********************
C          ===================================
C          =========== amptes  ===============
C          ===================================
***
* New test of amplitude with L3 cut, 
* Amplitude by Scott Yost and Ela Richter-Was
***
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.1415926535897932D0)
      COMMON / INOUT  / NINP,NOUT
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT
      COMMON / INPARM / CMSENE,AMFIN,KEYBRM,KEYZET
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)
! commons to erichter -------
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2,COMPFI,AMPLFI   
      SAVE   / WGTD2  /
! commons to syost ----------
      COMMON /trmkey/ DCSkey, LLkey, Zkey
      integer         DCSkey, LLkey, Zkey
      SAVE / INOUT  /,/ MOMSET /,/ INPARM  /,/ WGTALL /
      SAVE /trmkey/
! --------
      SAVE IEV,ievac,p1,p2
      DIMENSION ph1(4),ph2(4),sph(4),p1(4),p2(4)
      DATA iev,ievac /0,0/
 
      IF (mode.eq.-1) THEN
!     ====================
      DCSkey = 100
      LLkey  = 100
      Zkey   = KEYZET
      write(   6,*) " Amptes: Zkey   = KEYZET = ",Zkey
      write(nout,*) " Amptes: Zkey   = KEYZET = ",Zkey
 
      CALL AMPLIT(-1)
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
      iev=iev+1
      if(nphot.ne.2) RETURN
 
      DO 30 k=1,4
      ph1(k) = sphot(1,k)
      ph2(k) = sphot(2,k)
      sph(k) = ph1(k)+ph2(k)
 30   CONTINUE

! Richter-Was
      CALL amplit( 0)
      wterw = AMPLFI/COMPFI
! Yost
      CALL twopho(dcs, dcsll, dcss, p1, qf1, p2, qf2, ph1, ph2)
C.. cross section compact formula from YFS3 version january 1990
      CALL dfinap2(cmsene,    p1, qf1, p2, qf2, ph1, ph2, compf2)
      wtsyo = DCS/COMPF2
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Controll printouts
      sumeng  = sph(4)
      sum_xk  = 2*sumeng/cmsene
      xk_min  = min(sphot(1,4),sphot(2,4))/cmsene
      if(amph.lt. 50.d0 ) goto 800
ccc      if(sumeng.lt.0.2) goto 800
ccc      if(xk_min.lt.0.1) goto 800
ccc      if(sum_xk.gt.0.00001) goto 800
ccc      if(sum_xk.lt.0.001. or. sum_xk.gt.0.01) RETURN
      IF(ievac.gt.30) goto 800
      ievac=ievac+1
      write( 6,*) "--------------- iev,ievac",iev,ievac
      write(16,*) "--------------- iev,ievac",iev,ievac
ccc      CALL dumps( 6)
      CALL dumps(16)
      write( 6,"(a,23f15.10)") "sum_xk tgmin,amph=",sum_xk,tgmin,amph
      write(16,"(a,23f15.10)") "sum_xk tgmin,amph=",sum_xk,tgmin,amph
      write( 6,*) "ytes: DCS,      DCS/COMPF2=",DCS,      DCS/COMPF2
      write(16,*) "ytes: DCS,      DCS/COMPF2=",DCS,      DCS/COMPF2
c      write( 6,*) "ytes: DCSll  ,DCSll/COMPF2=",DCSll,  DCSll/COMPF2
c      write(16,*) "ytes: DCSll  ,DCSll/COMPF2=",DCSll,  DCSll/COMPF2
      write( 6,*) "ytes: AMPLFI,AMPLFI/COMPFI=",AMPLFI,AMPLFI/COMPFI
      write(16,*) "ytes: AMPLFI,AMPLFI/COMPFI=",AMPLFI,AMPLFI/COMPFI
      write( 6,*)  "ytes:        COMPF2/COMPFI=",COMPF2/COMPFI
      write(16,*)  "ytes:        COMPF2/COMPFI=",COMPF2/COMPFI
 800  continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      ENDIF
 
      END

      SUBROUTINE dfinap2(cmsene,p1,p2,q1,q2,pk1,pk2,xection)
!     ******************************************************
! Compact form for double final bremstr cross section from yfs3 M.C.
! Note BHLUMI notation for  p1,q1,p2,q2.
C     ***************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DCSkey, LLkey, Zkey
      PARAMETER(pi=3.1415926535897932d0, alfinv=137.03604d0)
      DIMENSION p1(4),q1(4),p2(4),q2(4),pk1(4),pk2(4),xx(4)
      COMMON /trmkey/ DCSkey, LLkey, Zkey
      SAVE /trmkey/
 
      alfa=1d0/alfinv
      S0=cmsene**2
      S1=(p2(4)+q2(4))**2-(p2(1)+q2(1))**2-(p2(2)+q2(2))**2
     $   -(p2(3)+q2(3))**2
      DO 10 i=1,4
      xx(i)=p1(i)+q1(i)
  10  CONTINUE
      CALL sfdist2(xx,p1,q1,p2,q2,pk1,pk2,dist2)
      xection= (S0/S1) * dist2
      END
