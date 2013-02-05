
      SUBROUTINE BOK3PH(MODE,XPAR,NPAR)
C     *********************************
C histograming distributions from TOPIK1
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /GENER2/ X1,X2,WTMOD
      COMMON /WEIGHTS/ WTINF, WTYFS, WTEXA
!      COMMON /MOMLAB2/ GLU1L(4),GLU2L(4),TOP1L(4),TOP2L(4),PHOT1L(4)
!     #                ,PHOT2L(4)
      COMMON /MOMLAB3/ 
     # GLU1L(4),GLU2L(4),P2L(4),Q2L(4),PHOT1L(4),PHOT2L(4),PHOT3L(4)

      DIMENSION XPAR(99),NPAR(99)
      DIMENSION PHX(4),PHY(4),PHZ(4)
      LOGICAL IVIS1,IVIS2,IVIS3
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)


      IF(MODE.EQ.-1) THEN
C     =======================
      IDENT = 7000
      NBIN  = 50
      CMSENE= XPAR(1)
      AMTOP2= XPAR(2)**2
      AMI2  = 4D0*AMTOP2/CMSENE**2
c histograms initialization for gamma energy
      TMAXA = 1D0-AMI2
      TMINA = 0.D0
       CALL GBOOK1(IDENT+85,'EPH  pierwszy photon $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+85,'ERRO')
       CALL GBOOK1(IDENT+86,'EPH  drugi photon    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+86,'ERRO')
       CALL GBOOK1(IDENT+87,'EPH  trzec photon    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+87,'ERRO')
c histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+91,' dsigma/dcos  2 photon $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+91,'ERRO')
      CALL GBOOK1(IDENT+97,'dsigma/dcos pierw photon$',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+97,'ERRO')
      CALL GBOOK1(IDENT+98,'dsigma/dcos drug photon $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+98,'ERRO')
      CALL GBOOK1(IDENT+99,'dsigma/dcos trzc photon $',NBIN,TMINB,TMAXB)
      CALL GIDOPT(IDENT+99,'ERRO')
c histograms initialization for gamma-gamma mass
      TMAXE =   200D0
      TMINE =   0.0D0
      CALL GBOOK1(IDENT+40,'gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+40,'ERRO')
      CALL GBOOK1(IDENT+41,'recoil   1  mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+41,'ERRO')
      CALL GBOOK1(IDENT+42,'recoil   2  mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+42,'ERRO')
c delta histograms initialization for gamma energy
      TMAXA = 1D0-AMI2
      TMINA = 0.D0
       CALL GBOOK1(IDENT+185,'d EPH  pierwszy ph $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+185,'ERRO')
       CALL GBOOK1(IDENT+186,'d EPH drugi ph    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+186,'ERRO')
c delta histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+191,'d dsig/dcos  2 ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+191,'ERRO')
      CALL GBOOK1(IDENT+197,'d dsig/dcos pierw ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+197,'ERRO')
      CALL GBOOK1(IDENT+198,'d dsig/dcos drug ph $',NBIN,TMINB,TMAXB) 
      CALL GIDOPT(IDENT+198,'ERRO')
c histograms initialization for gamma-gamma mass
      TMAXE =   200D0
      TMINE =   0.0D0
      CALL GBOOK1(IDENT+140,'d gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+140,'ERRO')

      CALL GBOOK1(IDENT+141,'d recoil  1   mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+141,'ERRO')

      CALL GBOOK1(IDENT+142,'d recoil  2   mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+142,'ERRO')

c delta histograms initialization for gamma energy
      TMAXA = 1D0-AMI2
      TMINA = 0.D0
       CALL GBOOK1(IDENT+385,'D EPH  pierwszy ph $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+385,'ERRO')
       CALL GBOOK1(IDENT+386,'D EPH drugi ph    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+386,'ERRO')
c delta histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+391,'D dsig/dcos  2 ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+391,'ERRO')
      CALL GBOOK1(IDENT+397,'D dsig/dcos pierw ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+397,'ERRO')
      CALL GBOOK1(IDENT+398,'D dsig/dcos drug ph $',NBIN,TMINB,TMAXB) 
      CALL GIDOPT(IDENT+398,'ERRO')
c histograms initialization for gamma-gamma mass
      TMAXE =   200D0
      TMINE =   0.0D0
      CALL GBOOK1(IDENT+340,'D gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+340,'ERRO')

C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
       EMINVIS=0.5
       ETRIG=3.0
       CVIS=  0.97  !!! 0.7
      DO K=1,4
       PHX(k)=phot1l(k)
       PHy(k)=phot2l(k)
       PHz(k)=phot3l(k)
      ENDDO
       NVIS=0
       IVIS1=((PHX(4).GT.EMINVIS).and.(abs(phx(3)/phx(4)).lt.cvis))
       IVIS2=((PHy(4).GT.EMINVIS).and.(abs(phy(3)/phy(4)).lt.cvis))
       IVIS3=((PHz(4).GT.EMINVIS).and.(abs(phz(3)/phz(4)).lt.cvis))
       IF (IVIS1) NVIS=NVIS+1
       IF (IVIS2) NVIS=NVIS+1
       IF (IVIS3) NVIS=NVIS+1


       DO k=1,4
       IF     (NVIS.eq.0) THEN
!        do nothing
       ELSEIF (NVIS.eq.1) THEN
        IF     (ivis1) then
!        do nothing
        elseif (ivis2) then
         phot1l(k)=phy(k)
         phot2l(k)=phx(k)
        else
         phot1l(k)=phz(k)
         phot3l(k)=phx(k)
        endif       
       ELSEIF (NVIS.eq.2) THEN
        IF     (ivis1.and.ivis2.and.phx(4).gt.phy(4)) then  
!        do nothing 
        elseIF (ivis1.and.ivis2.and.phx(4).lt.phy(4)) then  
         phot1l(k)=phy(k)
         phot2l(k)=phx(k)
         phot3l(k)=phx(k)
        elseif (ivis1.and.ivis3.and.phx(4).gt.phz(4)) then  
         phot1l(k)=phx(k)
         phot2l(k)=phz(k)
         phot3l(k)=phy(k)
        elseif (ivis1.and.ivis3.and.phx(4).lt.phz(4)) then  
         phot1l(k)=phz(k)
         phot2l(k)=phx(k)
         phot3l(k)=phy(k)
        elseif (ivis2.and.ivis3.and.phy(4).gt.phz(4)) then  
         phot1l(k)=phy(k)
         phot2l(k)=phz(k)
         phot3l(k)=phx(k)
        elseif (ivis2.and.ivis3.and.phy(4).lt.phz(4)) then  
         phot1l(k)=phz(k)
         phot2l(k)=phy(k)
         phot3l(k)=phx(k)
        endif
       ELSEIF (NVIS.eq.3) THEN
        if     (phx(4).gt.phy(4).and.phy(4).gt.phz(4)) then  
! do nothing
        elseif (phx(4).gt.phz(4).and.phz(4).gt.phy(4)) then  
         phot1l(k)=phx(k)
         phot2l(k)=phz(k)
         phot3l(k)=phy(k)
        elseif (phy(4).gt.phx(4).and.phx(4).gt.phz(4)) then  
         phot1l(k)=phy(k)
         phot2l(k)=phx(k)
         phot3l(k)=phz(k)
        elseif (phy(4).gt.phz(4).and.phz(4).gt.phx(4)) then  
         phot1l(k)=phy(k)
         phot2l(k)=phz(k)
         phot3l(k)=phx(k)
        elseif (phz(4).gt.phy(4).and.phy(4).gt.phx(4)) then  
         phot1l(k)=phz(k)
         phot2l(k)=phy(k)
         phot3l(k)=phx(k)
        elseif (phz(4).gt.phx(4).and.phx(4).gt.phy(4)) then  
         phot1l(k)=phz(k)
         phot2l(k)=phx(k)
         phot3l(k)=phy(k)
        endif

       ENDIF
       enddo


      ENE1   =  PHOT1L(4)/(CMSENE/2D0)
      COSTH1 =(PHOT1L(3)*GLU1L(3)+PHOT1L(2)*GLU1L(2)+PHOT1L(1)*GLU1L(1))
     %        /DSQRT(PHOT1L(3)**2+PHOT1L(2)**2+PHOT1L(1)**2)
     %        /DSQRT(GLU1L(3)**2+GLU1L(2)**2+GLU1L(1)**2)
       pt1=dsqrt(PHOT1L(2)**2+PHOT1L(1)**2)

       ENE2  = PHOT2L(4)/(CMSENE/2D0)
      COSTH2 =(PHOT2L(3)*GLU1L(3)+PHOT2L(2)*GLU1L(2)+PHOT2L(1)*GLU1L(1))
     %        /DSQRT(PHOT2L(3)**2+PHOT2L(2)**2+PHOT2L(1)**2)
     %        /DSQRT(GLU1L(3)**2+GLU1L(2)**2+GLU1L(1)**2)
       pt2=DSQRT(PHOT2L(2)**2+PHOT2L(1)**2)         
      XMASG2=((PHOT1L(4)+PHOT2L(4))**2-(PHOT1L(3)+PHOT2L(3))**2
     #            -(PHOT1L(2)+PHOT2L(2))**2-(PHOT1L(1)+PHOT2L(1))**2)
       ENE3  = PHOT3L(4)/(CMSENE/2D0)
      COSTH3 =(PHOT3L(3)*GLU1L(3)+PHOT3L(2)*GLU1L(2)+PHOT3L(1)*GLU1L(1))
     %        /DSQRT(PHOT3L(3)**2+PHOT3L(2)**2+PHOT3L(1)**2)
     %        /DSQRT(GLU1L(3)**2+GLU1L(2)**2+GLU1L(1)**2)
      IF(XMASG2.LE.0D0) THEN
        XMASGG=0D0
      ELSE
        XMASGG=DSQRT(XMASG2)
      ENDIF
      XMASR2=((PHOT3L(4)+P2L(4)+Q2L(4))**2
     #       -(PHOT3L(3)+P2L(3)+Q2L(3))**2
     #       -(PHOT3L(2)+P2L(2)+Q2L(2))**2
     #       -(PHOT3L(1)+P2L(1)+Q2L(1))**2)
      IF(XMASR2.LE.0D0) THEN
        XMASR2=0D0
      ELSE
        XMASR2=DSQRT(XMASR2)
      ENDIF
      XMASR1=((PHOT3L(4)+PHOT2L(4)+P2L(4)+Q2L(4))**2
     #       -(PHOT3L(3)+PHOT2L(3)+P2L(3)+Q2L(3))**2
     #       -(PHOT3L(2)+PHOT2L(2)+P2L(2)+Q2L(2))**2
     #       -(PHOT3L(1)+PHOT2L(1)+P2L(1)+Q2L(1))**2)
  
      IF(XMASR1.LE.0D0) THEN
        XMASR1=0D0
      ELSE
        XMASR1=DSQRT(XMASR1)
      ENDIF

            CALL GF1(IDENT+40,XMASGG,WTMOD)
            IF (ENE1.gt.etrig/(CMSENE/2D0)) then      
             if (NVIS.eq.1) CALL GF1(IDENT+41,XMASR1,WTMOD)      
             if (NVIS.eq.2) CALL GF1(IDENT+42,XMASR2,WTMOD)
            endif      
            CALL GF1(IDENT+85,ENE1,WTMOD)      
            CALL GF1(IDENT+86,ENE2,WTMOD) 
            CALL GF1(IDENT+87,ENE3,WTMOD) 
     
            CALL GF1(IDENT+91,COSTH1,WTMOD)      
            CALL GF1(IDENT+97,COSTH1,WTMOD)      

            CALL GF1(IDENT+91,COSTH2,WTMOD)      
            CALL GF1(IDENT+98,COSTH2,WTMOD)  
    
            CALL GF1(IDENT+99,COSTH3,WTMOD)      
            
            wtx=abs(wtyfs-wtexa)
            CALL GF1(IDENT+140,XMASGG,WTX)  
            IF (ENE1.gt.etrig/(CMSENE/2D0)) then    
             if (NVIS.eq.1) CALL GF1(IDENT+141,XMASR1,WTX)      
             if (NVIS.eq.2) CALL GF1(IDENT+142,XMASR2,WTX)      
            endif      
            CALL GF1(IDENT+185,ENE1,WTX)      
            CALL GF1(IDENT+186,ENE2,WTX) 
     
            CALL GF1(IDENT+191,COSTH1,WTX)      
            CALL GF1(IDENT+197,COSTH1,WTX)      


            CALL GF1(IDENT+191,COSTH2,WTX)      
            CALL GF1(IDENT+198,COSTH2,WTX)      
            wtx=(wtyfs-wtexa)
            CALL GF1(IDENT+340,XMASGG,WTX)      
            CALL GF1(IDENT+385,ENE1,WTX)      
            CALL GF1(IDENT+386,ENE2,WTX) 
     
            CALL GF1(IDENT+391,COSTH1,WTX)      
            CALL GF1(IDENT+397,COSTH1,WTX)      


            CALL GF1(IDENT+391,COSTH2,WTX)      
            CALL GF1(IDENT+398,COSTH2,WTX)      
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      XSECT =XPAR(10)
c...for each event two photons stored
      IEVENT=NPAR(10)*2
      FACT = NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 85,FACT,NBIN)
      CALL GPRINT(IDENT+85)
      CALL RENHIE(IDENT+ 86,FACT,NBIN)
      CALL GPRINT(IDENT+86)
      CALL RENHIE(IDENT+ 87,FACT,NBIN)
      CALL GPRINT(IDENT+87)
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 91,FACT,NBIN)
      CALL GPRINT(IDENT+91)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 97,FACT,NBIN)
      CALL GPRINT(IDENT+97)
      CALL RENHIE(IDENT+ 98,FACT,NBIN)
      CALL GPRINT(IDENT+98)
      CALL RENHIE(IDENT+ 99,FACT,NBIN)
      CALL GPRINT(IDENT+99)
      FACT = NBIN*XSECT/(TMAXe-TMINe)/NPAR(10)
      CALL RENHIE(IDENT+ 40,FACT,NBIN)
      CALL GPRINT(IDENT+40)
      CALL RENHIE(IDENT+ 41,FACT,NBIN)
      CALL GPRINT(IDENT+41)
      CALL RENHIE(IDENT+ 42,FACT,NBIN)
      CALL GPRINT(IDENT+42)
      IEVENT=NPAR(10)*2
      FACT = NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 185,FACT,NBIN)
!!!      CALL GPRINT(IDENT+185)
      CALL RENHIE(IDENT+ 186,FACT,NBIN)
!!!      CALL GPRINT(IDENT+186)
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 191,FACT,NBIN)
!!!      CALL GPRINT(IDENT+191)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 197,FACT,NBIN)
!!!      CALL GPRINT(IDENT+197)
      CALL RENHIE(IDENT+ 198,FACT,NBIN)
!!!      CALL GPRINT(IDENT+198)
      FACT = NBIN*XSECT/(TMAXe-TMINe)/NPAR(10)
      CALL RENHIE(IDENT+ 140,FACT,NBIN)
      CALL RENHIE(IDENT+ 141,FACT,NBIN)
      CALL RENHIE(IDENT+ 142,FACT,NBIN)
C
      IEVENT=NPAR(10)*2
      FACT = NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 385,FACT,NBIN)
!!!      CALL GPRINT(IDENT+385)
      CALL RENHIE(IDENT+ 386,FACT,NBIN)
!!!      CALL GPRINT(IDENT+386)
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 391,FACT,NBIN)
!!!      CALL GPRINT(IDENT+391)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 397,FACT,NBIN)
!!!      CALL GPRINT(IDENT+397)
      CALL RENHIE(IDENT+ 398,FACT,NBIN)
!!!      CALL GPRINT(IDENT+398)
      FACT = NBIN*XSECT/(TMAXe-TMINe)/NPAR(10)
      CALL RENHIE(IDENT+ 340,FACT,NBIN)

C
      IA=85
      CALL GBOOK1(ident+200+IA,'EPH  pierwszy ph $', NBIN,tmina,tmaxa)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=86
      CALL GBOOK1(ident+200+IA,'EPH  drugi ph    $', NBIN,tmina,tmaxa)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=97
      CALL GBOOK1(ident+200+IA,'cos pierw ph    $', NBIN,tminb,tmaxb)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=98
      CALL GBOOK1(ident+200+IA,'cos drug ph    $', NBIN,tminb,tmaxb)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=40
      CALL GBOOK1(ident+200+IA,'M gamma-gamma  $', NBIN,tmine,tmaxe)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=41
      CALL GBOOK1(ident+200+IA,'M recoil   1   $', NBIN,tmine,tmaxe)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
      IA=42
      CALL GBOOK1(ident+200+IA,'M recoil   2  $', NBIN,tmine,tmaxe)
      CALL GOPERA(ident+100+IA,'/',ident+IA,ident+200+IA,1.D0,1.D0)
      CALL GMINIM(ident+200+IA,-1.D0)
      CALL GMAXIM(ident+200+IA, 1.D0)
      CALL GPRINT(ident+200+IA)
C
      IA=85
      CALL GBOOK1(ident+400+IA,'EPH  pierwszy ph $', NBIN,tmina,tmaxa)
      CALL GOPERA(ident+300+IA,'/',ident+IA,ident+400+IA,1.D0,1.D0)
      CALL GMINIM(ident+400+IA,-1.D0)
      CALL GMAXIM(ident+400+IA, 1.D0)
      CALL GPRINT(ident+400+IA)
      IA=86
      CALL GBOOK1(ident+400+IA,'EPH  drugi ph    $', NBIN,tmina,tmaxa)
      CALL GOPERA(ident+300+IA,'/',ident+IA,ident+400+IA,1.D0,1.D0)
      CALL GMINIM(ident+400+IA,-1.D0)
      CALL GMAXIM(ident+400+IA, 1.D0)
      CALL GPRINT(ident+400+IA)
      IA=97
      CALL GBOOK1(ident+400+IA,'cos pierw ph    $', NBIN,tminb,tmaxb)
      CALL GOPERA(ident+300+IA,'/',ident+IA,ident+400+IA,1.D0,1.D0)
      CALL GMINIM(ident+400+IA,-1.D0)
      CALL GMAXIM(ident+400+IA, 1.D0)
      CALL GPRINT(ident+400+IA)
      IA=98
      CALL GBOOK1(ident+400+IA,'cos drug ph    $', NBIN,tminb,tmaxb)
      CALL GOPERA(ident+300+IA,'/',ident+IA,ident+400+IA,1.D0,1.D0)
      CALL GMINIM(ident+400+IA,-1.D0)
      CALL GMAXIM(ident+400+IA, 1.D0)
      CALL GPRINT(ident+400+IA)
      IA=40
      CALL GBOOK1(ident+400+IA,'M gamma-gamma  $', NBIN,tmine,tmaxe)
      CALL GOPERA(ident+300+IA,'/',ident+IA,ident+400+IA,1.D0,1.D0)
      CALL GMINIM(ident+400+IA,-1.D0)
      CALL GMAXIM(ident+400+IA, 1.D0)
      CALL GPRINT(ident+400+IA)

      ENDIF

      END

