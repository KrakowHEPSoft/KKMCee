      SUBROUTINE TRZYNUE(MODE,XPAR,NPAR)
C     *************************************
C histograming distributions from TOPIK1
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /WEIGHTS/ WTINF, WTYFS, WTEXA
      COMMON /GENER2/ X1,X2,WTMOD1
      COMMON /MOMLAB2/ GLU1L(4),GLU2L(4),TOP1L(4),TOP2L(4),PHOT1L(4)
     #                ,PHOT2L(4)
      DIMENSION XPAR(99),NPAR(99)

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
      TMAXA = 100.D0
      TMINA =   0.D0
      CALL GBOOK1(IDENT+81,'EPH  pierwszy photon $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+81,'ERRO')
      CALL GBOOK1(IDENT+82,'EPH  drugi photon    $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+82,'ERRO')
      CALL GBOOK1(IDENT+86,'PT  pierwszy photon $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+86,'ERRO')
      CALL GBOOK1(IDENT+87,'PT  drugi photon    $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+87,'ERRO')
c histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+83,'dsigma/dcos pierw photon$',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+83,'ERRO')
      CALL GBOOK1(IDENT+84,'dsigma/dcos drug photon $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+84,'ERRO')
c histograms initialization for recoling mass
      TMAXC = 200.D0
      TMINC =   0.D0
      CALL GBOOK1(IDENT+85,'recolling mass $',NBIN,TMINC,TMAXC)       
      CALL GIDOPT(IDENT+85,'ERRO')
c histograms initialization for gamma energy
      TMAXA = 100.D0
      TMINA =   0.D0
      CALL GBOOK1(IDENT+71,'EPH  pierwszy photon YFS$',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+71,'ERRO')
      CALL GBOOK1(IDENT+72,'EPH  drugi photon  YFS  $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+72,'ERRO')
      CALL GBOOK1(IDENT+76,'PT  pierwszy photon YFS$',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+76,'ERRO')
      CALL GBOOK1(IDENT+77,'PT  drugi photon  YFS  $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+77,'ERRO')
c histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+73,'dsigma/dcos pierw photon YFS$',
     #           NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+73,'ERRO')
      CALL GBOOK1(IDENT+74,'dsigma/dcos drug photon YFS $',
     #           NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+74,'ERRO')
c histograms initialization for recoling mass
      TMAXC = 200.D0
      TMINC =   0.D0
      CALL GBOOK1(IDENT+75,'recolling mass YFS  $',NBIN,TMINC,TMAXC)       
      CALL GIDOPT(IDENT+75,'ERRO')

C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
      ENE1=PHOT1L(4)
      COSTH1 = (PHOT1L(3)*GLU1L(3)+PHOT1L(2)*GLU1L(2)
     %        +PHOT1L(1)*GLU1L(1))
     %        /DSQRT(PHOT1L(3)**2+PHOT1L(2)**2+PHOT1L(1)**2)
     %        /DSQRT(GLU1L(3)**2+GLU1L(2)**2+GLU1L(1)**2)
      PTPHOT1=SQRT(PHOT1L(1)**2+PHOT1L(2)**2)
      ENE2=PHOT2L(4)
      COSTH2 = (PHOT2L(3)*GLU1L(3)+PHOT2L(2)*GLU1L(2)
     %        +PHOT2L(1)*GLU1L(1))
     %        /DSQRT(PHOT2L(3)**2+PHOT2L(2)**2+PHOT2L(1)**2)
     %        /DSQRT(GLU1L(3)**2+GLU1L(2)**2+GLU1L(1)**2)
      PTPHOT2=SQRT(PHOT2L(1)**2+PHOT2L(2)**2)
      IF(ENE1.GT.ENE2) THEN
        IF(ENE1.GT.0.5.AND.ABS(COSTH1).LT.0.96) THEN
            CALL GF1(IDENT+81,ENE1,WTMOD1)      
            CALL GF1(IDENT+83,COSTH1,WTMOD1)      
            CALL GF1(IDENT+86,PTPHOT1,WTMOD1)      
            CALL GF1(IDENT+71,ENE1,WTYFS)      
            CALL GF1(IDENT+73,COSTH1,WTYFS)      
            CALL GF1(IDENT+76,PTPHOT1,WTYFS)      
            IF(ENE2.GT.0.5.AND.ABS(COSTH2).LT.0.96) THEN
               CALL GF1(IDENT+82,ENE2,WTMOD1)      
               CALL GF1(IDENT+84,COSTH2,WTMOD1)
               CALL GF1(IDENT+87,PTPHOT2,WTMOD1)      
               CALL GF1(IDENT+72,ENE2,WTYFS)      
               CALL GF1(IDENT+74,COSTH2,WTYFS)
               CALL GF1(IDENT+77,PTPHOT2,WTYFS)      
            ENDIF
        ENDIF
      ELSE      
        IF(ENE2.GT.0.5.AND.ABS(COSTH2).LT.0.96) THEN
            CALL GF1(IDENT+81,ENE2,WTMOD1)      
            CALL GF1(IDENT+83,COSTH2,WTMOD1)      
            CALL GF1(IDENT+87,PTPHOT2,WTMOD1)      
            CALL GF1(IDENT+71,ENE2,WTYFS)      
            CALL GF1(IDENT+73,COSTH2,WTYFS)      
            CALL GF1(IDENT+77,PTPHOT2,WTYFS)      
            IF(ENE1.GT.0.5.AND.ABS(COSTH1).LT.0.96) THEN
               CALL GF1(IDENT+82,ENE1,WTMOD1)      
               CALL GF1(IDENT+84,COSTH1,WTMOD1)
               CALL GF1(IDENT+86,PTPHOT1,WTMOD1)      
               CALL GF1(IDENT+72,ENE1,WTYFS)      
               CALL GF1(IDENT+74,COSTH1,WTYFS)
               CALL GF1(IDENT+76,PTPHOT1,WTYFS)      
            ENDIF
        ENDIF
      ENDIF
      PXREC=TOP1L(1)+TOP2L(1)
      PYREC=TOP1L(2)+TOP2L(2)
      PZREC=TOP1L(3)+TOP2L(3)
      EEREC=TOP1L(4)+TOP2L(4)
      IF(ENE1.GT.0.5.AND.ABS(COSTH1).LT.0.96) THEN
        PXREC=PXREC+PHOT1L(1)
        PYREC=PYREC+PHOT1L(2)
        PZREC=PZREC+PHOT1L(3)
        EEREC=EEREC+PHOT1L(4)
        IF(ENE2.GT.0.5.AND.ABS(COSTH2).LT.0.96) THEN
          PXREC=PXREC+PHOT2L(1)
          PYREC=PYREC+PHOT2L(2)
          PZREC=PZREC+PHOT2L(3)
          EEREC=EEREC+PHOT2L(4)
        ENDIF            
        XMREC=EEREC**2-PXREC**2-PYREC**2-PZREC**2
        IF(XMREC.GT.0.001) XMREC=SQRT(XMREC)
        CALL GF1(IDENT+85,XMREC,WTMOD1)      
        CALL GF1(IDENT+75,XMREC,WTYFS)      
      ENDIF
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      XSECT =XPAR(10)
c...................
      IEVENT=NPAR(10)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      CALL RENHIE(IDENT+ 81,FACT,NBIN)
      CALL GPRINT(IDENT+81)
      CALL RENHIE(IDENT+ 82,FACT,NBIN)
      CALL GPRINT(IDENT+82)
      CALL RENHIE(IDENT+ 86,FACT,NBIN)
      CALL GPRINT(IDENT+86)
      CALL RENHIE(IDENT+ 87,FACT,NBIN)
      CALL GPRINT(IDENT+87)
      CALL RENHIE(IDENT+ 71,FACT,NBIN)
      CALL GPRINT(IDENT+71)
      CALL RENHIE(IDENT+ 72,FACT,NBIN)
      CALL GPRINT(IDENT+72)
      CALL RENHIE(IDENT+ 76,FACT,NBIN)
      CALL GPRINT(IDENT+76)
      CALL RENHIE(IDENT+ 77,FACT,NBIN)
      CALL GPRINT(IDENT+77)
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 83,FACT,NBIN)
      CALL GPRINT(IDENT+83)
      CALL RENHIE(IDENT+ 84,FACT,NBIN)
      CALL GPRINT(IDENT+84)
      CALL RENHIE(IDENT+ 73,FACT,NBIN)
      CALL GPRINT(IDENT+73)
      CALL RENHIE(IDENT+ 74,FACT,NBIN)
      CALL GPRINT(IDENT+74)
      FACT = NBIN*XSECT/(TMAXC-TMINC)/IEVENT
      CALL RENHIE(IDENT+ 85,FACT,NBIN)
      CALL GPRINT(IDENT+85)
      CALL RENHIE(IDENT+ 75,FACT,NBIN)
      CALL GPRINT(IDENT+75)
C
      ENDIF

      END


      SUBROUTINE BOK2PH(MODE,XPAR,NPAR)
C     *********************************
C histograming distributions from TOPIK1
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /GENER2/ X1,X2,WTMOD
      COMMON /WEIGHTS/ WTINF, WTYFS, WTEXA
      COMMON /MOMLAB2/ GLU1L(4),GLU2L(4),TOP1L(4),TOP2L(4),PHOT1L(4)
     #                ,PHOT2L(4)
      DIMENSION XPAR(99),NPAR(99)

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
c histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+91,' dsigma/dcos  2 photon $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+91,'ERRO')
      CALL GBOOK1(IDENT+97,'dsigma/dcos pierw photon$',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+97,'ERRO')
      CALL GBOOK1(IDENT+98,'dsigma/dcos drug photon $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+98,'ERRO')
c histograms initialization for gamma-gamma mass
      TMAXE =   200D0
      TMINE =   0.0D0
      CALL GBOOK1(IDENT+40,'gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+40,'ERRO')
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
      IF(XMASG2.LE.0D0) THEN
        XMASGG=0D0
      ELSE
        XMASGG=DSQRT(XMASG2)
      ENDIF
  
           if (ene2.gt.ene1) then
            xx=ene2
            ene2=ene1
            ene1=xx
            xx=COSTH2
            COSTH2=COSTH1
            COSTH1=xx
            xx=pt2
            pt2=pt1
            pt1=xx
           endif
           if (pt1.gt.3.and.pt2.gt.3) then
            CALL GF1(IDENT+40,XMASGG,WTMOD)      
            CALL GF1(IDENT+85,ENE1,WTMOD)      
            CALL GF1(IDENT+86,ENE2,WTMOD) 
     
            CALL GF1(IDENT+91,COSTH1,WTMOD)      
            CALL GF1(IDENT+97,COSTH1,WTMOD)      

            CALL GF1(IDENT+91,COSTH2,WTMOD)      
            CALL GF1(IDENT+98,COSTH2,WTMOD)      
            
            wtx=abs(wtyfs-wtexa)
            CALL GF1(IDENT+140,XMASGG,WTX)      
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
           endif
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
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 91,FACT,NBIN)
      CALL GPRINT(IDENT+91)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 97,FACT,NBIN)
      CALL GPRINT(IDENT+97)
      CALL RENHIE(IDENT+ 98,FACT,NBIN)
      CALL GPRINT(IDENT+98)
      FACT = NBIN*XSECT/(TMAXe-TMINe)/NPAR(10)
      CALL RENHIE(IDENT+ 40,FACT,NBIN)
      CALL GPRINT(IDENT+40)
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

      SUBROUTINE BOK3PH(MODE,XPAR,NPAR)
C     *********************************
C histograming distributions from TOPIK1
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /GENER2/ X1,X2,WTMODorig
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
c acolinearities of photons
      TMAXac =   180D0
      TMINac =   0.0D0

      CALL GBOOK1(IDENT+30,'gamma-gamma acopl ME3 $',NBIN,TMINac,TMAXac)
      CALL GIDOPT(IDENT+30,'ERRO')
      CALL GBOOK1(IDENT+31,'gamma-gamma acopl ME3 $',4   ,TMINac,TMAXac)       
      CALL GIDOPT(IDENT+31,'ERRO')
      CALL GBOOK1(IDENT+32,'gamma-gamma acopl YFS $',NBIN,TMINac,TMAXac)       
      CALL GIDOPT(IDENT+32,'ERRO')
      CALL GBOOK1(IDENT+33,'gamma-gamma acopl YFS $',4   ,TMINac,TMAXac)       
      CALL GIDOPT(IDENT+33,'ERRO')
      CALL GBOOK1(IDENT+34,'gamma-gamma acopl ME2 $',NBIN,TMINac,TMAXac)       
      CALL GIDOPT(IDENT+34,'ERRO')
      CALL GBOOK1(IDENT+35,'gamma-gamma acopl ME2 $',4   ,TMINac,TMAXac)       
      CALL GIDOPT(IDENT+35,'ERRO')

c delta histograms initialization for gamma energy
      TMAXA = 1D0-AMI2
      TMINA = 0.D0
       CALL GBOOK1(IDENT+185,'d EPH  pierwszy ph $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+185,'ERRO')
       CALL GBOOK1(IDENT+186,'d EPH drugi ph    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+186,'ERRO')
       CALL GBOOK1(IDENT+187,'d EPH trzec ph    $',NBIN,TMINA,TMAXA)       
       CALL GIDOPT(IDENT+187,'ERRO')
c delta histograms initialization for gamma angular distribution
      TMAXB =  1.D0
      TMINB = -1.D0
      CALL GBOOK1(IDENT+191,'d dsig/dcos  2 ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+191,'ERRO')
      CALL GBOOK1(IDENT+197,'d dsig/dcos pierw ph $',NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+197,'ERRO')
      CALL GBOOK1(IDENT+198,'d dsig/dcos drug ph $',NBIN,TMINB,TMAXB) 
      CALL GIDOPT(IDENT+198,'ERRO')
      CALL GBOOK1(IDENT+199,'d dsig/dcos trze ph $',NBIN,TMINB,TMAXB) 
      CALL GIDOPT(IDENT+199,'ERRO')
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
      uru=0
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

       goto 155
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
 155   continue
        pt1=phot1l(1)**2+phot1l(2)**2
        pt2=phot2l(1)**2+phot2l(2)**2
        pt3=phot3l(1)**2+phot3l(2)**2
        wtmod=wtmodorig
        if (pt3.gt.pt2) wtmod=0
        if (pt3.gt.pt1) wtmod=0

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
            xacopl=xacol(phot1l,phot2l,2)
!         if (xacopl.lt.3.6.and.wtmod.gt.0d0) uru=(uru+wtmod)
!         if (xacopl.lt.3.6.and.wtmod.gt.0d0) write(*,*) '##',wtmod,uru 
            if (wtmod.gt.10.d-9) then
            write(*,*) '##',wtmod,uru
            write(*,*) phot1l
            write(*,*) phot2l
            write(*,*) phot3l
            write(*,*) ' '
            endif


            if (wtmod.gt.0d0) CALL GF1(IDENT+30,Xacopl,WTMOD)
            if (wtmod.gt.0d0) CALL GF1(IDENT+31,Xacopl,WTMOD)
            if (wtyfs.gt.0d0) CALL GF1(IDENT+32,Xacopl,WTyfs)      
            if (wtyfs.gt.0d0) CALL GF1(IDENT+33,Xacopl,WTyfs)      
            if (wtinf.gt.0d0) CALL GF1(IDENT+34,Xacopl,WTINF)
            if (wtinf.gt.0d0) CALL GF1(IDENT+35,Xacopl,WTINF)
     
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
            CALL GF1(IDENT+187,ENE3,WTX) 
     
            CALL GF1(IDENT+191,COSTH1,WTX)      
            CALL GF1(IDENT+197,COSTH1,WTX)      


            CALL GF1(IDENT+191,COSTH2,WTX)      
            CALL GF1(IDENT+198,COSTH2,WTX)      

            CALL GF1(IDENT+199,COSTH2,WTX)      
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
      FACT = NBIN*XSECT/(TMAXac-TMINac)/NPAR(10)
      CALL GMINIM(ident+30,0.D0)
      CALL GMINIM(ident+31,0.D0)
      CALL GMINIM(ident+32,0.D0)
      CALL GMINIM(ident+33,0.D0)
      CALL GMINIM(ident+34,0.D0)
      CALL GMINIM(ident+35,0.D0)
!      CALL RENHIE(IDENT+ 30,FACT,NBIN)
      CALL GPRINT(IDENT+30)
!      CALL RENHIE(IDENT+ 31,FACT,NBIN)
      CALL GPRINT(IDENT+31)
!      CALL RENHIE(IDENT+ 32,FACT,NBIN)
      CALL GPRINT(IDENT+32)
!      CALL RENHIE(IDENT+ 33,FACT,NBIN)
      CALL GPRINT(IDENT+33)
!      CALL RENHIE(IDENT+ 34,FACT,NBIN)
      CALL GPRINT(IDENT+34)
!      CALL RENHIE(IDENT+ 35,FACT,NBIN)
      CALL GPRINT(IDENT+35)
      IEVENT=NPAR(10)*2
      FACT = NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 185,FACT,NBIN)
!!!      CALL GPRINT(IDENT+185)
      CALL RENHIE(IDENT+ 186,FACT,NBIN)
!!!      CALL GPRINT(IDENT+186)
      CALL RENHIE(IDENT+ 187,FACT,NBIN)
!!!      CALL GPRINT(IDENT+187)
      FACT = NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+ 191,FACT,NBIN)
!!!      CALL GPRINT(IDENT+191)
      FACT = NBIN*XSECT/(TMAXA-TMINA)/NPAR(10)
      CALL RENHIE(IDENT+ 197,FACT,NBIN)
!!!      CALL GPRINT(IDENT+197)
      CALL RENHIE(IDENT+ 198,FACT,NBIN)
!!!      CALL GPRINT(IDENT+198)
      CALL RENHIE(IDENT+ 199,FACT,NBIN)
!!!      CALL GPRINT(IDENT+199)
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
      IA=87
      CALL GBOOK1(ident+200+IA,'EPH  trzec ph    $', NBIN,tmina,tmaxa)
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
      IA=99
      CALL GBOOK1(ident+200+IA,'cos trzc ph    $', NBIN,tminb,tmaxb)
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


      SUBROUTINE BOKER7(MODE,XPAR,NPAR)
C     *********************************
C histograming distributions from TOPIK2
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /GENER2/ X1,X2,WTMOD1
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON /MOMLAB2/ GLU1L(4),GLU2L(4),TOP1L(4),TOP2L(4),PHOT1L(4)
     #                ,PHOT2L(4)
      DIMENSION XPAR(99),NPAR(99)


      IF(MODE.EQ.-1) THEN
C     =======================
      IDENT = 6000 
      NBIN  = 50
c histograms initialization for t tbar transverse momentum spectrum 
      TMAXA = 50.D0
      TMINA = 0.D0
      CALL GBOOK1(IDENT+10,'PT of t,tbar TOPIK2 $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+10,'ERRO')
      CALL GBOOK1(IDENT+110,'PT of t,tbar TOPIK2 $',NBIN,TMINA,TMAXA)       
      CALL GIDOPT(IDENT+110,'ERRO')
      TMAXD = 50.D0
      TMIND = 0.D0
      CALL GBOOK1(IDENT+20,'PT of photon TOPIK2 $',NBIN,TMIND,TMAXD)       
      CALL GIDOPT(IDENT+20,'ERRO')
      CALL GBOOK1(IDENT+120,'PT of photon TOPIK2 $',NBIN,TMIND,TMAXD)       
      CALL GIDOPT(IDENT+120,'ERRO')
      TMAXB =   3.D0
      TMINB =  -3.D0
      CALL GBOOK1(IDENT+50,'t eta TOPIK2 $', NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+50,'ERRO')
      CALL GBOOK1(IDENT+150,'t eta TOPIK2 $', NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+150,'ERRO')
      TMAXB =   3.D0
      TMINB =  -3.D0
      CALL GBOOK1(IDENT+60,'photon eta TOPIK2 $', NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+60,'ERRO')
      CALL GBOOK1(IDENT+160,'photon eta TOPIK2 $', NBIN,TMINB,TMAXB)       
      CALL GIDOPT(IDENT+160,'ERRO')
      TMAXC =   200D0
      TMINC =   0.0D0
      CALL GBOOK1(IDENT+30,'hard energy TOPIK2 $', NBIN,TMINC,TMAXC)       
      CALL GIDOPT(IDENT+30,'ERRO')
      CALL GBOOK1(IDENT+130,'hard energy TOPIK2 $', NBIN,TMINC,TMAXC)       
      CALL GIDOPT(IDENT+130,'ERRO')
      TMAXE =   200D0
      TMINE =   0.0D0
      CALL GBOOK1(IDENT+40,'gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+40,'ERRO')
      CALL GBOOK1(IDENT+140,'gamma-gamma mass $', NBIN,TMINE,TMAXE)       
      CALL GIDOPT(IDENT+140,'ERRO')
 
C
      ELSEIF(MODE.EQ.0) THEN
C     =======================
         CALL GF1(IDENT+30,YMSENE,WTMOD1)      
         CALL GF1(IDENT+130,YMSENE,WTMOD1)  
      PTTOP=DSQRT(TOP1L(1)**2+TOP1L(2)**2) 
         CALL GF1(IDENT+10,PTTOP,WTMOD1) 
         CALL GF1(IDENT+110,PTTOP,WTMOD1)
      PTGAM=DSQRT(PHOT1L(1)**2+PHOT1L(2)**2) 
         CALL GF1(IDENT+20,PTGAM,WTMOD1)
         CALL GF1(IDENT+120,PTGAM,WTMOD1) 
      PTGAM=DSQRT(PHOT2L(1)**2+PHOT2L(2)**2) 
         CALL GF1(IDENT+20,PTGAM,WTMOD1) 
         CALL GF1(IDENT+120,PTGAM,WTMOD1) 
      ETAUP = TOP1L(4)+TOP1L(3)
      ETADN = TOP1L(4)-TOP1L(3)
      IF(ABS(ETAUP).LT.1D-4.OR.ABS(ETADN).LT.1D-4) THEN
       ETOP=10.0
      ELSE  
       ETOP=0.5D0*  LOG(ABS(ETAUP/ETADN))
      ENDIF
         CALL GF1(IDENT+50,ETOP,WTMOD1)
         CALL GF1(IDENT+150,ETOP,WTMOD1)
      ETAUP = PHOT1L(4)+PHOT1L(3)
      ETADN = PHOT1L(4)-PHOT1L(3)
      IF(ABS(ETAUP).LT.1D-4.OR.ABS(ETADN).LT.1D-4) THEN
       EGAM=10.0
      ELSE  
       EGAM=0.5D0*  LOG(ABS(ETAUP/ETADN))
      ENDIF
         CALL GF1(IDENT+60,EGAM,WTMOD1)
         CALL GF1(IDENT+160,EGAM,WTMOD1)
      ETAUP = PHOT2L(4)+PHOT2L(3)
      ETADN = PHOT2L(4)-PHOT2L(3)
      IF(ABS(ETAUP).LT.1D-4.OR.ABS(ETADN).LT.1D-4) THEN
       EGAM=10.0
      ELSE  
       EGAM=0.5D0*  LOG(ABS(ETAUP/ETADN))
      ENDIF
         CALL GF1(IDENT+60,EGAM,WTMOD1)
         CALL GF1(IDENT+160,EGAM,WTMOD1)
      XMASG2=((PHOT1L(4)+PHOT2L(4))**2-(PHOT1L(3)+PHOT2L(3))**2
     #            -(PHOT1L(2)+PHOT2L(2))**2-(PHOT1L(1)+PHOT2L(1))**2)
      IF(XMASG2.LE.0D0) THEN
        XMASGG=0D0
      ELSE
        XMASGG=DSQRT(XMASG2)
      ENDIF
         CALL GF1(IDENT+40,XMASGG,WTMOD1)      
         CALL GF1(IDENT+140,XMASGG,WTMOD1) 
C
      ELSEIF(MODE.EQ.1) THEN
C     =========================
      XSECT=XPAR(10)
      IEVENT=NPAR(10)
      FACT=NBIN*XSECT/(TMAXA-TMINA)/IEVENT
      CALL RENHIE(IDENT+10,FACT,NBIN)
      FACT=NBIN*XSECT/(TMAXD-TMIND)/IEVENT/2
      CALL RENHIE(IDENT+20,FACT,NBIN)
      FACT=NBIN*XSECT/(TMAXB-TMINB)/IEVENT
      CALL RENHIE(IDENT+50,FACT,NBIN)
      FACT=NBIN*XSECT/(TMAXB-TMINB)/IEVENT/2
      CALL RENHIE(IDENT+60,FACT,NBIN)
      FACT=NBIN*XSECT/(TMAXC-TMINC)/IEVENT
      CALL RENHIE(IDENT+30,FACT,NBIN)
      FACT=NBIN*XSECT/(TMAXE-TMINE)/IEVENT
      CALL RENHIE(IDENT+40,FACT,NBIN)
      CALL GPRINT(IDENT+30)
      CALL GPRINT(IDENT+10)
      CALL GPRINT(IDENT+20)
      CALL GPRINT(IDENT+50)
      CALL GPRINT(IDENT+60)
      CALL GPRINT(IDENT+40)
      CALL CUMHI3(IDENT+110,NBIN,IEVENT,XSECT)
      CALL CUMHI3(IDENT+120,NBIN,IEVENT*2,XSECT)
      CALL CUMHI3(IDENT+150,NBIN,IEVENT,XSECT)
      CALL CUMHI3(IDENT+160,NBIN,IEVENT*2,XSECT)
      CALL CUMHI3(IDENT+130,NBIN,IEVENT,XSECT)
      CALL CUMHI3(IDENT+140,NBIN,IEVENT,XSECT)
C
      ENDIF

      END

      FUNCTION XACOL(X,Y,N)
C     ********************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8    X(*),Y(*)
      DIMENSION X1(4),Y1(4)
      DATA PI /3.1415926535897932D0/
      S=0.D0
      X2=0.D0
      Y2=0.D0
      DO 9  I=1,N
      X1(I)=X(I)
    9 Y1(I)=Y(I)
      DO 10 I=1,N
      S=S+X1(I)*Y1(I)
      X2=X2+X1(I)**2
   10 Y2=Y2+Y1(I)**2
      XACOL=ACOS(-S/SQRT(X2*Y2))*180.D0/PI
      RETURN
      END

