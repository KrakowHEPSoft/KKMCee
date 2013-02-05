C WARNING: FOR KEYDIS=410-422 DELB MODIFIED &&&&               
C WARNING NUMBERING FOR YFS KEYDIS CHANGED TO 301,302,303 
C------------------------------------------------------   
C Convention for KEYDIS      
C     KEYDIS   =  0 + R*100  Zero   Order     
C     KEYDIS   =  1 + R*100  First  Order     
C     KEYDIS   =  2 + R*100  Second Order     
C     KEYDIS   = 10 + R*100  Beta0         Zero   Order        
C     KEYDIS   = 11 + R*100  Beta0         First  Order        
C     KEYDIS   = 12 + R*100  Beta1            
C     KEYDIS   = 20 + R*100  Beta0         Second Order        
C     KEYDIS   = 21 + R*100  Beta1            
C     KEYDIS   = 22 + R*100  Beta2                Order        
C     R = 100 Kuraev Fadin   
C     R = 300 YFS            
C     R = 400 YFS single str. funct.          
C------------------------------------------------------   
      PROGRAM MAIN
C     ***********************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(20000)
      COMMON / INOUT  / NINP,NOUT  
     
      CALL GLIMIT(20000)                 
      NINP= 5           
      NOUT= 6           
      OPEN(NOUT,file="tbfig.output")
      CALL GOUTPU(NOUT)

      CALL SUSFIG          
C     ***********          
      CLOSE(NOUT)          
      END                  

      SUBROUTINE SUSFIG    
C     *****************    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT     
      COMMON / BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF              
      COMMON / DNAME / DNAME                
      CHARACTER*7 DNAME    
      SAVE   / INOUT  /,/ BEAMPM /,/ DNAME /

      WRITE(NOUT,*) '==============================================' 
      WRITE(NOUT,*) '=============S U S F I G======================' 
      WRITE(NOUT,*) '==============================================' 
                  
 3000 FORMAT(8I2)      
 3001 FORMAT(I10)   
 3002 FORMAT(F10.0)                  
 3010 FORMAT(A7)
C BEAM MASS (ELECTRON)               
      AMIN = 0.0005111               
      READ( NINP,3001)  NEVT
      READ( NINP,3010)  DNAME
      WRITE(   6,3010)  DNAME
      READ( NINP,3000)  KAT1,KAT2,KAT3,KAT4,KAT5,KAT6,KAT7,KAT8
      READ(NINP,3001) KEYBRM         
      READ(NINP,3002) CMSENE,SINW2,AMAZ,GAMMZ         
      READ(NINP,3002) AMFIN,EPS,VVMAX                 
      READ(NINP,3001) KEYRED,KEYWGT  
      ENE=CMSENE/2  
   
C ----restore histograms             
      NINPH=10      
      OPEN(NINPH,file="../" //dname// "/" //dname// ".hst")
      CALL GRFILE(NINPH,DNAME,' ')   
      CALL GRIN(   0,9999,0)         
C=====================================                
C ----initialize GPLOT               
      CALL GPLINT(0)                 
      NOUFIG=11
      OPEN(NOUFIG,file="./" //dname// ".tex")
      CALL GPLCAP(-NOUFIG)
      CALL GPLTIT(' This is SUSFIG   $')              
C=================================   
C Study of x-sections with L3 cut
      IF(KAT2.NE.0) CALL FIG2                                           
C D_2 compact / Exact for L3 cut
      IF(KAT3.NE.0) CALL FIG3                                           
C comparison of D_2 compact and spin amplit formulas for two 
C hard photons                                               
      IF(KAT4.NE.0) CALL FIG4                                           
C=================================              
      CALL GPLEND              
      END     

      SUBROUTINE FIG2        
C     ********************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
! This plots dsig/dM2 [nb] and tests of amplitudes
! for two hard acollinear photons

c============================================================
c============================================================
      CALL GPLTIT('dsigma/dE2  [nb], real cuts $')  
      CALL GPRINT(1000)
      CALL RENHST("NB  ",1000,1070)
      CALL GIDOPT(1070,'ERRO')         
      CALL GMAXIM(1070,  1.5D-3)            
      CALL GMINIM(1070,  0.0D0)            
      CALL GPRINT(1070)
**************
      CALL GPLSET('DMOD',2D0)            
      CALL GPLOT(1070,' ','*',0)       
c============================================================
      CALL GPLTIT('dsig/dMgg [nb], th(mu,g)min=10deg, cos(mu)max=0.8,
     $ cos(g)max=0.8 $')  
      CALL GPRINT(1002)
      CALL RENHST("NB  ",1002,1072)
      CALL GIDOPT(1072,'ERRO')         
      CALL GMAXIM(1072,  3.00D-4)            
      CALL GMINIM(1072,  0.0D0)            
      CALL GPRINT(1072)
**************
      CALL GPLSET('DMOD',2D0)            
      CALL GPLOT(1072,' ','*',0)
* larger scale
      CALL GMAXIM(1072,  4.0D-6)            
      CALL GMINIM(1072,  0.0D0)            
      CALL GPLSET('DMOD',2D0)            
      CALL GPLOT(1072,' ','*',0)       
c============================================================
      CALL GPLTIT(' /exact-amplit1/D2-compact., real cut $')  
      CALL GPRINT(1002)
      CALL GPRINT(1003)
      CALL GOPERA (1003,'/',1002,1054,1D0,1D0)    
      CALL GIDOPT( 1054,'ERRO')         
      CALL GMAXIM( 1054,  1.30D0)            
      CALL GMINIM( 1054,  0.50D0)            
      CALL GPRINT(1054)
**************
      CALL GPLSET('DMOD',7D0)            
      CALL GPLOT(1054,' ','*',0)       
c============================================================
      CALL GPLTIT(' exact-amplit2/D2-compact., real cut $')  
      CALL GPRINT(1002)
      CALL GPRINT(1004)
      CALL GOPERA (1004,'/',1002,1055,1D0,1D0)    
      CALL GIDOPT( 1055,'ERRO')         
      CALL GMAXIM( 1055,  1.30D0)            
      CALL GMINIM( 1055,  0.50D0)            
      CALL GPRINT(1055)
**************
      CALL GPLSET('DMOD',7D0)            
      CALL GPLOT(1055,' ','*',0)       

      END
      
      SUBROUTINE FIG3        
C     ********************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 

c============================================================
c============================================================
c      CALL GPLTIT(' D2-compact./exact-amplit, real cut $')  
c      CALL GPRINT(1000)
c      CALL GPRINT(1001)
c      CALL GOPERA (1001,'/',1000,1052,1D0,1D0)    
c      CALL GIDOPT( 1052,'ERRO')         
c      CALL GMAXIM( 1052,  1.30D0)            
c      CALL GMINIM( 1052,  0.50D0)            
c      CALL GPRINT(1052)
**************
c      CALL GPLSET('DMOD',7D0)            
c      CALL GPLOT(1052,' ','*',0)       
c============================================================
c============================================================
      CALL GPLTIT(' D2-compact./exact-amplit, real cut $')  
      CALL GPRINT(1002)
      CALL GPRINT(1003)
      CALL GOPERA (1003,'/',1002,1054,1D0,1D0)    
      CALL GIDOPT( 1054,'ERRO')         
      CALL GMAXIM( 1054,  1.30D0)            
      CALL GMINIM( 1054,  0.50D0)            
      CALL GPRINT(1054)
**************
      CALL GPLSET('DMOD',7D0)            
      CALL GPLOT(1054,' ','*',0)       
c============================================================
c============================================================
c      CALL GPLTIT('dsigma/dE2  [nb], real cuts $')  
c      CALL GPRINT(1000)
c      CALL RENHST("NB  ",1000,1070)
c      CALL GIDOPT(1070,'ERRO')         
c      CALL GMAXIM(1070,  2.5D-6)            
c      CALL GMINIM(1070,  0.0D0)            
c      CALL GPRINT(1070)
**************
c      CALL GPLSET('DMOD',2D0)            
c      CALL GPLOT(1070,' ','*',0)       
c============================================================
      CALL GPLTIT('dsigma/dM  [nb], real cuts $')  
      CALL GPRINT(1002)
      CALL RENHST("NB  ",1002,1072)
      CALL GIDOPT(1072,'ERRO')         
      CALL GMAXIM(1072,  2.50D-4)            
      CALL GMINIM(1072,  0.0D0)            
      CALL GPRINT(1072)
**************
      CALL GPLSET('DMOD',2D0)            
      CALL GPLOT(1072,' ','*',0)
* larger scale
      CALL GMAXIM(1072,  2.5D-6)            
      CALL GMINIM(1072,  0.0D0)            
      CALL GPLSET('DMOD',2D0)            
      CALL GPLOT(1072,' ','*',0)       

      END

      SUBROUTINE RENHST(CHAK,ID1,ID2)     
C     *******************************
*! This routine normalizes to x-section in NANOBARNS or to unity
*  CHAK = "NB  "    normal case [nb]
*  CHAK = "NB10"    log10 x-scale assumed [nb]
*  CHAK = "UNIT"    normalization to unity
*! id2.ne.id1 required
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      DIMENSION  XMEMO(100)
C
      IF(ID2.EQ.ID1) GOTO 900
      IDM = 1
      CALL GUNPAK(IDM,XMEMO,'  ',1)
      NEVT    = XMEMO( 1)
      XSCRNB  = XMEMO(20)   
      CALL GINBO1(ID1,TITLE,NBT,TMIN,TMAX)
      FACT=1D0
      IF( CHAK.EQ. "NB  ") THEN
        FACT = NBT*XSCRNB/(NEVT*(TMAX-TMIN))
      ELSEIF( CHAK.EQ. "NB10") THEN
        FLN10 = LOG(10.)
        FACT = NBT*XSCRNB/(NEVT*(TMAX-TMIN)*FLN10)
      ELSEIF( CHAK.EQ. "UNIT") THEN
        WTAVR  = XMEMO(10)/XMEMO(20)
        FACT = NBT/(NEVT*(TMAX-TMIN))/WTAVR 
      ELSE
        WRITE(6,*) '+++++ RENHST: wrong chak=',chak
      ENDIF 
C Multiply content
      CALL GOPERA(ID1,'+',ID1,ID2, FACT, 0D0)
      RETURN
 900  WRITE(6,*) '+++++ RENHST: ID1=ID2=',ID1
      END  

      SUBROUTINE FIG4        
C     ********************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C----------------------------------------------------------------------C
C D2 FUNCTION comparison of compact and spin amplitude formula         C
C in function of p_T/E_beam and p_T/E_beam                             C
C----------------------------------------------------------------------C
        
      DO 312 J=1,2              
      CALL GOPERA (170+J,'/',1170+J,170+J,1D0,1D0)    
      CALL GOPERA (180+J,'/',1180+J,180+J,1D0,1D0)   
      CALL GOPERA (190+J,'/',1190+J,190+J,1D0,1D0)   
      CALL GIDOPT( 170+J,'ERRO')         
      CALL GPRINT (170+J)      
      CALL GIDOPT( 180+J,'ERRO')         
      CALL GPRINT (180+J)      
      CALL GIDOPT( 190+J,'ERRO')         
  312 CALL GPRINT (190+J)  
    
c.....ratio of compact/spin amplit function of E/E_beam         
      CALL GMAXIM( 0,  1.30D0)            
      CALL GMINIM( 0,  0.50D0)            
c.....first photon      
      CALL GPLTIT(' FIG4, D2 comp./amplit vers  E/Ebeam firs phot  $')  
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(170+1,' ','*',0)       
C     CALL GPLOT(170+1,' ',' ',0)       
c.....second photon             
      CALL GPLTIT(' FIG4, D2 comp./amplit vers  E/Ebeam second phot  $')
        
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(170+2,' ','*',0)       
C     CALL GPLOT(170+2,' ',' ',0)       
c.....ratio of compact/spin amplit function of p_T/E_beam       
      CALL GMAXIM( 0,  1.30D0)            
      CALL GMINIM( 0,  0.50D0)            
c.....first photon      
      CALL GPLTIT(' FIG3, D2 comp./amplit vers pT/Ebeam first phot $')  
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(180+1,' ','*',0)       
C     CALL GPLOT(180+1,' ',' ',0)       
c.....second photon             
      CALL GPLTIT(' FIG3, D2 comp./amplit vers pT/Ebeam second phot $') 
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(180+2,' ','*',0)       
C     CALL GPLOT(180+2,' ',' ',0)       
c.....ratio of compact/spin amplit function of p_C/E_phot collinearity  
      CALL GMAXIM( 0,  1.30D0)            
      CALL GMINIM( 0,  0.50D0)            
c.....first photon      
      CALL GPLTIT(' FIG3, D2 comp./amplit vers pC/Ephot first phot $')  
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(190+1,' ','*',0)       
C     CALL GPLOT(190+1,' ',' ',0)       
c.....second photon             
      CALL GPLTIT(' FIG3, D2 comp./amplit vers pC/Ephot second phot $') 
      CALL GPLSET('DMOD',1D0)            
      CALL GPLOT(190+2,' ','*',0)       
C     CALL GPLOT(190+2,' ',' ',0)       
c..............         
      call gprint(171)          
      call gprint(172)          
      call gprint(181)          
      call gprint(182)          
      call gprint(191)          
      call gprint(192)     
      END
