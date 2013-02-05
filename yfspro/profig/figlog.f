
* to be activated /weking/ and /keyyfs/ has to be repalced with /xpar/
C WARNING: FOR KEYDIS=410-422 DELB MODIFIED &&&&               
C WARNING NUMBERING FOR YFS KEYDIS CHANGED TO 301,302,303 
C------------------------------------------------------   
C Convention for KEYDIS      
C Pedagogical exercises
C     KEYDIS   =  1     soft part YFS       First  Order
C     KEYDIS   =  2     soft part YFS       Second Order
C     KEYDIS   =  5     hard non-exp.       First  Order
C     KEYDIS   =  6     hard non-exp.       Second Order
C     KEYDIS   =  9     reference distr. of YFS paper
C Total results
C     KEYDIS   =  0 + R*100                 Zero   Order     
C     KEYDIS   =  1 + R*100                 First  Order     
C     KEYDIS   =  2 + R*100                 Second Order     
C Beta contributions
C     KEYDIS   = 10 + R*100     Beta0       Zero   Order        
C     KEYDIS   = 11 + R*100     Beta0       First  Order        
C     KEYDIS   = 12 + R*100     Beta1            
C     KEYDIS   = 20 + R*100     Beta0       Second Order        
C     KEYDIS   = 21 + R*100     Beta1            
C     KEYDIS   = 22 + R*100     Beta2        
C     R = 100 Kuraev-Fadin   
C     R = 300 YFS            
C     R = 400 YFS single electron LL str. funct.          
C------------------------------------------------------   
      PROGRAM MAIN
C     ***********************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib  / b(50000)
      COMMON / INOUT  / NINP,NOUT  
       
      CALL GLIMIT(50000)                 
      NINP=  5           
      NOUT= 16           
      CALL GOUTPU(NOUT)
      CALL YFSFIG          
C     ***********          
      END                  

      SUBROUTINE YFSFIG    
C     *****************    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT     
      COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF 
      COMMON / KEYYFS / KEYZET,KEYBRM,KEYFIX,KEYRED,KEYWGT     
      COMMON / VVREC  / VVMIN,VVMAX,VV,BETI  
      SAVE / INOUT  /,/ WEKING /,/ KEYYFS /,/ VVREC  /
      COMMON /CMSENE/ CMSENE
      SAVE   /CMSENE/
C-- this commeon transfers test name from data set to ploting 
      character*5 tname
      common / tname / tname 
      SAVE   / tname /
C
      WRITE( 6,*) '==============================================' 
      WRITE( 6,*) '=============Y F S F I G======================' 
      WRITE( 6,*) '==============================================' 
C
 3000 FORMAT(A5)
 3001 FORMAT(8I2)
 3002 FORMAT(I10)      
 3003 FORMAT(F10.0)    
      READ( NINP,3000) tname
      WRITE(   6,3000) tname
      READ(NINP,3001) KAT1,KAT2,KAT3,KAT4,KAT5,KAT6,KAT7,KAT8
      READ(NINP,3002) NEVTOT
      READ(NINP,3002) KEYBRM     
      READ(NINP,3003) CMSENE,SINW2,AMAZ,GAMMZ       
      READ(NINP,3003) AMFIN,VVMIN,VVMAX   
      READ(NINP,3002) KEYRED,KEYWGT,KEYZET
C
      ENE = CMSENE/2D0
      AMEL =  0.5111D-3   
      WRITE(   6,'(8A6/8I6)')
     $ 'KAT1','KAT2','KAT3','KAT4','KAT5','KAT6','KAT7','KAT8',
     $  KAT1 , KAT2 , KAT3 , KAT4 , KAT5 , KAT6 , KAT7 , KAT8

C=====================================================
C=====================================================
C
C Leading Logarithmic exercises
C
      CALL LLGFIG
C=================================              
      END     



      SUBROUTINE LLGFIG
C     *****************
c  ****************************************
c  *         l        ooo      ggg        *
c  *         l       o   o    g           *
c  *         l       o   o    g ggg       *
c  *         l       o   o    g   g       *
c  *         llll     ooo      ggg        *
c  ****************************************
      COMMON / INOUT  / NINP,NOUT  
      SAVE   / INOUT  /
      WRITE(NOUT,*) '==============================================' 
      WRITE(NOUT,*) '=============   LLGFIG   =====================' 
      WRITE(NOUT,*) '==============================================' 
C------------single fermion LL fagmentation-------------------   
C... differences (M.C.-anal.) for O(alf2,alf1) individual beta-s 
      CALL FIG6Y               
      CALL FIG6Z               
      CALL FIG6Z2              
      END





c  ****************************************
c  ****************************************
c  ****************************************
c  *         l        ooo      ggg        *
c  *         l       o   o    g           *
c  *         l       o   o    g ggg       *
c  *         l       o   o    g   g       *
c  *         llll     ooo      ggg        *
c  ****************************************
c  ****************************************
      SUBROUTINE FIG6Y          
C     ****************         
C--------------------------------------------------------------
C                S I N G L E    F R A G M E N A T I O N        
C v*rho(v) distribution        
C Contributions from various beta's in various orders          
C--------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      COMMON / KEYDIS / KEYDIS 
      SAVE   / INOUT  /,/ KEYDIS /
             
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    FIG6Y    =============='      
      WRITE(NOUT,*) ' ====================================='      
C BOOKS HISTO AND FILLS IT WITH A FUNCTION   
      TMIN =  -4D0             
      TMAX =  -0.001D0  
      NBIN = 80                
      KEYDIS=   410            
      CALL DFPLOT("VRHOS",-2010,' O(alf1)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   411            
      CALL DFPLOT("VRHOS",-2011,' O(alf1)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   420            
      CALL DFPLOT("VRHOS",-2020,' O(alf2)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   421            
      CALL DFPLOT("VRHOS",-2021,' O(alf2)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   422            
      CALL DFPLOT("VRHOS",-2022,' O(alf2)       $',NBIN,TMIN,TMAX) 
C--------------------------------------------------------------
C----------------------- O(alf1) ------------------------------
      CALL GPLTIT('FIG6Y  : Beta-s O(alf1)   $')
      CALL GMINIM( 0,-0.080D0)   
      CALL GMAXIM( 0, 0.130D0)   
      CALL GPRINT(310)         
      CALL GPRINT(311)         
C==== Monte Carlo              
      CALL GPLSET('DMOD',1D0)   
      CALL GPLOT( 310,' ',' ',0)             
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT( 311,'S',' ',0)             
C==== analitical               
      CALL GPLSET('DMOD',1D0)   
      CALL GPLOT(2010,'S','*',0)             
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT(2011,'S','*',0)             
C--------------------------------------------------------------
C----------------------- O(alf2) ------------------------------
      CALL GPLTIT('FIG6Y  : Beta-s O(alf2)   $')
      CALL GPRINT(320)         
      CALL GPRINT(321)         
      CALL GPRINT(322)         
C==== PLOTTING                 
C THE BETA2 CONTR. MULTIPLIED BY FACTOR 10   
      CALL GOPERA( 322,'+', 322, 9322,0D0,10D0)                  
      CALL GOPERA( 2022,'+',2022,2022,0D0,10D0)                  
C Monte Carlo                  
      CALL GPLSET('DMOD',1D0)   
      CALL GPLOT(320,' ',' ',0)              
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT(321,'S',' ',0)              
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT(9322,'S',' ',0)             
C Analytical 
      CALL GPLSET('DMOD',1D0)   
      CALL GPLOT(2020,'S','*',0)             
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT(2021,'S','*',0)             
      CALL GPLSET('DMOD',2D0)   
      CALL GPLOT(2022,'S','*',0)             
      CALL GDELET(2922)       
      CALL GDELET(2010)       
      CALL GDELET(2011)       
      CALL GDELET(2020)       
      CALL GDELET(2021)       
      CALL GDELET(2022)       
      END    
      SUBROUTINE FIG6Z         
C     ****************         
C--------------------------------------------------------------
C                S I N G L E    F R A G M E N A T I O N        
C Single fermion fragmentation LL model      
C V a r i o u s    b e t a s   
C differences among MC and analytical results                  
C--------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      COMMON / KEYDIS / KEYDIS 
      SAVE   / INOUT  /,/ KEYDIS /
             
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    FIG6Z    =============='      
      WRITE(NOUT,*) ' ====================================='      
C BOOKS HISTO AND FILLS IT WITH A FUNCTION   
      TMIN =  -4D0             
      TMAX =  -0.001D0  
      NBIN = 80                
      KEYDIS=   410            
      CALL DFPLOT("VRHOS",-2010,' O(alf1)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   411           
      CALL DFPLOT("VRHOS",-2011,' O(alf1)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   420            
      CALL DFPLOT("VRHOS",-2020,' O(alf2)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   421            
      CALL DFPLOT("VRHOS",-2021,' O(alf2)       $',NBIN,TMIN,TMAX) 
      KEYDIS=   422            
      CALL DFPLOT("VRHOS",-2022,' O(alf2)       $',NBIN,TMIN,TMAX) 
C unit line  
      CALL HISUNI(2010,2091) 
C----------------------- O(alf1) ------------------------------
      CALL GPLTIT('FIG6Z   : O(alf1) differences          $')
      CALL GOPERA( 310,'-',2010,2310,1D0,1D0)  
      CALL GOPERA( 311,'-',2011,2311,1D0,1D0)  
      CALL GPRINT(2310)        
      CALL GPRINT(2311)        
      CALL GOPERA( 320,'-',2020,2320,1D0,1D0)  
      CALL GOPERA( 321,'-',2021,2321,1D0,1D0)  
      CALL GOPERA( 322,'-',2022,2322,1D0,1D0)  
      CALL GPRINT(2320)        
      CALL GPRINT(2321)        
      CALL GPRINT(2322)        
      CALL GMINIM( 0,-0.010D0)   
      CALL GMAXIM( 0, 0.010D0)   
      CALL GPLSET('DMOD',2D0)      
      CALL GPLOT(2310,' ',' ',0)  
      CALL GPLSET('DMOD',3D0)      
      CALL GPLOT(2311,'S','*',0)  
C----------------------- O(alf2) ---------------------------------
      CALL GPLTIT('FIG6Z   : O(alf2) differences          $')
      CALL GPLSET('DMOD',2D0)      
      CALL GPLOT(2320,' ',' ',0)         
      CALL GPLSET('DMOD',3D0)             
      CALL GPLOT(2321,'S','*',0)         
      CALL GPLSET('DMOD',4D0)             
      CALL GPLOT(2322,'S','*',0)         
C---=====              
      CALL GDELET(2310)                  
      CALL GDELET(2311)                  
      CALL GDELET(2320)                  
      CALL GDELET(2321)                  
      CALL GDELET(2322)                  
      END              
      SUBROUTINE FIG6Z2                  
C     ***************** 
C------------------------------------------------------------
C Single fermion fragmentation L  model    
C I n f i n i t e   LL  r e s u l t        
C differences among MC and analytical results                
C------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      COMMON / KEYDIS / KEYDIS             
      SAVE   / INOUT  /,/ KEYDIS /
           
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    FIG6Z2   =============='      
      WRITE(NOUT,*) ' ====================================='      
C BOOKS HISTO AND FILLS IT WITH A FUNCTION 
      TMIN =  -4D0           
      TMAX =  -0.001D0
      NBIN = 80              
C--- O(alf2)  YFS analytical final state   
      KEYDIS=   302          
      CALL DFPLOT("VRHOS",-2002,' O(alf2)        ',NBIN,TMIN,TMAX)
C--- O(ALF2)  LL  ANALYTICAL FINAL STATE    
      KEYDIS=   502           
      CALL DFPLOT("VRHOS",-2502,' O(ALF2) LL     ',NBIN,TMIN,TMAX)
C unit line 
      CALL HISUNI(2002,2091)
C (((((((((((                 
      CALL GPRINT( 302)       
      CALL GPRINT( 330)       
      CALL GPRINT(2002)       
C----------------------- O(alf1,2) -----------------------------
      CALL GPLTIT('FIG6Z2 : O(alf2)anatyt. -O(alf.infin).MC.LL. $')
CC    CALL GOPERA( 302,'-',2002,2302,1D0,1D0) 
CC    CALL GOPERA( 330,'-',2502,2330,1D0,1D0) 
      CALL GOPERA(2002,'-', 302,2302,1D0,1D0) 
      CALL GOPERA(2502,'-', 330,2330,1D0,1D0) 
      CALL GOPERA(2502,'/', 330,2339,1D0,1D0) 
      CALL GPRINT(2302)       
      CALL GPRINT(2330)       
      CALL GMINIM( 0,-0.010D0)  
      CALL GMAXIM( 0, 0.010D0)  
C-------
c     CALL GPLSET('DMOD',2D0)    
C     CALL GPLOT(2302,' ',' ',0)              
C     CALL GPLSET('DMOD',3D0)    
C     CALL GPLOT(2330,' ','*',0)              
      CALL GMINIM( 0, 0.990D0)    
      CALL GMAXIM( 0, 1.010D0)    
      CALL GPLSET('DMOD',3D0)    
      CALL GPLOT(2339,' ','*',0)              
C---=====     
      CALL GDELET(2002)         
      CALL GDELET(2502)         
      CALL GDELET(2330)         
      END     


