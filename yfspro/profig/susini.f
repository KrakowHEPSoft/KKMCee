*////////////////////////////////////////////////////////////////////////////
*//   make susini
*////////////////////////////////////////////////////////////////////////////
*  ****************************************
*  *           sss   u   u   sss          *
*  *          s      u   u  s             *
*  *           ss    u   u   ss           *
*  *             s   u   u     s          *
*  *          sss     uuu   sss           *
*  ****************************************
* WARNING: FOR KeyDis=410-422 DELB MODIFIED &&&&               
* WARNING NUMBERING FOR YFS KeyDis CHANGED TO 301,302,303 
*------------------------------------------------------   
* Convention for KeyDis      
* Pedagogical exercises
*     KeyDis   =  1     soft part YFS       First  Order
*     KeyDis   =  2     soft part YFS       Second Order
*     KeyDis   =  5     hard non-exp.       First  Order
*     KeyDis   =  6     hard non-exp.       Second Order
*     KeyDis   =  9     reference distr. of YFS paper
* Total results
*     KeyDis   =  0 + R*100                 Zero   Order     
*     KeyDis   =  1 + R*100                 First  Order     
*     KeyDis   =  2 + R*100                 Second Order     
* Beta contributions
*     KeyDis   = 10 + R*100     Beta0       Zero   Order        
*     KeyDis   = 11 + R*100     Beta0       First  Order        
*     KeyDis   = 12 + R*100     Beta1            
*     KeyDis   = 20 + R*100     Beta0       Second Order        
*     KeyDis   = 21 + R*100     Beta1            
*     KeyDis   = 22 + R*100     Beta2        
*     R = 100 Kuraev-Fadin   
*     R = 300 YFS            
*     R = 400 YFS single electron LL str. funct.          
*------------------------------------------------------   
      PROGRAM MAIN
*  ****************************************
      IMPLICIT NONE
      SAVE
*---
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*---
      INTEGER imax
      PARAMETER(imax=1000)
      REAL*8          xpar
      COMMON  /xpar/  xpar(imax)
*---
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
      INTEGER lint
*----------------------------------------------------------
      ninp=  5           
      nout= 16           
      Tesnam    = 'susini'
      TeXfile   = 'susini.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

c      Dname  = '../ini40/ini40.data'    ! current
c      Hname  = '../ini40/pro.hst '      ! current

      Dname  = '../ini40/5M.ini40.data' ! standard KeyZet=2
      Hname  = '../ini40/5M.pro.hst '   ! standard KeyZet=2

*=====================================
* Read data, the same as in MC run
      CALL ReaDat(Dname,imax,xpar)
      CALL Semalib_Initialize(xpar)
* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
*=====================================
* Initialize GLK_Plot
      Lint=0
      CALL GLK_PlInitialize(Lint,TeXfile)
*=====================================================

*------------++++++++++
* Initial state !!!
* The distribution v*rho(v) as a function of log(v)   
* O(alf1), O(alf0) and O(alf1)-O(alf0)         
* O(alf1), O(alf0) and O(alf1)-O(alf0)         
      CALL RHTOTI      
*------------++++++++++
* The distribution v*rho(v) as function of log(v)         
* Contributions from various beta's in O(alf1) and Oalf(2)           
      CALL RHBETI      
*------------++++++++++
* The distribution v*rho(v) as function of log(v)         
* Differences (M.C. - analyt.) for O(alf2), O(alf1)
      CALL RHDIFI              
*----------------------
* Plots for Sussex talk  
*------------++++++++++              
* Pedagogical plots explaining the meaning of ad-hoc exponentiation
      CALL ADHOC1     
*----------------------
* The distribution,  rho(v)/reference,  like in YFS paper    
* O(alf2)  O(alf1) cases, seems to work!           
      CALL RHOCPC
*---------------------              
* Energies of fastest 3 photons
      CALL ORDENE
*----------------------
*===================================================================
* end GLK_Plot, close LaTeX file
      CALL GLK_PlEnd
*===================================================================
* Write all histograms into dump file, for control/debug
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=================================
      CLOSE(nout)              
      END



*  ****************************************
*  ****************************************
*  ****************************************
*  *           sss   u   u   sss          *
*  *          s      u   u  s             *
*  *           ss    u   u   ss           *
*  *             s   u   u     s          *
*  *          sss     uuu   sss           *
*  ****************************************
*  ****************************************
      SUBROUTINE RHTOTI             
*     *****************              
*-------------------------------------------------------------------
* The distribution v*rho(v) as a function of log(v)   
* O(alf1), O(alf0) and O(alf1)-O(alf0)         
* O(alf1), O(alf0) and O(alf1)-O(alf0)         
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON  /xpar/  xpar(1000)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
c[[[      CHARACTER*80 TITLE
    
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHTOTI   =============='      
      WRITE(NOUT,*) ' ====================================='      

      cmsene=xpar(1)
      born = Semalib_Born(cmsene**2)
* Renormalize M.C. results and divide by Born
      IDA=10000
      IDB=50000
      IDD=15000
c[[[      CALL GLK_hinbo1(IDA+300,TITLE,NBIN,TMIN,TMAX)
c[[[      CALL GLK_hinbo1(IDB+300,TITLE,NBIV,VMIN,VMAX)
      CALL RHONOR(IDA+300,IDD+300)
      CALL RHONOR(IDA+301,IDD+301)
      CALL RHONOR(IDA+302,IDD+302)
      CALL RHONOR(IDA+305,IDD+305)
      CALL RHONOR(IDA+306,IDD+306)
* BOOKS HISTO AND FILLS IT WITH A FUNCTION        
*--- O(alf0)      
      KeyDis=   300                 
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2300,' O(alf0)       $',IDA+300)
*--- O(alf1)      
      KeyDis=   301                 
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2301,' O(alf1)       $',IDA+300)
*--- O(alf2)      
      KeyDis=   302                 
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2302,' O(alf2)       $',IDA+300)
*--- Differences  
      CALL GLK_Operat(2301,'-',2300,2911,1D0,1D0)       
      CALL GLK_Operat(2302,'-',2301,2912,1D0,1D0)       
*--------------------------------   
*--- v*rho(v)     
*--- O(alf1), O(alf0) and O(alf1)-O(alf0)         
*--------------------------------   
      CALL GLK_PLTITLE('RHTOTI: O(alf1), O(alf0) and O(alf1)-O(alf0)$')
      CALL GLK_YMINIM( 0,-0.06D0)         
      CALL GLK_YMAXIM( 0, 0.12D0)         
      CALL GLK_PRINT(IDD+300)              
      CALL GLK_PRINT(IDD+301)              
      CALL GLK_PRINT(IDD+305)              
      CALL GLK_PRINT(2300)             
      CALL GLK_PRINT(2301)             
      CALL GLK_PRINT(2911)             
* Nonte Carlo     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+300,' ',' ',0)     
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(IDD+301,'S',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+305,'S',' ',0)     
* Analytical      
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2300,'S','*',0)    
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(2301,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2911,'S','*',0)    
*--------------------------------   
*--- v*rho(v)     
*--- O(alf2), O(alf1) and O(alf2)-O(alf1)         
*--------------------------------   
      CALL GLK_PLTITLE('RHTOTI: O(alf2) and (O(alf2)-O(alf1))X10 $')
      CALL GLK_YMINIM( 0,-0.01D0)         
      CALL GLK_YMAXIM( 0, 0.08D0)         
* THE DIFFERENCE O(ALF2)-O(ALF1) MULTIPLIED BY FACTOR 10            
      CALL GLK_Operat(   2912,'+',    2912,   2912,0D0,10D0)      
      CALL GLK_Operat(IDD+306,'+', IDD+306,IDD+306,0D0,10D0)      
*
      CALL GLK_PRINT(IDD+301)              
      CALL GLK_PRINT(IDD+302)              
      CALL GLK_PRINT(IDD+306)              
      CALL GLK_PRINT(2301)              
      CALL GLK_PRINT(2302)              
      CALL GLK_PRINT(2912)              
* Monte Carlo     
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(IDD+301,' ',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+302,'S',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+306,'S',' ',0)    
* Analytical      
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(2301,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2302,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2912,'S','*',0)    
*------         
      CALL GLK_DELET(2300)       
      CALL GLK_DELET(2301)       
      CALL GLK_DELET(2302)       
      CALL GLK_DELET(2911)       
      CALL GLK_DELET(2912)       
      CALL GLK_DELET(IDD+300)
      CALL GLK_DELET(IDD+301)
      CALL GLK_DELET(IDD+302)
      CALL GLK_DELET(IDD+305)
      CALL GLK_DELET(IDD+306)

*************************************************
****             NEW PART 1993              *****
*************************************************
****************** chi(log_xmax) ********************
      KeyDis=   301                 
      CALL Semalib_VVplot(KeyDis,"VCHI2",-100301,' O(alf2) ini+fin $',IDA+300)
      CALL GLK_Print(100301)
      KeyDis=   302                 
      CALL Semalib_VVplot(KeyDis,"VCHI2",-100302,' O(alf2) ini+fin $',IDA+300)
      CALL GLK_Print(100302)
      CALL GLK_Operat(100302,'-',100301,200306,1d0,1d0)       
****************** chi(xmax)     ********************
      KeyDis=   301                 
      CALL Semalib_VVplot(KeyDis,"VCHI2", 500301,' O(alf2) ini+fin $',IDB+300)
      CALL GLK_Print(500301)
      KeyDis=   302                 
      CALL Semalib_VVplot(KeyDis,"VCHI2", 500302,' O(alf2) ini+fin $',IDB+300)
      CALL GLK_Print(500302)
      CALL GLK_Operat(500302,'-',500301,600306,1D0,1D0)       
***************** rho(log_x) ************************
      KeyDis=   301
      CALL Semalib_VVplot(KeyDis,"VRHO2",-700301,' O(alf2) ini $',IDA+300)
      KeyDis=   302
      CALL Semalib_VVplot(KeyDis,"VRHO2",-700302,' O(alf2) ini $',IDA+300)
      CALL GLK_Operat(700302,'-',700301,800306,1D0,1D0)       

*** TECHNICAL/PHYSICAL PRECISION chi(log_xmax)
      YMIN = -0.010D0
      YMAX =  0.015D0
      XMAG =   1
      YMAG = 0.01
      CALL PLTECH("CUMU",
     $"RHTOTI: sig(log(vmax))/Born, Phys/Tech prec., YMAG = 0.01$",
     $ IDA+302, IDA+306, 100302, 200306,YMIN,YMAX,XMAG,YMAG,BORN)
*** TECHNICAL/PHYSICAL PRECISION chi(xmax)
      YMIN = -0.010D0
      YMAX =  0.015D0
      XMAG =   1
      YMAG = 0.01
      CALL PLTECH("CUMU",
     $"RHTOTI: sig(vmax)/Born, Phys/Tech prec., YMAG = 0.01$",
     $ IDB+302, IDB+306, 500302, 600306,YMIN,YMAX,XMAG,YMAG,BORN)
*** TECHNICAL/PHYSICAL PRECISION rho(log_x)
      YMIN = -0.005D0
      YMAX =  0.0200D0
      XMAG = 10
      YMAG = 1
      CALL PLTECH("NB10",
     $"RHTOTI: dsig/dlogx [nb] Phys/Tech. prec. , XMAG = 10$",
     $ IDA+302, IDA+306, 700302, 800306,YMIN,YMAX,XMAG,YMAG,1d0)
      END     

      SUBROUTINE RHBETI         
*     *****************          
*---------------------------------------------------------------
* The distribution v*rho(v) as function of log(v)         
* Contributions from various beta's in O(alf1) and Oalf(2)           
*---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
c[[[      CHARACTER*80 TITLE
              
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHBETI   =============='      
      WRITE(NOUT,*) ' ====================================='      
      
* Renormalize M.C. results and divide by Born
      IDA=10000
      IDD=15000
      CALL RHONOR(IDA+310,IDD+310)
      CALL RHONOR(IDA+311,IDD+311)
      CALL RHONOR(IDA+320,IDD+320)
      CALL RHONOR(IDA+321,IDD+321)
      CALL RHONOR(IDA+322,IDD+322)
* BOOKS HISTO AND FILLS IT WITH A FUNCTION    
c[[[      CALL GLK_hinbo1(IDA+310,TITLE,NBIN,TMIN,TMAX)
      KeyDis=   310             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2310,' O(alf1)       $',IDA+300)
      KeyDis=   311    
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2311,' O(alf1)       $',IDA+300)
      KeyDis=   320             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2320,' O(alf2)       $',IDA+300)
      KeyDis=   321             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2321,' O(alf2)       $',IDA+300)
      KeyDis=   322             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2322,' O(alf2)       $',IDA+300)
*---------------------------------------------------------------
*----------------------- O(alf1) -------------------------------
      CALL GLK_PLTITLE('RHBETI: beta-s,    O(alf1)    $')           
      CALL GLK_YMINIM( 0,-0.080D0)    
      CALL GLK_YMAXIM( 0, 0.130D0)    
      CALL GLK_PRINT(IDD+310)          
      CALL GLK_PRINT(IDD+311)          
      CALL GLK_PRINT(   2310)          
      CALL GLK_PRINT(   2311)          
*==== Monte Carlo               
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT( IDD+310,' ',' ',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT( IDD+311,'S',' ',0)              
*==== analitical                
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2310,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2311,'S','*',0)              
*---------------------------------------------------------------
*----------------------- O(alf2) -------------------------------
* THE BETA2 CONTR. MULTIPLIED BY FACTOR 10    
      CALL GLK_PLTITLE('RHBETI: beta-s,    O(alf2)    $')           
      CALL GLK_Operat( IDD+322,'+', IDD+322,IDD+322,0D0,10D0) 
      CALL GLK_Operat( 2322,   '+',    2322,   2322,0D0,10D0) 
      CALL GLK_PRINT(IDD+320)          
      CALL GLK_PRINT(IDD+321)          
      CALL GLK_PRINT(IDD+322)          
      CALL GLK_PRINT(2320)          
      CALL GLK_PRINT(2321)          
      CALL GLK_PRINT(2322)          
*==== PLOTTING                  
* Monte Carlo 
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(IDD+320,' ',' ',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(IDD+321,'S',' ',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(IDD+322,'S',' ',0)              
* Analytical  
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2320,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2321,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2322,'S','*',0)              
*
      CALL GLK_DELET(2310)       
      CALL GLK_DELET(2311)       
      CALL GLK_DELET(2320)       
      CALL GLK_DELET(2321)       
      CALL GLK_DELET(2322)       
      CALL GLK_DELET(IDD+310)          
      CALL GLK_DELET(IDD+311)          
      CALL GLK_DELET(IDD+320)          
      CALL GLK_DELET(IDD+321)          
      CALL GLK_DELET(IDD+322)          
      END     

      SUBROUTINE RHDIFI         
*     *****************          
*---------------------------------------------------------------
* The distribution v*rho(v) as function of log(v)         
* Differences (M.C. - analyt.) for O(alf2), O(alf1)
*---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
c[[[      CHARACTER*80 TITLE
              
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHDIFI   =============='      
      WRITE(NOUT,*) ' ====================================='      

* Renormalize M.C. results and divide by Born
      IDA=10000
      IDD=15000
      CALL RHONOR(IDA+310,IDD+310)
      CALL RHONOR(IDA+311,IDD+311)
      CALL RHONOR(IDA+320,IDD+320)
      CALL RHONOR(IDA+321,IDD+321)
      CALL RHONOR(IDA+322,IDD+322)
      CALL GLK_PRINT(IDD+310)
      CALL GLK_PRINT(IDD+311)
      CALL GLK_PRINT(IDD+320)
      CALL GLK_PRINT(IDD+321)
      CALL GLK_PRINT(IDD+322)
* BOOKS HISTO AND FILLS IT WITH A FUNCTION    
c[[[      CALL GLK_hinbo1(IDA+310,TITLE,NBIN,TMIN,TMAX)
      KeyDis=   310             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2310,' O(alf1)       $',IDA+300)
      KeyDis=   311             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2311,' O(alf1)       $',IDA+300)
      KeyDis=   320             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2320,' O(alf2)       $',IDA+300)
      KeyDis=   321             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2321,' O(alf2)       $',IDA+300)
      KeyDis=   322             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2322,' O(alf2)       $',IDA+300)
* unit line   
c[[[      CALL HISUNI(IDA+310,2091)
*----------------------- O(alf1) -------------------------------
      CALL GLK_Operat( IDD+310,'-',2310,5310,1D0,1D0)   
      CALL GLK_Operat( IDD+311,'-',2311,5311,1D0,1D0)   
      CALL GLK_PRINT(5310)         
      CALL GLK_PRINT(5311)         
      CALL GLK_Operat( IDD+320,'-',2320,5320,1D0,1D0)   
      CALL GLK_Operat( IDD+321,'-',2321,5321,1D0,1D0)   
      CALL GLK_Operat( IDD+322,'-',2322,5322,1D0,1D0)   
      CALL GLK_PRINT(5320)         
      CALL GLK_PRINT(5321)         
      CALL GLK_PRINT(5322)         
      CALL GLK_YMINIM( 0,-0.010D0)    
      CALL GLK_YMAXIM( 0, 0.010D0)    
*----------------------- O(alf1) -------------------------------
      CALL GLK_PLTITLE('RHDIFI:  (M.C. - analyt.),  O(alf1)  $') 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(5310,' ',' ',0)              
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(5311,'S','*',0)              
*----------------------- O(alf2) -------------------------------
      CALL GLK_PLTITLE('RHDIFI:  (M.C. - analyt.),  O(alf2)  $') 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(5320,' ',' ',0)              
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(5321,'S','*',0)              
      CALL GLK_PLSET('DMOD',4D0)    
      CALL GLK_PLOT(5322,'S','*',0)              
*----           
c[[[      CALL GLK_DELET(2091)       
      CALL GLK_DELET(5310)       
      CALL GLK_DELET(5311)       
      CALL GLK_DELET(5320)       
      CALL GLK_DELET(5321)       
      CALL GLK_DELET(5322)       
      CALL GLK_DELET(2310)       
      CALL GLK_DELET(2311)       
      CALL GLK_DELET(2320)       
      CALL GLK_DELET(2321)       
      CALL GLK_DELET(2322)       
      CALL GLK_DELET(IDD+310)
      CALL GLK_DELET(IDD+311)
      CALL GLK_DELET(IDD+320)
      CALL GLK_DELET(IDD+321)
      CALL GLK_DELET(IDD+322)
      END     

      SUBROUTINE ADHOC1        
*     *****************          
*------------------------------------------------------------------
* Pedagogical plots explaining the meaning of ad-hoc exponentiation
* rho(v) distribution, O(alf1) and O(alf2) cases: alanyt. + MC. 
*------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  


      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    ADHOC1   =============='      
      WRITE(NOUT,*) ' ====================================='      
              
* Renormalize M.C. results and divide by Born
      IDA=10000
      IDD=15000
      CALL RHONOR(IDA+301,IDD+301)
      CALL RHONOR(IDA+302,IDD+302)
* BOOKS HISTO AND FILLS IT WITH A FUNCTION      
      TMIN =  -4D0             
      TMAX =  -0.001D0           
      NBIN = 100               
*--- O(ALF1) 
      KeyDis=   101              
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5001,' SOFT          $',IDA+300) 
      KeyDis=   105          
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5005,' HARD          $',IDA+300) 
      KeyDis=   301          
      CALL Semalib_VVplot(KeyDis,"VRHO ",-5301,' O(ALF)  YFS   $',IDA+300) 
      KeyDis=   201          
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5101,' O(ALF)  K-F   $',IDA+300) 
*--- O(ALF2)               
      KeyDis=   102          
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5002,' SOFT          $',IDA+300) 
      KeyDis=   106          
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5006,' HARD          $',IDA+300) 
      KeyDis=   302         
      CALL Semalib_VVplot(KeyDis,"VRHO ",-5302,' O(ALF2) YFS   $',IDA+300) 
      KeyDis=   202          
      CALL Semalib_VVplot(KeyDis,"VRHOS",-5102,' O(ALF2) K-F   $',IDA+300) 
*-- O(alf1) First plot ----------               
      CALL GLK_PLTITLE('ADHOC1: initial O(alf1)         $')           
      CALL GLK_YMINIM( 0,0D0)  
      CALL GLK_YMAXIM( 0,.13D0)                  
      CALL GLK_PRINT(5001)      
      CALL GLK_PRINT(5005)      
      CALL GLK_PRINT(5301)     
      CALL GLK_PRINT(5101)        
      CALL GLK_PRINT(IDD+301)       
*==== PLOTTING O(ALF1)       
      CALL GLK_PLSET('DMOD',1D0) 
      CALL GLK_PLOT(IDD+301,' ',' ',0)            
      CALL GLK_PLSET('DMOD',1D0) 
      CALL GLK_PLOT(5001,'S','*',0)             
      CALL GLK_PLSET('DMOD',4D0) 
      CALL GLK_PLOT(5005,'S','*',0)             
      CALL GLK_PLSET('DMOD',2D0) 
      CALL GLK_PLOT(5301,'S','*',0)             
      CALL GLK_PLSET('DMOD',3D0) 
      CALL GLK_PLOT(5101,'S','*',0)             
*-- O(alf2)  Second plot ---------                 
      CALL GLK_PLTITLE('ADHOC2: initial O(alf2)          $')           
      CALL GLK_YMINIM( 0,0D0)    
      CALL GLK_YMAXIM( 0,.13D0)   
      CALL GLK_PRINT(5002)        
      CALL GLK_PRINT(5006)               
      CALL GLK_PRINT(5302)               
      CALL GLK_PRINT(IDD+302)              
*==== PLOTTING O(ALF2)              
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(IDD+302,' ',' ',0)     
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(5002,'S','*',0)      
      CALL GLK_PLSET('DMOD',4D0)        
      CALL GLK_PLOT(5006,'S','*',0)      
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(5302,'S','*',0)      
* cleaning
      CALL GLK_DELET(5001)               
      CALL GLK_DELET(5005)               
      CALL GLK_DELET(5301)               
      CALL GLK_DELET(5101)               
      CALL GLK_DELET(5002)               
      CALL GLK_DELET(5006)               
      CALL GLK_DELET(5302)               
      CALL GLK_DELET(5102)               
      CALL GLK_DELET(IDD+301)               
      CALL GLK_DELET(IDD+302)               
      END         

      SUBROUTINE RHOCPC         
*     *****************          
*---------------------------------------------------------------
* historical plot
* The distribution     rho(v)/reference       
* O(alf2)  O(alf1)    MC and analytical            
*---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  

      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHOCPC   =============='      
      WRITE(NOUT,*) ' ====================================='      
              
* Renormalize M.C. results and divide by Born
      IDA=10000
      IDD=15000
      CALL RHONOR(IDA+301,IDD+301)
      CALL RHONOR(IDA+302,IDD+302)
* BOOKS HISTO AND FILLS IT WITH A FUNCTION    
      TMIN =  -4D0              
      TMAX =  -0.001D0     
      NBIN = 200                
      KeyDis=   301             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2301,' O(alf1)       $',IDA+300)
      KeyDis=   302             
      CALL Semalib_VVplot(KeyDis,"VRHO ",-2302,' O(alf2)       $',IDA+300)
* reference for plots           
      KeyDis=     9             
      CALL Semalib_VVplot(KeyDis,"VRHOS",-2009,' reference     $',IDA+300)
* Reference for histos          
      TMIN = -4D0               
      TMAX = -0.000001D0          
      NBT  = 80                 
      KeyDis=     9             
      CALL Semalib_VVplot(KeyDis,"VRHOS",-2099,' reference     $',IDA+300)

c[[[      CALL HISUNI(2301,2090)
*----------------------------------           
* Divide by refrence            
      CALL GLK_Operat( 2301, '/' ,2009,   2301,1D0,1D0)       
      CALL GLK_Operat( 2302, '/' ,2009,   2302,1D0,1D0)       
      CALL GLK_Operat(  IDD+301, '/' ,2099,   5301,1D0,1D0)     
      CALL GLK_Operat(  IDD+302, '/' ,2099,   5302,1D0,1D0)     
*--- v*rho(v) 
*--- O(alf2), O(alf1) and O(alf2)-O(alf1)     
      CALL GLK_YMINIM( 0, 1.02D0)     
      CALL GLK_YMAXIM( 0, 1.14D0)     
      CALL GLK_PRINT(2301)          
      CALL GLK_PRINT(2302)          
      CALL GLK_PLTITLE('  RHOCPC           $')           
c[[[      CALL GLK_PLOT(2090,' ','B',0) 
* Monte Carlo 
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(5301,'S','B',0) 
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(5302,'S','B',0) 
* Analytical  
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2301,'S','*',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2302,'S','*',0) 
*leaning
      CALL GLK_DELET(2301)
      CALL GLK_DELET(2302)
      CALL GLK_DELET(2009)
      CALL GLK_DELET(2099)
c[[[      CALL GLK_DELET(2090)
      call GLK_Delet(5301)
      call GLK_Delet(5302)
      CALL GLK_DELET(IDD+301)               
      CALL GLK_DELET(IDD+302)               
      END     


      SUBROUTINE ORDENE         
*     *****************           
*---------------------------------------------------------------
* Historical Susex plot
* this looks like energies of photons ID=30+i , 50+i
*---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
      EXTERNAL TDIS             
              
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    ORDENE   =============='      
      WRITE(NOUT,*) ' ====================================='      

      IDE1 = 10050
      IDE  = 10750
      DO J= 1,3
         CALL GLK_RenHst("UN10",6,ide1+j,ide+j)
      ENDDO
      CALL GLK_PLTITLE('ORDENE : Energies of 3 fastest photons  $')  
* BOOKS HISTO AND FILLS IT WITH A FUNCTION    
      TMIN =  -14D0             
      TMAX =  -0.001D0      
      NBIN =   120                
      KeyDis=   1               
      CALL GLK_BOOKFUN1(2031,' SECOND ORDER  $',NBIN,TMIN,TMAX,TDIS)    
      KeyDis=   2               
      CALL GLK_BOOKFUN1(2032,'               $',NBIN,TMIN,TMAX,TDIS)    
      KeyDis=   3               
      CALL GLK_BOOKFUN1(2033,'               $',NBIN,TMIN,TMAX,TDIS)    
      CALL GLK_YMINIM( 0,0D0)       
      CALL GLK_YMAXIM( 0,.120D0)      
      CALL GLK_PRINT(2031)           
      CALL GLK_PRINT(2032)           
      CALL GLK_PRINT(2033)           
      CALL GLK_PRINT(IDE+1)          
      CALL GLK_PRINT(IDE+2)          
      CALL GLK_PRINT(IDE+3)          
*==== PLOTTING Monte Carlo                             
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(IDE +1 ,' ',' ',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(IDE +2 ,'S',' ',0) 
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(IDE +3 ,'S',' ',0) 
*---- Analytical
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(2031,'S','*',0)  
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2032,'S','*',0)  
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2033,'S','*',0)  
      CALL GLK_DELET(2031)
      CALL GLK_DELET(2032)
      CALL GLK_DELET(2033)
      END     


 
      FUNCTION TDIS(T)          
*     ****************          
* Called in SUS part only, in ORDENE
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER( pi=3.1415926535897932d0)     
      COMMON  /xpar/  xpar(1000)
      SAVE
*
      alfinv = xpar(30)
      cmsene=xpar(1)
      svar   = cmsene**2
      alfinv = xpar(30)
      alf1   = 1d0/pi/alfinv

      KFbeam = 11           ! KF=11 is electron
      ke = 500+10*KFbeam
      amel   = xpar(ke+6)

      bilg   = dlog(svar/amel**2)             
      beti   = 2d0*alf1*(bilg-1d0)            
*             
      ALPHA= BETI*LOG(10D0)      
      TT=ABS(T)                 
      NN=KeyDis-1               
      DIS=ALPHA*EXP(-ALPHA*TT)/LOG(10D0)       
      DO 10 I=1,NN              
   10 DIS=DIS*ALPHA*TT/I        
      TDIS =DIS                 
      END     



      SUBROUTINE rhonor(id1,id2)
*     **************************
* Called in SUS part only!
* Renormalize and DIVIDE OFF Born(s(1-v))
* WARNING!!! BORNV undefined !!!!
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.1415926535897932D0,ALFINV=137.03604D0)
      COMMON  /xpar/  xpar(1000)
      SAVE
*-------------------------------------------------------
      CHARACTER*80 TITLE
      DIMENSION  HSTG(300)
*-------------------------------------------------------
*
      cmsene=xpar(1)
      KeyZet = xpar(501)

      IdGen = 6
      CALL GLK_WtMon(   1,idgen,par1,par2,par3 )
      nevt    =  par3
      xscrnb  =  par1
*
      gnanob=389.385d-30*1.d33
      sig0nb =  4d0*pi/(alfinv**2*3d0*cmsene**2)*gnanob
      XSECR  =  xscrnb/sig0nb
*
      CALL GLK_UnPak(id1,hstg,'  ',1)
      CALL GLK_hinbo1(id1,title,nbt,tmin,tmax)
      fact = nbt*xsecr/(nevt*(tmax-tmin)*log(10.))
* Divide histogram by Born
      DO 500 IB=1,NBT
         T= TMIN+(TMAX-TMIN)*(IB-0.5)/NBT
         VV= 10.0**T
         SVAR = CMSENE**2
         BORV =Semalib_Borny(SVAR*(1D0-VV))/(1.-VV)
* constant x-section for tests
         IF(KEYZET.EQ.-2)  BORNV = BORNV* (1D0-VV) !!!!!??????
         HSTG(IB)= HSTG(IB)*FACT/BORV
  500 CONTINUE
* Store result into ID2
      IF(ID2.EQ.ID1) THEN
         CALL GLK_RESET(ID1,' ')
         CALL GLK_PAK(ID1,HSTG)
      ELSE
         CALL GLK_BOOK1(ID2,TITLE,NBT,TMIN,TMAX)
         CALL GLK_PAK(ID2,HSTG)
      ENDIF
      END
