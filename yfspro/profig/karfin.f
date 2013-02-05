*///////////////////////////////////////////////////////////////
*//
*//   make karfin
*//
*//
*///////////////////////////////////////////////////////////////
c  ****************************************
c  *          k  k    aaa    rrrr         *
c  *          k k    a   a   r   r        *
c  *          kk     a   a   rrrr         *
c  *          k k    aaaaa   r  r         *
c  *          k  k   a   a   r   r        *
c  ****************************************
C WARNING: FOR KeyDis=410-422 DELB MODIFIED &&&&               
C WARNING NUMBERING FOR YFS KeyDis CHANGED TO 301,302,303 
C------------------------------------------------------   
C Convention for KeyDis      
C Pedagogical exercises
C     KeyDis   =  1     soft part YFS       First  Order
C     KeyDis   =  2     soft part YFS       Second Order
C     KeyDis   =  5     hard non-exp.       First  Order
C     KeyDis   =  6     hard non-exp.       Second Order
C     KeyDis   =  9     reference distr. of YFS paper
C Total results
C     KeyDis   =  0 + R*100                 Zero   Order     
C     KeyDis   =  1 + R*100                 First  Order     
C     KeyDis   =  2 + R*100                 Second Order     
C Beta contributions
C     KeyDis   = 10 + R*100     Beta0       Zero   Order        
C     KeyDis   = 11 + R*100     Beta0       First  Order        
C     KeyDis   = 12 + R*100     Beta1            
C     KeyDis   = 20 + R*100     Beta0       Second Order        
C     KeyDis   = 21 + R*100     Beta1            
C     KeyDis   = 22 + R*100     Beta2        
C     R = 100 Kuraev-Fadin   
C     R = 300 YFS            
C     R = 400 YFS single electron LL str. funct.          
C------------------------------------------------------   

      PROGRAM MAIN
*     ***********************************      
      IMPLICIT NONE
      SAVE
*---
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*---
      INTEGER imax
      PARAMETER(imax=1000)
      REAL*8  xpar(imax)
*---
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
      INTEGER lint
*-------------------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'karfin'
      TeXfile   = 'karfin.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

c      Dname  = '../fin40/fin40.data'  ! current
c      Hname  = '../fin40/pro.hst '    ! current

      Dname  = '../fin40/5M.fin40.data'  ! standard KeyZet=2
      Hname  = '../fin40/5M.pro.hst '    ! standard KeyZet=2

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

C Final state !!!!
C The distribution v*rho(v) as a function of log(v)   
C O(alf1), O(alf0) and O(alf1)-O(alf0)         
C O(alf1), O(alf0) and O(alf1)-O(alf0)         
      CALL RHTOTF      
C------------++++++++++
C The distribution v*rho(v) as function of log(v)         
C Contributions from various beta's in O(alf1) and Oalf(2)           
      CALL RHBETF      
C------------++++++++++
C The distribution v*rho(v) as function of log(v)         
C Differences (M.C. - analyt.) for O(alf2), O(alf1)
      CALL RHDIFF              
C----------------------
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


      SUBROUTINE RHTOTF             
C     *****************              
C-------------------------------------------------------------------
C Final state!!!!
C The distribution v*rho(v) as a function of log(v)   
C O(alf1), O(alf0) and O(alf1)-O(alf0)         
C O(alf1), O(alf0) and O(alf1)-O(alf0)         
C-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF 
      COMMON / INOUT  / NINP,NOUT  
      SAVE
    
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHTOTF   =============='      
      WRITE(NOUT,*) ' ====================================='      
      IDA = 20000
      IDB = 50000
      IDD = 25000
      BORN = Semalib_Born(4*ENE**2)
      CALL RENHIF(IDA+300,IDD+300)
      CALL RENHIF(IDA+301,IDD+301)
      CALL RENHIF(IDA+302,IDD+302)
      CALL RENHIF(IDA+305,IDD+305)
      CALL RENHIF(IDA+306,IDD+306)
C BOOKS HISTO AND FILLS IT WITH A FUNCTION        
C--- O(alf0)      
      KeyDis=   300                 
      CALL Semalib_VVplot(KeyDis,"URHO ",-2300,' O(alf0) $',IDA+300)
C--- O(alf1)      
      KeyDis=   301                 
      CALL Semalib_VVplot(KeyDis,"URHO ",-2301,' O(alf1) $',IDA+300)
C--- O(alf2)
      KeyDis=   302                
      CALL Semalib_VVplot(KeyDis,"URHO ",-2302,' O(alf2) $',IDA+300)
C--- Differences  
      CALL GLK_OPERAT(2301,'-',2300,2911,1D0,1D0)       
      CALL GLK_OPERAT(2302,'-',2301,2912,1D0,1D0)       
C---------------------------------- 
C--- v*rho(v)     
C--- O(alf1), O(alf0) and O(alf1)-O(alf0)         
C--------------------------------   
      CALL GLK_PLTITLE('RHTOTF: O(alf1), O(alf0) and O(alf1)-O(alf0)$')
      CALL GLK_YMINIM( 0,-0.06D0)         
      CALL GLK_YMAXIM( 0, 0.12D0)         
      CALL GLK_PRINT(2300)             
      CALL GLK_PRINT(2301)             
      CALL GLK_PRINT(2911)             
      CALL GLK_PRINT(IDD+300)              
      CALL GLK_PRINT(IDD+301)              
      CALL GLK_PRINT(IDD+305)              
C Nonte Carlo     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+300,' ',' ',0)     
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(IDD+301,'S',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+305,'S',' ',0)     
C Analytical      
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2300,'S','*',0)    
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(2301,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2911,'S','*',0)    
C--------------------------------   
C--- v*rho(v)     
C--- O(alf2), O(alf1) and O(alf2)-O(alf1)         
C--------------------------------   
      CALL GLK_PLTITLE('RHTOTF: O(alf2) and (O(alf2)-O(alf1))*5 $')
      CALL GLK_YMINIM( 0,-0.01D0)         
      CALL GLK_YMAXIM( 0, 0.08D0)         
C THE DIFFERENCE O(ALF2)-O(ALF1) MULTIPLIED BY FACTOR XFAC
      XFAC=5D0           
      CALL GLK_OPERAT(2912,'+',2912,2912,0D0,XFAC)      
      CALL GLK_OPERAT( IDD+306,'+', IDD+306,IDD+306,0D0,XFAC)      
      CALL GLK_PRINT(IDD+301)              
      CALL GLK_PRINT(IDD+302)              
      CALL GLK_PRINT(IDD+306)              
      CALL GLK_PRINT(2301)             
      CALL GLK_PRINT(2302)             
      CALL GLK_PRINT(2912)             
C Monte Carlo     
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(IDD+301,' ',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+302,'S',' ',0)     
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(IDD+306,'S',' ',0)    
C Analytical      
      CALL GLK_PLSET('DMOD',1D0)        
      CALL GLK_PLOT(2301,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2302,'S','*',0)    
      CALL GLK_PLSET('DMOD',2D0)        
      CALL GLK_PLOT(2912,'S','*',0)    
C------         
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
      WRITE(6,*) ' OLD ended '
*************************************************
****             NEW PART 1993              *****
*************************************************
****************** chi(log_xmax) ********************
      KeyDis=   301
      CALL Semalib_VVplot(KeyDis,"UCHI2",-100301,' O(alf2) fin $',IDA+300)
      CALL GLK_Print(100301)
      KeyDis=   302
      CALL Semalib_VVplot(KeyDis,"UCHI2",-100302,' O(alf2) fin $',IDA+300)
      CALL GLK_Print(100302)
      CALL GLK_Operat(100302,'-',100301,200306,1D0,1D0)       
      CALL GLK_Print(200306)
****************** chi(xmax)     ********************
      KeyDis=   301
      CALL Semalib_VVplot(KeyDis,"UCHI2", 500301,' O(alf2) fin $',IDB+300)
      CALL GLK_Print(500301)
      KeyDis=   302
      CALL Semalib_VVplot(KeyDis,"UCHI2", 500302,' O(alf2) fin $',IDB+300)
      CALL GLK_Print(500302)
      CALL GLK_Operat(500302,'-',500301,600306,1D0,1D0)       
      CALL GLK_Print(600306)
***************** rho(log_x) ************************
      KeyDis=   301
      CALL Semalib_VVplot(KeyDis,"URHO2",-700301,' O(alf2) fin $',IDA+300)
      CALL GLK_Print(700301)
      KeyDis=   302
      CALL Semalib_VVplot(KeyDis,"URHO2",-700302,' O(alf2) fin $',IDA+300)
      CALL GLK_Print(700302)
      CALL GLK_Operat(700302,'-',700301,800306,1D0,1D0)       
      CALL GLK_Print(800306)

*** TECHNICAL/PHYSICAL PRECISION chi(log_xmax)
      YMIN = -0.010D0
      YMAX =  0.015D0
      XMAG =   1
      YMAG = 0.01
      CALL PLTECH("CUMU",
     $"RHTOTF: sig(log(vmax))/Born, Phys/Tech prec., YMAG = 0.01$",
     $ IDA+302, IDA+306, 100302, 200306,YMIN,YMAX,XMAG,YMAG,BORN)
*** TECHNICAL/PHYSICAL PRECISION chi(xmax)
      YMIN = -0.010D0
      YMAX =  0.015D0
      XMAG =   1
      YMAG = 0.01
      CALL PLTECH("CUMU",
     $"RHTOTF: sig(vmax)/Born, Phys/Tech prec., YMAG = 0.01$",
     $ IDB+302, IDB+306, 500302, 600306,YMIN,YMAX,XMAG,YMAG,BORN)
*** TECHNICAL/PHYSICAL PRECISION rho(log_x)
      YMIN = -0.002D0
      YMAX =  0.0080D0
      XMAG = 10
      YMAG = 1
      CALL PLTECH("NB10",
     $"RHTOTF: dsig/dlogx [nb] Phys/Tech. prec. , XMAG = 10$",
     $ IDA+302, IDA+306, 700302, 800306,YMIN,YMAX,XMAG,YMAG,1d0)
      WRITE(6,*) ' RHTOTF ended '
      END     

      SUBROUTINE RHBETF        
C     *****************          
C---------------------------------------------------------------
C Final state !!!
C The distribution v*rho(v) as function of log(v)         
C Contributions from various beta's in O(alf1) and Oalf(2)           
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
              
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHBETF   =============='      
      WRITE(NOUT,*) ' ====================================='

      IDA = 20000
      IDD = 25000
      CALL RENHIF(IDA+310,IDD+310)
      CALL RENHIF(IDA+311,IDD+311)
      CALL RENHIF(IDA+320,IDD+320)
      CALL RENHIF(IDA+321,IDD+321)
      CALL RENHIF(IDA+322,IDD+322)
C BOOKS HISTO AND FILLS IT WITH A FUNCTION    
      KeyDis=   310             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2310,' O(alf1)       $',IDA+300)
      KeyDis=   311    
      CALL Semalib_VVplot(KeyDis,"URHO ",-2311,' O(alf1)       $',IDA+300)
      KeyDis=   320             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2320,' O(alf2)       $',IDA+300)
      KeyDis=   321             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2321,' O(alf2)       $',IDA+300)
      KeyDis=   322             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2322,' O(alf2)       $',IDA+300)
C---------------------------------------------------------------
C----------------------- O(alf1) -------------------------------
      CALL GLK_PLTITLE('RHBETF: beta-s,    O(alf1)    $')           
      CALL GLK_YMINIM( 0,-0.080D0)    
      CALL GLK_YMAXIM( 0, 0.130D0)    
      CALL GLK_PRINT(IDD+310)          
      CALL GLK_PRINT(IDD+311)          
C==== Monte Carlo               
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT( IDD+310,' ',' ',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT( IDD+311,'S',' ',0)              
C==== analitical                
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2310,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2311,'S','*',0)              
C---------------------------------------------------------------
C----------------------- O(alf2) -------------------------------
C THE BETA2 CONTR. MULTIPLIED BY FACTOR XMAG
      XMAG=5D0    
      CALL GLK_PLTITLE('RHBETF: beta-s,    O(alf2), beta2*5    $')
      CALL GLK_OPERAT( IDD+322,'+', IDD+322, IDD+322, 0D0,XMAG) 
      CALL GLK_OPERAT( 2322,   '+',    2322,    2322, 0D0,XMAG) 
      CALL GLK_PRINT(IDD+320)          
      CALL GLK_PRINT(IDD+321)          
      CALL GLK_PRINT(IDD+322)          
      CALL GLK_PRINT(2320)          
      CALL GLK_PRINT(2321)          
      CALL GLK_PRINT(2322)          
C==== PLOTTING                  
C Monte Carlo 
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(IDD+320,' ',' ',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(IDD+321,'S',' ',0) 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(IDD+322,'S',' ',0)              
C Analytical  
      CALL GLK_PLSET('DMOD',1D0)    
      CALL GLK_PLOT(2320,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2321,'S','*',0)              
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(2322,'S','*',0)              
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

      SUBROUTINE RHDIFF         
C     *****************          
C---------------------------------------------------------------
C Final state !!!!
C The distribution v*rho(v) as function of log(v)         
C Differences (M.C. - analyt.) for O(alf2), O(alf1)
C---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / INOUT  / NINP,NOUT  
      SAVE
              
      WRITE(NOUT,*) ' ====================================='      
      WRITE(NOUT,*) ' ==========    RHDIFF   =============='      
      WRITE(NOUT,*) ' ====================================='      
      IDA = 20000
      IDD = 25000
      CALL RENHIF(IDA+310,IDD+310)
      CALL RENHIF(IDA+311,IDD+311)
      CALL RENHIF(IDA+320,IDD+320)
      CALL RENHIF(IDA+321,IDD+321)
      CALL RENHIF(IDA+322,IDD+322)
C BOOKS HISTO AND FILLS IT WITH A FUNCTION    
      KeyDis=   310             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2310,' O(alf1)       $',IDA+300)
      KeyDis=   311             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2311,' O(alf1)       $',IDA+300)
      KeyDis=   320             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2320,' O(alf2)       $',IDA+300)
      KeyDis=   321             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2321,' O(alf2)       $',IDA+300)
      KeyDis=   322             
      CALL Semalib_VVplot(KeyDis,"URHO ",-2322,' O(alf2)       $',IDA+300)
C----------------------- O(alf1) -------------------------------
      CALL GLK_OPERAT( IDD+310,'-',2310,5310,1D0,1D0)   
      CALL GLK_OPERAT( IDD+311,'-',2311,5311,1D0,1D0)   
      CALL GLK_PRINT(5310)         
      CALL GLK_PRINT(5311)         
      CALL GLK_OPERAT( IDD+320,'-',2320,5320,1D0,1D0)   
      CALL GLK_OPERAT( IDD+321,'-',2321,5321,1D0,1D0)   
      CALL GLK_OPERAT( IDD+322,'-',2322,5322,1D0,1D0)   
      CALL GLK_PRINT(5320)         
      CALL GLK_PRINT(5321)         
      CALL GLK_PRINT(5322)         
      CALL GLK_YMINIM( 0,-0.010D0)    
      CALL GLK_YMAXIM( 0, 0.010D0)    
C----------------------- O(alf1) -------------------------------
      CALL GLK_PLTITLE('RHDIFF:  (M.C. - analyt.),  O(alf1)  $') 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(5310,' ',' ',0)              
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(5311,'S','*',0)              
C----------------------- O(alf2) -------------------------------
      CALL GLK_PLTITLE('RHDIFF:  (M.C. - analyt.),  O(alf2)  $') 
      CALL GLK_PLSET('DMOD',2D0)    
      CALL GLK_PLOT(5320,' ',' ',0)              
      CALL GLK_PLSET('DMOD',3D0)    
      CALL GLK_PLOT(5321,'S','*',0)              
      CALL GLK_PLSET('DMOD',4D0)    
      CALL GLK_PLOT(5322,'S','*',0)              
C----           
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

      SUBROUTINE renhif(id1,id2)
*     **************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.1415926535897932D0,ALFINV=137.03604D0)
      CHARACTER*80 TITLE
      SAVE

      CALL Semalib_GetCMSene(CMSene)
*
      IdGen = 6
      CALL GLK_WtMon(   1,idgen,par1,par2,par3 )
      nevt    =  par3
      xscrnb  =  par1
*
      gnanob=389.385d-30*1.d33
      sig0nb =  4d0*pi/(alfinv**2*3d0*cmsene**2)*gnanob
      XSECR  =  xscrnb/sig0nb
*
      CALL GLK_hinbo1(id1,title,nbt,tmin,tmax)
      fact = nbt*xsecr/(nevt*(tmax-tmin)*log(10.))
      fact = fact/XSECR   !?????
*
      CALL GLK_Operat(id1,'+',id1,id2, fact, 0d0)

      END
