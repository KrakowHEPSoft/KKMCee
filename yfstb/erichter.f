

C          ---------------- LOGBOOK -----------------

**** CO TRZEBA ZROBIC ABY PODLACZYC MOJ KIT
**** 1. wziasc SUSAMPZ0.FOR jako bibloteke, wszystkie potrzebne twoje
**** funkcje sa tam zdublowane a wiec po podlaczeniu mozesz robic ewentualne
**** dochodzenie czy sie roznia, w tym sensie bibloteka nie korzysta z YFS3
C
C->> wkopiowalem SUSAMPZ0.FORTRAN na koniec YFS3LIB
C->> kompilacja bez klopotow (dwa niegrozne smiecie w kodzie)
C
**** 2. jezeli zmieniales COMMON's to jedyne miejsce gdzie trzeba zajrzec
**** to subr.DIALOG gdzie potrzebne common;s z YFS3 sa tlumaczone na
**** wewnetrzne tej bibloteki
C
C->> nie zmienialem common-ow ale zajrzalem...
C
**** 3. uwaga!!! do tej subr DIALOG musi byc przekazana wart KEYZET, ja
**** dodalam
****    COMMON /KEYS/ KEYZET ktory musi sie wypelniac w YFS3
C
C->> hm... dlaczego Ela nie wziela caly common /KEYYFS/ ?? 
C->> jeszce znalazlem taki komentarz:
*** as compare to previous version one need extra common /KYES/ to pass  
*** information of KEYZET from SPLOT to DIALOG. As a test you can compare...
C->> wyglada na to ze /KEYS/ jest niepotrzebny, wystepuje tylko w SPLOT
C->> i w DIALOG i nigdzie indziej (parallel spurious communication), 
C->> wstawiam /KEYYFS/ do subr. DIALOG i zobaczymy czy sie uda.
C
**** 4. ta bibloteka  jest
**** wolana przez CALL AMPLIT(MODE) wg standartowej konwencji
C
C->> gdzie jest aktualnie CALL AMPLIT ? krotkie poszukiwania: w ROBOL4.
C
**** 5. output w postaci monitorowania ilorazu jest wypisywany dla MODE>0
**** istotny common ktory ona wypelnia to
****      COMMON / WGTD2  / WTHARD,WT2,AWRWT2
C
C->> nie bardzo rozumiem ale zdaje sie trzeba przeczytac uwagi ponizej 
C
**** 
**** CO DALEJ............
****  1. wziasc SUSMAINZ0.for
****  2. wyjac z niego subr. ROBOL4(MODE)
****  3. wolac subr ROBOL4(MODE) wg zwyklej konwencji, musi byc
****  wolana po  EXPAND(MODE)
****  4. wyjac subr. ORDPTF, ORDPTC  ktore dodatkowo ona uzywa
C
C->> Aha! chodzi o to ze /WGTD2/ komunikuje AMPLIT i ROBOL4, tak jest. 
C->> Przekladam ROBOL4 do tbmain.f (ale heca! KAT4 juz byl)
C->> ORDPTC nie znalazlem, nie ma go w susmainz0.fortran wcale.
C->> pewno chodzi o ORDPCF. 
C->> No to kompilujemy! linkujemy! nic nie kopie!
C->> No to wykonujemy 100 przyp. Crash! poprawiam argumenty HOPERA,
C->> 200 przyp. idzie! 1000 przyp. idzie! 100 przyp. CRASH!
C->> Jakies niestabilnosci wziaz sa ale 100tys przyp. (70min.) przeszlo!
C
**** ABY ZROBIC RYSUNKI.....
****   1. wziasc ostatni SUSFIGZ0.for
****   2. wyjac subr ROB4FIG ktory sluzy do robienia rysunkow
C
C->> Przelozylem ROB4FIG do mojego ostatnego SUSFIG-a, przerobilem
C->> na double precision. Zpuscilem i poszlo! Dla 100k bledy na 
C->> caly rysunek. To by bylo na tyle!
C
* UWAGA!!! poniewaz brak mi konsekwenci i nie wiedzialam czy bedziesz
*  zmienial te bibloteki wersje ktore nie maja d. na koncu sa do tej
* nie pelnej dbl.prec. GLIBD i LPLOT.
* 
* Wyprodukowalam te rusunki dla 10**6 przypadkow z Z0 i na piku
* (dane susz0) . Iloraz energi I i II fotonu, p_T/E_beam I i II fotonu
* p_C/E_phot I i II fotonu ( to ma byc miara kolinearnosci) Prawidlowe
* granice to iloraz =1 dla E/E_beam=0,  p_T/E_beam = 0 , p_C/E_phot =1.
* Na odmiane wyslalam zbior .tex chociaz linia zdajsie przeklamuje.
* Zwroc uwage ze na piku z0 ten iloraz rozjezdza sie od jedynki troche
* inaczej niz poza pikiem i czyste QED (ela bledy sa za duze aby cos
* zdecydowanie powiedziec) ale ciagle jest dobrze!!!
C
C-> Ostatecznie przerzucilem wszystkie podprogramy Eli do glownego
C-> programu tbmain.f reszta jest tak jak w zwyklym programie yfs3
C-> Rysunki sa robione w tbfig.f ktory korzysta z glibl.f

C !!!!!!! Corrections since 24 Nov. 91 !!!!!
C
C  zbw/stj: SECTION=DREAL(SUM) marked with ==>>




C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C --------  01.XI.91  E. R-W        
C in this version initial bremstrahlung is kicked out and Z0 is properly
C as compare to previous version one need extra common /KYES/ to pass   
C information of KEYZET from SPLOT to DIALOG. As a test you can compare 
C previous and present version for Z0 switched off (they are equivalent)
C the procedure of switching on Z0 effect have been carefully tested on 
C the Born amlitude and compared with others programs.This two test     
C make the result "wiarygodny".     
C-- note important correction in BORNVV which concern Z0 propagator     
C-- the part proportional to FF1 have been omitted in previous version  
C --------  01.XI.91  E. R-W        
C---------------------------------------------------------------------- 
C---DIST2 - SPIN AMPLITUDE(E.R.-W.) AND COMPACT FORMULA(S.J.)-----------
C--ADDITIONAL BORNV FUNCTION WITH MASS CORRECTION---------------------- 
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE AMPLIT(MODE)       
C     ********************************      
C ONLY DOUBLE BREMSTRAHLUNG-PURE INITIAL AND PURE FINAL 
C THIS ROUTINE TEST PRECISION OF CROSS SECTION IN COMPACT FORMULA       
C AND SPIN AMPLITUDE FORMULA        
C Z0 STILL NOT PROPERLY INCLUDED    
C E.RICHTER-WAS NOVEMBER 1989       
C     **********************************    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON / MOMSET / QF1(4),QF2(4),SPHUM(4),SPHOT(100,4),NPHOT       
      COMMON / MOMINI / XF1(4),XF2(4),XPHUM(4),XPHOT(100,4),NPHOX       
      COMMON / MOMFIN / YF1(4),YF2(4),YPHUM(4),YPHOT(100,4),NPHOY       
      COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF     
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(100)  
      COMMON / INOUT  / NINP,NOUT     
      COMMON / WGTD2  / WTHARD,WT2,AWRWT2,COMPFI,AMPLFI   
      COMMON / BXFMTS / BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      SAVE / MOMSET /,/ MOMINI /,/ MOMFIN /,/ WEKING /
      SAVE / WGTALL /,/ INOUT  /,/ WGTD2  /,/ BXFMTS /
C  ...........  
        
      IF(MODE.EQ.-1) THEN   
C     ==================================================================
C     =====================INITIALIZATION===============================
C     ==================================================================
      CALL WMONIT(-1,41,DUMM1,DUMM2,DUMM3)  
      CALL WMONIT(-1,42,DUMM1,DUMM2,DUMM3)  
C...comunication between to cooled programs 
      CALL DIALOG(-1)       
        
      ELSEIF(MODE.EQ.0) THEN        
C     ==================================================================
C     ====DOUBLE BREMSTR. SPIN AMPLITUDE AND COMPACT FORMULA ===========
C     ==================================================================
        
      WTHARD=0D0    
C.....dealing with momenta  
      CALL DIALOG( 0 )      
C ...Return for events outside phase space  
      IF(WTCRU1*WTCRU2.EQ.0D0) GOTO 770     
C ...Return for events with NPHOT=(NPHOY+NPHOX) < >  2  
      IF((NPHOY+NPHOX).NE.2) GOTO 770       
CC...Return for events with no 2  hard photons  
C ...Return for events with no 1  hard photons  
      EHARD  =  0.1D0*ENE   
cc      IF(SPHOT(1,4).LT.EHARD.AND.SPHOT(2,4).LT.EHARD)  GOTO 770 
cc      IF(SPHOT(1,4).LT.EHARD. OR.SPHOT(2,4).LT.EHARD)  GOTO 770 
C.........................................  
      IF (NPHOX.EQ.0.AND.NPHOY.EQ.2) THEN   
C.. double bremstrahlung pure final state TWO HARDS PHOTONS     
C.. cross section compact formula from YFS version january 1990 
      WTHARD=1D0    
      CALL DFINAPR(COMPFI)  
C.. cross section spin amplitude formula E.W.   november 1989   
      CALL DUBLFIN(AMPLFI)  
      WT2 = COMPFI/AMPLFI     
      CALL WMONIT(0,42,WT2 ,1.D0 ,5D-4)     
      ENDIF     
*%%%%%%%%%%%%%%%%
c      write(6,*) "COMPFI,AMPLFI ",COMPFI,AMPLFI
*%%%%%%%%%%%%%%%%
  770 CONTINUE  
      ELSE      
C     ==================================================================
C     =====================FINAL WEIGHT REPORT==========================
C     ==================================================================
      CALL WMONIT(1,42,AWRWT2,ERR2 ,DUMM3)  
CCC      CALL WMONIT(2,42,AWRWT2,ERR2 ,DUMM3)   
        
C .....................Tests on double bremstrahlung formulas...........
C ......................................................................
        
      WRITE(NOUT,BXTXT) '*********************************'     
      WRITE(NOUT,BXTXT) '  TEST ON DOUBLE BREMSTRAHLUNG.    '   
      WRITE(NOUT,BXL2F) AWRWT2,AWRWT2*ERR2, 
     $                  'D2COM/D2AMPLI FIN  ','WT2   ','  ' 
      WRITE(NOUT,BXTXT) '*********************************'     
      ENDIF     
C     =====     
      END       
      SUBROUTINE DIALOG(MODE)       
C     ********************************      
C this routine transfers kinematics and parameters beetwen YFS3 
C and added spin amplitudes routines.This is a consequence of lack of   
C compactibility between this two packeges  
C     ********************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
c.....commons from YFS3     
      COMMON / MOMINI / XF1(4),XF2(4),XPHUM(4),XPHOT(100,4),NPHOX       
      COMMON / MOMFIN / QF1(4),QF2(4),YPHUM(4),YPHOT(100,4),NPHOY       
      COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF     
      COMMON / INOUT  / NINP,NOUT     
C->>> stj: wymieniam common /keys/ na standardowy /keyfs/
cc    COMMON /KEYS    / KEYZET  
      COMMON / KEYYFS / KEYZET,KEYBRM,KEYFIX,KEYRED,KEYWGT  
c.....commons from added SPIN AMPLITUD      
      COMMON / FOTON  / ARBIT(4)      
      COMMON / BHPAR2 / CMSENE,AAMEL        
      COMMON / BHPAR1 / CMS,AMMI    
      COMMON / WEKINP / AMAZZ,GAMMZZ,SINW22,IDEE  
      SAVE   / WEKINP /  
      COMMON / COEFF  / V,A 
      COMMON / MOMS   / TP1(4),TQ1(4),TP2(4),TQ2(4),PH1(4),PH2(4)       
      SAVE / MOMINI /,/ MOMFIN /,/ WEKING /,/ INOUT  /
      SAVE / KEYYFS /,/ FOTON/,/ BHPAR2 /,/ BHPAR1 /,/ COEFF/,/ MOMS /
      SAVE / WEKINP /
      DIMENSION AA(4),P1(4),P2(4)   
      SAVE AA,P1,P2
        
      IF(MODE.EQ.-1) THEN   
C     ==================================================================
C     =====================INITIALIZATION===============================
C     ==================================================================
        
c.....  
      CMSENE=2D0*ENE        
C... define beam momenta    
      CALL GIBEA(CMSENE,AMEL,P1,P2) 
C... fixing arbitrary gauge vector in common /FOTON/    
      AA(1)=0.0D0   
      AA(2)=0.0D0   
      AA(3)=0.5D0*CMSENE    
      AA(4)=0.5D0*CMSENE    
      ARBIT(4)=AA(4)        
      ARBIT(1)=AA(1)        
      ARBIT(3)=AA(3)*COS(0.5)-AA(2)*SIN(0.5)    
      ARBIT(2)=AA(2)*COS(0.5)+AA(3)*SIN(0.5)    
c.. masses and Zo couplings 
      AMAZZ  = AMAZ
      GAMMZZ = GAMMZ
      SINW22 = SINW2
      IDEE   = IDE
      AAMEL =AMEL    
      AMMI  =AMFIN   
      V=0D0     
      A=0D0     
c....Zo swichted off        
      IF(KEYZET.EQ.1) THEN  
        AV= 4D0*SQRT(SINW2*(1D0-SINW2))       
        V= (-1D0+4*SINW2)/AV  
        A= 1D0/AV     
      ENDIF 
      write(6,*) "     a,v=", a,v
      ELSEIF(MODE.EQ.0) THEN        
C     ======================        
C... in the spin amplit.packege fermions momenta are TP1,TQ1,TP2,TQ2    
C...this transform convension from YFS2 to BHLUMI in which spin 
C...amplitude formulas were written 
c.................................................................      
      IF (NPHOX.EQ.0.AND.NPHOY.EQ.2) THEN   
      DO 20 I=1,4   
      TP2(I)=QF1(I) 
      TQ2(I)=QF2(I) 
      TP1(I)=P1(I)  
      TQ1(I)=P2(I)  
      PH1(I) =YPHOT(1,I)    
 20   PH2(I) =YPHOT(2,I)    
      ENDIF     
        
      ENDIF     
C     =========     
      END       

C.. spin amplitudes subroutines E.W. november 1989      
C...................................................    
C.. PURE FINAL   STATE BREMSTRAHLUNG        
      SUBROUTINE DUBLFIN(SECTION)   
*     ***********************************   
      IMPLICIT REAL*8(A-H,O-Z)      
      PARAMETER(PI=3.1415926535897932D0,ALFINV= 137.03604D0)        
      COMPLEX*16 XMDFIN,S,SUM,C     
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / BHPAR2 / CMSENE,AMEL 
c it seems that NIC is not used !!!!!!!!!
      COMMON /NIC/ XNORM    
      SAVE / MOMS /,/ BHPAR2 /,/NIC/
        
      ALFA=1D0/ALFINV       
      S0=CMSENE**2  
ccc      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2     
ccc     $   -(P2(3)+Q2(3))**2  
      C=DCMPLX(1D0,0D0)     
      SUM=0D0*C     
!######      XNORM=ALFA**4/PI**4/S0/16D0   
      XNORM=4 
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS  
      DO 10 JL=1,3,2        
      LAM1=2-JL     
      LAM2=2-JL     
      DO 10 J=1,3,2 
      LAM3=2-J  
      DO 10 I=1,3,2 
      LAM4=2-I  
      DO 10 M=1,3,2 
      LEPS1=2-M     
      DO 10 N=1,3,2 
      LEPS2=2-N     
      S= XMDFIN( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2)       
      S=S*DCONJG(S)*XNORM   
c     IF((LAM1*LAM2.GT.0D0)) THEN   
c     SS=DBLE(S)    
c     PRINT *,SS,LAM1,LAM2,LEPS1,LEPS2      
c     ENDIF     
      SUM=SUM+S     
        
  10  CONTINUE  
*CROSS SECTION  
C==>>   SECTION=DBLE(SUM)     
      SECTION=DREAL(SUM)     
      END       

***********************************
** SPIN AMPLITUDE  FINAL  DOUBLE **
***********************************      
      FUNCTION XMDFIN(L1,L2,L3,L4,LE1,LE2)  
*     *******************************       
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX*16  XMDFIN,XP2P2,XP2Q2,XQ2Q2  
      COMPLEX *16 CR1,CR2,S,CPH1,CPH2,C,C1  
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)  
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / CURR / PH1(4),PH2(4),PH(4)   
      COMMON /BHPAR2/ CMSENE, AMEL  
      COMMON /FOTON/ ARBIT(4)       
      SAVE /HELP2/,/ MOMS /,/ CURR /,/BHPAR2/,/FOTON/
        
      S0= CMSENE**2 
      DO 15 I=1,4   
      PH(I)=PK1(I)+PK2(I)   
 15   ZER(I)=0D0    
      C=DCMPLX(1D0,0D0)     
      C1=C*DELTA(L1,L2)     
      S=(0D0,0D0)   
      DO 10 KK=1,2  
      IF (KK.EQ.1) THEN     
      DO 5 II=1,4   
      PH1(II)=PK1(II)       
  5   PH2(II)=PK2(II)       
      LL1=-LE1  
      LL2=-LE2  
      ELSE      
      DO 6 II=1,4   
      PH1(II)=PK2(II)       
  6   PH2(II)=PK1(II)       
      LL1=-LE2  
      LL2=-LE1  
      ENDIF     
      CPH1=DSQRT( PROPFIN(PH1,ARBIT))*C     
      CPH2=DSQRT( PROPFIN(PH2,ARBIT))*C     
*--------------------------------------------------------       
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE P1       
      CR1=PROPFIN(P2,PH1)*C 
      CR2=PROPFIN(P2,PH)*C  
      CALL FUBLEA(L1,L2,L3,L4,LL1,LL2,P2,PH1,CR1,P2,PH1,PH2,    
     $   CR2,XP2P2) 
*--------------------------------------------------------       
*SPIN AMPLITUDE - ONE BREMS FROM LINE P1,ONE FROM LINE Q1       
      CR1=PROPFIN(P2,PH1)*C 
      CR2=PROPFIN(Q2,PH2)*C 
      CALL FUBLEC(L1,L2,L3,L4,LL1,LL2,P2,PH1,CR1,Q2,PH2,CR2,XP2Q2)      
*--------------------------------------------------------       
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE Q1       
      CR1=PROPFIN(Q2,PH)*C  
      CR2=PROPFIN(Q2,PH2)*C 
      CALL FUBLEB(L1,L2,L3,L4,LL1,LL2,PH1,PH2,Q2,CR1,PH2,Q2,    
     %    CR2,XQ2Q2)        
*----------------------==--------------------------------       
*TOTAL SPIN AMPLITUDE       
      S=S+(XP2P2-XP2Q2+XQ2Q2)       
 10   CONTINUE  
      XMDFIN=S  
        
      END       

      SUBROUTINE FUBLEA(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,        
     %  R2A,R2B,R2C,CR2,X)  
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   CR1,CR2,X       
      COMPLEX *16 C1,C,CPH1,CPH2,CZ0P,CZ0M,HIS  
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SC     
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / CURR / PH1(4),PH2(4),PH(4)   
      COMMON /FOTON/ ARBIT(4)       
      COMMON / BHPAR2 / CMSENE,AMEL 
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)  
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),        
     $  SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2), 
     $  SS(2,2),SSA(2,2),SSB(2,2)   
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4)      
      SAVE / MOMS /,/ CURR /,/FOTON/,/ BHPAR2 /,/HELP2/

      S0= CMSENE**2 
*P1 LINE        
      CALL MULTI( L3,P2 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)     
*FIRST PROPAGATOR   
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1A,R1A,C , C,SPH1,CR1)      
      CALL MULTI(  1, R1A, R1A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)    
      CALL ILOCZ(SPH1,SP ,SSA)      
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1, CR1)     
      CALL MULTI(  1, R1B, R1B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)    
      CALL ILOCZ(SPH1,SP ,SSB)      
      CALL DODAJ(SSA,SSB,SS)        
*SECOND PROPAGATOR  
*Z0 coefficients    
      CZ0P=HIS(1,L1,S0)     
      CZ0M=HIS(-1,L1,S0)    
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2A,R2A,C , C,SPH2, CR2)      
      CALL MULTI(  1, R2A, R2A,CZ0P,CZ0M, L1,P1 ,Q1 ,C  , C,SFI,C1)     
      CALL ILOCZ(SPH2,SFI,SA)       
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2B,R2B,C , C,SPH2, CR2)      
      CALL MULTI(  1, R2B, R2B,CZ0P, CZ0M, L1,P1 ,Q1 ,C  , C,SFI,C1)    
      CALL ILOCZ(SPH2,SFI,SB)       
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2C,R2C,C , C,SPH2, CR2)      
      CALL MULTI(  1, R2C, R2C,CZ0P,CZ0M, L1,P1 ,Q1 ,C  , C,SFI,C1)     
      CALL ILOCZ(SPH2,SFI,SC)       
      CALL DODAJ(SA,SB,S)   
      CALL DODAJ(S,SC,S)    
*Q1 LINE        
      CALL MULTI( L2,Q1 ,P1 ,C, C, L4,Q2 ,ZER,C  , Z,SPO ,C  )  
      CALL ILOCZ(S,SPO,S)   
      CALL ADD3(SEL ,SS,S,X)        
        
      END       
      SUBROUTINE FUBLEC(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,R2A,R2B,CR2,X)  
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   CR1,CR2,X       
      COMPLEX *16 C1,C,CPH1,CPH2,CZ0P,CZ0M,HIS  
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB        
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / CURR / PH1(4),PH2(4),PH(4)   
      COMMON / BHPAR2 / CMSENE,AMEL 
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)  
      COMMON /FOTON/ ARBIT(4)       
      SAVE / MOMS /,/ CURR /,/ BHPAR2 /,/HELP2/,/FOTON/
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),        
     $  SP(2,2),SA(2,2),SB(2,2),S(2,2), 
     $  SS(2,2),SSA(2,2),SSB(2,2)   
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4) 

      S0= CMSENE**2 
*P1 LINE        
      CALL MULTI( L3,P2 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)     
*FIRST PROPAGATOR   
*Z0 coefficients    
      CZ0P=HIS(1,L1,S0)     
      CZ0M=HIS(-1,L1,S0)    
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1 , R1A,R1A,C , C,SPH1,CR1)       
      CALL MULTI(  1, R1A, R1A,CZ0P,CZ0M,L1 ,P1 ,Q1,C , C,SFI,C1)       
      CALL ILOCZ(SPH1,SFI,SSA)      
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1, CR1)     
      CALL MULTI(  1, R1B,R1B,CZ0P,CZ0M,L1 ,P1 ,Q1,C , C,SFI,C1)        
      CALL ILOCZ(SPH1,SFI,SSB)      
      CALL DODAJ(SSA,SSB,SS)        
*SECOND PROPAGATOR  
      CALL MULTI( L2 ,Q1 ,P1,C ,C ,1  , R2A,R2A,C , C,SP  , CR2)        
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)     
      CALL ILOCZ(SP  ,SPH2,SA)      
      CALL MULTI( L2 ,Q1 ,P1,C ,C, 1  , R2B,R2B,C , C,SP  , CR2)        
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)     
      CALL ILOCZ(SP  ,SPH2,SB)      
      CALL DODAJ(SA,SB,S)   
*Q1 LINE        
      CALL MULTI(LE2,PH2,ARBIT,C ,C,L4,Q2 ,ZER,C  , C,SPO,C )   
      CALL ADD4(SEL ,SS,S,SPO,X)    
      END       

      SUBROUTINE FUBLEB(L1,L2,L3,L4,LE1,LE2,R1A,R1B,R1C,CR1,    
     % R2A,R2B,CR2,X)       
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   CR1,CR2,X       
      COMPLEX *16 C1,C,CPH1,CPH2,CZ0P,CZ0M,HIS  
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SSC    
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / CURR / PH1(4),PH2(4),PH(4)   
      COMMON / BHPAR2 / CMSENE,AMEL 
      COMMON /FOTON/ ARBIT(4)       
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)  
      SAVE / MOMS /,/ CURR /,/ BHPAR2 /,/FOTON/,/HELP2/
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),        
     $  SP(2,2),SA(2,2),SB(2,2),S(2,2), 
     $  SS(2,2),SSA(2,2),SSB(2,2),SSC(2,2)      
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R1C(4)      

      S0= CMSENE**2 
*P1 LIN2        
*Z0 coefficients    
      CZ0P=HIS(L3,L1,S0)    
      CZ0M=HIS(-L3,L1,S0)   
      CALL MULTI( L3,P2 ,ZER,CZ0P,CZ0M,L1 ,P1 ,Q1,C , C ,SEL,C1)        
*FIRST PROPAGATOR   
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1A,R1A ,C , C ,SFI ,-CR1)        
      CALL MULTI(  1, R1A, R1A,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C ) 
      CALL ILOCZ(SFI ,SPH1,SSA)     
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1B ,R1B ,C , C ,SFI ,-CR1)       
      CALL MULTI(  1, R1B, R1B,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C ) 
      CALL ILOCZ(SFI ,SPH1,SSB)     
      CALL DODAJ(SSA,SSB,SS)        
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1C ,R1C ,C , C ,SFI ,-CR1)       
      CALL MULTI(  1, R1C, R1C,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C ) 
      CALL ILOCZ(SFI ,SPH1,SSC)     
      CALL DODAJ(SS,SSC,SS) 
*SECOND PROPAGATOR  
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2A,R2A,C , C,SP  ,-CR2)     
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)  
      CALL ILOCZ(SP  ,SPH2,SA)      
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2B,R2B,C , C,SP  ,-CR2)     
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)  
      CALL ILOCZ(SP  ,SPH2,SB)      
      CALL DODAJ(SA,SB,S)   
*Q1 LINE        
      CALL MULTI(LE2, PH2,ARBIT,C, C, L4,Q2 ,ZER,C  , C,SPO,C ) 
      CALL ADD4(SEL ,SS,S,SPO,X)    
      END       

*PROPAGATOR IN INITIAL STATE        
      FUNCTION PROPIN( P1,PH)       
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      DIMENSION P1(4),PH(4) 
        
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2  
      P1PH=P1(4)*PH(4)-P1(1)*PH(1)-P1(2)*PH(2)-P1(3)*PH(3)      
      PROPIN=1D0/(-2D0*P1PH+PHPH)   
      END       

*PROPAGATOR-FINAL STATE     
      FUNCTION PROPFIN(P1,PH)       
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      DIMENSION P1(4),PH(4) 
        
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2  
      P1PH=P1(4)*PH(4)-P1(1)*PH(1)-P1(2)*PH(2)-P1(3)*PH(3)      
      PROPFIN=1D0/( 2D0*P1PH+PHPH)  
      END       

      SUBROUTINE  DODAJ(SEL,SPO,S)  
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   SEL,SPO,S       
      DIMENSION SEL(2,2),SPO(2,2),S(2,2)    
        
* ADDING   MATRIX S =SEL+SPO        
      S(1,1)=SEL(1,1)+SPO(1,1)      
      S(1,2)=SEL(1,2)+SPO(1,2)      
      S(2,1)=SEL(2,1)+SPO(2,1)      
      S(2,2)=SEL(2,2)+SPO(2,2)      
      END       

      SUBROUTINE ILOCZ(SEL,SPO,S)   
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   SEL,SPO,S       
      DIMENSION SEL(2,2),SPO(2,2),S(2,2)    
        
* MULTIPLE MATRIX S =SEL*SPO        
      S(1,1)=SEL(1,1)*SPO(1,1)+SEL(1,2)*SPO(2,1)        
      S(1,2)=SEL(1,1)*SPO(1,2)+SEL(1,2)*SPO(2,2)        
      S(2,1)=SEL(2,1)*SPO(1,1)+SEL(2,2)*SPO(2,1)        
      S(2,2)=SEL(2,1)*SPO(1,2)+SEL(2,2)*SPO(2,2)        
      END       

      SUBROUTINE ADD4(SEL,SPH ,SFI,SPO ,X)  
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   SEL,S,SS,X,SPH,SFI ,SPO 
      DIMENSION SEL(2,2),S(2,2),SPO(2,2),SFI(2,2),SS(2,2),SPH(2,2)      
        
* MULTIPLE MATRIX S =SEL*SPH        
      CALL ILOCZ(SEL,SPH,S) 
* MULTIPLE MATRIX SS=S *SFI 
      CALL ILOCZ(S,SFI,SS)  
* MULTIPLE MATRIX S=SS*SPO  
      CALL ILOCZ(SS,SPO ,S) 
*CONSTRACTION TO THE SCALAR OBJECT  
      X=S (1,1)+S (1,2)+S (2,1)+S (2,2)     
        
      END       

      SUBROUTINE ADD3(SEL,SPH ,SFI,X)       
*     **************************************************        
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16   SEL,S,SS,X,SPH,SFI      
      DIMENSION SEL(2,2),S(2,2),SFI(2,2),SS(2,2),SPH(2,2)       
        
* MULTIPLE MATRIX S =SEL*SPH        
      CALL ILOCZ(SEL,SPH,S) 
* MULTIPLE MATRIX SS=S *SFI 
      CALL ILOCZ(S,SFI,SS)  
*CONSTRACTION TO THE SCALAR OBJECT  
      X=SS(1,1)+SS(1,2)+SS(2,1)+SS(2,2)     
      END       

*PRODUCT OF VB(L1,P1,Q1,A1,B1*V(L2,P2,Q2,A2,B2) 
      SUBROUTINE MULTI(L1,P1,Q1,A1,B1,L2,P2,Q2,A2,B2,SS,C)      
*     *******************************       
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX*16 A1,A2,SS,C,B1,B2,SSPM      
      DIMENSION P1(4),P2(4),Q1(4),Q2(4),SS(2,2) 
        
      A1=DCONJG(A1) 
      SS(1,1)=A1*A2*SSPM( L1, L2,P1,P2)*C   
      SS(1,2)=A1*B2*SSPM( L1,-L2,P1,Q2)*C   
      SS(2,1)=B1*A2*SSPM(-L1, L2,Q1,P2)*C   
      SS(2,2)=B1*B2*SSPM(-L1,-L2,Q1,Q2)*C   
      END       

*FUNCTION DIRAC DELTA FOR INTEGER ARGUMENTS 
      FUNCTION DELTA(L1,L2) 
*    ****************************************   
      IMPLICIT REAL*8(A-H,O-Z)      
        
      N=(1+L1*L2)   
      IF (N.EQ.2) DELTA=1D0 
        IF (N.EQ.0) DELTA=0D0       
      END       

*DEFINE LIGHT=LIKE MOMENTA FOR MASIVE FERMIONS  
      SUBROUTINE MOMENTA(P,P1,P2)   
*     *******************************       
      IMPLICIT REAL*8(A-H,O-Z)      
      DIMENSION P(4),P1(4),P2(4)    
      COMMON / BHPAR2 / CMSENE,AMEL 
      SAVE   / BHPAR2 /
        
      PM=DSQRT(P(1)**2+P(2)**2+P(3)**2)     
      HI1=1/2D0*(P(4)+PM)   
      HI2=1/2D0*(P(4)-PM)   
      DO 10 I=1,3   
      P1(I)= HI1*P(I)/PM    
 10   P2(I)=-HI2*P(I)/PM    
      P1(4)=HI1     
      P2(4)=HI2     
      END       

*FUNCTION OF SPINOR PRODUCT OF MASSIVE PARTICLES        
      FUNCTION SSPM(L1,L2,P1,P2)    
*     ********************  
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX*16 SPLUS,SMINS,SSPM   
      DIMENSION P1(4),P2(4) 
        
      IF (P1(4).LT.1D-15) THEN      
      SSPM=(0D0,0D0)        
      ELSE IF (P2(4).LT.1D-15) THEN 
      SSPM=(0D0,0D0)        
      ELSE      
      R1=P1(4)-P1(3)        
      R2=P2(4)-P2(3)        
      AM1=(P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2) 
      AM2=(P2(4)**2-P2(1)**2-P2(2)**2-P2(3)**2) 
      IF (AM1.GT.0D0) THEN  
      AM1=DSQRT(AM1)        
      IF (ABS(R1).LT.1D-3) AM1=-AM1 
      ELSE      
      AM1=0D0   
      ENDIF     
      IF (AM2.GT.0D0) THEN  
      AM2=DSQRT(AM2)        
      IF (ABS(R2).LT.1D-3) AM2=-AM2 
      ELSE      
      AM2=0D0   
      ENDIF     
      ETA2=DSQRT(P2(4)-P2(1))       
      ETA1=DSQRT(P1(4)-P1(1))       
      SSPM=DELTA(L1, 1)*DELTA(L2,-1)*SPLUS(P1,P2)       
     $   +DELTA(L1,-1)*DELTA(L2, 1)*SMINS(P1,P2)        
     $   +(DELTA(L1, 1)*DELTA(L2, 1)+DELTA(L1,-1)*DELTA(L2,-1)) 
     $   *(AM1*ETA2/ETA1+AM2*ETA1/ETA2)     
      ENDIF     
      END       

*SPINOR PRODUCT     
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)     
      FUNCTION SPLUS(P1,P2) 
*     ********************  
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX*16 SPLUS,X,Y  
      DIMENSION P1(4),P2(4) 
        
      X=DCMPLX(P1(2),P1(3)) 
      Y=DCMPLX(P2(2),P2(3)) 
      SPLUS=X*DSQRT(P2(4)-P2(1))/DSQRT(P1(4)-P1(1))     
     $     -Y*DSQRT(P1(4)-P1(1))/DSQRT(P2(4)-P2(1))     
      END       

*SPINOR PRODUCT     
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)     
      FUNCTION SMINS(P1,P2) 
*     ********************  
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX*16 SPLUS,SMINS        
      DIMENSION P1(4),P2(4) 
        
      SMINS=-DCONJG(SPLUS(P1,P2))   
      END       

*FROM BOHM,DENNIER,HOLLIK,NUCL.PHYS.B304(1988),687      
      FUNCTION HIS(NI1,NI2,X)       
*     ***************************   
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16 Z,HIS     
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE    
      COMMON / COEFF  / V,A    
      SAVE / WEKINP /,/COEFF/
        
      BZ0RE=(X-AMAZ**2)     
     $      /((X-AMAZ**2)**2+GAMMZ**2*AMAZ**2)  
      BZ0IM=-GAMMZ*AMAZ     
     $      /((X-AMAZ**2)**2+GAMMZ**2*AMAZ**2)  
      Z=DCMPLX(BZ0RE,BZ0IM) 
      HIS=(1D0/X+(V+NI1*A)*(V+NI2*A)*Z)     
      END       

!===================================================================
!===================================================================
!===================================================================
!===================================================================
!===================================================================


C.. compact formulas subroutines S.J. version january 1990      
C...................................................    
*COMPACT FORM FOR DOUBLE INI BREMSTR CROSS SECTION      
      SUBROUTINE DINIAPR(SECTION)   
C     ******************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.1415926535897932D0, ALFINV=137.03604D0)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / BHPAR2 / CMSENE,AMEL 
      SAVE / MOMS /,/ BHPAR2 /
      DIMENSION QQ(4)       
        
      ALFA=1D0/ALFINV       
      S0=CMSENE**2  
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2     
     $   -(P2(3)+Q2(3))**2  
      DO 10 I=1,4   
 10   QQ(I)=P1(I)+Q1(I)     
      CALL ENDIST2(QQ,P1,Q1,P2,Q2,PK1,PK2,DIST2)        
      SECTION=(ALFA/4D0/PI**2)**2*ALFA**2/4D0/S0*DIST2*S0/S1    
      END       

      SUBROUTINE ENDIST2(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)  
C     ***********************************************   
C Provides double bremsstrahlung distribution - INITIAL state brem.     
C INPUT:  P1,P2,Q1,Q2,PH1,PH2, four momenta 
C OUTPUT: DIST2     double bremsstrahlung distribution  
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)     
      DIMENSION PR1(4),PR2(4),PH1R(4),PH2R(4),QR1(4),QR2(4)     
        
      CALL EREDUZ2(QQ,P1,P2,PH1,PH2,PR1,PR2,PH1R,PH2R)  
      CALL EREDUZ0(QQ,Q1,Q2,QR1,QR2)        
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2       
C infrared factors from reduced momenta     
C double bremsstrahlung Xsect in next-to-leading log approx.    
      CALL EGSOFA2(P1,P2,PH1,PH2,GF1,GF2)   
      CALL EGTHET1(PR1,PR2,QR1,COSTH1,COSTH2)   
      ANDI11= BORNVV(SVAR1,COSTH1)  
      ANDI12= BORNVV(SVAR1,COSTH2)  
      DIST2 =   GF1*ANDI11+   GF2*ANDI12    
      END       

      SUBROUTINE EGSOFA2(P1,P2,PH1,PH2,F1,F2)   
C     **************************************    
C NEW VERSION BY ELA WAS    
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
      PP = P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)      
      AM2= P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2  
      AM = AM2/(2D0*PP)     
      B1 = (P2(4)*PH1(4)-P2(1)*PH1(1)-P2(2)*PH1(2)-P2(3)*PH1(3))/PP     
      A1 = (P1(4)*PH1(4)-P1(1)*PH1(1)-P1(2)*PH1(2)-P1(3)*PH1(3))/PP     
      B2 = (P2(4)*PH2(4)-P2(1)*PH2(1)-P2(2)*PH2(2)-P2(3)*PH2(3))/PP     
      A2 = (P1(4)*PH2(4)-P1(1)*PH2(1)-P1(2)*PH2(2)-P1(3)*PH2(3))/PP     
      SFAC1  =  2D0/(PP*A1*B1)*WWM(A1,B1)   
      SFAC2  =  2D0/(PP*A2*B2)*WWM(A2,B2)   
      A1P= A1/(1D0-A2)      
      B1P= B1/(1D0-B2)      
      A2P= A2/(1D0-A1)      
      B2P= B2/(1D0-B1)      
      IF((A1+B1).GT.(A2+B2)) THEN   
        X1=WM (A1   )*WMS(A2P,B2P) +WM (A1P    )*WMS(A2,B2)     
        X2=WM (   B1)*WMS(A2P,B2P) +WM (    B1P)*WMS(A2,B2)     
      ELSE      
        X1=WM (A2   )*WMS(A1P,B1P) +WM (A2P    )*WMS(A1,B1)     
        X2=WM (   B2)*WMS(A1P,B1P) +WM (    B2P)*WMS(A1,B1)     
      ENDIF     
      F1 = X1*SFAC1*SFAC2/8D0       
      F2 = X2*SFAC1*SFAC2/8D0       
C.. correction ELA WAS november 1989................................    
C.. this correction reconstructs properly double collinear limit        
C.. and affects below photon-fermion angle  <0.1 amel/ene       
      SFAC1  =  2D0/(PP*A1*B1)      
      SFAC2  =  2D0/(PP*A2*B2)      
      DELT=(AM2/(2D0*PP))**2*(B2**2*A1**2+A2**2*B1**2)* 
     #  ( B1*B2/(A1*A2)/(A1+A2)**2  
     #   +A1*A2/(B1*B2)/(B1+B2)**2  )       
      WMINF=2D0*DELT/(X1+X2)        
      WMM=WWM(A1,B1)*WWM(A2,B2)+WMINF       
      F1 = X1*SFAC1*SFAC2/8D0*WMM   
      F2 = X2*SFAC1*SFAC2/8D0*WMM   
C...end of correction............................................       
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
        
*COMPACT FORM FOR DOUBLE FINAL BREMSTR CROSS SECTION    
      SUBROUTINE DFINAPR(SECTION)   
C     ******************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      PARAMETER(PI=3.1415926535897932D0, ALFINV=137.03604D0)        
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)     
      COMMON / BHPAR2 / CMSENE,AMEL 
      SAVE / MOMS /,/ BHPAR2 /
      DIMENSION XX(4)       
        
      ALFA=1D0/ALFINV       
      S0=CMSENE**2  
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2     
     $   -(P2(3)+Q2(3))**2  
      DO 10 I=1,4   
      XX(I)=P1(I)+Q1(I)     
  10  CONTINUE  
        
      CALL EFDIST2(XX,P1,Q1,P2,Q2,PK1,PK2,DIST2)
! Normalization        
!####      SECTION=(ALFA/4D0/PI**2)**2*ALFA**2/4D0/S0 *DIST2*S0/S1    
      SECTION=  DIST2*S0/S1    
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
c      SFAC1  =  2D0/(PP*AA1*BB1)    
c      SFAC2  =  2D0/(PP*AA2*BB2)    
c      DELT=(AM2/(2D0*PP))**2*(B2**2*A1**2+A2**2*B1**2)* 
c     #  ( B1*B2/(A1*A2)/(A1+A2)**2  
c     #   +A1*A2/(B1*B2)/(B1+B2)**2  )       
c      WMINF=2D0*DELT/(X1+X2)        
c      WMM=WWM(A1,B1)*WWM(A2,B2)+WMINF       
c      F1 = X1*SFAC1*SFAC2/8D0*WMM   
c      F2 = X2*SFAC1*SFAC2/8D0*WMM   
C...end of correction............................................       
      END       

      FUNCTION BORNVV(SVARI,COSTHE) 
C     ***********************************   
C THIS ROUTINE PROVIDES BORN DIFFERENTIAL CROSS SECTION 
C a version without COMPLEX*16      
C     ***********************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON / KEYYFS / KEYZET,KEYBRM,KEYFIX,KEYRED,KEYWGT  
      COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF
      SAVE   / KEYYFS /,/ WEKING /
ccc      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE  
ccc      SAVE   / WEKINP /  
ccc      COMMON / BHPAR1 / CMS,AMFIN   
ccc      COMMON / BHPAR2 / CMSENE,AMEL 
ccc      SAVE   / BHPAR1 /,/ BHPAR2 /
ccc      COMMON / COEFF  / VE,AE  
ccc      SAVE   / COEFF  /


ccc      VF=VE     
ccc      AF=AE     
      QE= -1D0  
      QF= -1D0  
      AA= 4D0*SQRT(SINW2*(1D0-SINW2))       
      VE= (-1D0+4*SINW2)/AA 
      AE= 1D0/AA    
      VF= (-1D0+4*SINW2)/AA 
      AF= 1D0/AA    
      if(keyzet.eq.0) then
       VE=0
       AE=0
       VF=0
       AF=0
      endif
c%%%%%%
c      write(6,*) " SINW2,AMAZ,GAMMZ",   SINW2,AMAZ,GAMMZ
c      write(6,*) " VE,AE,VF,AF",   VE,AE,VF,AF 
c%%%%%%
      S = SVARI     
      CHI2 = S**2/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)      
      RECHI=(S-AMAZ**2)*S/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)      
      XE= VE**2 +AE**2      
      XF= VF**2 +AF**2      
      YE= 2*VE*AE   
      YF= 2*VF*AF   
      FF0= QE**2*QF**2 +2*RECHI*QE*QF*VE*VF +CHI2*XE*XF 
      FF1=     +2*RECHI*QE*QF*AE*AF +CHI2*YE*YF 
!####
       BORN    = (1D0+ COSTHE**2)*FF0 +2D0*COSTHE*FF1 
ccc       BORN    = 1D0+ COSTHE**2  
c      BORN    = (1D0+ COSTHE**2     
c     1 +4D0*(AMEL**2+AMFIN**2)/S*(1D0-COSTHE**2)        
c     1 +16*AMEL**2*AMFIN**2/S**2*COSTHE**2)*FF0 
c     1 +2D0*COSTHE*FF1   
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
        

