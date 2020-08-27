      PROGRAM TAUDEM
C     **************
C NOTE THAT THE ROUTINES ARE NOT LIKE IN CPC DECK THIS IS HISTORICAL !!
C=======================================================================
C====================== DECTES    : TEST OF TAU DECAY LIBRARY===========
C====================== KTORY = 1 : INTERFACE OF KORAL-Z TYPE ==========
C====================== KTORY = 2 : INTERFACE OF KORAL-B TYPE =========
C=======================================================================
C     COMMON  /PAWC/ BLAN(10000)
      COMMON  / / BLAN(10000)
      CHARACTER*7 DNAME
      COMMON / INOUT / INUT,IOUT
      DNAME='KKPI'
!      CALL GLIMIT(20000)
!      CALL GOUTPU(16)
      INUT=5
      IOUT=6
      OPEN(IOUT,FILE="./tauola.output")
       OPEN(INUT,FILE="./dane.dat")
      KTORY=1
      CALL DECTES(KTORY)
      KTORY=2
!      CALL DECTES(KTORY)
C      CALL testresu ! fine tune inputs: masses etc. 
      END
      SUBROUTINE DECTES(KTORY)
C     ************************
      REAL POL(4)
      DOUBLE PRECISION HH(4)
C SWITCHES FOR TAUOLA;
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
C I/O UNITS  NUMBERS
      COMMON / INOUT /  INUT,IOUT
C LUND TYPE IDENTIFIER FOR A1
      COMMON / IDPART / IA1
C /PTAU/ IS USED IN ROUTINE TRALO4
      COMMON /PTAU/ PTAU
      COMMON / TAURAD / XK0DEC,ITDKRC
      REAL*8            XK0DEC
      COMMON /TESTA1/ KEYA1
C special switch for tests of dGamma/dQ**2 in a1 decay
C KEYA1=1 constant width of a1 and rho
C KEYA1=2 free choice of rho propagator (defined in function FPIK)
C         and free choice of a1 mass and width. function g(Q**2)
C         (see formula 3.48 in Comp. Phys. Comm. 64 (1991) 275)
C         hard coded both in Monte Carlo and in testing distribution.
C KEYA1=3 function g(Q**2) hardcoded in the Monte Carlo
C         (it is timy to calculate!), but appropriately adjusted in
C         testing distribution.
C-----------------------------------------------------------------------
C          INITIALIZATION
C-----------------------------------------------------------------------
C======================================
      NINP=INUT
      NOUT=IOUT
 3000 FORMAT(A80)
 3001 FORMAT(8I2)
 3002 FORMAT(I10)
 3003 FORMAT(F10.0)
      IF (KTORY.EQ.1) THEN
      READ( NINP,3000) TESTIT
      WRITE(NOUT,3000) TESTIT
      READ( NINP,3001) KAT1,KAT2,KAT3,KAT4,KAT5,KAT6
      READ( NINP,3002) NEVT,JAK1,JAK2,ITDKRC
      READ( NINP,3003) PTAU,XK0DEC
      ENDIF
C======================================
C control output
      WRITE(NOUT,'(6A6/6I6)')
     $ 'KAT1','KAT2','KAT3','KAT4','KAT5','KAT6',
     $  KAT1 , KAT2 , KAT3 , KAT4 , KAT5 , KAT6
      WRITE(NOUT,'(4A12/4I12)')
     $  'NEVT','JAK1','JAK2','ITDKRC',
     $   NEVT,  JAK1 , JAK2 , ITDKRC
      WRITE(NOUT,'(2A12/2F12.6)')
     $ 'PTAU','XK0DEC',
     $  PTAU , XK0DEC
C======================================
      JAK=0
C      JAK1=5
C      JAK2=5
C LUND IDENTIFIER (FOR TAU+) -15
      IF (KTORY.EQ.1) THEN
        IDFF=-15
      ELSE
        IDFF= 15
      ENDIF
C KTO=1 DENOTES TAU DEFINED BY IDFF (I.E. TAU+)
C KTO=2 DENOTES THE OPPOSITE        (I.E. TAU-)
      KTO=2
      IF (KTO.NE.2) THEN
        PRINT *, 'for the sake of these tests KTO has to be 2'
        PRINT *, 'to change tau- to tau+ change IDFF from -15 to 15'
        STOP
      ENDIF
C TAU POLARIZATION IN ITS RESTFRAME;
      POL(1)=0.
      POL(2)=0.
      POL(3)=.9
C TAU MOMENTUM IN GEV;
C      PTAU=CMSENE/2.D0
C NUMBER OF EVENTS TO BE GENERATED;
      NEVTES=10
      NEVTES=NEVT
      PRINT *, 'NEVTES= ',NEVTES
      WRITE(IOUT,7011) KEYA1
C
      IF (KTORY.EQ.1) THEN
         WRITE(IOUT,7001) JAK,IDFF,POL(3),PTAU
      ELSE
         WRITE(IOUT,7004) JAK,IDFF,POL(3),PTAU
      ENDIF

C INITIALISATION OF TAU DECAY PACKAGE TAUOLA
C ******************************************

        CALL INIMAS
        CALL INITDK
        CALL INIPHY(0.1D0)
        CALL INISAMPL
C re initialization introduced from C wrappers
        CALL TauolaRedef ! register reinitialization function; could be invoked earlier.
        CALL INIofC      ! call reinitialization; has to be called at this place!
C -----------------

      IF (KTORY.EQ.1) THEN
         CALL DEXAY(-1,POL)
      ELSE 
         CALL DEKAY(-1,HH)
      ENDIF
C-----------------------------------------------------------------------
C          GENERATION
C-----------------------------------------------------------------------
      NEV=0
      DO 300 IEV=1,NEVTES
      NEV=NEV+1
C RESLU INITIALISE THE LUND RECORD


      CALL TAUFIL
C DECAY....
      IF (KTORY.EQ.1) THEN
         CALL DEXAY(KTO,POL)
      ELSE
         CALL DEKAY(KTO,HH)
         CALL DEKAY(KTO+10,HH)
      ENDIF
      CALL LUHEPC(2)
      IF(IEV.LE.44) THEN
       WRITE(IOUT,7002) IEV
       IF (KTORY.NE.1) THEN
         WRITE(IOUT,7003) HH
       ENDIF
C      CALL LULIST(11)
      CALL LULIST(2)
      ENDIF
      IPRI=MOD(NEV,1000)

      IF(IPRI.EQ.1) write(*,*) ' event no: ',NEV,' NEVTES: ',NEVTES
  300 CONTINUE
  301 CONTINUE
C-----------------------------------------------------------------------
C                     POSTGENERATION
C-----------------------------------------------------------------------
      IF (KTORY.EQ.1) THEN
         CALL DEXAY(100,POL)
      ELSE
         CALL DEKAY(100,HH)
      ENDIF
      RETURN
 7001 FORMAT(//4(/1X,15(5H=====))
     $ /,' ',     19X,'  NON INITIALIZED BBB-VERSION OF TAUOLA ',9X,1H ,
     $ /,' ',     19X,'    TESTS OF TAU DECAY ROUTINES         ',9X,1H ,
     $ /,' ',     19X,'    INTERFACE OF THE KORAL-Z TYPE       ',9X,1H ,
     $  2(/,1X,15(5H=====)),
     $ /,5X ,'JAK   =',I7  ,'  KEY DEFINING DECAY TYPE         ',9X,1H ,
     $ /,5X ,'IDFF  =',I7  ,'  LUND IDENTIFIER FOR FIRST TAU   ',9X,1H ,
     $ /,5X ,'POL(3)=',F7.2,'  THIRD COMPONENT OF TAU POLARIZ. ',9X,1H ,
     $ /,5X ,'PTAU  =',F7.2,'  THIRD COMPONENT OF TAU MOM. GEV ',9X,1H ,
     $  2(/,1X,15(5H=====))/)
 7002 FORMAT(///1X, '===== EVENT NO.',I4,1X,5H=====)
 7003 FORMAT(5X,'POLARIMETRIC VECTOR: ',
     $       7X,'HH(1)',7X,'HH(2)',7X,'HH(3)',7X,'HH(4)',
     $ /,    5X,'                     ', 4(1X,F11.8)   )
 7004 FORMAT(//4(/1X,15(5H=====))
     $ /,'  ',     19X,' NON INITIALIZED BBB-VERSION OF TAUOLA ',9X,1H ,
     $ /,'  ',     19X,'    TESTS OF TAU DECAY ROUTINES        ',9X,1H ,
     $ /,'  ',     19X,'    INTERFACE OF THE KORAL-B TYPE      ',9X,1H ,
     $  2(/,1X,15(5H=====)),
     $ /,5X ,'JAK   =',I7  ,'  KEY DEFINING DECAY TYPE         ',9X,1H ,
     $ /,5X ,'IDFF  =',I7  ,'  LUND IDENTIFIER FOR FIRST TAU   ',9X,1H ,
     $ /,5X ,'POL(3)=',F7.2,'  THIRD COMPONENT OF TAU POLARIZ. ',9X,1H ,
     $ /,5X ,'PTAU  =',F7.2,'  THIRD COMPONENT OF TAU MOM. GEV ',9X,1H ,
     $  2(/,1X,15(5H=====))/)
 7011 FORMAT(///1X, '===== TYPE OF CURRENT',I4,1X,5H=====)
      END
      
      SUBROUTINE INISAMPL
C Initialization of parameters used in optimalization of phase space generation
C
C only efficiency/speed of generation depend on the actual choice, 
C unless e.g. for a given channel probablilities (used for parallel presamplers) are negative 
C or their sum is bigger than 1.
C `widths and masses '  of presampling resonances which are unphysical and used to parametrize change
C of generation variables, have to be chosen reasonably as well.
      CALL INISAMPL2
      CALL INISAMPL3
      CALL INISAMPL4
      CALL INISAMPL5
      END

      SUBROUTINE INISAMPL2
      include '../TAUDCDsize.inc'
      COMMON /SAMPL2/ PROB1(NM2),PROB2(NM2),AM2(NM2),GAM2(NM2),AM3(NM2),GAM3(NM2)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           
C                                                                       
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           


C initialization for 2 scalars
      DO INUM=1,NM2
        AM2(INUM) =AMKST
        GAM2(INUM)=GAMKST
        AM3(INUM) =AMRO
        GAM3(INUM)=GAMRO

!         PROB(1) flat
!         PROB(2) K*
!         PROB(3) rho
        IF    (INUM.EQ.1) THEN
          PROB1(INUM)=0.0 !0.2
          PROB2(INUM)=0.0
        ELSEIF(INUM.EQ.2) THEN
          PROB1(INUM)=0.0 !0.2
          PROB2(INUM)=1.0
        ELSEIF(INUM.EQ.3) THEN
          PROB1(INUM)=0.0 !0.2
          PROB2(INUM)=1.0
        ELSEIF(INUM.EQ.4) THEN
          PROB1(INUM)=1.0 !0.2
          PROB2(INUM)=0.0
        ELSE
          PROB1(INUM)=1.0 !0.2
          PROB2(INUM)=0.0

        ENDIF
      ENDDO
      END

      SUBROUTINE INISAMPL3
      include '../TAUDCDsize.inc'
      COMMON /SAMPL3/ PROB1(NM3),PROB2(NM3),AMRX(NM3),GAMRX(NM3),AMRA(NM3),GAMRA(NM3),AMRB(NM3),GAMRB(NM3)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           
C                                                                       
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           

C initialization for 3 scalars
      AMROP=1.1
      GAMROP=0.36
      AMOM=.782
      GAMOM=0.0084


      DO MNUM=1,NM3
C     XXXXA CORRESPOND TO S2 CHANNEL !
       IF(MNUM.EQ.10) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =AMA1
        GAMRX(MNUM)=GAMA1
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
        AMRB(MNUM) =AMRO
        GAMRB(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.1) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =1.57
        GAMRX(MNUM)=0.9
        AMRB(MNUM) =AMKST
        GAMRB(MNUM)=GAMKST
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.2) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =1.57
        GAMRX(MNUM)=0.9
        AMRB(MNUM) =AMKST
        GAMRB(MNUM)=GAMKST
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.3) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =1.27
        GAMRX(MNUM)=0.3
        AMRA(MNUM) =AMKST
        GAMRA(MNUM)=GAMKST
        AMRB(MNUM) =AMKST
        GAMRB(MNUM)=GAMKST
       ELSEIF(MNUM.EQ.4) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =1.27
        GAMRX(MNUM)=0.3
        AMRA(MNUM) =AMKST
        GAMRA(MNUM)=GAMKST
        AMRB (MNUM)=AMKST
        GAMRB(MNUM)=GAMKST
       ELSEIF(MNUM.EQ.5) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =1.27
        GAMRX(MNUM)=0.3
        AMRA(MNUM) =AMKST
        GAMRA(MNUM)=GAMKST
        AMRB(MNUM) =AMRO
        GAMRB(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.6) THEN
        PROB1(MNUM)=0.4
        PROB2(MNUM)=0.4
        AMRX(MNUM) =1.27
        GAMRX(MNUM)=0.3
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
        AMRB(MNUM) =AMKST
        GAMRB(MNUM)=GAMKST
       ELSEIF(MNUM.EQ.7) THEN
        PROB1(MNUM)=0.0
        PROB2(MNUM)=1.0
        AMRX (MNUM)=1.27
        GAMRX(MNUM)=0.9
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
        AMRB(MNUM) =AMRO
        GAMRB(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.8) THEN
        PROB1(MNUM)=0.0
        PROB2(MNUM)=1.0
        AMRX(MNUM) =AMROP
        GAMRX(MNUM)=GAMROP
        AMRB(MNUM) =AMOM
        GAMRB(MNUM)=GAMOM
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
       ELSEIF(MNUM.EQ.9) THEN
        PROB1(MNUM)=0.5
        PROB2(MNUM)=0.5
        AMRX(MNUM) =AMA1
        GAMRX(MNUM)=GAMA1
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
        AMRB(MNUM) =AMRO
        GAMRB(MNUM)=GAMRO
       ELSE
        PROB1(MNUM)=0.0
        PROB2(MNUM)=0.0
        AMRX(MNUM) =AMA1
        GAMRX(MNUM)=GAMA1
        AMRA(MNUM) =AMRO
        GAMRA(MNUM)=GAMRO
        AMRB(MNUM) =AMRO
        GAMRB(MNUM)=GAMRO
       ENDIF

      ENDDO
      END

      SUBROUTINE INISAMPL4
      include '../TAUDCDsize.inc'
      COMMON /SAMPL4/ PROB1(NM4),PROB2(NM4),AMRX(NM4),GAMRX(NM4),AMRA(NM4),GAMRA(NM4)

      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           
C                                                                       
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           




C initialization for 4 scalars
      DO INUM=1,NM4
       AMOM=.782
       GAMOM=0.0084
       IF(INUM.EQ.1) THEN
        PROB1(INUM)=.35
        PROB2(INUM)=.35
        AMRX(INUM) =1.2
        GAMRX(INUM)=.46
        AMRA(INUM) =AMOM
        GAMRA(INUM)=GAMOM
       ELSEIF(INUM.EQ.2) THEN
        PROB1(INUM)=0.0
        PROB2(INUM)=0.0
        AMRX(INUM) =1.4
        GAMRX(INUM)=.6
        AMRA(INUM) =AMOM
        GAMRA(INUM)=GAMOM
       ELSEIF(INUM.GE.3.AND.INUM.LE.12) THEN
        PROB1(INUM)=0.0
        PROB2(INUM)=0.0
        AMRX(INUM) =1.4
        GAMRX(INUM)=.6
        AMRA(INUM) =AMOM
        GAMRA(INUM)=GAMOM
       ELSE
        PROB1(INUM)=0.0
        PROB2(INUM)=0.0
        AMRX(INUM) =AMA1
        GAMRX(INUM)=GAMA1
        AMRA(INUM) =AMRO
        GAMRA(INUM)=GAMRO
       ENDIF
      ENDDO
      END

      SUBROUTINE INISAMPL5
      include '../TAUDCDsize.inc'
      REAL*8          AMOM,GAMOM
      COMMON /SAMPL5/ PROBa2(NM5),PROBOM(NM5),ama2(NM5),gama2(NM5),AMOM(NM5),GAMOM(NM5)

      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           
C                                                                       
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU             
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1                
     *                 ,AMK,AMKZ,AMKST,GAMKST                           


C initialization for 5 scalars
      DO INUM=1,NM5
         PROBa2(INUM)=0.7
         PROBOM(INUM)=0.7
         ama2(INUM)=1.260
         gama2(INUM)=0.400
         AMOM(INUM)=.78257
         GAMOM(INUM)=.7
         IF (INUM.EQ.1.OR.INUM.EQ.2) GAMOM(INUM)= 0.00844
      ENDDO
      END


      SUBROUTINE INITDK
* ----------------------------------------------------------------------
*     INITIALISATION OF TAU DECAY PARAMETERS  and routines
*
*     called by : KORALZ
* ----------------------------------------------------------------------

      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBRA / GAMPRT(500),JLIST(500),NCHAN
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      include '../TAUDCDsize.inc'
      INTEGER JMAX
      PARAMETER (JMAX=9)
      COMMON / TAUDCD /IDFFIN(JMAX,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      integer KEY0,KEY1,KEY2,KEY3,KEY4,KEY5,KEY6
      COMMON /METYP/ KEY0(2),KEY1(NM1),KEY2(NM2),KEY3(NM3),
     $               KEY4(NM4),KEY5(NM5),KEY6(NM6)
      INTEGER k
c RChL+new currents switch
      INTEGER IVER
C
      CHARACTER OLDNAMES(NLT)*31
      CHARACTER*80 bxINIT
      PARAMETER (
     $  bxINIT ='(1x,1h*,g17.8,            16x, a31,a4,a4, 1x,1h*)'
     $ )
      REAL*4 PI,POL1(4)

C SWITCHES FOR OPTIONS IN INITIALIZATION:
C =======================================================================
! switch on initialization as in BaBar.
      INTEGER IFBABAR
      COMMON /SETINI/ IFBABAR
      DATA    IFBABAR /1/

C IFBABAR = 0 CLEO initialization
C IFBABAR = 1 BaBar initialization
C IFBABAR = 2 CLEO init. + some currents replaced with listed below
C parametrizations. Note that this may be unphysical 
C form the point of view of overall fit to global quantities.
C Affected channels are then:
C pi- pi0 -> Two options of Belle parametrization (PRD 78, 2008, 072006) 
C            decided by internal flag FF2PIRHO (by default = 2)
C pi- pi- pi+, pi- pi0 pi0 -> RChL initialization (PRD 88, 2013, 093012)
C pi- pi0 pi0 pi0, pi- pi- pi+ pi0 -> Two options steered by internal flag IFKARL,
C             IFKARL = 0 (default) CPC 146 2002, hep=ph/0201149
C             IFKARL = 1  hep-ph/9410260 (probably)

! switches 5pi mdes to description of Acta Phys Polon 2008
      INTEGER IF5PIAPP
      DATA    IF5PIAPP /0/ ! if 1 it overrules in part IFBABAR
C switch to initialize that all modes are represented equally.
      INTEGER IFEQUALBR
      DATA    IFEQUALBR /0/ ! if 1 it overrules both IFBABAR and IF5PIAPP

C==========================================================================
C      PARAMETER (NMODE=196,NM1=40,NM2=71,NM3=19,NM4=32,NM5=21,NM6=13)
*
*
* LIST OF BRANCHING RATIOS (from old times)
CAM normalised to e nu nutau channel
CAM                  enu   munu   pinu  rhonu   A1nu   Knu    K*nu   pi
CAM   DATA JLIST  /    1,     2,     3,     4,     5,     6,     7,

*AM   DATA GAMPRT /1.000,0.9730,0.6054,1.2432,0.8432,0.0432,O.O811,0.616
*AM
*AM  multipion decays
*
*    conventions of particles names
*                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
*                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
*                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
*                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
*                 ET,P-,P0 , P-,P0,GM  , P-,P0,P0
*                  9, 1, 2  , 1, 2, 8  ,  1, 2, 2
*

C
      DIMENSION NPIK(NMODE),NPIK6(NM6)
      DIMENSION NOPIK1(JMAX,NM1),NOPIK2(JMAX,NM2),NOPIK3(JMAX,NM3)
      DIMENSION NOPIK6(JMAX,NM6),NOPIK5(JMAX,NM5),NOPIK4(JMAX,NM4)
      DIMENSION KEYstrt0(2),KEYstrt1(NM1),KEYstrt2(NM2),KEYstrt3(NM3)
      DIMENSION             KEYstrt4(NM4),KEYstrt5(NM5),KEYstrt6(NM6)

*AM   outgoing multiplicity and flavors of multi-pion /multi-K modes    
      DATA   NPIK6  /                
     x                                                    5,    ! old npi starts here
     2                              6,                    6,
     a                              7,                    7,    ! new (may 2004)
     b                              8,                    6,    ! new (may 2004)
     c                              6,                    6,    ! new (may 2004)
     d                              6,                    6,    ! new (may 2004)
     e                              6,                    6    ! new (may 2004)
     o                                                    /    ! new (may 2004)          
C MatrEl key 0- channel not initialized,  1- constant ME flat phase space
C            2- default ME,               3- default ME, but one stable spin>0
C            4- default ME wrapped curr., 5- wrapped ME  
      DATA KEYstrt0 / 2, 2/
      DATA KEYstrt1 / 2, 2, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 0, 0/
      DATA KEYstrt2 / 2, 2, 2, 2, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 0, 0,
     a                0, 0, 0, 0, 0,   0, 0, 0, 0, 0,
     a                0/
      DATA KEYstrt3 / 2, 2, 2, 2, 2,   2, 2, 3, 2, 2,
     a                1, 1, 1, 1, 1,   1, 1, 1, 0   /
      DATA KEYstrt4 / 2, 2, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 1, 1, 1, 1,
     a                1, 1, 1, 1, 1,   1, 0, 0, 0, 0,
     a                0, 0/
      DATA KEYstrt5 / 2, 2, 2, 2, 2,   2, 1, 1, 2, 0,
     a                0, 0, 0, 0, 0,   0, 0, 0, 0, 0,
     a                0/
      DATA KEYstrt6 / 2, 2, 2, 2, 2,   2, 2, 0, 0, 0,
     a                0, 0, 0/

      DATA NOPIK4 / -1,-1, 1, 2, 0, 0,3*0,     2, 2, 2,-1, 0, 0,3*0,  
     a              12,-11,-11,11, 0, 0,3*0,  14,-13,-13,13, 0, 0,3*0,      ! new (may 2004)
     b              12,-11,-13,13, 0, 0,3*0,  14,-13,-11,11, 0, 0,3*0,      ! new (may 2004)
     c              -3, 2, 2, 2, 0, 0,3*0,     2, 2, 9,-3, 0, 0,3*0,     ! new (may 2004)
     d               2, 2, 4,-1, 0, 0,3*0,     2, 4, 9,-1, 0, 0,3*0,     ! new (may 2004)
     e               2,-1, 1,-3, 0, 0,3*0,     4,-1, 1,-1, 0, 0,3*0,     ! new (may 2004)
     a               2, 2, 9,-1, 0, 0,3*0,     4,-4, 9,-1, 0, 0,3*0,     ! new (sep 2004)
     b               4,-4, 2,-1, 0, 0,3*0,     4,-4, 4,-1, 0, 0,3*0,     ! new (sep 2004)
     c               4, 2, 2,-3, 0, 0,3*0,     4,-4, 2,-3, 0, 0,3*0,     ! new (sep 2004)
     d               2, 4, 9,-3, 0, 0,3*0,    -1, 1,-1, 9, 0, 0,3*0,     ! new (sep 2004)
     e              -1, 3,-3, 2, 0, 0,3*0,    -3, 3,-3, 2, 0, 0,3*0,     ! new (sep 2004)
     a              -3, 3,-3, 4, 0, 0,3*0,    -3, 1,-1, 4, 0, 0,3*0,     ! new (sep 2004)
     b              -3, 3,-1, 4, 0, 0,3*0,    -1, 1,-1,-223, 0, 0,3*0,     ! new (sep 2004)
     c               4, 2, 2,-1, 0, 0,3*0,     4, 2, 2,-1, 0, 0,3*0,     ! new (sep 2004)
     d               4, 2, 2,-1, 0, 0,3*0,     4, 2, 2,-1, 0, 0,3*0,     ! new (sep 2004)
     e               4, 2, 2,-1, 0, 0,3*0,     4, 2, 2,-1, 0, 0,3*0/     ! new (sep 2004)

      DATA NOPIK5 /
     1              -1,-1, 1, 2, 2, 0,3*0,  
     a              -1,-1, 1, 2, 2, 0,3*0,     2, 2, 2, 2, 2, 0,3*0,     ! new (may 2004)
     a               1,-1,-1, 2, 2, 0,3*0,    -1, 2, 2, 2, 2, 0,3*0,     ! new (may 2004)
     a              -1, 1, 1,-1,-1, 0,3*0,    -1,-1, 1, 1,-3, 0,3*0,     ! new (may 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1, 2, 2, 2, 2, 0,3*0,     ! new (may 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0,     ! new (may 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0,     ! new (sep 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0,     ! new (sep 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0,     ! new (sep 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0,     ! new (sep 2004)
     a              -1,-1, 1, 2, 4, 0,3*0,    -1,-1, 1, 2, 4, 0,3*0/     ! new (sep 2004)

      DATA NOPIK6 /
     x                                        -1,-1,-1, 1, 1, 0,3*0,     ! old npi starts here
     2              -1,-1,-1, 1, 1, 2,3*0,    -1,-1, 1, 2, 2, 2,3*0, 
     a              -1,-1,-1, 1, 1,2*2,2*0,   -1,-1,-1,-1, 1,2*1,2*0,     ! new (may 2004)  7PI 
     b              -1,-1,-1,-1,3*1, 2, 0,    -1,-1, 1, 1,-3, 2,3*0,     ! new (may 2004)
     c              -1,-1,-1, 1, 1, 1,3*0,    -1,-1, 1, 2, 2, 1,3*0,     ! new (may 2004)
     d              -1,-1,-1, 1, 1, 1,3*0,    -1,-1, 1, 2, 2, 1,3*0,     ! new (may 2004)
     e              -1,-1,-1, 1, 1, 1,3*0,    -1,-1, 1, 2, 2, 1,3*0/     ! new (may 2004)

      DATA NOPIK3 /
     3              -3,-1, 3, 0, 0, 0,3*0,    -4,-1, 4, 0, 0, 0,3*0,  
     4              -3, 2,-4, 0, 0, 0,3*0,     2, 2,-3, 0, 0, 0,3*0,  
     5              -3,-1, 1, 0, 0, 0,3*0,    -1, 4, 2, 0, 0, 0,3*0,  
     6               9,-1, 2, 0, 0, 0,3*0,    -1, 2, 8, 0, 0, 0,3*0,
C AJWMOD fix sign bug, 2/22/99
     7               2, 2,-1, 0, 0, 0,3*0,                           ! new (may 2004) but useful
     7              -1,-1, 1, 0, 0, 0,3*0,    -3,-3, 3, 0, 0, 0,3*0, ! new (may 2004)
     7              -3, 4, 4, 0, 0, 0,3*0,    -3, 9, 2, 0, 0, 0,3*0, ! new (may 2004)
     7               4, 9,-1, 0, 0, 0,3*0,    -3, 4,113, 0, 0, 0,3*0, ! new (may 2004)
     7              -1,333, 2, 0, 0, 0,3*0,    -3,333, 2, 0, 0, 0,3*0, ! new (may 2004)
     7               4, 9,-3, 0, 0, 0,3*0,     2, 2, 2, 0, 0, 0,3*0/ ! new (may 2004)  ostatnie 3

      DATA NOPIK2 /
     8                                        -1, 2, 0, 0, 0, 0,3*0,
     8               -1,-4, 0, 0, 0, 0,3*0,   -3, 2, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-4, 0, 0, 0, 0,3*0,   -13,-13,13, 0, 0, 0,3*0, ! new (may 2004)
     8               -13,-13, 11, 0, 0, 0,3*0, -13,-11, 13, 0, 0, 0,3*0, ! new (may 2004)
     8               -13,-11, 11, 0, 0, 0,3*0,  13,-11, 11, 0, 0, 0,3*0, ! new (may 2004)
     8               -11,-11, 11, 0, 0, 0,3*0, -1, 1, -11, 0, 0, 0,3*0, ! new (may 2004)
     8               -1, 1, -13, 0, 0, 0,3 *0,   1, -3, -11, 0, 0, 0,3*0, ! new (may 2004)
     8                1, -3, -13, 0, 0, 0,3*0,   -1, 3, -11, 0, 0, 0,3*0, ! new (may 2004)
     8               -1, 3, -13, 0, 0, 0,3*0,   -3, 3, -11, 0, 0, 0,3*0, ! new (may 2004)
     8               -3, 3, -13, 0, 0, 0,3*0,   4,4, -11, 0, 0, 0,3*0, ! new (may 2004)
     8                4, 4, -13, 0, 0, 0,3*0,   -1,-1, 11, 0, 0, 0,3*0, ! new (may 2004)
     8               -1,-1, 13, 0, 0, 0,3*0,   -1, -3, 11, 0, 0, 0,3*0, ! new (may 2004)
     8               -1,-3, 13, 0, 0, 0,3*0,   -3,-3, 11, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-3, 13, 0, 0, 0,3*0,   -13,-13,-2212, 0, 0, 0,3*0, ! new (may 2004)
     8               -13, 13, 2212, 0, 0, 0,3*0,   -11,-11, -2212, 0, 0, 0,3*0, ! new (may 2004)
     8               -11, 11, 2212, 0, 0, 0,3*0,   9,-3, 0, 0, 0, 0,3*0, ! new (may 2004)
     8                9, -1, 0, 0, 0, 0,3*0,   -1,333, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,333, 0, 0, 0, 0,3*0,   -1,223, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,223, 0, 0, 0, 0,3*0,   -1,331, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,331, 0, 0, 0, 0,3*0,    -11,13, 2212, 0, 0, 0,3*0, ! new (may 2004)
     8               11,-13, 2212, 0, 0, 0,3*0,   -11,-13, -2212, 0, 0, 0,3*0, ! new (may 2004)
     8               2, 2, -11, 0, 0, 0,3*0,   2,2, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               2, 9, -11, 0, 0, 0,3*0,   2, 9, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               2, 331, -11, 0, 0, 0,3*0,   2, 331, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               9, 9, -11, 0, 0, 0,3*0,   9,9, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               9, 331, -11, 0, 0, 0,3*0,  9, 331, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               2,310, -11, 0, 0, 0,3*0,   2,310, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               9,310, -11, 0, 0, 0,3*0,   9,310, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               331,310, -11, 0, 0, 0,3*0,  331,310, -13, 0, 0, 0,3*0, ! new (may 2004)
     8               1,-3, 2212, 0, 0, 0,3*0,   -1,-3, -2212, 0, 0, 0,3*0, ! new (may 2004)
     8               3,-1, 2212, 0, 0, 0,3*0,   2,2, 2212, 0, 0, 0,3*0, ! new (may 2004)
     8               2, 9, 2212, 0, 0, 0,3*0,   2,310, 2212, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-3, 0, 0, 0, 0,3*0,   -3,-3, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-3, 0, 0, 0, 0,3*0,   -3,-3, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-3, 0, 0, 0, 0,3*0,   -3,-3, 0, 0, 0, 0,3*0, ! new (may 2004)
     8               -3,-3, 0, 0, 0, 0,3*0,   -3,-3, 0, 0, 0, 0,3*0/ ! new (may 2004)  ostattnie 2

      DATA NOPIK1 /
     o               -1, 0, 0, 0, 0, 0,3*0,   -3, 0, 0, 0, 0, 0,3*0, ! new (may 2004)
     o                8, 0,-11,0, 0, 0,3*0,    8, 0,-13,0, 0, 0,3*0, ! new (may 2004)
     o                2, 0,-11,0, 0, 0,3*0,    2, 0,-13,0, 0, 0,3*0, ! new (may 2004)
     o                9, 0,-11, 0, 0, 0,3*0,   9, 0,-13, 0, 0, 0,3*0, ! new (may 2004)
     o                4, 0,-11, 0, 0, 0,3*0,   4, 0,-13, 0, 0, 0,3*0, ! new (may 2004)
     o               -11,0, 0,223, 0, 0,3*0,  -13, 0, 0,223, 0, 0,3*0, ! new (may 2004)
     o               -11,0, 0,333, 0, 0,3*0,  -13, 0, 0,333, 0, 0,3*0, ! new (may 2004)
     o               -11,0, 0,113, 0, 0,3*0,  -13, 0, 0,113, 0, 0,3*0, ! new (may 2004)
     o                10211, 0, 0, 0, 0, 0,3*0,  10213, 0, 0, 0, 0, 4*0, ! new (may 2004)
     o               -11, 0, 0, 311, 0, 0,3*0,   -13, 0, 0, 311, 0, 0,3*0, ! new (may 2004)
     o               2212, 0, 0, 22, 0, 0,3*0,   2212, 0, 0, 111, 0, 0,3*0, ! new (may 2004)
     o               2212, 0, 0, 221, 0, 0,3*0,   2212, 0, 0, 311, 0, 0,3*0, ! new (may 2004)
     o               -11, 0, 0, 331, 0, 0,3*0,   -13, 0, 0, 331, 0, 0,3*0, ! new (may 2004)
     o               -1, 0, 0, 3122, 0, 0,3*0,   -1, 0, 0, -3122, 0, 0,3*0, ! new (may 2004)
     o               -3, 0, 0, 3122, 0, 0,3*0,   -3, 0, 0, -3122, 0, 0,3*0, ! new (may 2004)
     o               -11, 0, 0, 313, 0, 0,3*0,   -11, 0, 0, -313, 0, 0,3*0, ! new (may 2004)
     o               -13, 0, 0, 313, 0, 0,3*0,   -13, 0, 0, -313, 0, 0,3*0, ! new (may 2004)
     o               -11, 0, 0, 10111, 0, 0,3*0,   -13, 0, 0, 10111, 0, 0,3*0, ! new (may 2004)
     o               -11, 0, 0, 10221, 0, 0,3*0,   -13, 0, 0, 10221, 0, 0,3*0, ! new (may 2004)
     o               -3, 0, 0, 0, 0, 0,3*0,   -3, 0, 0, 0, 0, 0,3*0 /! new (may 2004)

      NCHAN = NMODE + NLT

C Initialization of matrix elements, BRANCHING RATIOS, channel names and final state flavours
C it is grouped into segments accordingly to multiplicity for final state scalars
      if (ifbabar.eq.2) then
      call rchl_parameters
      endif

C=================================================
C LEPTONIC CHANNELS plus channels to be removed
C=================================================
      DO I=1,2
       KEY0(I)=KEYstrt0(I)
      ENDDO

      DO  I = 1,NLT
        JLIST(I) = I

        IF(I.EQ. 1) GAMPRT(I) =0.1800 ! E-  as in CLEO default
        IF(I.EQ. 2) GAMPRT(I) =0.1751 ! MU- as in CLEO default
        IF(I.EQ. 1) OLDNAMES(I)='  TAU-  -->   E-               '
        IF(I.EQ. 2) OLDNAMES(I)='  TAU-  -->  MU-               '

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(I.eq.1)  GAMPRT(I)=0.178651 ! E-
         if(I.eq.2)  GAMPRT(I)=0.173551 ! MU-
        endif
      ENDDO


C=================================================
C 4-scalar CHANNELS 
C=================================================

      DO II = 1,NM4
          I = II+NLT  ! position on the whole list of decay channels

        KEY4(II)=KEYstrt4(II)

        JLIST(I) = I
        NPIK(I-NLT)=4
        IF(II.EQ. 1) GAMPRT(I) =0.0450 ! CLEO br
        IF(II.EQ. 2) GAMPRT(I) =0.0100 ! CLEO br
        IF(II.EQ. 3) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ. 4) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ. 5) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ. 6) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ. 7) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ. 8) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ. 9) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.10) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ.11) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.12) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.13) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.14) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ.15) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.16) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.17) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.18) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ.19) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.20) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.21) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.22) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ.23) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.24) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.25) GAMPRT(I) =0.0100 *0 ! *1000000
        IF(II.EQ.26) GAMPRT(I) =0.0100 *0 ! *1000000

        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  --> 2PI-,  PI+,  PI0   '
        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  --> 3PI0,        PI-   '

        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  --> nu_e e- e- e+      '  !  (may 2004)
        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  --> nu_mu mu- mu- mu+  '  !  (may 2004)
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  --> nu_e e- mu- mu+    '  !  (may 2004)
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> nu_mu mu- e- e+    '  !  (may 2004)
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> K- 3PI0            '  !  (may 2004)
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> 2PI0 ETA K-        '  !  (may 2004)
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> 2PI0 K0  PI-       '  !  (may 2004)
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> PI0  K0  ETA PI-   '  !  (may 2004)
        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> PI0  PI- PI+ K-    '  !  (may 2004) 
        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> K0   PI- PI+ PI-   '  !  (may 2004)
        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> 2PI0 ETA PI-       '  !  (sep 2004)
        IF(II.EQ.14) NAMES(I-NLT)='  TAU-  --> K0 K0B ETA PI-     '  !  (sep 2004)
        IF(II.EQ.15) NAMES(I-NLT)='  TAU-  --> K0 K0B PI0 PI-     '  !  (sep 2004)
        IF(II.EQ.16) NAMES(I-NLT)='  TAU-  --> K0 K0B K0  PI-     '  !  (sep 2004)
        IF(II.EQ.17) NAMES(I-NLT)='  TAU-  --> K0 PI0 PI0 K-      '  !  (sep 2004)
        IF(II.EQ.18) NAMES(I-NLT)='  TAU-  --> K0 K0B PI0 K-      '  !  (sep 2004)
        IF(II.EQ.19) NAMES(I-NLT)='  TAU-  --> PI0  K0  ETA K-    '  !  (sep 2004)
        IF(II.EQ.20) NAMES(I-NLT)='  TAU-  --> PI-PI+PI-  ETA     '  !  (sep 2004)

        IF(II.EQ.21) NAMES(I-NLT)='  TAU-  --> PI-K+ K-   PI0     '  !  (sep 2004) 
        IF(II.EQ.22) NAMES(I-NLT)='  TAU-  --> K- K+ K-   PI0     '  !  (sep 2004)

        IF(II.EQ.23) NAMES(I-NLT)='  TAU-  --> K- K+ K-   K0      '  !  (sep 2004)
        IF(II.EQ.24) NAMES(I-NLT)='  TAU-  --> K- PI+PI-  K0      '  !  (sep 2004)
        IF(II.EQ.25) NAMES(I-NLT)='  TAU-  --> K- K+ PI-  K0      '  !  (sep 2004)
        IF(II.EQ.26) NAMES(I-NLT)='  TAU-  --> PI-PI+PI-  OMEGA   '  !  (sep 2004)
        IF(II.EQ.27) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)
        IF(II.EQ.28) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)
        IF(II.EQ.29) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)
        IF(II.EQ.30) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)
        IF(II.EQ.31) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)  
        IF(II.EQ.32) NAMES(I-NLT)='  TAU-  --> xxxxxxx4xxxxxxxx   '  !  (sep 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK4(J,I-NLT)
        ENDDO

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq.1)  GAMPRT(I)=0.043654 ! 2pi- pi+ pi0 
         if(II.eq.2)  GAMPRT(I)=0.012619 ! pi- 3pi0
        endif 
      ENDDO


C=================================================
C 5-scalar CHANNELS 
C=================================================

      DO II = 1,NM5
          I = II+NLT+NM4  ! position on the whole list of decay channels

        KEY5(II)=KEYstrt5(II)

        JLIST(I) = I
        NPIK(I-NLT)=5
        IF(II.EQ. 1) GAMPRT(I) =0.0009 ! CLEO br
        IF(II.EQ. 4) GAMPRT(I) =0.00
        IF(II.EQ. 5) GAMPRT(I) =0.00
        IF(II.EQ. 6) GAMPRT(I) =0.00
        IF(II.EQ. 7) GAMPRT(I) =0.001  *0 ! *1000000 

        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  --> 2PI-, PI+, 2PI0 old'

        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  --> a1 --> rho omega   '  !  (may 2004)
        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  --> benchmark curr     '  !  (may 2004)
        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  --> 2PI- PI+ 2PI0 app08'  !  (may 2004)
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  --> PI- 4PI0  app08    '  !  (may 2004)
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> 3PI- 2PI+ app08    '  !  (may 2004)
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> 2PI- 2PI+  K-      '  !  (may 2004)
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> 2PI- PI+ PI0 K0    '  !  (may 2004)
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> PI- 4PI0           '  !  (may 2004)
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (may 2004)
        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (may 2004)

        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.14) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.15) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.16) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.17) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.18) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.19) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.20) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)
        IF(II.EQ.21) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx5xxxxxx   '  !  (sep 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK5(J,I-NLT-NM4)
        ENDDO

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq. 1) GAMPRT(I)=0.005011 ! 2pi- pi+ 2pi0
        endif
      ENDDO


C=================================================
C multiple-scalar CHANNELS 
C=================================================

      DO II = 1,NM6
          I = II+NLT+NM4+NM5  ! position on the whole list of decay channels

        KEY6(II)=KEYstrt6(II)

        JLIST(I) = I
        NPIK(I-NLT)=NPIK6(I-NLT-NM4-NM5)
        IF(II.EQ. 1) GAMPRT(I) =0.0004 ! CLEO br
        IF(II.EQ. 2) GAMPRT(I) =0.0003 ! CLEO br
        IF(II.EQ. 3) GAMPRT(I) =0.0005 ! CLEO br
        IF(II.EQ. 4) GAMPRT(I) =0.0005 *0 ! *1000000   ! 7pi
        IF(II.EQ. 5) GAMPRT(I) =0.0005 *0 ! *1000000   ! 7pi
        IF(II.EQ. 6) GAMPRT(I) =0.0005 *0 ! *1000000   ! 8pi
        IF(II.EQ. 7) GAMPRT(I) =0.0005 *0 ! *1000000   ! 6 with K 

        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  --> 3PI-, 2PI+,        '
        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  --> 3PI-, 2PI+,  PI0   '
        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  --> 2PI-,  PI+, 3PI0   '
        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  --> 3pi- 2pi+ 2pi0     '  !  (may 2004)
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  --> 4PI- 3PI+          '  !  (may 2004)
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> 4PI- 3PI+  PI0     '  !  (may 2004)
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> 2PI- 2PI+ K- PI0   '  !  (may 2004)
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)
        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)
        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)
        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> xxxxxxxxxnxxxxxx   '  !  (may 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK6(J,I-NLT-NM4-NM5)
        ENDDO

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq. 1) GAMPRT(I)=0.000789 ! 3pi- 2pi+
         if(II.eq. 2) GAMPRT(I)=0.000183 ! 3pi- 2pi+ pi0
         if(II.eq. 3) GAMPRT(I)=0.000251 ! 2pi- pi+ 3pi0
        endif
      ENDDO


C=================================================
C 3-scalar CHANNELS 
C=================================================


      DO II = 1,NM3
          I = II+NLT+NM4+NM5+NM6  ! position on the whole list of decay channels

        KEY3(II)=KEYstrt3(II)

        JLIST(I) = I
        NPIK(I-NLT)=3
        IF(II.EQ. 1) GAMPRT(I) =0.0015 ! CLEO br
        IF(II.EQ. 2) GAMPRT(I) =0.0015 ! CLEO br
        IF(II.EQ. 3) GAMPRT(I) =0.0015 ! CLEO br
        IF(II.EQ. 4) GAMPRT(I) =0.0005 ! CLEO br
        IF(II.EQ. 5) GAMPRT(I) =0.0050 ! CLEO br
        IF(II.EQ. 6) GAMPRT(I) =0.0055 ! CLEO br
        IF(II.EQ. 7) GAMPRT(I) =0.0017 ! CLEO br
        IF(II.EQ. 8) GAMPRT(I) =0.0013 ! CLEO br
        IF(II.EQ. 9) GAMPRT(I) =0.1790 /2 ! CLEO a1 br/2
        IF(II.EQ.10) GAMPRT(I) =0.1790 /2 ! CLEO a1 br/2
        IF(II.EQ.11) GAMPRT(I) =0.0010 *0 ! *100000  ! K- K- K+
        IF(II.EQ.12) GAMPRT(I) =0.0010 *0 ! *100000  ! K- K0 K0

        IF(II.EQ.13) GAMPRT(I) =0.0010 *0 ! *100000  ! K-   ETA  PI0
        IF(II.EQ.14) GAMPRT(I) =0.0010 *0 ! *100000  ! K0   ETA  PI-
        IF(II.EQ.15) GAMPRT(I) =0.0010 *0 ! *100000  ! K-   K0   RHO0
        IF(II.EQ.16) GAMPRT(I) =0.0010 *0 ! *100000  ! PI-  PHI  PI0
        IF(II.EQ.17) GAMPRT(I) =0.0010 *0 ! *100000  ! K-   PHI  PI0
        IF(II.EQ.18) GAMPRT(I) =0.0010 *0 ! *100000  ! K0   ETA  K- 

        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  -->  K-, PI-,  K+      '
        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  -->  K0, PI-, K0B      '

        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  -->  K-,  PI0, K0      '

        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  --> PI0  PI0   K-      '
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  -->  K-  PI-  PI+      '
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> PI-  K0B  PI0      '
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> ETA  PI-  PI0      '
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> PI-  PI0  GAM      '
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> PI0  PI0  PI-      '
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> PI-  PI-  PI+      '  !  (may 2004)
        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> K-    K-   K+      '  !  (may 2004)
        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> K-    K0   K0      '  !  (may 2004)
        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> K-   ETA  PI0      '  !  (may 2004)
        IF(II.EQ.14) NAMES(I-NLT)='  TAU-  --> K0   ETA  PI-      '  !  (may 2004)
        IF(II.EQ.15) NAMES(I-NLT)='  TAU-  --> K-   K0   RHO0     '  !  (may 2004)
        IF(II.EQ.16) NAMES(I-NLT)='  TAU-  --> PI-  PHI  PI0      '  !  (may 2004)
        IF(II.EQ.17) NAMES(I-NLT)='  TAU-  --> K-   PHI  PI0      '  !  (may 2004)
        IF(II.EQ.18) NAMES(I-NLT)='  TAU-  --> K0   ETA  K-       '  !  (may 2004)
        IF(II.EQ.19) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx3xxxxxx   '  !  (may 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK3(J,I-NLT-NM4-NM5-NM6)
        ENDDO
      
        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq. 1) GAMPRT(I)=0.00159  ! K- pi- K+
         if(II.eq. 2) GAMPRT(I)=0.001672 ! pi- K0 K0
         if(II.eq. 3) GAMPRT(I)=0.001536 ! K- K0 pi0
         if(II.eq. 4) GAMPRT(I)=0.00068  ! K- pi0 pi0
         if(II.eq. 5) GAMPRT(I)=0.003009 ! K- pi- pi+
         if(II.eq. 6) GAMPRT(I)=0.003767 ! pi- K0 pi0
         if(II.eq. 7) GAMPRT(I)=0.00183  ! pi- pi0 ETA
         if(II.eq. 8) GAMPRT(I)=0.000802 ! pi- pi0 GAMMA
         if(II.eq. 9) GAMPRT(I)=0.091783 ! pi- 2pi0
         if(II.eq.10) GAMPRT(I)=0.091783 ! 2pi- pi+
        endif
      ENDDO


C=================================================
C 2-scalar or anomalous CHANNELS 
C=================================================
 
      DO II = 1,NM2
          I = II+NLT+NM4+NM5+NM6+NM3  ! position on the whole list of decay channels

        KEY2(II)=KEYstrt2(II)

        JLIST(I) = I
        NPIK(I-NLT)=2
        IF(II.EQ. 1) GAMPRT(I) =0.2515        ! CLEO br   ! pi-p0
        IF(II.EQ. 2) GAMPRT(I) =0.0134*0.6666 ! CLEO br   ! pi- K0
        IF(II.EQ. 3) GAMPRT(I) =0.0134*0.3334 ! CLEO br   ! K- pi0
        IF(II.EQ. 4) GAMPRT(I) =0.0010        ! CLEO br   ! K-K0
        IF(II.EQ. 5) GAMPRT(I) =0.0000 ! *100000  ! mu mu mu
        IF(II.EQ. 6) GAMPRT(I) =0.0000 ! *100000  ! mu mu e SS 
        IF(II.EQ. 7) GAMPRT(I) =0.0000 ! *100000  ! mu mu e OS  
        IF(II.EQ. 8) GAMPRT(I) =0.0000 ! *100000  ! mu e e OS
        IF(II.EQ. 9) GAMPRT(I) =0.0000 ! *100000  ! mu e e SS
        IF(II.EQ.10) GAMPRT(I) =0.0000 ! *100000  ! e e e 
        IF(II.EQ.11) GAMPRT(I) =0.0000 ! *100000  ! e pi pi  
        IF(II.EQ.12) GAMPRT(I) =0.0000 ! *100000  ! mu pi pi
        IF(II.EQ.13) GAMPRT(I) =0.0000 ! *100000  ! e pi+  K-
        IF(II.EQ.14) GAMPRT(I) =0.0000 ! *100000  ! mu pi+  K-
        IF(II.EQ.15) GAMPRT(I) =0.0000 ! *100000  ! e K+  pi-           
        IF(II.EQ.16) GAMPRT(I) =0.0000 ! *100000  ! mu K+  pi-
        IF(II.EQ.17) GAMPRT(I) =0.0000 ! *100000  ! e K+  K-  
        IF(II.EQ.18) GAMPRT(I) =0.0000 ! *100000  ! mu K+  K-   
        IF(II.EQ.19) GAMPRT(I) =0.0000 ! *100000  ! e K0  K0
        IF(II.EQ.20) GAMPRT(I) =0.0000 ! *100000  ! mu K0  K0   
        IF(II.EQ.21) GAMPRT(I) =0.0000 ! *100000  ! e+ pi- pi-                                    
        IF(II.EQ.22) GAMPRT(I) =0.0000 ! *100000  ! mu+ pi- pi-
        IF(II.EQ.23) GAMPRT(I) =0.0000 ! *100000  ! e+ K- pi- 
        IF(II.EQ.24) GAMPRT(I) =0.0000 ! *100000  ! mu+ K- pi-
        IF(II.EQ.25) GAMPRT(I) =0.0000 ! *100000  ! e+ K- K- 
        IF(II.EQ.26) GAMPRT(I) =0.0000 ! *100000  ! mu+ K- K-    
        IF(II.EQ.27) GAMPRT(I) =0.0000 ! *100000  ! mu- mu- p  SS      
        IF(II.EQ.28) GAMPRT(I) =0.0000 ! *100000  ! mu- mu+ p-  OS
        IF(II.EQ.29) GAMPRT(I) =0.0000 ! *100000  ! e- e- p  SS               
        IF(II.EQ.30) GAMPRT(I) =0.0000 ! *100000  ! e- e+ p-  OS    
        IF(II.EQ.31) GAMPRT(I) =0.0000 ! *100000
        IF(II.EQ.32) GAMPRT(I) =0.0000 ! *100000 
        IF(II.EQ.33) GAMPRT(I) =0.0000 ! *100000 
        IF(II.EQ.34) GAMPRT(I) =0.0000 ! *100000
        IF(II.EQ.35) GAMPRT(I) =0.0000 ! *100000  
        IF(II.EQ.36) GAMPRT(I) =0.0000 ! *100000
        IF(II.EQ.37) GAMPRT(I) =0.0000 ! *100000
        IF(II.EQ.38) GAMPRT(I) =0.0000 ! *100000
        IF(II.EQ.39) GAMPRT(I) =0.000  ! *100000  ! e- mu+ p-
        IF(II.EQ.40) GAMPRT(I) =0.000  ! *100000 ! e+ mu- p-
        IF(II.EQ.41) GAMPRT(I) =0.000  ! *100000 ! e- mu- p+   
        IF(II.EQ.42) GAMPRT(I) =0.000  ! *100000 ! e-  pi0 pi0
        IF(II.EQ.43) GAMPRT(I) =0.000  ! *100000 ! mu-  pi0 pi0                              
        IF(II.EQ.44) GAMPRT(I) =0.000  ! *100000 ! tau- -> e-  pi0 eta        
        IF(II.EQ.45) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu-  pi0 eta
        IF(II.EQ.46) GAMPRT(I) =0.000  ! *100000 ! tau- -> e-  pi0 eta' 
        IF(II.EQ.47) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  pi0 eta'
        IF(II.EQ.48) GAMPRT(I) =0.000  ! *100000 ! tau- -> e- eta eta  
        IF(II.EQ.49) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  eta eta
        IF(II.EQ.50) GAMPRT(I) =0.000  ! *100000 ! tau- -> e- eta eta'              
        IF(II.EQ.51) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  eta eta'             
        IF(II.EQ.52) GAMPRT(I) =0.000  ! *100000 ! tau- -> e- PI0 Ks  
        IF(II.EQ.53) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  PI0 Ks     
        IF(II.EQ.54) GAMPRT(I) =0.000  ! *100000 ! tau- -> e- eta  Ks                  
        IF(II.EQ.55) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  eta Ks         
        IF(II.EQ.56) GAMPRT(I) =0.000  ! *100000 ! tau- -> e- eta'  Ks 
        IF(II.EQ.57) GAMPRT(I) =0.000  ! *100000 ! tau- -> mu  eta' Ks 
        IF(II.EQ.58) GAMPRT(I) =0.000  ! *100000 ! tau- -> p-  pi+ K-      
        IF(II.EQ.59) GAMPRT(I) =0.000  ! *100000 ! tau- -> p+  pi- K-
        IF(II.EQ.60) GAMPRT(I) =0.000  ! *100000 ! tau- -> p-  K+ pi-
        IF(II.EQ.61) GAMPRT(I) =0.000  ! *100000 ! tau- -> p-  pi0 pi0
        IF(II.EQ.62) GAMPRT(I) =0.000  ! *100000 ! tau- -> p- pi0 eta
        IF(II.EQ.63) GAMPRT(I) =0.000  ! *100000 ! tau- -> p- pi0 Ks


        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  -->  PI- PI0           '
        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  -->  PI- K0            '  !  (may 2004)
        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  -->  K-  PI0           '  !  (may 2004)
        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  -->  K-  K0            '  !  (may 2004)
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  -->  mu-mu-mu+ !nu_tau '  !  (may 2004)
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> mu- mu- e+ !nu_tau '  !  (may 2004)
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> mu- e- mu+ !nu_tau '  !  (may 2004)
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> mu- e- e+  !nu_tau '  !  (may 2004)
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> mu+ e- e-  !nu_tau '  !  (may 2004)
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> e- e- e+   !nu_tau '  !  (may 2004)
        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> e-pi+pi-  !nu_tau  '  !  (may 2004)
        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> mu-pi+pi-  !nu_tau '  !  (may 2004)
        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> e-pi+K-  !nu_tau   '  !  (may 2004)
        IF(II.EQ.14) NAMES(I-NLT)='  TAU-  --> mu-pi+K-  !nu_tau  '  !  (may 2004)
        IF(II.EQ.15) NAMES(I-NLT)='  TAU-  --> e-pi-K+  !nu_tau   '  !  (may 2004)
        IF(II.EQ.16) NAMES(I-NLT)='  TAU-  --> mu-pi-K+  !nu_tau  '  !  (may 2004)
        IF(II.EQ.17) NAMES(I-NLT)='  TAU-  --> e-K-K+  !nu_tau    '  !  (may 2004)
        IF(II.EQ.18) NAMES(I-NLT)='  TAU-  --> mu-K-K+  !nu_tau   '  !  (may 2004)
        IF(II.EQ.19) NAMES(I-NLT)='  TAU-  --> e-K0K0  !nu_tau    '  !  (may 2004)
        IF(II.EQ.20) NAMES(I-NLT)='  TAU-  --> mu-K0K0  !nu_tau   '  !  (may 2004)
        IF(II.EQ.21) NAMES(I-NLT)='  TAU-  --> e+pi-pi-  !nu_tau  '  !  (may 2004)
        IF(II.EQ.22) NAMES(I-NLT)='  TAU-  --> mu+pi-pi-  !nu_tau '  !  (may 2004)
        IF(II.EQ.23) NAMES(I-NLT)='  TAU-  --> e+pi-K-  !nu_tau   '  !  (may 2004)
        IF(II.EQ.24) NAMES(I-NLT)='  TAU-  --> mu+pi-K-  !nu_tau  '  !  (may 2004)
        IF(II.EQ.25) NAMES(I-NLT)='  TAU-  --> e+K-K-  !nu_tau    '  !  (may 2004)
        IF(II.EQ.26) NAMES(I-NLT)='  TAU-  --> mu+K-K-  !nu_tau   '  !  (may 2004)
        IF(II.EQ.27) NAMES(I-NLT)='  TAU-  --> mu-mu- p+  !nu_tau '  !  (may 2004)
        IF(II.EQ.28) NAMES(I-NLT)='  TAU-  --> mu-mu+ p-  !nu_tau '  !  (may 2004)
        IF(II.EQ.29) NAMES(I-NLT)='  TAU-  -->  e - e- p+ !nu_tau '  !  (may 2004)
        IF(II.EQ.30) NAMES(I-NLT)='  TAU-  -->  e - e+ p- !nu_tau '  !  (may 2004)
        IF(II.EQ.31) NAMES(I-NLT)='  TAU-  --> eta k-             '  !  (may 2004)
        IF(II.EQ.32) NAMES(I-NLT)='  TAU-  --> eta pi-            '  !  (may 2004)
        IF(II.EQ.33) NAMES(I-NLT)='  TAU-  --> PI-  PHI           '  !  (may 2004)
        IF(II.EQ.34) NAMES(I-NLT)='  TAU-  -->  K-  PHI           '  !  (may 2004)
        IF(II.EQ.35) NAMES(I-NLT)='  TAU-  -->  PI- OMEGA         '  !  (may 2004)
        IF(II.EQ.36) NAMES(I-NLT)='  TAU-  -->  K-  OMEGA         '  !  (may 2004)
        IF(II.EQ.37) NAMES(I-NLT)='  TAU-  -->  PI- ETAprm        '  !  (may 2004)
        IF(II.EQ.38) NAMES(I-NLT)='  TAU-  -->  K-  ETAprm        '  !  (may 2004)
        IF(II.EQ.39) NAMES(I-NLT)='  TAU-  -->  e- mu+ p- !nu_tau '  !  (may 2004)                
        IF(II.EQ.40) NAMES(I-NLT)='  TAU-  -->  e+ mu- p- !nu_tau '  !  (may 2004)               
        IF(II.EQ.41) NAMES(I-NLT)='  TAU-  -->  e- mu- p+ !nu_tau '  !  (may 2004)                      
        IF(II.EQ.42) NAMES(I-NLT)='  TAU-  --> e- PI0 PI0  !nu_tau'  !  (may 2004)                    
        IF(II.EQ.43) NAMES(I-NLT)='  TAU-  --> mu- PI0 PI0 !nu_tau'  !  (may 2004)                      
        IF(II.EQ.44) NAMES(I-NLT)='  TAU-  --> e- PI0 eta !nu_tau '  !  (may 2004)     
        IF(II.EQ.45) NAMES(I-NLT)='  TAU-  --> mu- PI0 eta !nu_tau'  !  (may 2004)
        IF(II.EQ.46) NAMES(I-NLT)='  TAU-  -->  e- PI0 eta_p !nu_t'  !  (may 2004)
        IF(II.EQ.47) NAMES(I-NLT)='  TAU-  --> mu- PI0 eta_p !nu_t'  !  (may 2004)
        IF(II.EQ.48) NAMES(I-NLT)='  TAU-  --> e- eta eta  !nu_tau'  !  (may 2004)
        IF(II.EQ.49) NAMES(I-NLT)='  TAU-  --> mu- eta eta !nu_tau'  !  (may 2004)
        IF(II.EQ.50) NAMES(I-NLT)='  TAU-  --> e- eta eta_p !nu_t '  !  (may 2004)
        IF(II.EQ.51) NAMES(I-NLT)='  TAU-  --> mu- eta eta_p !nu_t'  !  (may 2004)
        IF(II.EQ.52) NAMES(I-NLT)='  TAU-  --> e- PI0 Ks  !nu_tau '  !  (may 2004)
        IF(II.EQ.53) NAMES(I-NLT)='  TAU-  --> mu- PI0 Ks !nu_tau '  !  (may 2004)
        IF(II.EQ.54) NAMES(I-NLT)='  TAU-  --> e- eta  Ks !nu_tau '  !  (may 2004)
        IF(II.EQ.55) NAMES(I-NLT)='  TAU-  --> mu- eta Ks !nu_tau '  !  (may 2004)
        IF(II.EQ.56) NAMES(I-NLT)='  TAU-  --> e- eta_p Ks !nu_tau'  !  (may 2004)
        IF(II.EQ.57) NAMES(I-NLT)='  TAU-  --> mu- eta_p Ks !nu_t '  !  (may 2004)
        IF(II.EQ.58) NAMES(I-NLT)='  TAU-  --> p- pi+ K- !nu_tau  '  !  (may 2004)
        IF(II.EQ.59) NAMES(I-NLT)='  TAU-  --> p+ pi- K-  !nu_tau '  !  (may 2004)
        IF(II.EQ.60) NAMES(I-NLT)='  TAU-  --> p- K+ pi- !nu_tau  '  !  (may 2004)
        IF(II.EQ.61) NAMES(I-NLT)='  TAU-  --> p- pi0 pi0 !nu_tau '  !  (may 2004)
        IF(II.EQ.62) NAMES(I-NLT)='  TAU-  --> p- pi0 eta !nu_tau '  !  (may 2004)
        IF(II.EQ.63) NAMES(I-NLT)='  TAU-  --> p- pi0 Ks !nu_tau  '  !  (may 2004)
        IF(II.EQ.64) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.65) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.66) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.67) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.68) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.69) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.70) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)
        IF(II.EQ.71) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx2xxxxxx   '  !  (may 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK2(J,I-NLT-NM4-NM5-NM6-NM3)
        ENDDO

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq. 1) GAMPRT(I)=0.253754 ! pi- pi0
         if(II.eq. 2) GAMPRT(I) =0.013641*0.6666 ! pi- K0
         if(II.eq. 3) GAMPRT(I) =0.013641*0.3334 ! K- pi0
         if(II.eq. 4) GAMPRT(I)=0.001651 ! K- K0
        endif
      ENDDO


C=================================================
C 1-scalar or anomalous CHANNELS  
C=================================================

      DO II = 1,NM1
          I = II+NLT+NM4+NM5+NM6+NM3+NM2  ! position on the whole list of decay channels

        KEY1(II)=KEYstrt1(II)

        JLIST(I) = I
        NPIK(I-NLT)=1
        IF(II.EQ. 1) GAMPRT(I) =0.1110 ! CLEO br ! pi-
        IF(II.EQ. 2) GAMPRT(I) =0.0071 ! CLEO br ! K- 
        IF(II.EQ. 3) GAMPRT(I) =0.0000 ! *100000 ! e gamma  
        IF(II.EQ. 4) GAMPRT(I) =0.0000 ! *100000 ! mu gamma 
        IF(II.EQ. 5) GAMPRT(I) =0.0000 ! *100000 ! Pi0 e 
        IF(II.EQ. 6) GAMPRT(I) =0.0000 ! *100000 ! PI0 mu
        IF(II.EQ. 7) GAMPRT(I) =0.0000 ! *100000 ! eta e   
        IF(II.EQ. 8) GAMPRT(I) =0.0000 ! *100000 ! eta mu        
        IF(II.EQ. 9) GAMPRT(I) =0.0000 ! *100000 ! e-  K0
        IF(II.EQ.10) GAMPRT(I) =0.0000 ! *100000 ! mu- K0
        IF(II.EQ.11) GAMPRT(I) =0.0000 ! *100000 ! e-  omega
        IF(II.EQ.12) GAMPRT(I) =0.0000 ! *100000 ! mu- omega
        IF(II.EQ.13) GAMPRT(I) =0.0000 ! *100000 ! e-  phi
        IF(II.EQ.14) GAMPRT(I) =0.0000 ! *100000 ! mu- phi
        IF(II.EQ.15) GAMPRT(I) =0.0000 ! *100000 ! e-  rho0 
        IF(II.EQ.16) GAMPRT(I) =0.0000 ! *100000 ! mu- rho0 
        IF(II.EQ.17) GAMPRT(I) =0.0000 ! *100000 ! a0-
        IF(II.EQ.18) GAMPRT(I) =0.0000 ! *100000 ! b1-
        IF(II.EQ.19) GAMPRT(I) =0.0000 ! *100000 ! e- K0
        IF(II.EQ.20) GAMPRT(I) =0.0000 ! *100000 ! mu- K0                                
        IF(II.EQ.21) GAMPRT(I) =0.0000 ! *100000 ! p gamma
        IF(II.EQ.22) GAMPRT(I) =0.0000 ! *100000 ! p pi0 
        IF(II.EQ.23) GAMPRT(I) =0.0000 !  *100000 ! p eta
        IF(II.EQ.24) GAMPRT(I) =0.0000 !  *100000 ! p K0
        IF(II.EQ.25) GAMPRT(I) =0.000  !  *100000 ! e- eta'                                       
        IF(II.EQ.26) GAMPRT(I) =0.000  !*100000 ! mu- eta'
        IF(II.EQ.27) GAMPRT(I) =0.000  !*100000 ! pi- lambda
        IF(II.EQ.28) GAMPRT(I) =0.000  !*100000 ! pi- lambda_bar 
        IF(II.EQ.29) GAMPRT(I) =0.000  !*100000 ! K- lambda           
        IF(II.EQ.30) GAMPRT(I) =0.000  !*100000 ! K- lambda_bar           
        IF(II.EQ.31) GAMPRT(I) =0.000  !  *100000 ! e K*
        IF(II.EQ.32) GAMPRT(I) =0.000  !  *100000 ! e K*_bar
        IF(II.EQ.33) GAMPRT(I) =0.000  !  *100000 ! mu K*  
        IF(II.EQ.34) GAMPRT(I) =0.000  !  *100000 ! mu K*_bar  
        IF(II.EQ.35) GAMPRT(I) =0.000  !  *100000 ! e a0(980)      
        IF(II.EQ.36) GAMPRT(I) =0.000  !  *100000 ! mu a0(980)
        IF(II.EQ.37) GAMPRT(I) =0.000  !  *100000 !  e f0(980) 
        IF(II.EQ.38) GAMPRT(I) =0.000  !  *100000 ! mu f0(980) 

        IF(II.EQ. 1) NAMES(I-NLT)='  TAU-  --> PI-                '  !  (may 2004)
        IF(II.EQ. 2) NAMES(I-NLT)='  TAU-  --> K-                 '  !  (may 2004)
        IF(II.EQ. 3) NAMES(I-NLT)='  TAU-  --> gamma e-   !nu_tau '  !  (may 2004)
        IF(II.EQ. 4) NAMES(I-NLT)='  TAU-  --> gamma mu-  !nu_tau '  !  (may 2004)
        IF(II.EQ. 5) NAMES(I-NLT)='  TAU-  --> PI0 e-     !nu_tau '  !  (may 2004)
        IF(II.EQ. 6) NAMES(I-NLT)='  TAU-  --> PI0 mu-    !nu_tau '  !  (may 2004)
        IF(II.EQ. 7) NAMES(I-NLT)='  TAU-  --> eta e-     !nu_tau '  !  (may 2004)
        IF(II.EQ. 8) NAMES(I-NLT)='  TAU-  --> eta mu-    !nu_tau '  !  (may 2004)
        IF(II.EQ. 9) NAMES(I-NLT)='  TAU-  --> e-  K0     !nu_tau '  !  (may 2004)
        IF(II.EQ.10) NAMES(I-NLT)='  TAU-  --> mu- K0     !nu_tau '  !  (may 2004)

        IF(II.EQ.11) NAMES(I-NLT)='  TAU-  --> e-  omega  !nu_tau '  !  (may 2004)
        IF(II.EQ.12) NAMES(I-NLT)='  TAU-  --> mu- omega  !nu_tau '  !  (may 2004)

        IF(II.EQ.13) NAMES(I-NLT)='  TAU-  --> e-  phi    !nu_tau '  !  (may 2004)
        IF(II.EQ.14) NAMES(I-NLT)='  TAU-  --> mu- phi    !nu_tau '  !  (may 2004)
        IF(II.EQ.15) NAMES(I-NLT)='  TAU-  --> e- rho0    !nu_tau '  !  (may 2004)
        IF(II.EQ.16) NAMES(I-NLT)='  TAU-  --> mu- rho0   !nu_tau '  !  (may 2004)
        IF(II.EQ.17) NAMES(I-NLT)='  TAU-  --> A0-                '  !  (may 2004)
        IF(II.EQ.18) NAMES(I-NLT)='  TAU-  --> B1-                '  !  (may 2004)
        IF(II.EQ.19) NAMES(I-NLT)='  TAU-  --> e- K0    !nu_tau   '  !  (may 2004)
        IF(II.EQ.20) NAMES(I-NLT)='  TAU-  --> mu- K0    !nu_tau  '  !  (may 2004)
        IF(II.EQ.21) NAMES(I-NLT)='  TAU-  -->  p gamma  !nu_tau  '  !  (may 2004)
        IF(II.EQ.22) NAMES(I-NLT)='  TAU-  --> p pi0     !nu_tau  '  !  (may 2004)
        IF(II.EQ.23) NAMES(I-NLT)='  TAU-  --> p eta    !nu_tau   '  !  (may 2004)
        IF(II.EQ.24) NAMES(I-NLT)='  TAU-  -->  p K0   !nu_tau    '  !  (may 2004)
        IF(II.EQ.25) NAMES(I-NLT)='  TAU-  --> e- eta_p  !nu_tau  '  !  (may 2004)
        IF(II.EQ.26) NAMES(I-NLT)='  TAU-  --> mu- eta_p !nu_tau  '  !  (may 2004)
        IF(II.EQ.27) NAMES(I-NLT)='  TAU-  --> pi- lambda !nu_tau '  !  (may 2004)
        IF(II.EQ.28) NAMES(I-NLT)='  TAU-  --> pi- lmb_br !nu_tau '  !  (may 2004)
        IF(II.EQ.29) NAMES(I-NLT)='  TAU-  --> K- lambda  !nu_tau '  !  (may 2004)
        IF(II.EQ.30) NAMES(I-NLT)='  TAU-  --> K- lmb_bar !nu_tau '  !  (may 2004)
        IF(II.EQ.31) NAMES(I-NLT)='  TAU-  --> e-  K*  !nu_tau    '  !  (may 2004)
        IF(II.EQ.32) NAMES(I-NLT)='  TAU-  --> e-  K*_bar !nu_tau '  !  (may 2004)
        IF(II.EQ.33) NAMES(I-NLT)='  TAU-  --> mu- K*_bar !nu_tau '  !  (may 2004)
        IF(II.EQ.34) NAMES(I-NLT)='  TAU-  --> mu-  K*  !nu_tau   '  !  (may 2004)
        IF(II.EQ.35) NAMES(I-NLT)='  TAU-  --> e- a0(980) !nu_tau '  !  (may 2004)
        IF(II.EQ.36) NAMES(I-NLT)='  TAU-  --> mu- a0(980) !nu_tau'  !  (may 2004)
        IF(II.EQ.37) NAMES(I-NLT)='  TAU-  --> e-  f0(980) !nu_tau'  !  (may 2004)
        IF(II.EQ.38) NAMES(I-NLT)='  TAU-  --> mu- f0(980) !nu_tau'  !  (may 2004)
        IF(II.EQ.39) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx1xxxxxx   '  !  (may 2004)
        IF(II.EQ.40) NAMES(I-NLT)='  TAU-  --> xxxxxxxxx1xxxxxx   '  !  (may 2004)

        MULPIK(I-NLT)=NPIK(I-NLT)
        DO J=1,JMAX
         IDFFIN(J,I-NLT)=NOPIK1(J,I-NLT-NM4-NM5-NM6-NM3-NM2)
        ENDDO

        if(IFBABAR.eq.1) then
C JKBZR 30.05.2014  BaBar gamprt for default channels
         gamprt(i)=0.0
         if(II.eq.1) GAMPRT(I)=0.110841 ! pi-
         if(II.eq.2) GAMPRT(I)=0.006946 ! K-
        endif
      ENDDO


c If IF5PIAPP=1 then 5pi channels will be generated with APP methods.
C appropriate gamprt-s are shifted.
      IF(IF5PIAPP.EQ.1) THEN
C       activated channels
        gamprt(NLT+NM4+4)=gamprt(NLT+NM4+1)    ! now 2PI- PI+ 2PI0 app08
        gamprt(NLT+NM4+5)=gamprt(NLT+NM4+9)    ! now PI- 4PI0  app08 
C       NOTE: in BaBar initialization BR for PI- 4PI0 channel is set to 0.0 therefore it won't be generated
        gamprt(NLT+NM4+6)=gamprt(NLT+NM4+NM5+1)! now 3PI- 2PI+ app08
C       de-activated channels
        gamprt(NLT+NM4+NM5+1)=0.0 ! 3PI-, 2PI+, (of multi-pions channel group)
        gamprt(NLT+NM4+9)=0.0   ! PI- 4PI0
        gamprt(NLT+NM4+1)=0.0   ! 2PI-, PI+, 2PI0 old
      ENDIF

      DO I = NLT+NM4+NM5+NM6+NM3+NM2+NM1+1,500  ! initialization of memory not allocated even to multiplicity.
        JLIST(I) = 0
        GAMPRT(I) = 0.
      ENDDO

*
*
* --- COEFFICIENTS TO FIX RATIO OF:
* --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
* --- PROBABILITY OF K0 TO BE KS
* --- PROBABILITY OF K0B TO BE KS
* --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
* --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
* --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
* --- NEGLECTS MASS-PHASE SPACE EFFECTS
      BRA1=1D0 ! 0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667
*

      GFERMI = 1.16637E-5
      CCABIB = 0.975
      GV     = 1.0
      GA     =-1.0

      

* ZW 13.04.89 HERE WAS AN ERROR
      SCABIB = SQRT(1.-CCABIB**2)
      PI =4.*ATAN(1.)
      GAMEL  = GFERMI**2*AMTAU**5/(192*PI**3)
*
*      CALL DEXAY(-1,pol1)
*

      IF(IFEQUALBR.EQ.1) THEN  ! uniform  initialization useful for some tests only
        FAC=0.001
        DO I=1,NCHAN
         GAMPRT(I) = 1D0/NCHAN
        ENDDO
      ENDIF

      RETURN
      END
      FUNCTION DCDMAS(IDENT)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      IF      (IDENT.EQ. 1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ.-1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ. 2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ.-2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ. 3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ.-3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ. 4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ.-4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ. 8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ.-8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ. 9) THEN
        APKMAS=0.5488
      ELSEIF  (IDENT.EQ.-9) THEN
        APKMAS=0.5488
      ELSEIF  (ABS(IDENT).EQ.13) THEN
        APKMAS=AMMU
      ELSEIF  (ABS(IDENT).EQ.11) THEN
        APKMAS=AMEL
      ELSEIF  (ABS(IDENT).EQ.12) THEN
        APKMAS=0.0
      ELSEIF  (ABS(IDENT).EQ.14) THEN
        APKMAS=0.0
      ELSEIF  (ABS(IDENT).EQ.22) THEN
        APKMAS=0.0
      ELSEIF  (ABS(IDENT).EQ.223) THEN
        APKMAS=0.7826
      ELSEIF  (ABS(IDENT).EQ.333) THEN
        APKMAS=1.019
      ELSEIF  (ABS(IDENT).EQ.113) THEN
        APKMAS=0.7755
      ELSEIF  (ABS(IDENT).EQ.331) THEN
        APKMAS=0.958
      ELSEIF  (ABS(IDENT).EQ.2212) THEN
         APKMAS=0.938272046
      ELSEIF  (ABS(IDENT).EQ.10211) THEN
         APKMAS=1.45
      ELSEIF  (ABS(IDENT).EQ.10213) THEN
         APKMAS=1.235
      ELSEIF  (ABS(IDENT).EQ.311) THEN
         APKMAS=AMPIZ
      ELSEIF  (ABS(IDENT).EQ.111) THEN
         APKMAS=AMPIZ
      ELSEIF  (ABS(IDENT).EQ. 221) THEN
         APKMAS=0.5488
      ELSEIF  (ABS(IDENT).EQ. 310) THEN
         APKMAS=0.497614
      ELSEIF  (ABS(IDENT).EQ. 313) THEN
         APKMAS=0.89166
      ELSEIF  (ABS(IDENT).EQ. 3122) THEN
         APKMAS=1.15683
      ELSEIF  (ABS(IDENT).EQ. 10111) THEN
         APKMAS=0.98  
      ELSEIF  (ABS(IDENT).EQ. 10221) THEN
         APKMAS=0.98
   

         
      ELSE
        PRINT *, 'STOP IN APKMAS, WRONG IDENT=',IDENT
        STOP
      ENDIF
      DCDMAS=APKMAS
      END
      FUNCTION LUNPIK(ID,ISGN)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      REAL*4 XIO(1)
      IDENT=ID*ISGN

      IF      (IDENT.EQ. 1) THEN
        IPKDEF=-211
      ELSEIF  (IDENT.EQ.-1) THEN
        IPKDEF= 211
      ELSEIF  (IDENT.EQ. 2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ.-2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ. 3) THEN
        IPKDEF=-321
      ELSEIF  (IDENT.EQ.-3) THEN
        IPKDEF= 321
        
      ELSEIF  (IDENT.EQ. 4) THEN
*
* K0 --> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ.-4) THEN
*
* K0B--> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0B) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ. 8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ.-8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ. 9) THEN
        IPKDEF= 221
      ELSEIF  (IDENT.EQ.-9) THEN
        IPKDEF= 221
      ELSEIF  (ABS(IDENT).EQ.13) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.11) THEN
        IPKDEF=IDENT 
      ELSEIF  (ABS(IDENT).EQ.12) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.14) THEN
        IPKDEF=IDENT 
      ELSEIF  (ABS(IDENT).EQ.223) THEN
        IPKDEF=ABS(IDENT)
      ELSEIF  (ABS(IDENT).EQ.333) THEN
        IPKDEF=ABS(IDENT)
      ELSEIF  (ABS(IDENT).EQ.113) THEN
        IPKDEF=ABS(IDENT)
       ELSEIF  (ABS(IDENT).EQ.331) THEN
        IPKDEF=ABS(IDENT)
      ELSEIF  ((IDENT).EQ.22) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.2212) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.10211) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.10213) THEN
        IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.311) THEN
         IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.311) THEN
         IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.221) THEN
         IPKDEF=IDENT
      ELSEIF  (ABS(IDENT).EQ.310) THEN
         IPKDEF=ABS(IDENT)
      ELSEIF  (ABS(IDENT).EQ.313) THEN
         IPKDEF=IDENT   
      ELSEIF  (ABS(IDENT).EQ.3122) THEN
         IPKDEF=IDENT   
      ELSEIF  (ABS(IDENT).EQ.10111 ) THEN
         IPKDEF=ABS(IDENT)   
      ELSEIF  (ABS(IDENT).EQ.10221 ) THEN
         IPKDEF=ABS(IDENT)


         
      ELSE
        PRINT *, 'STOP IN IPKDEF, WRONG IDENT=',IDENT
        STOP
      ENDIF
      LUNPIK=IPKDEF
      END




      SUBROUTINE TAURDF(KTO)
C THIS ROUTINE CAN BE CALLED BEFORE ANY TAU+ OR TAU- EVENT IS GENERATED
C IT CAN BE USED TO GENERATE TAU+ AND TAU- SAMPLES OF DIFFERENT
C CONTENTS
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / TAUBRA / GAMPRT(500),JLIST(500),NCHAN
C this is a template of routine for reinitialization of parameters
C to set them differently for tau+ and tau-
C it is prepared for cleo version of currents but is anyway inactive 
      IF (KTO.LT.1000) RETURN

      IF (KTO.EQ.1) THEN
C     ==================
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,1)
      BRKS = PKORB(4,3)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ELSE
C     ====
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,2)
      BRKS = PKORB(4,4)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ENDIF
C     =====

      END


      SUBROUTINE INIPHY(XK00)
* ----------------------------------------------------------------------
*     INITIALISATION OF PARAMETERS
*     USED IN QED and/or GSW ROUTINES
* ----------------------------------------------------------------------
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8 PI8,XK00
*
      PI8    = 4.D0*DATAN(1.D0)
      ALFINV = 137.03604D0
      ALFPI  = 1D0/(ALFINV*PI8)
      XK0=XK00
      END

      SUBROUTINE INIMAS
C ----------------------------------------------------------------------
C     INITIALISATION OF MASSES
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      include '../new-currents/RChL-currents/parameter.inc'
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C initialize rchl_parameters first
      CALL rchl_parameters(1)
C
C IN-COMING / OUT-GOING  FERMION MASSES
      AMTAU  = 1.7842
C --- let us update tau mass ...
      AMTAU  = 1.777
      AMNUTA = 0.0 ! 0.010
      AMEL   = 0.0005111
      AMNUE  = 0.0
      AMMU   = 0.105659 
      AMNUMU = 0.0
*
* MASSES USED IN TAU DECAYS

      AMPIZ  = 0.134976
      AMPI   = 0.139570
      AMRO   = 0.77590
      GAMRO  = 0.14790
*C    GAMRO  = 0.666
      AMA1   = 1.251
      GAMA1  = 0.599
      AMK    = 0.493677
      AMKZ   = 0.497672
      AMKST  = 0.89166
      GAMKST = 0.0508
C
C
C IN-COMING / OUT-GOING  FERMION MASSES
!!      AMNUTA = PKORB(1,2)
!!      AMNUE  = PKORB(1,4)
!!      AMNUMU = PKORB(1,6)
C
C MASSES USED IN TAU DECAYS  Cleo settings
!!      AMPIZ  = PKORB(1,7)
!!      AMPI   = PKORB(1,8)
!!      AMRO   = PKORB(1,9)
!!      GAMRO  = PKORB(2,9)
c      AMA1   = 1.275   !! PKORB(1,10)
c      GAMA1  = 0.615   !! PKORB(2,10)
!!      AMK    = PKORB(1,11)
!!      AMKZ  = PKORB(1,12)
!!      AMKST  = PKORB(1,13)
!!      GAMKST = PKORB(2,13)
C


      RETURN
      END
      SUBROUTINE TAUFIL
C     *****************
C SUBSITUTE OF tau PRODUCTION GENERATOR
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / IDFC  / IDFF
C positions of taus in the LUND common block
C it will be used by TAUOLA output routines.
      COMMON /TAUPOS / NPA,NPB
      DIMENSION XPB1(4),XPB2(4),AQF1(4),AQF2(4)
C
C --- DEFINING DUMMY EVENTS MOMENTA
      DO 4 K=1,3
        XPB1(K)=0.0
        XPB2(K)=0.0
        AQF1(K)=0.0
        AQF2(K)=0.0
  4   CONTINUE
        AQF1(4)=AMTAU
        AQF2(4)=AMTAU
C --- TAU MOMENTA
      CALL TRALO4(1,AQF1,AQF1,AM)
      CALL TRALO4(2,AQF2,AQF2,AM)
C --- BEAMS MOMENTA AND IDENTIFIERS
        KFB1= 11*IDFF/IABS(IDFF)
        KFB2=-11*IDFF/IABS(IDFF)
        XPB1(4)= AQF1(4)
        XPB1(3)= AQF1(4)
        IF(AQF1(3).NE.0.0)
     $  XPB1(3)= AQF1(4)*AQF1(3)/ABS(AQF1(3))
        XPB2(4)= AQF2(4)
        XPB2(3)=-AQF2(4)
        IF(AQF2(3).NE.0.0)
     $  XPB2(3)= AQF2(4)*AQF2(3)/ABS(AQF2(3))
C --- Position of first and second tau in LUND common
      NPA=3
      NPB=4
C --- FILL TO LUND COMMON
      CALL FILHEP(  1,3, KFB1,0,0,0,0,XPB1, AMEL,.TRUE.)
      CALL FILHEP(  2,3, KFB2,0,0,0,0,XPB2, AMEL,.TRUE.)
      CALL FILHEP(NPA,1, IDFF,1,2,0,0,AQF1,AMTAU,.TRUE.)
      CALL FILHEP(NPB,1,-IDFF,1,2,0,0,AQF2,AMTAU,.TRUE.)
      END
      SUBROUTINE TRALO4(KTO,P,Q,AM)
C     **************************
C SUBSITUTE OF TRALO4
      REAL  P(4),Q(4)
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /PTAU/ PTAU
      AM=AMAS4(P)
      ETAU=SQRT(PTAU**2+AMTAU**2)
      EXE=(ETAU+PTAU)/AMTAU
      IF(KTO.EQ.2) EXE=(ETAU-PTAU)/AMTAU
      CALL BOSTR3(EXE,P,Q)
C ======================================================================
C         END OF THE TEST JOB
C ======================================================================
      END
      SUBROUTINE FILHEP(N,IST,ID,JMO1,JMO2,JDA1,JDA2,P4,PINV,PHFLAG)
C ----------------------------------------------------------------------
C this subroutine fills one entry into the HEPEVT common
C and updates the information for affected mother entries
C
C written by Martin W. Gruenewald (91/01/28)
C
C     called by : ZTOHEP,BTOHEP,DWLUxy
C ----------------------------------------------------------------------
C
C this is the hepevt class in old style. No d_h_ class pre-name
      INTEGER NMXHEP
      PARAMETER (NMXHEP=4000)
      REAL*8  phep,  vhep ! to be real*4/ *8  depending on host
      INTEGER nevhep,nhep,isthep,idhep,jmohep,
     $        jdahep
      COMMON /hepevt/
     $      nevhep,               ! serial number
     $      nhep,                 ! number of particles
     $      isthep(nmxhep),   ! status code
     $      idhep(nmxhep),    ! particle ident KF
     $      jmohep(2,nmxhep), ! parent particles
     $      jdahep(2,nmxhep), ! childreen particles
     $      phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $      vhep(4,nmxhep)    ! vertex [mm]
* ----------------------------------------------------------------------
      LOGICAL qedrad
      COMMON /phoqed/ 
     $     qedrad(nmxhep)    ! Photos flag
* ----------------------------------------------------------------------
      SAVE hepevt,phoqed


C      PARAMETER (NMXHEP=2000)
C      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
C     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
C      SAVE  /HEPEVT/
C      COMMON/PHOQED/QEDRAD(NMXHEP)
C      LOGICAL QEDRAD
C      SAVE /PHOQED/
      LOGICAL PHFLAG
C
      REAL*4  P4(4)
C
C check address mode
      IF (N.EQ.0) THEN
C
C append mode
        IHEP=NHEP+1
      ELSE IF (N.GT.0) THEN
C
C absolute position
        IHEP=N
      ELSE
C
C relative position
        IHEP=NHEP+N
      END IF
C
C check on IHEP
      IF ((IHEP.LE.0).OR.(IHEP.GT.NMXHEP)) RETURN
C
C add entry
      NHEP=IHEP
      ISTHEP(IHEP)=IST
      IDHEP(IHEP)=ID
      JMOHEP(1,IHEP)=JMO1
      IF(JMO1.LT.0)JMOHEP(1,IHEP)=JMOHEP(1,IHEP)+IHEP
      JMOHEP(2,IHEP)=JMO2
      IF(JMO2.LT.0)JMOHEP(2,IHEP)=JMOHEP(2,IHEP)+IHEP
      JDAHEP(1,IHEP)=JDA1
      JDAHEP(2,IHEP)=JDA2
C
      DO I=1,4
        PHEP(I,IHEP)=P4(I)
C
C KORAL-B and KORAL-Z do not provide vertex and/or lifetime informations
        VHEP(I,IHEP)=0.0
      END DO
      PHEP(5,IHEP)=PINV
C FLAG FOR PHOTOS...
      QEDRAD(IHEP)=PHFLAG
C
C update process:
      DO IP=JMOHEP(1,IHEP),JMOHEP(2,IHEP)
        IF(IP.GT.0)THEN
C
C if there is a daughter at IHEP, mother entry at IP has decayed
          IF(ISTHEP(IP).EQ.1)ISTHEP(IP)=2
C
C and daughter pointers of mother entry must be updated
          IF(JDAHEP(1,IP).EQ.0)THEN
            JDAHEP(1,IP)=IHEP
            JDAHEP(2,IP)=IHEP
          ELSE
            JDAHEP(2,IP)=MAX(IHEP,JDAHEP(2,IP))
          END IF
        END IF
      END DO
C
      RETURN
      END
      include '../../randg/tauola-random.h'
