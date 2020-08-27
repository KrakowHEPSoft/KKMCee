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
      include 'TAUDCDsize.inc'
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
      include 'TAUDCDsize.inc'
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
      include 'TAUDCDsize.inc'
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
      include 'TAUDCDsize.inc'
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
      include 'TAUDCDsize.inc'
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
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
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
      SUBROUTINE ANGULU(PD1,PD2,Q1,Q2,COSTHE)
      REAL*8 PD1(4),PD2(4),Q1(4),Q2(4),COSTHE,P(4),QQ(4),QT(4)
C take effective beam which is less massive, it should be irrelevant
C but in case HEPEVT is particulary dirty may help.
C this routine calculate reduced system transver and cosine of scattering 
C angle.

      XM1=ABS(PD1(4)**2-PD1(3)**2-PD1(2)**2-PD1(1)**2)
      XM2=ABS(PD2(4)**2-PD2(3)**2-PD2(2)**2-PD2(1)**2)
      IF (XM1.LT.XM2) THEN
        SIGN=1D0
        DO K=1,4
          P(K)=PD1(K)
        ENDDO
      ELSE
        SIGN=-1D0
        DO K=1,4
          P(K)=PD2(K)
        ENDDO
      ENDIF
C calculate space like part of P (in Z restframe)
      DO K=1,4
       QQ(K)=Q1(k)+Q2(K)
       QT(K)=Q1(K)-Q2(K)
      ENDDO

       XMQQ=SQRT(QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2)

       QTXQQ=QT(4)*QQ(4)-QT(3)*QQ(3)-QT(2)*QQ(2)-QT(1)*QQ(1)
      DO K=1,4
       QT(K)=QT(K)-QQ(K)*QTXQQ/XMQQ**2
      ENDDO

       PXQQ=P(4)*QQ(4)-P(3)*QQ(3)-P(2)*QQ(2)-P(1)*QQ(1)
      DO K=1,4
       P(K)=P(K)-QQ(K)*PXQQ/XMQQ**2
      ENDDO
C calculate costhe
       PXP  =SQRT(p(1)**2+p(2)**2+p(3)**2-p(4)**2)
       QTXQT=SQRT(QT(3)**2+QT(2)**2+QT(1)**2-QT(4)**2)
       PXQT =P(3)*QT(3)+P(2)*QT(2)+P(1)*QT(1)-P(4)*QT(4)
       COSTHE=PXQT/PXP/QTXQT
       COSTHE=COSTHE*SIGN
      END

      FUNCTION PLZAP0(IDE,IDF,SVAR,COSTH0)
C this function calculates probability for the helicity +1 +1 configuration
C of taus for given Z/gamma transfer and COSTH0 cosine of scattering angle
      REAL*8 PLZAP0,SVAR,COSTHE,COSTH0,T_BORN

      COSTHE=COSTH0
C >>>>>      IF (IDE*IDF.LT.0) COSTHE=-COSTH0 ! this is probably not needed ID
C >>>>>      of first beam is used by T_GIVIZ0 including sign

      IF (IDF.GT.0) THEN
        CALL INITWK(IDE,IDF,SVAR)
      ELSE
        CALL INITWK(-IDE,-IDF,SVAR)
      ENDIF
      PLZAP0=T_BORN(0,SVAR,COSTHE,1D0,1D0)
     $  /(T_BORN(0,SVAR,COSTHE,1D0,1D0)+T_BORN(0,SVAR,COSTHE,-1D0,-1D0))


C      write(*,*) 'svar=  ',  svar
C      write(*,*) 'COSTHE=',  COSTHE
C      write(*,*) ide,'  ',idf
C      write(*,*) 'PLZAP0=', PLZAP0
C      COSTHE=0.999
C      write(*,*) 'TBORN+=', T_BORN(0,SVAR,COSTHE,1D0,1D0)

C      COSTHE=-0.999
C      write(*,*) 'TBORN-=', T_BORN(0,SVAR,COSTHE,-1D0,-1D0)
C 100  format (A,E8.4)

!      PLZAP0=0.5
      END
      FUNCTION T_BORN(MODE,SVAR,COSTHE,TA,TB)
C ----------------------------------------------------------------------
C THIS ROUTINE PROVIDES BORN CROSS SECTION. IT HAS THE SAME         
C STRUCTURE AS FUNTIS AND FUNTIH, THUS CAN BE USED AS SIMPLER       
C EXAMPLE OF THE METHOD APPLIED THERE                               
C INPUT PARAMETERS ARE: SVAR    -- transfer
C                       COSTHE  -- cosine of angle between tau+ and 1st beam
C                       TA,TB   -- helicity states of tau+ tau-
C
C     called by : BORNY, BORAS, BORNV, WAGA, WEIGHT
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8             SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      REAL*8            SEPS1,SEPS2
C=====================================================================
      COMMON / T_GSWPRM /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
      COMPLEX*16 ABORN(2,2),APHOT(2,2),AZETT(2,2)
      COMPLEX*16 XUPZFP(2),XUPZIP(2)
      COMPLEX*16 ABORNM(2,2),APHOTM(2,2),AZETTM(2,2)
      COMPLEX*16 PROPA,PROPZ
      COMPLEX*16 XR,XI
      COMPLEX*16 XUPF,XUPI
      COMPLEX*16 XTHING
      DATA XI/(0.D0,1.D0)/,XR/(1.D0,0.D0)/
      DATA MODE0 /-5/
      DATA IDE0 /-55/
      DATA SVAR0,COST0 /-5.D0,-6.D0/
      DATA PI /3.141592653589793238462643D0/
      DATA SEPS1,SEPS2 /0D0,0D0/
C
C MEMORIZATION =========================================================
      IF ( MODE.NE.MODE0.OR.SVAR.NE.SVAR0.OR.COSTHE.NE.COST0
     $    .OR.IDE0.NE.IDE)THEN
C
        KEYGSW=1
C ** PROPAGATORS
        IDE0=IDE
        MODE0=MODE
        SVAR0=SVAR
        COST0=COSTHE
        SINTHE=SQRT(1.D0-COSTHE**2)
        BETA=SQRT(MAX(0D0,1D0-4D0*AMFIN**2/SVAR))
C I MULTIPLY AXIAL COUPLING BY BETA FACTOR.
        XUPZFP(1)=0.5D0*(XUPZF(1)+XUPZF(2))+0.5*BETA*(XUPZF(1)-XUPZF(2))
        XUPZFP(2)=0.5D0*(XUPZF(1)+XUPZF(2))-0.5*BETA*(XUPZF(1)-XUPZF(2))
        XUPZIP(1)=0.5D0*(XUPZI(1)+XUPZI(2))+0.5*(XUPZI(1)-XUPZI(2))
        XUPZIP(2)=0.5D0*(XUPZI(1)+XUPZI(2))-0.5*(XUPZI(1)-XUPZI(2))
C FINAL STATE VECTOR COUPLING
        XUPF     =0.5D0*(XUPZF(1)+XUPZF(2))
        XUPI     =0.5D0*(XUPZI(1)+XUPZI(2))
        XTHING   =0D0

        PROPA =1D0/SVAR
        PROPZ =1D0/DCMPLX(SVAR-AMZ**2,SVAR/AMZ*GAMMZ)
        IF (KEYGSW.EQ.0) PROPZ=0.D0
        DO 50 I=1,2
         DO 50 J=1,2
          REGULA= (3-2*I)*(3-2*J) + COSTHE
          REGULM=-(3-2*I)*(3-2*J) * SINTHE *2.D0*AMFIN/SQRT(SVAR)
          APHOT(I,J)=PROPA*(XUPGI(I)*XUPGF(J)*REGULA)
          AZETT(I,J)=PROPZ*(XUPZIP(I)*XUPZFP(J)+XTHING)*REGULA
          ABORN(I,J)=APHOT(I,J)+AZETT(I,J)
          APHOTM(I,J)=PROPA*DCMPLX(0D0,1D0)*XUPGI(I)*XUPGF(J)*REGULM
          AZETTM(I,J)=PROPZ*DCMPLX(0D0,1D0)*(XUPZIP(I)*XUPF+XTHING)*REGULM
          ABORNM(I,J)=APHOTM(I,J)+AZETTM(I,J)
   50   CONTINUE
      ENDIF
C
C******************
C* IN CALCULATING CROSS SECTION ONLY DIAGONAL ELEMENTS
C* OF THE SPIN DENSITY MATRICES ENTER (LONGITUD. POL. ONLY.)
C* HELICITY CONSERVATION EXPLICITLY OBEYED
      POLAR1=  (SEPS1)
      POLAR2= (-SEPS2)
      BORN=0D0
      DO 150 I=1,2
       HELIC= 3-2*I
       DO 150 J=1,2
        HELIT=3-2*J
        FACTOR=KOLOR*(1D0+HELIC*POLAR1)*(1D0-HELIC*POLAR2)/4D0
        FACTOM=FACTOR*(1+HELIT*TA)*(1-HELIT*TB)
        FACTOR=FACTOR*(1+HELIT*TA)*(1+HELIT*TB)

        BORN=BORN+CDABS(ABORN(I,J))**2*FACTOR
C      MASS TERM IN BORN
        IF (MODE.GE.1) THEN
         BORN=BORN+CDABS(ABORNM(I,J))**2*FACTOM
        ENDIF

  150 CONTINUE
C************
      FUNT=BORN
      IF(FUNT.LT.0.D0)  FUNT=BORN

C
      IF (SVAR.GT.4D0*AMFIN**2) THEN
C PHASE SPACE THRESHOLD FACTOR
        THRESH=SQRT(1-4D0*AMFIN**2/SVAR)
        T_BORN= FUNT*SVAR**2*THRESH
      ELSE
        THRESH=0.D0
        T_BORN=0.D0
      ENDIF
C ZW HERE WAS AN ERROR 19. 05. 1989
!      write(*,*) 'KKKK ',PROPA,PROPZ,XUPGI,XUPGF,XUPZI,XUPZF
!      write(*,*) 'KKKK X',svar,costhe,TA,TB,T_BORN
      END

      SUBROUTINE INITWK(IDEX,IDFX,SVAR)
! initialization routine coupling masses etc.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8             SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      COMMON / T_GSWPRM /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
C
      ENE=SQRT(SVAR)/2
      AMIN=0.511D-3
      SWSQ=0.23147
      AMZ=91.1882
      GAMMZ=2.4952
      IF     (IDFX.EQ. 15) then       
        IDF=2  ! denotes tau +2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSEIF (IDFX.EQ.-15) then
        IDF=-2  ! denotes tau -2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSE
        WRITE(*,*) 'INITWK: WRONG IDFX'
        STOP
      ENDIF

      IF     (IDEX.EQ. 11) then      !electron
        IDE= 2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ.-11) then      !positron
        IDE=-2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ. 13) then      !mu+
        IDE= 2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.-13) then      !mu-
        IDE=-2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.  1) then      !d
        IDE= 4
        AMIN=0.05
      ELSEIF (IDEX.EQ.- 1) then      !d~
        IDE=-4
        AMIN=0.05
      ELSEIF (IDEX.EQ.  2) then      !u
        IDE= 3
        AMIN=0.02
      ELSEIF (IDEX.EQ.- 2) then      !u~
        IDE=-3
        AMIN=0.02
      ELSEIF (IDEX.EQ.  3) then      !s
        IDE= 4
        AMIN=0.3
      ELSEIF (IDEX.EQ.- 3) then      !s~
        IDE=-4
        AMIN=0.3
      ELSEIF (IDEX.EQ.  4) then      !c
        IDE= 3
        AMIN=1.3
      ELSEIF (IDEX.EQ.- 4) then      !c~
        IDE=-3
        AMIN=1.3
      ELSEIF (IDEX.EQ.  5) then      !b
        IDE= 4
        AMIN=4.5
      ELSEIF (IDEX.EQ.- 5) then      !b~
        IDE=-4
        AMIN=4.5
      ELSEIF (IDEX.EQ.  12) then     !nu_e
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 12) then     !nu_e~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  14) then     !nu_mu
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 14) then     !nu_mu~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  16) then     !nu_tau
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 16) then     !nu_tau~
        IDE=-1
        AMIN=0.1D-3

      ELSE
        WRITE(*,*) 'INITWK: WRONG IDEX'
        STOP
      ENDIF

C ----------------------------------------------------------------------
C
C     INITIALISATION OF COUPLING CONSTANTS AND FERMION-GAMMA / Z0 VERTEX
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      ITCE=IDE/IABS(IDE)
      JTCE=(1-ITCE)/2
      ITCF=IDF/IABS(IDF)
      JTCF=(1-ITCF)/2
      CALL T_GIVIZO( IDE, 1,AIZOR,QE,KDUMM)
      CALL T_GIVIZO( IDE,-1,AIZOL,QE,KDUMM)
      XUPGI(1)=QE
      XUPGI(2)=QE
      T3E    = AIZOL+AIZOR
      XUPZI(1)=(AIZOR-QE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      XUPZI(2)=(AIZOL-QE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      CALL T_GIVIZO( IDF, 1,AIZOR,QF,KOLOR)
      CALL T_GIVIZO( IDF,-1,AIZOL,QF,KOLOR)
      XUPGF(1)=QF
      XUPGF(2)=QF
      T3F    =  AIZOL+AIZOR
      XUPZF(1)=(AIZOR-QF*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      XUPZF(2)=(AIZOL-QF*SWSQ)/SQRT(SWSQ*(1-SWSQ))
C
      NDIAG0=2
      NDIAGA=11
      KEYA  = 1
      KEYZ  = 1
C
C
      RETURN
      END

      SUBROUTINE T_GIVIZO(IDFERM,IHELIC,SIZO3,CHARGE,KOLOR)
C ----------------------------------------------------------------------
C PROVIDES ELECTRIC CHARGE AND WEAK IZOSPIN OF A FAMILY FERMION
C IDFERM=1,2,3,4 DENOTES NEUTRINO, LEPTON, UP AND DOWN QUARK
C NEGATIVE IDFERM=-1,-2,-3,-4, DENOTES ANTIPARTICLE
C IHELIC=+1,-1 DENOTES RIGHT AND LEFT HANDEDNES ( CHIRALITY)
C SIZO3 IS THIRD PROJECTION OF WEAK IZOSPIN (PLUS MINUS HALF)
C AND CHARGE IS ELECTRIC CHARGE IN UNITS OF ELECTRON CHARGE
C KOLOR IS A QCD COLOUR, 1 FOR LEPTON, 3 FOR QUARKS
C
C     called by : EVENTE, EVENTM, FUNTIH, .....
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
C
      IF(IDFERM.EQ.0.OR.IABS(IDFERM).GT.4) GOTO 901
      IF(IABS(IHELIC).NE.1)                GOTO 901
      IH  =IHELIC
      IDTYPE =IABS(IDFERM)
      IC  =IDFERM/IDTYPE
      LEPQUA=INT(IDTYPE*0.4999999D0)
      IUPDOW=IDTYPE-2*LEPQUA-1
      CHARGE  =(-IUPDOW+2D0/3D0*LEPQUA)*IC
      SIZO3   =0.25D0*(IC-IH)*(1-2*IUPDOW)
      KOLOR=1+2*LEPQUA
C** NOTE THAT CONVENTIONALY Z0 COUPLING IS
C** XOUPZ=(SIZO3-CHARGE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      RETURN
 901  PRINT *,' STOP IN GIVIZO: WRONG PARAMS.'
      STOP
      END
      SUBROUTINE PHYFIX(NSTOP,NSTART)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      SAVE /LUJETS/ 
C NSTOP NSTART : when PHYTIA history ends and event starts.
      NSTOP=0
      NSTART=1
      DO I=1, N
       IF(K(I,1).NE.21) THEN
           NSTOP = I-1
           NSTART= I
           GOTO 500
       ENDIF
      ENDDO
 500  CONTINUE
      END
 

      SUBROUTINE TAUPI0(PI0,K)
C no initialization required. Must be called once after every:
C   1)    CALL DEKAY(1+10,...)
C   2)    CALL DEKAY(2+10,...)
C   3)    CALL DEXAY(1,...)
C   4)    CALL DEXAY(2,...)
C subroutine to decay originating from TAUOLA's taus: 
C 1) etas (with CALL TAUETA(JAK))
C 2) later pi0's from taus.
C 3) extensions to other applications possible. 
C this routine belongs to >tauola universal interface<, but uses 
C routines from >tauola< utilities as well.  25.08.2005      
C this is the hepevt class in old style. No d_h_ class pre-name
 
C position of taus, must be defined by host program:
      COMMON /TAUPOS/ NP1,NP2
c
      REAL  PHOT1(4),PHOT2(4)
      REAL*8  R,X(4),Y(4),PI0(4)
      INTEGER JEZELI(3),K
      DATA JEZELI /0,0,0/
      SAVE JEZELI

! random 3 vector on the sphere, masless
        R=SQRT(PI0(4)**2-PI0(3)**2-PI0(2)**2-PI0(1)**2)/2D0
        CALL SPHERD(R,X)
        X(4)=R
        Y(4)=R
        
        Y(1)=-X(1)
        Y(2)=-X(2)
        Y(3)=-X(3)
! boost to lab and to real*4
        CALL bostdq(-1,PI0,X,X)
        CALL bostdq(-1,PI0,Y,Y)
        DO L=1,4
         PHOT1(L)=X(L)
         PHOT2(L)=Y(L)
        ENDDO
C to hepevt
        CALL FILHEP(0,1,22,K,K,0,0,PHOT1,0.0,.TRUE.)
        CALL FILHEP(0,1,22,K,K,0,0,PHOT2,0.0,.TRUE.)

C
      END
      SUBROUTINE TAUETA(PETA,K)
C subroutine to decay etas's from taus. 
C this routine belongs to tauola universal interface, but uses 
C routines from tauola utilities. Just flat phase space, but 4 channels.
C it is called at the beginning of SUBR. TAUPI0(JAK)
C and as far as hepevt search it is basically the same as TAUPI0.  25.08.2005    
C this is the hepevt class in old style. No d_h_ class pre-name
 
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

C position of taus, must be defined by host program:
c
      REAL  RRR(1),BRSUM(3), RR(2)
      REAL  PHOT1(4),PHOT2(4),PHOT3(4)
      REAL*8    X(4),    Y(4),    Z(4)
      REAL                                YM1,YM2,YM3
      REAL*8  R,RU,PETA(4),XM1,XM2,XM3,XLAM
      REAL*8 a,b,c
      INTEGER K
      XLAM(a,b,c)=SQRT(ABS((a-b-c)**2-4.0*b*c))
C position of decaying particle:
C        DO L=1,4
C          PETA(L)= phep(L,K)  ! eta 4 momentum
C        ENDDO
C       eta cumulated branching ratios:
        BRSUM(1)=0.389  ! gamma gamma
        BRSUM(2)=BRSUM(1)+0.319  ! 3 pi0
        BRSUM(3)=BRSUM(2)+0.237  ! pi+ pi- pi0 rest is thus pi+pi-gamma
        CALL RANMAR(RRR,1) 
        
        IF (RRR(1).LT.BRSUM(1)) THEN ! gamma gamma channel exactly like pi0
! random 3 vector on the sphere, masless   
         R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)/2D0
         CALL SPHERD(R,X) 
         X(4)=R
         Y(4)=R
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
! boost to lab and to real*4
         CALL bostdq(-1,PETA,X,X)  
         CALL bostdq(-1,PETA,Y,Y)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
         ENDDO
C to hepevt
         CALL FILHEP(0,1,22,K,K,0,0,PHOT1,0.0,.TRUE.)
         CALL FILHEP(0,1,22,K,K,0,0,PHOT2,0.0,.TRUE.)
        ELSE ! 3 body channels
         IF(RRR(1).LT.BRSUM(2)) THEN  ! 3 pi0
          ID1= 111
          ID2= 111
          ID3= 111
          XM1=AMPIZ ! masses
          XM2=AMPIZ
          XM3=AMPIZ
         ELSEIF(RRR(1).LT.BRSUM(3)) THEN ! pi+ pi- pi0
          ID1= 211
          ID2=-211
          ID3= 111
          XM1=AMPI ! masses
          XM2=AMPI
          XM3=AMPIZ
         ELSE                            ! pi+ pi- gamma 
          ID1= 211
          ID2=-211
          ID3=  22
          XM1=AMPI ! masses
          XM2=AMPI
          XM3=0.0
         ENDIF
 7       CONTINUE  ! we generate mass of the first pair:
          CALL RANMAR(RR,2)
          R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)
          AMIN=XM1+XM2
          AMAX=R-XM3
          AM2=SQRT(AMIN**2+RR(1)*(AMAX**2-AMIN**2))
C         weight for flat phase space
          WT=XLAM(1D0*R**2,1D0*AM2**2,1D0*XM3**2)
     &      *XLAM(1D0*AM2**2,1D0*XM1**2,1D0*XM2**2)
     &           /R**2                    /AM2**2
         IF (RR(2).GT.WT) GOTO 7

         RU=XLAM(1D0*AM2**2,1D0*XM1**2,1D0*XM2**2)/AM2/2  ! momenta of the
                                              ! first two products
                                              ! in the rest frame of that pair
         CALL SPHERD(RU,X)
         X(4)=SQRT(RU**2+XM1**2)
         Y(4)=SQRT(RU**2+XM2**2)
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
C generate momentum of that pair in rest frame of eta:
         RU=XLAM(1D0*R**2,1D0*AM2**2,1D0*XM3**2)/R/2
         CALL SPHERD(RU,Z)
         Z(4)=SQRT(RU**2+AM2**2)
C and boost first two decay products to rest frame of eta.
         CALL bostdq(-1,Z,X,X)
         CALL bostdq(-1,Z,Y,Y)
C redefine Z(4) to 4-momentum of the last decay product: 
         Z(1)=-Z(1)
         Z(2)=-Z(2)
         Z(3)=-Z(3)
         Z(4)=SQRT(RU**2+XM3**2)
C boost all to lab and move to real*4; also masses
         CALL bostdq(-1,PETA,X,X)
         CALL bostdq(-1,PETA,Y,Y)
         CALL bostdq(-1,PETA,Z,Z)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
          PHOT3(L)=Z(L)
         ENDDO
         YM1=XM1
         YM2=XM2
         YM3=XM3
C to hepevt
         CALL FILHEP(0,1,ID1,K,K,0,0,PHOT1,YM1,.TRUE.)
         CALL FILHEP(0,1,ID2,K,K,0,0,PHOT2,YM2,.TRUE.)
         CALL FILHEP(0,1,ID3,K,K,0,0,PHOT3,YM3,.TRUE.)
        ENDIF


C
      END
      SUBROUTINE TAUK0S(PETA,K)
C subroutine to decay K0S's from taus. 
C this routine belongs to tauola universal interface, but uses 
C routines from tauola utilities. Just flat phase space, but 4 channels.
C it is called at the beginning of SUBR. TAUPI0(JAK)
C and as far as hepevt search it is basically the same as TAUPI0.  25.08.2005   

      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU 
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

C position of taus, must be defined by host program:
      COMMON /TAUPOS/ NP1,NP2
c
      REAL  RRR(1),BRSUM(3)
      REAL  PHOT1(4),PHOT2(4)
      REAL*8    X(4),    Y(4)
      REAL                                YM1,YM2
      REAL*8  R,PETA(4),XM1,XM2,XLAM
      REAL*8 a,b,c
      INTEGER K
      XLAM(a,b,c)=SQRT(ABS((a-b-c)**2-4.0*b*c))
C position of decaying particle:

      
!        DO L=1,4
!          PETA(L)= phep(L,K)  ! K0S 4 momentum  (this is cloned from eta decay)
!        ENDDO
C       K0S cumulated branching ratios:
        BRSUM(1)=0.313  ! 2 PI0
        BRSUM(2)=1.0 ! BRSUM(1)+0.319  ! Pi+ PI-
        BRSUM(3)=BRSUM(2)+0.237  ! pi+ pi- pi0 rest is thus pi+pi-gamma
        CALL RANMAR(RRR,1) 

         IF(RRR(1).LT.BRSUM(1)) THEN  ! 2 pi0
          ID1= 111
          ID2= 111
          XM1=AMPIZ ! masses
          XM2=AMPIZ
         ELSEIF(RRR(1).LT.BRSUM(2)) THEN ! pi+ pi- 
          ID1= 211
          ID2=-211
          XM1=AMPI ! masses
          XM2=AMPI
         ELSE                            ! gamma gamma unused !!!
          ID1= 22
          ID2= 22
          XM1= 0.0 ! masses
          XM2= 0.0
         ENDIF
        
! random 3 vector on the sphere, of equal mass !!  
         R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)/2D0
         R4=R
         R=SQRT(ABS(R**2-XM1**2))
         CALL SPHERD(R,X) 
         X(4)=R4
         Y(4)=R4
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
! boost to lab and to real*4
         CALL bostdq(-1,PETA,X,X)  
         CALL bostdq(-1,PETA,Y,Y)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
         ENDDO

         YM1=XM1
         YM2=XM2
C to hepevt
         CALL FILHEP(0,1,ID1,K,K,0,0,PHOT1,YM1,.TRUE.)
         CALL FILHEP(0,1,ID2,K,K,0,0,PHOT2,YM2,.TRUE.)

C
      END

      subroutine bostdq(idir,vv,pp,q)
*     *******************************
c Boost along arbitrary vector v (see eg. J.D. Jacson, Classical 
c Electrodynamics).
c Four-vector pp is boosted from an actual frame to the rest frame 
c of the four-vector v (for idir=1) or back (for idir=-1). 
c q is a resulting four-vector.
c Note: v must be time-like, pp may be arbitrary.
c
c Written by: Wieslaw Placzek            date: 22.07.1994
c Last update: 3/29/95                     by: M.S.
c 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (nout=6)
      DOUBLE PRECISION v(4),p(4),q(4),pp(4),vv(4)  
      save
!
      do 1 i=1,4
      v(i)=vv(i)
 1    p(i)=pp(i)
      amv=(v(4)**2-v(1)**2-v(2)**2-v(3)**2)
      if (amv.le.0d0) then
        write(6,*) 'bosstv: warning amv**2=',amv
      endif
      amv=sqrt(abs(amv))
      if (idir.eq.-1) then
        q(4)=( p(1)*v(1)+p(2)*v(2)+p(3)*v(3)+p(4)*v(4))/amv
        wsp =(q(4)+p(4))/(v(4)+amv)
      elseif (idir.eq.1) then
        q(4)=(-p(1)*v(1)-p(2)*v(2)-p(3)*v(3)+p(4)*v(4))/amv
        wsp =-(q(4)+p(4))/(v(4)+amv)
      else
        write(nout,*)' >>> boostv: wrong value of idir = ',idir
      endif
      q(1)=p(1)+wsp*v(1)
      q(2)=p(2)+wsp*v(2)
      q(3)=p(3)+wsp*v(3)
      end
 
C wrapper to the random number generator should be included here.
C the actual choice must depend on the user environment.
C this is important for parallel generations to have 
C independent seeds.  In this case it is env. of Tauolapp.
      include '../randg/tauola-random.h'

