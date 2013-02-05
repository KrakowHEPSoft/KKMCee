C---------------------------------------------------------------------------C
C --------    SUBROUTINES FOR KEYOPT=4   ------------------------------------C
C---------------------------------------------------------------------------C
      SUBROUTINE Z0DEC2(MODE,XPAR,NPAR)
C     ***********************************
C main generating routine of TOPIK
C process   q qbar ->l lbar photon photon
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      PARAMETER(GNANOB = 389385D0)
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /PARFUN/ XMIN,Z1,Z2
      COMMON /BREMSTR/ KEYBRE,KEYPRO
      COMMON /GENER2/  X1,X2,WTMOD
      COMMON /WEIGHTS/ WTINF, WTYFS, WTEXA
      COMMON /MOMLAB2/ P1L(4),Q1L(4),P2L(4),Q2L(4),PHOT1L(4),PHOT2L(4)
      COMMON /MOMCMS2/ P1(4) ,Q1(4) ,P2(4) ,Q2(4) ,PHOT1(4) ,PHOT2(4)
      dimension phot3(4),test(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      COMMON /FLAVOUR/ IFLEV
      DIMENSION NPAR(99),XPAR(99)
      REAL *4 XPPR(-6:6),ULALPS 
      DIMENSION RN(1),APHOT(4)
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      EXTERNAL FUNSKI

      IF(MODE.EQ.-1) THEN
C     =======================
C ...BX-formats for nice and flexible outputs
      BXOPE =  '(//1X,15(5H*****)    )'
      BXTXT =  '(1X,1H*,                  A48,25X,    1H*)'
      BXL1I =  '(1X,1H*,I17,                 16X, A20,A12,A7, 1X,1H*)'
      BXL1F =  '(1X,1H*,F17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2F =  '(1X,1H*,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXL1G =  '(1X,1H*,G17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2G =  '(1X,1H*,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXCLO =  '(1X,15(5H*****)/   )'

      SVAR  = XPAR(1)**2
      AMINI = XPAR(2)
      AMFIN = XPAR(3)
      XMIN  = XPAR(4)
      PHCUT = XPAR(5)
      PHMAX = XPAR(6)
      PHMIN = XPAR(7)
      PTL   = XPAR(8)
      XMIMAS= XPAR(9)
      KEYBRE= NPAR(1)
      KEYPRO=NPAR(2)
c------GSW parameters
      AMAZ    = XPAR(50)  
      GAMMZ   = XPAR(51)
      SINW2   = XPAR(52) 

C
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '*     *****************         *'
      WRITE(NOUT,BXTXT) '*     ***   Z0DEC2   ***        *'
      WRITE(NOUT,BXTXT) '*     *****************         *'
      WRITE(NOUT,BXTXT) '*    September      1993         *'
      WRITE(NOUT,BXTXT) '*         AUTHORS               *'
      WRITE(NOUT,BXTXT) '* .......... E. Richter-Was     *'
      WRITE(NOUT,BXTXT) '* ......................        *'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)
C
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '   ===== Z0DEC2          ======   '
      WRITE(NOUT,BXTXT) '   initialization starts....     '

      WRITE(NOUT,BXL1I) KEYRND,     ' VARRAN    switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYBRE,     ' bremsstr. switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYPRO,     ' propagator switch ','KEYPRO','  '
      WRITE(NOUT,BXL1F) XPAR(1),    ' CMSENE [GeV]      ','CMSENE','X1'
      WRITE(NOUT,BXL1F) AMFIN,      ' quark  mass [GeV] ','AMFIN ','X1'
      WRITE(NOUT,BXL1F) AMINI,      ' lepton mass [GeV] ','AMINI ','X1'
      WRITE(NOUT,BXL1F) AMAZ ,      ' Z_0 mass [GeV]    ','AMAZ  ','X1'
      WRITE(NOUT,BXL1F) GAMMZ,      ' Z_0 width   [GeV] ','GAMMZ ','X1'
      WRITE(NOUT,BXL1F) SINW2,      ' sin theta_W^2     ','SINW2 ','X1'
      WRITE(NOUT,BXL1F) XMIN ,      ' kinemat XMIN      ','XMIN  ','X1'
      WRITE(NOUT,BXL1F) PTMIN,      ' gener min   P_T   ','PTMIN ','X1'
      WRITE(NOUT,BXL1F) PHCUT,      ' photon min energy ','PHCUT ','X1'
      WRITE(NOUT,BXL1F) PHMAX,      ' photon max energy ','PHMAX ','X1'
      WRITE(NOUT,BXTXT) '   defining trigger             '
      WRITE(NOUT,BXL1F) PHMIN,      ' photon min   P_T  ','PTMIN ','X1'       
      WRITE(NOUT,BXL1F) PTL ,       ' l/l_bar minP_T    ','PTL   ','X1'
      WRITE(NOUT,BXL1F) XMIMAS,     ' l/l_bar min mass  ','XMIMAS','X1'
       
      WRITE(NOUT,BXCLO)

C****      CALL VESK2W(-1,FUNSKI,DUM1,DUM2,WT)
c weight monitoring initialization
      CALL WMONIT(-1,1,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,2,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,3,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,11,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,12,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,13,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,21,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,22,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,23,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,50,DUMM1,DUMM2,DUMM3)


      ELSEIF(MODE.EQ.0) THEN
C     =======================
C generating x1,x2 according to crude distribution (1/x1**2)*(1/x2**2) in range
C (x_min,1) and (svar*x1*x2) > 4*amtop**2
 
C****      CALL VESK2W(0,FUNSKI,DUM1,DUM2,WTVES)
C****      IF(WTVES.EQ.0D0) THEN
C****        WTMOD1=0D0
C****        WTMOD2=0D0
C****        WTMOD3=0D0
C****      ELSE
C****        X1 = Z1
C****        X2 = Z2
        WTVES = 1D0
        X1 = 1D0
        X2 = 1D0
C crude x1,x2 distribution
C***        DSCRUD = 1D0/X1**2/X2**2
        DSCRUD = 1D0
C reduced frame energy
        XMSENE = DSQRT(SVAR*X1*X2)
C calculating glu-glu vectors in reduced frame
        CALL BEAMAS(XMSENE,AMINI,P1,Q1)
C generating tree-body (top,top,photon) phase space and vectors in reduced frame
        AMPHOT=0D0
        IF(KEYBRE.EQ.3)
     #  CALL KINE4F(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
        IF(KEYBRE.EQ.2)
     #  CALL KINE4IFZ0(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
        IF(KEYBRE.EQ.1)
     #  CALL KINE4IIZ0(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
c.....control integral on the phase space generator
        CALL WMONIT(0,50,WTKIN,1D0,0D0)

        YMSENE=DSQRT( (P2(4)+Q2(4))**2-(P2(3)+Q2(3))**2
     #               -(P2(2)+Q2(2))**2-(P2(1)+Q2(1))**2   )
C momenta transformation to laboratory frame
        EXE =X1/DSQRT(X1*X2)
        CALL BOSTD3(EXE, P1, P1L)
        CALL BOSTD3(EXE, Q1, Q1L)
        CALL BOSTD3(EXE, P2, P2L)
        CALL BOSTD3(EXE, Q2, Q2L)
        CALL BOSTD3(EXE,PHOT1,PHOT1L)
        CALL BOSTD3(EXE,PHOT2,PHOT2L)
C rejecting infrared photons, but laboratory frame
        EPH1 = PHOT1L(4)
        EPH2 = PHOT2L(4)
        IF(EPH1.LT.PHCUT.OR.EPH1.GT.PHMAX) WTKIN = 0D0
        IF(EPH2.LT.PHCUT.OR.EPH2.GT.PHMAX) WTKIN = 0D0
C rejecting photons in proton-proton CMS frame
        PT1L=DSQRT(PHOT1L(1)**2+PHOT1L(2)**2)
        PT2L=DSQRT(PHOT2L(1)**2+PHOT2L(2)**2)
        IF(PT1L.LT.PHMIN.OR.PT2L.LT.PHMIN) WTKIN = 0D0
C rejecting top in proton-proton CMS frame
        PT=DSQRT(P2L(1)**2+P2L(2)**2)
        IF(PT.LT.PTL) WTKIN = 0D0
        PT=DSQRT(Q2L(1)**2+Q2L(2)**2)
        IF(PT.LT.PTL) WTKIN = 0D0
C checking lepton-lepton mass
        IF(YMSENE.LE.XMIMAS) WTKIN=0D0
C avoiding matrix element calculation
        IF(WTKIN.EQ.0D0) GOTO 100
C calculating Mandelstam variable in reduced frame
        SDOT=XMSENE**2
C calculating energy scale Q2 for str. func, alfa strong...
        QSQR = SDOT
C....functions  ULALPS ,ULALEM from PYTHIA library
ccc        ALPS = DBLE(ULALPS(REAL(QSQR)))
ccc        ALEM = DBLE(ULALEM(REAL(QSQR)))
C calculating gluon distributions from PYTHIA library
C***        CALL PYSTPR(REAL(X1),REAL(QSQR),XPPR)
C***        XQUA1=DBLE(XPPR(0))
C***        CALL PYSTPR(REAL(X2),REAL(QSQR),XPPR)
C***        XQUA2=DBLE(XPPR(0))
C  structure functions f(x_1,QSQR)*f(x_2,QSQR) distribution
C***         DISTRY = XQUA1*XQUA2/X1/X2
         DISTRY = 1D0
C   cross section  for q qbar --> Z0 --> l lbar photon
C   according to exact matrix element calculated numerically with help of spin
C   amplitude technique
C and compact formulas from YFS3
C quark up in initial state
       CALL INTEFB
       
       XNORM= ALFA**2*(4D0*PI*ALFA)**2 
       IFLEV=0 !
       write(*,*) 'iflev=',iflev,'; zero defns. beams charge 1'
C tutaj amplitudy w trzech wersjach 
C  DINIINF --podczerwony
C  DINIAPR --YFS3
C  DUBLINI --spinowe
       IF(KEYBRE.EQ.1) THEN
         CALL DINIINF(SECT1)
         XCROS1 =SECT1 *XNORM *(2D0/3D0)**4
         CALL DINIAPR(SECT2)
         XCROS2 =SECT2 *XNORM *(2D0/3D0)**4
         CALL DUBLINI(SECT3)
         XCROS3 =SECT3 *XNORM *(2D0/3D0)**4
 
       ELSEIF(KEYBRE.EQ.3) THEN
         CALL DFININF(SECT11)
         XCROS11=SECT11*XNORM
         CALL DFINAPR(SECT12)
         XCROS12=SECT12*XNORM
         CALL DUBLFIN(SECT13)
         XCROS13=SECT13*XNORM
      ELSEIF(KEYBRE.EQ.2) THEN
c extra factor *2D0 because the conbinatorial factor 1/2 from
c the phase space shoul be cancel, the are two photons
c but matrix element in principle shoul be averaged over them
c it is not easy because of the generation presampling
c is it O.K.   !!!!!!?????????!!!!!!!!!!!
         CALL DMIXINF(SECT21)
         XCROS21=SECT21*XNORM *(2D0/3D0)**2 *2D0
         CALL DMIXAPR(SECT22)
         XCROS22=SECT22*XNORM *(2D0/3D0)**2 *2D0
         CALL DUBLMIX(SECT23)
         XCROS23=SECT23*XNORM *(2D0/3D0)**2 *2D0
        
       ENDIF




C   normalization factor (comming from the Born cross section callibration
C   with PYTHIA)
         COMFAC=PI/SDOT * (8D0*PI) 
C...matrix element not calculated for WTKIN=0D0
 

 100     CONTINUE

C final weights
      IF(KEYBRE.EQ.1) THEN
C INITIAL STATE BREMSSTRAHLUNG
         WTMOD1 = WTVES*WTKIN*DISTRY*XCROS1*COMFAC/DSCRUD
         WTMOD2 = WTVES*WTKIN*DISTRY*XCROS2*COMFAC/DSCRUD
         WTMOD3 = WTVES*WTKIN*DISTRY*XCROS3*COMFAC/DSCRUD
      ELSEIF(KEYBRE.EQ.3) THEN
C FINAL STATE BREMSSTRAHLUNG
         WTMOD11= WTVES*WTKIN*DISTRY*XCROS11*COMFAC/DSCRUD
         WTMOD12= WTVES*WTKIN*DISTRY*XCROS12*COMFAC/DSCRUD
         WTMOD13= WTVES*WTKIN*DISTRY*XCROS13*COMFAC/DSCRUD
      ELSEIF(KEYBRE.EQ.2) THEN
C INITIAL-FINAL STATE BREMSSTRAHLUNG
         WTMOD21= WTVES*WTKIN*DISTRY*XCROS21*COMFAC/DSCRUD
         WTMOD22= WTVES*WTKIN*DISTRY*XCROS22*COMFAC/DSCRUD
         WTMOD23= WTVES*WTKIN*DISTRY*XCROS23*COMFAC/DSCRUD
      ENDIF


         IF(KEYBRE.EQ.1) THEN
             WTMOD=WTMOD3
             WTINF=WTMOD1
             WTYFS=WTMOD2
             WTEXA=WTMOD3
         ELSEIF(KEYBRE.EQ.2) THEN
             WTMOD=WTMOD23
             WTINF=WTMOD21
             WTYFS=WTMOD22
             WTEXA=WTMOD23
         ELSEIF(KEYBRE.EQ.3) THEN
             WTMOD=WTMOD13
             WTINF=WTMOD11
             WTYFS=WTMOD12
             WTEXA=WTMOD13
         ENDIF


         IF(KEYBRE.EQ.3) WTMOD=WTMOD13
         IF(KEYBRE.EQ.2) WTMOD=WTMOD23

C***      ENDIF
c weight monitoring 
      CALL WMONIT(0,1,WTMOD1,1D0,0D0)
      CALL WMONIT(0,2,WTMOD2,1D0,0D0)
      CALL WMONIT(0,3,WTMOD3,1D0,0D0)
      CALL WMONIT(0,11,WTMOD11,1D0,0D0)
      CALL WMONIT(0,12,WTMOD12,1D0,0D0)
      CALL WMONIT(0,13,WTMOD13,1D0,0D0)
      CALL WMONIT(0,21,WTMOD21,1D0,0D0)
      CALL WMONIT(0,22,WTMOD22,1D0,0D0)
      CALL WMONIT(0,23,WTMOD23,1D0,0D0)



      ELSEIF(MODE.EQ.1) THEN
C     ======================= 

      NEVGEN = NPAR(10)
C***      CALL VESK2W(1,FUNSKI,AWT,EREL,XCRU)
C***      DWT      = AWT*EREL
      XCRU = 1D0
      DWT  = 1D0
      AWT  = 1D0
      XPAR(10) = XCRU*GNANOB 
      CALL WMONIT(1,1,AWT1,DWT1,0D0)
      XSEC1=AWT1*XCRU*GNANOB
      DXSEC1=XSEC1*DWT1
      CALL WMONIT(1,2,AWT2,DWT2,0D0)
      XSEC2=AWT2*XCRU*GNANOB
      DXSEC2=XSEC2*DWT2
      CALL WMONIT(1,3,AWT3,DWT3,0D0)
      XSEC3=AWT3*XCRU*GNANOB
      DXSEC3=XSEC3*DWT3
      CALL WMONIT(1,11,AWT11,DWT11,0D0)
      XSEC11=AWT11*XCRU*GNANOB
      DXSEC11=XSEC11*DWT11
      CALL WMONIT(1,12,AWT12,DWT12,0D0)
      XSEC12=AWT12*XCRU*GNANOB
      DXSEC12=XSEC12*DWT12
      CALL WMONIT(1,13,AWT13,DWT13,0D0)
      XSEC13=AWT13*XCRU*GNANOB
      DXSEC13=XSEC13*DWT13
      CALL WMONIT(1,21,AWT21,DWT21,0D0)
      XSEC21=AWT21*XCRU*GNANOB
      DXSEC21=XSEC21*DWT21
      CALL WMONIT(1,22,AWT22,DWT22,0D0)
      XSEC22=AWT22*XCRU*GNANOB
      DXSEC22=XSEC22*DWT22
      CALL WMONIT(1,23,AWT23,DWT23,0D0)
      XSEC23=AWT23*XCRU*GNANOB
      DXSEC23=XSEC23*DWT23

      CALL WMONIT(1,50,AWT50,DWT50,DUM3)

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC2       : WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '   X.sect. in [nb] units         '
      WRITE(NOUT,BXTXT) '   for total generated sample    '
      WRITE(NOUT,BXL1I) NEVGEN,     'generated events   ','NEVGEN','A1'
      WRITE(NOUT,BXL2F) AWT50 ,DWT50   ,'phase space    ','AWT50 ','A2'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC2       : INITIAL        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) 'according to formula in infrared limit YFS3'
      WRITE(NOUT,BXL2F) XSEC1,DXSEC1,'xsection 1   [nb]  ','XSEC4','A2'   
      WRITE(NOUT,BXTXT) 'according to formula in YFS3'
      WRITE(NOUT,BXL2F) XSEC2,DXSEC2,'xsection     [nb]  ','XSEC2','A2'   
      WRITE(NOUT,BXTXT) 'according to formula from spin amplitudes'
      WRITE(NOUT,BXL2F) XSEC3,DXSEC3,'xsection     [nb]  ','XSEC3','A2'   
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC2       : FINAL          '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) 'according to formula in infrared limit YFS3'
      WRITE(NOUT,BXL2F) XSEC11,DXSEC11,'xsection    [nb] ','XSEC11','A2'   
      WRITE(NOUT,BXTXT) 'according to formula in YFS3'
      WRITE(NOUT,BXL2F) XSEC12,DXSEC12,'xsection    [nb] ','XSEC12','A2'   
      WRITE(NOUT,BXTXT) 'according to formula from spin  amplitudes'
      WRITE(NOUT,BXL2F) XSEC13,DXSEC13,'xsection    [nb] ','XSEC13','A2'   
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC2       : INITIAL- FINAL '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) 'according to formula in infrared limit YFS3'
      WRITE(NOUT,BXL2F) XSEC21,DXSEC21,'xsection    [nb] ','XSEC21','A2'   
      WRITE(NOUT,BXTXT) 'according to formula in YFS3'
      WRITE(NOUT,BXL2F) XSEC22,DXSEC22,'xsection    [nb] ','XSEC22','A2'   
      WRITE(NOUT,BXTXT) 'according to formula from spin  amplitudes'
      WRITE(NOUT,BXL2F) XSEC23,DXSEC23,'xsection    [nb] ','XSEC23','A2'   
      WRITE(NOUT,BXCLO)


      ENDIF
C     =====      
      END

       SUBROUTINE 
     #     KINE4C(AMTAU,AMP3,AMP2,AMP1,AMNUTA,PIM3,PIM2,PIM1,PN,WT)
C     *******************************************************************
C generator of 4 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA  - masses of particles
C         PIM1,PIM2,PIM3,PN  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C    factor 1/2 for two identical particles included
C         WT  - weight
C subroutine is based on subroutine DPHTRE from TAUOLA
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PR(4),PAA(4),PT(4)
      DIMENSION RRR(6)
      COMMON /POMOC2/ Y,Z,DELZ,DEL1,DEL2,BETA,C,SPHI,C1,S1,CJAC

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1.D0/2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAU
C
      CALL VARRAN(RRR,6)
C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C MASS OF (REAL/VIRTUAL) A1 -> ( RHO + PIM1) flat phase space
C>           AMS1=(AMP1+AMP2+AMP3)**2
C>           AMS2=(AMTAU-AMNUTA)**2
C>           AM3SQ=AMS1+   RRR(1)*(AMS2-AMS1)
C>           AM3 =SQRT(AM3SQ)
C>           PHSPAC = PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PN (y generated with density 1/y)
        AMS1=(AMP1+AMP2+AMP3)**2       
        Y1  = (AMTAU**2-AMS1+AMNUTA**2)/AMTAU**2
        Y0  = 0.001D0
        IF(Y1.LT.Y0) Y1=Y0
        YL1 = DLOG(Y1/Y0)
        YL0 = DLOG(Y0)
        Y   = DEXP(YL1*RRR(1)+YL0)
        AM3SQ=(1D0-Y)*AMTAU**2+AMNUTA**2
        AM3=SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMTAU**2*YL1*Y
        IF(PHSPAC.EQ.0D0) GOTO 900
C MASS OF (REAL/VIRTUAL) RHO -> (PIM2+PIM3) flat phase space
C>        AMS1=(AMP2+AMP3)**2
C>        AMS2=(AM3-AMP1)**2
C>        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
C>        AM2 =SQRT(AM2SQ)
C>        PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PIM1 (y generated with density 1/y)
        AMS1=(AMP2+AMP3)**2
        Z1  = (AM3**2-AMS1+AMP1**2)/AM3**2
        Z0  = 0.001D0
        IF(Z1.LT.Z0) Z1=Z0
        ZL1 = DLOG(Z1/Z0)
        ZL0 = DLOG(Z0)
        Z   = DEXP(ZL1*RRR(2)+ZL0)
        AM2SQ=(1D0-Z)*AM3**2+AMP1**2
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*AM3**2*ZL1*Z
        IF(PHSPAC.EQ.0D0) GOTO 900
* RHO RESTFRAME, DEFINE PIPL AND PIM1
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PIM3)
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
* A1 REST FRAME, DEFINE PIM1
*       RHO  MOMENTUM
        PR(1)=0.D0
        PR(2)=0.D0
        PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
        PIM1(1)=0.D0
        PIM1(2)=0.D0
        PIM1(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMP1**2)
        PIM1(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
* ALL PIONS AND RHO ROTATED IN THE A1 REST FRAME
      THET =ACOS(-1.D0+2*RRR(3))
      PHI = 2*PI*RRR(4)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0.D0
      PN(2)=0.D0
      PN(4)=1.D0/(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PR,PR)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(5))
      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END
      SUBROUTINE 
     #     KINE4F(AMTAU,AMP3,AMP2,AMP1,AMNUTA,PIM3,PIM2,PIM1,PN,WT)
C     *******************************************************************
C generator of 4 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA  - masses of particles
C         PIM1,PIM2,PIM3,PN  - generated momenta 
C         presampling on collinear and infrared singularity
C         for PN, PIM1 momenta
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on collinear singlularity in final state and infrared singularity
C subroutine is based on subroutine DPHTRE from TAUOLA
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PR(4),PAA(4),PT(4)
      DIMENSION RRR(20)
      DIMENSION PA(4), PB(4)
      common /ala/ca,xjaca,cb,xjacb,cc,c1,c2,f1,f2,g1,g2

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1.D0/2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAU
C
      CALL VARRAN(RRR,20)
C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C MASS OF (REAL/VIRTUAL) A1 -> ( RHO + PIM1) flat phase space
C>           AMS1=(AMP1+AMP2+AMP3)**2
C>           AMS2=(AMTAU-AMNUTA)**2
C>           AM3SQ=AMS1+   RRR(1)*(AMS2-AMS1)
C>           AM3 =SQRT(AM3SQ)
C>           PHSPAC = PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PN (y generated with density 1/y)
        AMS1=(AMP1+AMP2+AMP3)**2       
        Y1  = (AMTAU**2-AMS1+AMNUTA**2)/AMTAU**2
        Y0  = 0.01D0
        IF(Y1.LT.Y0) Y1=Y0
        YL1 = DLOG(Y1/Y0)
        YL0 = DLOG(Y0)
        Y   = DEXP(YL1*RRR(1)+YL0)
        AM3SQ=(1D0-Y)*AMTAU**2+AMNUTA**2
        AM3=SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMTAU**2*YL1*Y
        IF(PHSPAC.EQ.0D0) GOTO 900
C MASS OF (REAL/VIRTUAL) RHO -> (PIM2+PIM3) flat phase space
C>        AMS1=(AMP2+AMP3)**2
C>        AMS2=(AM3-AMP1)**2
C>        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
C>        AM2 =SQRT(AM2SQ)
C>        PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PIM1 (y generated with density 1/y)
        AMS1=(AMP2+AMP3)**2
        Z1  = (AM3**2-AMS1+AMP1**2)/AM3**2
        Z0  = 0.01D0
        IF(Z1.LT.Z0) Z1=Z0
        ZL1 = DLOG(Z1/Z0)
        ZL0 = DLOG(Z0)
        Z   = DEXP(ZL1*RRR(2)+ZL0)
        AM2SQ=(1D0-Z)*AM3**2+AMP1**2
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*AM3**2*ZL1*Z
        IF(PHSPAC.EQ.0D0) GOTO 900
* RHO RESTFRAME, DEFINE PIPL AND PIM1
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        PHI1=2*PI*RRR(7)
        EPS=4*(AMP2/AMTAU)**2
C>>>>
c        EPS=4*(AMP2/AM2)**2
C>>>>
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
        ETA  =EXP(XL1*RRR(6)+XL0)
        CTHET1=1+EPS-ETA
        IF (RRR(8).GT.0.5) THEN
          CTHET1=-CTHET1
          IFLAG1=2
        ELSE
          IFLAG1=1
        ENDIF
        CA=CTHET1
        XJACA=XL1/2*ETA
        THET =ACOS(CTHET)
C this will be recovered later
C        PHSPAC=PHSPAC*XL1/2*ETA
C                            CTHET1=-1.0+2.0*RRR(6)
        STHET1=SQRT(1-CTHET1**2)
        PIM3(1)=PPPI*STHET1*SIN(PHI1)
        PIM3(2)=PPPI*STHET1*COS(PHI1)
        PIM3(3)=PPPI*CTHET1
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
* A1 REST FRAME, DEFINE PIM1
*       RHO  MOMENTUM
        PR(1)=0.D0
        PR(2)=0.D0
        PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
        PIM1(1)=0.D0
        PIM1(2)=0.D0
        PIM1(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMP1**2)
        PIM1(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
* ALL PIONS AND RHO ROTATED IN THE A1 REST FRAME
      IF(RRR(9).GT.0.5D0) THEN
       IFLAG2=1
C      PIM3 along 3-rd axix
      GHI=-PHAN1E(PIM3(1),PIM3(2))
      CALL PHRO3E(GHI,PIM3)
      CALL PHRO3E(GHI,PIM1)
      CALL PHRO3E(GHI,PIM2)
      CALL PHRO3E(GHI,PR)
      GHET=PI-PHAN2E(PIM3(3),SQRT(PIM3(1)**2+PIM3(2)**2))
      CALL PHRO2E(GHET,PIM3)
      CALL PHRO2E(GHET,PIM1)
      CALL PHRO2E(GHET,PIM2)
      CALL PHRO2E(GHET,PR)
      FIK=2.*PI*RRR(10)
      CALL PHRO3E(FIK,PIM1)
      CALL PHRO3E(FIK,PIM2)
      ELSE
       IFLAG2=2
C      PIM2 along third axix
      GHI=-PHAN1E(PIM2(1),PIM2(2))
      CALL PHRO3E(GHI,PIM3)
      CALL PHRO3E(GHI,PIM1)
      CALL PHRO3E(GHI,PIM2)
      CALL PHRO3E(GHI,PR)
      GHET=PI-PHAN2E(PIM2(3),SQRT(PIM2(1)**2+PIM2(2)**2))
      CALL PHRO2E(GHET,PIM3)
      CALL PHRO2E(GHET,PIM1)
      CALL PHRO2E(GHET,PIM2)
      CALL PHRO2E(GHET,PR)
      FIK=2.*PI*RRR(10)
      CALL PHRO3E(FIK,PIM3)
      CALL PHRO3E(FIK,PIM1)
      ENDIF
      RR3=RRR(3)
      RR4=RRR(4)
CAM   THET =PI*RR3
        PHI=2*PI*RRR(4)
        EPS=4*(AMP2/AMTAU)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
        ETA  =EXP(XL1*RRR(3)+XL0)
        CTHET=1+EPS-ETA
        THET =ACOS(CTHET)
        CB=CTHET
        XJACB=XL1/2*ETA
c        write(6,*) 'uffi',CTHET,'  ',XL1/2*ETA
C this will be recovered later
C        PHSPAC=PHSPAC*XL1/2*ETA
C      THET =ACOS(-1.+2*RR3)
C      PHI = 2*PI*RR4
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0.D0
      PN(2)=0.D0
      PN(4)=1.D0/(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PR,PR)
      DO K=1,4
       PA(K)=PIM1(K)+PIM2(K)
       PB(K)=PIM3(K)+PIM1(K)
      ENDDO
        EPS=4*(AMP2/AMTAU)**2
        XL1=LOG((2+EPS)/EPS)
C>>>>>
c        EPSA=4*(AMP2/AM2)**2
c        XL1A=LOG((2+EPSA)/EPSA)
C>>>>>
      CALL COSIK(PIM3,PIM2,PIM1,CC)
      CALL COSIK(PIM3,PA,PN,C1)
      CALL COSIK(PIM2,PB,PN,C2)
        F1=XL1/2*(1+EPS-CC)
        F2=XL1/2*(1+EPS+CC)
C>>>>
c        F1=XL1A/2*(1+EPSA-CC)
c        F2=XL1A/2*(1+EPSA+CC)
C>>>>
C        F1=1.0
C        F2=1.0
        G1=XL1/2*(1+EPS+C1)
        G2=XL1/2*(1+EPS+C2)
C        G1=1.0
C        G2=1.0
        PHSPAC=PHSPAC*4.D0/(1.D0/F1+1.D0/F2)/(1.D0/G1+1.D0/G2)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(11))
      PHI = 2*PI*RRR(12)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)

C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END
      FUNCTION PHAN1E(X,Y)
C     *********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=3.14159265358979324D0
      TWOPI=6.28318530717958648D0
      IF (ABS(Y).LT.ABS(X)) THEN
        PHAN1E=ATAN(ABS(Y/X))
        IF (X.LE.0.D0) PHAN1E=PI-PHAN1E
      ELSE
        PHAN1E=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF (Y.LT.0.D0) PHAN1E=TWOPI-PHAN1E
      RETURN

      END
      FUNCTION PHAN2E(X,Y)
C     ********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=3.14159265358979324D0
      TWOPI=6.28318530717958648D0

      IF (ABS(Y).LT.ABS(X)) THEN
        PHAN2E=ATAN(ABS(Y/X))
        IF (X.LE.0.D0) PHAN2E=PI-PHAN2E
      ELSE
        PHAN2E=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      RETURN

      END
      SUBROUTINE PHRO2E(ANGLE,PVEC)
C     *****************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4)

      CS=COS(ANGLE)*PVEC(1)+SIN(ANGLE)*PVEC(3)
      SN=-SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(3)
      PVEC(1)=CS
      PVEC(3)=SN
      RETURN

      END
      SUBROUTINE PHRO3E(ANGLE,PVEC)
C     ******************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4)

      CS=COS(ANGLE)*PVEC(1)-SIN(ANGLE)*PVEC(2)
      SN=SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(2)
      PVEC(1)=CS
      PVEC(2)=SN
      RETURN

      END
      FUNCTION XMUL(PP,QQ)
C     *********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PP(4),QQ(4)

      XMUL=PP(4)*QQ(4)-PP(3)*QQ(3)-PP(2)*QQ(2)-PP(1)*QQ(1)

      END
      SUBROUTINE COSIK(PA,PB,PH,COSTHE)
C     ********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PA(4),PB(4),PH(4),PAMB(4),PAPB(4)
      XLAM(X,Y,Z)=SQRT(ABS((X-Y-Z)**2-4.0*Y*Z))

      DO K=1,4
        PAMB(K)=PA(K)-PB(K)
        PAPB(K)=PA(K)+PB(K)
      ENDDO
      S=XMUL(PAPB,PAPB)
      XMA2=XMUL(PA,PA)
      XMB2=XMUL(PB,PB)
      P=(XLAM(S,XMA2,XMB2))/SQRT(S)
      E1=SQRT(P**2/4+XMA2)
      E2=SQRT(P**2/4+XMB2)
      COSTHE=((XMUL(PAMB,PH) - (E1-E2)/SQRT(S)*XMUL(PAPB,PH)))
     $        /XMUL(PAPB,PH)/P*SQRT(S)
      END
      SUBROUTINE 
     #     KINE4IIZ0(AMTAU,AMP3,AMP2,AMP1,AMNUTA,PIM3,PIM2,PIM1,PN,WT)
C     *******************************************************************
C generator of 4 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA  - masses of particles
C         PIM1,PIM2,PIM3,PN  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(10)
      COMMON /POMOC2/ Y,Z,DELZ,DEL1,DEL2,BETA,C,SPHI,C1,S1,CJAC

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1.D0/2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAU
C
      CALL VARRAN(RRR,10)
C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in (PIM1+PIM2)**2 
C GENERATING MASS of l+l- pair
        AMS1=(AMP2+AMP3)**2       
c...        AMS2=AMTAU**2-(AMP1+AMNUTA)**2 
c...  xk0 defines min. photon energy in photon-photon rest frame
        XK0=0.01D0 
        AMS2=AMTAU**2-(2D0*XK0*AMTAU)**2       
        Y1  = 1D0-AMS1/AMTAU**2
        Y0  = 1D0-AMS2/AMTAU**2
        IF(Y1.LT.Y0) Y1=Y0
        AZ0=AMAZ/AMTAU
        BZ0=AMAZ*GAMMZ/AMTAU**2
        XP2=1D0/BZ0*
     #   (ATAN((1D0-Y0-AZ0**2)/BZ0)-ATAN((1D0-Y1-AZ0**2)/BZ0))
        XP3=DLOG((1D0-Y0)/(1D0-Y1))
        IF(RRR(9).LE.(XP2)/(XP2+XP3)) THEN
          Y=1D0-AZ0**2-BZ0*TAN(
     #            +RRR(1) *ATAN((1D0-Y1-AZ0**2)/BZ0)
     #      +(1D0 -RRR(1))*ATAN((1D0-Y0-AZ0**2)/BZ0)
     #                       )
        ELSEIF(RRR(9).GT.(XP2)/(XP2+XP3)) THEN
           Y=1D0-(1D0-Y0)*((1D0-Y1)/(1D0-Y0))**RRR(1)
        ENDIF
        XJAC2=ABS(BZ0*( ATAN((1D0-Y1-AZ0**2)/BZ0)
     #           -ATAN((1D0-Y0-AZ0**2)/BZ0)   )
     #         /COS(ATAN((1D0-Y-AZ0**2)/BZ0))**2)
        XJAC3=ABS((1D0-Y)*DLOG((1D0-Y0)/(1D0-Y1)))
        XJAC=(XP2+XP3)/(XP2/XJAC2+XP3/XJAC3)
        AM3SQ=(1D0-Y)*AMTAU**2
        AM3=SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMTAU**2*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900
* AM3 RESTFRAME, DEFINE PIM2 AND PIM3
        ENQ1=(AM3SQ-AMP2**2+AMP3**2)/(2*AM3)
        ENQ2=(AM3SQ+AMP2**2-AMP3**2)/(2*AM3)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM3)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PIM3)
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
C MASS OF gam-gam pair 
        AMS1=(AMP1+AMNUTA)**2
        AMS2=(AMTAU-AM3)**2
        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
        PHSPAC=PHSPAC*(AMS2-AMS1)
* AM2 RESTFRAME, DEFINE PIM1 AND PN
        ENQ1=(AM2SQ-AMP1**2+AMNUTA**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP1**2-AMNUTA**2)/(2*AM2)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERD(PPPI,PN)
        PN(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 40 I=1,3
 40     PIM1(I)=-PN(I)
        PIM1(4)= ENQ2
* NOW TO THE TAU REST FRAME, DEFINE AM2 AND AM3 MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AM2**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
      PBB(1)=0.D0
      PBB(2)=0.D0
      PBB(4)=1.D0/(2*AMTAU)*(AMTAU**2+AM2**2-AM3**2)
      PBB(3)=-PAA(3)
* ALL PIONS BOOSTED  TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO AM2 MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      EXE=(PBB(4)+PBB(3))/AM2
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(5))
      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END
      SUBROUTINE 
     #     KINE4IFZ0(AMTAU,AMP3,AMP2,AMP1,AMNUTA,PIM3,PIM2,PIM1,PN,WT)
C     *******************************************************************
C generator of 4 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA  - masses of particles
C         PIM1,PIM2,PIM3,PN  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on collinear singularity in initial/final state
C infrared singularity and resonance in (pim3+pim2+pim1)
C subroutine is based on subroutine DPHTRE from TAUOLA
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PR(4),PAA(4),PT(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(15)


C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1.D0/2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAU
C
      CALL VARRAN(RRR,15)
C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C MASS OF (REAL/VIRTUAL) A1 -> ( RHO + PIM1) flat phase space
C>           AMS1=(AMP1+AMP2+AMP3)**2
C>           AMS2=(AMTAU-AMNUTA)**2
C>           AM3SQ=AMS1+   RRR(1)*(AMS2-AMS1)
C>           AM3 =SQRT(AM3SQ)
C>           PHSPAC = PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH PRESAMPLING ON INFRARED SINGULARITY FOR PN 
C and resonance in AM3
C>>>        AMS1=(AMP1+AMP2+AMP3)**2 
c...for initial state radiation
        AMS1=MAX(4D0*AMFIN**2,4D0*AMINI**2)      
        Y1  = (AMTAU**2-AMS1+AMNUTA**2)/AMTAU**2
        Y0  = 0.01D0
        IF(Y1.LT.Y0) Y1=Y0
        AZ0=AMAZ/AMTAU
        BZ0=AMAZ*GAMMZ/AMTAU**2
        XP1=DLOG(Y1/Y0)
        XP2=1D0/BZ0*
     #   (ATAN((1D0-Y0-AZ0**2)/BZ0)-ATAN((1D0-Y1-AZ0**2)/BZ0))
        XP3=DLOG((1D0-Y0)/(1D0-Y1))
        IF(RRR(9).LT.(XP1/(XP1+XP2+XP3)) ) THEN
           Y=Y0*(Y1/Y0)**RRR(1)
        ELSEIF(RRR(9).LT.(XP1+XP2)/(XP1+XP2+XP3)) THEN
          Y=1D0-AZ0**2-BZ0*TAN(
     #            +RRR(1) *ATAN((1D0-Y1-AZ0**2)/BZ0)
     #      +(1D0 -RRR(1))*ATAN((1D0-Y0-AZ0**2)/BZ0)
     #                       )
        ELSEIF(RRR(9).GT.(XP1+XP2)/(XP1+XP2+XP3)) THEN
           Y=1D0-(1D0-Y0)*((1D0-Y1)/(1D0-Y0))**RRR(1)
        ENDIF
        XJAC1=ABS(Y*DLOG(Y1/Y0))
        XJAC2=ABS(BZ0*( ATAN((1D0-Y1-AZ0**2)/BZ0)
     #           -ATAN((1D0-Y0-AZ0**2)/BZ0)   )
     #         /COS(ATAN((1D0-Y-AZ0**2)/BZ0))**2)
        XJAC3=ABS((1D0-Y)*DLOG((1D0-Y0)/(1D0-Y1)))
        XJAC=(XP1+XP2+XP3)/(XP1/XJAC1+XP2/XJAC2+XP3/XJAC3)
        AM3SQ=(1D0-Y)*AMTAU**2+AMNUTA**2
        AM3=SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMTAU**2*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900
C MASS OF (REAL/VIRTUAL) RHO -> (PIM2+PIM3) flat phase space
C>        AMS1=(AMP2+AMP3)**2
C>        AMS2=(AM3-AMP1)**2
C>        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
C>        AM2 =SQRT(AM2SQ)
C>        PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PIM1 (y generated with density 1/y)
        AMS1=(AMP2+AMP3)**2
        Z1  = (AM3**2-AMS1+AMP1**2)/AM3**2
        Z0  = 0.01D0
        IF(Z1.LT.Z0) Z1=Z0
        ZL1 = DLOG(Z1/Z0)
        ZL0 = DLOG(Z0)
        Z   = DEXP(ZL1*RRR(2)+ZL0)
        AM2SQ=(1D0-Z)*AM3**2+AMP1**2
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*AM3**2*ZL1*Z
        IF(PHSPAC.EQ.0D0) GOTO 900
* RHO RESTFRAME, DEFINE PIPL AND PIM1
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
C>>>        CALL SPHERD(PPPI,PIM3)
C>>>        PIM3(4)=ENQ1
C ANGULAR DISTRIBUTION WITH PRESAMPLING ON COLLINEAR SINGULARITY
C IN FINAL STATE
      AMM2=4D0*AMP3**2/AM2SQ
      BETA=SQRT(1.D0-AMM2)                
      EPS=AMM2/(1.D0+SQRT(1.D0-AMM2))      
      DEL1=(2.D0-EPS)*(EPS/(2.D0-EPS))**RRR(3) 
      DELL = DEL1
      DEL2=2.D0-DEL1 
      XJAK=-1D0/2D0/BETA*DLOG(EPS/(2D0-EPS))*DELL*(2D0-DELL)
C SYMMETRIZATION                         
      IF(RRR(4).LE.0.5D0) THEN              
        A=DEL1                           
        DEL1=DEL2                        
        DEL2=A                           
      ENDIF                              
C CALCULATION OF SIN AND COS THETA FROM INTERNAL VARIABLES 
      COSTHG=(1.D0-DEL1)/BETA            
      SINTHG=DSQRT(DEL1*DEL2-AMM2)/BETA 
      PHSPAC=PHSPAC*XJAK
      PHI = 2D0*PI*RRR(5)
      PIM3(1)=PPPI*SINTHG*COS(PHI)
      PIM3(2)=PPPI*SINTHG*SIN(PHI)
      PIM3(3)=PPPI*COSTHG
      PIM3(4) = ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
* A1 REST FRAME, DEFINE PIM1
*       RHO  MOMENTUM
        PR(1)=0.D0
        PR(2)=0.D0
        PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
        PIM1(1)=0.D0
        PIM1(2)=0.D0
        PIM1(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMP1**2)
        PIM1(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
* ALL PIONS AND RHO ROTATED IN THE A1 REST FRAME
      THET =ACOS(-1.D0+2*RRR(6))
      PHI = 2*PI*RRR(7)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0.D0
      PN(2)=0.D0
      PN(4)=1.D0/(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PR,PR)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
C>>>      THET =ACOS(-1.D0+2*RRR(5))
C>>>      PHI = 2*PI*RRR(6)
C ANGULAR DISTRIBUTION WITH PRESAMPLING ON COLLINEAR SINGULARITY
C IN INITIAL STATE
      AM2=4D0*AMINI**2/AMTAU**2
      BETA=SQRT(1.D0-AM2)                
      EPS=AM2/(1.D0+SQRT(1.D0-AM2))      
      DEL1=(2.D0-EPS)*(EPS/(2.D0-EPS))**RRR(7) 
      DELL = DEL1
      DEL2=2.D0-DEL1 
      XJAK=-1D0/2D0/BETA*DLOG(EPS/(2D0-EPS))*DELL*(2D0-DELL)
C SYMMETRIZATION                         
      IF(RRR(8).LE.0.5D0) THEN              
        A=DEL1                           
        DEL1=DEL2                        
        DEL2=A                           
      ENDIF                              
C CALCULATION OF SIN AND COS THETA FROM INTERNAL VARIABLES 
      COSTHG=(1.D0-DEL1)/BETA            
      SINTHG=DSQRT(DEL1*DEL2-AM2)/BETA 
      PHSPAC=PHSPAC*XJAK
      THET =ACOS(COSTHG)
      PHI =2*PI*RRR(9)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END


