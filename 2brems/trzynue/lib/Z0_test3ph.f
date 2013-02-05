C---------------------------------------------------------------------------C
C --------    SUBROUTINES FOR KEYOPT=4  ------------------------------------C
C---------------------------------------------------------------------------C
      SUBROUTINE Z0DEC3(MODE,XPAR,NPAR)
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
      COMMON /MOMLAB2/ 
     #       PP1L(4),QQ1L(4),PP2L(4),QQ2L(4),PPHOT1L(4),PPHOT2L(4)
      COMMON /MOMCMS2/ 
     #       PP1(4) ,QQ1(4) ,PP2(4) ,QQ2(4) ,PPHOT1(4) ,PPHOT2(4)
      COMMON /MOMLAB3/ 
     #       P1L(4),Q1L(4),P2L(4),Q2L(4),PHOT1L(4),PHOT2L(4),PHOT3L(4)
      COMMON /MOMCMS3/ 
     #       P1(4) ,Q1(4) ,P2(4) ,Q2(4) ,PHOT1(4) ,PHOT2(4),PHOT3(4)

      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      COMMON /FLAVOUR/ IFLEV
      common /jakobiany/  CTHET,xccos1,xccos2,xj4,xj5,rrr
      common /jakobi1/ xjac1,xjac2,xjac3
      DIMENSION RRR(19)
      DIMENSION NPAR(99),XPAR(99)
      REAL *4 XPPR(-6:6),ULALPS 
      DIMENSION RN(1),APHOT(4)
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      EXTERNAL FUNSKI
      REAL*8 phsum(4)

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
      WRITE(NOUT,BXTXT) '*     ***   Z0DEC3   ***        *'
      WRITE(NOUT,BXTXT) '*     *****************         *'
      WRITE(NOUT,BXTXT) '*    September      1993         *'
      WRITE(NOUT,BXTXT) '*         AUTHORS               *'
      WRITE(NOUT,BXTXT) '* .......... E. Richter-Was     *'
      WRITE(NOUT,BXTXT) '* ......................        *'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)
C
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '   ===== Z0DEC3          ======   '
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
      CALL WMONIT(-1,50,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,1,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,2,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,3,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,4,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,5,DUMM1,DUMM2,DUMM3)
      CALL WMONIT(-1,6,DUMM1,DUMM2,DUMM3)
      xaus1=0
      xaus2=0

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
        AMPHOT=0.0
! niedomazka !!!!  ########################################
!        IF(KEYBRE.EQ.3)
!     #  CALL KINE4F(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
!     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
!        IF(KEYBRE.EQ.2)
!     #  CALL KINE4IFZ0(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
!     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
!        IF(KEYBRE.EQ.1)
!     #  CALL KINE4IIZ0(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,
!     #                                 P2,Q2,PHOT1,PHOT2,WTKIN)
        IF(KEYBRE.EQ.1)
     #  CALL KINE5IIZ0(XMSENE,AMFIN,AMFIN,AMPHOT,AMPHOT,AMPHOT,
     #                           P2,Q2,PHOT1,PHOT2,PHOT3,WTKIN)
      if (wtkin.lt.-10-xaus1) then
       xaus1=wtkin
       write(*,*) 'wtkin= ',wtkin
       write(*,*) keybre,wtkin
       write(*,*) ' '
       write(*,*) p2
       write(*,*) q2
       write(*,*) phot1
       write(*,*) phot2
       write(*,*) phot3
       do k=1,4
        phsum(k)=p2(k)+q2(k)+phot1(k)+phot2(k)+phot3(k)
       enddo
       write(*,*) phsum
      endif
!     stop
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
        CALL BOSTD3(EXE,PHOT3,PHOT3L)
C rejecting infrared photons, but laboratory frame
        EPH1 = PHOT1L(4)
        EPH2 = PHOT2L(4)
        EPH3 = PHOT3L(4)
        IF(EPH1.LT.PHCUT.OR.EPH1.GT.PHMAX) WTKIN = 0D0
        IF(EPH2.LT.PHCUT.OR.EPH2.GT.PHMAX) WTKIN = 0D0
        IF(EPH3.LT.PHCUT.OR.EPH3.GT.PHMAX) WTKIN = 0D0
        IF(EPH3.LT.PHCUT.OR.EPH3.GT.PHMAX) WTKIN = 0D0
C rejecting photons in proton-proton CMS frame
        PT1L=DSQRT(PHOT1L(1)**2+PHOT1L(2)**2)
        PT2L=DSQRT(PHOT2L(1)**2+PHOT2L(2)**2)
        PT3L=DSQRT(PHOT3L(1)**2+PHOT3L(2)**2)*1000000
        PT3Ll=DSQRT(PHOT3L(1)**2+PHOT3L(2)**2)
        IF(PT1L.LT.PHMIN.OR.PT2L.LT.PHMIN
     #      .OR.PT3L.LT.PHMIN) WTKIN = 0D0
C third photon must have smallest pt.
        IF(PT1L.LT.pt3ll) wtkin=0.d0
        IF(PT2L.LT.pt3ll) wtkin=0.d0
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

! we fill extra commons for mixed matrix element:  exact times soft factor.
       DO I=1,4
        PP1(I)=P1(I)
        PP2(I)=P2(I)
        QQ1(I)=Q1(I)
        QQ2(I)=Q2(I)
        PPHOT1(I)=PHOT1(I)
        PPHOT2(I)=PHOT2(I)
        PP1L(I)=P1L(I)
        PP2L(I)=P2L(I)
        QQ1L(I)=Q1L(I)
        QQ2L(I)=Q2L(I)
        PPHOT1L(I)=PHOT1L(I)
        PPHOT2L(I)=PHOT2L(I)
       ENDDO
! further common refillings, gauge fixing etc.
       CALL INTEFB
       CALL INTEFC       
       XNORM= ALFA**2*(4D0*PI*ALFA)**2 
       IFLEV=1
! m.e. caclulation normalization is `free and wild' but consitstent
!       IF(KEYBRE.EQ.1) THEN  me calculation is off !!!
       IF(KEYBRE.EQ.1) THEN
! soft as in yfs
         CALL TRZIYFSINF(SECT1)
         XCROS1 =SECT1 *XNORM *(2D0/3D0)**4
! soft as in yfs but calculated from spin ampl.
! xk=zero in numerators of propagators
         CALL TRZIINIINF(SECT2)
         XCROS2 =SECT2 *XNORM *(2D0/3D0)**4
! 1 soft 2 hard pragmatic as in yfs ph3 is soft !!!
         CALL TRZIYFS(SECT3)
         XCROS3 =SECT3 *XNORM *(2D0/3D0)**4
! 1 soft 2 hard exact spin amplitudes; ph3 is soft !!!
         CALL TRZIDBL(SECT4)
         XCROS4 =SECT4 *XNORM *(2D0/3D0)**4
! 3 hard exact
         CALL TRZYINI(SECT5)
         XCROS5 =SECT5 *XNORM *(2D0/3D0)**4
! 3 hard exact other gauge (done by INTEFD)
!         CALL INTEFD
!         CALL TRZYINI(SECT6)
          sect6=0d0
         XCROS6 =SECT6 *XNORM *(2D0/3D0)**4
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
         WTMOD4 = WTVES*WTKIN*DISTRY*XCROS4*COMFAC/DSCRUD
         WTMOD5 = WTVES*WTKIN*DISTRY*XCROS5*COMFAC/DSCRUD
         WTMOD6 = WTVES*WTKIN*DISTRY*XCROS6*COMFAC/DSCRUD
      ENDIF


         IF(KEYBRE.EQ.1) THEN
             WTMOD=WTMOD5
             WTINF=WTMOD4
             WTYFS=WTMOD3
             WTEXA=WTMOD5
         ENDIF


C***      ENDIF
c weight monitoring 
      CALL WMONIT(0,1,WTMOD1,1D0,0D0)
      CALL WMONIT(0,2,WTMOD2,1D0,0D0)
      CALL WMONIT(0,3,WTMOD3,1D0,0D0)
      CALL WMONIT(0,4,WTMOD4,1D0,0D0)
      CALL WMONIT(0,5,WTMOD5,1D0,0D0)
      CALL WMONIT(0,6,WTMOD6,1D0,0D0)
      if (wtmod.gt.xaus2) then
       xaus2=wtmod
       write(*,*) 'maximum wt: xcros1 =',wtmod,' $$$$$$$$$$$$$$$'
       write(*,*)  CTHET,'>>>',1d0/xccos1,1d0/xccos2,xj4,xj5 !,rrr
       write(*,*) xjac1,xjac2,xjac3
       write(*,*) 'keybre= ',keybre,'  wtkin= ',wtkin
       write(*,*) ' '
       write(*,*) ' p2=',p2
       write(*,*) ' q2=',q2
       write(*,*) ' '
       write(*,*) 'ph1=',phot1
       write(*,*) 'ph2=',phot2
       write(*,*) 'ph3=',phot3
       do k=1,4
        phsum(k)=p2(k)+q2(k)+phot1(k)+phot2(k)+phot3(k)
       enddo
       write(*,*) '----------------------------------'
       write(*,*) 'sum=',phsum
       write(*,*) '=================================='
      endif


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
      CALL WMONIT(1,4,AWT4,DWT4,0D0)
      XSEC4=AWT4*XCRU*GNANOB
      DXSEC4=XSEC4*DWT4
      CALL WMONIT(1,5,AWT5,DWT5,0D0)
      XSEC5=AWT5*XCRU*GNANOB
      DXSEC5=XSEC5*DWT5
      CALL WMONIT(1,6,AWT6,DWT6,0D0)
      XSEC6=AWT6*XCRU*GNANOB
      DXSEC6=XSEC6*DWT6

      CALL WMONIT(1,50,AWT50,DWT50,DUM3)

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC3       : WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '   X.sect. in [nb] units         '
      WRITE(NOUT,BXTXT) '   for total generated sample    '
      WRITE(NOUT,BXL1I) NEVGEN,     'generated events   ','NEVGEN','A1'
      WRITE(NOUT,BXL2F) AWT50 ,DWT50   ,'phase space    ','AWT50 ','A2'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   Z0DEC3       : INITIAL        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) 'according to formula in infrared limit YFS3'
      WRITE(NOUT,BXL2F) XSEC1,DXSEC1,'xsection 1   [nb]  ','XSEC4','A2'   
      WRITE(NOUT,BXTXT) 'according to spin amplitude infrared limit'
      WRITE(NOUT,BXL2F) XSEC2,DXSEC2,'xsection     [nb]  ','XSEC2','A2'   
      WRITE(NOUT,BXTXT) 'according to YFS 2 phot * infrared limit ph3'
      WRITE(NOUT,BXL2F) XSEC3,DXSEC3,'xsection     [nb]  ','XSEC3','A2'   
      WRITE(NOUT,BXTXT) 
     $             'according to spin ampl 2 phot * infrared limit ph3'
      WRITE(NOUT,BXL2F) XSEC4,DXSEC4,'xsection     [nb]  ','XSEC4','A2'   
      WRITE(NOUT,BXTXT) 'according to formula from spin amplitudes'
      WRITE(NOUT,BXL2F) XSEC5,DXSEC5,'xsection     [nb]  ','XSEC5','A2'   
      WRITE(NOUT,BXTXT) 'according to formula from spin amplitudes'
      WRITE(NOUT,BXTXT) ' different gauge: commented out !        '
      WRITE(NOUT,BXL2F) XSEC6,DXSEC6,'xsection     [nb]  ','XSEC6','A2'   
      WRITE(NOUT,BXCLO)


      ENDIF
C     =====      
      END

      SUBROUTINE INTEFD
C     *******************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FOTON/ ARBIT(4)
      DIMENSION AA(4)
      COMMON /MOMCMS3/ 
     #       P1(4),Q1(4),P2(4),Q2(4),PHOT1(4),PHOT2(4),PHOT3(4)


*===========================================INITIALIZATION FOTON
*ARBITRARY VECTOR FOR PHOTON POLARIZATION IS FIXED
*ARBITRARY VECTOR FOR PHOTON POLARIZATION IS FIXED
      CALL MOMENTA(P1,AA,X)
      ARBIT(4)=AA(4)
      ARBIT(1)=AA(1)
      ARBIT(3)=AA(3)*COS(75.)-AA(2)*SIN(75.)
      ARBIT(2)=AA(2)*COS(75.)+AA(3)*SIN(75.)

      END


      SUBROUTINE KINE5IIZ0(
     # AMTAx,AMP3,AMP2,AMP1,AMNUTA,amnut2,PIM3,PIM2,PIM1,PN,pn2,WT)
C     *******************************************************************
C generator of 5 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA,amnut2  - masses of particles
C         PIM1,PIM2,PIM3,PN,PN2  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C         pn2-flat new additional particle. Not yet integrated.
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),
     $          pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      common /jakobiany/  CTHET,xccos1,xccos2,xj4,xj5,rrr
      common /jakobi1/ xjac1,xjac2,xjac3
      DIMENSION RRR(19)

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
!      PHSPAC=1.D0/2**17/PI**8
      PHSPAC=1.D0/2**23/PI**11
! for collinear singularity pt1 first photon pt2 second photon
      pt1=2.5
      pt2=2.5
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAx
C
      CALL VARRAN(RRR,19)
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in (PIM1+PIM2)**2 
C GENERATING MASS of l+l- pair
        AMS1=(AMP2+AMP3)**2       
        AMS2=(AMTAx-AMP1-AMNUTA-amnut2)**2 
        ALP1=ATAN((AMS1-AMaz**2)/AMaz/GAMmz)
        ALP2=ATAN((AMS2-AMaz**2)/AMaz/GAMmz)
        prob1=0.5
        prob2=0.5
        IF(RRR(9).LE.prob1) THEN
         am2sq=ams1+(ams2-ams1)*rrr(1)
        ELSE
         ALP=ALP1+RRr(1)*(ALP2-ALP1)
         AM2sq=AMaz**2+AMaz*GAMmz*TAN(ALP) 
        ENDIF
        XJAC1=ams2-ams1
        XJAC2=((AM2sq-AMaz**2)**2+(AMaz*GAMmz)**2)
     $       /(AMaz*GAMmz)*(ALP2-ALP1)
        XJAC=1d0/(prob1/XJAC1+prob2/XJAC2)
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900

C MASS OF ll gam
          AMS1=(am2+amp1)**2
          AMS2=(AMTAx-AMnuta-amnut2)**2
          prob1=0.5
          prob2=0.5
          XK0=0.0051D0 
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0)
          XL0=LOG(2*XK0)
          if (rrr(8).lt.prob1) then
           RR1=RRR(2)
           XK=EXP(XL1*RR1+XL0)
           AM3SQ=(1-XK)*AMS2
           AM3 =SQRT(AM3SQ)
          else
           AM3SQ=AMS1+   RRR(2)*(AMS2-AMS1)
           AM3 =SQRT(AM3SQ)
          endif
          xjac2=(AMS2-AMS1)
          xk=1-am3sq/ams2
          xjac1=AMS2*XL1*XK
          XJAC=1d0/(prob1/XJAC1+prob2/XJAC2)
          PHSPAC=PHSPAC*XJAC
          if (xk.lt.xk0) phspac=0
        IF(PHSPAC.EQ.0D0) GOTO 900
C GENERATING MASS of l+l- gamam gamma 
c...  xk0 defines min. photon pn2 energy in lab 
        ams1=(AM3+AMNUTA)**2
        AMS2=AMTAx**2
        prob1=0.4
        prob2=0.4
        prob3=0.2
          XK0=0.0251D0 
!          XK0=0.0051D0 
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0)
          XL0=LOG(2*XK0)
          if (rrr(7).lt.prob1) then
           RR1=RRR(11)
           XK=EXP(XL1*RR1+XL0)
           AM4SQ=(1-XK)*AMS2
           AM4 =SQRT(AM4SQ)
          elseif (rrr(7).lt.prob1+prob2) then
           RR1=RRR(11)
           XK=EXP(XL1*RR1+XL0)
           AM4SQa=(1-XK)*AMS2
           am4sq=ams2-am4sqa+ams1
           AM4 =SQRT(AM4SQ)
          else
           AM4SQ=AMS1+   RRR(11)*(AMS2-AMS1)
           AM4 =SQRT(AM4SQ)
          endif
          xjac3=(AMS2-AMS1)
          xk=1-am4sq/ams2
          am4sqa=ams2-am4sq+ams1
          xka=1-am4sqa/ams2
          xjac1=AMS2*XL1*XK
          xjac2=AMS2*XL1*xka
          XJAC=1d0/(prob1/XJAC1+prob2/XJAC2+prob3/xjac3)
          PHSPAC=PHSPAC*XJAC
          if (xk.lt.  xk0) phspac=0
          if (xka.lt.xk0) phspac=0
        IF(PHSPAC.EQ.0D0) GOTO 900
       
      amtau=sqrt(am4sq)

* AM2 RESTFRAME, DEFINE PIM2 AND PIM3
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
* NOW boost TO THE am3 REST FRAME
      paa4=1./(2*am3)*(am3**2-amp1**2+am2**2)                    
      paa3= sqrt(abs(paa4**2-am2**2))                                            
      exe=(paa4+paa3)/am2
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)
      eee=0
      do k=1,3
       pim1(k)=-pim2(k)-pim3(k)
       eee=eee+pim1(k)**2
      enddo
      pim1(4)=sqrt(eee)
        PHSPAC=PHSPAC*(4*PI)*(2*pim1(4)/AM3)
* ALL ROTATED IN THE am3 REST FRAME
! here pt1 for mass and rrr(18) for branch)
      THET =ACOS(-1.D0+2*RRR(16))
        rr3=rrr(16)
        prev=0.3
        IF(RRR(18).lt.PREV) then
         EPS=min((pt1/am3)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         ETA  =EXP(XL1*RR3+XL0)
         CTHET=-(1+EPS-ETA)
         THET =ACOS(CTHET)
        elseIF(RRR(18).lt.2*PREV) then
         EPS=min((pt1/am3)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         ETA  =EXP(XL1*RR3+XL0)
         CTHET=(1+EPS-ETA)
         THET =ACOS(CTHET)
        else
         CTHET=-1+2*rr3
         THET =ACOS(CTHET)
        endif
         EPS=min((pt1/am3)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         eta1=1+eps+cthet
         eta2=1+eps-cthet
!       write(*,*) thet2, 1/(PREV/(XL1/2*ETA)+(1d0-prev)/1D0)


        PHSPAC=PHSPAC/( PREV/(XL1/2*ETA1)+PREV/(XL1/2*ETA2)
     $                +(1d0-2*prev)/1D0)


      PHI = 2*PI*RRR(17)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)

* AM4/amtau RESTFRAME, DEFINE PN
        ENQ1=(Amtau**2-AM3**2+AMNUTA**2)/(2*AMtau)
        ENQ2=(AMtau**2+AM3**2-AMNUTA**2)/(2*AMtau)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMtau)
* PI MINUS MOMENTUM IN RHO REST FRAME
* NOW boost TO THE am4 REST FRAME
      paa4=1./(2*amtau)*(amtau**2-amnuta**2+am3**2)                    
      paa3= sqrt(abs(paa4**2-am3**2))                                            
      exe=(paa4+paa3)/am3
      CALL BOSTd3(EXE,PIm1,PIm1)
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)

        pn(1)=0d0
        pn(2)=0d0
        pn(3)=-pppi
        PN(4)=ENQ1
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
! here pt2 for mass and rrr(19) for branch

        rr3=rrr(5)
        prev=0.3
        IF(RRR(19).lt.PREV) then
         EPS=min((pt2/amtau)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         ETA  =EXP(XL1*RR3+XL0)
         CTHET=-(1+EPS-ETA)
         THET =ACOS(CTHET)
        elseIF(RRR(19).lt.2*PREV) then
         EPS=min((pt2/amtau)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         ETA  =EXP(XL1*RR3+XL0)
         CTHET=(1+EPS-ETA)
         THET =ACOS(CTHET)
        else
         CTHET=-1+2*rr3
         THET =ACOS(CTHET)
        endif
         EPS=min((pt2/amtau)**2,0.8d0)
         XL1=LOG((2+EPS)/EPS)
         XL0=LOG(EPS)
         eta1=1+eps+cthet
         eta2=1+eps-cthet
!       write(*,*) thet2, 1/(PREV/(XL1/2*ETA)+(1d0-prev)/1D0)


        PHSPAC=PHSPAC/( PREV/(XL1/2*ETA1)+PREV/(XL1/2*ETA2)
     $                +(1d0-2*prev)/1D0)


      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)

* now to the tau rest frame, define paa and neutrino momenta            
* paa  momentum                                                         
      paa(1)=0                                                          
      paa(2)=0                                                          
      paa(4)=1./(2*amtax)*(amtax**2-amnut2**2+amtau**2)                    
      paa(3)= sqrt(abs(paa(4)**2-amtau**2))                                            
      phspac=phspac*(4*pi)*(2*paa(3)/amtax)                             
* tau-neutrino momentum                                                 
      pn2(1)=0                                                           
      pn2(2)=0                                                           
      pn2(4)=1./(2*amtax)*(amtax**2+amnut2**2-amtau**2)                    
      pn2(3)=-paa(3)  

      exe=(paa(4)+paa(3))/amtau                                           
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)

        prob1=.2
        prob2=.2
        prob3=.2
        prob4=0.2
        prob5=0.2
        EPS=(AMp3/AMtax)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
      IF    (RRR(15).lt.PROB1) then 
       THET =ACOS(-1.D0+2*RRR(12))
       CTHET=COS(THET)
      elseIF(RRR(15).lt.(PROB1+PROB2)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=-(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
!         CTHET=-cthet
        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)

        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3+PROB4)) then
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=AM2SQX-2D0+sqrt(1d0-eps)
         xx=eps
         beta=sqrt(1d0-eps)
         xlog=-log((1+beta)**2/xx)
         xlog1=-log(16D0/xx)
          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
         CTHET=-cthet
        THET =ACOS(CTHET)
      else
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=-AM2SQX+2D0-sqrt(1d0-eps)
         xx=eps
         beta=sqrt(1d0-eps)
         xlog=-log((1+beta)**2/xx)
         xlog1=-log(16D0/xx)
          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
        THET =ACOS(CTHET)
      endif
      if (cthet**2.gt.1d0) then
       cthet=cthet/cthet**2
       write(*,*) 'cthet error -- arbi action'
       write(*,*) cthet,rrr(12),rrr(15)
       write(*,*) ams1,ams2,am2sq
        THET =ACOS(CTHET)
      endif
      eta1=1+eps+cthet
      eta2=1+eps-cthet
      xx=eps
      beta=sqrt(1d0-eps)
      xx=eps
      xlog=-log((1+beta)**2/xx)
      xlog1=-log(16D0/xx)
      ct=-cthet

      xa=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))
      ct=cthet
      xb=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))

      xccos1=1d0/(XL1/2*ETA1)
      xccos2=1d0/(XL1/2*ETA2)

        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
        n=1
        AM2SQX= CTHET+2D0-sqrt(1d0-eps)
        xj4=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
        n=1
        AM2SQX=-CTHET+2D0-sqrt(1d0-eps)
        xj5=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
      xj4=1.d0/xa
      xj5=1.d0/xb
      FF=PROB1/1d0+PROB2*xccos1+PROB3*xccos2+PROB4/XJ4+PROB5/XJ5

      PHSPAC=PHSPAC/FF                                                   
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
!      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)


C FINAL WEIGHT
      WT = PHSPAC
C THE STATISTICAL FACTOR FOR IDENTICAL gammas 
C is replaced with ordering in p_T
      pt3=pn2(1)**2+pn2(2)**2
      pt2=pn(1)**2+pn(2)**2
      pt1=pim1(1)**2+pim1(2)**2
      if (pt1.lt.pt2) wt=0
      if (pt2.lt.pt3) wt=0
      RETURN
 900  WT=0D0

      END

      SUBROUTINE KINE5IIZ0a(
     # AMTAx,AMP3,AMP2,AMP1,AMNUTA,amnut2,PIM3,PIM2,PIM1,PN,pn2,WT)
C     *******************************************************************
C generator of 5 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA,amnut2  - masses of particles
C         PIM1,PIM2,PIM3,PN,PN2  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C         pn2-flat new additional particle. Not yet integrated.
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),
     $          pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(17)

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
!      PHSPAC=1.D0/2**17/PI**8
      PHSPAC=1.D0/2**23/PI**11

C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAx
C
      CALL VARRAN(RRR,17)
! flat phase space
      ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
      ams2=(amtax-amnut2)**2
!     amx4=ams1+rrr(11)*(ams2-ams1)
!      phspac=phspac*(AMS2-AMS1)
        if (amtax.gt.ams1) goto 5
        amro=100
        gamro=100
        prob1=.5
        prob2=.5
        PROB3=.0
        PROB4=.0
        rr2=rrr(11)
        ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
        ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
        IF (RRR(14).LT.PROB1) THEN
         AM2SQ=AMS1+   RR2*(AMS2-AMS1)
         AM2 =SQRT(AM2SQ)
        elseIF (RRR(14).LT.(PROB1+PROB2)) THEN  
         B=LOG(AMS1)
         A=LOG(AMS2)
         AM2SQ=AMS2*EXP((B-A)*RR2)
         am2sq=ams2+ams1-am2sq
         AM2 =SQRT(AM2SQ)     
        ELSEIF (RRR(14).LT.(PROB1+PROB2+PROB3)) THEN
         ALP=ALP1+RR2*(ALP2-ALP1)
         AM2SQ=AMRO**2+AMRO*GAMRO*TAN(ALP)
         AM2 =SQRT(AM2SQ)
        ELSE
         n=1
          if(n.eq.1) then
         AM2SQ=AMS1/(1D0-RR2*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQ=AMS1/sqrt(1D0-RR2*(1-(ams1/ams2)**n))
          else
         AM2SQ=AMS1*(1D0-RR2*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
         AM2 =SQRT(AM2SQ)
         if (am2sq.gt.ams2) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.gt.ams2) stop
         if (am2sq.lt.ams1) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.lt.ams1) stop

        ENDIF

        XJ1=(AMS2-AMS1)
         B=LOG(AMS1)
         A=LOG(AMS2)
         am2sqx=ams2+ams1-am2sq
        xj2=AM2SQx*(A-B)
        xj3=((AM2SQ-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
        xj3=xj3*(ALP2-ALP1)
        n=1
        xj4=am2SQ**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)
!        sum=Sum+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
!        sum2=Sum2+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)**2
!        enddo
!        sum=sum/nn
!        sum2=sum2/nn
!        err=sqrt((sum2-sum**2)/nn)
!        write(*,*) sum,'+-',err
!        write(*,*) '28761.2547837270613 +- 0'
!        stop
        PHSPAC=PHSPAC/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
 5     continue

C>>> 1 soft photon
C>>>          RRR(11)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in pn2(4)
C GENERATING MASS of l+l- gamam gamma 
c...  xk0 defines min. photon pn2 energy in lab 
        ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
        XK0=0.0051D0 
        AMS2=AMTAx**2

          RR1=RRR(11)
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0)
          XL0=LOG(2*XK0)
          XK=EXP(XL1*RR1+XL0)
          AM3SQ=(1-XK)*AMS2
          AM3 =SQRT(AM3SQ)
          PHSPAC=PHSPAC*AMS2*XL1*XK
        IF(PHSPAC.EQ.0D0) GOTO 900
       
      amtau=sqrt(am3sq)


C>>> do testow na konfiguracje podczerwone
C>>> 2 soft photons
C>>>      RRR(1)=0.999998D0
C>>>      RRR(2)=0.999998D0
C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in (PIM1+PIM2)**2 
C GENERATING MASS of l+l- pair
        AMS1=(AMP2+AMP3)**2       
        AMS2=(AMTAU-AMP1-AMNUTA)**2 
        ALP1=ATAN((AMS1-AMaz**2)/AMaz/GAMmz)
        ALP2=ATAN((AMS2-AMaz**2)/AMaz/GAMmz)
        prob1=0.5
        prob2=0.5
        IF(RRR(9).LE.prob1) THEN
         am2sq=ams1+(ams2-ams1)*rrr(1)
        ELSE
         ALP=ALP1+RRr(1)*(ALP2-ALP1)
         AM2sq=AMaz**2+AMaz*GAMmz*TAN(ALP) 
        ENDIF
        XJAC1=ams2-ams1
        XJAC2=((AM2sq-AMaz**2)**2+(AMaz*GAMmz)**2)
     $       /(AMaz*GAMmz)*(ALP2-ALP1)
        XJAC=1d0/(prob1/XJAC1+prob2/XJAC2)
        AM2=SQRT(AM2SQ)
        PHSPAC=PHSPAC*XJAC
        IF(PHSPAC.EQ.0D0) GOTO 900

C MASS OF ll gam
        AMS1=(am2+amp1)**2
        AMS2=(AMTAU-AMnuta)**2
        AM3SQ=AMS1+   RRR(2)*(AMS2-AMS1)
        AM3 =SQRT(AM3SQ)
        PHSPAC=PHSPAC*(AMS2-AMS1)

* AM2 RESTFRAME, DEFINE PIM2 AND PIM3
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
* NOW boost TO THE am3 REST FRAME
      paa4=1./(2*am3)*(am3**2-amp1**2+am2**2)                    
      paa3= sqrt(abs(paa4**2-am2**2))                                            
      exe=(paa4+paa3)/am2
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)
      eee=0
      do k=1,3
       pim1(k)=-pim2(k)-pim3(k)
       eee=eee+pim1(k)**2
      enddo
      pim1(4)=sqrt(eee)
        PHSPAC=PHSPAC*(4*PI)*(2*pim1(4)/AM3)
* ALL ROTATED IN THE am3 REST FRAME
      THET =ACOS(-1.D0+2*RRR(16))
      PHI = 2*PI*RRR(17)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)

* AM4/amtau RESTFRAME, DEFINE PN
        ENQ1=(Amtau**2-AM3**2+AMNUTA**2)/(2*AMtau)
        ENQ2=(AMtau**2+AM3**2-AMNUTA**2)/(2*AMtau)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMtau)
* PI MINUS MOMENTUM IN RHO REST FRAME
* NOW boost TO THE am4 REST FRAME
      paa4=1./(2*amtau)*(amtau**2-amnuta**2+am3**2)                    
      paa3= sqrt(abs(paa4**2-am3**2))                                            
      exe=(paa4+paa3)/am3
      CALL BOSTd3(EXE,PIm1,PIm1)
      CALL BOSTd3(EXE,PIm2,PIm2)
      CALL BOSTd3(EXE,PIm3,PIm3)

        pn(1)=0d0
        pn(2)=0d0
        pn(3)=-pppi
        PN(4)=ENQ1
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(5))
      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)

* now to the tau rest frame, define paa and neutrino momenta            
* paa  momentum                                                         
      paa(1)=0                                                          
      paa(2)=0                                                          
      paa(4)=1./(2*amtax)*(amtax**2-amnut2**2+amtau**2)                    
      paa(3)= sqrt(abs(paa(4)**2-amtau**2))                                            
      phspac=phspac*(4*pi)*(2*paa(3)/amtax)                             
* tau-neutrino momentum                                                 
      pn2(1)=0                                                           
      pn2(2)=0                                                           
      pn2(4)=1./(2*amtax)*(amtax**2+amnut2**2-amtau**2)                    
      pn2(3)=-paa(3)  

      exe=(paa(4)+paa(3))/amtau                                           
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)

        prob1=1.
        prob2=.0
        prob3=.0
        prob4=0
        prob5=0
        EPS=(AMp3/AMtax)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
      IF    (RRR(15).lt.PROB1) then 
       THET =ACOS(-1.D0+2*RRR(12))
       CTHET=COS(THET)
      elseIF(RRR(15).lt.(PROB1+PROB2)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=-(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
!         CTHET=-cthet
        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)

        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3+PROB4)) then
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=AM2SQX-2D0+sqrt(1d0-eps)
        THET =ACOS(CTHET)
      else
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=-AM2SQX+2D0-sqrt(1d0-eps)
        THET =ACOS(CTHET)
      endif
      if (cthet**2.gt.1d0) then
       cthet=cthet/cthet**2
       write(*,*) 'cthet error -- arbi action'
       write(*,*) cthet,rrr(12),rrr(15)
       write(*,*) ams1,ams2,am2sq
        THET =ACOS(CTHET)
      endif
      eta1=1+eps+cthet
      eta2=1+eps-cthet
      xx=eps
      beta=sqrt(1d0-eps)
      xx=eps
      xlog=-log((1+beta)**2/xx)
      xlog1=-log(16D0/xx)
      ct=-cthet

      xccos1=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))
      ct=cthet
      xccos2=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))

      xccos1=1d0/(XL1/2*ETA1)
      xccos2=1d0/(XL1/2*ETA2)

        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
        n=1
        AM2SQX= CTHET+2D0-sqrt(1d0-eps)
        xj4=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
        n=1
        AM2SQX=-CTHET+2D0-sqrt(1d0-eps)
        xj5=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
       
      FF=PROB1/1d0+PROB2*xccos1+PROB3*xccos2+PROB4/XJ4+PROB5/XJ5

      PHSPAC=PHSPAC/FF                                                   
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
!      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)


C FINAL WEIGHT
      WT = PHSPAC
C THE STATISTICAL FACTOR FOR IDENTICAL gammas 
C is replaced with ordering in p_T
      pt3=pn2(1)**2+pn2(2)**2
      pt2=pn(1)**2+pn(2)**2
      pt1=pim1(1)**2+pim1(2)**2
      if (pt1.lt.pt2) wt=0
      if (pt2.lt.pt3) wt=0
      RETURN
 900  WT=0D0

      END
      SUBROUTINE KINE5IIZ0okx(
     # AMTAx,AMP3,AMP2,AMP1,AMNUTA,amnut2,PIM3,PIM2,PIM1,PN,pn2,WT)
C     *******************************************************************
C generator of 5 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA,amnut2  - masses of particles
C         PIM1,PIM2,PIM3,PN,PN2  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C         pn2-flat new additional particle. Not yet integrated.
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),
     $          pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(15)

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
!      PHSPAC=1.D0/2**17/PI**8
      PHSPAC=1.D0/2**23/PI**11

C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAx
C
      CALL VARRAN(RRR,15)
! flat phase space
      ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
      ams2=(amtax-amnut2)**2
!     amx4=ams1+rrr(11)*(ams2-ams1)
!      phspac=phspac*(AMS2-AMS1)
        if (amtax.gt.ams1) goto 5
        amro=100
        gamro=100
        prob1=.5
        prob2=.5
        PROB3=.0
        PROB4=.0
        rr2=rrr(11)
        ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
        ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
        IF (RRR(14).LT.PROB1) THEN
         AM2SQ=AMS1+   RR2*(AMS2-AMS1)
         AM2 =SQRT(AM2SQ)
        elseIF (RRR(14).LT.(PROB1+PROB2)) THEN  
         B=LOG(AMS1)
         A=LOG(AMS2)
         AM2SQ=AMS2*EXP((B-A)*RR2)
         am2sq=ams2+ams1-am2sq
         AM2 =SQRT(AM2SQ)     
        ELSEIF (RRR(14).LT.(PROB1+PROB2+PROB3)) THEN
         ALP=ALP1+RR2*(ALP2-ALP1)
         AM2SQ=AMRO**2+AMRO*GAMRO*TAN(ALP)
         AM2 =SQRT(AM2SQ)
        ELSE
         n=1
          if(n.eq.1) then
         AM2SQ=AMS1/(1D0-RR2*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQ=AMS1/sqrt(1D0-RR2*(1-(ams1/ams2)**n))
          else
         AM2SQ=AMS1*(1D0-RR2*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
         AM2 =SQRT(AM2SQ)
         if (am2sq.gt.ams2) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.gt.ams2) stop
         if (am2sq.lt.ams1) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.lt.ams1) stop

        ENDIF

        XJ1=(AMS2-AMS1)
         B=LOG(AMS1)
         A=LOG(AMS2)
         am2sqx=ams2+ams1-am2sq
        xj2=AM2SQx*(A-B)
        xj3=((AM2SQ-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
        xj3=xj3*(ALP2-ALP1)
        n=1
        xj4=am2SQ**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)
!        sum=Sum+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
!        sum2=Sum2+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)**2
!        enddo
!        sum=sum/nn
!        sum2=sum2/nn
!        err=sqrt((sum2-sum**2)/nn)
!        write(*,*) sum,'+-',err
!        write(*,*) '28761.2547837270613 +- 0'
!        stop
        PHSPAC=PHSPAC/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
 5     continue

C>>> 1 soft photon
C>>>          RRR(11)=0.99999D0
C PHASE SPACE WITH PRESAMPLING ON RESONANCE and POLE in pn2(4)
C GENERATING MASS of l+l- gamam gamma 
c...  xk0 defines min. photon pn2 energy in lab 
        ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
        XK0=0.0051D0 
        AMS2=AMTAx**2

          RR1=RRR(11)
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0)
          XL0=LOG(2*XK0)
          XK=EXP(XL1*RR1+XL0)
          AM3SQ=(1-XK)*AMS2
          AM3 =SQRT(AM3SQ)
          PHSPAC=PHSPAC*AMS2*XL1*XK
        IF(PHSPAC.EQ.0D0) GOTO 900
       
      amtau=sqrt(am3sq)


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
        XK0=0.0051D0 
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




* now to the tau rest frame, define paa and neutrino momenta            
* paa  momentum                                                         
      paa(1)=0                                                          
      paa(2)=0                                                          
      paa(4)=1./(2*amtax)*(amtax**2-amnut2**2+amtau**2)                    
      paa(3)= sqrt(abs(paa(4)**2-amtau**2))                                            
      phspac=phspac*(4*pi)*(2*paa(3)/amtax)                             
* tau-neutrino momentum                                                 
      pn2(1)=0                                                           
      pn2(2)=0                                                           
      pn2(4)=1./(2*amtax)*(amtax**2+amnut2**2-amtau**2)                    
      pn2(3)=-paa(3)  

      exe=(paa(4)+paa(3))/amtau                                           
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)

        prob1=.2
        prob2=.4
        prob3=.4
        prob4=0
        prob5=0
        EPS=(AMp3/AMtax)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
      IF    (RRR(15).lt.PROB1) then 
       THET =ACOS(-1.D0+2*RRR(12))
       CTHET=COS(THET)
      elseIF(RRR(15).lt.(PROB1+PROB2)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=-(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
!         CTHET=-cthet
        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)

        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3+PROB4)) then
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=AM2SQX-2D0+sqrt(1d0-eps)
        THET =ACOS(CTHET)
      else
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=-AM2SQX+2D0-sqrt(1d0-eps)
        THET =ACOS(CTHET)
      endif
      if (cthet**2.gt.1d0) then
       cthet=cthet/cthet**2
       write(*,*) 'cthet error -- arbi action'
       write(*,*) cthet,rrr(12),rrr(15)
       write(*,*) ams1,ams2,am2sq
        THET =ACOS(CTHET)
      endif
      eta1=1+eps+cthet
      eta2=1+eps-cthet
      xx=eps
      beta=sqrt(1d0-eps)
      xx=eps
      xlog=-log((1+beta)**2/xx)
      xlog1=-log(16D0/xx)
      ct=-cthet

      xccos1=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))
      ct=cthet
      xccos2=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))

      xccos1=1d0/(XL1/2*ETA1)
      xccos2=1d0/(XL1/2*ETA2)

        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
        n=1
        AM2SQX= CTHET+2D0-sqrt(1d0-eps)
        xj4=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
        n=1
        AM2SQX=-CTHET+2D0-sqrt(1d0-eps)
        xj5=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
       
      FF=PROB1/1d0+PROB2*xccos1+PROB3*xccos2+PROB4/XJ4+PROB5/XJ5

      PHSPAC=PHSPAC/FF                                                   
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
!      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)


C FINAL WEIGHT
      WT = PHSPAC
C THE STATISTICAL FACTOR FOR IDENTICAL gammas 
C is replaced with ordering in p_T
      pt3=pn2(1)**2+pn2(2)**2
      pt2=pn(1)**2+pn(2)**2
      pt1=pim1(1)**2+pim1(2)**2
      if (pt1.lt.pt2) wt=0
      if (pt2.lt.pt3) wt=0
      RETURN
 900  WT=0D0

      END

      SUBROUTINE KINE5IIZ0ok(
     # AMTAx,AMP3,AMP2,AMP1,AMNUTA,amnut2,PIM3,PIM2,PIM1,PN,pn2,WT)
C     *******************************************************************
C generator of 5 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA,amnut2  - masses of particles
C         PIM1,PIM2,PIM3,PN,PN2  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C         pn2-flat new additional particle. Not yet integrated.
C    factor 1/2 for two identical particles included
C         WT  - weight
C presampling on singularity and resonance in (PIM2+PIM3) mass
C subroutine is based on subroutine DPHTRE from TAUOLA
C but algorithm of generation is sleightly different
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PAA(4),PBB(4),PT(4),
     $          pn2(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      DIMENSION RRR(15)

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
!      PHSPAC=1.D0/2**17/PI**8
      PHSPAC=1.D0/2**23/PI**11

C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAx
C
      CALL VARRAN(RRR,15)
! flat phase space
      ams1=(AMP3+AMP2+AMP1+AMNUTA)**2
      ams2=(amtax-amnut2)**2
!     amx4=ams1+rrr(11)*(ams2-ams1)
!      phspac=phspac*(AMS2-AMS1)
        amro=100
        gamro=100
        prob1=.5
        prob2=.5
        PROB3=.0
        PROB4=.0
        rr2=rrr(11)
        ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
        ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
        IF (RRR(14).LT.PROB1) THEN
         AM2SQ=AMS1+   RR2*(AMS2-AMS1)
         AM2 =SQRT(AM2SQ)
        elseIF (RRR(14).LT.(PROB1+PROB2)) THEN  
         B=LOG(AMS1)
         A=LOG(AMS2)
         AM2SQ=AMS2*EXP((B-A)*RR2)
         am2sq=ams2+ams1-am2sq
         AM2 =SQRT(AM2SQ)     
        ELSEIF (RRR(14).LT.(PROB1+PROB2+PROB3)) THEN
         ALP=ALP1+RR2*(ALP2-ALP1)
         AM2SQ=AMRO**2+AMRO*GAMRO*TAN(ALP)
         AM2 =SQRT(AM2SQ)
        ELSE
         n=1
          if(n.eq.1) then
         AM2SQ=AMS1/(1D0-RR2*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQ=AMS1/sqrt(1D0-RR2*(1-(ams1/ams2)**n))
          else
         AM2SQ=AMS1*(1D0-RR2*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
         AM2 =SQRT(AM2SQ)
         if (am2sq.gt.ams2) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.gt.ams2) stop
         if (am2sq.lt.ams1) WRITE(*,*) 'am2sq',am2sq,ams1,ams2,rr2
         if (am2sq.lt.ams1) stop

        ENDIF

        XJ1=(AMS2-AMS1)
         B=LOG(AMS1)
         A=LOG(AMS2)
         am2sqx=ams2+ams1-am2sq
        xj2=AM2SQx*(A-B)
        xj3=((AM2SQ-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
        xj3=xj3*(ALP2-ALP1)
        n=1
        xj4=am2SQ**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)
!        sum=Sum+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
!        sum2=Sum2+1d0/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)**2
!        enddo
!        sum=sum/nn
!        sum2=sum2/nn
!        err=sqrt((sum2-sum**2)/nn)
!        write(*,*) sum,'+-',err
!        write(*,*) '28761.2547837270613 +- 0'
!        stop
        PHSPAC=PHSPAC/(PROB1/XJ1+PROB2/XJ2+PROB3/XJ3+PROB4/XJ4)
      amtau=sqrt(am2sq)


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
        XK0=0.0051D0 
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




* now to the tau rest frame, define paa and neutrino momenta            
* paa  momentum                                                         
      paa(1)=0                                                          
      paa(2)=0                                                          
      paa(4)=1./(2*amtax)*(amtax**2-amnut2**2+amtau**2)                    
      paa(3)= sqrt(abs(paa(4)**2-amtau**2))                                            
      phspac=phspac*(4*pi)*(2*paa(3)/amtax)                             
* tau-neutrino momentum                                                 
      pn2(1)=0                                                           
      pn2(2)=0                                                           
      pn2(4)=1./(2*amtax)*(amtax**2+amnut2**2-amtau**2)                    
      pn2(3)=-paa(3)  

      exe=(paa(4)+paa(3))/amtau                                           
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PN,PN)

        prob1=.4
        prob2=.3
        prob3=.3
        prob4=0
        prob5=0
        EPS=(AMp3/AMtax)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)

      IF    (RRR(15).lt.PROB1) then 
       THET =ACOS(-1.D0+2*RRR(12))
       CTHET=COS(THET)
      elseIF(RRR(15).lt.(PROB1+PROB2)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=-(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!         cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)
!         CTHET=-cthet
        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3)) then
        ETA  =EXP(XL1*RRR(12)+XL0)
        CTHET=(1+EPS-ETA)
!         xx=eps
!         beta=sqrt(1d0-eps)
!         xlog=-log((1+beta)**2/xx)
!         xlog1=-log(16D0/xx)
!          u=(log((1D0+beta)/4D0))**2 +xlog*xlog1*rrr(12)
!          cthet=-1D0/beta*(4D0*EXP(-SQRT(u))-1)

        THET =ACOS(CTHET)
      elseIF    (RRR(15).lt.(PROB1+PROB2+PROB3+PROB4)) then
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=AM2SQX-2D0+sqrt(1d0-eps)
        THET =ACOS(CTHET)
      else
        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
          n=1
          if(n.eq.1) then
         AM2SQX=AMS1/(1D0-RRr(12)*(1-(ams1/ams2)**n))
          elseif(n.eq.2) then
         AM2SQX=AMS1/sqrt(1D0-RRr(12)*(1-(ams1/ams2)**n))
          else
         AM2SQX=AMS1*(1D0-RRr(12)*(1-(ams1/ams2)**n))**(-1d0/n)
          endif
        CTHET=-AM2SQX+2D0-sqrt(1d0-eps)
        THET =ACOS(CTHET)
      endif
      if (cthet**2.gt.1d0) then
       cthet=cthet/cthet**2
       write(*,*) 'cthet error -- arbi action'
       write(*,*) cthet,rrr(12),rrr(15)
       write(*,*) ams1,ams2,am2sq
        THET =ACOS(CTHET)
      endif
      eta1=1+eps+cthet
      eta2=1+eps-cthet
      xx=eps
      beta=sqrt(1d0-eps)
      xx=eps
      xlog=-log((1+beta)**2/xx)
      xlog1=-log(16D0/xx)
      ct=-cthet

      xccos1=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))
      ct=cthet
      xccos2=beta/(xlog*xlog1
     $     /log(4d0/(xx/(1d0+beta)+beta*(1D0-ct)))
     $     /(4d0/(xx/(1d0+beta)+beta*(1D0-ct))))!!! +1d0/(1+beta*costhe))

      xccos1=1d0/(XL1/2*ETA1)
      xccos2=1d0/(XL1/2*ETA2)

        ams1=1-sqrt(1d0-eps)
        ams2=3-sqrt(1d0-eps)
        n=1
        AM2SQX= CTHET+2D0-sqrt(1d0-eps)
        xj4=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
        n=1
        AM2SQX=-CTHET+2D0-sqrt(1d0-eps)
        xj5=am2SQX**(n+1)*n*(1D0/ams1**n-1D0/ams2**n)/2
       
      FF=PROB1/1d0+PROB2*xccos1+PROB3*xccos2+PROB4/XJ4+PROB5/XJ5

      PHSPAC=PHSPAC/FF                                                   
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
!      THET =ACOS(-1.D0+2*RRR(12))
      PHI = 2*PI*RRR(13)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
      CALL ROTPOD(THET,PHI,PN2)

C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
!##########        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END








