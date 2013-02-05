      PROGRAM MAIN
!     ***********************************
! To execute:  make BeMaFigs-dvi
!              make BeMaFigs-ps
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      Tesnam    = 'BeMaFigs'
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms 
!--------------------------------------------------------
      lendan = 0
      lendan = lendan+1
      Dname(lendan)  = '../prod2/prod2.data.2037M'  ! April 96
      Hname(lendan)  = '../prod2/prod2.hst.2037M'   ! April 96
      Dname(lendan)  = '../prod2/prod2.data'        ! Current
      Hname(lendan)  = '../prod2/bhl.hst'           ! Current
      lendan = lendan+1
      Dname(lendan)  = '../obis2/obis2.data.2118M'  ! April 96
      Hname(lendan)  = '../obis2/obis2.hst.2118M'   ! April 96
      Dname(lendan)  = '../obis2/obis2.data'        ! Current
      Hname(lendan)  = '../obis2/bhl.hst'           ! Current
      lendan = lendan+1
      Dname(lendan)  = '../llog2/llog2.data.2114M'  ! April 96
      Hname(lendan)  = '../llog2/llog2.hst.2114M'   ! April 96
      Dname(lendan)  = '../llog2/llog2.data'        ! Current
      Hname(lendan)  = '../llog2/bhl.hst'           ! Current
!==========================================================
! Table for workshop
      CALL wshop
!==========================================================
!--------------------------------------------------------
!  dumping histogram for control
!--------------------------------------------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL grfile(nouth,' ','N')
      CALL grout( 0,ICY,' ')
      CALL grend(tname)
      CLOSE(nout)
      END


      SUBROUTINE wshop
!     ************************
! tables for wshop proceedings
!     *********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI     =  3.1415926535897932D0)
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )
      SAVE
      COMMON / INOUT  / NINP,NOUT
!-----------------------------------------------------------------------
! Communicates with MAIN
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!-----------------------------------------------------------------------
! Communicates with READAT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
!-----------------------------------------------------------------------
      DIMENSION cut(100)
      DIMENSION bin1(100),err1(100)
      DIMENSION bin2(100),err2(100)
      DIMENSION bin3(100),err3(100)
      CHARACTER*80 title
      DIMENSION BorW(10),BorN(10)
!-------------------------
! Parameters for tables
      DIMENSION    idl(5)
      CHARACTER*16 capt(6)
      CHARACTER*8  fmt(3), fmtx,fmty
!-----------------------------------
!-----------------------------------
      CHARACTER*64 cpsical1(50)
      DATA cpsical1 /
     $'This is from the output file 4.x-cpc/figs/BeMaFigs.tex.',
     $'The difference $d_{3}=',
     $'({\\rm OLDBIS-LUMLOG}_{h.o.}-{\\rm BHLUMI}.4x)/{\\rm Born}$',
     $'for the SICAL event selection',
     $'as a function of the energy cut $1-U_{\min}$.',
     $'Dash-box marks 1.5 permile precision limits.',
     $'\\label{fig:BeMaFigs-OL}',
     $'% end-of-caption'/
!-----------------------------------
      CHARACTER*64 cpsical2(50)
      DATA cpsical2 /
     $'This is from the output file 4.x-cpc/figs/BeMaFigs.tex.',
     $'The difference $d_{3}=',
     $'({\\rm BHLPV}-{\\rm BHLUMI}.4x)/{\\rm Born}$',
     $'for the SICAL event selection',
     $'as a function of the energy cut $1-U_{\min}$.',
     $'BHLP is an emulation of SBASPV using OLDBIS and LUMLOG.',
     $'Dash-box marks 1.0 permile precision limits.',
     $'\\label{fig:BeMaFigs-PV}',
     $'% end-of-caption'/
!-----------------------------------
!-------------------------
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (dot     ='\\circle*{10}')
!-----------------------------------
!-----------------------------------
      CHARACTER*64 labex1(40)
      DATA labex1 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put( 700, 1170){\\makebox(0,0)[t]{\\LARGE ',
     $'  ${\\rm{(OLDBIS+LUMLOG_{h.o.})} ',
     $'   -\\rm{BHLUMI 4.04} \\over\\rm{Born} }$ }}',
     $'%%%%%%%%',
     $' \\put(100,880){\\begin{picture}( 300,400)',
!
     $' \\put(  0,120){\\circle*{10}}',
     $' \\put( 50,120){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,120){\\makebox(0,0)[l]{\\large WW}}',
!
     $' \\put(  0, 60){\\circle*{20}}',
     $' \\put( 50, 60){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280, 60){\\makebox(0,0)[l]{\\large NW}}',
!
     $' \\put(  0,  0){\\circle{30}}',
     $' \\put( 50,  0){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,  0){\\makebox(0,0)[l]{\\large NN}}',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 300, 300){\\dashbox{15}(600,600){ }} ',
     $' \\put( 900,  40){\\makebox(0,0)[b]{\\LARGE $ 1-U_{\\min} $}}',
!>     $' \\put( 673,0){\\line(0,1){1200}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!-----------------------------------
!-----------------------------------
      CHARACTER*64 labex2(40)
      DATA labex2 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put( 700, 1170){\\makebox(0,0)[t]{\\LARGE ',
     $'  ${\\rm{(BHLPV)} ',
     $'   -\\rm{BHLUMI 4.04} \\over\\rm{Born} }$ }}',
     $'%%%%%%%%',
     $' \\put(150,150){\\begin{picture}( 300,400)',
!
     $' \\put(  0,120){\\circle*{10}}',
     $' \\put( 50,120){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,120){\\makebox(0,0)[l]{\\large WW}}',
!
     $' \\put(  0, 60){\\circle*{20}}',
     $' \\put( 50, 60){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280, 60){\\makebox(0,0)[l]{\\large NW}}',
!
     $' \\put(  0,  0){\\circle{30}}',
     $' \\put( 50,  0){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,  0){\\makebox(0,0)[l]{\\large NN}}',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 300, 400){\\dashbox{15}(600,400){ }} ',
     $' \\put( 900,  40){\\makebox(0,0)[b]{\\LARGE $ 1-U_{\\min} $}}',
!>     $' \\put( 673,0){\\line(0,1){1200}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!-----------------------------------
!==================================================================

! ------------------------------------------------
! ----------------- BHLUMI -----------------------
! ------------------------------------------------
      jene=1
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
      JBHL = 30000 +2000
      kbhl = JBHL  +10000000
!------- normalize in nanobarns
      KeyGen = 3
      DO itr=1,3
! BHLUMI O(alf2)exp
        CALL cumhis(KeyGen,JBHL+340+itr,kbhl+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBHL+350+itr,kbhl+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBHL+360+itr,kbhl+360+itr) ! SICAL2
        CALL cumhis(KeyGen,JBHL+300+itr,kbhl+300+itr) ! SICAL
      ENDDO
! ------------------------------------------------
! ----------------- OLDBIS -----------------------
! ------------------------------------------------
      jene=2
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
      JBIS= 10000 +2000
      kbis= JBis +10000000
!------- normalize in nanobarns
      KeyGen = 1
      DO itr=1,3
        CALL cumhis(KeyGen,JBIS+300+itr,kbis+300+itr) ! SICAL
        CALL cumhis(KeyGen,JBIS+340+itr,kbis+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBIS+350+itr,kbis+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBIS+360+itr,kbis+360+itr) ! SICAL2
      ENDDO
! ------------------------------------------------
! ----------------- LUMLOG -----------------------
! ------------------------------------------------
      jene=3
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
!------- normalize in nanobarns
      JLOG= 20000 +4000
      klog= JLOG +10000000
      KeyGen = 2
      DO itr=1,3
! O(alf3-alf1) !!!
        CALL cumhis(KeyGen,JLOG+340+itr,klog+340+itr) ! BARE1
        CALL cumhis(KeyGen,JLOG+350+itr,klog+350+itr) ! CALO2
        CALL cumhis(KeyGen,JLOG+360+itr,klog+360+itr) ! SICAL2
        CALL cumhis(KeyGen,JLOG+300+itr,klog+300+itr) ! SICAL
      ENDDO
      LLOG= 20000 +2000     ! from JDA  O(alf1) for SABSPV
      mlog= LLOG +10000000
      DO itr=1,3
        CALL cumhis(KeyGen,LLOG+340+itr,mlog+340+itr) ! BARE1
        CALL cumhis(KeyGen,LLOG+350+itr,mlog+350+itr) ! CALO2
        CALL cumhis(KeyGen,LLOG+360+itr,mlog+360+itr) ! SICAL2
        CALL cumhis(KeyGen,LLOG+300+itr,mlog+300+itr) ! SICAL
      ENDDO
!==================================================================
!==================================================================
!             B-(O+L)
!==================================================================
! Born for Sical, wide and narrow
      THsic1 = .024d0
      THsic2 = .058d0
      Nthe = 16
      PAD = (THsic2-THsic1)/Nthe
      BorW(2)  = BORNB(CMSENE,THsic1   +PAD, THsic2   -PAD)
      BorN(2)  = BORNB(CMSENE,THsic1 +2*PAD, THsic2 -4*PAD)
!-------------------
      khyb= 50000 +2000
      kbpv= 60000 +2000
      kdif= khyb +10000000
      nlog= llog +20000000
      kdip= llog +30000000
      DO itr=1,3
         fac = 1/BorN(2)
         IF(itr .EQ. 1) fac = 1/BorW(2)
! (B-(O+L))/Born
         CALL gopera(kbis+300+itr,'+',klog+300+itr,khyb+300+itr,1d0,1d0)
         CALL gopera(khyb+300+itr,'-',kbhl+300+itr,kdif+300+itr,fac,fac)
!     (O(alf3)LL-O(alf1))LL +O(alf1)LL => O(alf3)LL
!     (Oldbis   -O(alf1)LL) /Born      => NLL/Born
!     O(alf3)LLexp * (NLL/Born)        => Corr
!     O(alf3)LLexp + Corr              => BhPV
!     (BhPV  - Bhlum4) /Born           => BhPdiff
         CALL gopera(klog+300+itr,'+',mlog+300+itr,nlog+300+itr,1d0,1d0)
         CALL gopera(kbis+300+itr,'-',mlog+300+itr,kbpv+300+itr,fac,fac)
         CALL gopera(nlog+300+itr,'*',kbpv+300+itr,kbpv+300+itr,1d0,1d0)
         CALL gopera(nlog+300+itr,'+',kbpv+300+itr,kbpv+300+itr,1d0,1d0)
         CALL gopera(kbpv+300+itr,'-',kbhl+300+itr,kdip+300+itr,fac,fac)
      ENDDO
!=========================================================
!  PLOTS PLOTS PLOTS PLOTS PLOTS
!=========================================================
!--------------------------------------------------------
      ymin=-3d-3
      ymax= 3d-3
      fmtx='f10.2'
      fmty='f10.3'
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT( 0)
      NOUFIG=11
      TeXfile   = 'BeMaFigs.tex'
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!=========================================================
!=========================================================
!-------------------------------------------------------
!   (BHLPV-BHLUMI)/Born  SICAL WW, NN, NW            !
!-------------------------------------------------------
      DO itr=1,3
      CALL gplcapt(cpsical2)
         CALL gidopt(kdip+300+itr,'ERRO')
         CALL gmimax(kdip+300+itr,ymin,ymax)
      ENDDO
      CALL gplot2(kdip+301,' ','*',dot    ,fmtx,fmty)
      CALL gplot2(kdip+302,'S','*',circle ,fmtx,fmty)
      CALL gplot2(kdip+303,'S','*',disc   ,fmtx,fmty)
      CALL gplabel(labex2)
!-----------------------------------------------!
!   ((O+L)-BHLUMI)/Born  SICAL WW, NN, NW            !
!-----------------------------------------------!
      CALL gpltit(' BeMaFigs B-(O+L) $')
      DO itr=1,3
         CALL gidopt(kdif+300+itr,'ERRO')
         CALL gmimax(kdif+300+itr,ymin,ymax)
      ENDDO
!-------------------------------------------------------
      CALL gplcapt(cpsical1)
      CALL gplot2(kdif+301,' ','*',dot    ,fmtx,fmty)
      CALL gplot2(kdif+302,'S','*',circle ,fmtx,fmty)
      CALL gplot2(kdif+303,'S','*',disc   ,fmtx,fmty)
      CALL gplabel(labex1)
!=========================================================
!   The end of plots writing
      CALL gplend
!=========================================================
      WRITE(   6,BXCLO)
      END

      SUBROUTINE subtra(id1,id2,id3)
!     ******************************
! subrats id2 from id1 and divides by id2
! Errors are only divided, not combined as in gopera!!!
!     *********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION bin1(100),bin2(100)
      DIMENSION err1(100),err2(100)

      CALL gunpak(id1,bin1,' ',0)
      CALL gunpak(id1,err1,'ERRO',0)
      CALL gunpak(id2,bin2,' ',0)
      DO i=1,100
         bin1(i)=(bin1(i)-bin2(i))/bin2(i)
         err1(i)=          err1(i)/bin2(i)
      ENDDO
      CALL gpak (id3,bin1)
      CALL gpake(id3,err1)
      END
