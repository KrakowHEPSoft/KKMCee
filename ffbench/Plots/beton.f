      PROGRAM MAIN
!***********************************
! To execute: make -f beton.makefile beton-eps
!***************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!---------------------------------------------------------------------- 
      Tesnam    = 'beton'
      TeXfile   = 'beton.tex'
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms
!--------------------------------------------------------
      lendan = 0
! LUMLOG
      lendan = lendan+1
!
      Dname(lendan)  = '../drun2/drun2.data.1122M'  ! 1995 Oct
      Hname(lendan)  = '../drun2/drun2.hst.1122M'   ! 1995 Oct

!==========================================================
      CALL Plbeta
!==========================================================
! ------------dumping histogram for control -------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL GRFILE(NOUTH,DNAME,'N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)

      CLOSE(NOUT)
      END




      SUBROUTINE Plbeta
!     *****************
! Basic plots BLUM2, OLDBIS, LUMLOG separately
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
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / TMINW,RAXIW,TMINN,RAXIN,VMAXE
      COMMON / TRANSR / TRANS,TRMIN,TRMAX
      COMMON / KEYDIS / KEYDIS
      CHARACTER*80 TITLE
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
! Parameters for tables/figures
      DIMENSION    idl(6)
      CHARACTER*16 capt(7)
      CHARACTER*8  fmt(3),fmtx,fmty
      LOGICAL gexist
!---------------------------------------------------------------------- 
      DIMENSION BorW(10),BorN(10)
!---------------------------------------------------------------------- 
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,plus,box,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (plus    ='\\makebox(0,0){\\LARGE $+$}')
      PARAMETER (box     ='\\makebox(0,0){\\Large $\\Box$}')
      PARAMETER (dot     ='\\circle*{10}')
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Latot2(60)
      DATA Latot2 /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(450,750){\\begin{picture}( 300,400)',
     $' \\put(    90,250){\\makebox(0,0)[l]{\\huge ',
     $'                     TOTAL }}',
     $' \\put(   -15,160){\\line(1,0){140}}',
     $' \\put(   200,160){\\makebox(0,0)[l]{\\LARGE ',
     $'                   MC $\\times 10^{-3}$}} ',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'                  ANL $\\times 10^{-3}$}}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE MC$-$ANL}} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R^{(2)}(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Lab0(60)
      DATA Lab0 /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(450,750){\\begin{picture}( 300,400)',
     $' \\put(    90,250){\\makebox(0,0)[l]{\\huge ',
     $'                     $\\bar{\\beta}_0^{(0)}$ }}',
     $' \\put(   -15,160){\\line(1,0){140}}',
     $' \\put(   200,160){\\makebox(0,0)[l]{\\LARGE ',
     $'                   MC $\\times 10^{-3}$}} ',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'                  ANL $\\times 10^{-3}$}}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE MC$-$ANL}} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Lab1(60)
      DATA Lab1 /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(450,750){\\begin{picture}( 300,400)',
     $' \\put(    90,250){\\makebox(0,0)[l]{\\huge ',
     $'                     $\\bar{\\beta}_1^{(2)}$ }}',
     $' \\put(   -15,160){\\line(1,0){140}}',
     $' \\put(   200,160){\\makebox(0,0)[l]{\\LARGE ',
     $'                   MC $\\times 10^{-2}$}} ',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'                  ANL $\\times 10^{-2}$}}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE MC$-$ANL}} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Lab2ul(60)
      DATA Lab2ul /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(450,750){\\begin{picture}( 300,400)',
     $' \\put(    90,280){\\makebox(0,0)[l]{\\huge ',
     $'                     $\\bar{\\beta}_{2UL}^{(2)}$ }}',
     $' \\put(   -15,160){\\line(1,0){140}}',
     $' \\put(   200,160){\\makebox(0,0)[l]{\\LARGE ',
     $'                   MC $\\times 10^{-1}$}} ',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'                  ANL $\\times 10^{-1}$}}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE MC$-$ANL}} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Lab2uu(60)
      DATA Lab2uu /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(450,750){\\begin{picture}( 300,400)',
     $' \\put(    90,280){\\makebox(0,0)[l]{\\huge ',
     $'                     $\\bar{\\beta}_{2UU}^{(2)}$ }}',
     $' \\put(   -15,160){\\line(1,0){140}}',
     $' \\put(   200,160){\\makebox(0,0)[l]{\\LARGE ',
     $'                   MC }} ',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'                  ANL }}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE MC$-$ANL}} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 Ladifba(60)
      DATA Ladifba /
     $'\\put(300,250){\\begin{picture}( 1200,1200)',
     $'%%%%%%%%',
     $'\\put(400,950){\\begin{picture}( 300,400)',
     $' \\put(     0,180){\\makebox(0,0)[l]{\\huge ',
     $'                     TOTAL (B)-(A)}}',
!
     $' \\multiput(0, 80)(60,0){3}{\\circle{30}}',
     $' \\put(   200, 80){\\makebox(0,0)[l]{\\LARGE ',
     $'   ${\\cal O}(\\alpha^{(2)})$ }}',
!
     $' \\multiput(0,  0)(60,0){3}{\\circle*{10}}',
     $' \\put(   200,  0){\\makebox(0,0)[l]{\\LARGE ',
     $'  ${\\cal O}(\\alpha^{(1)})$  $\\times 10^{-1}$ }} ',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $'\\put( -10, 1000){\\makebox(0,0)[r]{\\LARGE ',
     $'    $R^{(B-A)}(t,V_{_{\\max}})$ }}',
     $'\\put( 600,  40){\\makebox(0,0)[b]{\\huge $ V_{\\max} $}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------

!==================================================================
!     Restore data and histo (1)
!==================================================================
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
!==========================================================
! Normalization
!... Calculate minimum/maximum angles
      TMAXW  = ACOS(1-(1-COS(TMINW*PI/180))*RAXIW)*180/PI
      TMAXN  = ACOS(1-(1-COS(TMINN*PI/180))*RAXIN)*180/PI
      TH1W=TMINW*PI/180
      TH2W=TMAXW*PI/180
      TH1N=TMINN*PI/180
      TH2N=TMAXN*PI/180
!... Calculate narrow and wide Born
      BORNW  = BORNB(CMSENE,TH1W,TH2W)
      BORNN  = BORNB(CMSENE,TH1N,TH2N)
      write(NOUT,*) ' bornN =',bornN
      write(NOUT,*) ' bornW =',bornW
      write(   6,*) ' bornN =',bornN
      write(   6,*) ' bornW =',bornW
*** narrow and wide limits on transfer (born level)
      TRMINW =CMSENE**2*      (1D0-COS(TMINW*PI/180))/2D0
      TRMAXW =CMSENE**2*      (1D0-COS(TMAXW*PI/180))/2D0
      DTR = TRMAXW-TRMINW
      TRANS = (TRMAXW+TRMINW)/2
      write(NOUT,*) ' TRMAXW,TRMINW=',TRMAXW,TRMINW
      write(NOUT,*) ' TRANS =       ',TRANS
      write(   6,*) ' TRMAXW,TRMINW=',TRMAXW,TRMINW
      write(   6,*) ' TRANS =       ',TRANS
      DBORN = DSIG0(TRMINW,TRMAXW)
      BORN  = DBORN/DTR
      write(NOUT,*) ' DBorn,Born  = ',DBORN,BORN
!---------------------------------------------------------------
!---------------------------------------------------------------
      KeyGen = 3
      JDAX =  10000*KEYGEN +2000
      KDAX = 100000 +JDAX
      LDAX = 200000 +JDAX
      IDAX = 300000 +JDAX
!=========================================================
! Two lines
!=========================================================
      DO i=140,154
! betas
        CALL cumhis(KeyGen,JDAX+i,KDAX+i)
      ENDDO
      CALL cumhis(KeyGen,JDAX+ 42,KDAX+ 42)
      CALL cumhis(KeyGen,JDAX+ 41,KDAX+ 41)
!=========================================================
      fact = 1/dBorn
      CALL GINBO1(JDAX+140,TITLE,NBV,VMIN,VMAX)
!---------------------------------------------------------------
!     Analytical:  already divided by Born
!---------------------------------------------------------------
! [142] total  O(alf2)
      iMCtot2b =KDAX+142
      iANtot2b =LDAX+142
      iDItot2b =IDAX+142
      iMCtot2a =KDAX+ 42
      iDItot2  =IDAX+642
      KEYDIS=   2142
      CALL DFPLOT('CUGS',iANtot2b,'Total O(alf0) $',NBV,VMIN,VMAX)
      ymin=-0.0011d0
      ymax= 0.0011d0
      XMAG =1d0
      YMAG =0.001d0
      CALL gopera(iMCtot2b,'+',iMCtot2b,iMCtot2b,0d0,fact)  ! MC/Born
      CALL gopera(iMCtot2a,'+',iMCtot2a,iMCtot2a,0d0,fact)  ! MC/Born
      CALL gopera(iMCtot2b,'-',iMCtot2a,iDItot2 ,XMAG,XMAG) ! MC-MC
      CALL gopera(iMCtot2b,'-',iANtot2b,iDItot2b,XMAG,XMAG) ! MC-AN
      CALL gopera(iMCtot2b,'+',iMCtot2b,iMCtot2b ,0d0,YMAG) ! magnify MC
      CALL gopera(iANtot2b,'+',iANtot2b,iANtot2b ,0d0,YMAG) ! magnify AN
      CALL gmimax(iMCtot2b ,ymin,ymax)
      CALL gmimax(iANtot2b ,ymin,ymax)
      CALL gmimax(iDItot2b ,ymin,ymax)
      CALL gmimax(iDItot2  ,ymin,ymax)
      CALL gprint(iMCtot2b)
      CALL gprint(iANtot2b)
      CALL gprint(iDItot2b)
      CALL gprint(iDItot2)
!---------------------------------------------------------------
! [141] total  O(alf1), difference B-A only
      iMCtot1b =KDAX+141
      iMCtot1a =KDAX+ 41
      iDItot1  =IDAX+641
      CALL gopera(iMCtot1b,'+',iMCtot1b,iMCtot1b,0d0,fact)  ! MC/Born
      CALL gopera(iMCtot1a,'+',iMCtot1a,iMCtot1a,0d0,fact)  ! MC/Born
      YMAG =0.1d0
      CALL gopera(iMCtot1b,'-',iMCtot1a,iDItot1 ,YMAG,YMAG) ! MC-MC
      CALL gmimax(iDItot1  ,ymin,ymax)
      CALL gprint(iDItot1)
!---------------------------------------------------------------
! [140] beta0, 2-lines O(alf0) (= beta_zero)
      iMCbet0 =KDAX+140
      iANbet0 =LDAX+140
      iDIbet0 =IDAX+140
      KEYDIS=   2140
      CALL DFPLOT('CUGS',iANbet0,'Total O(alf0) $',NBV,VMIN,VMAX)
      ymin=-0.0011d0
      ymax= 0.0011d0
      XMAG =1d0
      YMAG =0.001d0
      CALL gopera(iMCbet0,'+',iMCbet0,iMCbet0,0d0,fact)  ! MC/Born
      CALL gopera(iMCbet0,'-',iANbet0,iDIbet0,XMAG,XMAG) ! MC-AN
      CALL gopera(iMCbet0,'+',iMCbet0,iMCbet0 ,0d0,YMAG) ! magnify MC
      CALL gopera(iANbet0,'+',iANbet0,iANbet0 ,0d0,YMAG) ! magnify AN
      CALL gmimax(iMCbet0 ,ymin,ymax)
      CALL gmimax(iANbet0 ,ymin,ymax)
      CALL gmimax(iDIbet0 ,ymin,ymax)
      CALL gprint(iMCbet0)
      CALL gprint(iANbet0)
      CALL gprint(iDIbet0)
!---------------------------------------------------------------
! [150] beta1u, 2-lines O(alf2)
      iMCbet1 =KDAX+150
      iANbet1 =LDAX+150
      iDIbet1 =IDAX+150
      KEYDIS=   2150
      CALL DFPLOT('CUGS',iANbet1,'Total O(alf0) $',NBV,VMIN,VMAX)
      ymin=-0.0011d0
      ymax= 0.0011d0
      XMAG =1d0
      YMAG =0.01d0
      CALL gopera(iMCbet1,'+',iMCbet1,iMCbet1,0d0,fact)  ! MC/Born
      CALL gopera(iMCbet1,'-',iANbet1,iDIbet1,XMAG,XMAG) ! MC-AN
      CALL gopera(iMCbet1,'+',iMCbet1,iMCbet1 ,0d0,YMAG) ! magnify MC
      CALL gopera(iANbet1,'+',iANbet1,iANbet1 ,0d0,YMAG) ! magnify AN
      CALL gmimax(iMCbet1 ,ymin,ymax)
      CALL gmimax(iANbet1 ,ymin,ymax)
      CALL gmimax(iDIbet1 ,ymin,ymax)
      CALL gprint(iMCbet1)
      CALL gprint(iANbet1)
      CALL gprint(iDIbet1)
!---------------------------------------------------------------
! [152] beta2ul, 2-lines O(alf2)
      iMCbt2ul =KDAX+152
      iANbt2ul =LDAX+152
      iDIbt2ul =IDAX+152
      KEYDIS=   2152
      CALL DFPLOT('CUGS',iANbt2ul,'Total O(alf0) $',NBV,VMIN,VMAX)
      ymin=-0.0011d0
      ymax= 0.0011d0
      XMAG =1d0
      YMAG =0.1d0
      CALL gopera(iMCbt2ul,'+',iMCbt2ul,iMCbt2ul,0d0,fact)  ! MC/Born
      CALL gopera(iMCbt2ul,'-',iANbt2ul,iDIbt2ul,XMAG,XMAG) ! MC-AN
      CALL gopera(iMCbt2ul,'+',iMCbt2ul,iMCbt2ul ,0d0,YMAG) ! magnify MC
      CALL gopera(iANbt2ul,'+',iANbt2ul,iANbt2ul ,0d0,YMAG) ! magnify AN
      CALL gmimax(iMCbt2ul ,ymin,ymax)
      CALL gmimax(iANbt2ul ,ymin,ymax)
      CALL gmimax(iDIbt2ul ,ymin,ymax)
      CALL gprint(iMCbt2ul)
      CALL gprint(iANbt2ul)
      CALL gprint(iDIbt2ul)
!---------------------------------------------------------------
! [153] beta2uu, 2-lines O(alf2)
      iMCbt2uu =KDAX+153
      iANbt2uu =LDAX+153
      iDIbt2uu =IDAX+153
      KEYDIS=   2153
      CALL DFPLOT('CUGS',iANbt2uu,'Total O(alf0) $',NBV,VMIN,VMAX)
      ymin=-0.0011d0
      ymax= 0.0011d0
      XMAG =1d0
      YMAG =1d0
      CALL gopera(iMCbt2uu,'+',iMCbt2uu,iMCbt2uu,0d0,fact)  ! MC/Born
      CALL gopera(iMCbt2uu,'-',iANbt2uu,iDIbt2uu,XMAG,XMAG) ! MC-AN
      CALL gopera(iMCbt2uu,'+',iMCbt2uu,iMCbt2uu ,0d0,YMAG) ! magnify MC
      CALL gopera(iANbt2uu,'+',iANbt2uu,iANbt2uu ,0d0,YMAG) ! magnify AN
      CALL gmimax(iMCbt2uu ,ymin,ymax)
      CALL gmimax(iANbt2uu ,ymin,ymax)
      CALL gmimax(iDIbt2uu ,ymin,ymax)
      CALL gprint(iMCbt2uu)
      CALL gprint(iANbt2uu)
      CALL gprint(iDIbt2uu)
!---------------------------------------------------------------
      CALL gprint(iDItot2)
      CALL gprint(iDItot1)

!=========================================================
!  Separate PLOTS for translation into eps files
!=========================================================
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-tot2.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDItot2b ,' ','*',dot     ,fmtx,fmty) ! N-N
      CALL gplot2(iMCtot2b ,'S',' ',' '     ,fmtx,fmty) ! N-N
      CALL gplot2(iANtot2b ,'S','*',circle  ,fmtx,fmty) ! N-N
      CALL gplabel(Latot2)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-beta0.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDIbet0 ,' ','*',dot     ,fmtx,fmty) ! N-N
      CALL gplot2(iMCbet0 ,'S',' ',' '     ,fmtx,fmty) ! N-N
      CALL gplot2(iANbet0 ,'S','*',circle  ,fmtx,fmty) ! N-N
      CALL gplabel(Lab0)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-beta1.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDIbet1 ,' ','*',dot     ,fmtx,fmty) ! N-N
      CALL gplot2(iMCbet1 ,'S',' ',' '     ,fmtx,fmty) ! N-N
      CALL gplot2(iANbet1 ,'S','*',circle  ,fmtx,fmty) ! N-N
      CALL gplabel(Lab1)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-bt2ul.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDIbt2ul ,' ','*',dot     ,fmtx,fmty) ! N-N
      CALL gplot2(iMCbt2ul ,'S',' ',' '     ,fmtx,fmty) ! N-N
      CALL gplot2(iANbt2ul ,'S','*',circle  ,fmtx,fmty) ! N-N
      CALL gplabel(Lab2ul)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-bt2uu.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDIbt2uu ,' ','*',dot     ,fmtx,fmty) ! N-N
      CALL gplot2(iMCbt2uu ,'S',' ',' '     ,fmtx,fmty) ! N-N
      CALL gplot2(iANbt2uu ,'S','*',circle  ,fmtx,fmty) ! N-N
      CALL gplabel(Lab2uu)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './beton-difba.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.4'
      CALL gplot2(iDItot1 ,' ','*',dot     ,fmtx,fmty) ! B-A
      CALL gplot2(iDItot2 ,'S','*',circle  ,fmtx,fmty) ! B-A
      CALL gplabel(Ladifba)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
      END


