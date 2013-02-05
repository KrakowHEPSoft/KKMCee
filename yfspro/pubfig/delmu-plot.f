      PROGRAM MAIN
*****************************
*     gmake delta-ps
*     gmake delmu-sigma1.eps
*     gmake delmu-sigma2.eps
*****************************
     IMPLICIT NONE
*---------------------------------------------------------------------------
      CHARACTER*60  TeXfile
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 star,diamond,circle,ring,times,disc,plus,box,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (ring    ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{30}')
      PARAMETER (plus    ='\\makebox(0,0){\\LARGE $+$}')
      PARAMETER (dot     ='\\circle*{15}')
*---------------------------------------------------------------------------
      CHARACTER*80 LabelBasic(2)
      DATA LabelBasic /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
ccc     $'\\Pave{',
ccc     $'\\PaveLt{600}{1170}{\\color{green}\\Large $\\cal KK$ Monte Carlo 1999}',
ccc     $'} % -- end Pave',
     $'% end-of-label'/
*---------------------------------------------------------------------------
      CHARACTER*80 Label1(13)
      DATA Label1 /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{  50}{1100}{\\Huge $\\sigma_{\\rm int.}$ [pb] }',
     $'\\PaveLb{900}{  40}{\\huge $\\sqrt{s}-M_{Z}$}',
*
     $'\\PaveLr{200}{ 280}{$\\circle*{15}\\quad$ }',
     $'\\PaveL{ 200}{ 280}{\\large $q$ CEEX $\\cal KK$  4.01}',
     $'\\PaveLr{200}{ 220}{\\color{green} $\\circle{30}\\quad$ }',
     $'\\PaveL{ 200}{ 220}{\\color{green}\\large  $\\mu$ CEEX $\\cal KK$ 4.01}',
     $'\\PaveLr{200}{ 160}{\\color[named]{Purple}\\LARGE$\\diamond\\;\\;$ }',
     $'\\PaveL{ 200}{ 160}{\\color[named]{Purple}\\large $\\mu$ KORALZ 1-st ord.}',
     $'\\PaveLr{200}{ 100}{\\color{red} $\\times\\;\\;$ }',
     $'\\PaveL{ 200}{ 100}{\\color{red}\\large   $\\mu$ ZFITTER/TOPAZ0 }',
     $'} % -- End Pave',
     $'% end-of-label'/
*---------------------------------------------------------------------------
      CHARACTER*80 Label2(13)
      DATA Label2 /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveLr{ -10}{1100}{\\Huge $\\sigma_{\\rm int}\\over\\sigma_{\\rm tot}$}',
     $'\\PaveLb{600}{  40}{\\huge $\\sqrt{s}-M_{Z}$}',
* 
     $'\\PaveLr{360}{1140}{$\\circle*{15}\\quad$ }',
     $'\\PaveL{ 360}{1140}{\\large $q$ CEEX $\\cal KK$ $\\times10$ }',
     $'\\PaveLr{360}{1080}{\\color{green} $\\circle{30}\\quad$ }',
     $'\\PaveL{ 360}{1080}{\\color{green}\\large $\\mu$ CEEX $\\cal KK$ }',
     $'\\PaveLr{360}{1020}{\\color[named]{Purple}\\LARGE$\\diamond\\;\\;$ }',
     $'\\PaveL{ 360}{1020}{\\color[named]{Purple}\\large $\\mu$ KORALZ 1-st ord.}',
     $'\\PaveLr{360}{ 960}{\\color{red} $\\times\\;\\;$ }',
     $'\\PaveL{ 360}{ 960}{\\color{red}\\large  $\\mu$ ZFITTER/TOPAZ0 }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*---------------------------------------------------------------------------
      INTEGER     nout
      INTEGER     npoint, i,j,k,ibin,nbin
      PARAMETER (nbin = 1000)
      DOUBLE PRECISION xmin,xmax
      DOUBLE PRECISION dCMSene(10),  SigG1(10), SigG1Int(10), dSigG1Int(10), SigZfInt(10)
      DOUBLE PRECISION hCMSene(10),  SihG1(10), SihG1Int(10), dSihG1Int(10)
      DOUBLE PRECISION SigMsInt(10), dSigMsInt(10)
      DOUBLE PRECISION zfit(nbin) 
      INTEGER          izfit
      DOUBLE PRECISION ceex(nbin), cint(nbin),dcint(nbin)
      INTEGER          iceex, icint
      DOUBLE PRECISION heex(nbin), hint(nbin),dhint(nbin)
      INTEGER          iheex, ihint
      DOUBLE PRECISION mint(nbin),dmint(nbin) ! Koralz/mustral
      INTEGER          imint
*----------------------------------------

      nout= 16
      OPEN( nout, FILE='output-delmu-plot')

      DO j=1,nbin
         zfit(j) = -1d19
         ceex(j) = -1d9
         cint(j) = -1d29
         dcint(j) = 0d0
         heex(j) = -1d9
         hint(j) = -1d29
         dhint(j) = 0d0
         mint(j) = -1d29
         dmint(j) = 0d0
      ENDDO

      xmin = -4d0
      xmax = +4d0
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                  read input from delmu program
      OPEN( 10, FILE='delmu.data')
      READ( 10,*)
      READ( 10,*)
      READ( 10,*) npoint
      WRITE(*,*)  npoint
      DO k=1, npoint
         READ ( 10,*)           
     $    dCMSene(k), SigG1(k), SigZfInt(k), SigG1Int(k), dSigG1Int(k), SigMsInt(k),dSigMsInt(k)
         WRITE(  *,'(20g12.4)') 
     $    dCMSene(k), SigG1(k), SigZfInt(k), SigG1Int(k), dSigG1Int(k), SigMsInt(k),dSigMsInt(k)
         ibin = (dCMSene(k)-xmin)/(xmax-xmin)*nbin +1
         IF(ibin.GE.1 .AND. ibin.LE.nbin) THEN
            ceex( ibin) =  SigG1(k)
            cint( ibin) =  SigG1Int(k)
            dcint(ibin) = dSigG1Int(k)
            zfit( ibin) =  SigZfInt(k)
            mint( ibin) =  SigMsInt(k)
            dmint(ibin) = dSigMsInt(k)
         ENDIF
      ENDDO
      CLOSE(10)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                  read input from delha program
      OPEN( 10, FILE='delhad.data')
      READ( 10,*)
      READ( 10,*)
      READ( 10,*) npoint
      WRITE(*,*)  npoint
      DO k=1, npoint
         READ ( 10,*)           hCMSene(k), SihG1(k), SihG1Int(k), dSihG1Int(k)
         WRITE(  *,'(20f12.4)') hCMSene(k), SihG1(k), SihG1Int(k), dSihG1Int(k)
         ibin = (hCMSene(k)-xmin)/(xmax-xmin)*nbin +1
         IF(ibin.GE.1 .AND. ibin.LE.nbin) THEN
            heex( ibin) = SihG1(k)
            hint( ibin) = SihG1Int(k)
            dhint(ibin) = dSihG1Int(k)
         ENDIF
      ENDDO
      CLOSE(10)
*/////////////////////////////////////////////////////////////////////////////////////////////
      izfit  = 1001
      iceex  = 1002
      icint  = 1003
      iheex  = 1004
      ihint  = 1005
      imint  = 1006
      CALL GLK_Book1( iceex,'mu  Sigma total     CEEX   $',nbin,xmin,xmax)
      CALL GLK_Book1( icint,'mu  Sigma interf.   CEEX   $',nbin,xmin,xmax)
      CALL GLK_Book1( izfit,'mu  Sigma interf. Zfiter   $',nbin,xmin,xmax)
      CALL GLK_Book1( iheex,'had Sigma total     CEEX   $',nbin,xmin,xmax)
      CALL GLK_Book1( ihint,'had Sigma interf.   CEEX   $',nbin,xmin,xmax)
      CALL GLK_Book1( imint,'had Sigma interf. KorZ/Must$',nbin,xmin,xmax)

      CALL GLK_Pak(  iceex, ceex)
      CALL GLK_Pak(  izfit, zfit)
      CALL GLK_Pak(  icint, cint)
      CALL GLK_Pake( icint,dcint)

      CALL GLK_Pak(  iheex, heex)
      CALL GLK_Pak(  ihint, hint)
      CALL GLK_Pake( ihint,dhint)

      CALL GLK_Pak(  imint, mint)
      CALL GLK_Pake( imint,dmint)

      fmtx='f10.2'
      fmty='f10.2'
      TeXfile   = 'delmu-sigma1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetYminYmax(icint,-3d0, 3d0)
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  icint, ' ','*',ring     ,fmtx,fmty)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  ihint, 'S','*',dot      ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  izfit, 'S','*',times     ,fmtx,fmty)
      CALL GLK_SetColor('\\color[named]{Purple}$')
      CALL GLK_plot2(  imint, 'S','*',diamond   ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label1)
      CALL GLK_PlEnd

      CALL GLK_Operat(izfit,'/',iceex, izfit, 1d0, 1d0) ! relative
      CALL GLK_Operat(icint,'/',iceex, icint, 1d0, 1d0) ! relative
      CALL GLK_Operat(ihint,'/',iheex, ihint, 1d1, 1d0) ! relative*10
      CALL GLK_Operat(imint,'/',iceex, imint, 1d0, 1d0) ! relative
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'delmu-sigma2.txp'
      CALL GLK_PlInitialize(2,TeXfile)
ccc      CALL GLK_SetYminYmax(icint,-3d-3, 3d-3) !!! for ibin=99
      CALL GLK_SetYminYmax(icint,-9d-3, 9d-3) !!! for ibin=10-90
ccc      CALL GLK_SetYminYmax(icint,-3d-2, 3d-2) !!! for ibin=1
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  icint, ' ','*',ring      ,fmtx,fmty)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  ihint, 'S','*',dot       ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  izfit, 'S','*',times     ,fmtx,fmty)
      CALL GLK_SetColor('\\color[named]{Purple}$')
      CALL GLK_plot2(  imint, 'S','*',diamond   ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label2)
      CALL GLK_PlEnd

      END

