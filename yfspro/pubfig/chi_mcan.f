*****************************************************************
* alias kmake='make -f KKMakefile'
*     kmake chi_mcan-dvi
*     kmake chi_mcan-ps
*****************************************************************
*     kmake chi_mcan-O2mca.eps
*     kmake chi_mcan-O2dif.eps
*     kmake chi_mcan-O2mO1.eps
*     kmake chi_mcan-O3dO2.eps
*     kmake chi_mcan-O0dif.eps
*     kmake chi_mcan-O3dan.eps
*     kmake chi_mcan-O0tech.eps
*----------------------------------------------------------------
      PROGRAM MAIN
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INCLUDE 'chi.h'
*
      CHARACTER*60  Tesnam, Dname
      CHARACTER*60  Hname, DumpFile
*--------------------------------------------------------------------
      ninp=  5
      nout= 16

      Tesnam    = 'chi_mcan'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

      Dname  = '../mix200/mix200.input.vmax=1_51M' ! Dec 97
*=====================================
* Read data, the same as in MC run
******CALL ReaDat(Dname,imax,xpar)
******CALL Semalib_Initialize(xpar)
* Read histograms from MC and prepare
      Hname  = './chi.hst'
      CALL GLK_ReadFile(Hname)
*=====================================================
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL Plot_O1mca           ! obsolete
      CALL Plot_O2mca
      CALL Plot_O1dif
      CALL Plot_O2dif
      CALL Plot_O3dan
*=====================================================
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL Plot_O2mO1
      CALL Plot_O3dO2
ccc      CALL Plot_O2dO1 ! not ploted
*=====================================================
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL Plot_O0dif
      CALL Plot_O0tech
ccc      CALL Plot_O2dO1 ! not ploted
ccc      CALL Plot_O3dan ! not ploted
*=========================================================================
* Write all histograms into dump file, for control/debug
******DumpFile = './chi.hst'
******CALL GLK_WriteFile(DumpFile)
*=================================
      CLOSE(nout)
      END


      SUBROUTINE Plot_O2mca
*------------------------------------
*     kmake chi_mcan-O2mca.eps
*------------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(11)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveLr{200}{700}{',
     $'             \\color{red}\\circle{30}\\; \\circle{30}\\; \\circle{30}}',
     $'\\PaveL{ 230}{700}{\\color{red}\\LARGE EEX2, Semi-Analyt.}',
     $'\\PaveLr{230}{600}{\\color{green}\\line(1,0){120}}',
     $'\\PaveL{ 230}{600}{\\color{green}\\LARGE EEX2, ${\\cal KK}$ M.C.}',
* Labeling style 1
     $'\\PaveL{ 40}{1000}{\\Huge ',
     $' ${\\sigma\\left(s^\\prime_{\\min} \\right) \\over\\sigma_{_{\\rm{Born}}}}$}',
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*              $_________|_________|_________|_________|_________|_________|_________|_________|
*                Ploting
      fmtx='f10.2'
      fmty='f10.2'
*===================================================================
*  First plot: O(alf2) MC
      TeXfile   = 'chi_mcan-O2mca.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2( mO2tot,' ',' ',dot    ,fmtx,fmty) ! MC
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2( iO2tot,'S','*',circle ,fmtx,fmty) ! analytical
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
*-----------------------------------------
      END

      SUBROUTINE Plot_O2dif
*-----------------------------------
*     kmake chi_mcan-O2dif.eps
*-----------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(7)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
c     $'\\PaveL{ 40}{1000}{\\LARGE ',
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma_{_{\\rm MC}}-\\sigma_{_{\\rm Sem}}\\over\\sigma_{_{\\rm Sem}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
* description label
cc     $'\\PaveLr{300}{800}{\\color{blue}\\line(1,0){200}}',
     $'\\PaveL{ 400}{800}{\\color{blue} EEX2: M.C. $-$ Sem.An. }',
*
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
*===================================================================
*                Ploting
      fmtx='f10.2'
      fmty='f10.2'
*===================================================================
*  Second plot: O(alf2) MC-SAN
      TeXfile   = 'chi_mcan-O2dif.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.025d0
      ymax =  0.025d0
      CALL GLK_Ymimax( iO2dif,ymin,ymax)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  iO2dif,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_O1mca
*--------------------------------
*     kmake chi_mcan-O1mca.eps
*---------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma_{_{\\rm MC}}-\\sigma_{_{\\rm Sem}}\\over\\sigma_{_{\\rm Sem}} }$}',!
*
     $'\\PaveLr{200}{700}{',
     $'             \\color{red}\\circle{30}\\; \\circle{30}\\; \\circle{30}}',
     $'\\PaveL{ 230}{700}{\\color{red}\\LARGE EEX1, Semi-Analyt.}',
     $'\\PaveLr{230}{600}{\\color{green}\\line(1,0){120}}',
     $'\\PaveL{ 230}{600}{\\color{green}\\LARGE EEX1, ${\\cal KK}$ M.C.}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
*===================================================================
*                Ploting
      fmtx='f10.2'
      fmty='f10.2'
*===================================================================
*  Third plot: O(alf1) MC
      TeXfile   = 'chi_mcan-O1mca.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2( mO1tot,' ',' ',dot    ,fmtx,fmty) ! MC
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2( iO1tot,'S','*',circle ,fmtx,fmty) ! analytical
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_O1dif
*------------------------------
*     kmake chi_mcan-O1dif.eps
*------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(5)
      DATA Label /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put(500,850){\\makebox(0,0)[l]{\\LARGE ${\\cal O}(\\alpha^1)_{exp}$ }}',
     $'\\put(500,750){\\makebox(0,0)[l]{\\LARGE M.C. $-$ S-Analyt. }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
*===================================================================
*                Ploting
      fmtx='f10.2'
      fmty='f10.2'
*===================================================================
*  Forth plot: O(alf1) MC-SAN
      TeXfile   = 'chi_mcan-O1dif.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.070d0
      ymax =  0.060d0
      CALL GLK_Ymimax( iO1dif,ymin,ymax)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  iO1dif,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_O3dan
*----------------------------
*     kmake chi_mcan-O3dan.eps
*----------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(7)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma_{_{\\rm MC}}-\\sigma_{_{\\rm Sem}}\\over\\sigma_{_{\\rm Sem}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveL{ 350}{800}{\\color[named]{Magenta} EEX3: M.C. $-$ Sem.An. }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
      fmtx='f10.2'
      fmty='f10.3'
*===================================================================
*  eighth plot: O(alf3) MC-SAN
      TeXfile   = 'chi_mcan-O3dan.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.020d0
      ymax =  0.020d0
      CALL GLK_Ymimax(mO3dan,ymin,ymax)
      CALL GLK_SetColor('\\color[named]{Magenta}$')
      CALL GLK_plot2( mO3dan,' ',' ',times ,fmtx,fmty)
      CALL GLK_print( mO3dan)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_O2mO1
*-----------------------------------
*     kmake chi_mcan-O2mO1.eps
*-----------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(11)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma^{^{(2)}}-\\sigma^{^{(1)}}\\over\\sigma^{^{(2)}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLr{200}{700}{',
     $'             \\color{red}\\circle{30}\\; \\circle{30}\\; \\circle{30}}',
     $'\\PaveL{ 230}{700}{\\color{red} EEX2$-$EEX1, Semi-Analyt.}',
     $'\\PaveLr{230}{600}{\\color{green}\\line(1,0){120}}',
     $'\\PaveL{ 230}{600}{\\color{green} EEX2$-$EEX1, ${\\cal KK}$ M.C.}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
      fmtx='f10.2'
      fmty='f10.3'
*===================================================================
*  Fifth plot: O(alf2)-O(alf1) MC and SAN
      TeXfile   = 'chi_mcan-O2mO1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_Ymimax( iO2mO1,-0.02d0, 0.07d0)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2( iO2mO1,' ','*',circle ,fmtx,fmty) ! analytical
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2( mO2mO1,'S',' ',dot    ,fmtx,fmty) ! MC direct
***   x-check only
***   CALL GLK_plot2( mO2vO1,'S','*',times  ,fmtx,fmty) ! MC diff. of wts
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_O2dO1
*------------------------------
*     kmake chi_mcan-O2dO1.eps
*------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(6)
      DATA Label /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put(500,850){\\makebox(0,0)[l]{\\LARGE ',
     $'                ${\\cal O}(\\alpha^2)-{\\cal O}(\\alpha^1)$ }}',
     $'\\put(500,750){\\makebox(0,0)[l]{\\LARGE M.C. $-$ S-Analyt. }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
      fmtx='f10.2'
      fmty='f10.3'
*-----------------------------------------
*  Sixth plot: O(alf2)-O(alf1), MC-SAN
      TeXfile   = 'chi_mcan-O2dO1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.060d0
      ymax =  0.060d0
      CALL GLK_Ymimax(mO2dO1,ymin,ymax)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2( mO2dO1,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_O3dO2
*-------------------------------
*     kmake chi_mcan-O3dO2.eps
*-------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(7)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma^{^{(3)}}-\\sigma^{^{(2)}}\\over\\sigma^{^{(2)}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveL{ 300}{800}{\\color{blue} EEX3$-$EEX2: ${\\cal KK}$  M.C. }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
      fmtx='f10.2'
      fmty='f10.4'
*-----------------------------------------
*  Seventh plot: O(alf3)-O(alf2), MC-SAN
      TeXfile   = 'chi_mcan-O3dO2.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_Ymimax(mO3vO2,-0.0011d0,+0.0011d0)
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2( mO3vO2,' ',' ',times ,fmtx,fmty)
      CALL GLK_print( mO3vO2)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END




      SUBROUTINE Plot_O0dif
*------------------------------
*     kmake chi_mcan-O0dif.eps
*------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(7)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma_{_{\\rm MC}}-\\sigma_{_{\\rm Sem}}\\over\\sigma_{_{\\rm Sem}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 600}{750}{\\color{red} EEX0 $-$ Sem.An. }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
*===================================================================
*                Ploting
      fmtx='f10.2'
      fmty='f10.2'
*===================================================================
*  O(alf1) MC-SAN
      TeXfile   = 'chi_mcan-O0dif.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.020d0
      ymax =  0.020d0
      CALL GLK_Ymimax( iO0dif,ymin,ymax)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  iO0dif,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_O0tech
*--------------------------------
*     kmake chi_mcan-O0tech.eps
*--------------------------------
      IMPLICIT NONE
      INCLUDE 'chi.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(7)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* y-axis label
     $'\\PaveL{ 40}{1000}{\\huge ',
     $'${\\sigma_{_{\\rm MC}}-\\sigma_{_{\\rm Sem}}\\over\\sigma_{_{\\rm Sem}} }$}',!
* x-axis label
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 600}{750}{\\color{red} EEX0 $-$ Sem.An.Best }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
*===================================================================
*                Ploting
      fmtx='f10.2'
      fmty='f10.3'
*===================================================================
*  O(alf0) MC-SAN BEST!!!
      TeXfile   = 'chi_mcan-O0tech.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      ymin = -0.005d0
      ymax =  0.005d0
      ymin = -0.02d0            ! for vmax=1
      ymax =  0.02d0            ! for vmax=1
      CALL GLK_Ymimax( iO0tech,ymin,ymax)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  iO0tech,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

