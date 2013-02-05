*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*------------------------------------------------------------------------------
      INTEGER    imax
      PARAMETER (imax=10000)
      REAL*8     xpar(imax)
      INTEGER    KeyDis
      REAL*8     ymin,ymax
*------------------------------------------------------------------------------
      INTEGER idb,idc,nbiv
      INTEGER iO0nll,iO0tech
      INTEGER iO0tot,mO0tot,iO0dif
      INTEGER iO1tot,mO1tot,iO1dif
      INTEGER iO2tot,mO2tot,iO2dif
      INTEGER iO2mO1,mO2mO1,mO2dO1
      INTEGER mO3tot
      INTEGER mO2vO1,mO3vO2
      INTEGER iO3best,mO3dan
*------------------------------------------------------------------------------
* basic histos 
      PARAMETER(
     $  idb =  50000,
     $  idc = 150000,
     $  mO0tot  = idc+71,            ! MC O(alf0) beta0=total
     $  mO1tot  = idc+72,            ! MC O(alf1) total
     $  mO2tot  = idc+73,            ! MC O(alf2) total
     $  mO3tot  = idc+74,            ! MC O(alf3) total
     $  mO2vO1  = idc+76,            ! MC O(alf2-alf1) from wt differences
     $  mO3vO2  = idc+77,            ! MC O(alf3-alf2) from wt differences
     $  iO0nll  = 9000000 + 400300,  ! SAN O(alf0) O(alf3)pragm
     $  iO0tot  = 9000000 + 300300,  ! SAN O(alf0) standard
     $  iO1tot  = 9000000 + 301301,  ! SAN O(alf1) standard
     $  iO2tot  = 9000000 + 302302,  ! SAN O(alf2) standard
     $  iO3best = 5000000 + 305302,  ! SAN O(alf3) the best 
* differences
     $  iO0tech = iO0tot  +2000000,  ! MC_O(alf0)-SAN_O(alf3)prag
     $  iO0dif  = iO0tot  +1000000,  ! MC_O(alf0)-SAN_O(alf0)
     $  iO1dif  = iO1tot  +1000000,  ! MC_O(alf1)-SAN_O(alf1)
     $  iO2dif  = iO2tot  +1000000,  ! MC_O(alf2)-SAN_O(alf2)
     $  iO2mO1  = 9000000 + 306306,  ! SAN O(alf2)-O(alf1) standard
     $  mO2mO1  = 9000000 + mO2tot,  ! MC_O(alf2) - MC_O(alf1) normal total
     $  mO2dO1  = 9000000 + mO2mO1,  ! (MC-SAN) O(alf2-alf1); diff. of differences!!
     $  mO3dan  = 5000000 + iO3best  ! MC_O(alf3) - SAN_best
     $ )
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,plus,box,dot
      PARAMETER (diamond ='\\makebox(0,0){\\Large $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\Large $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\Large $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (plus    ='\\makebox(0,0){\\Large $+$}')
      PARAMETER (box     ='\\makebox(0,0){\\Large $\\Box$}')
      PARAMETER (dot     ='\\circle*{10}')
*---------------------------------------------------------------------------
      CHARACTER*60  TeXfile
      CHARACTER*16  fmtx,fmty
*---------------------------------------------------------------------------
      CHARACTER*80 LabelBasic(8)
      DATA LabelBasic /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
* General KK logo
!!!! $'\\PaveLt{600}{1170}{\\color{green}\\Large $\\cal KK$ Monte Carlo 1999}',
     $'\\PaveLt{600}{1170}{\\color{green}',
     $'        \\large $\\cal KK$ MC 1999, S.Jadach, Z. W\\c{a}s, B.F.L. Ward}',
* Additional comment on x-axis
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
     $'} % -- end Pave',
     $'% end-of-label'/
* Labeling style 1
c     $'\\PaveL{ 40}{1000}{\\Huge ',
c     $' ${\\sigma\\left(s^\\prime_{\\min} \\right) \\over\\sigma_{_{\\rm{Born}}}}$}',
c     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
* Labeling style 2
c     $'\\PaveL{ 40}{1000}{\\Huge ',
c     $' ${\\sigma\\left(v<v_{\\rm{\\max}}\\right) \\over\\sigma_{_{\\rm{Born}}} }$ }',
c     $'\\PaveLb{700}{ 40}{\\huge $v_{\\rm{\\max}}$}',
*---------------------------------------------------------------------------

