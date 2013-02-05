*------------------------------------------------------
*     make Sli-IFImu-ps
*------------------------------------------------------
      PROGRAM MAIN
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  Dname
      CHARACTER*60  Hname, DumpFile

      m_FlagSem=0               ! desabling KKsem in SABANG
      m_FlagSem=1               ! enabling  KKsem in SABANG
*-------------------------------------------------------------------------------
c      Dname  = '../189GeV/189GeV.input'
c      Hname  = '../189GeV/pro.hst' ! current
*
      Dname  = '../189GeV/MuIFIon189GeV.input'
      Hname  = '../189GeV/MuIFIon189GeVreal.hst.39M' ! Jun 29 vmax =0.99
c      Hname  = '../189GeV/MuIFIon189GeV.hst.17M' ! May 26 vmax =0.99
c      Hname  = '../189GeV/MuIFIon189GeVreal.hst.23M' ! June24  vmax =0.999, the same as 17M
*------------
c      Dname  = '../200GeV/MuIFIon200GeV.input'
c      Hname  = '../200GeV/MuIFIon200GeV.hst.12M'
*------------
c      Dname  = '../206GeV/MuIFIon206GeV.input'
c      Hname  = '../206GeV/MuIFIon206GeV.hst.30M'      ! June 29 vmax =0.99
*-------------------------------------------------------------------------------
      m_out = 16
      CALL GLK_SetNout(m_out)

      OPEN(m_out,FILE='./PlotISR.output')
*===================================================================================================
      CALL GLK_ReadFile(Hname)  ! Read histograms from MC run
      CALL GLK_ListPrint(m_out) ! debug
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,m_imax,m_xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,m_imax,m_xpar)  ! Read user input
      CALL KK2f_Initialize( m_xpar)                ! Initialize generator with the production data
      CALL KKsem_Initialize(m_xpar)                ! Initialize semianalytical package
*=========================================================================
      CALL ISRprepare
      CALL TabIFI1
      CALL IFIfigSig
      CALL IFIfigAfb
*=========================================================================
* Write all histograms into dump file, for control
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=========================================================================
      CLOSE(m_out)
      END


      SUBROUTINE ISRprepare
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER  iMod
*-------------
      CALL SAB_MakeRen
* Semi-analytical
      CALL SAB_FilKKsem
* Data from ZFitter
***** CALL SAB_RDatZF1('./zf_2f_ForSJ_mu_Minv.log.B',m_kSABsem) ! old data file with bug
c      iMod = 2   ! IFI off
c      CALL SAB_RDatZF2('./zf_ZUATSM_muons_IFIoff.out',m_kSABsem,iMod) ! old data file with bug
c      iMod = 1   ! IFI ON
c      CALL SAB_RDatZF2('./zf_ZUATSM_muons_IFIon_.dat',m_kSABsem,iMod) ! old data file with bug
* Data from ZFitter, IFI with bug
c      CALL SAB_RDatZF1a('./zf_2f_ForSJ_prima_Test_B_muons.out.bug',m_kSABsem) ! new 13 june
* Data from ZFitter, IFI corrected
      CALL SAB_RDatZF1a('./zf_2f_ForSJ_prima_Test_B_muos.out',m_kSABsem) ! new 24 june
      END

      SUBROUTINE TabIFI1
*////////////////////////////////////////////////////////////////////////////////                             
*     gmake TabIFI-1.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  TeXfile
      INTEGER       nCols
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(8)
      CHARACTER*80  capt(9)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
      CHARACTER*6   Energy
      INTEGER       i,j,k,ikFf,iThe,iMod,nSel,iSel1,iSel2
      DOUBLE PRECISION CMSene, xObs(50),dObs(50),sysErr(50),PhysPrec
      DOUBLE PRECISION GLK_hie, GLK_hi
      DOUBLE PRECISION xSel(m_iSel1,m_iSel2),eSel(m_iSel1,m_iSel2)
      CHARACTER*32 TabLab1(21)
      DATA TabLab1 /  ' $ v<0.01$ ', ' $ v<0.05$ ',
     $ ' $ v<0.10$ ', ' $ v<0.15$ ', ' $ v<0.20$ ', ' $ v<0.25$ ',
     $ ' $ v<0.30$ ', ' $ v<0.35$ ', ' $ v<0.40$ ', ' $ v<0.45$ ',
     $ ' $ v<0.50$ ', ' $ v<0.55$ ', ' $ v<0.60$ ', ' $ v<0.65$ ',
     $ ' $ v<0.70$ ', ' $ v<0.75$ ', ' $ v<0.80$ ', ' $ v<0.85$ ',
     $ ' $ v<0.90$ ', ' $ v<0.95$ ', ' $ v<0.99$ ' /
      CHARACTER*32 TabLab1s(11)
      DATA TabLab1s / ' $ v<0.01$ ',
     $ ' $ v<0.10$ ', ' $ v<0.20$ ',
     $ ' $ v<0.30$ ', ' $ v<0.40$ ',
     $ ' $ v<0.50$ ', ' $ v<0.60$ ',
     $ ' $ v<0.70$ ', ' $ v<0.80$ ',
     $ ' $ v<0.90$ ', ' $ v<0.99$ ' /
*-------------------------------------------------------------------------------
      CMSene = m_xpar(1)
      IF( ABS(CMSene-189d0).LT.001) Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) Energy = '206GeV'
*-------------------------------------------------------------------------------m_nVmax
      ikFf =13
      iThe = 1  ! xsections
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  1,m_nVmax, 1d3, ix_Best)  ! KKMC IFI-off
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  1,m_nVmax, 1d3, ix_NoInt) ! KKMC IFI-on
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 7,  1,m_nVmax, 1d3, ix_Int)   ! IFI from wt-difs
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe,10,  1,m_nVmax, 1d3, ix_Semi)     ! KKsem
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 2,  1,m_nVmax, 1d0, ix_ZFNoInt)  ! ZF IFI-off
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 1,  1,m_nVmax, 1d0, ix_ZFIntExp) ! ZF IFI-exp
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 3,  1,m_nVmax, 1d0, ix_ZFInt)    ! ZF IFI-on
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,    1d0,1d0)   ! IFI
      CALL GLK_Operat(ix_NoInt,'-', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_ZFter,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      iThe = 2 ! asymetries
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  1,m_nVmax, 1d0, ia_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  1,m_nVmax, 1d0, ia_NoInt) !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 7,  1,m_nVmax, 1d0, ia_Int)   ! IFI from wt-difs
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe,10,  1,m_nVmax, 1d0, ia_Semi)  ! KKsem
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 2,  1,m_nVmax, 1d0, ia_ZFNoInt)  ! ZF IFI-off
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 1,  1,m_nVmax, 1d0, ia_ZFIntExp) ! ZF IFI-exp
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 3,  1,m_nVmax, 1d0, ia_ZFInt)    ! ZF IFI-on
      CALL GLK_Operat(ia_Best, '-', ia_NoInt, ia_IFI,    1d0,1d0)   ! IFI
* gain in error less then factor 2, not worth realy
      CALL GLK_Operat(ia_Int,  '*', ix_Int,   ia_Int,    1d0,1d0)   ! correct it!!!!!!
      CALL GLK_Operat(ia_Int,  '/', ix_NoInt, ia_Int,    1d0,1d0)   ! correct it!!!!!!
      CALL GLK_Operat(ia_NoInt, '*',ix_Int,   ia_Work,   1d0,1d0)   ! correct it!!!!!!
      CALL GLK_Operat(ia_Work,  '/',ix_NoInt, ia_Work,   1d0,1d0)   ! correct it!!!!!!
      CALL GLK_Operat(ia_Int,   '-',ia_Work,  ia_Int,    1d0,1d0)   ! correct it!!!!!!
*
      CALL GLK_Operat(ia_NoInt,'-', ia_Semi,  ia_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ia_ZFter,'-', ix_Semi,  ia_Diff2,  1d0,1d0)   ! difference1
*-------------------------------------------
      TeXfile   = 'TabIFI-1.txp'
cc      TeXfile   = 'YRtabMu'//Energy//'.tex'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
*===================================================================
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='{\\color{blue}$v_{\\max}$}'
      capt(2)='{\\color{blue}(a) ${\\cal KK}$sem Refer.}'
      capt(3)='{\\color{blue}(b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm IFIoff}$}'!
      capt(4)='{\\color{blue}(c) ZF IFIoff }'
      capt(5)='{\\color{blue}(d) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm IFIon}$}'!
      capt(6)='{\\color{blue}(e) ZF IFIon }'
      capt(7)='{\\color{blue}(f) ZF IFIexp }'
c[[[[
c      capt(7)='{\\color{blue}(f) IFI KKMC }' ! just for control !!!
c      capt(8)='{\\color{blue}(f) IFI KKMC }' ! just for control !!!
c]]]]
!----------------------------------------------------------
*                  SIGMA
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\mu^+\\mu^-)$, PRIMITIVE, at '//Energy
      idl(1) = ix_Semi
      idl(2) = ix_NoInt
      idl(3) = ix_ZFNoInt
      idl(4) = ix_Best
      idl(5) = ix_ZFInt
      idl(6) = ix_ZFIntExp
      nCols  = 6
c[[[[
c      idl(6) = ix_IFI
c      idl(7) = ix_Int
c      nCols  = 7
c]]]]
!----------------------------------------------------------
cc      CALL GLK_SetTabRan(1,21,1)
cc      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_SetTabRan(1,21,2)
      CALL GLK_SetTabLab(11,TabLab1s)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
*                  AFB
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\mu^+\\mu^-)$, PRIMITIVE, at '//Energy
      idl(1) = ia_Semi
      idl(2) = ia_NoInt
      idl(3) = ia_ZFNoInt
      idl(4) = ia_Best
      idl(5) = ia_ZFInt
      idl(6) = ia_ZFIntExp
c[[[[
c      idl(6) = ia_IFI
c      idl(7) = ia_Int
c      nCols  = 7
c]]]]
      CALL GLK_SetTabRan(1,21,2)
      CALL GLK_SetTabLab(11,TabLab1s)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_Int)
      CALL GLK_Delet(ix_ZFNoInt)
      CALL GLK_Delet(ix_ZFInt)
      CALL GLK_Delet(ix_ZFIntExp)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Int)
      CALL GLK_Delet(ia_ZFNoInt)
      CALL GLK_Delet(ia_ZFInt)
      CALL GLK_Delet(ia_ZFIntExp)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      END


      SUBROUTINE IFIfigAfb
*////////////////////////////////////////////////////////////////////////////////
*     make FigFSR-Mu.eps
*     make FigIFI-Mu.eps
*     make FigIFI-Mu1.eps
*     make FigIFI-Mu2.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
      DOUBLE PRECISION ymin,ymax
* relative to KKsem, IFI OFF, make FigFSR-MuAfb.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigFSR-MuAfb.txp'
      ymin =-0.04d0
      ymax = 0.04d0
      CALL FigFSRAfb(ikFf,Labx,TeXfile, ymin,ymax)
* KKMC-ZF, IFI OFF, make FigFSR-MuAfb1.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigFSR-MuAfb1.txp'
      ymin =-0.012d0
      ymax = 0.012d0
      CALL FigFSRAfb1(ikFf,Labx,TeXfile, ymin,ymax)
* IFI with respect ZFexp, make FigIFI-MuAfb1.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigIFI-MuAfb1.txp'
      ymin =-0.020d0
      ymax = 0.050d0
      CALL FigIFIAfb1(ikFf,Labx,TeXfile, ymin,ymax)
* IFI from each separately, make FigIFI-MuAfb2.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigIFI-MuAfb2.txp'
      CALL FigIFIAfb2(ikFf,Labx,TeXfile,-0.01d0,0.085d0)
      END


      SUBROUTINE FigFSRAfb(ikFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(12)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large ref = KKsem}',!
     $'\\PaveL{350}{1060}{\\large\\color{blue} $\\diamond$ KKMC,  IFIoff }',!
     $'\\PaveL{350}{1000}{\\large\\color{red}  $\\times$ ZF, IFIoff}',!
*
     $'\\PaveL{20}{1100}{\\large $A_{FB} - A_{FB}^{\\rm ref} $ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 2  ! AFB
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d0, ia_NoInt) ! KKMC IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d0, ia_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ia_ZFter) ! ZF IFIoff
      CALL GLK_Operat(ia_NoInt,'-', ia_Semi,  ia_Diff1,  1d0,1d0)   ! KKMC IFIoff
      CALL GLK_Operat(ia_ZFter,'-', ia_Semi,  ia_Diff2,  1d0,1d0)   ! ZF IFIoff
      CALL GLK_SetYminYmax(ia_Diff1,ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff2,ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  ia_Diff1, ' ','*', m_diamond  ,fmtx,fmty) ! KKMC
      CALL GLK_SetColor('\\color{red}\\thicklines$')
      CALL GLK_plot2(  ia_Diff2, 'S','*', m_times    ,fmtx,fmty) ! ZFitter
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Semi)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      END

      SUBROUTINE FigFSRAfb1(ikFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1000}{\\large\\color{red} $\\times$\\;\\;  ZF$-$KKMC, IFIoff}',!
*
     $'\\PaveL{20}{1100}{\\large $A_{FB}^{\\rm ZF} - A_{FB}^{\\cal KK} $ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 2  ! AFB
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d0, ia_NoInt) ! KKMC IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d0, ia_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ia_ZFter) ! ZF IFIoff
      CALL GLK_Operat(ia_ZFter,'-', ia_NoInt, ia_Diff1,  1d0,1d0)   ! ZF IFIoff
      CALL GLK_SetYminYmax(ia_Diff1,ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{red}\\thicklines$')
      CALL GLK_plot2(  ia_Diff1, ' ','*', m_times    ,fmtx,fmty) ! ZFitter
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Semi)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_Diff1)
      END


      SUBROUTINE FigIFIAfb1(iKFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe, idRef
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(12)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{20}{1100}{\\large $A_{FB}^{\\rm ZF} - A_{FB}^{\\rm KKMC} $ }',!
*
     $'%%\\PaveL{390}{1120}{\\large\\color{blue}  $\\diamond$ KKMC}',!
     $'\\PaveL{390}{1060}{\\large\\color{red}   $\\times$ ZF IFIon}',!
     $'\\PaveL{350}{1000}{\\large\\color{green} $\\circle{30}$\\;\\; ZF, IFIexp}',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 2  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 1, 1d0, ia_Best)  ! KKMC IFI-on
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d0, ia_NoInt) ! KKMC IFI-off
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ia_ZFter)    ! ZF IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 1, 1d0, ia_ZFIntExp) ! ZF IFI-exp
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 3, 1d0, ia_ZFInt)    ! ZF IFI-on
      idRef = ia_Best
      CALL GLK_Operat(ia_Best,    '-', idRef,  ia_Diff1,  1d0,1d0)   ! KKMC IFI
      CALL GLK_Operat(ia_ZFInt,   '-', idRef,  ia_Diff2,  1d0,1d0)   ! ZF IFI
      CALL GLK_Operat(ia_ZFIntExp,'-', idRef,  ia_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_SetYminYmax(ia_Diff1, ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff2, ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff3, ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  ia_Diff2, ' ','*', m_times    ,fmtx,fmty) ! ZF IFI
****
**      CALL GLK_SetColor('\\color{blue}\\thicklines$')
**      CALL GLK_plot2(  ia_Diff1, 'S','*', m_diamond  ,fmtx,fmty) ! KKMC IFI
****
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  ia_Diff3, 'S','*', m_circle   ,fmtx,fmty) ! ZF IFI-exp
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   ! Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Int)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_ZFInt)
      CALL GLK_Delet(ia_ZFIntExp)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      CALL GLK_Delet(ia_Diff3)
      END


      SUBROUTINE FigIFIAfb2(iKFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(12)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large\\color{blue} $\\diamond$ KKMC, IFI }',!
     $'\\PaveL{350}{1060}{\\large\\color{red}  $\\times$ ZF, IFI}',!
     $'\\PaveL{350}{1000}{\\large\\color{green} $\\circle{30}$\\; ZF, IFIexp}',! ZFexp
*
     $'\\PaveL{20}{1100}{\\Large $ \\delta A_{FB}^{\\rm IFI} $ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 2  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 1, 1d0, ia_Best)  ! KKMC IFI-on
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d0, ia_NoInt) ! KKMC IFI-off
cc      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 7, 1d0, ia_Int)   ! KKMC IFI-off
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ia_ZFter)    ! ZF IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 1, 1d0, ia_ZFIntExp) ! ZF IFI-exp
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 3, 1d0, ia_ZFInt)    ! ZF IFI-on
      CALL GLK_Operat(ia_Best,    '-', ia_NoInt,  ia_Diff1,  1d0,1d0)   ! KKMC IFI
      CALL GLK_Operat(ia_ZFInt,   '-', ia_ZFter,  ia_Diff2,  1d0,1d0)   ! ZF IFI
      CALL GLK_Operat(ia_ZFIntExp,'-', ia_ZFter,  ia_Diff3,  1d0,1d0)   ! ZF IFI-exp
cc      CALL GLK_SetYminYmax(ia_Int, ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff1, ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff2, ymin, ymax) !
      CALL GLK_SetYminYmax(ia_Diff3, ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  ia_Diff1, ' ','*', m_diamond  ,fmtx,fmty) ! KKMC IFI
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  ia_Diff2, 'S','*', m_times    ,fmtx,fmty) ! ZF IFI
****
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  ia_Diff3, 'S','*', m_circle   ,fmtx,fmty) ! ZF IFI-exp
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  i_Misc2,  'S','*', m_dot      ,fmtx,fmty) ! Misc
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   ! Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Int)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_ZFInt)
      CALL GLK_Delet(ia_ZFIntExp)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      CALL GLK_Delet(ia_Diff3)
      END



      SUBROUTINE IFIfigSig
*////////////////////////////////////////////////////////////////////////////////
*     make FigFSR-Mu.eps
*     make FigIFI-Mu.eps
*     make FigIFI-Mu1.eps
*     make FigIFI-Mu2.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
      DOUBLE PRECISION ymin,ymax
* relative to KKsem, IFI OFF, make FigFSR-Mu.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigFSR-Mu.txp'
      CALL FigFSR(ikFf,Labx,TeXfile, -0.012d0, 0.012d0)
* relative to KKsem, IFI ON, make FigIFI-Mu.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigIFI-Mu.txp'
      CALL FigIFI( ikFf,Labx,TeXfile, -0.020d0, 0.120d0)
* relative to ZFexp, IFI ON, make FigIFI-Mu1.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigIFI-Mu1.txp'
      CALL FigIFI1( ikFf,Labx,TeXfile, -0.02d0, 0.05d0)
* IFI alone from each calculation intependently, make FigIFI-Mu2.eps
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigIFI-Mu2.txp'
      CALL FigIFI2( ikFf,Labx,TeXfile,-0.015d0, 0.110d0,)
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      END


      SUBROUTINE FigFSR(ikFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(13)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large ref = KKsem}',!
     $'\\PaveL{350}{1060}{\\large\\color{blue} $\\diamond$ KKMC, IFIoff }',!
     $'\\PaveL{350}{1000}{\\large\\color{red} $\\times$ ZF, IFIoff}',!
*
     $'\\PaveL{40}{1100}{\\Huge ',
     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d3, ix_NoInt) ! KKMC IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ix_ZFter) ! ZF IFIoff
      CALL GLK_Operat(ix_NoInt,'-', ix_Semi,  ix_Diff1,  1d0,1d0)   ! KKMC IFIoff
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! KKMC IFIoff
      CALL GLK_Operat(ix_ZFter,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   ! ZF IFIoff
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! ZF IFIoff
      CALL GLK_SetYminYmax(ix_Diff1,ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff2,ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  ix_Diff1, ' ','*', m_diamond  ,fmtx,fmty) ! KKMC
      CALL GLK_SetColor('\\color{red}\\thicklines$')
      CALL GLK_plot2(  ix_Diff2, 'S','*', m_times    ,fmtx,fmty) ! ZFitter
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFter)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      END



      SUBROUTINE FigIFI(iKFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(14)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large ref = KKsem}',!
     $'\\PaveL{350}{1060}{\\large\\color{blue} $\\diamond$ KKMC, IFIon }',!
     $'\\PaveL{350}{1000}{\\large\\color{red} $\\times$ ZF, IFIon}',!
     $'\\PaveL{350}{ 940}{\\large\\color{green} $\\circle{30}$\\; ZF, IFIexp}',!
*
     $'\\PaveL{40}{1100}{\\Huge ',
     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 1, 1d3, ix_Best)  ! KKMC IFI-on
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 7, 1d3, ix_Int)   ! IFI
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 1, 1d0, ix_ZFIntExp)  ! ZF IFI-exp
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 3, 1d0, ix_ZFInt)     ! ZF IFI-on
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,    1d0,1d0)   ! IFI alone
      CALL GLK_Operat(ix_Best, '-', ix_Semi,  ix_Diff1,  1d0,1d0)   ! KKMC IFI-on
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! KKMC IFI-on
      CALL GLK_Operat(ix_ZFInt,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   ! ZF IFI-on
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! ZF IFI-on
      CALL GLK_Operat(ix_ZFIntExp,'-', ix_Semi,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_Operat(ix_Diff3,   '/', ix_Semi,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_SetYminYmax(ix_Diff1, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff2, ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  ix_Diff1, ' ','*', m_diamond  ,fmtx,fmty) ! KKMC
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  ix_Diff2, 'S','*', m_times    ,fmtx,fmty) ! ZF IFI-on
** ZFexp
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  ix_Diff3, 'S','*', m_circle   ,fmtx,fmty) ! ZF IFI-exp
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  i_Misc1,  'S','*', m_dot      ,fmtx,fmty) ! Misc
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   ! Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_IFI)
      CALL GLK_Delet(ix_Int)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFInt)
      CALL GLK_Delet(ix_ZFIntExp)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      CALL GLK_Delet(ix_Diff3)
      END



      SUBROUTINE FigIFI1(iKFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe, idRef
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(13)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large ref = KKMC}',!
     $'\\PaveL{350}{1060}{\\large\\color{red}  $\\times$ ZF, IFIon }',!
ccc     $'\\PaveL{350}{1000}{\\large\\color{blue} $\\diamond$ KKMC,  IFIon }',!
     $'\\PaveL{350}{1000}{\\large\\color{green} $\\circle{30}$\\; ZF, IFIexp}',!
*
     $'\\PaveL{40}{1100}{\\Huge ',
     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 1, 1d3, ix_Best)  ! KKMC IFI-on
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 7, 1d3, ix_Int)   ! IFI
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 1, 1d0, ix_ZFIntExp)  ! ZF IFI-exp
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 3, 1d0, ix_ZFInt)     ! ZF IFI-on
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,    1d0,1d0)   ! IFI alone
*
      idRef = ix_ZFIntExp
      idRef = ix_Best
      CALL GLK_Operat(ix_Best,    '-', idRef,  ix_Diff1,  1d0,1d0)   ! KKMC IFI-on
      CALL GLK_Operat(ix_Diff1,   '/', idRef,  ix_Diff1,  1d0,1d0)   ! KKMC IFI-on
      CALL GLK_Operat(ix_ZFInt,   '-', idRef,  ix_Diff2,  1d0,1d0)   ! ZF IFI-on
      CALL GLK_Operat(ix_Diff2,   '/', idRef,  ix_Diff2,  1d0,1d0)   ! ZF IFI-on
      CALL GLK_Operat(ix_ZFIntExp,'-', idRef,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_Operat(ix_Diff3,   '/', idRef,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_SetYminYmax(ix_Diff1, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff2, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff3, ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  ix_Diff2, ' ','*', m_times    ,fmtx,fmty) ! ZF IFI-on
**** KKMC is now reference
**      CALL GLK_SetColor('\\color{blue}\\thicklines$')
**      CALL GLK_plot2(  ix_Diff1, 'S','*', m_diamond  ,fmtx,fmty) ! KKMC
* ZFexp
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  ix_Diff3, 'S','*', m_circle   ,fmtx,fmty) ! ZF IFI-exp
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   ! Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_IFI)
      CALL GLK_Delet(ix_Int)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFInt)
      CALL GLK_Delet(ix_ZFIntExp)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      CALL GLK_Delet(ix_Diff3)
      END


      SUBROUTINE FigIFI2(iKFf,Labx,TeXfile,ymin,ymax)
*////////////////////////////////////////////////////////////////////////////////                             
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*---------------------------------------------------------------------------
      CHARACTER*16     fmtx,fmty
      INTEGER          iThe
      DOUBLE PRECISION ymin,ymax
*---------------------------------------------------------------------------
      CHARACTER*80 Label(13)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $' ???  ',
     $'\\PaveL{350}{1120}{\\large\\color{blue} $\\diamond$ KKMC, IFI }',!
     $'\\PaveL{350}{1060}{\\large\\color{red} $\\times$ ZF, IFI}',!
     $'\\PaveL{350}{1000}{\\large\\color{green} $\\circle{30}$\\; ZF, IFIexp}',! ZFexp
*
     $'\\PaveL{30}{1100}{\\Huge ',
     $'     ${\\Delta\\sigma_{_{\\rm IFI}} \\over \\sigma}$ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',!
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      WRITE(*,*) '===>>> FigFSR: ikFf=',ikFf
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 1, 1d3, ix_Best)  ! KKMC IFI-on
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d3, ix_NoInt) ! KKMC IFI-off
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 7, 1d3, ix_Int)   ! KKMC IFI-off
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ix_ZFter)    ! ZF IFIoff
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 1, 1d0, ix_ZFIntExp) ! ZF IFI-exp
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 3, 1d0, ix_ZFInt)    ! ZF IFI-on
      CALL GLK_Operat(ix_Best, '-', ix_NoInt,  ix_Diff1,  1d0,1d0)   ! KKMC IFI
      CALL GLK_Operat(ix_Diff1,'/', ix_NoInt,  ix_Diff1,  1d0,1d0)   ! KKMC IFI
      CALL GLK_Operat(ix_Int,  '/', ix_NoInt,  ix_Int,    1d0,1d0)   ! KKMC IFI
      CALL GLK_Operat(ix_ZFInt,   '-', ix_ZFter,  ix_Diff2,  1d0,1d0)   ! ZF IFI
      CALL GLK_Operat(ix_Diff2,   '/', ix_ZFter,  ix_Diff2,  1d0,1d0)   ! ZF IFI
      CALL GLK_Operat(ix_ZFIntExp,'-', ix_ZFter,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_Operat(ix_Diff3,   '/', ix_ZFter,  ix_Diff3,  1d0,1d0)   ! ZF IFI-exp
      CALL GLK_SetYminYmax(ix_Int, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff1, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff2, ymin, ymax) !
      CALL GLK_SetYminYmax(ix_Diff3, ymin, ymax) !
*--------------------------------------------------------------------------
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  ix_Int, ' ','*', m_diamond  ,fmtx,fmty) ! KKMC IFI
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  ix_Diff2, 'S','*', m_times    ,fmtx,fmty) ! ZF IFI
**** ZFexp not shown
      CALL GLK_SetColor('\\color{green}$')
      CALL GLK_plot2(  ix_Diff3, 'S','*', m_circle   ,fmtx,fmty) ! ZF IFI-exp
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  i_Misc1,  'S','*', m_dot      ,fmtx,fmty) ! Misc
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd   ! Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_Int)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFter)
      CALL GLK_Delet(ix_ZFInt)
      CALL GLK_Delet(ix_ZFIntExp)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      CALL GLK_Delet(ix_Diff3)
      END
