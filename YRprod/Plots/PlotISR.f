*------------------------------------------------------
*     gmake Sli-ISRquks-ps
*     gmake TabISR-1.eps
*     gmake TabISR-2.eps
*     gmake FigISR-Qs.eps
*     gmake FigISR-Mu.eps
*     gmake FigISR-Down.eps
*     gmake FigISR-Up.eps
*     gmake FigISR-Stran.eps
*     gmake FigISR-Charm.eps
*     gmake FigISR-Bottom.eps
*  New!!!!        AFB,     gmake TabISR-2c.eps
*  New!!!!        IFI,     gmake TabISR-2d.eps
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
*---
      Dname  = '../189GeV/IncISRonly189GeV.input.140M'     !(N) Xcheck, Sept.2002
      Hname  = '../189GeV/IncISRonly189GeV.hst.140M'       !(N) quarks+mu
*---
c      Dname  = '../189GeV/IncISRonly189GeV.input.173M'     !(N) Xcheck, Oct.2001
c      Hname  = '../189GeV/IncISRonly189GeV.hst.173M'       !(N) quarks+mu
*---
c      Dname  = '../189GeV/IncISRonly189GeV.input'
c      Hname  = '../189GeV/IncISRonly189GeV.hst.217M'
*---
c      Dname  = '../200GeV/IncISRonly200GeV.input'
c      Hname  = '../200GeV/IncISRonly200GeV.hst.151M'
*---
c      Dname  = '../206GeV/IncISRonly206GeV.input'
c      Hname  = '../206GeV/IncISRonly206GeV.hst.239M'
*-------------------------------------------------------------------------------
c      Dname  = '../189GeV/IncIFIon189GeV.input'   ! SPECIAL IFI ON
c      Hname  = '../189GeV/IncIFIon189GeV.hst.20M' ! SPECIAL IFI ON
*---
*      Dname  = '../206GeV/IncIFIon206GeV.input'   !!! SPECIAL IFI ON
*      Hname  = '../206GeV/IncIFIon206GeV.hst.15M' !!! SPECIAL IFI ON new
c      Hname  = '../206GeV/IncIFIon206GeV.hst.17M' ! SPECIAL IFI ON old, obsolete
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
      CALL TabISR1
      CALL ISRfig1
      CALL TabISR2
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
      INTEGER iMod,iBox
*-------------
      CALL SAB_MakeRen
* Semi-analytical
      CALL SAB_FilKKsem
*------------------------------------------------------------------------
* Data from ZFitter, old incorrect version
****  CALL SAB_RDatZF1('./zf_2f_ForSJ_1.log.A',m_kSABsem) !
c      iMod = 2   ! IFI off
c      CALL SAB_RDatZF2('./zf_For_SJ_quarks_IFIoff_Test_A.out',m_kSABsem,iMod)  !
c      CALL SAB_RDatZF2('./zf_For_SJ_hadrons_IFIoff_Test_A.out',m_kSABsem,iMod) ! total
c      CALL SAB_RDatZF2('./zf_For_SJ_muons_IFIoff_Test_A.out',m_kSABsem,iMod)   ! muons
c      CALL SAB_RDatZF2b('./zf_For_SJ_allcha_IFIoff_Test_A.out',m_kSABsem,iMod)  ! all, ibox=0,2, obsolete
*------------------------------------------------------------------------
* Dima says it is important plot, shows that running of couplings is important for ZRR.
      iMod = 2   ! IFI off
      iBox = 0
**      CALL SAB_RDatZF2b('./zf_2f_ForSJ_prima_Test_A_leptns.out',m_kSABsem,iMod,iBox)  ! all, ibox=0,2, obsolete
**      CALL SAB_RDatZF2b('./zf_2f_ForSJ_prima_Test_A_quarks.out',m_kSABsem,iMod,iBox)  ! all, ibox=0,2, obsolete
*------------------------------------------------------------------------
* This is finally correct... all, ibox=0,2, new
      iMod = 2   ! IFI off
      iBox=0
      iBox=2
      CALL SAB_RDatZF2b('./zf_2f_ForSJ_prima_Test_A.out.new',m_kSABsem,iMod,iBox)     ! all, ibox=0,2, new
c[[[[[[[[[[[
C   see also [[[ ]]] in SabAng.f
c      CALL SAB_RDatZF1b('./zf_2f_ForSJ_primBORN_CF3.log.dima',m_kSABsem)   ! Simplistic
c]]]]]]]]]]]
      END


      SUBROUTINE TabISR1
*////////////////////////////////////////////////////////////////////////////////                                
*     gmake TabISR-1.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
*------------------------------------------------------------------
! Parameters for tables
      CHARACTER*60  TeXfile
      INTEGER       nCols
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
*------------------------------------------------------------------
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
*-------------------------------------------------------------------------------
      ikFf =7
      iThe = 1  ! xsections
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d3, ix_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d3, ix_NoInt) !
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe,10,  m_iSel1,m_iSel2, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ix_ZFter) ! ZFter
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,    1d0,1d0)   ! IFI
      CALL GLK_Operat(ix_NoInt,'-', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_ZFter,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      iThe = 2 ! asymetries
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d0, ia_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ia_NoInt) !
      CALL GLK_Operat(ia_Best, '-', ia_NoInt, ia_IFI,    1d0,1d0)   ! IFI
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe,10,  m_iSel1,m_iSel2, 1d0, ia_Semi)  ! KKsem
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ia_ZFter) ! ZFter
      CALL GLK_Operat(ia_NoInt,'-', ia_Semi,  ia_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ia_ZFter,'-', ia_Semi,  ia_Diff2,  1d0,1d0)   ! difference1
*-------------------------------------------
      TeXfile   = 'TabISR-1.txp'
cc      TeXfile   = 'YRtabMu'//Energy//'.tex'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
*===================================================================
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
!----------------------------------------------------------
*         SIGMA
!----------------------------------------------------------
      capt(1)='{\\color{blue}$v_{\\max}$}'
      capt(2)='{\\color{blue}(a) ${\\cal KK}$sem Refer.}'
      capt(3)='{\\color{blue}(b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm intOFF}$}'!
      capt(4)='{\\color{blue}(c) ZFter }'
      capt(5)='{\\color{blue}(b-a)/a }'
      capt(6)='{\\color{blue}(c-a)/a }'
      idl(1) = ix_Semi
      idl(2) = ix_NoInt
      idl(3) = ix_ZFter
      idl(4) = ix_Diff1
      idl(5) = ix_Diff2
      nCols  = 5
      Mcapt = '$\\sigma(q\\bar{q})$, PRIMITIVE, at '//Energy
***      CALL GLK_SetTabRan(1,21,1)
***      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_SetTabRan(1,21,2)
      CALL GLK_SetTabLab(11,TabLab1s)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
*          AFB
!----------------------------------------------------------
      idl(1) = ia_Semi
      idl(2) = ia_NoInt
      idl(3) = ia_ZFter
      idl(4) = ia_Diff1
      idl(5) = ia_Diff2
      nCols  = 5
      Mcapt = '$A_{FB}(q\\bar{q})$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,2)
      CALL GLK_SetTabLab(11,TabLab1s)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      END

      SUBROUTINE ISRfig1
*////////////////////////////////////////////////////////////////////////////////                                
*     gmake FigISR-Qs.eps
*     gmake FigISR-Mu.eps
*     gmake FigISR-Down.eps
*     gmake FigISR-Up.eps
*     gmake FigISR-Stran.eps
*     gmake FigISR-Charm.eps
*     gmake FigISR-Bottom.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          ikFf
      CHARACTER*60     TeXfile
      CHARACTER*80     Labx
*
      iKFf =7
      Labx = '\\PaveL{800}{1120}{\\Large Quarks, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Qs.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =13
      Labx = '\\PaveL{800}{1120}{\\Large $f=\\mu$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Mu.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =1
      Labx = '\\PaveL{800}{1120}{\\Large $f=d$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Down.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =2
      Labx = '\\PaveL{800}{1120}{\\Large $f=u$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Up.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =3
      Labx = '\\PaveL{800}{1120}{\\Large $f=s$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Stran.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =4
      Labx = '\\PaveL{800}{1120}{\\Large $f=c$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Charm.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
      iKFf =5
      Labx = '\\PaveL{800}{1120}{\\Large $f=b$, '//m_Energy//'}'!
      TeXfile   = 'FigISR-Bottom.txp'
      CALL ISRfigOne(ikFf,Labx,TeXfile)
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      END


      SUBROUTINE ISRfigOne(ikFf,Labx,TeXfile)
*////////////////////////////////////////////////////////////////////////////////                                
*     gmake FigISR-Qs.eps
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
      DOUBLE PRECISION ymax,ymin
*---------------------------------------------------------------------------
      CHARACTER*80 Label(13)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $'\\PaveL{750}{1120}{\\Large Quarks}',!
     $'\\PaveL{350}{1120}{\\large ref = KKsem}',!
     $'\\PaveL{350}{1060}{\\large\\color{blue} $\\diamond$ KKMC, IFIoff }',!
     $'\\PaveL{350}{1000}{\\large\\color{red} $\\times$ Zfitter, IFIoff}',!
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
*
      Label(3) = Labx
*-------------------------------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetHistV(m_kSABren,ikFf,iThe, 2, 1d3, ix_NoInt) !
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe,10, 1d3, ix_Semi)  ! KKsem
      CALL SAB_GetHistV(m_kSABsem,ikFf,iThe, 2, 1d0, ix_ZFter) ! ZFter
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,    1d0,1d0)   ! IFI
      CALL GLK_Operat(ix_NoInt,'-', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_ZFter,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference1
      ymax = 0.030d0
      ymin =-0.030d0
      ymax = 0.011d0
      ymin =-0.011d0
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
      CALL GLK_Delet(ix_IFI)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFter)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      END


      SUBROUTINE TabISR2
*////////////////////////////////////////////////////////////////////////////////
*     gmake TabISR-2.eps
*     gmake TabISR-2b.eps
*     gmake TabISR-2c.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  TeXfile
*------------------------------------------------------------------
      CHARACTER*80  mcapt
      CHARACTER*6   Energy
      INTEGER       iSel
      DOUBLE PRECISION CMSene
*-------------------------------------------------------------------------------
      CMSene = m_xpar(1)
      IF( ABS(CMSene-189d0).LT.001) Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) Energy = '206GeV'
*=================================================================================
      TeXfile   = 'TabISR-2.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
!----------------------------------------------------------
      iSel = 3  ! v<0.10
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.10$, }'//Energy
      CALL TabISR2sigma(iSel,' ',Mcapt)
!----------------------------------------------------------
      iSel = 5  ! v<0.20
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.20$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 15  ! v<0.70
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.70$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 19  ! v<0.90
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.90$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
      TeXfile   = 'TabISR-2b.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
!----------------------------------------------------------
      iSel = 1  ! v<0.01
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.01$, }'//Energy
      CALL TabISR2sigma(iSel,' ',Mcapt)
!----------------------------------------------------------
      iSel = 7  ! v<0.30
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.30$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 17  ! v<0.80
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.80$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 21  ! v<0.99
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.99$, }'//Energy
      CALL TabISR2sigma(iSel,'S',Mcapt)
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
*         New!!!!        AFB,     gmake TabISR-2c.eps
*=================================================================================
      TeXfile   = 'TabISR-2c.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
!----------------------------------------------------------
      iSel = 1  ! v<0.01
      Mcapt ='{\\color{red}$A_{FB}(v_{\\max})$[pb], $v_{\\max}=0.01$, }'//Energy
      CALL TabISR2afb(iSel,' ',Mcapt)
!----------------------------------------------------------
      iSel = 3  ! v<0.10
      Mcapt ='{\\color{red}$A_{FB}(v_{\\max})$[pb], $v_{\\max}=0.10$, }'//Energy
      CALL TabISR2afb(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 5  ! v<0.20
      Mcapt ='{\\color{red}$A_{FB}(v_{\\max})$[pb], $v_{\\max}=0.20$, }'//Energy
      CALL TabISR2afb(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 7  ! v<0.30
      Mcapt ='{\\color{red}$A_{FB}(v_{\\max})$[pb], $v_{\\max}=0.30$, }'//Energy
      CALL TabISR2afb(iSel,'S',Mcapt)
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
*         New!!!!        IFI,     gmake TabISR-2d.eps
*=================================================================================
      TeXfile   = 'TabISR-2d.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
!----------------------------------------------------------
      iSel = 1  ! v<0.01
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$[pb], $v_{\\max}=0.01$, }'//Energy
ccc      iSel = 3  ! v<0.10
ccc      Mcapt ='{\\color{red}ISR$\\otimes$FSR, $v_{\\max}=0.10$, }'//Energy
      CALL TabISR2ifi(iSel,' ',Mcapt)
!----------------------------------------------------------
      iSel = 5  ! v<0.20
      Mcapt ='{\\color{red}ISR$\\otimes$FSR, $v_{\\max}=0.20$, }'//Energy
      CALL TabISR2ifi(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 7  ! v<0.30
      Mcapt ='{\\color{red}ISR$\\otimes$FSR, $v_{\\max}=0.30$, }'//Energy
      CALL TabISR2ifi(iSel,'S',Mcapt)
!----------------------------------------------------------
      iSel = 19  ! v<0.90
      Mcapt ='{\\color{red}ISR$\\otimes$FSR, $v_{\\max}=0.90$, }'//Energy
      CALL TabISR2ifi(iSel,'S',Mcapt)
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
      END

      SUBROUTINE TabISR2sigma(iSel,Char,Mcapt)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER       iSel
      CHARACTER*1 Char
      CHARACTER*80  Mcapt
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      INTEGER       iThe,nCols
      CHARACTER*32  TabLab1(5), TabLab2(1), TabLab3(1)
      DATA TabLab1 / ' $d$ ', ' $u$ ', ' $s$ ', ' $c$ ', ' $b$ '/ !
      DATA TabLab2 / ' $all$' /
      DATA TabLab3 / ' $\\mu$' /
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='{\\color{blue} $f$}'
      capt(2)='{\\color{blue} (a) ${\\cal KK}$sem }'
      capt(3)='{\\color{blue} (b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm intOFF}$}' !
      capt(4)='{\\color{blue} (c) Zfitter 6.x }'
      capt(5)='{\\color{blue} (b-a)/a }'
      capt(6)='{\\color{blue} (c-a)/a }'
      idl(1) = ix_Semi
      idl(2) = ix_Best
      idl(3) = ix_ZFter
      idl(4) = ix_Diff1
      idl(5) = ix_Diff2
      nCols  = 5
!----------------------------------------------------------
      iThe = 1  ! xsections
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  1,iSel, 1d3, ix_Best) ! KKMC wtbest
      CALL SAB_GetKFHist(m_kSABsem,m_iKFf1,m_iKFf2, iThe, 10,iSel, 1d3, ix_Semi) ! KKsem semianalytical
      CALL SAB_GetKFHist(m_kSABsem,m_iKFf1,m_iKFf2, iThe,  2,iSel, 1d0, ix_ZFter) ! ZFter
      CALL GLK_Operat(ix_ZFter,'-', ix_Semi,  ix_Diff2,  1d0,1d0)   !
      CALL GLK_Operat(ix_Diff2,'/', ix_Semi,  ix_Diff2,  1d0,1d0)   ! difference2
      CALL GLK_Operat(ix_Best, '-', ix_Semi,  ix_Diff1,  1d0,1d0)   !
      CALL GLK_Operat(ix_Diff1,'/', ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
      CALL GLK_SetTabRan(1,5,1)
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,Char,'R',' ')
      CALL GLK_SetTabRan(7,7,1)
      CALL GLK_SetTabLab(1,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_SetTabRan(13,13,1)
      CALL GLK_SetTabLab(1,TabLab3)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_Semi)
      CALL GLK_Delet(ix_ZFter)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      END

      SUBROUTINE TabISR2afb(iSel,Char,Mcapt)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER       iSel
      CHARACTER*1 Char
      CHARACTER*80  Mcapt
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      INTEGER       iThe,nCols
      CHARACTER*32  TabLab1(5), TabLab2(1), TabLab3(1)
      DATA TabLab1 / ' $d$ ', ' $u$ ', ' $s$ ', ' $c$ ', ' $b$ '/ !
      DATA TabLab2 / ' $all$' /
      DATA TabLab3 / ' $\\mu$' /
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='{\\color{blue} $f$}'
      capt(2)='{\\color{blue} (a) ${\\cal KK}$sem }'
      capt(3)='{\\color{blue} (b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm intOFF}$}' !
      capt(4)='{\\color{blue} (c) Zfitter 6.x }'
      capt(5)='{\\color{blue} (b-a)/a }'
      capt(6)='{\\color{blue} (c-a)/a }'
      idl(1) = ia_Semi
      idl(2) = ia_Best
      idl(3) = ia_ZFter
      idl(4) = ia_Diff1
      idl(5) = ia_Diff2
      nCols  = 5
!----------------------------------------------------------
      iThe = 2  ! afb
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  1,iSel, 1d0, ia_Best) ! KKMC wtbest
      CALL SAB_GetKFHist(m_kSABsem,m_iKFf1,m_iKFf2, iThe, 10,iSel, 1d0, ia_Semi) ! KKsem semianalytical
      CALL SAB_GetKFHist(m_kSABsem,m_iKFf1,m_iKFf2, iThe,  2,iSel, 1d0, ia_ZFter) ! ZFter
      CALL GLK_Operat(ia_ZFter,'-', ia_Semi,  ia_Diff2,  1d0,1d0)   !
      CALL GLK_Operat(ia_Best, '-', ia_Semi,  ia_Diff1,  1d0,1d0)   !
      CALL GLK_SetTabRan(1,5,1)
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,Char,'R',' ')
      CALL GLK_SetTabRan(7,7,1)
      CALL GLK_SetTabLab(1,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_SetTabRan(13,13,1)
      CALL GLK_SetTabLab(1,TabLab3)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_Semi)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      END

      SUBROUTINE TabISR2ifi(iSel,Char,Mcapt)
*////////////////////////////////////////////////////////////////////////////////
*//               gmake TabISR-2d.eps
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER       iSel
      CHARACTER*1 Char
      CHARACTER*80  Mcapt
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      INTEGER       iThe,nCols
      CHARACTER*32  TabLab1(5), TabLab2(1), TabLab3(1)
      DATA TabLab1 / ' $d$ ', ' $u$ ', ' $s$ ', ' $c$ ', ' $b$ '/ !
      DATA TabLab2 / ' $all$' /
      DATA TabLab3 / ' $\\mu$' /
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='{\\color{blue} $f$}'
      capt(2)='{\\color{blue} $\\sigma^{IFI}/\\sigma$ }'
      capt(3)='{\\color{blue} $A^{IFI}_{FB}$  }' !
      idl(1) = ix_IFI
      idl(2) = ia_Int
      nCols  = 2
c      idl(1) = ix_Best
c      idl(2) = ia_Best
!----------------------------------------------------------
      iThe = 1  ! sigma
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  1,iSel, 1d3, ix_Best)  ! KKMC wtbest
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  7,iSel, 1d3, ix_IFI)   ! KKMC IFI
      CALL GLK_Operat(ix_IFI, '/', ix_Best,  ix_IFI,  1d0,1d0)   !
      iThe = 2  ! AFB
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  1,iSel, 1d0, ia_Best)   ! KKMC IFI
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  2,iSel, 1d0, ia_NoInt)  ! KKMC IFI
      CALL SAB_GetKFHist(m_kSABren,m_iKFf1,m_iKFf2, iThe,  6,iSel, 1d0, ia_IFI)    ! Not good!!!
      CALL GLK_Operat(ia_Best, '-', ia_NoInt,  ia_Int,  1d0,1d0)   !

      CALL GLK_SetTabRan(1,5,1)
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,Char,'R',' ')
      CALL GLK_SetTabRan(7,7,1)
      CALL GLK_SetTabLab(1,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_SetTabRan(13,13,1)
      CALL GLK_SetTabLab(1,TabLab3)
      CALL GLK_PlTable2(nCols,idl,capt,'   ',fmt,'S','R',' ')
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_Semi)
      CALL GLK_Delet(ia_ZFter)
      CALL GLK_Delet(ia_Diff1)
      CALL GLK_Delet(ia_Diff2)
      END
