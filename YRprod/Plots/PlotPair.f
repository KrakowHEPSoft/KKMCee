*------------------------------------------------------
*     gmake Sli-Pair.ps
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
*-------------------------------------------------------------------------------
c      Dname  = '../189GeV/PairsISRonly.input'
c      Hname  = '../189GeV/PairsISRonly.hst'
c      Dname  = '../189GeV/PairsFSR.input'
c      Hname  = '../189GeV/PairsFSR.hst' ! current
      Dname  = '../189GeV/PairsFSR.input'
      Hname  = '../189GeV/pro.hst' ! current
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
*-------------
      CALL SAB_MakeRen
* Semi-analytical
      CALL SAB_FilKKsem
* Data from ZFitter
c      CALL SAB_RDatZF1('./zf_2f_ForSJ_mu_Minv.log.B',m_kSABsem) !
      END

      SUBROUTINE TabIFI1
*////////////////////////////////////////////////////////////////////////////////                             
*     gmake TabPair-1.eps
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
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 6,  1,m_nVmax, 1d3, ix_Pair)  ! KKMC Pairs
      CALL SAB_GetHistSel(m_kSABsem,ikFf,iThe,10,  1,m_nVmax, 1d3, ix_Semi)  ! KKsem
      CALL GLK_Operat(ix_Pair, '/',  ix_Semi,  ix_Diff1,  1d0,1d0)   ! difference1
*-------------------------------------------
      TeXfile   = 'TabPair-1.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
*===================================================================
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='{\\color{blue}$v_{\\max}$}'
      capt(2)='{\\color{blue}(a) ${\\cal KK}$sem Refer.}'
      capt(3)='{\\color{blue}(b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm IFIon}$}'!
      capt(4)='{\\color{blue}(c) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm IFIoff}$}'!
      capt(5)='{\\color{blue}(c) Pair }'
      capt(6)='{\\color{blue}(c) Pair/tot }'
      idl(1) = ix_Semi
      idl(2) = ix_Best
      idl(3) = ix_NoInt
      idl(4) = ix_Pair
      idl(5) = ix_Diff1
      nCols  = 5
!----------------------------------------------------------
      Mcapt = '$\\sigma(q\\bar{q})$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_ZFNoInt)
      CALL GLK_Delet(ix_ZFInt)
      CALL GLK_Delet(ix_ZFIntExp)
      CALL GLK_Delet(ix_Diff1)
      CALL GLK_Delet(ix_Diff2)
      END
