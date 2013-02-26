*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*!!!!! for files marked (O) before 6-th June, MUST correct RobAll.h !!!!
*------------------------------------------------------
*     make YRtabQuark-ps
*     make YRtabMu-ps
*     make YRtabTau-ps
*------------------------------------------------------
      PROGRAM MAIN
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      CHARACTER*60  Dname
      CHARACTER*60  Hname, DumpFile

c      Dname  = '../189GeV/189GeV.input'
c      Hname  = '../189GeV/pro.hst' ! current
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*   for files marked (O) before 6-th June 2000, MUST correct RobAll.h
*   for nungam exercise there also might be corrections in   RobAll.h
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*-------------------------------------------------------------------------------
* New muon with realistic and photonic
c     make YRtabMu-ps
c       Dname  = '../206GeV/MuIFIon206GeVreal.input.12M'  !(N) XCHECK, Oct.2001
c       Hname  = '../206GeV/MuIFIon206GeVreal.hst.12M'    !(N) vmax =0.999
c       Dname  = '../206GeV/MuIFIon206GeVreal.input.11M'  !(N) XCHECK, Nov.2000
c       Hname  = '../206GeV/MuIFIon206GeVreal.hst.11M'    !(N) vmax =0.999
**      Dname  = '../189GeV/MuIFIon189GeVreal.input'    !    Mu for YR
**      Hname  = '../189GeV/MuIFIon189GeVreal.hst.15M'  !(N) vmax =0.999
**      Hname  = '../189GeV/MuIFIon189GeV.hst.17M'      !(O) vmax =0.99
**      Dname  = '../206GeV/MuIFIon206GeVreal.input'    !    Mu for YR
**      Hname  = '../206GeV/MuIFIon206GeVreal.hst.30M'  !(N) vmax =0.999
**      Hname  = '../206GeV/MuIFIon206GeV.hst.9M'       !(O) vmax =0.99
*-------------------------------------------------------------------------------
*     make YRtabQuark-ps
c      Dname  = '../189GeV/IncISRonly189GeV.input.140M'     !(N) Xcheck, Sept.2002
c      Hname  = '../189GeV/IncISRonly189GeV.hst.140M'       !(N) quarks+mu
c      Dname  = '../189GeV/IncISRonly189GeV.input.173M'     !(N) Xcheck, Oct.2001
c      Hname  = '../189GeV/IncISRonly189GeV.hst.173M'       !(N) quarks+mu
c      Dname  = '../189GeV/IncISRonly189GeV.input.157M'     !(N) Xcheck, May.2001
c      Hname  = '../189GeV/IncISRonly189GeV.hst.157M'       !(N) quarks+mu
c      Dname  = '../189GeV/IncISRonly189GeV.input.200M'     !(N) Xcheck, Nov.2000
c      Hname  = '../189GeV/IncISRonly189GeV.hst.200M'       !(N) quarks+mu
**      Dname  = '../189GeV/IncISRonly189GeV.input'     ! quarks+mu
**      Hname  = '../189GeV/IncISRonly189GeV.hst.217M'  !(O) quarks+mu
**      Dname  = '../200GeV/IncISRonly200GeV.input'     !    quarks+mu
**      Hname  = '../200GeV/IncISRonly200GeV.hst.151M'  !(O) quarks+mu
**      Dname  = '../206GeV/IncISRonly206GeV.input'     !    quarks+mu
**      Hname  = '../206GeV/IncISRonly206GeV.hst.239M'  !(O) quarks+mu
*-------------------------------------------------------------------------------
      Dname  = '../189GeV/MuIFIon189GeV.input'        !    Mu for YR
      Hname  = '../189GeV/MuIFIon189GeV.hst.17M'      !(O) Mu for YR
c      Dname  = '../200GeV/MuIFIon200GeV.input'        !    Mu for YR
c      Hname  = '../200GeV/MuIFIon200GeV.hst.12M'      !(O) Mu for YR
c      Dname  = '../206GeV/MuIFIon206GeV.input'        !    Mu for YR
c      Hname  = '../206GeV/MuIFIon206GeV.hst.9M'       !(O) Mu for YR
*-------------------------------------------------------------------------------
c      Dname  = '../189GeV/TauIFIon189GeV.input'       !    Tau for YR
c      Hname  = '../189GeV/TauIFIon189GeV.hst.27M'     !(O) Tau for YR
c      Dname  = '../200GeV/TauIFIon200GeV.input'       !    Tau for YR
c      Hname  = '../200GeV/TauIFIon200GeV.hst.26M'     !(O) Tau for YR
c      Dname  = '../206GeV/TauIFIon206GeV.input'       !    Tau for YR
c      Hname  = '../206GeV/TauIFIon206GeV.hst.16M'     !(O) Tau for YR
*-------------------------------------------------------------------------------
      m_out = 16
      CALL GLK_SetNout(m_out)

      OPEN(m_out,FILE='./PlotYR.output')
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,m_imax,m_xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(                 Dname, 0,m_imax,m_xpar)  ! reading actual user input

* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
      CALL GLK_ListPrint(m_out)     ! debug
*=========================================================================
      IF( Nint(m_xpar(401)).EQ.1)  THEN
         CALL YRtabQuark
      ELSE
         IF( Nint(m_xpar(413)).EQ.1)  CALL YRtabMuon
         IF( Nint(m_xpar(415)).EQ.1)  CALL YRtabTau
      ENDIF
*=========================================================================
* Write all histograms into dump file, for control
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=========================================================================
      CLOSE(m_out)
      END



      SUBROUTINE PrtLnSig(ind1,ind2,fLabel,CMSene,xObs,dObs,sysErr)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          i,ind1,ind2
      DOUBLE PRECISION CMSene, xObs(*),dObs(*),sysErr(*)
      CHARACTER*13   fLabel(*)
*
      DO i = ind1,ind2
         WRITE(30, m_GrandFormat)
     $        '=',    fLabel(i),     ! observable label position    2-14
     $        '=',    'KKMC 4.14  ', ! program name     position   16-29
     $        '=',    'S',           ! Observable type  S or A     31
     $        '=',    CMSene,        ! CMS energy       position   33-39
     $        '=',    xObs(i),       ! prediction       position   41-54
     $        '=',    dObs(i),       ! Statistical err. position   56-69
     $        '=',    sysErr(i),     ! Systematic error position   71-84
     $        '=',    'pb       ',      ! Prediction units position   86-94
     $        '=',    'pb       ',      ! Error      units position   96-104
     $        '=',    m_Comment,     ! Comment          position  106-131
     $        '='                    ! Separator        position  105
      ENDDO
      END

      SUBROUTINE PrtLnAFB(ind1,ind2,fLabel,CMSene,xObs,dObs,sysErr)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          i,ind1,ind2
      DOUBLE PRECISION CMSene, xObs(*),dObs(*),sysErr(*)
      CHARACTER*13   fLabel(*)
*
      DO i = ind1,ind2
         WRITE(30, m_GrandFormat)
     $        '=',    fLabel(i),     ! observable label position    2-14
     $        '=',    'KKMC 4.14  ', ! program name     position   16-29
     $        '=',    'A',           ! Observable type  S or A     31
     $        '=',    CMSene,        ! CMS energy       position   33-39
     $        '=',    xObs(i),       ! prediction       position   41-54
     $        '=',    dObs(i),       ! Statistical err. position   56-69
     $        '=',    sysErr(i),     ! Systematic error position   71-84
     $        '=',    '      ',      ! Prediction units position   86-94
     $        '=',    '      ',      ! Error      units position   96-104
     $        '=',    m_Comment,     ! Comment          position  106-131
     $        '='                    ! Separator        position  105
      ENDDO
      END

      SUBROUTINE PrtDatSig(idHis, n1, n2, fLabel, CMSene, PhysPrec)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          idHis, n1, n2
      DOUBLE PRECISION PhysPrec
      DOUBLE PRECISION CMSene, xObs(1000),dObs(1000),sysErr(1000)
      DOUBLE PRECISION GLK_hie, GLK_hi
      CHARACTER*13   fLabel(*)
      INTEGER i,k1,k2
*////////////////////////////////
      DO i= n1, n2
         xObs(i)  = GLK_hi( idHis,m_nYR+ i) !  Sigma
         dObs(i)  = GLK_hie(idHis,m_nYR+ i) ! dSigma statistical
         sysErr(i)= GLK_hi( idHis,m_nYR+ i)*PhysPrec ! pSigma physical
      ENDDO
      CALL PrtLnSig(n1, n2,fLabel,CMSene,xObs,dObs,sysErr)
      END


      SUBROUTINE PrtDatAFB(idHis, n1, n2, fLabel, CMSene, PhysPrec)
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER          idHis, n1, n2
      DOUBLE PRECISION PhysPrec
      DOUBLE PRECISION CMSene, xObs(1000),dObs(1000),sysErr(1000)
      DOUBLE PRECISION GLK_hie, GLK_hi
      CHARACTER*13   fLabel(*)
      INTEGER i,k1,k2
*////////////////////////////////
      DO i= n1, n2
         xObs(i)  = GLK_hi( idHis,m_nYR+ i) !  AFB
         dObs(i)  = GLK_hie(idHis,m_nYR+ i) ! dAFB statistical
         sysErr(i)= PhysPrec                ! dAFB physical
      ENDDO
      CALL PrtLnAFB(n1, n2,fLabel,CMSene,xObs,dObs,sysErr)
      END


      SUBROUTINE YRtabMuon
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  TeXfile
      INTEGER       Lint, nCols
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*32  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
      CHARACTER*6   Energy
      INTEGER       i,j,k,ikFf,iThe,iMod,nSel,iSel1,iSel2
      DOUBLE PRECISION CMSene,PhysPrec
      DOUBLE PRECISION xSel(m_iSel1,m_iSel2),eSel(m_iSel1,m_iSel2)
      CHARACTER*32 TabLab1(21)
      DATA TabLab1 /  ' $ v<0.01$ ', ' $ v<0.05$ ',
     $ ' $ v<0.10$ ', ' $ v<0.15$ ', ' $ v<0.20$ ', ' $ v<0.25$ ',
     $ ' $ v<0.30$ ', ' $ v<0.35$ ', ' $ v<0.40$ ', ' $ v<0.45$ ',
     $ ' $ v<0.50$ ', ' $ v<0.55$ ', ' $ v<0.60$ ', ' $ v<0.65$ ',
     $ ' $ v<0.70$ ', ' $ v<0.75$ ', ' $ v<0.80$ ', ' $ v<0.85$ ',
     $ ' $ v<0.90$ ', ' $ v<0.95$ ', ' $ v<0.99$ ' /
      CHARACTER*32 TabLab2i(11) ! Idealized
      DATA TabLab2i /
     $ 'IAleph5','IAleph6','IDelphi5','IDelphi6','ILT9','ILT10','ILT11','IOpal6','IOpal7','IOpal8','IOpal9'/!
      CHARACTER*32 TabLab2r(8) ! Realistic
      DATA TabLab2r /
     $ 'Aleph5', 'Aleph6', 'Delphi4','Delphi5','LT9','LT10','Opal6','Opal7' /!
      CHARACTER*32 TabLab2p(10) ! Photonic
      DATA TabLab2p /
     $ 'Aleph12', 'Aleph15', 'Delphi9','Delphi12','LT15','LT18','Opal14','Opal15','Opal20','Opal21'/!
*-------------------------------------------------------------------------------
      WRITE(*,*) '======================YRtabMuon Entered ===================='
      CMSene = m_xpar(1)
      IF( ABS(CMSene-189d0).LT.001) Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) Energy = '206GeV'
*-------------------------------------------------------------------------------
      CALL SAB_MakeRen
      ikFf =13
      iThe = 1  ! xsections
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d3, ix_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d3, ix_NoInt) !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 4,  m_iSel1,m_iSel2, 1d3, ix_EEX2) !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 5,  m_iSel1,m_iSel2, 1d3, ix_EEX3) !
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,   1d0,1d0)   ! IFI
      iThe = 2 ! asymetries
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d0, ia_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ia_NoInt) !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 4,  m_iSel1,m_iSel2, 1d0, ia_EEX2) !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 5,  m_iSel1,m_iSel2, 1d0, ia_EEX3) !
      CALL GLK_Operat(ia_Best, '-', ia_NoInt, ia_IFI,   1d0,1d0)   ! IFI
*-------------------------------------------
      TeXfile   = 'YRtabMu.tex'
cc      TeXfile   = 'YRtabMu'//Energy//'.tex'
      Lint =0
      CALL GLK_PlInitialize(Lint,TeXfile) !Initialize GLK_Plot
cc      CALL GLK_PlCapt(CapTab)
*===================================================================
*                Tables sigma
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='Select.'
      capt(2)='IFI on'
      capt(3)='IFI off'
      capt(4)='EEX3'
      capt(5)='EEX2'
      idl(1) = ix_Best
      idl(2) = ix_NoInt
      idl(3) = ix_EEX3
      idl(4) = ix_EEX2
      nCols  = 4
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\mu^-\\mu^+)$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\mu^-\\mu^+)$, YR Idealized, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+11 ,1)
      CALL GLK_SetTabLab(11,TabLab2i)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\mu^-\\mu^+)$, YR Realistic, at '//Energy
      CALL GLK_SetTabRan(m_nYR+12, m_nYR+19 ,1)
      CALL GLK_SetTabLab(8,TabLab2r)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      fmt(2)='F10.5'
      fmt(3)= 'F8.5'
      Mcapt = '$\\sigma(\\mu^-\\mu^+\\gamma)$, $\\sigma(\\mu^-\\mu^+2\\gamma)$, '//Energy !
      CALL GLK_SetTabRan(m_nYR+20, m_nYR+29 ,1)
      CALL GLK_SetTabLab(10,TabLab2p)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
*===================================================================
*                Tables AFB
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      idl(1) = ia_Best
      idl(2) = ia_NoInt
      idl(3) = ia_EEX3
      idl(4) = ia_EEX2
      nCols  = 4
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\mu^-\\mu^+)$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\mu^-\\mu^+)$, YR Idealized, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+11, 1)
      CALL GLK_SetTabLab(11,TabLab2i)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\mu^-\\mu^+)$, YR Realistic, at '//Energy
      CALL GLK_SetTabRan(m_nYR+12, m_nYR+19, 1)
      CALL GLK_SetTabLab(8,TabLab2r)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
*=================================================================================
*     WRITING Observable data card for Sigma nad AFB
*=================================================================================
*=================================================================================
*                        XSECTIONS
*=================================================================================
      OPEN(30,File='MuonPred'//Energy//'.data')
* Idealized: IAleph5,IAleph6,IDelphi5,IDelphi6,ILT9,ILT10,ILT11
      PhysPrec = 0.02d0
      m_Comment = ' Minv IFIon   '
      CALL PrtDatSig(ix_Best,  1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv IFIoff  '
      CALL PrtDatSig(ix_NoInt, 1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv EEX3    '
      CALL PrtDatSig(ix_EEX3,  1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv EEX2    '
      CALL PrtDatSig(ix_EEX2,  1, 7,m_MuLabel,CMSene,PhysPrec)
* Idealized: IOpal6','IOpal7','IOpal8','IOpal9'
      PhysPrec = 0.02d0
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtDatSig(ix_NoInt, 8,11, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Mprop-, EEX3  '
      CALL PrtDatSig(ix_EEX3 , 8,11, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Mprop-, EEX2  '
      CALL PrtDatSig(ix_EEX2 , 8,11, m_MuLabel,CMSene,PhysPrec)
* Realistic Aleph5, Aleph6,Delphi4,Delphi5,LT9,LT10,Opal6,Opal7
      PhysPrec = 0.02d0
      m_Comment = ' CEEX2         '
      CALL PrtDatSig(ix_Best,  12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' CEEX2   IFIoff'
      CALL PrtDatSig(ix_NoInt, 12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX2          '
      CALL PrtDatSig(ix_EEX2 , 12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX3          '
      CALL PrtDatSig(ix_EEX3 , 12,19, m_MuLabel,CMSene,PhysPrec)
* Photonic  Aleph12, Aleph15, Delphi9, Delphi12, LT15, LT18, Opal14, Opal15, Opal20, Opal21
      PhysPrec = 0.00d0
      m_Comment = ' CEEX2         '
      CALL PrtDatSig(ix_Best,  20,29, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' CEEX2   IFIoff'
      CALL PrtDatSig(ix_NoInt, 20,29, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX2          '
      CALL PrtDatSig(ix_EEX2 , 20,29, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX3          '
      CALL PrtDatSig(ix_EEX3 , 20,29, m_MuLabel,CMSene,PhysPrec)
*=================================================================================
*                        AFB
*=================================================================================
* Idealized: IAleph5,IAleph6,IDelphi5,IDelphi6,ILT9,ILT10,ILT11
      PhysPrec = 0.02d0
      m_Comment = ' Minv, IFIon   '
      CALL PrtDatAFB(ia_Best, 1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv, IFIoff  '
      CALL PrtDatAFB(ia_NoInt, 1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv, EEX3    '
      CALL PrtDatAFB(ia_EEX3, 1, 7,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Minv, EEX2    '
      CALL PrtDatAFB(ia_EEX2, 1, 7,m_MuLabel,CMSene,PhysPrec)
* Idealized: IOpal6','IOpal7','IOpal8','IOpal9'
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtDatAFB(ia_NoInt, 8, 11,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Mprop-, EEX3  '
      CALL PrtDatAFB(ia_EEX3, 8, 11,m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' Mprop-, EEX2  '
      CALL PrtDatAFB(ia_EEX2, 8, 11,m_MuLabel,CMSene,PhysPrec)
* Realistic Aleph5, Aleph6,Delphi4,Delphi5,LT9,LT10,Opal6,Opal7
      PhysPrec = 0.02d0
      m_Comment = ' CEEX2         '
      CALL PrtDatAFB(ia_NoInt, 12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' CEEX2   IFIoff'
      CALL PrtDatAFB(ia_NoInt, 12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX2          '
      CALL PrtDatAFB(ia_EEX2 , 12,19, m_MuLabel,CMSene,PhysPrec)
      m_Comment = ' EEX3          '
      CALL PrtDatAFB(ia_EEX3 , 12,19, m_MuLabel,CMSene,PhysPrec)
      CLOSE(30)
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ix_EEX2)
      CALL GLK_Delet(ix_EEX3)
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      CALL GLK_Delet(ia_EEX2)
      CALL GLK_Delet(ia_EEX3)
      WRITE(*,*) '======================YRtabMuon Ended ===================='
      END


      SUBROUTINE YRtabTau
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  TeXfile
      INTEGER       Lint, nCols
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*32  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
      CHARACTER*6   Energy
      INTEGER       i,j,k,ikFf,iThe,iMod,nSel,iSel1,iSel2
      DOUBLE PRECISION CMSene, xObs(50),dObs(50),sysErr(50),PhysPrec
      DOUBLE PRECISION GLK_hie, GLK_hi
      DOUBLE PRECISION xSel(m_iSel1,m_iSel2),eSel(m_iSel1,m_iSel2)
      CHARACTER*32 TabLab1(21), TabLab2(11)
      DATA TabLab1 /  ' $ v<0.01$ ', ' $ v<0.05$ ',
     $ ' $ v<0.10$ ', ' $ v<0.15$ ', ' $ v<0.20$ ', ' $ v<0.25$ ',
     $ ' $ v<0.30$ ', ' $ v<0.35$ ', ' $ v<0.40$ ', ' $ v<0.45$ ',
     $ ' $ v<0.50$ ', ' $ v<0.55$ ', ' $ v<0.60$ ', ' $ v<0.65$ ',
     $ ' $ v<0.70$ ', ' $ v<0.75$ ', ' $ v<0.80$ ', ' $ v<0.85$ ',
     $ ' $ v<0.90$ ', ' $ v<0.95$ ', ' $ v<0.99$ ' /
      DATA TabLab2 /
     $ 'IAleph7','IAleph8','IDelphi7','IDelphi8','ILT12','ILT13','ILT14','IOpal10','IOpal11','IOpal12','IOpal13' / !
*-------------------------------------------------------------------------------
      WRITE(*,*) '======================YRtabTau Eneterd ===================='
      CMSene = m_xpar(1)
      IF( ABS(CMSene-189d0).LT.001) Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) Energy = '206GeV'
*-------------------------------------------------------------------------------
      CALL SAB_MakeRen
      ikFf =15
      iThe = 1  ! xsections
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d3, ix_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d3, ix_NoInt) !
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,   1d0,1d0)   ! IFI
      iThe = 2 ! asymetries
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d0, ia_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ia_NoInt) !
      CALL GLK_Operat(ia_Best, '-', ia_NoInt, ia_IFI,   1d0,1d0)   ! IFI
*-------------------------------------------
      TeXfile   = 'YRtabTau.tex'
cc      TeXfile   = 'YRtabTau'//Energy//'.tex'
      Lint =0
      CALL GLK_PlInitialize(Lint,TeXfile) !Initialize GLK_Plot
cc      CALL GLK_PlCapt(CapTab)
*===================================================================
*                Tables
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='Select.'
      capt(2)='IFI on'
      capt(3)='IFI off'
      capt(4)='IFI alone'
      idl(1) = ix_Best
      idl(2) = ix_NoInt
      idl(3) = ix_IFI
      nCols  = 3
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\tau^-\\tau^+)$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$\\sigma(\\tau^-\\tau^+)$, YR selections, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+11 ,1)
      CALL GLK_SetTabLab(11,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      idl(1) = ia_Best
      idl(2) = ia_NoInt
      idl(3) = ia_IFI
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\tau^-\\tau^+)$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\tau^-\\tau^+)$, YR selections, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+11, 1)
      CALL GLK_SetTabLab(11,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
*  WRITING Observable data card fro Sigma nad AFB
*                        XSECTIONS
      OPEN(30,File='TauPred'//Energy//'.data')
      DO i=1,7
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ix_Best,m_nYR+ i)            !  Sigma
         dObs(i)  = GLK_hie(ix_Best,m_nYR+ i)            ! dSigma statistical
         sysErr(i)= GLK_hi( ix_Best,m_nYR+ i)*PhysPrec   ! pSigma physical
      ENDDO
      m_Comment = ' Minv IFI on  '
      CALL PrtLnSig(1, 7,m_TauLabel,CMSene,xObs,dObs,sysErr)
      DO i=8,11
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ix_NoInt,m_nYR+ i)            !  Sigma
         dObs(i)  = GLK_hie(ix_NoInt,m_nYR+ i)            ! dSigma statistical
         sysErr(i)= GLK_hi( ix_NoInt,m_nYR+ i)*PhysPrec   ! pSigma physical
      ENDDO
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtLnSig(8,11,m_TauLabel,CMSene,xObs,dObs,sysErr)
*=================================================================================
*                        AFB
      DO i=1,7
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ia_Best,m_nYR+ i)   !  AFB
         dObs(i)  = GLK_hie(ia_Best,m_nYR+ i)   ! dAFB statistical
         sysErr(i)= PhysPrec                   ! pAFB physical
      ENDDO
      m_Comment = ' Minv, IFI on  '
      CALL PrtLnAFB(1, 7,m_TauLabel,CMSene,xObs,dObs,sysErr)
      DO i=8,11
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ia_NoInt,m_nYR+ i)   !  AFB
         dObs(i)  = GLK_hie(ia_NoInt,m_nYR+ i)   ! dAFB statistical
         sysErr(i)= PhysPrec                 ! pAFB physical
      ENDDO
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtLnAFB(8,11,m_TauLabel,CMSene,xObs,dObs,sysErr)
      CLOSE(30)
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      WRITE(*,*) '======================YRtabTau Ended ===================='
      END


      SUBROUTINE YRtabQuark
*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*60  TeXfile
      INTEGER       Lint, nCols
*------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*32  capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
      CHARACTER*6   Energy
      INTEGER       i,j,k,ikFf,iThe,iMod,nSel,iSel1,iSel2
      DOUBLE PRECISION CMSene, xObs(50),dObs(50),sysErr(50),PhysPrec
      DOUBLE PRECISION GLK_hie, GLK_hi
      DOUBLE PRECISION xSel(m_iSel1,m_iSel2),eSel(m_iSel1,m_iSel2)
      CHARACTER*32 TabLab1(21), TabLab2(9)
      DATA TabLab1 /  ' $ v<0.01$ ', ' $ v<0.05$ ',
     $ ' $ v<0.10$ ', ' $ v<0.15$ ', ' $ v<0.20$ ', ' $ v<0.25$ ',
     $ ' $ v<0.30$ ', ' $ v<0.35$ ', ' $ v<0.40$ ', ' $ v<0.45$ ',
     $ ' $ v<0.50$ ', ' $ v<0.55$ ', ' $ v<0.60$ ', ' $ v<0.65$ ',
     $ ' $ v<0.70$ ', ' $ v<0.75$ ', ' $ v<0.80$ ', ' $ v<0.85$ ',
     $ ' $ v<0.90$ ', ' $ v<0.95$ ', ' $ v<0.99$ ' /
      DATA TabLab2 /
     $ 'IAleph1$\\dag$','IAleph2$\\dag$','IDelphi1$\\dag$','IDelphi2$\\dag$','ILT1','ILT2','ILT3','IOpal1','IOpal2' / !
*-------------------------------------------------------------------------------
      WRITE(*,*) '======================YRtabQuark Eneterd ===================='
      CMSene = m_xpar(1)
      IF( ABS(CMSene-189d0).LT.001) Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) Energy = '206GeV'
*-------------------------------------------------------------------------------
      CALL SAB_MakeRen
      ikFf = 1  ! d-quark
      ikFf = 2  ! u-quark
      ikFf = 3  ! s-quark
      ikFf = 4  ! c-quark
      ikFf = 5  ! b-quark
      ikFf = 7  ! the sum over all quarks
      iThe = 1  ! xsections
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d3, ix_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d3, ix_NoInt) !
      CALL GLK_Operat(ix_Best, '-', ix_NoInt, ix_IFI,   1d0,1d0)   ! IFI
      iThe = 2 ! asymetries
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 1,  m_iSel1,m_iSel2, 1d0, ia_Best)  !
      CALL SAB_GetHistSel(m_kSABren,ikFf,iThe, 2,  m_iSel1,m_iSel2, 1d0, ia_NoInt) !
      CALL GLK_Operat(ia_Best, '-', ia_NoInt, ia_IFI,   1d0,1d0)   ! IFI
*-------------------------------------------
      TeXfile   = 'YRtabQuark.tex'
cc      TeXfile   = 'YRtabQuark'//Energy//'.tex'
      Lint =0
      CALL GLK_PlInitialize(Lint,TeXfile) !Initialize GLK_Plot
cc      CALL GLK_PlCapt(CapTab)
*===================================================================
*                Tables
!----------------------------------------------------------
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.3'
      capt(1)='Select.'
      capt(2)='IFI on'
      capt(3)='IFI off'
      capt(4)='IFI alone'
      idl(1) = ix_Best
      idl(2) = ix_NoInt
      idl(3) = ix_IFI
      nCols  = 3
!----------------------------------------------------------
      Mcapt = '$\\sigma(q\\bar{q})$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$\\sigma(q\\bar{q})$, YR with Mprop, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+9 ,1)
      CALL GLK_SetTabLab(9,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      idl(1) = ia_Best
      idl(2) = ia_NoInt
      idl(3) = ia_IFI
!----------------------------------------------------------
      Mcapt = '$A_{FB}(q\\bar{q})$, PRIMITIVE, at '//Energy
      CALL GLK_SetTabRan(1,21,1)
      CALL GLK_SetTabLab(21,TabLab1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
!----------------------------------------------------------
      Mcapt = '$A_{FB}(\\tau^-\\tau^+)$, YR with Mprop, at '//Energy
      CALL GLK_SetTabRan(m_nYR+1, m_nYR+9, 1)
      CALL GLK_SetTabLab(9,TabLab2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
!----------------------------------------------------------
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*=================================================================================
*  WRITING Observable data card fro Sigma nad AFB
*                        XSECTIONS
      OPEN(30,File='QuarkPred'//Energy//'.data')
      CMSene = m_xpar(1)
      DO i=1,9
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ix_NoInt,m_nYR+ i)            !  Sigma
         dObs(i)  = GLK_hie(ix_NoInt,m_nYR+ i)            ! dSigma statistical
         sysErr(i)= GLK_hi( ix_NoInt,m_nYR+ i)*PhysPrec   ! pSigma physical
      ENDDO
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtLnSig(5,9,m_QuLabel,CMSene,xObs,dObs,sysErr)
*=================================================================================
*                        AFB
      CMSene = m_xpar(1)
      DO i=1,9
         PhysPrec = 0.02d0
         xObs(i)  = GLK_hi( ia_NoInt,m_nYR+ i)   !  AFB
         dObs(i)  = GLK_hie(ia_NoInt,m_nYR+ i)   ! dAFB statistical
         sysErr(i)= PhysPrec                 ! pAFB physical
      ENDDO
      m_Comment = ' Mprop-, IFIoff'
      CALL PrtLnAFB(5,9,m_QuLabel,CMSene,xObs,dObs,sysErr)
      CLOSE(30)
* cleaning
      CALL GLK_Delet(ix_Best)
      CALL GLK_Delet(ix_NoInt)
      CALL GLK_Delet(ia_Best)
      CALL GLK_Delet(ia_NoInt)
      WRITE(*,*) '======================YRtabQuark Ended ===================='
      END

