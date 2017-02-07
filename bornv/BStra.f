*//////////////////////////////////////////////////////////////////////////////////////
*//                                                                                  //
*//          Pseudoclass BStra                                                       //
*//                                                                                  //
*//                                                                                  //
*//   Foam is now the basic MC sampler for beamstrahlung and ISR.                    //
*//   Notes: BStra has internal rejection procedure,                                 //
*//   consequently, normalization is determined by average weight,                   //
*//   which is monitored by MBrB                                                     //
*//   Brancher MBrB has only one branch, its purpose to monitoring weight of FoamC   //
*//////////////////////////////////////////////////////////////////////////////////////

      SUBROUTINE BStra_Initialize(KeyGrid,KeyWgt,Xcrude)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Initialization phase                                                           //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      INCLUDE 'BXformat.h'
      INTEGER            KeyGrid, KeyWgt
      DOUBLE PRECISION   XCrude
      INTEGER  k,j
      DOUBLE PRECISION   XsCru(10), WMList(10)
      DOUBLE PRECISION   XsectA, XsectB,  XsectC
      DOUBLE PRECISION   BornV_RhoFoamC
      EXTERNAL           BornV_RhoFoamC
      INTEGER            nCallsA,  nCallsB,  nCallsC
      INTEGER            IterMaxA, IterMaxB, IterMaxC,  Idyfs,  IdBra, Nbin
*-------------------------------
      m_out     = 16
      m_Nevgen  =  0
*//////////////////////////////////////////////////////////////////////////////////////
      m_KeyWgt =   KeyWgt
      m_Nev    = 0d0
      m_SWT    = 0d0
      m_SSWT   = 0d0
*//////////////////////////////////////////////////////////////////////////////////////
*//   3-dimensional case                                                             //
*//////////////////////////////////////////////////////////////////////////////////////
      WRITE(*,*) '*****************************************************************'
      WRITE(*,*) '****** BE PATIENT FoamC CREATING GRID FOR BEAMSTRAHLUNG *********'
      WRITE(*,*) '*****************************************************************'
      CALL FoamC_SetKdim(       3) ! No. of dimensions<5
      CALL FoamC_SetnBuf(   40000) ! Length of buffer<5000,  =Maximum No. of cells
      CALL FoamC_SetnSampl(100000) ! No. of MC sampling inside single cell, default=100
      CALL FoamC_SetnBin(       4) ! No of bins for edge explorations
      CALL FoamC_SetEvPerBin(  25) ! No. of equiv. MC events per bin
      CALL FoamC_SetOptEdge(    0) ! OptEdge excludes vertices
      CALL FoamC_SetOptDrive(   2) ! 0,1,2 =True,Sigma,WtMax !!!
      CALL FoamC_SetOptRanIni(  0) ! No internal initialization of rand.num.gen
      CALL FoamC_SetOptRanLux( -1) ! Ranmar choosen
      CALL FoamC_SetChat(       1) ! printout level =0,1,2
      CALL FoamC_Initialize(BornV_RhoFoamC)
*//////////////////////////////////////////////////////////////////////////////////////
      CALL FoamC_GetTotPrim( XsectC) ! Crude from Initialization
      m_XCrude =  XsectC ! is it really used ?????
*
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  BStra  Initializator                '
      WRITE(m_out,bxtxt) '  Foam initialization finished        '
      WRITE(m_out,bxl1g) XsectC ,    'XsectC  3-dimen.  ','XsectC','**'
      WRITE(m_out,bxl1g) m_XCrude,   'XCrude  total.    ','XCrude','**'
      WRITE(m_out,bxclo)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Initialization of MBrB, the own copy of a brancher                             //
*//   Now brancher has only one branch, its purpose to monitoring weight of FoamC
      CALL KK2f_GetIdyfs(Idyfs)
      IdBra = Idyfs+200
      CALL MBrB_Initialize(m_out,IdBra,50, 1d0, 'MBrB: Bstra main weight$')
      Nbin      = 500
      WMList(1) = 1.00d0
      XsCru(1)  = XsectC
      CALL MBrB_AddBranch(1, Nbin, WMList(1), 'MBrB: branch for FoamC !!$')
      CALL MBrB_SetXSList(XsCru)
      CALL MBrB_GetXCrude(m_XCrude)
*//////////////////////////////////////////////////////////////////////////////////////
*// Because in Bstra we have internal rejection loop we send to Karlud and KK2f      //
*// the best estimator of integral we have at this moment                            // 
*// It will be used for histogram normalization (sometimes)                          // 
*// Note that xsection from KK2f_finalize uses m_XCrude*<wt)> or m_XCrude            //
*//////////////////////////////////////////////////////////////////////////////////////
      XCrude   = m_XCrude
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  BStra  Initializator, PreGeneration '
      WRITE(m_out,bxl1g) XsCru(1) ,   'XsCru(1)  3-dimen.  ','XsCru(1) ','**'
      WRITE(m_out,bxl1f) WMlist(1) ,  'WMlist(1) 3-dimen.  ','WMlist(1)','**'
      WRITE(m_out,bxl1g) m_XCrude ,   'XCrude   total      ','XCrude   ','**'
      WRITE(m_out,bxclo)
      END

      SUBROUTINE BStra_Make(vv, x1, x2, MCwt)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Genearete set of 3 ISR variables for beamsstrahlung ISR                        //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      DOUBLE PRECISION  vv, x1, x2, MCwt, x, Wt_KF
      DOUBLE PRECISION  BornV_RhoFoamC
      EXTERNAL          BornV_RhoFoamC
      REAL              Qrand(10)        ! for PseuMar
      INTEGER           Itype
      DOUBLE PRECISION  rand
*-------------------------------------------------------------------------------
      m_NevCru = 0
 100  CONTINUE
      m_Nevgen  =  m_Nevgen +1
      CALL MBrB_GenKF(Itype, Wt_KF)  ! Needed to define KF_last
      Wt_KF = 1

      CALL FoamC_MakeEvent(BornV_RhoFoamC) ! generate MC event
      CALL FoamC_GetMCwt(  MCwt) ! get MC weight

      CALL BornV_GetVXX(vv,x1,x2)

* Rejection or wted events
*[[[[  m_KeYWgt=0 option does it work properly!!!!???
      IF(   m_KeYWgt .EQ. 0) THEN
        CALL PseuMar_MakeVec(Qrand,1)
        rand = Qrand(1)
        CALL MBrB_Fill(MCwt   ,rand)
        m_NevCru = m_NevCru+1
        IF(rand .GT. MCwt) GOTO 100
        MCwt = 1d0
      ELSE
        m_NevCru = 1
        MCwt = MCwt *Wt_KF
      ENDIF

      CALL MBrB_Fill(MCwt   ,0d0)
      END                       ! BStra_Make

      SUBROUTINE BStra_GetXCrude(XCrude)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Get TRUE crude integraml                                                       //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      DOUBLE PRECISION   XCrude
      XCrude   = m_XCrude
      END                       ! BStra_GetXCrude

      SUBROUTINE BStra_Finalize(Integ,Errel)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Provides Crude integral at the end of MC generation based on <wt>      //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      INCLUDE 'BXformat.h'
      DOUBLE PRECISION     Integ,Errel
      DOUBLE PRECISION     IntegMC,ErrelMC
      DOUBLE PRECISION     AverWt, WtSup
      DOUBLE PRECISION     AwtNew
*-----------------------------------------------------------------------------
      CALL MBrB_MgetAve(AverWt, ErRelMC, WtSup)
      IntegMC= m_XCrude*AverWt
      Integ  = IntegMC
      ErRel  = ErRelMC
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  BStra  Finalize MC results     '
      WRITE(m_out,bxl1g) IntegMC,   'MC integral   ','IntegMC','**'
      WRITE(m_out,bxl1f) ErRelMC,   'relat. error  ','ErRelMC','**'
      WRITE(m_out,bxl1f) AverWt,    'average wt    ','AverWt ','**'
      WRITE(m_out,bxl1f) WtSup,     'maximum wt    ','WtSup  ','**'
      WRITE(m_out,bxtxt) '  From grid building (initializ.)'  
      WRITE(m_out,bxl1g) m_XCrude,  'XCrude  total.','XCrude','**'
      WRITE(m_out,bxclo)
* Print more on the main weight
      CALL MBrB_Print0
* Print even more on the weight in each branch!
      CALL MBrB_Print1
      END       ! BStra_Finalize

      SUBROUTINE BStra_GetSumCru(XCrude,NevCru)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Get TRUE crude local sums                                                     //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      DOUBLE PRECISION     XCrude
      INTEGER              NevCru
      XCrude   = m_XCrude
      NevCru   = m_NevCru
      END                       ! BStra_GetXGridB

      SUBROUTINE BStra_GetIntegMC(IntegMC,ErRelMC)
*//////////////////////////////////////////////////////////////////////////////////////
*//   Get TRUE Monte Carlo run Integral and errors                                   //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      DOUBLE PRECISION     IntegMC,ErRelMC
      DOUBLE PRECISION     AverWt,WtSup
      CALL MBrB_MgetAve(AverWt, ErRelMC, WtSup)
      IntegMC= m_XCrude*AverWt
      END                       ! BStra_GetIntegMC

      SUBROUTINE BStra_GetAveWt(AveWt,RatWt)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BStra.h'
      DOUBLE PRECISION     AveWt,RatWt
      DOUBLE PRECISION     AverWt, ErRela, WtSup
*-----------------------------------------------------------------------------
      CALL MBrB_MgetAve(AverWt, ErRela, WtSup)
      AveWt = AverWt
      RatWt = AverWt/WtSup
      END
*//////////////////////////////////////////////////////////////////////////////////////
*//                                                                                  //
*//          END of Pseudoclass BStra                                                //
*//                                                                                  //
*//////////////////////////////////////////////////////////////////////////////////////
