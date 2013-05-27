*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb-sig2-ps
*/////////////////////////////////////////////////////////////////////////////////
*----------------------------------------------------------------
      PROGRAM MAIN
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INCLUDE 'afb.h'
*
      CHARACTER*60  Tesnam, Dname, Hname
      CHARACTER*60  Hname2, DumpFile
      REAL*8        sigma55
*--------------------------------------------------------------------
      ninp=  5
      nout= 16

      Tesnam    = 'afb_int2'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

*=====================================
* Read data, the same as in MC run
* Read histograms from MC and prepare
      Hname2  = './afb.hst'
      CALL GLK_ReadFile(Hname2)
*=====================================================
* Read names of production input files
      OPEN( 23, file='afb_prepare.files')
      READ(23,'(a)') Dname
      READ(23,'(a)') Hname
      CLOSE(23)
c      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
c      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
c      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
****  CALL KK2f_GetCMSene( m_CMSene)                    ! get CMS energy
      CALL KKsem_Initialize(xpar)                     ! Initialize semianalytical package
      CALL KKsem_GetCMSene(m_CMSene) ! get CMS energy 
*=====================================================
      CALL GLK_GetBin(isigO3 ,55, sigma55)
      m_iscale = 0
      IF( sigma55 .GT. 0.2d0) m_iscale = 1 !!! {{{{{ big scale }}}}} close to MZ
      IF( sigma55 .GT. 2.0d0) m_iscale = 2 !!! {{{{{ BIG scale }}}}} hadronic close to MZ
      WRITE(*,*) 'sigma55,m_iscale  = ',sigma55,m_iscale
*=====================================================
      CALL Plot_tab1
      CALL Plot_Gsig
      CALL Plot_Gafb
*
      CALL Plot_afbHO
      CALL Plot_sigHO
*=========================================================================
* Write all histograms into dump file, for control/debug
*=================================
      CLOSE(nout)
      END



      SUBROUTINE Plot_tab1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int2-tab1.eps
*//   Table of everything
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
! Parameters for tables
      INTEGER       idl(7),nColumn
      CHARACTER*80  capt(8)
      CHARACTER*16  fmt(3)
      CHARACTER*80  mcapt
*===================================================================
*                Writing table in LaTex format
      TeXfile   = 'afb_int2-tab1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      capt(1)='{\\color{blue}$v_{\\max}$}'
      capt(2)='{\\color{blue} ${\\cal KK}$sem Refer.}'
      capt(3)='{\\color{blue}${\\cal O}(\\alpha^3)_{\\rm EEX3}$ }'
      capt(4)='{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ intOFF}'
      capt(5)='{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }'
      capt(6)='{\\color{blue} KORALZ }' ! LEP2
      capt(7)='{\\color{blue} KORALZ Interf.}'
      fmt(1)='F10.2'
      fmt(2)='F10.4'
      fmt(3)='F8.4'
      CALL GLK_Operat(isigO3,    '+', isigO3,     isigO3,    1d3, 0d0) ! [nb]-->[pb]
      CALL GLK_Operat(isigG1,    '+', isigG1,     isigG1,    1d3, 0d0) ! [nb]-->[pb]
      CALL GLK_Operat(isigG1nin, '+', isigG1nin,  isigG1nin, 1d3, 0d0) ! [nb]-->[pb]
      CALL GLK_Operat(isigG2,    '+', isigG2,     isigG2,    1d3, 0d0) ! [nb]-->[pb]
      CALL GLK_Operat(isigG2nin, '+', isigG2nin,  isigG2nin, 1d3, 0d0) ! [nb]-->[pb]
!------------------------------------------------
      idl(1)= iO3best
      idl(2)= isigO3
      idl(3)= isigG2nin
      idl(4)= isigG2
      idl(5)= iMustSig
      idl(6)= iMustSigInt
      nColumn=6
c\\\\\\\\\\\
      nColumn=4   ! KORALZ eliminated
c      nColumn=3 !!!!!!!!!!!!!
c\\\\\\\\\\\
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], ${\\cal KK}$ M.C. and KORALZ 1-st order}'!
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb]}'!
      CALL GLK_SetTabRan(1,1,1)
      CALL GLK_PlTable2(nColumn,idl,capt,Mcapt,fmt,' ','R',' ')

      CALL GLK_SetTabRan(10,90,20)
cc      CALL GLK_SetTabRan(10,90,10) !!!!!!!!!!!!!
      CALL GLK_PlTable2(nColumn,idl,capt,' ',fmt,'S','R',' ')

      CALL GLK_SetTabRan(99,99,1)
      CALL GLK_PlTable2(nColumn,idl,capt,' ',fmt,'S','R',' ')
*------------------------------------------------
      idl(1)= iO3bestAFB
      idl(2)= iafbO3
      idl(3)= iafbG2nin
      idl(4)= iafbG2
      idl(5)= iMustAfb
      idl(6)= iMustAfbInt
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}$A_{\\rm FB}(v_{\\max})$, ${\\cal KK}$ M.C. and KORALZ 1-st order}'!
      Mcapt ='{\\color{red}$A_{\\rm FB}(v_{\\max})$}'!
      CALL GLK_SetTabRan(1,1,1)
      CALL GLK_PlTable2(nColumn,idl,capt,Mcapt,fmt,'S','R',' ')

      CALL GLK_SetTabRan(10,90,20)
cc      CALL GLK_SetTabRan(10,90,10) !!!!!!!!!!!!!
      CALL GLK_PlTable2(nColumn,idl,capt,' ',fmt,'S','R',' ')

      CALL GLK_SetTabRan(99,99,1)
      CALL GLK_PlTable2(nColumn,idl,capt,' ',fmt,'S','R',' ')

      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_Gsig
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int2-Gsig.eps
*//
*//   Plot of SigTot(vMax) versus reference EEX3(vMax)
*//   (SigTot-SigRef)/SigRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER                iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(16)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
*  upper close to Z
c     $'\\PaveLr{ 0}{1100}{\\Huge ',
c     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',
*  lower
c     $'\\PaveLr{ 0}{ 900}{\\Huge ',
c     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',
*  inside
     $'\\PaveL{150}{1000}{\\Huge ',
     $'     ${\\sigma-\\sigma_{_{\\rm ref}} \\over \\sigma_{_{\\rm ref}}}$ }',!
*
     $'\\PaveLb{ 600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
*
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{550}{1060}{\\color{red} ',
     $' $\\circle{20}\\;\\circle{20}\\;\\circle{20}\\;\\circle{20}\\;\\circle{20}$}',!
     $'\\PaveL{700}{1060}{\\large\\color{red} CEEX2, Int.ON}',
     $'\\PaveL{500}{1000}{\\color{blue}\\line(1,0){150}}',
     $'\\PaveL{700}{1000}{\\large\\color{blue} CEEX2, Int.OFF }',
*
     $'\\PaveL{550}{ 940}{\\color{red} ',
     $'         $\\circle*{10}\\;\\circle*{10}\\;\\circle*{10}\\;\\circle*{10}$}',!
     $'\\PaveL{700}{ 940}{\\large\\color{black} EEX3, NO Int.}',
*
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
ccc      iRef = isigO3
      iRef =iO3best
      fmtx='f10.2'
      fmty='f10.3'

      CALL GLK_Operat(isigO3,     '-', iRef,     isigO3,       1d0, 1d0)
      CALL GLK_Operat(isigO3,     '/', iRef,     isigO3,       1d0, 1d0)

      CALL GLK_Operat(isigG1nin,  '-', iRef,     isigG1nin,    1d0, 1d0)
      CALL GLK_Operat(isigG1nin,  '/', iRef,     isigG1nin,    1d0, 1d0)

      CALL GLK_Operat(isigG2nin,  '-', iRef,     isigG2nin,    1d0, 1d0)
      CALL GLK_Operat(isigG2nin,  '/', iRef,     isigG2nin,    1d0, 1d0)
      IF( m_CMSene .LE. 200d0) THEN
         CALL GLK_Operat(isigComb,    '-', iRef,       isigComb,     1d0, 1d0) !
         CALL GLK_Operat(isigComb,    '/', iRef,       isigComb,     1d0, 1d0) !
         CALL GLK_Operat(iHyb1Sig,    '-', iRef,       iHyb1Sig,     1d0, 1d0) !
         CALL GLK_Operat(iHyb1Sig,    '/', iRef,       iHyb1Sig,     1d0, 1d0) !
         CALL GLK_Operat(iHyb2Sig,    '-', iRef,       iHyb2Sig,     1d0, 1d0) !
         CALL GLK_Operat(iHyb2Sig,    '/', iRef,       iHyb2Sig,     1d0, 1d0) !
      ENDIF

      ymin= -0.050d0
      ymax=  0.075d0
      CALL GLK_SetYminYmax(isigG1nin,ymin, ymax) !!!! LEP2, leptons
      CALL GLK_SetYminYmax(isigG2nin,ymin, ymax) !!!! LEP2, leptons
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(isigG1nin,-0.004d0, 0.004d0) !!!! {{{{ 91GeV leptons
         CALL GLK_SetYminYmax(isigG2nin,-0.004d0, 0.004d0) !!!! {{{{ 91GeV leptons
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(isigG1nin,-0.009d0, 0.009d0) !!!! {{{{ 91GeV hadrons
         CALL GLK_SetYminYmax(isigG2nin,-0.009d0, 0.009d0) !!!! {{{{ 91GeV hadrons
      ENDIF

      CALL GLK_Operat(isigG1,    '-',  iRef,     isigG1,       1d0, 1d0)
      CALL GLK_Operat(isigG1,    '/',  iRef,     isigG1,       1d0, 1d0)

      CALL GLK_Operat(isigG2,    '-',  iRef,     isigG2,       1d0, 1d0)
      CALL GLK_Operat(isigG2,    '/',  iRef,     isigG2,       1d0, 1d0)

      CALL GLK_Operat(iKorzSig,   '-', iRef,     iKorzSig,      1d0, 1d0)
      CALL GLK_Operat(iKorzSig,   '/', iRef,     iKorzSig,      1d0, 1d0)
      CALL GLK_idopt( iKorzSig,'ERRO')
*--------------------------------------------------------------------------
      TeXfile   = 'afb_int2-Gsig.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
***      CALL GLK_plot2(  isigG1nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_plot2(  isigG2nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
***      CALL GLK_plot2(  isigG1,    'S','*',ring          ,fmtx,fmty)
      CALL GLK_plot2(  isigG2,    'S','*',ring          ,fmtx,fmty)

      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  isigO3,   'S','*',dot        ,fmtx,fmty)

      IF( m_CMSene .LE. 200d0) THEN
      ENDIF
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
*--------------------------------------------------------------------------
      END



      SUBROUTINE Plot_Gafb
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int2-Gafb.eps
*//
*//   Plot of Afb(vMax) versus reference AfbRef of EEX3
*//   (Afb-AfbRef)/AfbRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER  iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(17)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
* outside
*     $'\\PaveLr{-15}{900}{ $A_{_{\\rm FB}}\\!\\!-\\!A^{\\rm ref}_{_{\\rm FB}}$}',!
* inside
     $'\\PaveL{150}{1050}{ $A_{_{\\rm FB}}\\!\\!-\\!A^{\\rm ref}_{_{\\rm FB}}$}',!
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',!
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{550}{1060}{\\color{red} ',
     $' $\\circle{20}\\;\\circle{20}\\;\\circle{20}\\;\\circle{20}\\;\\circle{20}$}',!
     $'\\PaveL{700}{1060}{\\large\\color{red} CEEX2, Int.ON}',
     $'\\PaveL{500}{1000}{\\color{blue}\\line(1,0){150}}',
     $'\\PaveL{700}{1000}{\\large\\color{blue}     CEEX2, Int.OFF}',
*
     $'\\PaveL{550}{ 940}{\\color{red} ',
     $'         $\\circle*{10}\\;\\circle*{10}\\;\\circle*{10}\\;\\circle*{10}$}',!
     $'\\PaveL{700}{ 940}{\\large\\color{black} EEX3, NO Int.}',
*
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      iRef = iO3bestAFB
**      iRef = iafbO3

      CALL GLK_Operat(iafbG1nin,  '-', iRef,     iafbG1nin,    1d0, 1d0)
      CALL GLK_Operat(iafbG2nin,  '-', iRef,     iafbG2nin,    1d0, 1d0)
      CALL GLK_Operat(iafbO3,     '-', iRef,     iafbO3,       1d0, 1d0)
      IF( m_CMSene .LE. 200d0) THEN
      CALL GLK_Operat(iafbComb,   '-', iRef,     iafbComb,     1d0, 1d0)
      CALL GLK_Operat(iHyb1Afb,   '-', iRef,     iHyb1Afb,     1d0, 1d0)
      CALL GLK_Operat(iHyb2Afb,   '-', iRef,     iHyb2Afb,     1d0, 1d0)
      ENDIF

      ymin= -0.035d0
      ymax=  0.050d0
c      ymin= -0.015d0
c      ymax=  0.150d0
      CALL GLK_SetYminYmax(iafbG1nin,ymin,ymax)
      CALL GLK_SetYminYmax(iafbG2nin,ymin,ymax)
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(iafbG1nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
         CALL GLK_SetYminYmax(iafbG2nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(iafbG1nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
         CALL GLK_SetYminYmax(iafbG2nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
      ENDIF

      CALL GLK_Operat(iafbG1,    '-',  iRef,     iafbG1,       1d0, 1d0)
      CALL GLK_Operat(iafbG2,    '-',  iRef,     iafbG2,       1d0, 1d0)

      CALL GLK_Operat(iKorzAfb,   '-', iRef,     iKorzAfb,      1d0, 1d0)
      CALL GLK_idopt( iKorzAfb,'ERRO')
* Plotting starts here
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int2-Gafb.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
cc      CALL GLK_plot2(  iafbG1nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_plot2(  iafbG2nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
cc      CALL GLK_plot2(  iafbG1,    'S','*',ring          ,fmtx,fmty)
      CALL GLK_plot2(  iafbG2,    'S','*',ring          ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  iafbO3,'S','*',dot               ,fmtx,fmty)

      IF( m_CMSene .LE. 200d0) THEN
      ENDIF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_sigHO
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int2-sigHO.eps
*//
*//   Plot of  dSigma/dvMax, Difference O(alf2)-O(alf1)
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(11)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1000}{\\Huge${\\sigma^{(2)}-\\sigma^{(1)}\\over\\sigma^{(1)}}$}',!
     $'\\PaveLb{600}{ 40}{\\huge $v_{\\max}$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveLr{800}{1000}{\\color{blue}\\line(1,0){150}\\; }',
     $'\\PaveL{ 800}{1000}{\\color{blue} IFI ON }',
     $'\\PaveLr{800}{ 900}{\\color{black}\\line(1,0){150}\\; }',
     $'\\PaveL{ 800}{ 900}{\\color{black} IFI OFF }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int2-sigHO.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetColor('\\color{green}\\thicklines$')

      CALL GLK_SetYminYmax(isigG2mG1,-0.012d0, 0.012d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(isigG2mG1,-0.003d0, 0.003d0)
      ENDIF
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      isigG2mG1,' ',' ',circle    ,fmtx,fmty) ! O(alf2-alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(   isigG2mG1nin,'S',' ',star      ,fmtx,fmty) ! O(alf1-alf1) intOFF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_afbHO
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int2-afbHO.eps
*//
*//   Plot of  AFB(vMax), Difference O(alf2)-O(alf1)
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(11)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{100}{1000}{\\Huge $A_{_{\\rm FB}}^{(2)}-A_{_{\\rm FB}}^{(1)}$ }',!
     $'\\PaveLb{600}{ 40}{\\huge $v_{\\max}$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveLr{800}{1000}{\\color{blue}\\line(1,0){150}\\; }',
     $'\\PaveL{ 800}{1000}{\\color{blue} IFI ON }',
     $'\\PaveLr{800}{ 900}{\\color{black}\\line(1,0){150}\\; }',
     $'\\PaveL{ 800}{ 900}{\\color{black} IFI OFF }',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.4'
      TeXfile   = 'afb_int2-afbHO.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetColor('\\color{green}\\thicklines$')
*
      CALL GLK_SetYminYmax(iafbG2mG1,-0.004d0, 0.004d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(iafbG2mG1,-0.003d0, 0.003d0)
      ENDIF
*
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      iafbG2mG1,' ',' ',circle    ,fmtx,fmty) ! O(alf2-alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(   iafbG2mG1nin,'S',' ',star      ,fmtx,fmty) ! O(alf1-alf1) intOFF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END

