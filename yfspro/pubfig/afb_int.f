* With older *.hst you may need to comment out *[[[[[ ....c]]]]
*/////////////////////////////////////////////////////////////////////////////////
*//   make all first part;       gmake afb-sig-ps
*//   make all second part;      gmake afb_int-ps
*//   make all second part;      gmake afb-ang-ps
*/////////////////////////////////////////////////////////////////////////////////
*// 
*//   make all first part;       gmake afb-sig-ps
*//   Table of everything;                               gmake afb_int-tab1.eps
*//   Plot of SigTot(vMax) versus reference EEX3(vMax);  gmake afb_int-Gsig.eps
*//   Plot of Afb(vMax) versus reference AfbRef of EEX3; gmake afb_int-Gafb.eps
*//
*//   make all second part;      gmake afb_int-ps
*//   Plot of SigmaInt/SigmaTot(vMax);                   gmake afb_int-sig1.eps
*//   Plot of AFB_Int(vmax);                             gmake afb_int-afb1.eps
*//   Plot of AFB(v), bin per bin dependence on v;       gmake afb_int-afb2.eps
*//
*//   make all second part;      gmake afb-ang-ps
*//   Plot of dSigma/dCosTheta  for v<0.90;              gmake afb_int-G1.eps
*//   Plot of dSigma/dCosTheta  for v<0.10;              gmake afb_int-G1x.eps
*//   Plot of dSigma/dCosTheta  for vp<0.10;             gmake afb_int-G1xxx.eps
*//
*//    gmake afb_int-sig1S.eps
*//    gmake afb_int-afb1S.eps
*/////////////////////////////////////////////////////////////////////////////////

* UNCOMMENT marked {{{{{ }}}}} for 91GeV
* Uncomment c\\\\\\ for less columns in table1 (no comparison with KORALZ)
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

      Tesnam    = 'afb_int'
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
c[[[[[
      CALL Plot_tab2
      CALL Plot_tab3
c]]]]]
      CALL Plot_G1
      CALL Plot_G1x
      CALL Plot_G1xxx
*
      CALL Plot_afb1
      CALL Plot_afb2
*
      CALL Plot_sig1
      CALL Plot_Gsig
      CALL Plot_GsigZF
      CALL Plot_Gafb
      CALL Plot_GafbZF
*
      CALL Plot_afb1S
      CALL Plot_sig1S
*
      CALL Plot_com1
      CALL Plot_com1x
      CALL Plot_com1xxx
*
      CALL Plot_AngMx
      CALL Plot_comMx
*
      CALL Plot_tabEWG1
      CALL Plot_tabEWG2
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
*//   gmake afb_int-tab1.eps
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
      TeXfile   = 'afb_int-tab1.txp'
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
c      nColumn=4   ! KORALZ eliminated
cc      nColumn=3 !!!!!!!!!!!!!
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


      SUBROUTINE Plot_tab2
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-tab2.eps
*//   Table of everything
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
! Parameters for tables
      INTEGER       idl(7)
      CHARACTER*80  capt(8)
      CHARACTER*16  fmt(3)
      CHARACTER*80  mcapt
      INTEGER       nColumn,ibin1,ibin2,iKf1,iKf2,iKfs,idum
      CHARACTER*32 TabLab1(5), TabLab2(1)
      DATA TabLab1 / ' $d$ ', ' $u$ ', ' $c$ ', ' $s$ ', ' $b$ '/ !
      DATA TabLab2 / ' $all$' /
*===================================================================
*                Writing table in LaTex format
      TeXfile   = 'afb_int-tab2.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      capt(1)='{\\color{blue} $f$}'
      capt(2)='{\\color{blue} (a) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm intOFF}$}'!
      capt(3)='{\\color{blue} (b) KORALZ 4.02 }'
      capt(4)='{\\color{blue} (c) KORALZ 4.04 }'
      capt(5)='{\\color{blue} (b-a)/a }'
      capt(6)='{\\color{blue} (c-a)/a }'
      fmt(1)='F5.0'
      fmt(2)='F10.4'
      fmt(3)='F8.4'
*------------------------------------------------
      iKFs = 7                  ! poiter of sum over flafours
      CALL Plot_SumFlav(iCeex2,    iKFs)
      CALL Plot_SumFlav(iKor402sig,iKFs)
      CALL Plot_SumFlav(iKor404sig,iKFs)
!------------------------------------------------
      idl(1)= iCeex2
      idl(2)= iKor402sig
      idl(3)= iKor404sig
      idl(4)= iKor402dif
      idl(5)= iKor404dif
      nColumn = 3
      nColumn = 5
      iKf1  = 1
      iKf2  = 5 
      CALL Plot_DifFlav(iCeex2,iKor402sig,iKor402dif,iKF1,iKFs)
      CALL Plot_DifFlav(iCeex2,iKor404sig,iKor404dif,iKF1,iKFs)
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.1$ }'!
      iCut = 2
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,' ','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*---
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.2$}'!
      iCut = 3
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*---
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.3$}'!
      iCut = 4
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*---
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.9$}'!
      iCut = 10
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*---
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.99$}'!
      iCut = 11
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*------------------------------------------------
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_tab3
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-tab3.eps
*//   All flavours, Table of KKMC versus Dizet 
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
! Parameters for tables
      INTEGER       idl(7)
      CHARACTER*80  capt(8)
      CHARACTER*16  fmt(3)
      CHARACTER*80  mcapt
      INTEGER       nColumn,ibin1,ibin2,iKf1,iKf2,iKfs,iKfm,idum
      CHARACTER*32  TabLab1(5), TabLab2(1), TabLab3(1)
      DATA TabLab1 / ' $d$ ', ' $u$ ', ' $s$ ', ' $c$ ', ' $b$ '/ !
      DATA TabLab2 / ' $all$' /
      DATA TabLab3 / ' $\\mu$' /
*===================================================================
*                Writing table in LaTex format
      TeXfile   = 'afb_int-tab3.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      capt(1)='{\\color{blue} $f$}'
      capt(2)='{\\color{blue} (a) ${\\cal KK}$sem }'
      capt(3)='{\\color{blue} (b) ${\\cal O}(\\alpha^2)_{\\rm CEEX}^{\\rm intOFF}$}' !
      capt(4)='{\\color{blue} (c) Zfitter 6.x }'
      capt(5)='{\\color{blue} (b-a)/a }'
      capt(6)='{\\color{blue} (c-b)/b }'
      fmt(1)='F5.0'
      fmt(2)='F10.4'
      fmt(3)='F8.4'
*------------------------------------------------
      iKFs = 7                  ! poiter of sum over flafours
      CALL Plot_SumFlav(iKKsem,     iKFs)
      CALL Plot_SumFlav(iCeex2,     iKFs)
      CALL Plot_SumFlav(iZft620sig, iKFs)
!------------------------------------------------
      idl(1)= iKKsem
      idl(2)= iCeex2
      idl(3)= iZft620sig
      idl(4)= iCeex2dif
      idl(5)= iZft620dif
      nColumn = 5
      iKf1  = 1
      iKf2  = 5 
      iKfm  = 13
      CALL Plot_DifFlav(iKKsem,iCeex2,    iCeex2dif ,iKF1,iKFm)
      CALL Plot_DifFlav(iCeex2,iZft620sig,iZft620dif,iKF1,iKFm)
*            $_________|_________|_________|_________|_________|_________|_________|_________|
c      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.01$ }'!
c      iCut = 1
c      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab1)
c      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,' ','L',' ')
c      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab2)
c      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
c      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab3)
c      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*--- 0.10
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.1$ }'!
      iCut = 2
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,' ','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab3)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*--- 0.20
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.2$}'!
      iCut = 3
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab3)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*--- 0.30
c      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.3$}'!
c      iCut = 4
c      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab1)
c      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
c      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab2)
c      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
c      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
c      CALL GLK_SetTabLab(5,TabLab3)
c      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*--- 0.90
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.9$}'!
      iCut = 10
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab3)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*--- 0.99
      Mcapt ='{\\color{red}$\\sigma(v_{\\max})$ [pb], $v_{\\max}=0.99$}'!
      iCut = 11
      CALL GLK_SetTabRan(  nCut*(iKF1-1)+iCut,  nCut*(iKF2-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab1)
      CALL GLK_PlTable2(nColumn, idl,capt,Mcapt,fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFs-1)+iCut,  nCut*(iKFs-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab2)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
      CALL GLK_SetTabRan(  nCut*(iKFm-1)+iCut,  nCut*(iKFm-1)+iCut,  nCut) !
      CALL GLK_SetTabLab(5,TabLab3)
      CALL GLK_PlTable2(nColumn, idl,capt,'   ',fmt,'S','L',' ')
*------------------------------------------------
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_SumFlav(iHist,iSumKF)
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER           iHist,iSumKF
      DOUBLE PRECISION  xbin(400),xerr(400),sum1,sum2
      INTEGER           iKf1,iKf2,idum
!------------------------------------------------
      CALL GLK_UnPak(iHist,xbin,'    ',idum)
      CALL GLK_UnPak(iHist,xerr,'ERRO',idum)
      DO iCut =1,nCut
         sum1=0d0
         sum2=0d0
         DO iKF =1,5
            ibin =  nCut*(iKF-1)+iCut
            sum1=sum1 +xbin(ibin)
            sum2=sum2 +xerr(ibin)**2
         ENDDO
cc         iSumKF =  6
         ibin =  nCut*(iSumKF-1)+iCut
         xbin(ibin)=sum1
         xerr(ibin)=SQRT(sum2)
      ENDDO
      CALL GLK_Pak( iHist,xbin)
      CALL GLK_Pake(iHist,xerr)
      END

      SUBROUTINE Plot_DifFlav(iHis1,iHis2,iHis3,iKF1,iKF2)
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER           iHis1,iHis2,iHis3
      INTEGER           iKf1,iKf2,idum
      DOUBLE PRECISION  xbin1(400),xerr1(400),xbin2(400),xerr2(400),xbin3(400),xerr3(400) !
      DOUBLE PRECISION  sum1,sum2
!------------------------------------------------
      CALL GLK_UnPak(iHis1,xbin1,'    ',idum)
      CALL GLK_UnPak(iHis1,xerr1,'ERRO',idum)
      CALL GLK_UnPak(iHis2,xbin2,'    ',idum)
      CALL GLK_UnPak(iHis2,xerr2,'ERRO',idum)
      CALL GLK_Clone1(iHis1,iHis3,   'KORALZ 4.x difference with KKMC $')
      DO iCut =1,nCut
         DO iKF =iKF1,iKF2
            ibin =  nCut*(iKF-1)+iCut
            xbin3(ibin)=0d0
            xerr3(ibin)=0d0
            IF(xbin1(ibin).NE.0d0) THEN
               xbin3(ibin) =  (xbin2(ibin)-xbin1(ibin)) /xbin1(ibin)
               xerr3(ibin) =  SQRT(xerr1(ibin)**2+xerr2(ibin)**2)/xbin1(ibin) !
            ENDIF
         ENDDO
      ENDDO
      CALL GLK_Pak( iHis3,xbin3)
      CALL GLK_Pake(iHis3,xerr3)
      END

      SUBROUTINE Plot_Gsig
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int-Gsig.eps
*//
*//   Plot of SigTot(vMax) versus reference EEX3(vMax)
*//   (SigTot-SigRef)/SigRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER                iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(19)
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
     $'\\PaveL{550}{ 880}{\\large\\color{black} $\\times$ EEX2 of KORALZ 4.03}',!
     $'\\PaveL{550}{ 820}{\\large\\color{red} ',!
     $'         $\\star$\\; CEEX2 $+$ IFI at ${\\cal O}(\\alpha^1)$  }',!
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

      ymin= -0.025d0
      ymax=  0.100d0
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
c[[[
      CALL GLK_SetYminYmax(isigG1nin,-0.020d0, 0.020d0) !!!! {{{{ 91GeV hadrons
      CALL GLK_SetYminYmax(isigG2nin,-0.020d0, 0.020d0) !!!! {{{{ 91GeV hadrons
      CALL GLK_SetYminYmax(isigG1nin,-0.060d0, 0.060d0) !!!! {{{{ 91GeV hadrons
      CALL GLK_SetYminYmax(isigG2nin,-0.060d0, 0.060d0) !!!! {{{{ 91GeV hadrons
      CALL GLK_SetYminYmax(isigG1nin,-0.200d0, 0.200d0) !!!! {{{{ 91GeV hadrons
      CALL GLK_SetYminYmax(isigG2nin,-0.200d0, 0.200d0) !!!! {{{{ 91GeV hadrons
c]]]

      CALL GLK_Operat(isigG1,    '-',  iRef,     isigG1,       1d0, 1d0)
      CALL GLK_Operat(isigG1,    '/',  iRef,     isigG1,       1d0, 1d0)

      CALL GLK_Operat(isigG2,    '-',  iRef,     isigG2,       1d0, 1d0)
      CALL GLK_Operat(isigG2,    '/',  iRef,     isigG2,       1d0, 1d0)

      CALL GLK_Operat(iKorzSig,   '-', iRef,     iKorzSig,      1d0, 1d0)
      CALL GLK_Operat(iKorzSig,   '/', iRef,     iKorzSig,      1d0, 1d0)
      CALL GLK_idopt( iKorzSig,'ERRO')
*--------------------------------------------------------------------------
      TeXfile   = 'afb_int-Gsig.txp'
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
* KORALZ
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  iKorzSig,   'S','*',times        ,fmtx,fmty)
* Combined
      CALL GLK_SetColor('\\color{red}$')
c      CALL GLK_plot2(  isigComb,   'S','*',star        ,fmtx,fmty)
      CALL GLK_plot2(  iHyb1Sig,   'S','*',star        ,fmtx,fmty)
c      CALL GLK_plot2(  iHyb2Sig,   'S','*',star        ,fmtx,fmty)
      ENDIF
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
*--------------------------------------------------------------------------
      END


      SUBROUTINE Plot_GsigZF
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int-GsigZF.eps
*//
*//   Plot of SigTot(vMax) versus reference EEX3(vMax)
*//   (SigTot-SigRef)/SigRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER                iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(19)
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
     $'\\PaveL{550}{ 920}{\\large\\color{red}  $\\star$\\;    ZFITTER 6.x Int.ON}',!
     $'\\PaveL{550}{ 860}{\\large\\color{blue} $\\diamond$\\; ZFITTER 6.x Int.OFF}',!
     $'} % -- End Pave',
     $'% ',!
     $'% ',!
     $'% ',!
     $'% ',!
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*
ccc      iRef = isigO3
      iRef =iO3best
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-GsigZF.txp'
      CALL GLK_PlInitialize(2,TeXfile)

      ymin= -0.025d0
      ymax=  0.050d0
      ymax=  0.110d0
      CALL GLK_SetYminYmax(isigG2nin,ymin,ymax) !!!! LEP2, leptons
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(isigG2nin,-0.004d0, 0.004d0) !!!! {{{{ 91GeV leptons
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(isigG1nin,-0.009d0, 0.009d0) !!!! {{{{ 91GeV hadrons
         CALL GLK_SetYminYmax(isigG2nin,-0.009d0, 0.009d0) !!!! {{{{ 91GeV hadrons
      ENDIF

      CALL GLK_Operat(iDizSigNin, '-', iRef,     iDizSigNin,   1d0, 1d0)
      CALL GLK_Operat(iDizSigNin, '/', iRef,     iDizSigNin,   1d0, 1d0)
      CALL GLK_idopt( iDizSigNin,'ERRO')

      CALL GLK_Operat(iDizSig,    '-', iRef,     iDizSig,      1d0, 1d0)
      CALL GLK_Operat(iDizSig,    '/', iRef,     iDizSig,      1d0, 1d0)
      CALL GLK_idopt( iDizSig,'ERRO')

      CALL GLK_SetColor('\\color{blue}\\thicklines$')
      CALL GLK_plot2(  isigG2nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  isigG2,    'S','*',ring          ,fmtx,fmty)

      IF( m_CMSene .LE. 200d0) THEN
* ZFITTER
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  iDizSigNin,'S','*',diamond       ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  iDizSig,   'S','*',star          ,fmtx,fmty)
      ENDIF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd

      END



      SUBROUTINE Plot_Gafb
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int-Gafb.eps
*//
*//   Plot of Afb(vMax) versus reference AfbRef of EEX3
*//   (Afb-AfbRef)/AfbRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER  iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(20)
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
     $'\\PaveL{550}{ 880}{\\large\\color{black}    $\\times$  EEX2, KORALZ 4.03}',!
     $'\\PaveL{550}{ 820}{\\large\\color{red} ',!
     $'         $\\star$\\; CEEX2 $+$ IFI at ${\\cal O}(\\alpha^1)$  }',!
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

      ymin= -0.025d0
      ymax=  0.060d0
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
      TeXfile   = 'afb_int-Gafb.txp'
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
* KORALZ 
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(  iKorzAfb,   'S','*',times        ,fmtx,fmty)
* Combined
      CALL GLK_SetColor('\\color{red}$')
c      CALL GLK_plot2(  iafbComb,   'S','*',star         ,fmtx,fmty)
      CALL GLK_plot2(  iHyb1Afb,   'S','*',star         ,fmtx,fmty)
c      CALL GLK_plot2(  iHyb2Afb,   'S','*',star         ,fmtx,fmty)
      ENDIF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_GafbZF
*//////////////////////////////////////////////////////////////////
*//   gmake afb_int-GafbZF.eps
*//
*//   Plot of Afb(vMax) versus reference AfbRef of EEX3
*//   (Afb-AfbRef)/AfbRef as function of vMax
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER  iRef
*---------------------------------------------------------------------------
      CHARACTER*80 Label(20)
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
     $'% ',!
     $'% ',!
     $'% ',!
     $'% ',!
     $'\\PaveL{550}{ 920}{\\large\\color{red}  $\\star$\\;    ZFITTER 6.x Int.ON }',!
     $'\\PaveL{550}{ 860}{\\large\\color{blue} $\\diamond$\\; ZFITTER 6.x Int.OFF}',!
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
*

      iRef = iO3bestAFB
      ymin= -0.035d0
      ymax=  0.090d0
cc      ymax=  0.150d0
      CALL GLK_SetYminYmax(iafbG1nin,ymin, ymax)
      CALL GLK_SetYminYmax(iafbG2nin,ymin, ymax)
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(iafbG1nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
         CALL GLK_SetYminYmax(iafbG2nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(iafbG1nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
         CALL GLK_SetYminYmax(iafbG2nin,-0.005d0, 0.007d0) !!!! {{{{ 91GeV }}}}
      ENDIF

      CALL GLK_Operat(iDizAfbNin, '-', iRef,     iDizAfbNin,   1d0, 1d0)
      CALL GLK_idopt( iDizAfbNin,'ERRO')
      CALL GLK_Operat(iDizAfb,    '-', iRef,     iDizAfb,      1d0, 1d0)
      CALL GLK_idopt( iDizAfb,'ERRO')
* Plotting
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-GafbZF.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{blue}\\thicklines$')
cc      CALL GLK_plot2(  iafbG1nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_plot2(  iafbG2nin, ' ',' ',dot           ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
cc      CALL GLK_plot2(  iafbG1,    'S','*',ring          ,fmtx,fmty)
      CALL GLK_plot2(  iafbG2,    'S','*',ring          ,fmtx,fmty)
* Zfitter
      IF( m_CMSene .LE. 200d0) THEN
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(  iDizAfbNin,'S','*',diamond       ,fmtx,fmty)
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  iDizAfb,   'S','*',star          ,fmtx,fmty)
      ENDIF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_sig1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-sig1.eps
*//
*//   Plot of  SigmaInt/SigmaTot(vMax)
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(20)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
ccc     $'\\PaveL{ 60}{1100}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$ }',
     $'\\PaveL{250}{1000}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$ }',
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{ 700}{1070}{\\color{green} $\\;{\\cal O}(\\alpha^2)_{\\rm CEEX}$ }',
     $'\\PaveLr{700}{1070}{\\color{green} \\line(1,0){150}} ',
     $'\\PaveLr{700}{ 970}{\\color{red} $\\circle{20}$\\;\\;$\\circle{20}$\\;\\; }',
     $'\\PaveL{ 700}{ 970}{\\color{red}\\large KORALZ 1-st ord. }',
     $'\\PaveLr{700}{ 870}{\\color{red} $\\star$\\;$\\star$\\;$\\star$\\; }',
     $'\\PaveL{ 700}{ 870}{\\color{red}\\large ZFITTER 6.x}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|

* We may divide O(alf1) interf. by O(alf1) or O(alf3) sigma, small difference 
cc      CALL GLK_Operat(iMustSigInt,'/',iMustSig   , idWork, 1d0, 1d0)
      CALL GLK_Operat(iMustSigInt,'/',isigG2nin     , idWork, 1d0, 1d0) !
cc      CALL GLK_Operat(iMustSigInt,'/',isigO3     , idWork, 1d0, 1d0)

      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-sig1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYminYmax(isigG1int,-0.02d0, 0.10d0)
      CALL GLK_SetYminYmax(isigG2int,-0.02d0, 0.10d0)
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(isigG1int,-0.010d0, 0.010d0) !!!! {{{{{  91GeV }}}}}
         CALL GLK_SetYminYmax(isigG2int,-0.010d0, 0.010d0) !!!! {{{{{  91GeV }}}}}
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(isigG1int,-0.002d0, 0.002d0) !!!! hadrons
         CALL GLK_SetYminYmax(isigG2int,-0.002d0, 0.002d0) !!!! hadrons
      ENDIF


      CALL GLK_SetColor('\\color{green}\\thicklines$')
cc      CALL GLK_plot2(  isigG1int,' ',' ',dot       ,fmtx,fmty)
      CALL GLK_plot2(  isigG2int,' ',' ',dot       ,fmtx,fmty)

ccc      CALL GLK_SetYminYmax(isigS2int,-0.02d0, 0.10d0)
ccc      CALL GLK_plot2(  isigS2int,' ',' ',dot       ,fmtx,fmty)

c      CALL GLK_SetColor('\\color{blue}\\thicklines$')
c      CALL GLK_plot2(  isigG0int,' ',' ',circle    ,fmtx,fmty)

      IF( m_CMSene .LE. 200d0) THEN
      IF(m_iscale .NE. 2) THEN
         CALL GLK_SetColor('\\color{red}$')
         CALL GLK_plot2(  idWork,'S','*',ring     ,fmtx,fmty)
         CALL GLK_SetColor('\\color{red}$')
         CALL GLK_plot2(  iDizSigInt,'S','*',star     ,fmtx,fmty)
      ENDIF
      ENDIF
      CALL GLK_Delet(idWork)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_afb1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-afb1.eps
*//
*//   Plot of AFB_Int(vmax)
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(20)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
ccc     $'\\PaveLr{-10}{1100}{\\Huge $A_{_{\\rm FB}}^{\\rm int}$ }',
     $'\\PaveL{100}{1000}{\\Huge $A_{_{\\rm FB}}^{\\rm int}$ }',
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime_{\\min}/s$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{ 700}{1000}{\\color{green} ${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }',
     $'\\PaveLr{700}{1000}{\\color{green}\\line(1,0){150}\\; }',
     $'\\PaveLr{700}{ 900}{\\color{red} $\\circle{20}$\\;\\;$\\circle{20}$\\;\\; }',
     $'\\PaveL{ 700}{ 900}{\\color{red}\\large KORALZ 1-st ord. }',
     $'\\PaveLr{700}{ 800}{\\color{red} $\\star$\\;$\\star$\\; }',
     $'\\PaveL{ 700}{ 800}{\\color{red}\\large ZFITTER 6.x}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-afb1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYminYmax(iafbG0int,-0.02d0, 0.10d0)
      CALL GLK_SetYminYmax(iafbG1int,-0.02d0, 0.10d0)
      CALL GLK_SetYminYmax(iafbG2int,-0.02d0, 0.10d0)
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(iafbG0int,-0.009d0, 0.009d0) !!!! {{{{{  91GeV }}}}}
         CALL GLK_SetYminYmax(iafbG1int,-0.009d0, 0.009d0) !!!! {{{{{  91GeV }}}}}
         CALL GLK_SetYminYmax(iafbG2int,-0.009d0, 0.009d0) !!!! {{{{{  91GeV }}}}}
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(iafbG0int,-0.002d0, 0.002d0) !!!! {{{{{  91GeV }}}}}
         CALL GLK_SetYminYmax(iafbG1int,-0.002d0, 0.002d0) !!!! {{{{{  91GeV }}}}}
         CALL GLK_SetYminYmax(iafbG2int,-0.002d0, 0.002d0) !!!! {{{{{  91GeV }}}}}
      ENDIF
      CALL GLK_SetColor('\\color{green}\\thicklines$')
cc      CALL GLK_plot2(   iafbG0int,' ',' ',dot       ,fmtx,fmty) ! O(alf0)
cc      CALL GLK_plot2(   iafbG1int,' ',' ',dot       ,fmtx,fmty) ! O(alf1)
      CALL GLK_plot2(   iafbG2int,' ',' ',dot       ,fmtx,fmty) ! O(alf2)

c[[[      CALL GLK_SetYminYmax(iafbS2int,-0.02d0, 0.10d0)
c[[[      CALL GLK_plot2(   iafbS2int,' ',' ',dot       ,fmtx,fmty) ! O(alf1)

      IF( m_CMSene .LE. 200d0) THEN
      IF( m_iscale .NE. 2) THEN
         CALL GLK_SetColor('\\color{red}$')
         CALL GLK_plot2(  iMustAfbInt,'S','*',ring     ,fmtx,fmty)
         CALL GLK_SetColor('\\color{red}$')
         CALL GLK_plot2(  iDizAfbInt, 'S','*',star     ,fmtx,fmty)
      ENDIF
      ENDIF

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_afb2
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-afb2.eps
*//
*//   Plot of AFB(v), bin per bin dependence on v
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{ 40}{1000}{\\Huge $A_{_{\\rm FB}}^{\\rm int}$ }',
     $'\\PaveLb{600}{ 40}{\\huge $1-s^\\prime/s$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{ 700}{1000}{\\color{green} ${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }',
     $'\\PaveLr{700}{1000}{\\color{green}\\line(1,0){150}\\; }',
****     $'\\PaveLr{700}{ 900}{\\color{red} $\\star$\\;$\\star$\\;$\\star$\\; }',
****     $'\\PaveL{ 700}{ 900}{\\color{red}\\large ZFITTER 6.x}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-afb2.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYminYmax(kafbG2int,-0.10d0, 0.10d0)
      IF( m_iscale .EQ. 1) THEN
         CALL GLK_SetYminYmax(kafbG2int,-0.08d0, 0.05d0) !!!! {{{{{  91GeV }}}}}
      ENDIF
      IF( m_iscale .EQ. 2) THEN
         CALL GLK_SetYminYmax(kafbG2int,-0.08d0, 0.05d0) !!!! {{{{{  91GeV }}}}}
      ENDIF
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_plot2(   kafbG2int,' ',' ',dot       ,fmtx,fmty) ! O(alf0)

***      CALL GLK_SetColor('\\color{red}$')
***      CALL GLK_plot2(  iDizAfbInt,'S','*',star      ,fmtx,fmty)

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_G1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-G1.eps
*//
*//   plot of dSigma/dCosTheta  for v<0.90
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{810}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $s^\\prime/s>0.1$ }',
     $'\\PaveL{ 50}{720}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{720}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{650}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{650}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelZ(10)   ! upper position
      DATA LabelZ /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{310}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $s^\\prime/s>0.1$ }',
     $'\\PaveL{ 50}{220}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{220}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{150}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{150}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-G1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYmin(iangG2,0d0)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_plot2(  iangG2,' ',' ',dot       ,fmtx,fmty)  !! CEEX intON
      CALL GLK_SetColor('\\color{red}\\thinlines$')
      CALL GLK_plot2(  iangN2,'S',' ',circle    ,fmtx,fmty)  !!! CEEX intOFF
*
cc      CALL GLK_SetColor('\\color{blue}$')
cc      CALL GLK_plot2(  iangO2,'S',' ',circle    ,fmtx,fmty)  !! EEX ofsolete
*
      CALL GLK_PlLabel(LabelBasic)
      IF( m_iscale .EQ. 0) THEN
         CALL GLK_PlLabel(Label)
      ELSE
         CALL GLK_PlLabel(LabelZ)
      ENDIF
      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_G1x
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-G1x.eps
*//
*//   plot of dSigma/dCosTheta  for v<0.10
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{810}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $s^\\prime/s>0.81$ }',
     $'\\PaveL{ 50}{720}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{720}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{650}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{650}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelZ(10)   ! upper position
      DATA LabelZ /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{310}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $s^\\prime/s>0.81$ }',
     $'\\PaveL{ 50}{220}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{220}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{150}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{150}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-G1x.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetYmin(iangG2x,0d0)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_plot2(  iangG2x,' ',' ',dot       ,fmtx,fmty)  !! CEEX intON
      CALL GLK_SetColor('\\color{red}$')
      CALL GLK_plot2(  iangN2x,'S',' ',circle    ,fmtx,fmty)  !! CEEX intOFF
ccccc      CALL GLK_SetColor('\\color{blue}$')
ccccc      CALL GLK_plot2(  iangO2x,'S',' ',circle    ,fmtx,fmty) !! EEX, obsolete
      CALL GLK_PlLabel(LabelBasic)
      IF( m_iscale .EQ. 0) THEN
         CALL GLK_PlLabel(Label)
      ELSE
         CALL GLK_PlLabel(LabelZ)
      ENDIF
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_G1xxx
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-G1xxx.eps
*//
*//   plot of dSigma/dCosTheta  for vp<0.2775
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{810}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $v_p<0.2775$ }',
     $'\\PaveL{ 50}{720}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{720}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{650}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{650}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelZ(10)   ! upper position
      DATA LabelZ /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{310}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $v_p<0.2775$ }',
     $'\\PaveL{ 50}{220}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{220}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{150}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{150}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-G1xxx.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYmin(iangG2xxx,0d0)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_plot2(  iangG2xxx,' ',' ',dot       ,fmtx,fmty)  !! CEEX intON
      CALL GLK_SetColor('\\color{red}\\thicklines$')
      CALL GLK_plot2(  iangN2xxx,'S',' ',circle    ,fmtx,fmty)  !! CEEX intOFF
*
      CALL GLK_PlLabel(LabelBasic)
      IF( m_iscale .EQ. 0) THEN
         CALL GLK_PlLabel(Label)
      ELSE
         CALL GLK_PlLabel(LabelZ)
      ENDIF
      CALL GLK_PlEnd
      END




      SUBROUTINE Plot_sig1S
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-sig1S.eps
*//
*//   Plot of  SigmaInt/SigmaTot(vMax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
ccc     $'\\PaveL{ 60}{1100}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$ }',
     $'\\PaveL{250}{1000}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$ }',
     $'\\PaveLb{600}{ 40}{\\huge $v_p$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{ 700}{1000}{\\color{green} $\\;{\\cal O}(\\alpha^2)_{\\rm CEEX}$ }',
     $'\\PaveLr{700}{1000}{\\color{green} \\line(1,0){150}} ',
cc     $'\\PaveLr{700}{ 900}{\\color{red} $\\circle{20}$\\;\\;$\\circle{20}$\\;\\; }',
cc     $'\\PaveL{ 700}{ 900}{\\color{red}\\large KORALZ 1-st ord. }',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-sig1S.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetColor('\\color{green}\\thicklines$')

      CALL GLK_SetYminYmax(isigS2int,-0.02d0, 0.10d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(isigS2int,-0.003d0, 0.010d0)
      ENDIF
      CALL GLK_plot2(      isigS2int,' ',' ',dot       ,fmtx,fmty)

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_afb1S
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-afb1S.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
ccc     $'\\PaveLr{-10}{1100}{\\Huge $A_{_{\\rm FB}}^{\\rm int}$ }',
     $'\\PaveL{100}{1000}{\\Huge $A_{_{\\rm FB}}^{\\rm int}$ }',
     $'\\PaveLb{600}{ 40}{\\huge $v_p$}',
     $'\\PaveLb{ 200}{ 20}{\\large\\color{blue} $\\leftarrow$ Strong Cut}',
     $'\\PaveLb{1000}{ 20}{\\large\\color{blue} No Cut $\\rightarrow$}',
*
     $'\\PaveL{ 700}{1000}{\\color{green} ${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }',
     $'\\PaveLr{700}{1000}{\\color{green}\\line(1,0){150}\\; }',
cc     $'\\PaveLr{700}{ 900}{\\color{red} $\\circle{20}$\\;\\;$\\circle{20}$\\;\\; }',
cc     $'\\PaveL{ 700}{ 900}{\\color{red}\\large KORALZ 1-st ord. }',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-afb1S.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetColor('\\color{green}\\thicklines$')

      CALL GLK_SetYminYmax(iafbS2int,-0.02d0, 0.10d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(iafbS2int,-0.003d0, 0.010d0)
      ENDIF
      CALL GLK_plot2(      iafbS2int,' ',' ',dot       ,fmtx,fmty) ! O(alf1)

      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_com1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-com1.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveLr{150}{1000}{\\Large\\color{blue}$\\circle{30}\\;\\;\\;\\circle{30}$}',
     $'\\PaveL{ 200}{1000}{\\Huge\\color{blue} $A_{_{\\rm FB}}^{\\rm int}$  }',
*
     $'\\PaveLr{550}{1000}{\\Large  $\\star$  $\\star$ }',
     $'\\PaveL{ 600}{1000}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$}',
*
     $'\\PaveLb{600}{ 40}{\\huge $\\cos\\theta_{\\max}$}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-com1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_SetYminYmax(igAfbG2,-0.02d0, 0.04d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(igAfbG2,-0.004d0, 0.004d0)
      ENDIF
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      igAfbG2,' ','*',circle    ,fmtx,fmty) ! O(alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(      igSigG2,'S','*',star      ,fmtx,fmty) ! O(alf1)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_com1x
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-com1x.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveLr{150}{1000}{\\Large\\color{blue}$\\circle{30}\\;\\;\\;\\circle{30}$}',
     $'\\PaveL{ 200}{1000}{\\Huge\\color{blue} $A_{_{\\rm FB}}^{\\rm int}$  }',
*
     $'\\PaveLr{550}{1000}{\\Large  $\\star$  $\\star$ }',
     $'\\PaveL{ 600}{1000}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$}',
*
     $'\\PaveLb{600}{ 40}{\\huge $\\cos\\theta_{\\max}$}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              ',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-com1x.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_SetYminYmax(igAfbG2x,-0.02d0, 0.04d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(igAfbG2x,-0.004d0, 0.004d0)
      ENDIF
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      igAfbG2x,' ','*',circle    ,fmtx,fmty) ! O(alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(      igSigG2x,'S','*',star      ,fmtx,fmty) ! O(alf1)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_com1xxx
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-com1xxx.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\Pave{',
     $'\\PaveLr{150}{1000}{\\Large\\color{blue}$\\circle{30}\\;\\;\\;\\circle{30}$}',
     $'\\PaveL{ 200}{1000}{\\Huge\\color{blue} $A_{_{\\rm FB}}^{\\rm int}$  }',
     $'\\PaveL{ 50}{810}{${\\cal O}(\\alpha^2)_{\\rm CEEX},$ $v_p<0.2775$ }',
*
     $'\\PaveLr{550}{1000}{\\Large $\\star$  $\\star$ }',
     $'\\PaveL{ 600}{1000}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$}',
*
     $'\\PaveLb{600}{ 40}{\\huge $\\cos\\theta_{\\max}$}',
     $'} % -- End Pave',
     $'% end-of-label',
     $'              '/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-com1xxx.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_SetYminYmax(igAfbG2xxx,-0.02d0, 0.04d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(igAfbG2xxx,-0.004d0, 0.004d0)
      ENDIF
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      igAfbG2xxx,' ','*',circle    ,fmtx,fmty) ! O(alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(      igSigG2xxx,'S','*',star      ,fmtx,fmty) ! O(alf1)
      CALL GLK_PlLabel(LabelBasic)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_AngMx
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-AngMx.eps
*//
*//   plot of dSigma/dCosTheta  for v<0.1
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(11)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{810}{${\\cal O}(\\alpha^1)$ KORALZ $s^\\prime/s>0.9$ }',
     $'\\PaveL{ 50}{720}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{720}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{650}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{650}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelZ(11)   ! upper position
      DATA LabelZ /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $'\\PaveL{150}{1020}{\\Huge ${d\\sigma \\over d\\cos\\theta }$ }',
     $'\\PaveL{ 50}{310}{${\\cal O}(\\alpha^1)$ KORALZ $s^\\prime/s>0.9$ }',
     $'\\PaveL{ 50}{220}{\\color{green}\\line(1,0){150}}',
     $'\\PaveL{250}{220}{\\color{green} ISR*FSR  ON }',
     $'\\PaveL{ 50}{150}{\\color{red}\\line(1,0){150}}',
     $'\\PaveL{250}{150}{\\color{red} ISR*FSR  OFF }',
     $'\\PaveLb{1000}{ 40}{\\Huge $\\cos\\theta$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-AngMx.txp'
      CALL GLK_PlInitialize(2,TeXfile)
*
      CALL GLK_SetYmin(iMustAng,0d0)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_plot2(  iMustAng,    ' ',' ',dot       ,fmtx,fmty)  !! KORALZ intON
      CALL GLK_SetColor('\\color{red}\\thicklines$')
      CALL GLK_plot2(  iMustAngNint,'S',' ',circle    ,fmtx,fmty)  !! KORALZ intOFF
*
      IF( m_iscale .EQ. 0) THEN
         CALL GLK_PlLabel(Label)
      ELSE
         CALL GLK_PlLabel(LabelZ)
      ENDIF
      CALL GLK_PlEnd
      END



      SUBROUTINE Plot_comMx
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-comMx.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
      CHARACTER*80 Label(10)
      DATA Label /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $'\\PaveL{ 50}{1100}{${\\cal O}(\\alpha^1)$ KORALZ $s^\\prime/s>0.9$ }',
*
     $'\\PaveLr{150}{950}{\\Large\\color{blue}$\\circle{30}\\;\\;\\;\\circle{30}$}',
     $'\\PaveL{ 200}{950}{\\Huge\\color{blue} $A_{_{\\rm FB}}^{\\rm int}$  }',
*
     $'\\PaveLr{550}{950}{\\Large $\\star$  $\\star$ }',
     $'\\PaveL{ 600}{950}{\\Huge ${\\sigma^{\\rm int}\\over\\sigma^{\\rm tot.}}$}',
*
     $'\\PaveLb{600}{ 40}{\\huge $\\cos\\theta_{\\max}$}',
     $'} % -- End Pave',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      fmtx='f10.2'
      fmty='f10.3'
      TeXfile   = 'afb_int-comMx.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL GLK_SetColor('\\color{green}\\thicklines$')
      CALL GLK_SetYminYmax(jMustAfb,-0.02d0, 0.05d0)
      IF( m_iscale .GT. 0) THEN
         CALL GLK_SetYminYmax(jMustAfb,-0.004d0, 0.004d0)
      ENDIF
      CALL GLK_SetColor('\\color{blue}$')
      CALL GLK_plot2(      jMustAfb,' ','*',circle    ,fmtx,fmty) ! O(alf1)
      CALL GLK_SetColor('\\color{black}$')
      CALL GLK_plot2(      jMustSig,'S','*',star      ,fmtx,fmty) ! O(alf1)
      CALL GLK_PlLabel(Label)
      CALL GLK_PlEnd
      END


      SUBROUTINE Plot_tabEWG1
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-tabEWG1.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3)
      CHARACTER*80  mcapt
*---------------------------------------------------------------------------
      fmt(1)='F10.2'
      fmt(2)='F10.4'
      fmt(3)='F8.4'

      TeXfile   = 'afb_int-tabEWG1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      capt(1)='{\\color{blue}\\bf $\\cos\\vartheta_{\\max}$ }'
      capt(2)='{\\color{blue}\\bf (5-1) }'
      capt(3)='{\\color{blue}\\bf (3-5) }'
      capt(4)='{\\color{blue}\\bf (3-1) }'
      capt(5)='{\\color{blue}\\bf (4-3) }'
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}  Differences (relative) in cross section, ${\\cal KK}$M.C.       }' !
      idl(1)= igSig51
      idl(2)= igSigG2x
      idl(3)= igSig31
      idl(4)= igSig43
      CALL GLK_SetTabRan(18,20,1)
      CALL GLK_PlTable2(4,idl,capt,Mcapt,fmt,' ','R',' ')

*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}  Differences in charge asymmetry, ${\\cal KK}$M.C.               }' !
      idl(1)= igAfb51
      idl(2)= igAfbG2x
      idl(3)= igAfb31
      idl(4)= igAfb43
      CALL GLK_SetTabRan(18,20,1)
      CALL GLK_PlTable2(4,idl,capt,Mcapt,fmt,'S','R',' ')

      CALL GLK_PlEnd
      END

      SUBROUTINE Plot_tabEWG2
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-tabEWG2.eps
*//
*//   Plot of AFB_Int(vmax) where vMax is for ISR only
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
*---------------------------------------------------------------------------
! Parameters for tables
      INTEGER       idl(5)
      CHARACTER*80  capt(6)
      CHARACTER*16  fmt(3)
      CHARACTER*80  mcapt
*---------------------------------------------------------------------------
      fmt(1)='F10.2'
      fmt(2)='F10.4'
      fmt(3)='F8.4'

      TeXfile   = 'afb_int-tabEWG2.txp'
      CALL GLK_PlInitialize(2,TeXfile)

      CALL MakeAng(iangG2x,   igAfbG2x,   igSigG2x)  !(3) recycling
      CALL AngCut(            igSigG2x,   igAfbG2x)  !(3) recycling
      CALL MakeAng(iangG2xx,  igAfbG2xx,  igSigG2xx) !(4) recycling
      CALL AngCut(            igSigG2xx,  igAfbG2xx) !(4) recycling
      CALL MakeAng(iangN2x,   igAfb51,    igSig51)   !(5) recycling
      CALL AngCut(            igSig51,    igAfb51)   !(5) recycling
      CALL MakeAng(iangN2prop,igAfb43,    igSig43)   !(1) recycling
      CALL AngCut(            igSig43,    igAfb43)   !(1) recycling

      capt(1)='{\\color{blue}\\bf $\\cos\\vartheta_{\\max}$ }'
      capt(2)='{\\color{blue}\\bf (1) }'
      capt(3)='{\\color{blue}\\bf (5) }'
      capt(4)='{\\color{blue}\\bf (3) }'
      capt(5)='{\\color{blue}\\bf (4) }'
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}  Angular cut in cross section, ${\\cal KK}$M.C.       }' !
      idl(1)= igSig43           !(1)
      idl(2)= igSig51           !(5)
      idl(3)= igSigG2x          !(3)
      idl(4)= igSigG2xx         !(4)
      CALL GLK_SetTabRan(18,20,1)
      CALL GLK_PlTable2(4,idl,capt,Mcapt,fmt,' ','R',' ')
*            $_________|_________|_________|_________|_________|_________|_________|_________|
      Mcapt ='{\\color{red}  Angular cut in charge asymmetry, ${\\cal KK}$M.C.       }' !
      idl(1)= igAfb43           !(1)
      idl(2)= igAfb51           !(5)
      idl(3)= igAfbG2x          !(3)
      idl(4)= igAfbG2xx         !(4)
      CALL GLK_SetTabRan(18,20,1)
      CALL GLK_PlTable2(4,idl,capt,Mcapt,fmt,'S','R',' ')

      CALL GLK_PlEnd

      END


      SUBROUTINE MakeAng(idAng,idAfb,idSig)
*/////////////////////////////////////////////////////////////////////////////////
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER          idAng,idAfb,idSig
      INTEGER          Nb1,Nb2,Nb,ib
      DOUBLE PRECISION XsForw,XeForw,Xsback,Xeback !! ON
      DOUBLE PRECISION Afb(100),AfbErr(100)
      DOUBLE PRECISION Sig(100),SigErr(100)
      DOUBLE PRECISION GLK_hi,GLK_hie
      LOGICAL GLK_Exist
*     -----------------------------------------------------------
      IF( .NOT.GLK_Exist(idAng)) THEN
         WRITE(*,*) ' ##### MakeAng, histos do not exist idAng=',idAng
         STOP
      ENDIF
      CALL GLK_GetNb(idAng,Nb1)
      IF( MOD(NB1,2).NE.0 ) GOTO 900
      Nb=Nb1/2
      CALL GLK_Book1(idAfb,  'Afb $',Nb,0d0,1d0)
      CALL GLK_Book1(idSig,  'Sig $',Nb,0d0,1d0)
      XsForw =0d0
      XeForw =0d0
      Xsback =0d0
      Xeback =0d0
      DO ib=1,Nb
         XsForw =  XsForw+ GLK_hi( idAng,NB+ib)
         XeForw =  XeForw+ GLK_hie(idAng,NB+ib)**2
         Xsback =  Xsback+ GLK_hi( idAng,NB-ib+1)
         Xeback =  Xeback+ GLK_hie(idAng,NB-ib+1)**2
cc         WRITE(*,*) '    ib,YsForw,Ysback= ',ib,XsForw,Xsback
         Afb(ib)    = (XsForw-Xsback)/(XsForw+Xsback)          !! ON
         AfbErr(ib) = SQRT(XeForw + Xeback)/(XsForw+Xsback)
         Sig(ib)    = XsForw+Xsback
         SigErr(ib) = SQRT(XeForw + Xeback)
      ENDDO
      CALL GLK_Pak(  idAfb,Afb)
      CALL GLK_Pake( idAfb,AfbErr)
      CALL GLK_Pak(  idSig,Sig)
      CALL GLK_Pake( idSig,SigErr)
*
cc      WRITE(6,*) '################## this is MakeAng ##################'
cc      CALL GLK_SetNout(6)
      CALL GLK_Print(idAfb )
      CALL GLK_Print(idSig )
cc      CALL GLK_SetNout(16)
cc      WRITE(6,*) '#####################################################'
*      
      RETURN
      WRITE(*,*) '+++++++++ STOP in MakeAng, nb1,nb2= ',nb1,nb2
 900  STOP
      END



      SUBROUTINE AngCut(idSig,idAfb)
*/////////////////////////////////////////////////////////////////////////////////
*//
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER          idSig,idAfb
      INTEGER          Nb1,Nb2,Nb,ib
      DOUBLE PRECISION SigTot,AfbTot,SigTotErr,AfbTotErr,Xs,Xe,Ac,Ae
      DOUBLE PRECISION Afb(100),AfbErr(100)
      DOUBLE PRECISION Sig(100),SigErr(100)
      DOUBLE PRECISION GLK_hi,GLK_hie
      LOGICAL GLK_Exist
*     -----------------------------------------------------------
      IF( (.NOT.GLK_Exist(idSig)) .OR. (.NOT.GLK_Exist(idSig))) THEN
         WRITE(*,*) ' ##### makeAng, histos do not exist idSig,idAfb=',idSig,idAfb
         STOP
      ENDIF
      CALL GLK_GetNb(idSig,Nb1)
      CALL GLK_GetNb(idAfb,Nb2)
      IF( NB1 .NE. NB2 )    GOTO 900
      NB=NB1
      SigTot    = GLK_hi(  idSig,NB)
      SigTotErr = GLK_hie( idSig,NB)**2
      AfbTot    = GLK_hi(  idAfb,NB)
      AfbTotErr = GLK_hie( idAfb,NB)**2
      DO ib=1,Nb
         Xs =  GLK_hi( idSig,ib)
         Xe =  GLK_hie(idSig,ib)**2
         Ac =  GLK_hi( idAfb,ib)
         Ae =  GLK_hie(idAfb,ib)**2
***         WRITE(*,*) '=== ib,Xs,Afb= ',ib,Xs,Xe,Ac,Ae
         Afb   (ib)= Ac-AfbTot
         AfbErr(ib)= SQRT(Ae+AfbTotErr)
         Sig   (ib)= (Xs-SigTot)/SigTot
         SigErr(ib)= SQRT(Xe+SigTotErr)/SigTot ! approximation
***         WRITE(*,*) '=== ib,Afb,Xse= ',ib,Afb(ib),AfbErr(ib), Sig(ib),SigErr(ib)
      ENDDO
      CALL GLK_Pak(  idSig,Sig)
      CALL GLK_Pake( idSig,SigErr)
      CALL GLK_Pak(  idAfb,Afb)
      CALL GLK_Pake( idAfb,AfbErr)
*
cc      WRITE(6,*) ' ################## this is AngCut ##################'
cc      CALL GLK_SetNout(6)
      CALL GLK_Print(idSig )
      CALL GLK_Print(idAfb )
cc      CALL GLK_SetNout(16)
cc      WRITE(6,*) ' ####################################################'
*      
      RETURN
      WRITE(*,*) '+++++++++ STOP in afb_int, nb1,nb2= ',nb1,nb2
 900  STOP
      END



      SUBROUTINE Plot_sigHO
*/////////////////////////////////////////////////////////////////////////////////
*//   gmake afb_int-sigHO.eps
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
      TeXfile   = 'afb_int-sigHO.txp'
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
*//   gmake afb_int-afbHO.eps
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
      TeXfile   = 'afb_int-afbHO.txp'
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

