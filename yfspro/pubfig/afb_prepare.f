*/////////////////////////////////////////////////////////////////////////////////
*//   alias kmake='make -f KKMakefile'
*//   dont forget about (cd ../; kmake makflag;)
*//   the order of executing is essential!
*//                     kmake clean   !!!! <-- do it when data changed !!!!
*//   first part;       kmake afb_sig-ps
*//   second part;      kmake afb_int-ps
*//   second part;      kmake afb_ang-ps
*//   Moreover
*//                     kmake afb_LepEWG-ps
*//                     kmake afb_sig2-ps  clone for 10GeV
*//                     kmake delta-ps
*//                     kmake bornex
*****************************************************************
*   kmake afb_int-ps
*   kmake afb_sig-ps
*   kmake afb_sig2-ps
*   kmake afb_ang-ps
*   kmake afb_LepEWG-ps
*   kmake bornex
* One by one, for example:
*     kmake afb.hst
*     kmake afb_sig-ps
*     kmake afb_sig2-ps  clone for 10GeV
*     kmake afb_int-ps
*     kmake afb_ang-ps
*     kmake afb_int-tab1.eps
*     kmake afb_int-Gsig.eps
* With older *.hst you may need to comment out *[[[[[ ....c]]]] also in afb-int.f
*----------------------------------------------------------------------------------
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      INCLUDE 'afb.h'
*
      CHARACTER*60      Tesnam, Dname
      CHARACTER*60      Hname, DumpFile
*--------------------------------------------------------------------
      ninp=  5
      nout= 16
 
      Tesnam    = 'afb_prepare'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

************************** 189 GeV **************************
      Hname  = '../E189GeV/pro.hst'            ! current
      Dname  = '../E189GeV/pro.input'          ! current
c
      Dname  = '../E189GeV/E189GeV_17M.input' ! 2013
      Hname  = '../E189GeV/E189GeV_17M.hst'   ! 2013
c
c      Dname  = '../E189GeV/E189GeV_8M.input' ! 2013
c      Hname  = '../E189GeV/E189GeV_8M.hst'   ! 2013
c
c      Dname  = '../E189GeV/E189GeV_PRD63.input' ! PRD63, 2000
c      Hname  = '../E189GeV/E189GeV_PRD63.hst'   ! PRD63, 2000
c
c      Dname  = '../E189GeV/E189GeV.input.ori'      ! latest PRD63
c      Hname  = '../E189GeV/E189GeV.hst'            ! latest PRD63
c      Hname  = '../E189GeV/E189GeV,theta1.hst' ! another theta!!!
c((((((
c      Dname  = '../E189GeV/QuarksFSRoffYR.input'
c      Hname  = '../E189GeV/QuarksFSRoffYR.hst.440M'
c      Dname  = '../E189GeV/InclusiveFSRoffYR.input'
c      Hname  = '../E189GeV/InclusiveFSRoffYR.hst.250M'
c      Dname  = '../E189GeV/MuFSRoff.input'
c      Hname  = '../E189GeV/MuFSRoff.hst'
c      Dname  = '../E189GeV/DownFSRoff.input'
c      Hname  = '../E189GeV/DownFSRoff.hst'
c      Dname  = '../E189GeV/UpFSRoff.input'
c      Hname  = '../E189GeV/UpFSRoff.hst'
c      Dname  = '../E189GeV/QuarksFSRoff.input'
c      Hname  = '../E189GeV/QuarksFSRoff.hst'
c      Dname  = '../E189GeV/QuarksFSRoffWT=1.input'
c      Hname  = '../E189GeV/QuarksFSRoffWT=1.hst'
c      Hname  = '../E189GeV/pro.hst'            ! current 
c      Hname  = '../E189GeV/QuarksFSRoffYR.hst'
c      Dname  = '../E189GeV/QuarksFSRoffWTed.input'
c      Hname  = '../E189GeV/QuarksFSRoffWTed.hst'
c))))))
c      Dname  = '../E189GeV/H189GeV.input'      ! hadrons
c      Hname  = '../E189GeV/H189GeV.hst'        ! hadrons
c      Hname  = '../E189GeV/E189GeV,new.hst'    ! july. 99, 21M
* <<<<<<<<<<<<<<<<<<<<<<< 91 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E91GeV/E91+0GeV.input'     ! latest
c      Hname  = '../E91GeV/E91+0GeV.hst'       ! latest
* <<<<<<<<<<<<<<<<<<<<<<< 120 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E120GeV/E120GeV.input'     ! may. 99, 22M
c      Hname  = '../E120GeV/E120GeV.hst'       ! may. 99, 22M
c      Hname  = '../E120GeV/E120GeV.hst.10M'   ! may. 99, 22M NEW!!!!
* <<<<<<<<<<<<<<<<<<<<<<< 500 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E500GeV/E500GeV.input'     ! Apr. 99, 5M
c      Hname  = '../E500GeV/E500GeV.hst'       ! Apr. 99, 5M
*********************** around resonance
c      Dname  = '../E91GeV/E91-3GeV.input'     ! jan. 99
c      Hname  = '../E91GeV/E91-3GeV.hst'       ! jan. 99
*
c      Dname  = '../E91GeV/E91-1.8GeV.input'   ! jan. 99
c      Hname  = '../E91GeV/E91-1.8GeV.hst'     ! jan. 99
*
c      Dname  = '../E91GeV/E91+1.8GeV.input'   ! jan. 99
c      Hname  = '../E91GeV/E91+1.8GeV.hst'     ! jan. 99
*
c      Dname  = '../E91GeV/E91+3GeV.input'     ! jan. 99
c      Hname  = '../E91GeV/E91+3GeV.hst'       ! jan. 99
**** specials ********
c      Dname  = '../E91GeV/E91GeV,Zonly.input'     ! feb. 99
c      Hname  = '../E91GeV/E91GeV,Zonly.hst'       ! feb. 99
c      Dname  = '../E91GeV/E91+0GeV,GamCon.input'   ! jan. 99
c      Hname  = '../E91GeV/E91+0GeV,GamCon.hst'     ! jan. 99
********
c      Dname  = '../E91GeV/U91+0GeV,noFSR.input'  ! agrees with KORZ
c      Hname  = '../E91GeV/U91+0GeV,noFSR.hst'    ! agrees with KORZ
c      Dname  = '../E91GeV/B91+0GeV,noFSR.input'  ! agrees with KORZ
c      Hname  = '../E91GeV/B91+0GeV,noFSR.hst'    ! agrees with KORZ
c      Dname  = '../E91GeV/U91+0GeV.input'   ! looks OK
c      Hname  = '../E91GeV/U91+0GeV.hst'     ! looks OK
c      Dname  = '../E91GeV/B91+0GeV.input'   ! looks OK
c      Hname  = '../E91GeV/B91+0GeV.hst'     ! looks OK
*
c      Dname  = '../E91GeV/H91-1.8GeV.input' ! 37M 23feb99
c      Hname  = '../E91GeV/H91-1.8GeV.hst'   ! 37M 23feb99
c      Dname  = '../E91GeV/H91+0GeV.input'   ! 37M 23feb99
c      Hname  = '../E91GeV/H91+0GeV.hst'     ! 37M 23feb99
c      Dname  = '../E91GeV/H91+1.8GeV.input' ! 46M 23feb99
c      Hname  = '../E91GeV/H91+1.8GeV.hst'   ! 46M 23feb99
****
c      Dname  = '../E91GeV/E91GeV.input'   !!! Actual
c      Hname  = '../E91GeV/pro.hst'        !!! Actual
***********************

*===================================================================================================
      CALL GLK_ReadFile(Hname)  ! Read histograms from MC run
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
* Write energy to special input tex file
****  CALL KK2f_GetCMSene( m_CMSene)                   ! get CMS energy
      CALL KKsem_GetCMSene(m_CMSene)
      OPEN( 23, file='Energy.tex')
      IF(m_CMSene.LT.91.187d0+10d0) THEN
         WRITE(23,'(a,f8.3,a)') '\\def\\Energy{',m_CMSene,'GeV} '
      ELSE
         WRITE(23,'(a,i4,a)')   '\\def\\Energy{',NINT(m_CMSene),'GeV} '
      ENDIF
      CLOSE(23)
* dump input file names for further use
      OPEN( 23, file='afb_prepare.files')
      WRITE(23,'(a)') Dname
      WRITE(23,'(a)') Hname
      CLOSE(23)
*===================================================================================================
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL Prepare
*=========================================================================
* Write all histograms into dump file, for control/debug
      DumpFile = './afb.hst'
      CALL GLK_WriteFile(DumpFile)
*=================================
      CLOSE(nout)
      END

      SUBROUTINE Prepare
*     ******************
      IMPLICIT NONE
      INCLUDE 'afb.h'
*------------
      CHARACTER*80 title
      CHARACTER*4  CHAK
      DOUBLE PRECISION       Born,BornPb, QCDcor
      DOUBLE PRECISION       xbin(400),xerr(400),ybin(400),yerr(400),xdel(400),ydel(400) !
      CHARACTER*80           DataFile1, DataFile2, DataFile3, DataFile4 
      CHARACTER*80           DataFile5, DataFile6, DataFile7, DataFile8, DataFile9 !
      CHARACTER*5            chType
      INTEGER                i,j,k, Nb, KeyFSR, KeyDis
      LOGICAL    GLK_Exist
*------------------------------------------------------------------------------
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,*) 'figchi:: BornPb= ',BornPb
      WRITE(nout,*) 'figchi:: BornPb= ',BornPb
*===================================================================
*               Monte Carlo
*===================================================================
      chak = "NB  "
*////////////////////////////////////////////////
*//  dSigma/dCosTheta for EEX and CEEX M.E.    //
*////////////////////////////////////////////////
      CALL GLK_RenHst(chak,IdGenYFS3, idc+4073, iangO2)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+3073, iangO2x)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+2073, iangO2xx)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+1073, iangO2xxx)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+4203, iangG2)     !     vBARE <0.90
      CALL GLK_RenHst(chak,IdGenYFS3, idc+3203, iangG2x)    ! (3) vBARE <0.19
      CALL GLK_RenHst(chak,IdGenYFS3, idc+2203, iangG2xx)   ! (4)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+1203, iangG2xxx)  ! Aleph v_p<0.2775
      CALL GLK_RenHst(chak,IdGenYFS3, idc+4253, iangN2)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+3253, iangN2x)    ! (5)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+2253, iangN2xx)
      CALL GLK_RenHst(chak,IdGenYFS3, idc+1253, iangN2xxx)
      CALL GLK_Delet(idc+4073)
      CALL GLK_Delet(idc+3073)
      CALL GLK_Delet(idc+2073)
      CALL GLK_Delet(idc+1073)
      CALL GLK_Delet(idc+4203)
      CALL GLK_Delet(idc+3203)
      CALL GLK_Delet(idc+2203)
      CALL GLK_Delet(idc+1203)
      CALL GLK_Delet(idc+4253)
      CALL GLK_Delet(idc+3253)
      CALL GLK_Delet(idc+2253)
      CALL GLK_Delet(idc+1253)
      CALL GLK_Delet(idc+4072)
      CALL GLK_Delet(idc+3072)
      CALL GLK_Delet(idc+2072)
      CALL GLK_Delet(idc+1072)
      CALL GLK_Delet(idc+4202)
      CALL GLK_Delet(idc+3202)
      CALL GLK_Delet(idc+2202)
      CALL GLK_Delet(idc+1202)
      CALL GLK_Delet(idc+4252)
      CALL GLK_Delet(idc+3252)
      CALL GLK_Delet(idc+2252)
      CALL GLK_Delet(idc+1252)
      CALL MakeAngInt(iangG2,   iangN2,    igAfbG2,   igSigG2)
      CALL MakeAngInt(iangG2x,  iangN2x,   igAfbG2x,  igSigG2x)  ! (3-5), BAREintON - BAREintOFF
      CALL MakeAngInt(iangG2xx, iangN2xx,  igAfbG2xx, igSigG2xx)
      CALL MakeAngInt(iangG2xxx,iangN2xxx, igAfbG2xxx,igSigG2xxx)! Aleph v_p<0.2775
* WG LEP2
      CALL GLK_RenHst(chak,IdGenYFS3, idc+5253, iangN2prop) ! (1)
      CALL GLK_Delet(idc+5253)
      CALL GLK_Delet(idc+5252)
      CALL MakeAngInt(iangN2x, iangN2prop, igAfb51, igSig51)     ! (5-1), BAREintOFF - PROPintOFF
      CALL MakeAngInt(iangG2x, iangN2prop, igAfb31, igSig31)     ! (3-1), BAREintON  - PROPintOFF
      CALL MakeAngInt(iangG2xx,iangG2x,    igAfb43, igSig43)     ! (4-3), BAREintON  - BAREintON
*////////////////////////////////////////////////
*//      CEEX  for sigma(vmax) total           //
*////////////////////////////////////////////////
      CALL GLK_CumHis(     IdGenYFS3, idb+  74, isigO3)       ! EEX3
      CALL GLK_CumHis(     IdGenYFS3, idb+ 201, isigG0)       ! CEEX0 full
      CALL GLK_CumHis(     IdGenYFS3, idb+ 251, isigG0nin)    ! CEEX0 intOFF
      CALL GLK_CumHis(     IdGenYFS3, idb+ 202, isigG1)       ! CEEX1 full
      CALL GLK_CumHis(     IdGenYFS3, idb+ 252, isigG1nin)    ! CEEX1 intOFF
      CALL GLK_CumHis(     IdGenYFS3, idb+ 203, isigG2)       ! CEEX2 full
      CALL GLK_CumHis(     IdGenYFS3, idb+ 253, isigG2nin)    ! CEEX2 intOFF
      CALL GLK_RenHst(chak,IdGenYFS3, idb+ 203, ksigG2)       ! per bins CEEX2
      CALL GLK_RenHst(chak,IdGenYFS3, idb+ 253, ksigG2nin)    ! per bins CEEX2 intOFF
      CALL GLK_CumHis(     IdGenYFS3, idb+1203, isigG2mG1)    ! CEEX O(alf2-alf1)         !!!new
      CALL GLK_CumHis(     IdGenYFS3, idb+1253, isigG2mG1nin) ! CEEX O(alf2-alf1) intOFF  !!!new
      CALL GLK_Delet(idb+ 74)
      CALL GLK_Delet(idb+201)
      CALL GLK_Delet(idb+202)
      CALL GLK_Delet(idb+251)
      CALL GLK_Delet(idb+252)
      CALL GLK_Delet(idb+203)
      CALL GLK_Delet(idb+253)
      CALL GLK_Delet(idb+1203)
      CALL GLK_Delet(idb+1253)
***      CALL GLK_Operat(isigG2,    '-',   isigG1,      isigG2mG1,    1d0, 1d0) ! directly
***      CALL GLK_Operat(isigG2nin, '-',   isigG1nin,   isigG2mG1nin, 1d0, 1d0) ! directly
      CALL GLK_Operat(isigG2mG1,   '/',isigG2, isigG2mG1,    1d0, 1d0) ! CEEX O(alf2-alf1) 
      CALL GLK_Operat(isigG2mG1nin,'/',isigG2, isigG2mG1nin, 1d0, 1d0) ! CEEX O(alf2-alf1) intOFF
c[[[
c      CALL GLK_SetNout(6)
c      CALL GLK_Print(isigG2mG1 )
c      CALL GLK_Print(isigG2mG1nin )
c      CALL GLK_SetNout(16)
c]]]
* SEMIREALISTIC
      CALL GLK_CumHis(     IdGenYFS3, ids+ 203, isigS2)       ! CEEX O(alf0)
      CALL GLK_CumHis(     IdGenYFS3, ids+ 253, isigS2nin)    ! CEEX O(alf0) int OFF
      CALL GLK_Print(isigS2)
      CALL GLK_Print(isigS2nin)
      CALL GLK_Delet(ids+202)
      CALL GLK_Delet(ids+252)
      CALL GLK_Delet(ids+203)
      CALL GLK_Delet(ids+253)
*////////////////////////////////////////////////
*//       CEEX         AFB(vmax)               //
*////////////////////////////////////////////////
      CALL GLK_RenHst(chak,IdGenYFS3, idb+2203, kafbG2)       ! CEEX O(alf2)
      CALL GLK_RenHst(chak,IdGenYFS3, idb+2253, kafbG2nin)    ! CEEX O(alf2) int OFF
*
      CALL GLK_CumHis(     IdGenYFS3, idb+2074, iafbO3)       ! EEX O(alf3)
      CALL GLK_CumHis(     IdGenYFS3, idb+2201, iafbG0)       ! CEEX O(alf0)
      CALL GLK_CumHis(     IdGenYFS3, idb+2251, iafbG0nin)    ! CEEX O(alf0) int OFF
      CALL GLK_CumHis(     IdGenYFS3, idb+2202, iafbG1)       ! CEEX O(alf1)
      CALL GLK_CumHis(     IdGenYFS3, idb+2252, iafbG1nin)    ! CEEX O(alf1) int OFF
      CALL GLK_CumHis(     IdGenYFS3, idb+2203, iafbG2)       ! CEEX O(alf2)
      CALL GLK_CumHis(     IdGenYFS3, idb+2253, iafbG2nin)    ! CEEX O(alf2) int OFF
      CALL GLK_CumHis(     IdGenYFS3, idb+3203, iafbG2mG1)    ! CEEX O(alf2-alf1)  !!!new
      CALL GLK_CumHis(     IdGenYFS3, idb+3253, iafbG2mG1nin) ! CEEX O(alf2-alf1) int OFF !!!new
*
      CALL GLK_Operat(iafbO3,   '/',isigO3,    iafbO3,    1d0, 1d0) ! EEX O(alf3)
      CALL GLK_Operat(iafbG0,   '/',isigG0,    iafbG0,    1d0, 1d0) ! CEEX O(alf0)
      CALL GLK_Operat(iafbG0nin,'/',isigG0nin, iafbG0nin, 1d0, 1d0) ! CEEX O(alf0) int OFF
      CALL GLK_Operat(iafbG1,   '/',isigG1,    iafbG1,    1d0, 1d0) ! CEEX O(alf1)
      CALL GLK_Operat(iafbG1nin,'/',isigG1nin, iafbG1nin, 1d0, 1d0) ! CEEX O(alf1) int OFF
      CALL GLK_Operat(iafbG2,   '/',isigG2,    iafbG2,    1d0, 1d0) ! CEEX O(alf2)
      CALL GLK_Operat(iafbG2nin,'/',isigG2nin, iafbG2nin, 1d0, 1d0) ! CEEX O(alf2) int OFF
*
      CALL GLK_CumHis(     IdGenYFS3, ids+2203, iafbS2)       ! CEEX O(alf1)
      CALL GLK_CumHis(     IdGenYFS3, ids+2253, iafbS2nin)    ! CEEX O(alf1) int OFF
      CALL GLK_Operat(iafbS2,   '/',isigS2,    iafbS2,    1d0, 1d0) ! CEEX O(alf1)
      CALL GLK_Operat(iafbS2nin,'/',isigS2nin, iafbS2nin, 1d0, 1d0) ! CEEX O(alf1) int OFF
      CALL GLK_Print(iafbS2)
      CALL GLK_Print(iafbS2nin)
      CALL GLK_Delet(idb+2203)
      CALL GLK_Delet(idb+2253)
      CALL GLK_Delet(idb+2074)
      CALL GLK_Delet(idb+2251)
      CALL GLK_Delet(idb+2202)
      CALL GLK_Delet(idb+2252)
      CALL GLK_Delet(idb+2201)
      CALL GLK_Delet(idb+3203) !!!new
      CALL GLK_Delet(idb+3253) !!!new
      CALL GLK_Delet(ids+2202)
      CALL GLK_Delet(ids+2252)
      CALL GLK_Delet(ids+2203)
      CALL GLK_Delet(ids+2253)
*////////////////////////////////////////////////
*//     CEEX      sigma(vmax) Interf.          //
*////////////////////////////////////////////////
      CALL GLK_CumHis(     IdGenYFS3, idb+ 281, isigG0int)    ! CEEX O(alf0)
      CALL GLK_CumHis(     IdGenYFS3, idb+ 282, isigG1int)    ! CEEX O(alf1)
      CALL GLK_CumHis(     IdGenYFS3, idb+ 283, isigG2int)    ! CEEX O(alf2)
      CALL GLK_CumHis(     IdGenYFS3, ids+ 283, isigS2int)    ! CEEX O(alf1)
      CALL GLK_Delet(idb+281)
      CALL GLK_Delet(idb+282)
      CALL GLK_Delet(idb+283)
c      CALL GLK_Delet(ids+281)
      CALL GLK_Delet(ids+282)
      CALL GLK_Delet(ids+283)
      CALL GLK_Operat(isigG0int,'/',isigG0,  isigG0int, 1d0, 1d0)
      CALL GLK_Operat(isigG1int,'/',isigG1,  isigG1int, 1d0, 1d0)
      CALL GLK_Operat(isigG2int,'/',isigG2,  isigG2int, 1d0, 1d0)
      CALL GLK_Operat(isigS2int,'/',isigS2,  isigS2int, 1d0, 1d0)
      CALL GLK_Print( isigG0int)
      CALL GLK_Print( isigG1int)
      CALL GLK_Print( isigS2int)
*////////////////////////////////////////////////
*//       CEEX        AFB(vmax) Interf.        //
*////////////////////////////////////////////////
      CALL GLK_CumHis(     IdGenYFS3, idb+2281, iafbG0int)           ! cumulate
      CALL GLK_Operat(iafbG0int, '/', isigG0,   iafbG0int, 1d0, 1d0) ! CEEX O(alf0)
      CALL GLK_Operat(isigG0int, '*', iafbG0,   idWork   , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbG0int, '-', idWork,   iafbG0int, 1d0, 1d0) ! CEEX O(alf1)

      CALL GLK_CumHis(     IdGenYFS3, idb+2282, iafbG1int)           ! cumulate
      CALL GLK_Operat(iafbG1int, '/', isigG1,   iafbG1int, 1d0, 1d0) ! CEEX O(alf1)approx.
      CALL GLK_Operat(isigG1int, '*', iafbG1,   idWork   , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbG1int, '-', idWork,   iafbG1int, 1d0, 1d0) ! CEEX O(alf1)

      CALL GLK_CumHis(     IdGenYFS3, idb+2283, iafbG2int)           ! cumulate
      CALL GLK_Operat(iafbG2int, '/', isigG2,   iafbG2int, 1d0, 1d0) ! CEEX O(alf1)approx.
      CALL GLK_Operat(isigG2int, '*', iafbG2,   idWork   , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbG2int, '-', idWork,   iafbG2int, 1d0, 1d0) ! CEEX O(alf1)

      CALL GLK_Delet(idWork)
      CALL GLK_CumHis(     IdGenYFS3, ids+2283, iafbS2int)           ! cumulate
      CALL GLK_Operat(iafbS2int, '/', isigS2,   iafbS2int, 1d0, 1d0) ! CEEX O(alf1)approx.
      CALL GLK_Operat(isigS2int, '*', iafbS2,   idWork   , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbS2int, '-', idWork,   iafbS2int, 1d0, 1d0) ! CEEX O(alf1)
*/////////////////////////////////////////////////////////////////////////////////////////
* direct way ----- in fact I do not that it is any worse????
*/////////////////////////////////////////////////////////////////////////////////////////
      CALL GLK_Operat(iafbG0,   '-',iafbG0nin,  iafbG0int, 1d0, 1d0) ! CEEX O(alf1)
      CALL GLK_Operat(iafbG1,   '-',iafbG1nin,  iafbG1int, 1d0, 1d0) ! CEEX O(alf1)
      CALL GLK_Operat(iafbS2,   '-',iafbS2nin,  iafbS2int, 1d0, 1d0) ! CEEX O(alf1)
      CALL GLK_Print( iafbG0int)
      CALL GLK_Print( iafbG1int)
      CALL GLK_Print( iafbS2int)

*////////////////////////////////////////////////
*//      CEEX         AFB(v) Interf.           //
*////////////////////////////////////////////////
      CALL GLK_RenHst(chak, IdGenYFS3, idb+2282, kafbG2int)        ! CEEX O(alf2)
      CALL GLK_Operat(kafbG2int, '/', ksigG2, kafbG2int, 1d0, 1d0) ! CEEX O(alf1)
*////////////////////////////////////////////////
*//     CEEX   O(alf2-alf1)                    //
*////////////////////////////////////////////////
* tricky way
      CALL GLK_Delet(idWork)
      CALL GLK_Operat(iafbG2mG1, '/', isigG1,   iafbG2mG1, 1d0, 1d0) ! CEEX O(alf0)
      CALL GLK_Operat(isigG2mG1, '*', iafbG1,   idWork   , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbG2mG1, '-', idWork,   iafbG2mG1, 1d0, 1d0) ! CEEX O(alf1)
* tricky way
      CALL GLK_Operat(iafbG2mG1nin, '/', isigG1nin,   iafbG2mG1nin, 1d0, 1d0) ! CEEX O(alf0)
      CALL GLK_Operat(isigG2mG1nin, '*', iafbG1nin,   idWork      , 1d0, 1d0) ! Correction
      CALL GLK_Operat(iafbG2mG1nin, '-', idWork,      iafbG2mG1nin, 1d0, 1d0) ! CEEX O(alf1)
*/////////////////////////////////////////////////////////////////////////////////////////
* direct way ---- OK!!!  here loose a lot without the trick
*/////////////////////////////////////////////////////////////////////////////////////////
c      CALL GLK_Operat(iafbG2,   '-',iafbG1,    iafbG2mG1,    1d0, 1d0) ! CEEX O(alf2-alf1)
c      CALL GLK_Operat(iafbG2nin,'-',iafbG1nin, iafbG2mG1nin, 1d0, 1d0) ! CEEX O(alf2-alf1) int OFF
*----------------------------------------
c[[[
c      CALL GLK_SetNout(6)
c      CALL GLK_Print(iafbG2mG1 )
c      CALL GLK_Print(iafbG2mG1nin )
c      CALL GLK_SetNout(16)
c]]]
*===================================================================
*            Input Dizet and KoralZ results
*===================================================================
      CALL GLK_Clone1(isigO3,iDizSig,   'Sigma Dizet interf. ON$')
      CALL GLK_Clone1(isigO3,iDizSigNin,'Sigma Dizet interf. OFF$')
      CALL GLK_Clone1(isigO3,iDizAfb,   'AFB   Dizet interf. ON$')
      CALL GLK_Clone1(isigO3,iDizAfbNin,'AFB   Dizet interf. OFF$')
      CALL GLK_Clone1(isigO3,iKorzSig,  'sig   KORALZ $')
      CALL GLK_Clone1(isigO3,iKorzAfb,  'sig   KORALZ $')
      CALL GLK_Clone1(isigO3,iMustSig,  'sig   KORALZ O(alf1) interf. ON$')
      CALL GLK_Clone1(isigO3,iMustAfb,  'AFB   KORALZ O(alf1) interf. ON$')
      CALL GLK_Clone1(isigO3,iMustSigInt,  'sig   KORALZ O(alf1) interf.$')
      CALL GLK_Clone1(isigO3,iMustAfbInt,  'AFB   KORALZ O(alf1) interf.$')

      CALL GLK_GetNb(iDizSig,Nb)
*---------------------------------------------------------
      DataFile1 = 'Zfitter189GeV'
      DataFile2 = 'Zfitter189GeVnoint'
      DataFile3 = 'Koralz189GeV'
      DataFile4 = 'Mustra189GeV'
      DataFile5 = 'MusAng189GeV'
      DataFile6 = 'MusAng189GeVnoint'
      DataFile7 = 'Kor402-189GeV' !! ficticious default ??
      DataFile8 = 'Kor404-189GeV' !! ficticious default ??
      DataFile9 = 'Zft620-189GeV' !! ficticious default ??
 !!!!!! ficticious data for 500GeV !!!!!!
      IF( ABS(m_CMSene-500d0) .LT. 1d-5 ) THEN
         DataFile1 = 'Zfitter189GeV'
         DataFile2 = 'Zfitter189GeVnoint'
         DataFile3 = 'Koralz189GeV'
         DataFile4 = 'Mustra189GeV'
         DataFile5 = 'MusAng189GeV'
         DataFile6 = 'MusAng189GeVnoint'
      ENDIF
      IF( ABS(m_CMSene-189d0) .LT. 1d-5 ) THEN
         DataFile1 = 'Zfitter189GeV'
         DataFile2 = 'Zfitter189GeVnoint'
         DataFile3 = 'Koralz189GeV'
         DataFile4 = 'Mustra189GeV'
         DataFile5 = 'MusAng189GeV'
         DataFile6 = 'MusAng189GeVnoint'
         DataFile7 = 'Kor402-189GeV'
         DataFile8 = 'Kor404-189GeV'
         DataFile9 = 'Zft620-189GeV'
      ENDIF
*---------------------------------------------------------
      IF( ABS(m_CMSene-120d0) .LT. 1d-5 ) THEN
         DataFile1 = 'Zfitter120GeV'
         DataFile2 = 'Zfitter120GeVnoint'
         DataFile3 = 'Koralz120GeV'
         DataFile4 = 'Mustra120GeV'
         DataFile5 = 'MusAng120GeV'
         DataFile6 = 'MusAng120GeVnoint'
      ENDIF
*---------------------------------------------------------
      IF( ABS(m_CMSene-91.187d0) .LT. 0.010d0 ) THEN
         DataFile1 = 'ZfitterMZ'
         DataFile2 = 'ZfitterMZnoint'
         DataFile3 = 'KoralzMZ'
         DataFile4 = 'MustraMZ'
         DataFile5 = 'MusAng120GeV'      ! ficticious!!!!
         DataFile6 = 'MusAng120GeVnoint' ! ficticious!!!!
      ENDIF
      IF( ABS(m_CMSene-92.987d0) .LT. 0.010d0 ) THEN
         DataFile1 = 'ZfitterMZ+1.8'
         DataFile2 = 'ZfitterMZ+1.8noint'
         DataFile3 = 'KoralzMZ+1.8'
         DataFile4 = 'MustraMZ+1.8'
      ENDIF
      IF( ABS(m_CMSene-89.387d0) .LT. 0.010d0 ) THEN
         DataFile1 = 'ZfitterMZ-1.8'
         DataFile2 = 'ZfitterMZ-1.8noint'
         DataFile3 = 'KoralzMZ-1.8'
         DataFile4 = 'MustraMZ-1.8'
      ENDIF
      IF( ABS(m_CMSene-94.187d0) .LT. 0.010d0 ) THEN
         DataFile1 = 'ZfitterMZ+3'
         DataFile2 = 'ZfitterMZ+3noint'
         DataFile3 = 'KoralzMZ+3'
         DataFile4 = 'MustraMZ+3'
      ENDIF
      IF( ABS(m_CMSene-88.187d0) .LT. 0.010d0 ) THEN
         DataFile1 = 'ZfitterMZ-3'
         DataFile2 = 'ZfitterMZ-3noint'
         DataFile3 = 'KoralzMZ-3'
         DataFile4 = 'MustraMZ-3'
      ENDIF
*///////////////////////////////////////////////////////////////
*//             Zfitter/TopaZ0                                //
*///////////////////////////////////////////////////////////////
      CALL RData1(DataFile1,Nb,xbin,xerr,ybin,yerr)
      CALL GLK_Pak(  iDizSig,xbin)
      CALL GLK_Pake( iDizSig,xerr)
      CALL GLK_Pak(  iDizAfb,ybin)
      CALL GLK_Pake( iDizAfb,yerr)
      CALL RData2(DataFile2,Nb,xbin,xerr,ybin,yerr)
      CALL GLK_Pak(  iDizSigNin,xbin)
      CALL GLK_Pake( iDizSigNin,xerr)
      CALL GLK_Pak(  iDizAfbNin,ybin)
      CALL GLK_Pake( iDizAfbNin,yerr)
      CALL GLK_Operat(iDizSig, '-',iDizSigNin, iDizSigInt, 1d0, 1d0) ! Dizet iterf.
      CALL GLK_Operat(iDizSigInt, '/',iDizSig, iDizSigInt, 1d0, 1d0)
      CALL GLK_Operat(iDizAfb, '-',iDizAfbNin, iDizAfbInt, 1d0, 1d0) ! Dizet iterf.
*///////////////////////////////////////////////////////////////
*//                KORALZ/YFS3                                //
*///////////////////////////////////////////////////////////////     
      CALL RData3(DataFile3,Nb,xbin,xerr,ybin,yerr)
      CALL GLK_Pak(  iKorzSig,xbin)
      CALL GLK_Pake( iKorzSig,xerr)
      CALL GLK_Pak(  iKorzAfb,ybin)
      CALL GLK_Pake( iKorzAfb,yerr)
      CALL GLK_Operat(iKorzSig, '+',iKorzSig, iKorzSig, 1d3, 0d0) ! nb-->pb
*///////////////////////////////////////////////////////////////
*//           KoralZ/Mustraal O(alf^1) results                //
*//  Absolute values are very much off-scale!!!               //
*///////////////////////////////////////////////////////////////
      CALL  RData4(DataFile4,Nb,xbin,xerr,ybin,yerr,xdel,ydel)
      CALL GLK_Pak(   iMustSig,xbin)    ! Sigma
      CALL GLK_Pake(  iMustSig,xerr)    ! Sigma
      CALL GLK_Pak(   iMustAfb,ybin)    ! Afb
      CALL GLK_Pake(  iMustAfb,yerr)    ! Afb
      CALL GLK_Pak(   iMustSigInt,xdel) ! delta-Sigma Interf. contr. alone
c      CALL GLK_Pake(  iMustSigInt,xerr) ! delta-Sigma Interf. contr. alone (xerr is dummy)
      CALL GLK_Pak(   iMustAfbInt,ydel) ! delta-Afb   Interf. contr. alone
c      CALL GLK_Pake(  iMustAfbInt,yerr) ! delta-Afb   Interf. contr. alone (yerr is dummy)
      CALL GLK_Operat(iMustSig,   '+',iMustSig,    iMustSig,    1d3, 0d0) ! [nb]-->[pb]
      CALL GLK_Operat(iMustSigInt,'+',iMustSigInt, iMustSigInt, 1d3, 0d0) ! [nb]-->[pb]
* Correction factor
      CALL GLK_Operat(iMustSig,   '-',iMustSigInt, iMustCor,    1d0, 1d0) ! sigMust without IFI
      CALL GLK_Operat(isigG2nin,  '/',iMustCor ,   iMustCor,    1d3, 1d0) ! sigMustInt renormalized
* Combined IFI and CEEX2
      IF( m_CMSene .LE. 200d0) THEN
* recrealte sig and Afb Without IFI
      CALL GLK_Operat(iMustSig,  '-',iMustSigInt , iMustSigNin,  1d0, 1d0) !
      CALL GLK_Operat(iMustAfb,  '-',iMustAfbInt , iMustAfbNin,  1d0, 1d0) !
* Forw/back IFI on
      CALL GLK_Operat(iMustSig,  '*',iMustAfb,  idWork,     1d0, 1d0) !
      CALL GLK_Operat(iMustSig,  '+',idWork,    iMustSigF,  0.5d0, 0.5d0) !
      CALL GLK_Operat(iMustSig,  '-',idWork,    iMustSigB,  0.5d0, 0.5d0) !
* Forw/back IFI off
      CALL GLK_Operat(iMustSigNin,  '*',iMustAfbNin,  idWork,     1d0, 1d0) !
      CALL GLK_Operat(iMustSigNin,  '+',idWork,    iMustSigNinF,  0.5d0, 0.5d0) !
      CALL GLK_Operat(iMustSigNin,  '-',idWork,    iMustSigNinB,  0.5d0, 0.5d0) !
* delta-Sigma Forw/Back due to IFI
      CALL GLK_Operat(iMustSigF,  '-',iMustSigNinF, iMustDelSigF, 1d0, 1d0) !
      CALL GLK_Operat(iMustSigB,  '-',iMustSigNinB, iMustDelSigB, 1d0, 1d0) !
* Correction factor Forw/Back due to IFI
      CALL GLK_Operat(iMustSigF, '/',iMustSigNinF, iMustRatSigF, 1d0, 1d0) !
      CALL GLK_Operat(iMustSigB, '/',iMustSigNinB, iMustRatSigB, 1d0, 1d0) !
* Forw/back CEEX2 IFI off
      CALL GLK_Operat(isigG2nin, '*',iafbG2nin,  idWork,     1d0, 1d0) !
      CALL GLK_Operat(isigG2nin,  '+',idWork,    isigG2ninF, 0.5d0, 0.5d0) !
      CALL GLK_Operat(isigG2nin,  '-',idWork,    isigG2ninB, 0.5d0, 0.5d0) !
* Hybrid1
      CALL GLK_Operat(isigG2ninF, '+',iMustDelSigF, iHyb1SigF ,  1d3, 1d0) ! F+delta_F
      CALL GLK_Operat(isigG2ninB, '+',iMustDelSigB, iHyb1SigB ,  1d3, 1d0) ! B+delta_B
      CALL GLK_Operat(iHyb1SigF,  '+',iHyb1SigB,    iHyb1Sig ,   1d0, 1d0) ! F+B
      CALL GLK_Operat(iHyb1SigF,  '-',iHyb1SigB,    iHyb1Afb ,   1d0, 1d0) ! F-B
      CALL GLK_Operat(iHyb1Afb,   '/',iHyb1Sig,     iHyb1Afb ,   1d0, 1d0) ! (F-B)/(F+B)
* Hybrid2
      CALL GLK_Operat(isigG2ninF, '*',iMustRatSigF, iHyb2SigF ,  1d3, 1d0) ! F+delta_F
      CALL GLK_Operat(isigG2ninB, '*',iMustRatSigB, iHyb2SigB ,  1d3, 1d0) ! B+delta_B
      CALL GLK_Operat(iHyb2SigF,  '+',iHyb2SigB,    iHyb2Sig ,   1d0, 1d0) ! F+B
      CALL GLK_Operat(iHyb2SigF,  '-',iHyb2SigB,    iHyb2Afb ,   1d0, 1d0) ! F-B
      CALL GLK_Operat(iHyb2Afb,   '/',iHyb2Sig,     iHyb2Afb ,   1d0, 1d0) ! (F-B)/(F+B)
*
      CALL GLK_Operat(isigG2nin,   '+',iMustSigInt ,  isigComb,  1d3, 1d0) ! Or Renormalized
      CALL GLK_Operat(iMustAFBInt, '+', iafbG2nin,    iafbComb,  1d0, 1d0) ! Direct O(alf1)
      CALL GLK_idopt( isigComb,'ERRO')
      CALL GLK_idopt( iafbComb,'ERRO')
      CALL GLK_Delet(idWork)
      ENDIF
* We may divide O(alf1) interf. by O(alf1) or O(alf3) sigma
cc      CALL GLK_Operat(iMustSigInt,'/',isigG1     , iMustSigInt, 1d0, 1d3) !
cc      CALL GLK_Operat(iMustSigInt,'/',iMustSig   , iMustSigInt, 1d0, 1d0)
cc      CALL GLK_Operat(iMustSigInt,'/',isigO3     , iMustSigInt, 1d0, 1d3)
      WRITE(nout,*) '[[[[[[[[[[[[[[[[[['
      CALL GLK_Print(iMustRatSigF)
      CALL GLK_Print(iMustRatSigB)
      CALL GLK_Print(iHyb2Afb)
      CALL GLK_Print(iHyb2Sig)
      WRITE(nout,*) ']]]]]]]]]]]]]]]]]]'
*///////////////////////////////////////////////////////////////
*//           KoralZ/Mustraal O(alf^1) results                //
*//               Angular distribution                        //
*///////////////////////////////////////////////////////////////
      Nb=40
      CALL  RData5(DataFile5,Nb,xbin,xerr) ! Int ON
      CALL  RData5(DataFile6,Nb,ybin,yerr) ! Int OFF
      CALL GLK_Book1(iMustAng,      'O(alf1) dSigma/dTheta intOFF$',Nb,-1d0,1d0)
      CALL GLK_Book1(iMustAngNint,  'O(alf1) dSigma/dTheta intON $',Nb,-1d0,1d0)
      CALL GLK_Pak(   iMustAng,     xbin) ! Interf. contr. alone
      CALL GLK_Pake(  iMustAng,     xerr) ! Interf. contr. alone
      CALL GLK_Pak(   iMustAngNint, ybin) ! Interf. contr. alone
      CALL GLK_Pake(  iMustAngNint, yerr) ! Interf. contr. alone
      CALL MakeAngInt(iMustAng, iMustAngNint,   jMustAfb,   jMustSig)
cc      CALL GLK_SetNout(6)
      CALL GLK_Print(iMustAng )
      CALL GLK_Print(iMustAngNint )
      CALL GLK_Print(jMustAfb )
      CALL GLK_Print(jMustSig )
cc      CALL GLK_SetNout(16)
*///////////////////////////////////////////////////////////////
*//           KoralZ 4.02  results                            //
*//           All flavours                                    //
*///////////////////////////////////////////////////////////////
      IF( GLK_Exist(idf) ) THEN
         CALL GLK_RenHst(chak,IdGenYFS3, idf, iCeex2) ! table KF versus cut
         CALL GLK_Operat(iCeex2,'+',iCeex2, iCeex2, 1d3, 0d0)
         CALL GLK_idopt(iCeex2,    'ERRO')
         CALL GLK_Print(iCeex2 )
*     
         CALL GLK_Clone1(iCeex2,iKKsem,      'sigma   KKsem $')
         CALL GLK_idopt( iKKsem,'ERRO')
         CALL GLK_Clone1(iCeex2,iZft620sig,  'sigma   Zfitter 6.2x $')
         CALL RData7(DataFile9,xbin,xerr,ybin,yerr)
         CALL GLK_Pak(  iZft620sig,xbin)
         CALL GLK_Pake( iZft620sig,xerr)
         CALL GLK_Print(iZft620sig )
         CALL GLK_idopt(iZft620sig,'ERRO')
*     
         CALL GLK_Clone1(iCeex2,iKor402sig,  'sigma   KORALZ 4.02 $')
         CALL RData7(DataFile7,xbin,xerr,ybin,yerr)
         CALL GLK_Pak(  iKor402sig,xbin)
         CALL GLK_Pake( iKor402sig,xerr)
         CALL GLK_Print(iKor402sig )
         CALL GLK_idopt(iKor402sig,'ERRO')
*     
         CALL GLK_Clone1(iCeex2,iKor404sig,  'sigma   KORALZ 4.04 $')
         CALL RData7(DataFile8,xbin,xerr,ybin,yerr)
         CALL GLK_Pak(  iKor404sig,xbin)
         CALL GLK_Pake( iKor404sig,xerr)
         CALL GLK_Print(iKor404sig )
         CALL GLK_idopt(iKor404sig,'ERRO')
      ENDIF
*=================================================================================================
*============================ Semi-analytical ====================================================
*=================================================================================================
      CALL KKsem_GetKeyFSR(KeyFSR)
      IF( KeyFSR .NE. 0) THEN
         KeyDis =  303302     ! O(alf2)LL
         KeyDis =  304302     ! ISR with O(alf3)LL exp(3/4*gam)
         KeyDis =  305302     ! The best ISR,  O(alf3)LL + O(alf2)NLL
         chType = "XCHI2"     ! ISR+FSR
      ELSE
         KeyDis=   305
         chType = "VCHI2"       ! FSR only
      ENDIF
      CALL KKsem_SetKeyFoB( 1)  ! forward
      CALL KKsem_VVplot(KeyDis,chType, iO3bestF,' O(alf3) ini+fin best Forward$' ,isigO3) !
      CALL KKsem_SetKeyFoB(-1)  ! backward
      CALL KKsem_VVplot(KeyDis,chType, iO3bestB,' O(alf3) ini+fin best Backward$',isigO3) !
      CALL KKsem_SetKeyFoB( 0)  ! back to normal
*[[[[[[ this is for Inclusive run
ccc      CALL Plot_RefFill(iKKsem,KeyDis,chType)
*]]]]]]
      CALL GLK_Operat(iO3bestF,  '+', iO3bestB ,    iO3best,     1d0, 1d0) ! F+B

      CALL GLK_Operat(iO3bestF,  '-', iO3bestB,     iO3bestAFB,  1d0, 1d0)
      CALL GLK_Operat(iO3bestAFB,'/', iO3best ,     iO3bestAFB,  1d0, 1d0)

      CALL GLK_Operat(iO3best,   '+', iO3best,      iO3best,     1d3, 0d0) ! --> picobarns
*(((((((((
*      CALL GLK_SetNout(6)
      CALL GLK_Print(iO3best)
      CALL GLK_Print(iafbO3)
      CALL GLK_Print(iO3bestAFB)
*      CALL GLK_SetNout(16)
*)))))))))
      CALL GLK_Delet(iO3bestF)
      CALL GLK_Delet(iO3bestB)
*=================================================================================================
*=================================================================================================
*=================================================================================================
      CALL KKsem_GetBorn(Born)  ! This re-sets QCDcor
      CALL BornV_GetQCDcor(QCDcor)
      WRITE(   *,'(a,f15.7)') 'AFB_prepare:: QCDcor= ',QCDcor
      WRITE(nout,'(a,f15.7)') 'AFB_prepare:: QCDcor= ',QCDcor
      CALL KKsem_SetKeyQCD(0d0)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f15.7)') 'afb_prepare:: Born [Pb], QCDoff= ',BornPb
      WRITE(nout,'(a,f15.7)') 'afb_prepare:: Born [Pb], QCDoff= ',BornPb
      END


      SUBROUTINE   RData7(DataFile1,xbin,xerr,ybin,yerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, KORALZ results         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*)
      INTEGER           Nb,k,i,npoint,nacc
      DOUBLE PRECISION  vmax,sig,dsig,afb,dafb
      CHARACTER*5       chDummy
      CHARACTER*5       Marker

      WRITE(*,*) '======== RData7 will read diskfile ',DataFile1
      Nb = nCut*nKF
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d9
         yerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      DO i=1,20
         READ(20,*) Marker
         IF(  Marker .EQ. '<get>' ) THEN
            READ(20,*) iKF
            READ(20,*) chDummy
            READ(20,*) chDummy
            READ(20,*) nPoint
            IF( (iKF.GT.nKF) .OR. (iKF.LE.0) ) GOTO 900
            IF( nPoint .GT. nCut) GOTO 900
            DO iCut=1,nPoint
               READ(20,*) vmax, nacc,sig,dsig, afb,dafb
               ibin = nCut*(iKF-1) +iCut
               IF(ibin.LT.1 .OR. ibin.GT.Nb) GOTO 900
               xbin(ibin) =  sig
               xerr(ibin) = dsig
               ybin(ibin) =  afb
               yerr(ibin) = dafb
               WRITE(*,*) vmax,sig,dsig,afb,dafb
            ENDDO
         ELSEIF( Marker .EQ. '<end>' ) THEN
            GOTO 200
         ELSE
            GOTO 900
         ENDIF
      ENDDO
 200  CONTINUE
      CLOSE(20)
      RETURN
 900  WRITE(*,*) '############# STOP in  RData7 !!!!'
      STOP
      END


      SUBROUTINE   RData1(DataFile1,Nb,xbin,xerr,ybin,yerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  from diskfile, error is zero                         //
*// data from Zfiter                                                        //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,afb
      CHARACTER*4       ch

      WRITE(*,*) '======== RData1 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d9
         yerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax,sig,afb
         ibin = vmax*100.0001d0
         IF(ibin.LT.1 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) = sig
         ybin(ibin) = afb
         WRITE(*,*) ibin,xbin(ibin),ybin(ibin)
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

      SUBROUTINE   RData2(DataFile1,Nb,xbin,xerr,ybin,yerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  from diskfile, error is zero                         //
*// data from Zfiter, clone of  RData1 olny different empty bins            //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,afb
      CHARACTER*4       ch

      WRITE(*,*) '======== RData2 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  -1d9
         xerr(  k) =  0d0
         ybin(  k) =  -1d19
         yerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax,sig,afb
         ibin = vmax*100.0001d0
         IF(ibin.LT.1 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) = sig
         ybin(ibin) = afb
         WRITE(*,*) ibin,xbin(ibin),ybin(ibin)
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

      SUBROUTINE   RData3(DataFile1,Nb,xbin,xerr,ybin,yerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, KORALZ results         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,dsig,afb,dafb
      CHARACTER*26      ch

      WRITE(*,*) '======== RData3 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d9
         yerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax, sig,dsig, afb,dafb
         ibin = vmax*100.0001d0
         IF(ibin.LT.1 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) =  sig
         xerr(ibin) = dsig
         ybin(ibin) =  afb
         yerr(ibin) = dafb
         WRITE(*,*) vmax,sig,dsig,afb,dafb
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

      SUBROUTINE   RData4(DataFile1,Nb,xbin,xerr,ybin,yerr,xdel,ydel)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and AFB from diskfile, data from KoralZ/Mustraal          //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*),xdel(*),ydel(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,dsig,afb,dafb,sigint,afbint
      CHARACTER*26      ch

      WRITE(*,*) '======== RData4 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d19
         yerr(  k) =  0d0
         xdel(  k) =  1d29
         ydel(  k) =  1d29
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax,  sig,sigint,dsig,  afb,afbint,dafb
         ibin = vmax*100.0001d0
         IF(ibin.LT.1 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) =  sig
         xerr(ibin) = dsig
         ybin(ibin) =  afb
         yerr(ibin) = dafb
         xdel(ibin) = sigint
         ydel(ibin) = afbint
         WRITE(*,*) vmax,sig,dsig,afb,dafb,sigint,afbint
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

      SUBROUTINE   RData5(DataFile1,Nb,xbin,xerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading angular distribution                                            //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*)
      INTEGER           Nb,k,i,ibin
      CHARACTER*26      ch
      DOUBLE PRECISION  x,sig,dsig
*
      WRITE(*,*) '======== RData5 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      DO ibin=1,Nb
         READ(20,*) x,sig,dsig
         xbin(ibin) =  sig
         xerr(ibin) = dsig
cc         WRITE(*,*) '== RData5 x,sig,dsig=', ibin, x,sig,dsig
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData5 !!!!'
      STOP
      END

      SUBROUTINE MakeAngInt(idAngI,idAngN,idAfbInt,idSigInt)
*/////////////////////////////////////////////////////////////////////////////////
*// Interference correction/contribution to sigma and Afb as a function of costheta_max
*// Input: two histograms of dsigma/dcostheta
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER          idAngI,idAngN,idAfbInt,idSigInt
      INTEGER          Nb1,Nb2,Nb,ib
      DOUBLE PRECISION XsForw,XeForw,Xsback,Xeback !! IntON
      DOUBLE PRECISION YsForw,YeForw,Ysback,Yeback !! IntOFF
      DOUBLE PRECISION AfbInt(100),AfbIntErr(100)
      DOUBLE PRECISION SigInt(100),SigIntErr(100)
      DOUBLE PRECISION GLK_hi,GLK_hie
      LOGICAL GLK_Exist
*     -----------------------------------------------------------
      IF( (.NOT.GLK_Exist(idAngI)) .OR. (.NOT.GLK_Exist(idAngI))) THEN
         WRITE(*,*) ' ##### MakeAngInt, histos do not exist idAngI,idAngN=',idAngI,idAngN
         STOP
      ENDIF
      CALL GLK_GetNb(idAngI,Nb1)
      CALL GLK_GetNb(idAngN,Nb2)
      IF( NB1 .NE. NB2 )    GOTO 900
      IF( MOD(NB1,2).NE.0 ) GOTO 900
      Nb=Nb1/2
c      idAfbInt = idAngI+1000000
c      idSigInt = idAngI+2000000
      CALL GLK_Book1(idAfbInt,  'AfbInt $',Nb,0d0,1d0)
      CALL GLK_Book1(idSigInt,  'SigInt $',Nb,0d0,1d0)
      XsForw =0d0
      XeForw =0d0
      Xsback =0d0
      Xeback =0d0
      YsForw =0d0
      YeForw =0d0
      Ysback =0d0
      Yeback =0d0
      DO ib=1,Nb
         XsForw =  XsForw+ GLK_hi( idAngI,NB+ib)
         XeForw =  XeForw+ GLK_hie(idAngI,NB+ib)**2
         Xsback =  Xsback+ GLK_hi( idAngI,NB-ib+1)
         Xeback =  Xeback+ GLK_hie(idAngI,NB-ib+1)**2
***         WRITE(*,*) '=== ib,XsForw,Xsback= ',ib,XsForw,Xsback
         YsForw =  YsForw+ GLK_hi( idAngN,NB+ib)
         YeForw =  YeForw+ GLK_hie(idAngN,NB+ib)**2
         Ysback =  Ysback+ GLK_hi( idAngN,NB-ib+1)
         Yeback =  Yeback+ GLK_hie(idAngN,NB-ib+1)**2
***         WRITE(*,*) '    ib,YsForw,Ysback= ',ib,YsForw,Ysback
***         AfbInt(ib)    = (XsForw-Xsback)/(XsForw+Xsback)          !! IntON
***         AfbIntErr(ib) = SQRT(XeForw + Xeback)/(XsForw+Xsback)
***         AfbInt(ib)    = (YsForw-Ysback)/(YsForw+Ysback)          !! IntOFF
***         AfbIntErr(ib) = SQRT(YeForw + Yeback)/(YsForw+Ysback)
***         AfbInt(ib)    = ((XsForw-Xsback)-(YsForw-Ysback))/(YsForw+Ysback)               !! Incorrect

         AfbInt(ib)    = (XsForw-Xsback)/(XsForw+Xsback)-(YsForw-Ysback)/(YsForw+Ysback) !! IntON-IntOFF
         AfbIntErr(ib) = SQRT(XeForw + Xeback +YeForw + Yeback)/(YsForw+Ysback)

         SigInt(ib)    = ((XsForw+Xsback)-(YsForw+Ysback))/(XsForw+Xsback)               !! (IntON-IntOFF)
         SigIntErr(ib) = SQRT(XeForw + Xeback +YeForw + Yeback)/(YsForw+Ysback)

      ENDDO
      CALL GLK_Pak(  idAfbInt,AfbInt)
      CALL GLK_Pake( idAfbInt,AfbIntErr)
      CALL GLK_Pak(  idSigInt,SigInt)
      CALL GLK_Pake( idSigInt,SigIntErr)
*
***      CALL GLK_SetNout(6)
***      CALL GLK_Print(idAfbInt )
***      CALL GLK_Print(idSigInt )
***      CALL GLK_SetNout(16)
*      
      RETURN
      WRITE(*,*) '+++++++++ STOP in afb-int, nb1,nb2= ',nb1,nb2
 900  STOP
      END

      SUBROUTINE Plot_RefFill(iHist,KeyDis,chType)
*//////////////////////////////////////////////////////////////////
* Fill reference spreadsheet with semianalytical results
*//////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'afb.h'
      INTEGER           iHist,KeyDis
      DOUBLE PRECISION  xbin(1000),xbin2(1000)
      INTEGER           i,j,k,idum
      CHARACTER*5       chType
!------------------------------------------------
      DO iKF =1,16
         IF(iKF.LE.5 .OR. iKF.EQ.13) THEN ! Quarks and muon
            CALL BornV_SetKF(iKF)
            CALL KKsem_VVplot(KeyDis, chType, iSemRef, 'O(alf3) ini+fin reference$',isigO3) !
            CALL GLK_UnPak(iSemRef,xbin,'    ',idum)
            DO iCut =1,nCut
               ibin =  nCut*(iKF-1)+iCut
               j= 10*(iCut-1)
               IF( iCut .EQ. 1  ) j=1
               IF( iCut .EQ. 11 ) j=99
               xbin2(ibin)= xbin(j)*1d3
ccc   write(*,*) ibin,j, xbin(j)
            ENDDO
         ENDIF
      ENDDO
      CALL GLK_Pak( iHist,xbin2)
      END

