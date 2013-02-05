*//////////////////////////////////////////////////////////////////////////
*//                                                                      //
*//                                                                      //
*//   make figmix-dvi
*//   make figmix-ps
*//                                                                      //
*//   MC results from yfspro/robol4                                      //
*//                                                                      //
*//                                                                      //
*//////////////////////////////////////////////////////////////////////////
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*---
      INTEGER   imax
      PARAMETER(imax=10000)
      REAL*8  xpar(imax)
*
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
      INTEGER lint
*---------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'figmix'
      TeXfile   = 'figmix.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

* PRINCIPAL  200GeV
c      Dname  = '../mix200/mix200.input.375M.xenph=1.25' ! Jan 99
c      Hname  = '../mix200/pro.hst.375M.xenph=1.25'      ! Jan 98
c      Dname  = '../mix200/mix200.input.320M' ! June 98
c      Hname  = '../mix200/pro.hst.320M'      ! June 98

      Dname  = '../mix200/mix200.input' ! current
      Hname  = '../mix200/pro.hst'      ! current

* GPS TESLA/NLC
c      Dname  = '../mix1000/mix1000.input.GPS'  ! current
c      Hname  = '../mix1000/pro.hst.GPS__6M'    ! current

* GPS CLEO
c      Dname  = '../mix10/mix10.input.GPS'  ! current
c      Hname  = '../mix10/pro.hst.GPS__7M'  ! current

* GPS Z-peak
cc      Dname  = '../mix200/mix200.input.GPS'
cc      Hname  = '../mix200/pro.hst.GPS_vmin-5__29M'   ! vmin=1d-5
cc      Hname  = '../mix200/pro.hst.GPS_vmin-9__8M'    ! vmin=1d-9

* GPS 200GeV
**      Dname  = '../mix200/mix200.input.GPS__38M' !!!!
**      Hname  = '../mix200/pro.hst.GPS__38M'      !!!! default vmin=1d-5

cc      Dname  = '../mix200/mix200.input.11M' ! Beam measuerement, for FigRap
cc      Hname  = '../mix200/pro.hst.11M'      ! Beam measuerement, for FigRap

* <<<<<<<<<<<<<<<<<<<<<<< 189 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E189GeV/E189GeV.input'
c      Hname  = '../E189GeV/E189GeV.hst'  ! last production
c      Hname  = '../E189GeV/pro.hst'       ! current
c      Hname  = '../E189GeV/E189GeV.new.hst' ! temporary
* <<<<<<<<<<<<<<<<<<<<<<< 120 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E120GeV/E120GeV.input.bullet__24M'
c      Hname  = '../E120GeV/E120GeV.hst.bullet__24M'
* <<<<<<<<<<<<<<<<<<<<<<< 91 GeV >>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Dname  = '../E91GeV/E91GeV.input.KeyZet=9__31M' ! jan. 99
c      Hname  = '../E91GeV/E91GeV.hst.KeyZet=9__31M'   ! jan. 99
c      Dname  = '../E91GeV/E91+1.8GeV.input'   ! jan. 99
c      Hname  = '../E91GeV/E91+1.8GeV.hst'     ! jan. 99
c      Dname  = '../E91GeV/E91-1.8GeV.input'   ! jan. 99
c      Hname  = '../E91GeV/E91-1.8GeV.hst'     ! jan. 99
c      Dname  = '../E91GeV/E91+3GeV.input'     ! jan. 99
c      Hname  = '../E91GeV/E91+3GeV.hst'       ! jan. 99
c      Dname  = '../E91GeV/E91-3GeV.input'     ! jan. 99
c      Hname  = '../E91GeV/E91-3GeV.hst'       ! jan. 99
c      Dname  = '../E91GeV/E91+0GeV.input'     ! last production
c      Hname  = '../E91GeV/E91+0GeV.hst'       ! last production
c      Hname  = '../E91GeV/pro.hst'            ! current

******************/////////////****************
c      Dname  = '../E91GeV/E91GeV.input' ! current
c      Hname  = '../E91GeV/pro.hst'      ! current

c      Dname  = '../mix91/mix91.input' ! current
c      Hname  = '../mix91/pro.hst'      ! current

*=====================================
* Read data, the same as in MC run
      CALL ReaDat(Dname,imax,xpar)
      CALL Semalib_Initialize(xpar)
* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
*=====================================
* Initialize GLK_Plot
      Lint=0
      CALL GLK_PlInitialize(Lint,TeXfile)
*=====================================================
c      CALL GPSfig               !! sig(vmax), AFB((vmax)
c      write(*,*) ' enter figAFB '
c      CALL figAFB               !! costheta
c      write(*,*) ' exit  figAFB '
*----------------------------------
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL figchi
*----------------------------------
* MC<-->SAN, O(alf2), sigma(xi_max), where xi=log(v/(1-v))
c      CALL figzet
*----------------------------------
* ISR+FSR<-->ISR, Pure M.C., sigma(vmax) and d(sigma)/(vmax)
c      CALL ini_fin
*----------------------------------
      CALL FigReconS
*----------------------------------
* end GLK_Plot, close LaTeX file
c      CALL GLK_PlEnd
*----------------------------------
* Write all histograms into dump file, for control/debug
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=====================================================
      CLOSE(nout)
      END


      SUBROUTINE GPSfig
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Study of ISR/FSR interference close to Z peak: sigma(vmax) AFB(vmax)          //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)

      COMMON / inout  / ninp,nout
      INTEGER           ninp,nout
      DOUBLE PRECISION Born,ymin,ymax,aa,bb
      INTEGER          idb,idz

      idb = 50000
      idz = 70000
      CALL Semalib_GetBorn(Born)

      aa = 0.010                ! 91GeV
      bb = 0.010                ! 91GeV
      aa = 0.050                ! 200GeV
      bb = 0.050                ! 200GeV

      CALL GLK_CumHis(     IdGenYFS3, idb+  74, idb+9074)
      CALL GLK_CumHis(     IdGenYFS3, idb+ 271, idb+9271)
      CALL GLK_CumHis(     IdGenYFS3, idb+ 272, idb+9272)
      CALL GLK_CumHis(     IdGenYFS3, idb+ 273, idb+9273)
      CALL GLK_Operat(idb+9271,  '/', idb+9074,   idb+9271,    1d0, 1d0)
      CALL GLK_Operat(idb+9272,  '/', idb+9074,   idb+9272,    1d0, 1d0)
      CALL GLK_Operat(idb+9273,  '/', idb+9074,   idb+9273,    1d0, 1d0)
      CALL plt1mc('    ','R(vmax)=(CEEX0-EEX0)/EEX3 IntOff $',idb+9271, -aa, bb,1d0)
      CALL plt1mc('    ','R(vmax)=(CEEX1-EEX1)/EEX3 IntOff $',idb+9272, -aa, bb,1d0)
      CALL plt1mc('    ','R(vmax)=(CEEX2-EEX2)/EEX3 IntOff $',idb+9273, -aa, bb,1d0)

      CALL plt1mc('CUMU','R(vmax)=(CEEX0-EEX0)/born: IntOff $',idb+ 271, -aa*2d0, bb*2d0,born)
      CALL plt1mc('CUMU','R(vmax)=(CEEX1-EEX1)/born: IntOff $',idb+ 272, -aa*2d0, bb*2d0,born)
      CALL plt1mc('CUMU','R(vmax)=(CEEX2-EEX2)/born: IntOff $',idb+ 273, -aa*2d0, bb*2d0,born)

      CALL plt1mc('CUMU','R(vmax)=(CEEX1-CEEX0)/born: IntOn $',idb+1202, -aa*30d0, bb*30d0,born)
      CALL plt1mc('CUMU','R(vmax)=(CEEX1-CEEX0)/born: NoInt $',idb+1252, -aa*30d0, bb*30d0,born)

      CALL plt1mc('CUMU','R(vmax)=(CEEX2-CEEX1)/born: IntOn $',idb+1203, -0.3d0*aa, 0.3d0*bb,born)
      CALL plt1mc('CUMU','R(vmax)=(CEEX2-CEEX1)/born: NoInt $',idb+1253, -0.3d0*aa, 0.3d0*bb,born)
      RETURN

      ymax =  0.010d0
      ymin = -0.010d0
      CALL plt1mc('CUMU','R(vmax)=CEEX0/born, Interf. ON--OFF $',idb+281, ymin, ymax, born)
      CALL plt1mc('CUMU','R(vmax)=CEEX1/born, Interf. ON--OFF $',idb+282, ymin, ymax, born)
      CALL plt1mc('CUMU','R(vmax)=CEEX2/born, Interf. ON--OFF $',idb+283, ymin, ymax, born)

      CALL plt1mc('CUMU','R(vmax)=CEEX0/born, Interf. ON--OFF $',idz+281, ymin, ymax, born)
      CALL plt1mc('CUMU','R(vmax)=CEEX1/born, Interf. ON--OFF $',idz+282, ymin, ymax, born)
      CALL plt1mc('CUMU','R(vmax)=CEEX2/born, Interf. ON--OFF $',idz+282, ymin, ymax, born)

      ymax =  0.030d0
      ymin = -0.030d0
      CALL plotAFB('CUMU','AFB(vmax): CEEX0, Interf. ON--OFF $', idb+201, idb+2281,ymin,ymax)
      CALL plotAFB('CUMU','AFB(vmax): CEEX1, Interf. ON--OFF $', idb+202, idb+2282,ymin,ymax)
      CALL plotAFB('CUMU','AFB(vmax): CEEX2, Interf. ON--OFF $', idb+203, idb+2282,ymin,ymax)

      ymax =  0.150d0
      ymin = -0.150d0
      CALL plotAFB('CUMU','AFB(zmax): CEEX0, Interf. ON--OFF $', idz+201, idz+2281,ymin,ymax)
      CALL plotAFB('CUMU','AFB(zmax): CEEX1, Interf. ON--OFF $', idz+202, idz+2282,ymin,ymax)
      CALL plotAFB('CUMU','AFB(zmax): CEEX2, Interf. ON--OFF $', idz+203, idz+2282,ymin,ymax)

      ymax =  0.150d0
      ymin = -0.150d0
      CALL plotAFB('CUMU','AFB(vmax): CEEX1, Interf. ON  $', idb+202, idb+2202,ymin,ymax)
      CALL plotAFB('CUMU','AFB(vmax): CEEX1, Interf. OFF $', idb+252, idb+2252,ymin,ymax)

      ymax =  0.300d0
      ymin = -0.300d0
      CALL plotAFB('CUMU','AFB(zmax): CEEX1, Interf. ON  $', idz+202, idz+2202,ymin,ymax)
      CALL plotAFB('CUMU','AFB(zmax): CEEX1, Interf. OFF $', idz+252, idz+2252,ymin,ymax)
      END


 

      SUBROUTINE figAFB
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   CosTheta distributions and                                                    //
*//   Forward - backward charge asymmetries                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      CHARACTER*4 CHAK
      SAVE
*----------------------
      CALL Semalib_GetBorn(Born)

*/////////////////////////////////////////////////////////////////////////////////////
*//   CosTheta distributions for several cuts                                       //
*/////////////////////////////////////////////////////////////////////////////////////
      idc  = 10000
      chak = 'NB  '

      aa = 1.4d0*Born
      CALL plotANG2(chak,'CosTh: GPS/OLD, O(alf1), vmax=0.2 $',idc+3072,idc+3202, 0d0,aa)

      chak = 'UNIT'
      ymin = 0.0d0
      ymax = 1.7d0
c      CALL plotANG1(chak,'CosTheta: OLD, O(alf1), vmax=1d-4$', idc+1072, ymin,ymax)
c      CALL plotANG1(chak,'CosTheta: GPS, O(alf1), vmax=1d-4$', idc+1202, ymin,ymax)
c      CALL plotANG1(chak,'CosTheta: OLD, O(alf1), vmax=1d-2$', idc+2072, ymin,ymax)
c      CALL plotANG1(chak,'CosTheta: GPS, O(alf1), vmax=1d-2$', idc+2202, ymin,ymax)
      CALL plotANG1(chak,'CosTheta: OLD, O(alf1), vmax=0.2 $', idc+3072, ymin,ymax)
      CALL plotANG1(chak,'CosTheta: GPS, O(alf1), vmax=0.2 $', idc+3202, ymin,ymax)
c      CALL plotANG1(chak,'CosTheta: OLD, O(alf1), vmax=0.9 $', idc+4072, ymin,ymax)
c      CALL plotANG1(chak,'CosTheta: GPS, O(alf1), vmax=0.9 $', idc+4202, ymin,ymax)

*/////////////////////////////////////////////////////////////////////////////////////
*//   AFB depenedence on vmax                                                       //
*/////////////////////////////////////////////////////////////////////////////////////
      CALL Semalib_GetCMSene(CMSene)
      idb = 50000
      aa = 1.00d0
      bb = 0.10d0
      IF( ABS(CMSene-91d0) .LT. 2.d0 ) THEN
         aa =  4.00d0
         bb = 0.025d0
      ENDIF
      chak = 'CUMU'
***      chak = '    '     !!!<-- asymetry bin-per-bin
*---GPS: O(alf0)
      CALL plotAFB(chak,'AFB(vmax): CEEX0), Interf. ON  $', idb+201, idb+2201, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): CEEX0), Interf. OFF $', idb+251, idb+2251, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): CEEX0),ON minus OFF $', idb+201, idb+2281, -bb, bb)
*---GPS: O(alf1)
      CALL plotAFB(chak,'AFB(vmax): CEEX1), Interf. ON  $', idb+202, idb+2202, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): CEEX1), Interf. OFF $', idb+252, idb+2252, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): CEEX1),ON minus OFF $', idb+202, idb+2282, -bb, bb)
*---OLD: O(alf1) and compare with GPS
      CALL plotAFB(chak,'AFB(vmax): OLD O(alf0), Interf. OFF $', idb+ 71, idb+2071, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): OLD O(alf1), Interf. OFF $', idb+ 72, idb+2072, -aa, aa)
      aa = 0.02
      CALL plotAFB(chak,'AFB(vmax): O(alf0), GPS--OLD, NoInt $', idb+ 71, idb+2271, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): O(alf1), GPS--OLD, NoInt $', idb+ 72, idb+2272, -aa, aa)
      aa = 0.25
      CALL plotAFB(chak,'AFB(vmax): GPS  O(alf1-alf0), NoInt $', idb+252, idb+3252, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): GPS  O(alf1-alf0), IntOn $', idb+202, idb+3202, -aa, aa)
c[[[
      aa = 1.30
      bb = 0.15
      CALL plt1mc('CUMU','sig(vmax)/born: GPS  O(alf1-alf0), IntOn$',idb+1202, -aa, bb,born)
      CALL plt1mc('CUMU','sig(vmax)/born: GPS  O(alf1-alf0), NoInt$',idb+1252, -aa, bb,born)
c]]] 
*/////////////////////////////////////////////////////////////////////////////////////
*//   AFB dependence on zmax, z=ln(z/(1-z))                                         //
*/////////////////////////////////////////////////////////////////////////////////////
      idz = 70000
*---GPS: O(alf0)
      CALL plotAFB(chak,'AFB(zmax): CEEX0), Interf. ON  $', idz+201, idz+2201, -aa, aa)
      CALL plotAFB(chak,'AFB(zmax): CEEX0), Interf. OFF $', idz+251, idz+2251, -aa, aa)
      CALL plotAFB(chak,'AFB(zmax): CEEX0), ON minus OFF$', idz+201, idz+2281, -bb, bb)
*---GPS: O(alf1)
      CALL plotAFB(chak,'AFB(zmax): CEEX1), Interf. ON  $', idz+202, idz+2202, -aa, aa)
      CALL plotAFB(chak,'AFB(zmax): CEEX1), Interf. OFF $', idz+252, idz+2252, -aa, aa)
      CALL plotAFB(chak,'AFB(zmax): CEEX1), ON minus OFF$', idz+202, idz+2282, -bb, bb)
*---OLD: O(alf1) and compare with GPS
      CALL plotAFB(chak,'AFB(zmax): OLD O(alf0), Interf. OFF $', idz+ 71, idz+2071, -aa, aa)
      CALL plotAFB(chak,'AFB(zmax): OLD O(alf1), Interf. OFF $', idz+ 72, idz+2072, -aa, aa)
      aa = 0.02
      CALL plotAFB(chak,'AFB(vmax): O(alf0), GPS--OLD, NoInt $', idz+ 71, idz+2271, -aa, aa)
      CALL plotAFB(chak,'AFB(vmax): O(alf1), GPS--OLD, NoInt $', idz+ 72, idz+2272, -aa, aa)
*--------------
      END

      SUBROUTINE figchi
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*// MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      CHARACTER*4 CHAK
      SAVE
*------------------------------------
      WRITE(NOUT,*) ' ====================================='
      WRITE(NOUT,*) ' ==========    figchi   =============='
      WRITE(NOUT,*) ' ====================================='

      CALL Semalib_GetBorn(Born)

      WRITE(   *,*) 'figchi:: Born= ',born
      WRITE(nout,*) 'figchi:: Born= ',born
      ida = 30000
      idb = 50000

*===================================================================
*===================================================================
*===================================================================
      KeyDis=     400300
      iNLLPhsP = 7000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'XCHI2', iNLLPhsP,' O(alf1) ini+fin $',idb+300)
      CALL GLK_Print(iNLLPhsP)
*
      KeyDis=     300300
      iO0LL = 7000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'XCHI2', iO0LL,' O(alf1) ini+fin $',idb+300)
      CALL GLK_Print(iO0LL)
*
      KeyDis=     301301
      iO1LL = 7000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'XCHI2', iO1LL,' O(alf1) ini+fin $',idb+300)
      CALL GLK_Print(iO1LL)
*
      KeyDis=     302302
      iO2LL = 7000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'XCHI2', iO2LL,' O(alf2) ini+fin $',idb+300)
      CALL GLK_Print(iO2LL)
*
      KeyDis=     303303        ! ISR with O(alf3)LL
      KeyDis=     304302        ! ISR with O(alf3)LL exp(3/4*gam)
      KeyDis=     305302        ! The best ISR,  O(alf3)LL + O(alf2)LL
      iO3LL = 7000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'XCHI2', iO3LL,' O(alf2) ini+fin $',idb+300)
      CALL GLK_Print(iO3LL)
*
* O(alf2)-O(alf1)
      iO2vO1 = 7306306
      CALL GLK_Operat(iO2LL,'-',iO1LL,iO2vO1,1d0,1d0)
      CALL GLK_Print( iO2vO1)
* O(alf3)-O(alf2)
      iO3vO2 = iO3LL +500000
      CALL GLK_Operat(iO3LL,'-',iO2LL,iO3vO2,1d0,1d0)
*===================================================================
      ymax = 0.090d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01

      CALL pldiff('CUMU','[400*300] Technical test!!! sig(xmax)/Born ,YMAG=.01$',
     $     idb+71, iNLLPhsP,   ymin,ymax, xmag,ymag, born)
      CALL pldiff('CUMU','[300*300] O(alf0), sig(xmax)/Born  ,YMAG=.01$',
     $     idb+71, iO0LL,   ymin,ymax, xmag,ymag, born)
**********GPS**********(((((((((((((((((((((
      CALL pldiff('CUMU','GPS: O(alf0) vs [300*300], Interf. ON  ,YMAG=.01$',
     $     idb+201, iO0LL,   ymin,ymax, xmag,ymag, born)
      CALL pldiff('CUMU','GPS: O(alf0) vs [300*300], Interf. OFF ,YMAG=.01$',
     $     idb+251, iO0LL,   ymin,ymax, xmag,ymag, born)
      CALL plt1mc('CUMU','GPS: O(alf0), Interf. ON minus OFF $',
     $     idb+281, ymin, ymax, born)
***********************)))))))))))))))))))))
      CALL pldiff('CUMU','O(alf3)MC - O(alf3)Best , sig(xmax)/Born  ,YMAG=.01$',
     $     idb+74, iO3LL,   ymin,ymax, xmag,ymag, born)
      CALL pldiff('CUMU','O(alf2)MC - O(alf2)SAN , sig(xmax)/Born  ,YMAG=.01$',
     $     idb+73, iO2LL,   ymin,ymax, xmag,ymag, born)
      CALL pldiff('CUMU','O(alf1)MC - O(alf1)SAN, sig(xmax)/Born  ,YMAG=.01$',
     $     idb+72, iO1LL,   ymin,ymax, xmag,ymag, born)
**********GPS**********(((((((((((((((((((((
      CALL pldiff('CUMU','GPS: O(alf1) vs [301*301], Interf. ON   ,YMAG=.01$',
     $     idb+202, iO1LL,   ymin,ymax, xmag,ymag, born)
      CALL pldiff('CUMU','GPS: O(alf1) vs [301*301], Interf. OFF   ,YMAG=.01$',
     $     idb+252, iO1LL,   ymin,ymax, xmag,ymag, born)
      CALL plt1mc('CUMU','GPS: O(alf1), Interf. ON minus OFF !!! $',
     $     idb+282, ymin, ymax, born)
***********************)))))))))))))))))))))
*===================================================================
      CALL Semalib_GetCMSene(CMSene)
      ymax =  0.006d0
      IF( cmsene .GT. 150 )  ymax = 0.10d0
      ymin = -ymax
      xmag =   1
      ymag =   1
      CALL pldiff('CUMU','O(alf2)-O(alf1) total ,YMAG=1$',
     $     idb+76, iO2vO1,   ymin,ymax, xmag,ymag, born)
*---
c      ymax=0.10d0
c      ymin = -ymax
c      xmag =   1
c      ymag =   0.05
c      CALL pld2mc('CUMU','O(alf2) vs. O(alf1), sig(vmax)/Born ,YMAG=0.01$',
c     $     idb+73, idb+72,   ymin,ymax, xmag,ymag, born)
*===================================================================
      ymax=   0.02d0
      ymin = -ymax
      xmag =   1
      ymag =   0.005
c      CALL pld2mc('CUMU','O(alf3) vs. O(alf2), sig(vmax)/Born ,YMAG=0.005$',
c     $     idb+74, idb+73,   ymin,ymax, xmag,ymag, born)
      CALL plt1mc('CUMU','O(alf3)-O(alf2) MONTE CARLO, sig(vmax)/Born $',
     $     idb+77,ymin,ymax,born)
*-------
      kO3vO2 = iO3vO2 +10000
      CALL plt1an('O(alf3)LL-O(alf2) PURE SAN, sig(vmax)/Born $',
     $     iO3vO2,kO3vO2,ymin,ymax,born)
*===================================================================
      CALL GLK_Delet(iO1LL)
      CALL GLK_Delet(iO2LL)
      CALL GLK_Delet(iO2vO1)
      CALL GLK_Delet(iO3vO2)
      CALL GLK_Delet(kO3vO2)

      END

      SUBROUTINE figzet
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//  MC<-->SAN, O(alf2), sigma(xi_max), where xi=log(v/(1-v))                       //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE
      CHARACTER*4 CHAK
*-------------------------------------------------------------------

      CALL Semalib_GetBorn(Born)
      idz = 70000

*===================================================================
* O(alf2) ini+fin
      KeyDis=     400300
      jNLLPhsP = 9000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'ZCHI2', jNLLPhsP,' O(alf0) ini+fin$',idz+73)
      CALL GLK_Print(jNLLPhsP)

* O(alf2) ini+fin
      KeyDis=     302302
      jO2LL = 9000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'ZCHI2', jO2LL,' O(alf2) ini+fin$',idz+73)
      CALL GLK_Print(jO2LL)

* O(alf3)ini*O(alf2)fib
      KeyDis=     305302  ! The best ISR,  O(alf3)NLL * O(alf2)LL
      jO3NLL = 9000000 +KeyDis
      CALL Semalib_VVplot(KeyDis,'ZCHI2', jO3NLL,' O(alf3) ini+fin$',idz+73)
      CALL GLK_Print(jO3NLL)
*===================================================================

* O(alf2) MC vs. SAN, sig(zmax)/Born
      ymax =   0.06              ! 200GeV
      ymin =  -0.06
      xmag =  1
      ymag =  .01

c      iNLLPhsP = 7000000 +400300
c      CALL pld2an('[400*300] test of GAUSS integration, YMAG=0.01$',
c     $     iNLLPhsP,jNLLPhsP,ymin,ymax,xmag,ymag,born)

      CALL pldiff('CUMU','[400*300] Technical test!!! sig(zmax)/Born,YMAG=.01$',
     $     idz+71, jNLLPhsP,   ymin,ymax, xmag,ymag, born)

      CALL pldiff('CUMU','O(alf2)MC - [302302], sig(zmax)/Born ,YMAG=.01$',
     $     idz+73, jO2LL,   ymin,ymax, xmag,ymag, born)

      CALL pldiff('CUMU','O(alf3)MC - [305302]Best!, sig(zmax)/Born ,YMAG=.01$',
     $     idz+74, jO3NLL,   ymin,ymax, xmag,ymag, born)

*--------------------------------------------------------------------
* O(alf2) vs. O(alf1)
      ymax =  0.20              ! 200GeV
      ymin = -0.5*ymax
      xmag =  1
      ymag =  .06
      CALL pld2mc('CUMU','O(alf2) vs. O(alf1), sig(zmax)/Born ,YMAG=.02$',
     $     idz+73, idz+72,   ymin,ymax, xmag,ymag, born)

      ymax =  0.08              ! 200GeV
      ymin = -0.5*ymax
      xmag =  1
      ymag =  .02
      CALL pld2mc('CUMU','O(alf3) vs. O(alf2), sig(zmax)/Born ,YMAG=.02$',
     $     idz+74, idz+73,   ymin,ymax, xmag,ymag, born)

      END


      SUBROUTINE ini_fin
*     ***********************************
* ISR+FSR <--> ISR, Pure M.C.
* sigma(vmax)     and sigma(log10(vmax)), 
* d(sigma)/(vmax) and d(sigma)/d(log10(vmax))
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE
*---------------------------------------------
      CALL Semalib_GetBorn(Born)
      ida = 30000
      idb = 50000
*===================================================================
* dsig/dlog(x)/Born
      CALL Semalib_GetCMSene(CMSene)
      ymax =  .4                       ! 100GeV
      IF( CMSene .GT. 150 ) ymax =  25 ! 200GeV
      ymin = -.25*ymax
      xmag =  1
      ymag =  1
      CALL pld2mc('NB10','O(alf2) init+fin and init, dsig/d(x)/Born  ,YMAG=1$',
     $     idb+302, idb+73,   ymin,ymax, xmag,ymag, born)

* Cumulative dsig(log(xmax))/Born
      ymax =  1.5                       ! 100GeV
      IF( CMSene .GT. 150 ) ymax =  3.0 ! 200GeV
      ymin = -0.25*ymax
      xmag =  1
      ymag =  1
      CALL pld2mc('CUMU','O(alf2) init+fin and init, sig(xmax)/Born ,YMAG=1$',
     $     idb+302, idb+73,   ymin,ymax, xmag,ymag, born)

      END




      SUBROUTINE FigReconS
*///////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout

      SAVE

      WRITE(NOUT,*) ' ====================================='
      WRITE(NOUT,*) ' ==========    FigReconS   ==========='
      WRITE(NOUT,*) ' ====================================='
      CALL Semalib_GetBorn(Born)

      WRITE(   *,*) 'figchi:: Born= ',born
      WRITE(nout,*) 'figchi:: Born= ',born
      ida = 30000
      ids = 40000
      idb = 50000

      CALL GLK_Print(ids +74)
      CALL GLK_Print(ids+203)
* O(alf2-alf1)
      CALL GLK_Operat(ids+90,'/',ids+203,ids+90, 1d0, 1d0)
      CALL GLK_Print( ids+90)
* IFIon-IFIoff
      CALL GLK_Operat(ids+283,'/',ids+203,ids+283, 1d0, 1d0)
      CALL GLK_Print( ids+283)
      CALL plt1mc('NB  ', 'O(alf3) true dsigma/dv [nb] $',ids +74,     0d0,0.150d0,1d0)
      CALL plt1mc('    ', 'O(alf2)-O(alf1)             $',ids +90,-0.100d0,0.100d0,1d0)
      CALL plt1mc('    ', 'ISR*FSR                     $',ids+283,-0.100d0,0.100d0,1d0)
*     ******
      RETURN
*     ******

      CALL GLK_Print(ida+73)
      CALL GLK_Print(idb+73)
c[      CALL GLK_Print(ids+ 202)

      ymin = 0d0
      ymax =0.150d0
      CALL plt1mc('NB  ', 'O(alf2) true v [nb]         $',idb+73 ,ymin,ymax,1d0)
      CALL plt1mc('NB  ', 'O(alf2)  reconstructed v [nb]$',ids+73 ,ymin,ymax,1d0)
      CALL plt1mc('NB  ', 'CEEX1) reconstruc v [nb]$',ids+2202,ymin,ymax,1d0)
c[      CALL plt1mc('NB  ', 'CEEX1) reconstruc v [nb]$',ids+202,ymin,ymax,1d0)
      ymax =0.020d0
      CALL plt1mc('NB  ', 'O(alf2) rapidity [nb]$',ida+73 ,ymin,ymax,1d0)

      ymax = 0.006d0
      CALL plt1mc('CUMU', 'O(alf2) rapidity cumulative [nb]$',ida+73 ,ymin,ymax,1d0)

      CALL GLK_Operat(ids+76,'/',ids+72,ids+76, 1d0, 1d0)
      CALL GLK_Operat(ids+77,'/',ids+72,ids+77, 1d0, 1d0)
      CALL GLK_Operat(ids+78,'/',ids+72,ids+78, 1d0, 1d0)


      CALL GLK_Operat(ida+76,'/',ida+72,ida+76, 1d0, 1d0)
      CALL GLK_Operat(ida+77,'/',ida+72,ida+77, 1d0, 1d0)
      CALL GLK_Operat(ida+78,'/',ida+72,ida+78, 1d0, 1d0)

      ymin = -.50d0
      ymax =  .50d0
      CALL GLK_Operat(ids+2202,'/',ids+72,ids+2202, 1d0, 1d0)
      CALL plt1mc('    ','Int On/Off, KK2f MC, v-reconst. $',ids+2202,ymin,ymax,1d0)

      ymin = -0.05d0
      ymax =  0.05d0
      CALL plt1mc('    ','delMZ=+10MeV, KK2f MC, v-recon. O(alf1) NoInt$',ids +78,ymin,ymax,1d0)
      CALL plt1mc('    ','O(alf2-alf1)/O(alf1) v-reconst. $',ids+76,ymin,ymax,1d0)
      CALL plt1mc('    ','O(alf3-alf2)/O(alf1) v-reconst. $',ids+77,ymin,ymax,1d0)

      CALL plt1mc('    ','deltaMZ= +10MeV  MC  rapidity   $',ida+78,ymin,ymax,1d0)
      CALL plt1mc('    ','O(alf2-alf1)/O(alf1) rapidity   $',ida+76,ymin,ymax,1d0)
      CALL plt1mc('    ','O(alf3-alf2)/O(alf1) rapidity   $',ida+77,ymin,ymax,1d0)

      END
