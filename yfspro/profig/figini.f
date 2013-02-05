*///////////////////////////////////////////////////////////////////////
*//   Initial state emission
*//   Technical test of beta's MC versus SAnalytical
*//   MC results: yfspro/robol1
*//
*//   make figini-dvi
*//   make figini-ps
*//
*///////////////////////////////////////////////////////////////////////
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      SAVE
*---
      INTEGER imax
      PARAMETER(imax=10000)
      REAL*8  xpar(imax)
*---
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*---
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
      INTEGER   lint
      INTEGER   KeyZet,KeyISR
*----------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'figini'
      TeXfile   = 'figini.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

* Current run
c      Dname  = '../ini200/ini200.input' ! current
c      Hname  = '../ini200/pro.hst'      ! current

* Basic normal xenph=1.25 with fleps=1d-45, jlim2=4096
      Dname  = '../ini200/ini200.input.xenph=1.25__30M'  ! test
      Hname  = '../ini200/pro.hst.xenph=1.25__30M'       ! test

* Basic normal xenph=1.0, with fleps=1d-30, jlim2=1024 TO BE REPEATED!!!
c      Dname  = '../ini200/ini200.input.vmax=1_ewl=0_84M' ! jan 99
c      Hname  = '../ini200/pro.hst.vmax=1_ewl=0_84M'      ! jan 99
* Normal run
c       Dname  = '../ini200/ini200.input.vmax=1_ewl=0_78M' ! new jan 98
c       Hname  = '../ini200/pro.hst.vmax=1_ewl=0_78M'      ! new jan 98

* -------- Tests ----------------------------------------------------
c      Dname  = '../ini200/ini200.input.KeyZet-2__1M' !
c      Hname  = '../ini200/pro.hst.KeyZet-2__1M'      !

* Special runs with alternative ISR generator (no dilatation, flat xsection)
c      Dname  = '../ini140/ini140.input.keyisr=2.vmin-4.402M' ! vmin eff. seen
c      Hname  =      '../ini140/pro.hst.keyisr=2.vmin-4.402M' ! vmin eff. seen

c      Dname  = '../ini140/ini140.input.keyisr=2.vmin-6.738M'
c      Hname  =      '../ini140/pro.hst.keyisr=2.vmin-6.738M'

* Exercise with flat x-section
c      Dname  = '../ini140/ini140.input.keyzet-2.vmin-8.103M' ! flat xsec
c      Hname  =      '../ini140/pro.hst.keyzet-2.vmin-8.103M' ! flat xsec


*******************oldies******************
* Normal run
c      Dname  = '../ini140/ini140.data.vvmax=.9999.5M'  !
c      Hname  = '../ini140/pro.hst.vvmax=.9999.5M'      !


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
c[[[      CALL figtech3

      KeyISR = xpar(20)
      KeyZet = xpar(501)
* test of all beta's for flat x-section
      IF( KeyZet .EQ. -2  .OR. KeyISR .EQ. 2) THEN
         CALL figchi_bet
      ELSE
* normal run
*****    CALL figMasZ
         CALL figchi_lin
         CALL figchi_z
      ENDIF
* technical test
      IF( KeyISR .EQ. 1) THEN
         CALL figtech1
      ELSEIF( KeyISR .EQ. 2) THEN
         CALL figtech2
      ENDIF
*=================================
* end GLK_Plot, close LaTeX file
      CALL GLK_PlEnd
*===================================================================
* Write all histograms into dump file, for control/debug
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=================================
      CLOSE(nout)
      END

      SUBROUTINE figchi_bet
*-------------------------------------------------------------------
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_bet   =============='
      WRITE(nout,*) ' ====================================='

      WRITE(*,*) ' here we are !!!!!!'

      CALL Semalib_GetBorn(Born)

      idb  = 50000
*==============================================================
*=========================O(alf1)==============================
*==============================================================
*--- O(alf0) total
      KEYDIS=     300
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8300,' O(alf0) sig(vmax) $',IDB+1)
      CALL GLK_Print(8300)
*--------------------------------------------------------------
* GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS
      ymax =  0.013d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
c      CALL pldiff("CUMU","ISR GPS O(alf0) sig(vmax)/born ,YMAG=1$",
c     $     idb +201, 8300,   ymin,ymax, xmag,ymag, born)
*==============================================================
*=========================O(alf1)==============================
*==============================================================
*--- O(alf1) total
      KEYDIS=     301
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8301,' O(alf1) sig(ln(vmax)) $',IDB+1)
      CALL GLK_Print(8301)
*--- O(alf1) bt0xbt0 etc...
      KEYDIS=     310
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8310,' O(alf1) bt0xbt0  $',IDB+1)
      CALL GLK_Print(8310)
      KEYDIS=     311
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8311,' O(alf1) bt1xbt0  $',IDB+1)
      CALL GLK_Print(8311)
*==============================================================
* MC-SAN
      ymax =  0.013d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
*--------------------------------------------------------------
* GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS
c      CALL pldiff("CUMU","ISR GPS O(alf1) sig(vmax)/born ,YMAG=1$",
c     $     idb +202, 8301,   ymin,ymax, xmag,ymag, born)
*---------------------
      CALL pldiff("CUMU","ISR O(alf1) TOTAL, sig(vmax)/born ,YMAG=.01$",
     $     idb +2, 8301,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","init O(alf1), sig(vmax)/Born bt01, YMAG=.01$",
     $     idb+10, 8310,   ymin,ymax, xmag,ymag, born)

      ymag =  0.1
      CALL pldiff("CUMU","init O(alf1), sig(vmax)/Born bt11, YMAG=.1$",
     $     idb+11, 8311,   ymin,ymax, xmag,ymag, born)

*==============================================================
*=====================O(alf2)==================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8302,'O(alf2) sig(vmax) TOTAL$',IDB+1)
*--- O(alf2)  bt0xbt0 etc...
      KEYDIS=     320
      CALL Semalib_VVplot(KeyDis,"VCHI2",8320,' O(alf2) bt0xbt0  $',IDB+1)
      CALL GLK_Print(8320)
      KEYDIS=     321
      CALL Semalib_VVplot(KeyDis,"VCHI2",8321,' O(alf2) bt1xbt0  $',IDB+1)
      CALL GLK_Print(8321)
      KEYDIS=     322
      CALL Semalib_VVplot(KeyDis,"VCHI2",8322,' O(alf2) bt2xbt0  $',IDB+1)
      CALL GLK_Print(8322)
*==============================================================
* MC-SAN
      ymax =  0.015d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
      CALL pldiff("CUMU","init O(alf2) TOTAL, sig(vmax)/Born  ,YMAG=.01$",
     $     idb+3, 8302,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","init O(alf2), sig(vmax)/born  bt02, YMAG=.01$",
     $     idb+20, 8320,   ymin,ymax, xmag,ymag, born)

      ymag =  0.1
      CALL pldiff("CUMU","init O(alf2), sig(vmax)/born  bt1, YMAG=.1$",
     $     idb+21, 8321,   ymin,ymax, xmag,ymag, born)

      ymag =  1.
      CALL pldiff("CUMU","init O(alf2), sig(vmax)/born  bt2, YMAG=1$",
     $     idb+22, 8322,   ymin,ymax, xmag,ymag, born)

*==============================================================
*-----------------------O(alf2)-O(alf1) -----------------------
      ymax =  0.012d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
      CALL pld2mc("CUMU",'TOTAL O(alf2) vs O(alf1) ymag=0.01 $',
     $     idb+3,idb+2, ymin,ymax,xmag,ymag,born)
      ymag =  0.01
      CALL pld2mc("CUMU",'bt0 O(alf2) vs O(alf1) ymag=0.01 $',
     $     idb+20,idb+10, ymin,ymax,xmag,ymag,born)
      ymag =  0.1
      CALL pld2mc("CUMU",'bt1 O(alf2) vs O(alf1) ymag=0.1 $',
     $     idb+21,idb+11, ymin,ymax,xmag,ymag,born)
      CALL plt1mc("CUMU",  'bt2 O(alf2) only $',
     $     idb+22, ymin,ymax,born)
*-----------------------O(alf3)-O(alf2) -----------------------
cc      ymax =  0.0010d0
cc      ymin = -ymax
cc      xmag =  1
cc      ymag = .001
cc      CALL pld2mc("CUMU",'!!!TOTAL O(alf3) vs O(alf2) ymag=0.001 $',
cc     $     idb+4,idb+3, ymin,ymax,xmag,ymag,born)
cc      CALL pld2mc("CUMU",'!!!bt0 O(alf3) vs O(alf2) ymag=0.001 $',
cc     $     idb+30,idb+20, ymin,ymax,xmag,ymag,born)
cc      ymag = .01
cc      CALL pld2mc("CUMU",'!!!bt1 O(alf3) vs O(alf2) ymag=0.01 $',
cc     $     idb+31,idb+21, ymin,ymax,xmag,ymag,born)
cc      ymag = 0.1
cc      CALL pld2mc("CUMU",'!!!bt2 O(alf3) vs O(alf2) ymag=0.1 $',
cc     $     idb+32,idb+22, ymin,ymax,xmag,ymag,born)
cc      ymax =  0.00002d0
cc      ymin = -ymax
cc      CALL plt1mc("CUMU", '!!!bt3 O(alf3) only $',
cc     $     idb+33, ymin,ymax,born)
*-----------------------O(alf3)-O(alf2) -----------------------
* the case of histograms with differences
      ymax = 0.00020d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'TOTAL O(alf3)-O(alf2)!!!$',idb+ 7,ymin,ymax,born)
      CALL plt1mc("CUMU", 'bt0 O(alf3)-O(alf2)!!!  $',idb+30,ymin,ymax,born)
      ymax = 0.00020d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'bt1 O(alf3)-O(alf2)!!! $',idb+31, ymin,ymax,born)
      CALL plt1mc("CUMU", 'bt2 O(alf3)-O(alf2)!!! $',idb+32, ymin,ymax,born)
      ymax =  0.00002d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'bt3 O(alf3)!!!         $',idb+33, ymin,ymax,born)
*==============================================================
      CALL GLK_Delet(8301)
      CALL GLK_Delet(8310)
      CALL GLK_Delet(8311)
      CALL GLK_Delet(8302)
      CALL GLK_Delet(8320)
      CALL GLK_Delet(8321)
      CALL GLK_Delet(8322)
      END


      SUBROUTINE figMasZ
*//////////////////////////////////////////////////////////////////////////
*//                                                                      //
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      REAL*8  ymin,ymax,Born,xmag,ymag,Zmass
      INTEGER idm,ida
*-------------------------------------------------------------------
      ida =  30000

      CALL GLK_Operat(ida+6,'/',ida+2,ida+6, 1d0, 1d0)
      CALL GLK_Operat(ida+7,'/',ida+2,ida+7, 1d0, 1d0)
      CALL GLK_Operat(ida+8,'/',ida+2,ida+8, 1d0, 1d0)
      ymin = -0.05d0
      ymax =  0.05d0
      CALL plt1mc("    ",'O(alf2-alf1)/O(alf1)  $',ida+6 ,ymin,ymax,1d0)
      CALL plt1mc("    ",'O(alf3-alf2)/O(alf1)  $',ida+7 ,ymin,ymax,1d0)
      CALL plt1mc("    ",'deltaMZ= +10MeV !! MC $',ida+8 ,ymin,ymax,1d0)

      ymin = 0d0
      CALL GLK_Print(ida+3)
      ymax =0.006d0
      CALL plt1mc("CUMU", 'O(alf2)  [nb]$',ida+3 ,ymin,ymax,1d0)
      ymax =0.020d0
      CALL plt1mc("NB  ", 'O(alf2)  [nb]$',ida+3 ,ymin,ymax,1d0)

*/////////////////////////

      idm =  10000

      CALL Semalib_MZplot( 301,"VRHO2", 8301,'O(alf1) $',idm+2)
      CALL Semalib_MZplot( 302,"VRHO2", 8302,'O(alf2) $',idm+2)

* Difference O(alf2)-O(alf1) ANALITICALLY 
      CALL GLK_Operat(8302,'-',8301,5302, 1d0, 1d0)
      CALL GLK_Operat(5302,'/',8301,5302, 1d0, 1d0)

* calculate ANALITICALLY change due to MZ shift id=9302
      CALL Semalib_GetZmass(Zmass)
      WRITE(*,*) " original Zmass = ",Zmass
      Zmass=Zmass -0.005d0
      WRITE(*,*) " changed  Zmass = ",Zmass
      CALL Semalib_SetZmass(Zmass)
      CALL Semalib_MZplot( 302,"VRHO2", 9302,'O(alf2) $',idm+2)
      CALL GLK_Print(9302)
      CALL GLK_Operat(9302,'-',8302,9302, 1d0, 1d0)
      CALL GLK_Operat(9302,'/',8302,9302, 1d0, 1d0)

      ymin = 0d0
      ymax =0.006d0
      CALL plt1mc("CUMU", 'O(alf1)  [nb]$',idm+2 ,ymin,ymax,1d0)
      ymax =0.001d0
      CALL GLK_Print(idm+1)
      CALL plt1mc("NB  ", 'O(alf1)  [nb]$',idm+2 ,ymin,ymax,1d0)
c[[      CALL plt1mc("NB  ", 'O(alf2)  [nb]$',idm+3 ,ymin,ymax,1d0)
c[[      CALL plt1mc("NB  ", 'O(alf3)  [nb]$',idm+4 ,ymin,ymax,1d0)

      CALL GLK_Operat(idm+6,'/',idm+2,idm+6, 1d0, 1d0)
      CALL GLK_Operat(idm+7,'/',idm+2,idm+7, 1d0, 1d0)
      CALL GLK_Operat(idm+8,'/',idm+2,idm+8, 1d0, 1d0)

      ymin = -0.05d0
      ymax =  0.05d0
      CALL plt1an("O(alf2) relat. difference due to MZ shift of -5MeV$",
     $                     9302, 9802, ymin, ymax, 1d0)
      CALL plt1an("O(alf2-alf1)/O(alf1) SAN$",
     $                     5302, 5302, ymin, ymax, 1d0)
      CALL plt1mc("    ", 'O(alf2-alf1)/O(alf1) $',idm+6 ,ymin,ymax,1d0)
      ymin = -0.0005d0
      ymax =  0.0005d0
      CALL plt1mc("    ", 'O(alf3-alf2)/O(alf1) $',idm+7 ,ymin,ymax,1d0)
      ymin = -0.05d0
      ymax =  0.05d0
      CALL plt1mc("    ", 'deltaMZ= +10MeV !! MC$',idm+8 ,ymin,ymax,1d0)

      ymax =  0.0007d0
      ymin = -0.0001d0
      xmag =  1
      ymag =  1

      CALL pldiff("NB  ","[301] O(alf1) MC.vs.SAN, YMAG=.01$",
     $     idm+2, 8301,   ymin,ymax, xmag,ymag, 1d0)
      CALL pldiff("NB  ","[302] O(alf2) MC.vs.SAN, YMAG=.01$",
     $     idm+3, 8302,   ymin,ymax, xmag,ymag, 1d0)


      CALL GLK_Delet(8301)
      CALL GLK_Delet(8302)
      CALL GLK_Delet(9302)
      CALL GLK_Delet(9802)
      CALL GLK_Delet(5302)

      END

      SUBROUTINE figchi_lin
*//////////////////////////////////////////////////////////////////////////
*//                                                                      //
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_lin   =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      idb  = 50000
      idz  = 70000
*==============================================================
*     Basics
*==============================================================
****  CALL Semalib_VVplot(  51,"VCHI2", 8051,' [51]  WTCRUDE           $',idb+1)
      CALL Semalib_VVplot( 300,"VCHI2", 8300,' [300] bt00(vv) Standard $',idb+1)
      CALL Semalib_VVplot( 400,"VCHI2", 8400,' [400] bt00(vv) NLLphsp  $',idb+1)
* !!!!!!! debug !!!!!!!
      CALL Semalib_VVplot( 400,"ZCHI1", 9400,' [400] bt00(zz) NLLphsp  $',idz+1)
*==============================================================
* MC-SAN
      ymax =  0.075d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
cc      CALL pldiff("CUMU","[51] WTCRUDE, sig(vmax)/born ,YMAG=.01$",
cc     $     idb+51, 8051,   ymin,ymax, xmag,ymag, born)

cc      CALL pld2an('[400] test of gauss integration, YMAG=0.01$',
cc     $     8400,9400,ymin,ymax,xmag,ymag,born)

cc      CALL pld2mc("CUMU",'[400] Test of MC normalization, YMAG=0.01$',
cc     $     idb+1,idz+1, ymin,ymax,xmag,ymag,born)

      CALL pldiff("CUMU","[400] O(alf0) beta00, NLLphsp!!! sig(vmax)/born ,YMAG=.01$",
     $     idb+ 1, 8400,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[300] O(alf0) beta00 Standard, sig(vmax)/born ,YMAG=.01$",
     $     idb+ 1, 8300,   ymin,ymax, xmag,ymag, born)
*--------------------------------------------------------------
* GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS
c      CALL pldiff("CUMU","GPS  O(alf0) vs [300], sig(vmax)/born ,YMAG=.01$",
c     $     idb+ 201, 8300,   ymin,ymax, xmag,ymag, born)
*==============================================================
*--- O(alf1) total
      KEYDIS=     301
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8301,' O(alf1) sig(ln(vmax)) $',idb+1)
*--- O(alf1) bt0
      KEYDIS=     310
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8310,' O(alf1) bt0  $',idb+1)
*--- O(alf1) bt1
      KEYDIS=     311
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8311,' O(alf1) bt1xbt0  $',IDB+1)
*==============================================================
* MC-SAN
      ymax =  0.075d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01

      CALL pldiff("CUMU","[301] O(alf1) TOTAL, sig(vmax)/born ,YMAG=.01$",
     $     idb+ 2, 8301,   ymin,ymax, xmag,ymag, born)
*--------------------------------------------------------------
* GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS
c      CALL pldiff("CUMU","GPS O(alf1) vs. [301], sig(vmax)/born ,YMAG=.01$",
c     $     idb+ 202, 8301,   ymin,ymax, xmag,ymag, born)

c      CALL pldiff("CUMU","[310] O(alf1), sig(vmax)/Born bt0, YMAG=.01$",
c     $     idb+10, 8310,   ymin,ymax, xmag,ymag, born)

c      CALL pldiff("CUMU","[311] O(alf1), sig(vmax)/Born bt1, YMAG=.01$",
c     $     idb+11, 8311,   ymin,ymax, xmag,ymag, born)

*==============================================================
*=====================O(alf2)==================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8302,'O(alf2) sig(vmax) TOTAL$',IDB+1)
*--- O(alf2)
      KEYDIS=     320
      CALL Semalib_VVplot(KeyDis,"VCHI2",8320,' O(alf2) bt0  $',IDB+1)
      KEYDIS=     321
      CALL Semalib_VVplot(KeyDis,"VCHI2",8321,' O(alf2) bt1  $',IDB+1)
      KEYDIS=     322
      CALL Semalib_VVplot(KeyDis,"VCHI2",8322,' O(alf2) bt2  $',IDB+1)
*--- O(alf3) BEST!!!
      KEYDIS=     305
      CALL Semalib_VVplot(KeyDis,"VCHI2", 8305,'O(alf3) sig(vmax) BEST $',IDB+1)
*==============================================================
* MC-SAN
      ymax =  0.075d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
      CALL pldiff("CUMU","[302] O(alf2) TOTAL, sig(vmax)/Born  ,YMAG=.01$",
     $     idb +3, 8302,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[305]Best O(alf3) TOTAL, sig(vmax)/Born  ,YMAG=.01$",
     $     idb +4, 8305,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[320] O(alf2), sig(vmax)/born  bt0, YMAG=.01$",
     $     idb+20, 8320,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[321] O(alf2), sig(vmax)/born  bt1, YMAG=.01$",
     $     idb+21, 8321,   YMIN,YMAX, XMAG,YMAG, BORN)
      ymag =  1
      CALL pldiff("CUMU","[320] O(alf2), sig(vmax)/born  bt2, YMAG=1$",
     $     idb+22, 8322,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[322] O(alf3), sig(vmax)/born  bt3 !!!!, YMAG=1$",
     $     idb+33, 8322,   ymin,ymax, xmag,ymag, born)

*==============================================================
      CALL GLK_Delet(8300)
      CALL GLK_Delet(8301)
      CALL GLK_Delet(8310)
      CALL GLK_Delet(8311)
      CALL GLK_Delet(8302)
      CALL GLK_Delet(8305)
      CALL GLK_Delet(8320)
      CALL GLK_Delet(8321)
      CALL GLK_Delet(8322)
      CALL GLK_Delet(8400)
      CALL GLK_Delet(9400)

      END  !figchi_lin 


      SUBROUTINE figchi_z
*-------------------------------------------------------------------
* sig(ln(vmax))
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_z     =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      idz  = 70000
*==============================================================
*--- O(alf1) total
      KEYDIS=     301
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7301,'ini O(alf1) sig(z) $',idz+1)
      CALL GLK_Print(7301)
*--- O(alf1) bt0xbt0 etc...
      KEYDIS=     310
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7310,' O(alf1) bt0xbt0  $',idz+1)
      CALL GLK_Print(7310)
      KEYDIS=     311
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7311,' O(alf1) bt1xbt0  $',idz+1)
      CALL GLK_Print(7311)
*==============================================================
* MC-SAN
      ymax =  .075
      ymin = -ymax
      xmag =  1
      ymag =  .01

      CALL pldiff("CUMU",
     $     "init O(alf1) TOTAL, sig(zmax)/born ,YMAG=.01$",
     $     idz+2, 7301,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU",
     $     "init O(alf1), sig(zmax)/Born bt0, YMAG=.01 [310]$",
     $     idz+10, 7310,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU",
     $     "init O(alf1), sig(zmax)/Born bt1, YMAG=.01$",
     $     idz+11, 7311,   ymin,ymax, xmag,ymag, born)
*==============================================================
*--- O(alf2) total
      KEYDIS=     302
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7302,'O(alf2) sig(zmax) TOTAL$',idz+1)
      CALL GLK_Print(7302)
*--- O(alf2)  bt0xbt0 etc...
      KEYDIS=     320
      CALL Semalib_VVplot(KeyDis,"ZCHI1",7320,' O(alf2) bt0  $',idz+1)
      CALL GLK_Print(7320)
      KEYDIS=     321
      CALL Semalib_VVplot(KeyDis,"ZCHI1",7321,' O(alf2) bt1  $',idz+1)
      CALL GLK_Print(7321)
      KEYDIS=     322
      CALL Semalib_VVplot(KeyDis,"ZCHI1",7322,' O(alf2) bt2  $',idz+1)
      CALL GLK_Print(7322)
*==============================================================
* MC-SAN
      ymax =  .075
      ymin = -ymax
      xmag =  1
      ymag =  .01
      CALL pldiff("CUMU",
     $     "init O(alf2) TOTAL, sig(zmax)/born ,YMAG=.01$",
     $     idz+3, 7302,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU",
     $     "init O(alf2), sig(zmax)/born  bt0, YMAG=.01$",
     $     idz+20, 7320,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU",
     $     "init O(alf2), sig(zmax)/born  bt1, YMAG=.01$",
     $     idz+21, 7321,   ymin,ymax, xmag,ymag, born)

      ymag =  1
      CALL pldiff("CUMU",
     $     "init O(alf2), sig(zmax)/born  bt2, YMAG=1$",
     $     idz+22, 7322,   ymin,ymax, xmag,ymag, born)

*-----------------------O(alf3)-O(alf2) -----------------------
      ymax =  0.020d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'TOTAL O(alf3)-O(alf2)!!! $',idz+7 ,ymin,ymax,born)
      CALL plt1mc("CUMU", 'bt0   O(alf3)-O(alf2)!!! $',idz+30,ymin,ymax,born)
      ymax =  0.020d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'bt1 O(alf3)-O(alf2)!!!   $',idz+31,ymin,ymax,born)
      CALL plt1mc("CUMU", 'bt2 O(alf3)-O(alf2)!!!   $',idz+32,ymin,ymax,born)
      ymax =  0.002d0
      ymin = -ymax
      CALL plt1mc("CUMU", 'bt3 O(alf3)!!!           $',idz+33,ymin,ymax,born)
*==============================================================
      CALL GLK_Delet(7301)
      CALL GLK_Delet(7302)

      END  !figchi_z


      SUBROUTINE figtech1
*-------------------------------------------------------------------
* technical test
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figtech      =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      idz  = 70000
      idb  = 50000
*==============================================================
*--- O(alf0) bt0xbt0 etc...
      KEYDIS=     400
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7400,' O(alf0) bt0 NLLphsp  $',idz+1)
      CALL GLK_Print(7400)
*--- wtcrude
      KEYDIS=      51
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7051,' [51]  $',idz+1)
      CALL GLK_Print(7051)
*----
      KEYDIS=      52
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7052,' [52]  $',idz+1)
      CALL GLK_Print(7052)
*----
      KEYDIS=      53
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7053,' [53]  $',idz+1)
      CALL GLK_Print(7053)
*----
      KEYDIS=      54
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7054,' [54]  $',idz+1)
      CALL GLK_Print(7054)
*----
      KEYDIS=      55
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7055,' [55]  $',idz+1)
      CALL GLK_Print(7055)
*==============================================================
* MC-SAN
cc      ymax =  .250   ! for keyzet = 1 !!!! alfinv =65
cc      ymax =  .0011  ! for keyzet = -2
cc      ymag =  .001   ! for keyzet = -2

      ymax =  .050 ! for keyzet = 1
      ymin = -ymax
      ymag =  .0002  ! for keyzet = 1
      xmag =  1
ccc      CALL GLK_idopt(idz+52,   'ERRO')
ccc      CALL GLK_Print(idz+52)
      CALL pldiff("CUMU","[52] pure vesko, YMAG=.002 $",
     $     idz+52, 7052,   ymin,ymax, xmag,ymag, born)

      ymax =  .150 ! for keyzet = 1
      ymax =  .050 ! for keyzet = 1
      ymin = -ymax
      ymag =  .01  ! for keyzet = 1
      xmag =  1

      CALL pldiff("CUMU","[400] beta0 NLL phsp !!! , YMAG=.01$",
     $     idz+ 1, 7400,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[51] WTcrude, YMAG=.01 $",
     $     idz+51, 7051,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[53] WTcrude primitive, YMAG=.01 $",
     $     idz+53, 7053,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[54] WTcrude primitive, YMAG=.01 $",
     $     idz+54, 7054,   ymin,ymax, xmag,ymag, born)

***      CALL plt1mc("CUMU", ' [55] wt1-wt2 $',
***     $     idz+55, ymin,ymax,born)

      ymax =  .05               ! for keyzet = 1
      ymag =  .2                ! for keyzet = 1
cc      ymax =  .0011 ! for keyzet = -2
cc      ymag =   .10  ! for keyzet = -2
      ymin = -ymax
      xmag=1
      CALL pldiff("CUMU","[55] , wt1-wt2, YMAG=0.2 $",
     $     idz+55, 7055,   ymin,ymax, xmag,ymag, born)

*--------------------------------------------------*
*                      UNEXP                       *
*--------------------------------------------------*
*---
      KEYDIS=      660
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7160,' [660]  $',idz+1)
      CALL GLK_Print(7160)
*---
      KEYDIS=      661
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7161,' [661]  $',idz+1)
      CALL GLK_Print(7161)
*---
      KEYDIS=      662
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7162,' [662]  $',idz+1)
      CALL GLK_Print(7162)
*---
      KEYDIS=      663
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7163,' [663]  $',idz+1)
      CALL GLK_Print(7163)

      ymax =  .070 ! for keyzet = 1
      ymag =  .01  ! for keyzet = 1
cc      ymax =  .011 ! for keyzet = -2
cc      ymag =  .010 ! for keyzet = -2
ccccc      ymax =  .250 ! for keyzet = 1 !!!! alfinv =65
      ymin = -ymax
      xmag =  1
ccc      CALL plt1mc("CUMU", ' [161] $',$     idz+161, ymin,ymax,born)
*
      CALL pldiff("CUMU","[160] , YMAG=1 $",
     $     idz+160, 7160,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[161] , YMAG=1 $",
     $     idz+161, 7161,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[162] , YMAG=1 $",
     $     idz+162, 7162,   ymin,ymax, xmag,ymag, born)

      ymax =  .20   ! for keyzet = 1
      ymag =    1   ! for keyzet = 1
cc      ymax =  .0011  ! for keyzet = -2
cc      ymag =   .10  ! for keyzet = -2
cccccc      ymax =  1.0 ! for keyzet = 1 !!!!!!! alfinv =65
      ymin = -ymax
      xmag = 1
      CALL pldiff("CUMU","[163] , YMAG=1 dilatation effect $",
     $     idz+163, 7163,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[164] , YMAG=1 approx. NLL dilatation $",
     $     idz+164, 7163,   ymin,ymax, xmag,ymag, born)

      END  !figtech1



      SUBROUTINE figtech2
*-------------------------------------------------------------------
* technical test
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figtech      =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      idz  = 70000
*==============================================================
*--- O(alf0) bt0
      KEYDIS=     300
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7300,' O(alf0) bt0xbt0  $',idz+1)
      CALL GLK_Print(7300)
*--- wtcrude
      KEYDIS=      51
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7051,' [51]  $',idz+1)
      CALL GLK_Print(7051)
*==============================================================
* MC-SAN
      ymax =  .0011  ! for keyzet = -2
      ymag =  .001  ! for keyzet = -2
      ymin = -ymax
      xmag =  1

      CALL pldiff("CUMU","[300] beta0, YMAG=.01$",
     $     idz+ 1, 7300,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[51] WTcrude, YMAG=.01 $",
     $     idz+51, 7051,   ymin,ymax, xmag,ymag, born)

*--------------------------------------------------*
*                      UNEXP                       *
*--------------------------------------------------*
*---
      KEYDIS=      660
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7160,' [660]  $',idz+1)
      CALL GLK_Print(7160)
*---
      KEYDIS=      661
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7161,' [661]  $',idz+1)
      CALL GLK_Print(7161)
*---
      KEYDIS=      662
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7162,' [662]  $',idz+1)
      CALL GLK_Print(7162)
*---
      KEYDIS=      663
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7163,' [663]  $',idz+1)
      CALL GLK_Print(7163)

      ymax =  .0011  ! for keyzet = -2
      ymag =  .001  ! for keyzet = -2
      ymin = -ymax
      xmag =  1

      CALL pldiff("CUMU","[160] , YMAG=1 $",
     $     idz+160, 7160,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[161] , YMAG=1 $",
     $     idz+161, 7161,   ymin,ymax, xmag,ymag, born)
      CALL pldiff("CUMU","[162] , YMAG=1 $",
     $     idz+162, 7162,   ymin,ymax, xmag,ymag, born)

      CALL pldiff("CUMU","[164] , YMAG=1 vvnll !!!$",
     $     idz+164, 7162,   ymin,ymax, xmag,ymag, born)

***********************************************************
* very special tests for forward-forward, forward-backward
* photon emission
      ymax =  .0011  ! for keyzet = -2
      ymax =  .0022  !!!!!!!!
      ymag =   .10  ! for keyzet = -2
      ymin = -ymax
      xmag =  1

*      CALL pldiff("CUMU","[165] , YMAG=0.1 z10 versus zref $",
*     $     idz+165, 7163,   ymin,ymax, xmag,ymag, born)

*      CALL plt1mc("CUMU", ' [166] z10 versus znll $',
*     $     idz+166, ymin,ymax,born)

      END  !figtech2

      SUBROUTINE figtech3
*-------------------------------------------------------------------
* technical test
*-------------------------------------------------------------------
      IMPLICIT NONE
      COMMON / inout  / ninp,nout
      INTEGER           ninp,nout

      INTEGER           idz,keydis,idb
      DOUBLE PRECISION  ymin,ymax,xmag,ymag,Born
*--------------------------------------------------
      CALL Semalib_GetBorn(Born)
      idz  = 70000
      idb  = 50000
*----
      KEYDIS=      52
      CALL Semalib_VVplot(KeyDis,"ZCHI1", 7052,' [52]  $',idz+1)
      CALL GLK_Print(7052)

      ymax =  .100 ! for keyzet = 1
      ymin = -ymax
      ymag =  .0002  ! for keyzet = 1
      xmag =  1

      write(nout,*) '[[[[[[[[[[[[[[[[['
      CALL GLK_Print(idz+52)
      CALL pldiff("CUMU","[52] pure vesko, YMAG=.002 $",
     $     idz+52, 7052,   ymin,ymax, xmag,ymag, born)
      write(nout,*) ']]]]]]]]]]]]]]]]]'

      END  !figtech2
