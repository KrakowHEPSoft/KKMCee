*////////////////////////////////////////////////////////////////////
*//   Initial + Final state plots
*//   Technical test of beta's MC versus SAnalytical
*//   make figbet-dvi
*//   make figbet-ps
*//   make figbet-pubs <==== THIS IS MAIN RESULT
*////////////////////////////////////////////////////////////////////
      PROGRAM MAIN
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON / inout  / ninp,nout
      SAVE
*---
      INTEGER   imax
      PARAMETER(imax=10000)
      REAL*8  xpar(imax)
*
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
*-------------------------------------------------------------
*
      ninp=  5
      nout= 16
      Tesnam    = 'figbet'
      TeXfile   = 'figbet.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

* Exercise with flat x-section
ccc      Dname  = '../mix200/mix200.input' ! actual
ccc      Hname  = '../mix200/pro.hst'      ! actual

      Dname  = '../mix200/mix200_flat.input'     ! KeyZet=-2 ****
      Hname  = '../mix200/pro.hst.flat.16M'      ! KeyZet=-2 ****

*********************oldies******************************
cc      Dname  = '../mix2000/mix2000.input' ! KeyZet=-2
cc      Hname  = '../mix2000/pro.hst'       ! KeyZet=-2
* Exercise with flat x-section
*      Dname  = '../mixflat/5M.mix40.data'   ! KeyZet=-2 new
*      Hname  = '../mixflat/5M.pro40.hst'    ! KeyZet=-2 new
* Exercise with flat x-section
c      Dname  = '../mixflat/5M.mix200.data'  ! KeyZet=-2 new
c      Hname  = '../mixflat/5M.pro200.hst'   ! KeyZet=-2 new
*

*=====================================
* Read data, the same as in MC run
      CALL ReaDat(Dname,imax,xpar)
      CALL Semalib_Initialize(xpar)
* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
*=====================================
* Initialize GLK_Plot
      Lint=0                    ! Lint=0 several pages
      CALL GLK_PlInitialize(Lint,TeXfile)
*=====================================================
cccc      CALL figrho_log
cccc      CALL figchi_log
cccc      CALL figchi_lin  !  <==== THIS IS MAIN RESULT
cccc      CALL third_ord
*=========================================================================
* end GLK_Plot, close LaTeX file
      CALL GLK_PlEnd
*=========================================================================
      CALL fig_pubs  !  <==== THIS IS for publication
*=========================================================================
* Write all histograms into dump file, for control/debug
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=====================================================
      CLOSE(nout)
      END


      SUBROUTINE figrho_log
*-------------------------------------------------------------------
*  rho(ln(v)) 
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figrho_log   =============='
      WRITE(nout,*) ' ====================================='


      CALL Semalib_GetBorn(Born)
      ida = 30000
      idd = 35000
      idf = 40000
*==============================================================
*=======================O(alf1)================================
*==============================================================
* Analytical curves  **** dsigma/dx ****
*--- O(alf1) total
      KEYDIS=     301301
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,'O(alf1) rho(ln(v))$',ida+71)
*--- O(alf1) bt0xbt0
      KEYDIS=     310310
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf1) bt0xbt0   $',ida+71)
*--- O(alf1) bt0xbt1
      KEYDIS=     310311
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf1) bt0xbt1   $',ida+71)
*--- O(alf1) bt1xbt0
      KEYDIS=     311310
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf1) bt1xbt0   $',ida+71)
*==============================================================
      ymax =  0.012d0
      ymin = -0.012d0
      xmag =  1
      ymag =  0.05

      CALL PLDIFF("NB10","O(alf1) TOTAL,  rho(ln(v)) , YMAG=.05$",
     $     IDA+72,  301301,   YMIN,YMAX, XMAG,YMAG, BORN)

      CALL PLDIFF("NB10","O(alf1) ,rho(ln(x)) bt0xbt0, YMAG=.05$",
     $     IDA+80, 310310,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("NB10","O(alf1),rho(ln(x)) bt1*bt0, YMAG=.1$",
     $     IDA+82, 311310,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("NB10","O(alf1), rho(ln(x)) bt0*bt1, YMAG=.1$",
     $     IDA+83, 310311,   YMIN,YMAX, XMAG,YMAG, BORN)

*==============================================================
*====================== O(alf2)================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302302
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,'O(alf2) rho(ln(v))$',ida+71)
*--- O(alf2) bt0xbt0 etc...
      KEYDIS=     320320
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt0xbt0   $',ida+71)
      KEYDIS=     321320
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt1xbt0   $',ida+71)
      KEYDIS=     320321
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt0xbt1   $',ida+71)
      KEYDIS=     322320
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt2xbt0   $',ida+71)
      KEYDIS=     320322
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt0xbt2   $',ida+71)
* finally O(alf2) beta1*beta1
      KEYDIS=     321321
      CALL Semalib_VVplot(KeyDis,"XRHO2",-KEYDIS,' O(alf2) bt1xbt1   $',ida+71)
*==========================================================
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.05
      CALL PLDIFF("NB10","O(alf2) TOTAL, rho(ln(v)), YMAG=.05$",
     $     IDA+73,  302302,   YMIN,YMAX, XMAG,YMAG, BORN)

      CALL PLDIFF("NB10","O(alf2), rho(ln(x))  bt0xbt0, YMAG=.05$",
     $     IDA+90, 320320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("NB10","O(alf2), rho(ln(x))  bt1*bt0, YMAG=.1$",
     $     IDA+93, 321320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("NB10","O(alf2), rho(ln(x))  bt0*bt1, YMAG=.1$",
     $     IDA+94, 320321,   YMIN,YMAX, XMAG,YMAG, born)

      YMAG =  1.
      CALL PLDIFF("NB10","O(alf2), rho(ln(x))  bt2*bt0, YMAG=1$",
     $     IDA+95, 322320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  1
      CALL PLDIFF("NB10","O(alf2), rho(ln(x))  bt0*bt2, YMAG=1$",
     $     IDA+97, 320322,   YMIN,YMAX, XMAG,YMAG, born)

      YMAG =  1
      CALL PLDIFF("NB10","O(alf2) rho(ln(x)),  bt1*bt1, YMAG=1$",
     $     IDA+96, 321321,   YMIN,YMAX, XMAG,YMAG, born)
*==============================================================
*==============================================================
*==============================================================
      END ! figrho_log


      SUBROUTINE figchi_log
*-------------------------------------------------------------------
* sig(ln(vmax))
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ==========    figchi_log   =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      IDA = 30000
      IDD = 35000
      IDF = 40000
*==============================================================
*=========================O(alf1)==============================
*==============================================================
*--- O(alf1) total
      KEYDIS=     301301
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7301301,
     $     ' O(alf1) sig(ln(vmax)) $',ida+71)
*--- O(alf1)  bt0xbt0 etc...
      KEYDIS=     310310
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7310310,' O(alf1) bt0xbt0  $',ida+71)
      KEYDIS=     310311
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7310311,' O(alf1) bt0xbt1  $',ida+71)
      KEYDIS=     311310
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7311310,' O(alf1) bt1xbt0  $',ida+71)
*==============================================================
* MC-SAN
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01

      CALL PLDIFF("CUMU","O(alf1) TOTAL, sig(ln(vmax)) ,YMAG=.01$",
     $     IDA+72, 7301301,   YMIN,YMAX, XMAG,YMAG, BORN)

      CALL PLDIFF("CUMU","O(alf1), sig(ln(xmax)) bt0xbt0, YMAG=.01$",
     $     IDA+80, 7310310,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf1), sig(ln(xmax)) bt1*bt0, YMAG=.1$",
     $     IDA+82, 7311310,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf1), sig(ln(xmax)) bt0*bt1, YMAG=.1$",
     $     IDA+83, 7310311,   YMIN,YMAX, XMAG,YMAG, BORN)

*==============================================================
*=====================O(alf2)==================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302302
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7302302,'O(alf2) sig(ln(vmax)) TOTAL$',ida+71)
*---  O(alf2) bt0xbt0 etc...
      KEYDIS=     320320
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7320320,' O(alf2) bt0xbt0  $',ida+71)
      KEYDIS=     321320
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7321320,' O(alf2) bt1xbt0  $',ida+71)
      KEYDIS=     320321
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7320321,' O(alf2) bt0xbt1  $',ida+71)
      KEYDIS=     322320
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7322320,' O(alf2) bt2xbt0  $',ida+71)
      KEYDIS=     320322
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7320322,' O(alf2) bt0xbt2  $',ida+71)
* finally beta1*beta1
      KEYDIS=     321321
      CALL Semalib_VVplot(KeyDis,"XCHI2",-7321321,' O(alf2) bt1xbt1  $',ida+71)
********* sigma(xmax) ******************
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01
      CALL PLDIFF("CUMU","O(alf2) TOTAL, sig(ln(vmax))  ,YMAG=.01$",
     $     IDA+73, 7302302,   YMIN,YMAX, XMAG,YMAG, BORN)

      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt0xbt0, YMAG=.01$",
     $     IDA+90, 7320320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt1*bt0, YMAG=.1$",
     $     IDA+93, 7321320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt0*bt1, YMAG=.1$",
     $     IDA+94, 7320321,   YMIN,YMAX, XMAG,YMAG, born)

      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt2*bt0, YMAG=1$",
     $     IDA+95, 7322320,   YMIN,YMAX, XMAG,YMAG, BORN)

      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt0*bt2, YMAG=1$",
     $     IDA+97, 7320322,   YMIN,YMAX, XMAG,YMAG, born)

      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(ln(x))/born  bt1*bt1, YMAG=1$",
     $     IDA+96, 7321321,   YMIN,YMAX, XMAG,YMAG, born)

*==============================================================
*==============================================================
*==============================================================
      END ! figchi_lin



      SUBROUTINE figchi_lin
*-------------------------------------------------------------------
* sig(vmax)
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_lin   =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      IDA = 30000
      IDD = 35000
      IDF = 40000
      idb  = 50000
*==============================================================
*=========================O(alf1)==============================
*==============================================================
*--- O(alf1) total
      KEYDIS=     301301
      CALL Semalib_VVplot(KeyDis,"XCHI2", 8301301,' O(alf1) sig(ln(vmax)) $',IDB+72)
*--- O(alf1) bt0xbt0 etc...
      KEYDIS=     310310
      CALL Semalib_VVplot(KeyDis,"XCHI2", 8310310,' O(alf1) bt0xbt0  $',IDB+72)
      KEYDIS=     310311
      CALL Semalib_VVplot(KeyDis,"XCHI2", 8310311,' O(alf1) bt0xbt1  $',IDB+72)
      KEYDIS=     311310
      CALL Semalib_VVplot(KeyDis,"XCHI2", 8311310,' O(alf1) bt1xbt0  $',IDB+72)
*==============================================================
* MC-SAN
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01
      CALL PLDIFF("CUMU","O(alf1) TOTAL, sig(vmax)/born  ,YMAG=.01$",
     $     IDB+72, 8301301,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL PLDIFF("CUMU","O(alf1), sig(xmax)/Born bt0xbt0,YMAG=.01$",
     $     IDB+80, 8310310,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf1), sig(xmax)/Born bt1*bt0,YMAG=.1$",
     $     IDB+82, 8311310,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf1), sig(xmax)/Born bt0*bt1,YMAG=.1$",
     $     IDB+83, 8310311,   YMIN,YMAX, XMAG,YMAG, BORN)
*==============================================================
*=====================O(alf2)==================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302302
      CALL Semalib_VVplot(KeyDis,"XCHI2",8302302,' O(alf2) sig(vmax) TOTAL$',IDB+72) !
*--- O(alf2)  bt0xbt0 etc...
      KEYDIS=     320320
      CALL Semalib_VVplot(KeyDis,"XCHI2",8320320,' O(alf2) bt0xbt0  $',IDB+72) !
      KEYDIS=     321320
      CALL Semalib_VVplot(KeyDis,"XCHI2",8321320,' O(alf2) bt1xbt0  $',IDB+72) !
      KEYDIS=     320321
      CALL Semalib_VVplot(KeyDis,"XCHI2",8320321,' O(alf2) bt0xbt1  $',IDB+72) !
      KEYDIS=     322320
      CALL Semalib_VVplot(KeyDis,"XCHI2",8322320,' O(alf2) bt2xbt0  $',IDB+72) !
      KEYDIS=     320322
      CALL Semalib_VVplot(KeyDis,"XCHI2",8320322,' O(alf2) bt0xbt2  $',IDB+72) !
* finally beta1*beta1
      KEYDIS=     321321
      CALL Semalib_VVplot(KeyDis,"XCHI2",8321321,' O(alf2) bt1xbt1  $',IDB+72) !
********* sigma(xmax) ******************
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01
      CALL PLDIFF("CUMU","O(alf2) TOTAL, sig(xmax)/Born   , YMAG=.01$",
     $     IDB+73, 8302302,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0xbt0, YMAG=.01$",
     $     IDB+90, 8320320,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt1*bt0, YMAG=.1$",
     $     IDB+93, 8321320,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0*bt1, YMAG=.1$",
     $     IDB+94, 8320321,   YMIN,YMAX, XMAG,YMAG, born)
      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt2*bt0, YMAG=1$",
     $     IDB+95, 8322320,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0*bt2, YMAG=1$",
     $     IDB+97, 8320322,   YMIN,YMAX, XMAG,YMAG, born)
      YMAG =  1.
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt1*bt1, YMAG=1$",
     $     IDB+96, 8321321,   YMIN,YMAX, XMAG,YMAG, born)
*==============================================================
      END  !figchi_lin 

      SUBROUTINE third_ord
*-------------------------------------------------------------------
* third_ord: sig(vmax)
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout

      CALL Semalib_GetBorn(Born)
      idb  = 50000
*==============================================================
*-----------------------O(alf2)-O(alf1) -----------------------
      ymax =  0.012d0
      ymin = -ymax
      xmag =  1
      ymag =  0.01
      CALL pld2mc("CUMU",'TOTAL O(alf2) vs O(alf1) ymag=0.01 $',
     $     idb+73,idb+72, ymin,ymax,xmag,ymag,born)
      ymag =  0.01
      CALL pld2mc("CUMU",'bt0 O(alf2) vs O(alf1) ymag=0.01 $',
     $     idb+90,idb+80, ymin,ymax,xmag,ymag,born)
      ymag =  0.1
      CALL pld2mc("CUMU",'bt1 O(alf2) vs O(alf1) ymag=0.1 $',
     $     idb+91,idb+81, ymin,ymax,xmag,ymag,born)
      CALL plt1mc("CUMU",'bt2 O(alf2) only $',
     $     idb+92, ymin,ymax,born)
*==============================================================
*-----------------------O(alf3)-O(alf2) -----------------------
      ymax =  0.0012d0
      ymin = -ymax
      xmag =  1
      ymag = .001
      CALL pld2mc("CUMU",'!!!TOTAL O(alf3) vs O(alf2) ymag=0.001 $',
     $     idb+74,idb+73, ymin,ymax,xmag,ymag,born)
      CALL plt1mc("CUMU", 'TOTAL O(alf3)-O(alf2)!!!$',
     $     idb+77, ymin,ymax,born)
      CALL pld2mc("CUMU",'!!!bt0 O(alf3) vs O(alf2) ymag=0.001 $',
     $     idb+100,idb+90, ymin,ymax,xmag,ymag,born)
      ymag = .01
      CALL pld2mc("CUMU",'!!!bt1 O(alf3) vs O(alf2) ymag=0.01 $',
     $     idb+101,idb+91, ymin,ymax,xmag,ymag,born)
      ymag = .2
      CALL pld2mc("CUMU",'!!!beta21-beta20, O(alf3) vs O(alf2) ymag=0.2 $',
     $     idb+102,idb+92, ymin,ymax,xmag,ymag,born)
      CALL pld2mc("CUMU",'!!!xbet21-xbet20, O(alf3) vs O(alf2) ymag=0.2 $',
     $     idb+105,idb+95, ymin,ymax,xmag,ymag,born)
      CALL pld2mc("CUMU",'!!!xybet21-xybet20, O(alf3) vs O(alf2) ymag=0.2 $',
     $     idb+106,idb+96, ymin,ymax,xmag,ymag,born)
      CALL pld2mc("CUMU",'!!!ybet21-ybet20, O(alf3) vs O(alf2) ymag=0.2 $',
     $     idb+107,idb+97, ymin,ymax,xmag,ymag,born)
      CALL plt1mc("CUMU", '!!!bet30 all, O(alf3) only $',idb+108, ymin,ymax,born)
      ymax =  0.00020d0
      ymin = -ymax
      CALL plt1mc("CUMU", '!!!xxxbet30, O(alf3) only $',idb+109, ymin,ymax,born)
      CALL plt1mc("CUMU", '!!!xxybet30, O(alf3) only $',idb+110, ymin,ymax,born)
      CALL plt1mc("CUMU", '!!!xyybet30, O(alf3) only $',idb+111, ymin,ymax,born)
*==============================================================
      END  !third_ord



      SUBROUTINE fig_pubs
*-------------------------------------------------------------------
* gmage figbet-pubs
* sig(vmax) O(alf2)EEX, publication eps-figs
*-------------------------------------------------------------------
      IMPLICIT NONE
      COMMON / inout  / ninp,nout
      INTEGER           ninp,nout
      DOUBLE PRECISION Born,ymax,ymin,xmag,ymag
      INTEGER          idb,KeyDis
      CHARACTER*60  TeXfile

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_lin   =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      idb = 50000
*==============================================================
*=====================O(alf2)==================================
*==============================================================
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01
      YMAG =  0.01

      TeXfile   = 'flat_total.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(302302,"XCHI2",8302302,' O(alf2) sig(vmax) TOTAL$',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2) TOTAL, sig(xmax)/Born   , YMAG=.01$",
     $     IDB+73, 8302302,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt0xbt0.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(320320,"XCHI2",8320320,' O(alf2) bt0xbt0  $',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0xbt0, YMAG=.01$",
     $     IDB+90, 8320320,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt1xbt0.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(321320,"XCHI2",8321320,' O(alf2) bt1xbt0  $',IDB+72) !
      YMAG =  0.1
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt1*bt0, YMAG=.1$",
     $     IDB+93, 8321320,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt0xbt1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(320321,"XCHI2",8320321,' O(alf2) bt0xbt1  $',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0*bt1, YMAG=.1$",
     $     IDB+94, 8320321,   YMIN,YMAX, XMAG,YMAG, born)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt2xbt0.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      YMAG =  1.
      CALL Semalib_VVplot(322320,"XCHI2",8322320,' O(alf2) bt2xbt0  $',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt2*bt0, YMAG=1$",
     $     IDB+95, 8322320,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt0xbt2.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(320322,"XCHI2",8320322,' O(alf2) bt0xbt2  $',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt0*bt2, YMAG=1$",
     $     IDB+97, 8320322,   YMIN,YMAX, XMAG,YMAG, born)
      CALL GLK_PlEnd

      TeXfile   = 'flat_bt1xbt1.txp'
      CALL GLK_PlInitialize(2,TeXfile)
      CALL Semalib_VVplot(321321,"XCHI2",8321321,' O(alf2) bt1xbt1  $',IDB+72) !
      CALL PLDIFF("CUMU","O(alf2), sig(xmax)/born  bt1*bt1, YMAG=1$",
     $     IDB+96, 8321321,   YMIN,YMAX, XMAG,YMAG, born)
      CALL GLK_PlEnd

*==============================================================
      END  !fig_publ
