*/////////////////////////////////////////////////////////////////////
*//   Final state plots
*//   Technical test of beta's MC versus SAnalytical
*//
*//   make figfin-dvi
*//   make figfin-ps
*//
*/////////////////////////////////////////////////////////////////////
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      SAVE
*---
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*---
      INTEGER   imax
      PARAMETER(imax=10000)
      REAL*8    xpar(imax)
*---
      INTEGER    lint
      CHARACTER*60  Tesnam, TeXfile, Dname
      CHARACTER*60  Hname, DumpFile
*------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'figfin'
      TeXfile   = 'figfin.tex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

* current, running
*      Dname  = '../fin140/fin140.input'  ! current
*      Hname  = '../fin140/pro.hst'       ! current

* Exercise with flat x-section
      Dname  = '../fin140/fin140.input.keyzet-2.5M'  !
      Hname  = '../fin140/pro.hst.keyzet-2.5M'       !

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
      CALL figchi_lin
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


      SUBROUTINE figchi_lin
*-------------------------------------------------------------------
* sig(ln(vmax))
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      SAVE

      WRITE(nout,*) ' ====================================='
      WRITE(nout,*) ' ========= figchi_lin   =============='
      WRITE(nout,*) ' ====================================='

      CALL Semalib_GetBorn(Born)
      ida  = 30000
      idb  = 50000
*==============================================================
*=========================O(alf0)==============================
*==============================================================
*--- O(alf1) total
      KEYDIS=     300
      CALL Semalib_VVplot(KeyDis,"UCHI2", 8300,' O(alf1) sig(ln(vmax)) $',IDB+301)
      ymax =  0.012d0
      ymax =  2                 !!!
      ymin = -ymax
      xmag =  1
      ymag =  0.01
      ymag =  1                 !!!
c      CALL PLDIFF("CUMU","GPS final O(alf0), sig(vmax)/born ,YMAG=.01$",
c     $     idb+201, 8300,   ymin,ymax, xmag,ymag, born)
      CALL PLDIFF("CUMU","final O(alf0) TOTAL, sig(vmax)/born ,YMAG=.01$",
     $     idb+300, 8300,   ymin,ymax, xmag,ymag, born)
*==============================================================
*=========================O(alf1)==============================
*==============================================================
*--- O(alf1) total
      KEYDIS=     301
      CALL Semalib_VVplot(KeyDis,"UCHI2", 8301,' O(alf1) sig(ln(vmax)) $',IDB+301)
*--- O(alf1) bt0xbt0 etc...
      KEYDIS=     310
      CALL Semalib_VVplot(KeyDis,"UCHI2", 8310,' O(alf1) bt0xbt0  $',IDB+301)
      KEYDIS=     311
      CALL Semalib_VVplot(KeyDis,"UCHI2", 8311,' O(alf1) bt1xbt0  $',IDB+301)
*==============================================================
* MC-SAN
      ymax =  0.012d0
      ymin = -0.012d0
      xmag =  1
      ymag =  0.01
c      CALL PLDIFF("CUMU","GPS final O(alf1) , sig(vmax)/born ,YMAG=.01$",
c     $     idb+202, 8301,   ymin,ymax, xmag,ymag, born)
      CALL PLDIFF("CUMU","final O(alf1) TOTAL, sig(vmax)/born ,YMAG=.01$",
     $     idb+301, 8301,   ymin,ymax, xmag,ymag, born)
      CALL PLDIFF("CUMU","final O(alf1), sig(xmax)/Born bt0, YMAG=.01$",
     $     idb+310, 8310,   ymin,ymax, xmag,ymag, born)
      ymag =  0.1
      CALL PLDIFF("CUMU","final O(alf1), sig(xmax)/Born bt1, YMAG=.1$",
     $     idb+311, 8311,   ymin,ymax, xmag,ymag, born)
*==============================================================
*=====================O(alf2)==================================
*==============================================================
*--- O(alf2) total
      KEYDIS=     302
      CALL Semalib_VVplot(KeyDis,"UCHI2", 8302,'O(alf2) sig(vmax) TOTAL$',IDB+301)
*--- O(alf2)  bt0xbt0 etc...
      KEYDIS=     320
      CALL Semalib_VVplot(KeyDis,"UCHI2",8320,' O(alf2) bt0xbt0  $',IDB+301)
      KEYDIS=     321
      CALL Semalib_VVplot(KeyDis,"UCHI2",8321,' O(alf2) bt1xbt0  $',IDB+301)
      KEYDIS=     322
      CALL Semalib_VVplot(KeyDis,"UCHI2",8322,' O(alf2) bt2xbt0  $',IDB+301)
*==============================================================
* MC-SAN
      YMAX =  0.012D0
      YMIN = -0.012D0
      XMAG =  1
      YMAG =  0.01
      CALL PLDIFF("CUMU","final O(alf2) TOTAL, sig(xmax)/Born  ,YMAG=.01$",
     $     IDB+302, 8302,   YMIN,YMAX, XMAG,YMAG, BORN)
      CALL PLDIFF("CUMU","final O(alf2), sig(xmax)/born  bt0, YMAG=.01$",
     $     IDB+320, 8320,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  0.1
      CALL PLDIFF("CUMU","final O(alf2), sig(xmax)/born  bt1, YMAG=.1$",
     $     IDB+321, 8321,   YMIN,YMAX, XMAG,YMAG, BORN)
      YMAG =  1.
      CALL PLDIFF("CUMU","final O(alf2), sig(xmax)/born  bt2, YMAG=1$",
     $     IDB+322, 8322,   YMIN,YMAX, XMAG,YMAG, BORN)
*==============================================================
      END  !figchi_lin 
