*****************************************************************
*//   alias kmake='make -f KKMakefile'
*
* kmake chi_mcan-ps
*
**----------------------------------------------------------------
      PROGRAM MAIN
*     ***********************************
      IMPLICIT NONE
*
      INCLUDE 'chi.h'
*
      CHARACTER*60  Tesnam, Dname
      CHARACTER*60  Hname, DumpFile
*--------------------------------------------------------------------
      ninp=  5
      nout= 16

      Tesnam    = 'figchi'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

*** current 200GeV
**      Dname  = '../mix200/mix200.input' ! current
**      Hname  = '../mix200/pro.hst'      ! current
*
      Dname  = '../E189GeV/mix189.input' ! latest
      Hname  = '../E189GeV/mix189.hst'   ! vmax=1
cc      Hname  = '../E189GeV/mix189.hst.999'   ! vmax=0.999 

* newest
c      Dname  = '../mix200/mix200.input.375M.xenph=1.25' ! Jan 99
c      Hname  = '../mix200/pro.hst.375M.xenph=1.25'      ! Jan 99
*
c      Dname  = '../mix200/mix200.input__321M' ! Jan 99
c      Hname  = '../mix200/pro.hst__321M'      ! Jan 99

c      Dname  = '../mix200/mix200.input.320M' ! June 98
c      Hname  = '../mix200/pro.hst.320M'      ! June 98

c      Dname  = '../mix200/mix200.input.vmax=1_70M' ! Dec 97
c      Hname  = '../mix200/pro.hst.vmax=1_70M'      ! Dec 97

* principal  200GeV
c      Dname  = '../mix200/mix200.input.106M' ! Nov 97
c      Hname  = '../mix200/pro.hst.106M'      ! Nov 97

*=====================================
* Read data, the same as in MC run
      CALL ReaDat(Dname,imax,xpar)
      CALL Semalib_Initialize(xpar)
* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
*=====================================================
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
      CALL Prepare
*=========================================================================
* Write all histograms into dump file, for control/debug
      DumpFile = './chi.hst'
      CALL GLK_WriteFile(DumpFile)
*=================================
      CLOSE(nout)
      END

      SUBROUTINE Prepare
*     ******************
      IMPLICIT NONE
      SAVE
      INCLUDE 'chi.h'
*------------
      CHARACTER*80 title
      REAL*8    born,vmin,vmax,bornan
*------------------------------------------------------------------------------

      CALL Semalib_GetBorn(Born)
      WRITE(   *,*) 'figchi:: Born= ',born
      WRITE(nout,*) 'figchi:: Born= ',born
*===================================================================
*              Analytical
*===================================================================
      KeyDis=     400300
      CALL Semalib_VVplot(KeyDis,"XCHI2", iO0nll,' O(alf0) ini+fin $',idb+73)
*
      KeyDis=     300300
      CALL Semalib_VVplot(KeyDis,"XCHI2", iO0tot,' O(alf0) ini+fin $',idb+73)
*
      KeyDis=     301301
      CALL Semalib_VVplot(KeyDis,"XCHI2", iO1tot,' O(alf1) ini+fin $',idb+73)
*
      KeyDis=     302302
      CALL Semalib_VVplot(KeyDis,"XCHI2", iO2tot,' O(alf2) ini+fin $',idb+73)
*
      KeyDis=     304303        ! ISR with O(alf3)LL
      KeyDis=     304302        ! ISR with O(alf3)LL exp(3/4*gam)
      KeyDis=     305302        ! The best ISR,  O(alf3)LL + O(alf2)NLL
      CALL Semalib_VVplot(KeyDis,"XCHI2", iO3best,' O(alf3) ini+fin best $',idb+73)
*
      CALL GLK_Operat( iO0nll,'+',iO0nll,   iO0nll,  0d0, 1/born)
      CALL GLK_Operat( iO0tot,'+',iO0tot,   iO0tot,  0d0, 1/born)
      CALL GLK_Operat( iO1tot,'+',iO1tot,   iO1tot,  0d0, 1/born)
      CALL GLK_Operat( iO2tot,'+',iO2tot,   iO2tot,  0d0, 1/born)
      CALL GLK_Operat(iO3best,'+',iO3best,  iO3best, 0d0, 1/born)
*===================================================================
*               Monte Carlo
*===================================================================
      CALL GLK_CumHis(IdGenYFS3, idb+71,mO0tot)
      CALL GLK_CumHis(IdGenYFS3, idb+72,mO1tot)
      CALL GLK_CumHis(IdGenYFS3, idb+73,mO2tot)
      CALL GLK_CumHis(IdGenYFS3, idb+74,mO3tot)
      CALL GLK_CumHis(IdGenYFS3, idb+76,mO2vO1)
      CALL GLK_CumHis(IdGenYFS3, idb+77,mO3vO2)
*
      CALL GLK_Operat( mO0tot,'+',mO0tot, mO0tot, 0d0, 1/born)
      CALL GLK_Operat( mO1tot,'+',mO1tot, mO1tot, 0d0, 1/born)
      CALL GLK_Operat( mO2tot,'+',mO2tot, mO2tot, 0d0, 1/born)
      CALL GLK_Operat( mO3tot,'+',mO3tot, mO3tot, 0d0, 1/born)
      CALL GLK_Operat( mO2vO1,'+',mO2vO1, mO2vO1, 0d0, 1/born)
      CALL GLK_Operat( mO3vO2,'+',mO3vO2, mO3vO2, 0d0, 1/born)
* Divide diffrences by the total
      CALL GLK_Operat( mO2vO1,'/',mO2tot, mO2vO1, 1d0, 1d0)
      CALL GLK_Operat( mO3vO2,'/',mO3tot, mO3vO2, 1d0, 1d0)
*===================================================================
*               Differences
*===================================================================
* MC-SAN
      CALL GLK_Operat( mO0tot,'-',iO0nll, iO0tech, 1d0, 1d0)
      CALL GLK_Operat( mO0tot,'-',iO0tot, iO0dif,  1d0, 1d0)
      CALL GLK_Operat( mO1tot,'-',iO1tot, iO1dif,  1d0, 1d0)
      CALL GLK_Operat( mO2tot,'-',iO2tot, iO2dif,  1d0, 1d0)
*
      CALL GLK_Operat(iO0tech,'/',iO0tot, iO0tech, 1d0, 1d0)
      CALL GLK_Operat(iO0dif ,'/',iO0tot, iO0dif,  1d0, 1d0)
      CALL GLK_Operat(iO1dif ,'/',iO1tot, iO1dif,  1d0, 1d0)
      CALL GLK_Operat(iO2dif ,'/',iO2tot, iO2dif,  1d0, 1d0)

* O(alf2)-O(alf1)
      CALL GLK_Operat(iO2tot, '-',iO1tot, iO2mO1,  1d0,1d0)
      CALL GLK_Operat(mO2tot, '-',mO1tot, mO2mO1,  1d0,1d0)
*
      CALL GLK_Operat(iO2mO1, '/',iO2tot, iO2mO1,  1d0, 1d0)
      CALL GLK_Operat(mO2mO1, '/',mO2tot, mO2mO1,  1d0, 1d0)
      CALL GLK_Print(iO2mO1)
      CALL GLK_Print(mO2mO1)
* MC-SAN for the difference
      CALL GLK_Operat(mO2mO1,'-',iO2mO1,mO2dO1,1d0,1d0)
      CALL GLK_Operat(mO3tot,'-',iO3best,mO3dan,1d0,1d0) ! O(alf3) MC-San
* 
      CALL GLK_Operat(mO2mO1,'/',iO2tot, mO2dO1,  1d0,1d0)
      CALL GLK_Operat(mO3dan,'/',iO3best,mO3dan,  1d0,1d0)
*===================================================================
      END

