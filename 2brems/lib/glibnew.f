
C======================================================================
C======================= G L I B K  ===================================
C==================General Library of utilities========================
C===========It is imilar but not identical to HBOOK and HPLOT==========
C======================================================================
C   
C                      Version:    1.05
C              Last correction:    JULY 1992
C
C
C  Installation remarks: 
C  (1) printing backslash character depends on F77 compilator,
C      user may need to modify definition of BS variable in HPLCAP
C
C  Usage of the program:
C  (1) In most cases names and meanings of programs and their 
C      parameters is the same as in original CERN libraries HBOOK
C  (2) Unlike to original HBOOK and HPLOT, all floating parameters 
C      of the programs are in double precision!
C  (3) GLIBK stores histograms in double precision and always with
C      errors. REAL*8 storage is essential for 10**7 events statistics!
C  (4) Output from GLIBK is a picture recorded as regular a LaTeX file 
C      with frame and curves/histograms, it is easy to change fonts
C      add captions, merge plots, etc. by normal ediding. Finally,
C      picture may be inserted in any place into LaTeX source of the
C      article.
C
C  ********************************************************************
C  *  History of the program:                                         *
C  *  MINI-HBOOK writen by S. Jadach, Rutherford Lab. 1976            *
C  *  Rewritten December 1989 (S.J.)                                  *
C  *  Version with DOUBLE PRECISION ARGUMENTS ONLY!  and SAVE         *
C  *  Subrogram names start with G instead of H letter!               *
C  *  Entries:   Obligatory:  GLIMIT                                  *
C  *             Optional: see table below                            *
C  *  non-user subprograms in brackets                                *
C  ********************************************************************
C    SUBR/FUNC  1 PAR. 2 PAR. 3 PAR. 4 PAR. 5 PAR. 6 PAR.       
C  ====================================================================
*     (GINIT)   ----   ----    ----   ----   ----   ----        
*      GI       INT    INT     ----   ----   ----   ----        
*      GIE      INT    INT     ----   ----   ----   ----        
*      GF1      INT    DBL     DBL    ----   ----   ----        
*      GFILL    INT    DBL     DBL    DBL    ----   ----        
*      GBOOK1   INT    CHR*80  INT    DBL    DBL    ----  
*     (GOPTOU)  INT    INT     INT    INT    INT     INT
* (L.F. GEXIST) INT    -----  ------  ----   ----   ----        
*      GIDOPT   INT    CHR*4   -----  ----   ----   ----        
*      GBFUN1   INT    CHR*80   INT   DBL    DBL  DP-FUNC       
*      GIDOPT   INT    CHR*4   -----  ----   ----   ----        
*      GBOOK2   INT    CHR*80   INT   DBL    DBL     INT   DBL   DBL
*      GISTDO     ---   ----   ----   ----   ----   ----        
*      GOUTPU   INT     ----   ----   ----   ----   ----        
*      GPRINT   INT     ----   ----   ----   ----   ----        
*      GOPERA   INT    CHR*1   INT    INT    DBL    DBL         
*      GINBO1   INT    CHR*8   INT    DBL    DBL    ----        
*      GUNPAK   INT    DBL(*) CHR*(*) INT    ---    ----        
*      GPAK     INT    DBL(*)  ----   ----   ---    ----        
*      GPAKE    INT    DBL(*)  ----   ----   ---    ----       
*      GRANG1   INT    DBL     DBL    ----   ---    ----        
*      GINBO2   INT    INT     DBL    DBL    INT    DBL   DBL      
*      GMAXIM   INT    DBL     ----   ----   ---    ----        
*      GMINIM   INT    DBL     ----   ----   ---    ----        
*      GRESET   INT   CHR*(*)  ----   ----   ---    ----        
*      GDELET   INT     ----   ----   ----   ----   ----        
*      GLIMIT   INT     ----   ----   ----   ----   ----        
*     (COPCH)   CHR*80 CHR*80  ----   ----   ----   ----        
* (F. JADRES)   INT     ----   ----   ----   ----   ----        
*      GRFILE   INT   CHR*(*) CHR*(*) ----   ----   ----        
*      GROUT    INT    INT    CHR*8   ----   ----   ----        
*      GRIN     INT    INT     INT    ----   ----   ----        
*      GREND   CHR*(*) ----    ----   ----   ----   ----        
C  *******************  HPLOT entries ******************
*      GPLINT   INT    ----    ----   ----   ----   ----        
*      GPLCAP   INT    ----    ----   ----   ----   ----        
*      GPLEND   ----   ----    ----   ----   ----   ----        
*      GPLOT    INT    CHR*1   CHR*1   INT   ----   ----        
*     (LFRAM1)  INT      INT     INT  ----   ----   ----        
*     (SAXIX)   INT      DBL     DBL   INT    DBL   ----        
*     (SAXIY)   INT      DBL     DBL   INT    DBL   ----        
*     (PLHIST)  INT      INT     DBL   DBL    INT    INT        
*     (PLHIS2)  INT      INT     DBL   DBL    INT    INT        
*     (PLCIRC)  INT      INT     INT   DBL    DBL    DBL        
*     (APROF)   DBL      INT     DBL  ----   ----   ----        
*      GPLSET   INT      DBL    ----  ----   ----   ----        
*      GPLTIT   INT    CHR*80   ----  ----   ----   ----        
C  *******************  WMONIT entries ******************
*      GMONIT   INT ???
C  *******************************************************************
C                         END OF TABLE        
C  *******************************************************************
*          Map of memory for single histogram
*          ----------------------------------
*  (1-7) Header
*  ist +1   mark      9999999999999
*  ist +2   mark      9d12 + id*10 + 9
*  ist +3   iflag1    9d12 + iflag1*10 +9
*  ist +4   iflag2    9d12 + iflag2*10 +9
*  ist +5   scamin    minimum y-scale
*  ist +6   scamax    maximum y-scale
*  ist +7   reserve   9999999999999
*  One dimensional histogram            Two dimensional histog.
*  -------------------------            ----------------------
*  (8-11) Binning information           (8-15) Binning information
*  ist2 +1    NCHX                          ist2 +5   NCHY
*  ist2 +2      XL                          ist2 +6     YL
*  ist2 +3      XU                          ist2 +7     YU
*  ist2 +4   FACTX                          ist2 +8  FACTY
*
*  (12-24) Under/over-flow average x    (16-24)
*  ist3 +1   Underflow                     All nine combinations
*  ist3 +2   Normal                        (U,N,O) x (U,N,O)
*  ist3 +3   Overerflow                    sum wt only (no errors)
*  ist3 +4   U  sum w**2
*  ist3 +5   N  sum w**2
*  ist3 +6   O  sum w**2
*  ist3 +7   Sum 1
*  ist3 +8   Sum wt*x
*  ist3 +9   Sum wt*x*x
*  ist3 +10  nevzer    (gmonit)
*  ist3 +11  nevove    (gmonit)
*  ist3 +12  nevacc    (gmonit)
*  ist3 +13  maxwt     (gmonit)
*  -----------------------Bin content
*  (25 to 24+2*nchx)                     (25 to 24 +nchx*nchy)
*     sum wt and sum wt**2            sum wt only (no errors)
*  ----------------------------------------------------------------
      subroutine ginit
*     ****************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      data init/0/
      save /gind/
      save init
*
      if(init.ne.0) return
      init=1
c this is version version number
      nvrs=105
c default output unit
      nout=16
      lenmax=0
      length=0
      do 100 i=1,idmx
      do 40 k=1,3
   40 index(i,k)=0
      do 50 k=1,80
   50 titlc(i)(k:k)=' '
  100 continue
      end

      function gi(id,ib)
*     ******************
C getting out bin content
C S.J. 18-Nov. 90
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      save idmem,nch,lact,ist,ist2,ist3
      data idmem / -1256765/
c
      if(id.eq.idmem) goto 100
      idmem=id
c some checks, not repeated if id the same as previously
      lact=jadres(id)
      if(lact.eq.0) then
        write(nout,*) ' gi: nonexisting histo id=',id
        gi= 0d0
        return
      endif
      ist  = index(lact,2)
      ist2 = ist+7
      ist3 = ist+11
c checking if histo is of proper type
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.ne.1) then
        write(nout,*) ' gi: 1-dim histos only !!! id=',id
        gi= 0d0
        return
      endif
  100 continue
      nch  = nint(b(ist2+1))
      if(ib.eq.0) then
c underflow
         gi=   b(ist3 +1)
      elseif(ib.ge.1.and.ib.le.nch) then
c normal bin
         gi=   b(ist +nbuf+ib)
      elseif(ib.eq.nch+1) then
c overflow
         gi=   b(ist3 +3)
      else
c abnormal exit
         write(nout,*) ' gi: wrong binning id,ib=',id,ib
         gi=0d0
         return
      endif
      end

      function  gie(id,ib)
*     ********************
c getting out error of the bin
c s.j. 18-nov. 90
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      save idmem,nch,lact,ist,ist2,ist3
      data idmem / -1256765/
c
      if(id.eq.idmem) goto 100
      idmem=id
c some checks, not repeated if id the same as previously
      lact=jadres(id)
      if(lact.eq.0) then
        write(nout,*) ' gie: nonexisting histo id=',id
        gie= 0d0
        return
      endif
      ist  = index(lact,2)
      ist2 = ist+7
      ist3 = ist+11
c checking if histo is of proper type
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.ne.1) then
        write(nout,*) ' gie: 1-dim histos only !!! id=',id
        gie= 0d0
        return
      endif
  100 continue
      nch  = b(ist2+1)
      if(ib.eq.0) then
c underflow
         gie=   dsqrt( dabs(b(ist3 +4)))
      elseif(ib.ge.1.and.ib.le.nch) then
c...normal bin, error content
         gie=   dsqrt( dabs(b(ist+nbuf+nch+ib)) )
      elseif(ib.eq.nch+1) then
c overflow
         gie=   dsqrt( dabs(b(ist3 +6)))
      else
c abnormal exit
         write(nout,*) ' gie: wrong binning id, ib=',id,ib
         gie=0d0
         return
      endif
      end

      subroutine gf1(id,x,wtw)
*     ************************
c recommended fast filling 1-dim. histogram
c s.j. 18 nov. 90
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
c exit for non-existig histo
      if(lact.eq.0)  return
      ist  = index(lact,2)
      ist2 = ist+7
      ist3 = ist+11
c one-dim. histo only
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.ne.1) return
      xx= x
      wt= wtw
      index(lact,3)=index(lact,3)+1
c for average x
      b(ist3 +7)  =b(ist3 +7)   +1
      b(ist3 +8)  =b(ist3 +8)  +wt*xx
      b(ist3 +9)  =b(ist3 +9)  +wt*xx*xx
c filling bins
      nchx  =b(ist2 +1)
      xl    =b(ist2 +2)
      factx =b(ist2 +4)
      kx = (xx-xl)*factx+1d0
      if(kx.lt.1) then
c underflow
         b(ist3 +1)    = b(ist3 +1)         +wt
         b(ist3 +4)    = b(ist3 +4)         +wt*wt
      elseif(kx.gt.nchx) then
c overflow
         b(ist3 +3)    = b(ist3 +3)         +wt
         b(ist3 +6)    = b(ist3 +6)         +wt*wt
      else
c normal bin
         b(ist3 +2)    = b(ist3 +2)         +wt
         b(ist +nbuf+kx) = b(ist+nbuf+kx)   +wt
c error bin
         b(ist3 +5)    = b(ist3 +5)         +wt*wt
         b(ist +nbuf+nchx+kx) = b(ist+nbuf+nchx+kx)   +wt**2
      endif
      end

      subroutine gfill(id,x,y,wtw)
*     ****************************
c this routine not finished, 1-dim only!
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0)  return
      ist  = index(lact,2)
c one-dim. histo 
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.1) then
c...one-dim. histogram
        call gf1(id,x,wtw)
        return
      endif
c...two-dim. scattergram, no errors!
      ist2 = ist+7
      ist3 = ist+15
      xx= x
      yy= y
      wt= wtw
      index(lact,3)=index(lact,3)+1
c x-axis
      nchx  =b(ist2 +1)
      xl    =b(ist2 +2)
      factx =b(ist2 +4)
      kx=(xx-xl)*factx+1d0
      lx=2
      if(kx.lt.1)     lx=1
      if(kx.gt.nchx)  lx=3
      l     = ist+34  +lx
      b(l)  = b(l)    +wt
      k     = ist+nbuf2  +kx
      if(lx.eq.2) b(k)  =b(k)  +wt
      k2    = ist+nbuf2  +nchx+kx
      if(lx.eq.2) b(k2) =b(k2) +wt**2
c y-axix
      nchy  =b(ist2 +5)
      yl    =b(ist2 +6)
      facty =b(ist2 +8)
      ky=(yy-yl)*facty+1d0
      ly=2
      if(ky.lt.1)    ly=1
      if(ky.gt.nchy) ly=3
c under/over-flow
      l = ist3  +lx +3*(ly-1)
      b(l) =b(l)+wt
c regular bin
      k = ist+nbuf2 +kx +nchx*(ky-1)
      if(lx.eq.2.and.ly.eq.2) b(k)=b(k)+wt
      end

      subroutine gbook1(id,title,nnchx,xxl,xxu)
*     *****************************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      character*80 title
      logical gexist
c
      call ginit
      if(gexist(id)) goto 900
      ist=length
      lact=jadres(0)
c the case of no free entry in the index
      if(lact.eq.0) goto 901
      index(lact,1)=id
      index(lact,2)=length
      index(lact,3)=0
c -------
      call copch(title,titlc(lact))
      nchx =nnchx
      xl   =xxl
      xu   =xxu
c ---------- title and bin content ----------
      lengt2 = length +2*nchx +nbuf+1
      if(lengt2.ge.lenmax) goto 902
      do 10 j=length+1,lengt2+1
  10  b(j) = 0d0
      length=lengt2
c... default flags
      ioplog   = 1
      iopsla   = 1
      ioperb   = 1
      iopsc1   = 1
      iopsc2   = 1
      iflag1   = 
     $ ioplog+10*iopsla+100*ioperb+1000*iopsc1+10000*iopsc2
      ityphi   = 1
      iflag2   = ityphi
C examples of decoding flags 
c      iflag1   = nint(b(ist+3)-9d0-9d12)/10
c      ioplog = mod(iflag1,10)
c      iopsla = mod(iflag1,100)/10
c      ioperb = mod(iflag1,1000)/100
c      iopsc1 = mod(iflag1,10000)/1000
c      iopsc2 = mod(iflag1,100000)/10000
c      iflag2   = nint(b(ist+4)-9d0-9d12)/10
c      ityphi = mod(iflag2,10)
c--------- buffer -----------------
c header
      b(ist +1)  = 9999999999999d0
      b(ist +2)  = 9d12 +id*10 +9d0
      b(ist +3)  = 9d12 + iflag1*10 +9d0
      b(ist +4)  = 9d12 + iflag2*10 +9d0
c dummy vertical scale
      b(ist +5)  =  -100d0
      b(ist +6)  =   100d0
c reserve
      b(ist +7)  = 9999999999999d0
c information on binning
      ist2       = ist+7
      b(ist2 +1) = nchx
      b(ist2 +2) = xl
      b(ist2 +3) = xu
      ddx = xu-xl
      if(ddx.eq.0d0) goto 903
      b(ist2 +4) = float(nchx)/ddx
c under/over-flow etc.
      ist3       = ist+11
      do 100  j=1,13
 100  b(ist3 +j)=0d0
c
      return
 900  write(nout,*) ' gbook1: histo already exists!!!! id=',id
      return
 901  write(nout,*) ' gbook1: to many histos !!!!!, id=',id
      stop
 902  write(nout,*) ' gbook1: to litle storage!!!!, lenmax=',lenmax
      stop
 903  write(nout,*) ' gbook1: xl=xu, id=',id
      stop
      end

      logical function gexist(id)
c     ***************************
c this function is true when id  exists !!!! 
c     ***************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      gexist = lact.ne.0
      end

      subroutine goptou(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
c     ********************************************************
c decoding option flags
c     **********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/

      lact=jadres(id)
      if(lact.eq.0) return
      ist=index(lact,2)
c decoding flags 
      iflag1   = nint(b(ist+3)-9d0-9d12)/10
      ioplog = mod(iflag1,10)
      iopsla = mod(iflag1,100)/10
      ioperb = mod(iflag1,1000)/100
      iopsc1 = mod(iflag1,10000)/1000
      iopsc2 = mod(iflag1,100000)/10000
      end

      subroutine gidopt(id,ch)
c     ************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      character*4 ch
c
      lact=jadres(id)
      if(lact.eq.0) return
      ist=index(lact,2)
C decoding flags 
      call goptou(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      if(ch.eq.      'LOGY'  ) then
c log scale for print
        ioplog = 2 
      elseif(ch.eq.  'ERRO'  ) then
C errors in printing/plotting
       ioperb  = 2
      elseif(ch.eq.  'SLAN'  ) then
c slanted line in plotting
       iopsla  = 2
      elseif(ch.eq.  'YMIN'  ) then
       iopsc1  = 2
      elseif(ch.eq.  'YMAX'  ) then
       iopsc2  = 2
      ENDIF
c encoding back
      iflag1   = 
     $ ioplog+10*iopsla+100*ioperb+1000*iopsc1+10000*iopsc2
      b(ist+3) = 9d12 + iflag1*10 +9d0
      end

      subroutine gbfun1(id,title,nchx,xmin,xmax,func)
c     ***********************************************
c ...fills histogram with function func(x)
c     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /gind/
      dimension yy(200)
      external func
      character*80 title
      logical gexist
c
      call ginit
      if(gexist(id)) goto 900
 15   xl=xmin
      xu=xmax
      call gbook1(id,title,nchx,xl,xu)
c...slanted line in plotting
      call gidopt(id,'SLAN')
      if(nchx.gt.200) goto 901
      do 20 ib=1,nchx
      x= xmin +(xmax-xmin)/nchx*(ib-0.5d0)
      yy(ib) = func(x)
   20 continue
      call gpak(id,yy)
      return
 900  write(nout,*) ' +++gbfun1: already exists id=',id
      write(6   ,*) ' +++gbfun1: already exists id=',id      
      call gdelet(id)
      go to 15
 901  write(nout,*) ' +++gbfun1: to many bins'
      end

      SUBROUTINE GBOOK2(ID,TITLE,NCHX,XL,XU,NCHY,YL,YU)
*     *************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER( IDMX=400,NBUF=24,NBUF2=24)
      COMMON / Cglib / B(20000)
      COMMON /GIND/ NVRS,NOUT,LENMAX,LENGTH,INDEX(IDMX,3),TITLC(IDMX)
      CHARACTER*80 TITLC
      save /cglib/,/gind/
      CHARACTER*80 TITLE
      LOGICAL GEXIST
c
      CALL GINIT
      IF(GEXIST(ID)) GOTO 900
      ist=length
      LACT=JADRES(0)
      IF(LACT.EQ.0) GOTO 901
      index(LACT,1)=ID
      index(LACT,2)=length
      CALL COPCH(TITLE,TITLC(LACT))
      nnchx=NCHX
      nnchy=NCHY
      LENGT2 = LENGTH  +44+nnchx*nnchy
      IF(LENGT2.GE.LENMAX) GOTO 902
      DO 10 J=LENGTH+1,LENGT2+1
   10 B(J) = 0D0
      LENGTH=LENGT2
      B(ist+1)=nnchx
      B(ist+2)=XL
      B(ist+3)=XU
      B(ist+4)=float(nnchx)/(b(ist+3)-b(ist+2))
      B(ist+5)=nnchy
      B(ist+6)=YL
      B(ist+7)=YU
      B(ist+8)=float(nnchy)/(b(ist+7)-b(ist+6))
      RETURN
  900 WRITE(NOUT,*) ' GBOOK2: HISTO ALREADY EXISTS!!!! ID=',ID
      RETURN
  901 WRITE(NOUT,*) ' GBOOK2: TO MANY HISTOS !!!!!',LACT
      STOP
  902 WRITE(NOUT,*) ' GBOOK2: TO LITLE STORAGE!!!!',LENMAX
      STOP
      END

      subroutine gistdo
*     *****************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /gind/
      do 10 i=1,idmx
      id=index(i,1)
      if(id.gt.0) call gprint(id)
   10 continue
      end

      subroutine goutpu(ilun)
*     ***********************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /gind/
      call ginit
      nout=ilun
      end


      subroutine gprint(id)
*     *********************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      character*1 line(105),lchr(22),lb,lx,li,l0
      logical llg
      data lb,lx,li,l0 /' ','X','I','0'/
      data lchr/' ','1','2','3','4','5','6','7','8','9',
     $      'A','B','C','D','E','F','G','H','I','J','K','*'/
      save lb,lx,li,l0,lchr

      lact=jadres(id)
      if(lact.eq.0) goto 900
      ist  = index(lact,2)
      ist2 = ist+7
      ist3 = ist+11

      call goptou(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      ker    =  ioperb-1
      lmx = 67
      if(ker.eq.1) lmx=54
      nent=index(lact,3)
      if(nent.eq.0) goto 901
      write(nout,1000) id,titlc(lact)
 1000 FORMAT('1',/,1X,I6,10X,A)
c
c one-dim. histo 
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.ne.1) goto 200
      nchx =   b(ist2 +1)
      xl   =   b(ist2 +2)
      dx   =  (  b(ist2 +3)-b(ist2 +2)  )/float(nchx)
c fixing vertical scale
      istr=ist+nbuf+1
      bmax = b(istr)
      bmin = b(istr)
      do 15 ibn=istr,istr+nchx-1
      bmax = max(bmax,b(ibn))
      bmin = min(bmin,b(ibn))
  15  continue
      if(bmin.eq.bmax) goto 901
      if(iopsc1.eq.2) bmin=b(ist +5)
      if(iopsc2.eq.2) bmax=b(ist +6)
c
      llg=ioplog.eq.2
      if(llg.and.bmin.le.0d0) bmin=bmax/10000.d0
c
      deltb = bmax-bmin
      if(deltb.eq.0d0) goto 902
      fact  = (lmx-1)/deltb
      kzer  = -bmin*fact+1.00001d0
      if(llg) fact=(lmx-1)/(log(bmax)-log(bmin))
      if(llg) kzer=-log(bmin)*fact+1.00001d0
c
      undf = b(ist3 +1)
      ovef = b(ist3 +3)
      avex = 0d0
      sum  = b(ist3 +8)
      if(nent.ne.0) avex = sum/nent
      write(nout,'(4a15      )')  'nent','sum','bmin','bmax'
      write(nout,'(i15,3e15.5)')   nent,  sum,  bmin,  bmax
      write(nout,'(4a15  )')      'undf','ovef','avex'
      write(nout,'(4e15.5)')       undf,  ovef,  avex
c
      if(llg) write(nout,1105)
 1105 format(35x,17hlogarithmic scale)
c
      kzer=max0(kzer,0)
      kzer=min0(kzer,lmx)
      xlow=xl
      do 100 k=1,nchx
c first fill with blanks
      do  45 j=1,105
   45 line(j)  =lb
c then fill upper and lower boundry
      line(1)  =li
      line(lmx)=li
      ind=istr+k-1
      bind=b(ind)
      bind= max(bind,bmin)
      bind= min(bind,bmax)
      kros=(bind-bmin)*fact+1.0001d0
      if(llg) kros=log(bind/bmin)*fact+1.0001d0
      k2=max0(kros,kzer)
      k2=min0(lmx,max0(1,k2))
      k1=min0(kros,kzer)
      k1=min0(lmx,max0(1,k1))
      do 50 j=k1,k2
   50 line(j)=lx
      line(kzer)=l0
      z=b(ind)
      if(ker.ne.1) then
        write(nout,"(a, f7.4,  a, d14.6,  132a1)") 
     $             " ", xlow," ",     z," ",(line(i),i=1,lmx)
      else
        er=dsqrt(dabs(b(ind+nchx)))
        write(nout,"(a,f7.4,  a,d14.6,  a,d14.6, 132a1 )") 
     $             " ",xlow," ",    z," ",   er," ",(line(i),i=1,lmx)
      endif
      xlow=xlow+dx
  100 continue
      return
C------------- two dimensional requires complete restoration!!!----------------
  200 continue
      nchx=B(ist+1)
      nchy=B(ist+5)
      write(nout,2000) (lx,i=1,nchy)
 2000 format(1h ,10x,2hxx,100a1)
      do 300 kx=1,nchx
      do 250 ky=1,nchy
      k=ist +NBUF2 +kx+nchx*(ky-1)
      N=B(K)+1.99999D0
      n=max0(n,1)
      n=min0(n,22)
      if(DABS(b(k)).lt.1D-20) n=1
      line(ky)=lchr(n)
  250 continue
      line(nchy+1)=lx
      i1=nchy+1
      write(nout,2100) (line(i),i=1,i1)
 2100 format(1h ,10x,1hx,100a1)
  300 continue
      write(nout,2000) (lx,i=1,nchy)
      RETURN
  900 WRITE(NOUT,*) ' +++GPRINT: NONEXISTING HISTO',ID
      WRITE(6   ,*) ' +++GPRINT: NONEXISTING HISTO',ID
      RETURN
 901  WRITE(NOUT,*) ' +++GPRINT: NO ENTRIES  HISTO',ID
      WRITE(   6,*) ' +++GPRINT: NO ENTRIES  HISTO',ID
      RETURN
 902  WRITE(NOUT,*) ' +++GPRINT: wrong plotting limits',ID
      WRITE(   6,*) ' +++GPRINT: wrong plotting limits',ID
      END

      subroutine gopera(ida,chr,idb,idc,coef1,coef2)
*     **********************************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      character*80 title
      character*1  chr
c
      lacta=jadres(ida)
      if(lacta.eq.0) return
      ista  = index(lacta,2)
      ista2 = ista+7
      ncha  = b(ista2+1)
c
      lactb =jadres(idb)
      if(lactb.eq.0) return
      istb  = index(lactb,2)
      istb2 = istb+7
      nchb  = b(istb2+1)
      if(nchb.ne.ncha) goto 900
c
      lactc=jadres(idc)
      if(lactc.eq.0) then
c ...if nonexistent, histo idc is here defined
        call ginbo1(ida,title,nchx,xl,xu)
        call gbook1(idc,title,nchx,xl,xu)
        lactc = jadres(idc)
        istc  = index(lactc,2)
c...option copied from ida
        b(istc+ 3)= b(ista +3)
      endif
c...one nominal entry recorded
      index(lactc,3) = 1
c
      istc  =  index(lactc,2)
      istc2 =  istc+7
      nchc  =  b(istc2+1)
c
      if(nchc.ne.ncha) goto 900
      if(ncha.ne.nchb.or.nchb.ne.nchc) goto 900
      do 30 k=1,ncha
      i1 = ista+nbuf+k
      i2 = istb+nbuf+k
      i3 = istc+nbuf+k
      j1 = ista+nbuf+ncha+k
      j2 = istb+nbuf+ncha+k
      j3 = istc+nbuf+ncha+k
      if    (chr.eq.'+')   then
        b(i3) =    coef1*b(i1) +    coef2*b(i2)
        b(j3) = coef1**2*b(j1) + coef2**2*b(j2)
      elseif(chr.eq.'-')   then
        b(i3) = coef1*b(i1) - coef2*b(i2)
        b(j3) = coef1**2*b(j1) + coef2**2*b(j2)
      elseif(chr.eq.'*')   then
        b(j3) = (coef1*coef2)**2
     $          *(b(j1)*b(i2)**2 + b(j2)*b(i1)**2)
        b(i3) = coef1*b(i1) * coef2*b(i2)
      elseif(chr.eq.'/')   then
        if(b(i2).eq.0d0) then
          b(i3) = 0d0
          b(j3) = 0d0
        else
          b(j3) = (coef1/coef2)**2/b(i2)**4
     $          *(b(j1)*b(i2)**2 + b(j2)*b(i1)**2)
          b(i3) = (coef1*b(i1) )/( coef2*b(i2))
        endif
      else
        goto 901
      endif
   30 continue
      return
  900 write(nout,*) '+++++ gopera: non-equal no. bins ',ida,idb,idc
      write(   6,*) '+++++ gopera: non-equal no. bins ',ida,idb,idc
      return
  901 write(nout,*) '+++++ gopera: wrong chr=',chr
      end

      subroutine ginbo1(id,title,nchx,xl,xu)
*     **************************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      character*80 title
c
      lact=jadres(id)
      if(lact.eq.0) return
      ist=index(lact,2)
      ist2   = ist+7
      nchx   = b(ist2 +1)
      xl     = b(ist2 +2)
      xu     = b(ist2 +3)
      title  = titlc(lact)
      end

      subroutine gunpak(id,a,chd1,idum)
*     *********************************
c getting out histogram content (and error)
c chd1= 'ERRO' is nonstandard option (unpack errors)
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      character*(*) chd1
      dimension a(*)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0) goto 900
      ist   = index(lact,2)
      ist2  = ist+7
      nch   = b(ist2 +1)
      local = ist +nbuf
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.2) then
        nchy  = b(ist2+5)
        nch   = nch*nchy
        local = ist+ nbuf2
      endif
      do 10 ib=1,nch
      if(chd1.ne.'ERRO') then
c normal bin
        a(ib) = b(local+ib)
      else
c error content
        if(ityphi.eq.2) goto 901
        a(ib) = dsqrt( dabs(b(local+nch+ib) ))
      endif
   10 continue
      return
 900  write(nout,*) '+++gunpak: nonexisting id=',id
      write(6   ,*) '+++gunpak: nonexisting id=',id
      return
 901  write(nout,*) '+++gunpak: no errors, two-dim, id=',id
      end

      subroutine gpak(id,a)
*     *********************
c getting in histogram content
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      dimension  a(*)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0) goto 900
      ist  = index(lact,2)
      ist2 = ist+7
      nch=b(ist2 +1)
      local = ist+nbuf
c 2-dimens histo alowed
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.2) then
        nchy  = b(ist2+5)
        nch   = nch*nchy
        local = ist+nbuf2
      endif
      do 10 ib=1,nch
   10 b(local +ib) = a(ib)
c one nominal entry recorded
      index(lact,3)  = 1
      return
  900 write(nout,*) '+++gpak: nonexisting id=',id
      write(6   ,*) '+++gpak: nonexisting id=',id
      end

      subroutine gpake(id,a)
*     **********************
c getting in error content
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      dimension  a(*)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0) goto 901
      ist  = index(lact,2)
      ist2 = ist+7
      nch=b(ist2+1)
c 2-dimens histo NOT alowed
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.2) goto 900
      do 10 ib=1,nch
   10 b(ist+nbuf+nch+ib) = a(ib)**2
      return
  900 write(nout,*) ' +++++ gpake: only for one-dim histos'
      return
  901 write(nout,*) '+++ gpake: nonexisting id=',id
      write(6   ,*) '+++ gpake: nonexisting id=',id
      end


      subroutine grang1(id,ylr,yur)
*     *****************************
c provides y-scale for 1-dim plots
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0) return
      ist  = index(lact,2)
      ist2 = ist+7
      nch  = b(ist2 +1)
      yl   = b(ist+nbuf+1)
      yu   = b(ist+nbuf+1)
      do 10 ib=1,nch
      yl = min(yl,b(ist+nbuf+ib))
      yu = max(yu,b(ist+nbuf+ib))
   10 continue
      call goptou(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      if(iopsc1.eq.2) yl= b( ist +5)
      if(iopsc2.eq.2) yu= b( ist +6)
      ylr = yl
      yur = yu
      end

      subroutine ginbo2(id,nchx,xl,xu,nchy,yl,yu)
*     *******************************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.eq.0) goto 900
      ist  = index(lact,2)
      ist2 = ist+7
      nchx = b(ist2 +1)
      xl   = b(ist2 +2)
      xu   = b(ist2 +3)
      nchy = b(ist2 +5)
      yl   = b(ist2 +6)
      yu   = b(ist2 +7)
      return
  900 write(nout,*) ' +++ginbo2: nonexisting histo id= ',id 
      write(   6,*) ' +++ginbo2: nonexisting histo id= ',id
      end


      subroutine gmaxim(id,wmax)
*     **************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      if(id.ne.0) then
        lact=jadres(id)
        if(lact.eq.0) return
        ist= index(lact,2)
        b(ist+6) =wmax
        call gidopt(id,'YMAX')
      else
        do 20 k=1,idmx
        if(index(k,1).eq.0) goto 20
        ist=index(k,2)
        jd =index(k,1)
        b(ist+6) =wmax
        call gidopt(jd,'YMAX')
   20   continue
      endif
      end

      subroutine gminim(id,wmin)
*     **************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      if(id.ne.0) then
        lact=jadres(id)
        if(lact.eq.0) return
        ist =index(lact,2)
        b(ist+5) =wmin
        call gidopt(id,'YMIN')
      else
        do 20 k=1,idmx
        if(index(k,1).eq.0) goto 20
        ist=index(k,2)
        jd =index(k,1)
        b(ist+5) =wmin
        call gidopt(jd,'YMIN')
   20   continue
      endif
      end

      subroutine greset(id,chd1)
*     **************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      character*(*) chd1
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
c
      lact=jadres(id)
      if(lact.le.0) return
      ist  =index(lact,2)
      ist2 = ist+7
c 
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.1) then
c one-dim.
        ist3  = ist+11
        nchx  = b(ist2 +1)
        nch   = 2*nchx
        local = ist + nbuf
      elseif(ityphi.eq.2) then
c two-dim.
        ist3  = ist+15
        nchx  = b(ist2 +1)
        nchy  = b(ist2 +5)
        nch   = nchx*nchy
        local = ist +nbuf2
      else
         write(nout,*) '+++greset: wrong type id=',id
         write(6   ,*) '+++greset: wrong type id=',id
        return
      endif
c reset miscaelaneous entries and bins
      do 10 j=ist3+1,local +nch
  10  b(j)    = 0d0
c and no. of entries in index
      index(lact,3) = 0
      end

      SUBROUTINE GDELET(ID)
*     *********************
C Now it should work (stj Nov. 91) but watch out!
C should works for 2-dim histos, please check this!
*     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      logical gexist
c
      if(id.eq.0) goto 300
      if(.not.gexist(id)) goto 900
      lact = jadres(id)
      ist  = index(lact,2)
      ist2 = ist+7
      nch  = b(ist2 +1)
      iflag2   = nint(b(ist+4)-9d0-9d12)/10
      ityphi   = mod(iflag2,10)
      if(ityphi.eq.1) then
c one-dim.
        nchx  = b(ist2 +1)
        nch   = 2*nchx
c lenght of local histo to be removed
        local = nch+nbuf+1
      elseif(ityphi.eq.2) then
c two-dim.
        nchx  = b(ist2 +1)
        nchy  = b(ist2 +5)
        nch   = nchx*nchy
c lenght of local histo to be removed
        local = nch+nbuf2+1
      else
         write(nout,*) '+++gdelet: wrong type id=',id
         write(6   ,*) '+++gdelet: wrong type id=',id
        return
      endif
c starting position of next histo in storage b
      next = ist+1 +local
c move down all histos above this one 
      do 15 k =next,length
      b(k-local)=b(k)
   15 continue  
c define new end of storage
      length=length-local
c clean free space at the end of storage b
      do 20 k=length+1, length+local
   20 b(k)=0d0 
c shift adresses of all displaced histos 
      do 25 l=lact+1,idmx
      if(index(l,1).ne.0) index(l,2)=index(l,2)-local
   25 continue
c move entries in index down by one and remove id=lact entry
      do 30 l=lact+1,idmx
      index(l-1,1)=index(l,1)
      index(l-1,2)=index(l,2)
      index(l-1,3)=index(l,3)
      titlc(l-1)=titlc(l)
   30 continue
c last entry should be always empty
      index(idmx,1)=0
      index(idmx,2)=0
      index(idmx,3)=0 
      do 50 k=1,80
   50 titlc(idmx)(k:k)=' '
      return
C -----------------------------------
C Deleting all histos at once!!!
  300 length=0
      do 400 i=1,idmx
      do 340 k=1,3
  340 index(i,k)=0
      do 350 k=1,80
  350 titlc(i)(k:k)=' '
  400 continue
      return
  900 write(nout,*) ' +++gdelet: nonexisting histo id= ',id 
      end


      subroutine glimit(lenmx)
*     ************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /cglib/,/gind/
      call ginit
      if(lenmx.ge.lenmax) then
         lenmax=lenmx
      else
         write(nout,*) ' +++++ glimit: you cannot decrease storage'
         stop
      endif
      end

      subroutine copch(ch1,ch2)
*     *************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
* copies character*80 ch1 into ch2 up to a first $ sign
      character*80 ch1,ch2
      logical met
      met = .false.
      do 10 i=1,80
      if( ch1(i:i).eq.'$' .or. met )   then
        ch2(i:i)=' '
        met=.true.
      else
        ch2(i:i)=ch1(i:i)
      endif
  10  continue
      end

      function jadres(id)
*     *********************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      save /gind/
      jadres=0
      do 1 i=1,idmx
      if(index(i,1).eq.id) goto 2
    1 continue
      return
    2 jadres=i
      end


C--------------------------------------------------------------
C ----------- storing histograms in the disk file -------------
C--------------------------------------------------------------
      subroutine grfile(nhruni,dname,chd2)
c     ***********************************
      implicit double precision (a-h,o-z)
      character*(*) chd2, dname
      common / hruni / nhist
      save /hruni/
      nhist=nhruni
      end

      subroutine grout(idum1,idum2,chdum)
c     ***********************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      character*8 chdum
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      common / hruni / nhist
      character*80 titlc
      save /cglib/,/gind/, /hruni/
c
      call ginit
      nouth=nhist
      write(nouth,'(6i10)')   nvrs,nout,lenmax,length
      write(nouth,'(6i10)')   ((index(i,k),k=1,3),i=1,idmx)
      write(nouth,'(a80)')    titlc
      write(nouth,'(3d24.16)') (b(i),i=1,length)
      end


      SUBROUTINE GRIN(IDUM1,IDUM2,IDUM3)
C     **********************************
C New version which has possibility to ADD histograms
C identical ID is changed by adding 100000 !!!!
C compatible with GDELET which makes holes in INDEX
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
      common / hruni / nhist
      save /cglib/,/gind/, /hruni/
c copy from disk
      dimension lndex(idmx,3),titld(idmx)
      character*80 titld
      logical gexist
c
      call ginit 
      nouth=nhist
      read(nouth,'(6i10)')   nvrs3,nout3,lenma3,lengt3
      if(length+lengt3.ge.lenmax) goto 900
      if(nvrs.ne.nvrs3) write(nout,*)
     $ "  +++++ warning (grin): histos produced by older version",nvrs3
      read(nouth,'(6i10)')  ((lndex(i,k),k=1,3),i=1,idmx)
      read(nouth,'(a80)')    titld
c modify index and titlc
      lact = jadres(0)
      do 100 l=1,idmx
      if(lact.eq.0) goto 901
      idn= lndex(l,1)
      if(idn.eq.0) goto 100
c identical id is changed by adding 100000 !!!!
      if( gexist(idn) ) idn = idn +100000
      index(lact,1)=idn
      index(lact,2)=lndex(l,2)+length
      index(lact,3)=lndex(l,3)
      titlc(lact)  =titld(l)
      lact=lact+1
  100 continue  
c add content of new histos
      lenold=length
      length=length+lengt3
      read(nouth,'(3d24.16)') (b(i),i=lenold+1,length)
      return
  900 write(nout,*) ' ++++ stop in grin: to litle space '
      stop
  901 write(nout,*) ' ++++ stop in grin: to many histos '
      stop
      end

      subroutine grend(chdum)
c     ***********************
      implicit double precision (a-h,o-z)
      common / hruni / nhist
      save   /hruni/
      character*(*) chdum
      close(nhist)
c======================================================================
c======================end of gbook====================================
c======================================================================
      end

C======================================================================
C======================Mini-GPLOT======================================
C======================================================================
C... Plotting using LATeX
      SUBROUTINE GPLINT(IDUM)
C     ***********************
C ...dummy routine
      END
      SUBROUTINE GPLCAP(IFILE)
C     ***********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / LPLTIT / TITCH,KEYTIT
      CHARACTER*80 TITCH
C Note that backslash definition is varying from one 
C instalation/compiler to another, you have to figure out by yourself 
C how to fill backslash code into BS
      COMMON / BSLASH / BS
      CHARACTER*1 BS,BBS
      save /LPLDAT/, /LPLTIT/, /BSLASH/
C     DATA BBS / 1H\ /
      DATA BBS / "\\" /
      BS = BBS
cc      BS = "\\"
C---
      KEYTIT= 0
      ILINE = 1
      NOUH1=IABS(IFILE)
      NOUH2=NOUH1+1
      WRITE(NOUH1,'(A,A)') BS,'voffset =  1.0cm'
      WRITE(NOUH1,'(A,A)') BS,'hoffset = -1cm'
      WRITE(NOUH1,'(A,A)') BS,'documentstyle[12pt]{article}'
      WRITE(NOUH1,'(A,A)') BS,'textwidth  = 16cm'
      WRITE(NOUH1,'(A,A)') BS,'textheight = 24cm'
      WRITE(NOUH1,'(A,A)') BS,'begin{document}'
      WRITE(NOUH1,'(A)') '  '
      WRITE(NOUH1,'(A)') '  '
      END

      SUBROUTINE GPLEND
C     *****************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      save /LPLDAT/, /BSLASH/
      CHARACTER*1 BS
      WRITE(NOUH1,'(2A)') BS,'end{document}'
      CLOSE(NOUH1)
      END

      SUBROUTINE GPLOT(ID,CH1,CH2,KDUM)
C     *********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(200),YER(200)
      CHARACTER CH1,CH2,CHR
      CHARACTER*80 TITLE
      LOGICAL GEXIST
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/
      DATA CHR /' '/
C return if histo non-existing
      IF(.NOT.GEXIST(ID)) GOTO 900
C ...unpack histogram
      CALL GUNPAK(ID,YY ,'    ',IDUM)
      CALL GUNPAK(ID,YER,'ERRO',IDUM)
      CALL GINBO1(ID,TITLE,NCHX,DXL,DXU)
      XL = DXL
      XU = DXU
      CALL GRANG1(ID,YL,YU)
      KAX=1200
      KAY=1200
      IF(CH1.EQ.'S') THEN
C ...superimpose plot
        BACKSPACE(NOUH1)
      ELSE
C ...new frame only
        CHR=CH1
        CALL LFRAM1(ID,KAX,KAY)
      ENDIF
      WRITE(NOUH1,'(A)')    '%========== next plot (line) =========='
      WRITE(NOUH1,'(A,I6)') '%==== HISTOGRAM ID=',ID
      WRITE(NOUH1,'(A,A70 )') '% ',TITLE
C...cont. line for functions
      call goptou(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      ker = ioperb-1
      IF (iopsla.eq.2)  CHR='C'
C...suppress GPLOT assignments
      IF (CH2.EQ.'B')   CHR=' '
      IF (CH2.EQ.'*')   CHR='*'
      IF (CH2.EQ.'C')   CHR='C'
C...various types of lines
      IF     (CHR.EQ.' ') THEN
C...contour line used for histogram
          CALL PLHIST(KAX,KAY,NCHX,YL,YU,YY,KER,YER)
      ELSE IF(CHR.EQ.'*') THEN
C...marks in the midle of the bin
          CALL PLHIS2(KAX,KAY,NCHX,YL,YU,YY,KER,YER)
      ELSE IF(CHR.EQ.'C') THEN
C...slanted (dotted) line in plotting non-MC functions
          CALL PLCIRC(KAX,KAY,NCHX,YL,YU,YY)
      ENDIF
      WRITE(NOUH1,'(2A)') BS,'end{picture} % close entire picture '
      RETURN
  900 WRITE(*,*) ' ++++ GPLOT: NONEXISTIG HISTO ' ,ID
      END
      SUBROUTINE LFRAM1(ID,KAX,KAY)
C     *****************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      COMMON / LPLTIT / TITCH,KEYTIT
      CHARACTER*80 TITCH
      DIMENSION TIPSY(20),TIPSX(20)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      DOUBLE PRECISION DXL,DXU
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /LPLTIT/, /BSLASH/
      DATA ICONT/0/
      ICONT=ICONT+1
      CALL GINBO1(ID,TITLE,NCHX,DXL,DXU)
      XL = DXL
      XU = DXU
      CALL GRANG1(ID,YL,YU)
      IF(ICONT.GT.1) WRITE(NOUH1,'(2A)') BS,'newpage'
      WRITE(NOUH1,'(A)') '% =========== big frame, title etc. ======='
      WRITE(NOUH1,'(2A)') BS,'begin{center}'
      IF(KEYTIT.EQ.0) THEN
        WRITE(NOUH1,'(A)')     TITLE
      ELSE
        WRITE(NOUH1,'(A)')     TITCH
      ENDIF
      WRITE(NOUH1,'(2A)') BS,'end{center}'
      WRITE(NOUH1,'(4A)') BS,'setlength{',BS,'unitlength}{0.1mm}'
      WRITE(NOUH1,'(2A)') BS,'begin{picture}(1600,1500)'
      WRITE(NOUH1,'(4A)') BS,'put(0,0){',BS,'framebox(1600,1500){ }}'
      WRITE(NOUH1,'(A)') '% =========== small frame, labeled axis ==='
      WRITE(NOUH1,'(4A,I4,A,I4,A)')
     $    BS,'put(300,250){',BS,'begin{picture}( ',KAX,',',KAY,')'
      WRITE(NOUH1,'(4A,I4,A,I4,A)')
     $    BS,'put(0,0){',BS,'framebox( ',KAX,',',KAY,'){ }}'
      WRITE(NOUH1,'(A)') '% =========== x and y axis ================'
      CALL SAXIX(KAX,XL,XU,NTIPX,TIPSX)
      CALL SAXIY(KAY,YL,YU,NTIPY,TIPSY)
      WRITE(NOUH1,'(3A)') BS,'end{picture}}'
     $                ,'% end of plotting labeled axis'
      END
      SUBROUTINE SAXIX(KAY,YL,YU,NLT,TIPSY)
C     ***************************************
C plotting x-axis with long and short tips
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TIPSY(20)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/

      DY= ABS(YU-YL)
      LY = NINT( LOG10(DY) -0.49999 )
      JY = NINT(DY/10.**LY)
      DDYL = DY*10.**(-LY)
      IF( JY.EQ.1)             DDYL = 10.**LY*0.25
      IF( JY.GE.2.AND.JY.LE.3) DDYL = 10.**LY*0.5
      IF( JY.GE.4.AND.JY.LE.6) DDYL = 10.**LY*1.0
      IF( JY.GE.7)             DDYL = 10.**LY*2.0
      WRITE(NOUH1,'(A)') '% .......SAXIX........ '
      WRITE(NOUH1,'(A,I4)') '%  JY= ',JY
C-------
      NLT = INT(DY/DDYL)
      NLT = MAX0(MIN0(NLT,20),1)+1
      YY0L = NINT(YL/DDYL+0.5)*DDYL
      DDYS = DDYL/10.
      YY0S = NINT(YL/DDYS+0.49999)*DDYS
      P0L = KAY*(YY0L-YL)/(YU-YL)
      PDL = KAY*DDYL/(YU-YL)
      P0S = KAY*(YY0S-YL)/(YU-YL)
      PDS = KAY*DDYS/(YU-YL)
      NLT = INT(ABS(YU-YY0L)/DDYL+0.00001)+1
      NTS = INT(ABS(YU-YY0S)/DDYS+0.00001)+1
      DO 41 N=1,NLT
      TIPSY(N) =YY0L+ DDYL*(N-1)
  41  CONTINUE
      WRITE(NOUH1,1000)
     $ BS,'multiput('  ,P0L,  ',0)('  ,PDL,  ',0){'  ,NLT,  '}{',
     $ BS,'line(0,1){25}}',
     $ BS,'multiput('  ,P0S,  ',0)('  ,PDS,  ',0){'  ,NTS,  '}{',
     $ BS,'line(0,1){10}}'
      WRITE(NOUH1,1001)
     $ BS,'multiput('  ,P0L,  ','  ,KAY,  ')('  ,PDL,  ',0){'  ,NLT,
     $ '}{'  ,BS,  'line(0,-1){25}}',
     $ BS,'multiput('  ,P0S,  ','  ,KAY,  ')('  ,PDS,  ',0){'  ,NTS,
     $ '}{'  ,BS,  'line(0,-1){10}}'
 1000 FORMAT(2A,F8.2,A,F8.2,A,I4,3A)
 1001 FORMAT(2A,F8.2,A,I4,A,F8.2,A,I4,3A)
C ...labeling of axis
      SCMX = DMAX1(DABS(YL),DABS(YU))
      LEX  = NINT( LOG10(SCMX) -0.50001)
      DO 45 N=1,NLT
      K = KAY*(TIPSY(N)-YL)/(YU-YL)
      IF(LEX.LT.2.AND.LEX.GT.-1) THEN
C ...without exponent
      WRITE(NOUH1,'(2A,I4,5A,F8.3,A)')
     $ BS,'put(',K,',-25){',BS,'makebox(0,0)[t]{',BS,'large $ ',
     $ TIPSY(N), ' $}}'
      ELSE
C ...with exponent
      WRITE(NOUH1,'(2A,I4,5A,F8.3,2A,I4,A)')
     $ BS,'put('  ,K,  ',-25){',BS,'makebox(0,0)[t]{',BS,'large $ ',
     $ TIPSY(N)/(10.**LEX),BS,'cdot 10^{',LEX,'} $}}'
      ENDIF
  45  CONTINUE
      END
      SUBROUTINE SAXIY(KAY,YL,YU,NLT,TIPSY)
C     ***************************************
C plotting y-axis with long and short tips
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TIPSY(20)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/

      DY= ABS(YU-YL)
      LY = NINT( LOG10(DY) -0.49999 )
      JY = NINT(DY/10.**LY)
      DDYL = DY*10.**(-LY)
      IF( JY.EQ.1)             DDYL = 10.**LY*0.25
      IF( JY.GE.2.AND.JY.LE.3) DDYL = 10.**LY*0.5
      IF( JY.GE.4.AND.JY.LE.6) DDYL = 10.**LY*1.0
      IF( JY.GE.7)             DDYL = 10.**LY*2.0
      WRITE(NOUH1,'(A)') '% .......SAXIY........ '
      WRITE(NOUH1,'(A,I4)') '%  JY= ',JY
C-------
      NLT = INT(DY/DDYL)
      NLT = MAX0(MIN0(NLT,20),1)+1
      YY0L = NINT(YL/DDYL+0.49999)*DDYL
      DDYS = DDYL/10.
      YY0S = NINT(YL/DDYS+0.5)*DDYS
      P0L = KAY*(YY0L-YL)/(YU-YL)
      PDL = KAY*DDYL/(YU-YL)
      P0S = KAY*(YY0S-YL)/(YU-YL)
      PDS = KAY*DDYS/(YU-YL)
      NLT= INT(ABS(YU-YY0L)/DDYL+0.00001) +1
      NTS= INT(ABS(YU-YY0S)/DDYS+0.00001) +1
      DO 41 N=1,NLT
      TIPSY(N) =YY0L+ DDYL*(N-1)
  41  CONTINUE
C plotting tics on vertical axis
      WRITE(NOUH1,1000)
     $ BS,'multiput(0,'  ,P0L,  ')(0,'  ,PDL  ,'){'  ,NLT,  '}{',
     $ BS,'line(1,0){25}}',
     $ BS,'multiput(0,'  ,P0S,  ')(0,'  ,PDS,  '){'  ,NTS,  '}{',
     $ BS,'line(1,0){10}}'
      WRITE(NOUH1,1001)
     $ BS,'multiput('  ,KAY,  ','  ,P0L,  ')(0,'  ,PDL,  '){'  ,NLT,
     $ '}{',BS,'line(-1,0){25}}',
     $ BS,'multiput('  ,KAY,  ','  ,P0S,  ')(0,'  ,PDS,  '){'  ,NTS,
     $ '}{',BS,'line(-1,0){10}}'
 1000 FORMAT(2A,F8.2,A,F8.2,A,I4,3A)
 1001 FORMAT(2A,I4,A,F8.2,A,F8.2,A,I4,3A)
C ...Zero line if necessary
      Z0L = KAY*(-YL)/(YU-YL)
      IF(Z0L.GT.0D0.AND.Z0L.LT.FLOAT(KAY))
     $      WRITE(NOUH1,'(2A,F8.2,3A,I4,A)')
     $       BS,'put(0,'  ,Z0L,  '){',BS,'line(1,0){'  ,KAY,  '}}'
C ...labeling of axis
      SCMX = DMAX1(DABS(YL),DABS(YU))
      LEX  = NINT( LOG10(SCMX) -0.50001)
      DO 45 N=1,NLT
      K = KAY*(TIPSY(N)-YL)/(YU-YL)
      IF(LEX.LT.2.AND.LEX.GT.-1) THEN
C ...without exponent
      WRITE(NOUH1,'(2A,I4,5A,F8.3,A)')
     $  BS,'put(-25,'  ,K,  '){',BS,'makebox(0,0)[r]{',
     $  BS,'large $ '  ,TIPSY(N),  ' $}}'
      ELSE
C ...with exponent
      WRITE(NOUH1,'(2A,I4,5A,F8.3,2A,I4,A)')
     $ BS,'put(-25,'  ,K,  '){',BS,'makebox(0,0)[r]{',
     $ BS,'large $ '
     $ ,TIPSY(N)/(10.**LEX),  BS,'cdot 10^{'  ,LEX,  '} $}}'
      ENDIF
  45  CONTINUE
      END
      SUBROUTINE PLHIST(KAX,KAY,NCHX,YL,YU,YY,KER,YER)
C     ************************************************
C plotting contour line for histogram
C     ***********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(*),YER(*)
      CHARACTER*80 FMT1
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/
      WRITE(NOUH1,'(4A,I4,A,I4,A)')
     $  BS,'put(300,250){',BS,'begin{picture}( ',KAX,',',KAY,')'
      WRITE(NOUH1,'(A)') '% ========== plotting primitives =========='
C...various types of line
      IF(ILINE.EQ.1) THEN
         WRITE(NOUH1,'(2A)') BS,'thicklines '
      ELSE
         WRITE(NOUH1,'(2A)') BS,'thinlines '
      ENDIF
C...short macros for vertical/horizontal straight lines
      WRITE(NOUH1,'(8A)')
     $ BS,'newcommand{',BS,'x}[3]{',BS,'put(#1,#2){',
     $ BS,'line(1,0){#3}}}'
      WRITE(NOUH1,'(8A)')
     $ BS,'newcommand{',BS,'y}[3]{',BS,'put(#1,#2){',
     $ BS,'line(0,1){#3}}}'
      WRITE(NOUH1,'(8A)')
     $ BS,'newcommand{',BS,'z}[3]{',BS,'put(#1,#2){',
     $ BS,'line(0,-1){#3}}}'
C   error bars
      WRITE(NOUH1,'(8A)')
     $   BS,'newcommand{',BS,'e}[3]{',
     $   BS,'put(#1,#2){',BS,'line(0,1){#3}}}'
      IX0=0
      IY0=0
      DO 100 IB=1,NCHX
      IX1 = KAX*(IB-0.00001)/NCHX
      IY1 = KAY*(YY(IB)-YL)/(YU-YL)
      IDY = IY1-IY0
      IDX = IX1-IX0
      FMT1 = '(2(2A,I4,A,I4,A,I4,A))'
      IF( IDY.GE.0) THEN  
         IF(IY1.GE.0.AND.IY1.LE.KAY)
     $   WRITE(NOUH1,FMT1) BS,'y{',IX0,'}{',IY0,'}{',IDY,'}',
     $                     BS,'x{',IX0,'}{',IY1,'}{',IDX,'}'
      ELSE
         IF(IY1.GE.0.AND.IY1.LE.KAY)
     $   WRITE(NOUH1,FMT1) BS,'z{',IX0,'}{',IY0,'}{',-IDY,'}',
     $                     BS,'x{',IX0,'}{',IY1,'}{',IDX,'}'
      ENDIF
      IX0=IX1
      IY0=IY1
      IF(KER.EQ.1) THEN
        IX2 = KAX*(IB-0.5000)/NCHX
        IERR = KAY*((YY(IB)-YER(IB))-YL)/(YU-YL)
        IE = KAY*YER(IB)/(YU-YL)
        IF(IY1.GE.0.AND.IY1.LE.KAY) WRITE(NOUH1,8000) BS,IX2,IERR,IE*2
      ENDIF
 100  CONTINUE
8000  FORMAT(4(A1,2He{,I4,2H}{,I4,2H}{,I4,1H}:1X ))
      WRITE(NOUH1,'(3A)') BS,'end{picture}}',
     $       ' % end of plotting histogram'
C change line-style
      ILINE= ILINE+1
      IF(ILINE.GT.2) ILINE=1
      END
      SUBROUTINE PLHIS2(KAX,KAY,NCHX,YL,YU,YY,KER,YER)
C     ************************************************
C marks in the midle of the bin
C     **********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(*),YER(*)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/

      WRITE(NOUH1,'(4A,I4,A,I4,A)')
     $ BS,'put(300,250){',BS,'begin{picture}( ',KAX,',',KAY,')'
      WRITE(NOUH1,'(A)') '% ========== plotting primitives =========='
C...various types of mark
      IRAD1= 6
      IRAD2=10
      IF(ILINE.EQ.1) THEN
C   small filled circle
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle*{',IRAD1,'}}}'
      ELSEIF(ILINE.EQ.2) THEN
C   small open circle
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle{',IRAD1,'}}}'
      ELSEIF(ILINE.EQ.3) THEN
C   big filled circle
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle*{',IRAD2,'}}}'
      ELSEIF(ILINE.EQ.4) THEN
C   big open circle
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle{',IRAD2,'}}}'
C Other symbols
      ELSEIF(ILINE.EQ.5) THEN
       WRITE(NOUH1,'(10A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'makebox(0,0){$',BS,'diamond$}}}'
      ELSE
       WRITE(NOUH1,'(10A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'makebox(0,0){$',BS,'star$}}}'
      ENDIF
C   error bars
      WRITE(NOUH1,'(8A)')
     $   BS,'newcommand{',BS,'e}[3]{',
     $   BS,'put(#1,#2){',BS,'line(0,1){#3}}}'
      DO 100 IB=1,NCHX
      IX1 = KAX*(IB-0.5000)/NCHX
      IY1 = KAY*(YY(IB)-YL)/(YU-YL)
      IF(IY1.GE.0.AND.IY1.LE.KAY) WRITE(NOUH1,7000) BS,IX1,IY1
      IF(KER.EQ.1) THEN
        IERR = KAY*((YY(IB)-YER(IB))-YL)/(YU-YL)
        IE = KAY*YER(IB)/(YU-YL)
        IF(IY1.GE.0.AND.IY1.LE.KAY) WRITE(NOUH1,8000) BS,IX1,IERR,IE*2
      ENDIF
 100  CONTINUE
7000  FORMAT(4(A1,2Hr{,I4,2H}{,I4,1H}:1X ))
8000  FORMAT(4(A1,2He{,I4,2H}{,I4,2H}{,I4,1H}:1X ))
      WRITE(NOUH1,'(3A)') BS,'end{picture}}',
     $    ' % end of plotting histogram'
C change line-style
      ILINE= ILINE+1
      IF(ILINE.GT.6) ILINE=1
      END
      SUBROUTINE PLCIRC(KAX,KAY,NCHX,YL,YU,YY)
C     ****************************************
C plots equidistant points, four-point interpolation,
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(*),IX(3000),IY(3000)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      COMMON / BSLASH / BS
      CHARACTER*1 BS
      save /LPLDAT/, /BSLASH/
      SAVE DS

C ...various types of line
C ...distance between points is DS, radius of a point is IRAD
      IRAD2=6
      IRAD1=3
C .............
      WRITE(NOUH1,'(4A,I4,A,I4,A)')
     $  BS,'put(300,250){',BS,'begin{picture}( ',KAX,',',KAY,')'
      WRITE(NOUH1,'(A)') '% ========== plotting primitives =========='
      IF(ILINE.EQ.1) THEN
C   small filled circle
       DS = 10
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle*{',IRAD1,'}}}'
      ELSEIF(ILINE.EQ.2) THEN
C   small open circle
       DS = 10
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle{',IRAD1,'}}}'
      ELSEIF(ILINE.EQ.3) THEN
C   big filled circle
       DS = 20
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle*{',IRAD2,'}}}'
      ELSEIF(ILINE.EQ.4) THEN
C   big open circle
       DS = 20
       WRITE(NOUH1,'(8A,I3,A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'circle{',IRAD2,'}}}'
C Other symbols
      ELSEIF(ILINE.EQ.5) THEN
       DS = 20
       WRITE(NOUH1,'(10A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'makebox(0,0){$',BS,'diamond$}}}'
      ELSE
       DS = 20
       WRITE(NOUH1,'(10A)')
     $   BS,'newcommand{',BS,'r}[2]{',
     $   BS,'put(#1,#2){',BS,'makebox(0,0){$',BS,'star$}}}'
      ENDIF
      FACY = KAY/(YU-YL)
C plot first point
      AI  = 0.
      AJ  = (APROF((AI/KAX)*NCHX+0.5, NCHX, YY) -YL)*FACY
      IPNT =1
      IX(IPNT) = INT(AI)
      IY(IPNT) = INT(AJ)
      DX =  DS
      AI0 = AI
      AJ0 = AJ
C plot next points
      DO 100 IPOIN=2,3000
C iteration to get (approximately) equal distance among ploted points
      DO  50 ITER=1,3
      AI  = AI0+DX
      AJ  = (APROF((AI/KAX)*NCHX+0.5, NCHX, YY) -YL)*FACY
      DX  = DX *DS/SQRT(DX**2 + (AJ-AJ0)**2)
  50  CONTINUE
      IF(INT(AJ).GE.0.AND.INT(AJ).LE.KAY.AND.INT(AI).LE.KAX) THEN
         IPNT = IPNT+1
         IX(IPNT) = INT(AI)
         IY(IPNT) = INT(AJ)
      ENDIF
      AI0 = AI
      AJ0 = AJ
      IF(INT(AI).GT.KAX) GOTO 101
 100  CONTINUE
 101  CONTINUE
      WRITE(NOUH1,7000) (BS,IX(I),IY(I), I=1,IPNT)
7000  FORMAT(4(A1,2Hr{,I4,2H}{,I4,1H}:1X ))
      WRITE(NOUH1,'(2A)') BS,'end{picture}} % end of plotting line'
C change line-style
      ILINE= ILINE+1
      IF(ILINE.GT.2) ILINE=1
      END
      FUNCTION APROF(PX,NCH,YY)
C     *************************
C PX is a continuous extension of the index in array YY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(*)
      X=PX
      IF(X.LT.0.0.OR.X.GT.FLOAT(NCH+1)) THEN
        APROF= -1E-20
        RETURN
      ENDIF
      IP=INT(X)
      IF(IP.LT.2)     IP=2
      IF(IP.GT.NCH-2) IP=NCH-2
      P=X-IP
      APROF = -(1./6.)*P*(P-1)*(P-2)  *YY(IP-1)
     $        +(1./2.)*(P*P-1)*(P-2)  *YY(IP  )
     $        -(1./2.)*P*(P+1)*(P-2)  *YY(IP+1)
     $        +(1./6.)*P*(P*P-1)      *YY(IP+2)
      END
      SUBROUTINE GPLSET(CH,XX)
*     ************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / LPLDAT / NOUH1,NOUH2,ILINE
      save / LPLDAT /
      CHARACTER*4 CH
      KTY=NINT(XX)
      IF(CH.EQ.'DMOD') THEN
        ILINE=KTY
      ENDIF
      END
      SUBROUTINE GPLTIT(TITLE)
*     ************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      COMMON / LPLTIT / TITCH,KEYTIT
      CHARACTER*80 TITCH
      save / LPLTIT /
      KEYTIT=1
      DO 50 K=1,80
   50 TITCH(K:K)=' '
      CALL COPCH(TITLE,TITCH)
      END



      subroutine gmonit(mode,id,wt,wtmax,rn)
c     **************************************
c Utility program for monitoring m.c. rejection weights.
c ---------------------------------------------------------
C It is backward compatible with WMONIT except:
c  (1) for id=-1 one  should call as follows:
c      call(-1,id,0d0,1d0,1d0) or skip initialisation completely!
c  (2) maximum absolute weight is looked for,
c  (3) gprint(-id) prints weight distribution, net profit!
c  (4) no restriction id<100 any more!
c ---------------------------------------------------------
c wt is weight, wtmax is maximum weight and rn is random number.
c if(mode.eq.-1) then
c          initalization if entry id, 
c        - wtmax is maximum weight used for couting overweighted
c          other arguments are ignored
c elseif(mode.eq.0) then
c          summing up weights etc. for a given event for entry id
c        - wt is current weight.
c        - wtmax is maximum weight used for couting overweighted
c          events with wt>wtmax.
c        - rn is random number used in rejection, it is used to
c          count no. of accepted (rn<wt/wtmax) and rejected
c          (wt>wt/wtmax) events,
c          if ro rejection then put rn=0d0.
c elseif(mode.eq.1) then
c          in this mode wmonit repports on accumulated statistics
c          and the information is stored in common /cmonit/
c        - averwt= average weight wt counting all event
c        - errela= relative error of averwt
c        - nevtot= total nimber of accounted events
c        - nevacc= no. of accepted events (rn<wt\wtmax)
c        - nevneg= no. of events with negative weight (wt<0)
c        - nevzer= no. of events with zero weight (wt.eq.0d0)
c        - nevove= no. of overweghted events (wt>wtmax)
c          and if you do not want to use cmonit then the value
c          the value of averwt is assigned to wt,
c          the value of errela is assigned to wtmax and
c          the value of wtmax  is assigned to rn in this mode.
c elseif(mode.eq.2) then
c          all information defined for entry id defined above
c          for mode=2 is just printed of unit nout
c endif
c note that output repport (mode=1,2) is done dynamically just for a
c given entry id only and it may be repeated many times for one id and
c for various id's as well.
c     ************************
      implicit double precision (a-h,o-z)
      parameter( idmx=400,nbuf=24,nbuf2=24)
      common / cglib / b(20000)
      common /gind/ nvrs,nout,lenmax,length,index(idmx,3),titlc(idmx)
      character*80 titlc
c special gmonit common
      common / cmonit/ averwt,errela,nevtot,nevacc,nevneg,nevove,nevzer
      save / cglib /,/gind/, /cmonit/
c
      idg = -id
      if(id.le.0) then
           write(nout,*) ' =====wmonit: wrong id',id
           stop
      endif
      if(mode.eq.-1) then
c     *******************
           nbin = nint(dabs(rn))
           if(nbin.gt.100) nbin =100 
           if(nbin.eq.0)   nbin =1
           xl   =  wt
           xu   =  wtmax
           if(xu.le.xl) then
             xl = 0d0
             xu = 1d0
           endif
           lact=jadres(idg)
           if(lact.eq.0) then
              call gbook1(idg,' gmonit $',nbin,xl,xu)
           else
              call greset(idg,'  ')
           endif
      elseif(mode.eq.0) then
c     **********************
           lact=jadres(idg)
           if(lact.eq.0) then
              write(nout,*) ' ***** Gmonit initialized, id=',id
              call gbook1(idg,' gmonit $',1,0d0,1d0)
           endif
c     standard entries
           call gf1(idg,wt,1d0)
c     additional goodies
           ist  = index(lact,2)
           ist2 = ist+7
           ist3 = ist+11
c    maximum weight -- maximum by absolute value but keeping sign
           b(ist3+13)    = max( dabs(b(ist3+13)) ,dabs(wt))
           if(wt.ne.0d0) b(ist3+13)=b(ist3+13) *wt/dabs(wt)
c    nevzer,nevove,nevacc
           if(wt.eq.0d0)        b(ist3+10) =b(ist3+10) +1d0
           if(wt.gt.wtmax)      b(ist3+11) =b(ist3+11) +1d0
           if(rn*wtmax.le.wt)   b(ist3+12) =b(ist3+12) +1d0
      elseif(mode.ge.1.or.mode.le.3) then
c     ***********************************
           lact=jadres(idg)
           if(lact.eq.0) then
              write(nout,*) ' ==== warning from wmonit: '
              write(nout,*) ' lack of initialization, id=',id
              return
           endif
           ist    = index(lact,2)
           ist2   = ist+7
           ist3   = ist+11
           ntot = nint(b(ist3 +7))
           swt    =    b(ist3 +8)
           sswt   =    b(ist3 +9)
           if(ntot.le.0 .or. swt .eq. 0d0 )  then
              averwt=0d0
              errela=0d0
           else
              averwt=swt/float(ntot)
              errela=sqrt(abs(sswt/swt**2-1d0/float(ntot)))
           endif
c   output through common
           nevtot = ntot
           nevacc = b(ist3 +12)
           nevneg = b(ist3  +1)
           nevzer = b(ist3 +10)
           nevove = b(ist3 +11)
           wwmax  = b(ist3 +13)
c   output through parameters
           wt     = averwt
           wtmax  = errela
           rn     = wwmax
c  no printout for mode > 1
c  ************************
           if(mode.eq.1) return
           write(nout,1003) id, averwt, errela, wwmax
           write(nout,1004) nevtot,nevacc,nevneg,nevove,nevzer
           if(mode.eq.2) return
           call gprint(idg)
      else
c     ****
           write(nout,*) ' =====wmonit: wrong mode',mode
           stop
      endif
c     *****
 1003 format(
     $  ' =======================gmonit========================'
     $/,'   id           averwt         errela            wwmax'
     $/,    i5,           e17.7,         f15.9,           e17.7)
 1004 format(
     $  ' -----------------------------------------------------------'
     $/,'      nevtot      nevacc      nevneg      nevove      nevzer'
     $/,   5i12)
      end
