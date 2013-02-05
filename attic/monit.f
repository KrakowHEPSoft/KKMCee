*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------

      SUBROUTINE monin(mode)
*     **********************
*  monitornig weights and x-sections for the initial state bremss.
*  auxiliary weights on initial state only
*  second, first and zero order x-sections
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / cmonit/ averwt,errela,nevtot,nevacc,nevneg,nevove,nevzer
      COMMON / wgtall / wtmod,wtcru1,wtcru2,wtset(100)
      COMMON / inout  / ninp,nout
      COMMON / bxfmts / bxope,bxclo,bxtxt,bxl1i,bxl1f,bxl2f,bxl1g,bxl2g
      CHARACTER*80      bxope,bxclo,bxtxt,bxl1i,bxl1f,bxl2f,bxl1g,bxl2g
      SAVE / cmonit/,/ wgtall /,/ inout  /,/ bxfmts /
      SAVE idyfs

      IF(mode .EQ. -1) THEN
*     ===================
      idyfs = 0
      CALL gmonit(-1,idyfs+ 1,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+ 2,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+ 3,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+ 4,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+ 5,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+10,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+11,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+20,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+21,0d0,1d0,1d0)
      CALL gmonit(-1,idyfs+22,0d0,1d0,1d0)
      ELSEIF(mode .EQ. 0) THEN
*     ======================
      wtmax=2.5d0
      wtcrud = wtcru1*wtcru2
      CALL gmonit(0,idyfs+ 3,wtcrud*wtset( 3),wtmax,0d0)
      CALL gmonit(0,idyfs+ 2,wtcrud*wtset( 2),wtmax,0d0)
      CALL gmonit(0,idyfs+ 1,wtcrud*wtset( 1),wtmax,0d0)
* ...and dIFferences
      CALL gmonit(0,idyfs+ 5,wtcrud*(wtset(3)-wtset(2)),wtmax,0d0)
      CALL gmonit(0,idyfs+ 4,wtcrud*(wtset(2)-wtset(1)),wtmax,0d0)
* ...second order beta0,1,2
      CALL gmonit(0,idyfs+20,wtcrud*wtset(20),wtmax,0d0)
      CALL gmonit(0,idyfs+21,wtcrud*wtset(21),wtmax,0d0)
      CALL gmonit(0,idyfs+22,wtcrud*wtset(22),wtmax,0d0)
* ...first order beta0,1
      CALL gmonit(0,idyfs+10,wtcrud*wtset(10),wtmax,0d0)
      CALL gmonit(0,idyfs+11,wtcrud*wtset(11),wtmax,0d0)
      ELSEIF(mode .EQ. 1) THEN
*     ======================
      xkarl  = wtcru1
      erkarl = wtcru2
* ........................output winDOw a...............................
      WRITE(nout,bxope)
      WRITE(nout,bxtxt) '        model  output - winDOw a '
      WRITE(nout,bxtxt) '  monte carlo initial state only '
      WRITE(nout,bxtxt) '           x-sections in r-units '
      CALL gmonit(1,idyfs+ 3,dumm1,dumm2,dumm3)
*)))))CALL gmonit(2,idyfs+ 3,dumm1,dumm2,dumm3)
*     WRITE(2,*) ' yfsmod, xkarl,averwt',xkarl,averwt
*))))))))
      xs03   =  xkarl*averwt
      dxs03  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+ 2,dumm1,dumm2,dumm3)
      xs02   =  xkarl*averwt
      dxs02  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+ 1,dumm1,dumm2,dumm3)
      xs01   =  xkarl*averwt
      dxs01  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+ 5,dumm1,dumm2,dumm3)
      xs05   =  xkarl*averwt
      dxs05  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+ 4,dumm1,dumm2,dumm3)
      xs04   =  xkarl*averwt
      dxs04  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      WRITE(nout,bxl2f) xs03,dxs03,'x-section','o(alf2)',  'a01'
      WRITE(nout,bxl2f) xs02,dxs02,'x-section','o(alf1)',  'a02'
      WRITE(nout,bxl2f) xs01,dxs01,'x-section','o(alf0)',  'a03'
      IF(xs02.ne.0d0) WRITE(nout,bxl2f)
     $ xs05/xs02,dxs05/xs02,'(o(alf2)-o(alf1))','/o(alf1)','a04'
      IF(xs01.ne.0d0) WRITE(nout,bxl2f)
     $ xs04/xs01,dxs04/xs01,'(o(alf1)-o(alf0))','/o(alf0)','a05'
* ...beta contributions absolute
      CALL gmonit(1,idyfs+20,dumm1,dumm2,dumm3)
      xs20   =  xkarl*averwt
      dxs20  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+21,dumm1,dumm2,dumm3)
      xs21   =  xkarl*averwt
      dxs21  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+22,dumm1,dumm2,dumm3)
      xs22   =  xkarl*averwt
      dxs22  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+10,dumm1,dumm2,dumm3)
      xs10   =  xkarl*averwt
      dxs10  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+11,dumm1,dumm2,dumm3)
      xs11   =  xkarl*averwt
      dxs11  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      WRITE(nout,bxtxt) 'beta contributions in r-units'
      WRITE(nout,bxl2f) xs20,dxs20,'beta0          ','o(alf2)','a06'
      WRITE(nout,bxl2f) xs21,dxs21,'     beta1     ','       ','a07'
      WRITE(nout,bxl2f) xs22,dxs22,'          beta2','       ','a08'
      WRITE(nout,bxl2f) xs10,dxs10,'beta0          ','o(alf1)','a09'
      WRITE(nout,bxl2f) xs11,dxs11,'     beta1     ','       ','a10'
      WRITE(nout,bxl2f) xs01,dxs01,'beta0          ','o(alf0)','a11'
* ...beta contributions relative
      WRITE(nout,bxtxt) 'beta contributions - relative'
      IF(xs03.ne.0d0) WRITE(nout,bxl2f)
     $ xs20/xs03,dxs20/xs03,'bet0/(bet0+1+2)','o(alf2)','a12'
      IF(xs03.ne.0d0) WRITE(nout,bxl2f)
     $ xs21/xs03,dxs21/xs03,'bet1/(bet0+1+2)','       ','a13'
      IF(xs03.ne.0d0) WRITE(nout,bxl2f)
     $ xs22/xs03,dxs22/xs03,'bet2/(bet0+1+2)','       ','a14'
      IF(xs02.ne.0d0) WRITE(nout,bxl2f)
     $ xs10/xs02,dxs10/xs02,'bet0/(bet0+1)','o(alf1)'  ,'a15'
      IF(xs02.ne.0d0) WRITE(nout,bxl2f)
     $ xs11/xs02,dxs11/xs02,'bet1/(bet0+1)','       '  ,'a16'
      WRITE(nout,bxclo)
* --------------------------------------------------------------------
* ......................output winDOw b...............................
      WRITE(nout,bxope)
      WRITE(nout,bxtxt) '          model  output - winDOw b'
      WRITE(nout,bxtxt) '          initial state only cont.'
      WRITE(nout,bxtxt) 'analytical estimates of x-sections'
      prec = 1d-6
      ys00 =      bremkf( 300,prec)
      ys01 =      bremkf( 301,prec)
      ys02 =      bremkf( 302,prec)
      WRITE(nout,bxl1f) ys02,'x-section','o(alf2)','b1'
      WRITE(nout,bxl1f) ys01,'x-section','o(alf1)','b2'
      WRITE(nout,bxl1f) ys00,'x-section','o(alf0)','b3'
      WRITE(nout,bxl1f) ys02-ys01,'(o(alf2)-o(alf1))','/o(alf1)','b4'
      WRITE(nout,bxl1f) ys01-ys00,'(o(alf1)-o(alf0))','/o(alf0)','b5'
      ys20 = bremkf( 320,prec)
      ys21 = bremkf( 321,prec)
      ys22 = bremkf( 322,prec)
      ys10 = bremkf( 310,prec)
      ys11 = bremkf( 311,prec)
      WRITE(nout,bxl1f) ys20,'beta0          ','o(alf2)','b06'
      WRITE(nout,bxl1f) ys21,'     beta1     ','       ','b07'
      WRITE(nout,bxl1f) ys22,'          beta2','       ','b08'
      WRITE(nout,bxl1f) ys10,'beta0          ','o(alf1)','b09'
      WRITE(nout,bxl1f) ys11,'     beta1     ','       ','b10'
      WRITE(nout,bxclo)
* --------------------------------------------------------------------
* ......................output winDOw c...............................
      WRITE(nout,bxope)
      WRITE(nout,bxtxt) '             model  output - winDOw c'
      WRITE(nout,bxtxt) '                  initial state  only'
      WRITE(nout,bxtxt) 'comparison of mc and analytical calc.'
      WRITE(nout,bxtxt) '(montecarlo  - analytical)/analytical'
      rs01 = xs01/ys00 -1
      rs02 = xs02/ys01 -1
      rs03 = xs03/ys02 -1
      WRITE(nout,bxl2f) rs03,dxs03/ys02,'x-section','o(alf2)','c1'
      WRITE(nout,bxl2f) rs02,dxs02/ys01,'x-section','o(alf1)','c2'
      WRITE(nout,bxl2f) rs01,dxs01/ys00,'x-section','o(alf0)','c3'
      WRITE(nout,bxtxt) 'beta contributions'
      rs20 = xs20/ys20 -1
      rs21 = xs21/ys21 -1
      rs22 = xs22/ys22 -1
      rs10 = xs10/ys10 -1
      rs11 = xs11/ys11 -1
      drs20 = dxs20/ys20
      drs21 = dxs21/ys21
      drs22 = dxs22/ys22
      drs10 = dxs10/ys10
      drs11 = dxs11/ys11
      WRITE(nout,bxl2f) rs20,drs20,'beta0          ','o(alf2)','c04'
      WRITE(nout,bxl2f) rs21,drs21,'     beta1     ','       ','c05'
      WRITE(nout,bxl2f) rs22,drs22,'          beta2','       ','c06'
      WRITE(nout,bxl2f) rs10,drs10,'beta0          ','o(alf1)','c07'
      WRITE(nout,bxl2f) rs11,drs11,'     beta1     ','       ','c08'
      WRITE(nout,bxclo)
* ...the END of the initial-state report..............................
      ELSE
      WRITE(nout,*) ' +++++ wrong mode in monin'
      ENDIF
*     =====
      END


      SUBROUTINE monif
*     ****************
*  monitornig weights and x-sections for the final state bremss.
*  second, first and zero order x-sections
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / cmonit/ averwt,errela,nevtot,nevacc,nevneg,nevove,nevzer
      COMMON / wgtall / wtmod,wtcru1,wtcru2,wtset(100)
      COMMON / inout  / ninp,nout
      COMMON / bxfmts / bxope,bxclo,bxtxt,bxl1i,bxl1f,bxl2f,bxl1g,bxl2g
      CHARACTER*80      bxope,bxclo,bxtxt,bxl1i,bxl1f,bxl2f,bxl1g,bxl2g
      SAVE / cmonit/,/ wgtall /,/ inout  /,/ bxfmts /
      SAVE idyfs
* --------------------------------------------------------------------
* x-sctions in born units (for the moment)
      xkarl = 1
      erkarl= 0
      idyfs = 0
      CALL gmonit(1,idyfs+73,dumm1,dumm2,dumm3)
      xs03   =  xkarl*averwt
      dxs03  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+72,dumm1,dumm2,dumm3)
      xs02   =  xkarl*averwt
      dxs02  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+71,dumm1,dumm2,dumm3)
      xs01   =  xkarl*averwt
      dxs01  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+75,dumm1,dumm2,dumm3)
      xs05   =  xkarl*averwt
      dxs05  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+74,dumm1,dumm2,dumm3)
      xs04   =  xkarl*averwt
      dxs04  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
* ...beta contributions
      CALL gmonit(1,idyfs+90,dumm1,dumm2,dumm3)
      xs20   =  xkarl*averwt
      dxs20  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+91,dumm1,dumm2,dumm3)
      xs21   =  xkarl*averwt
      dxs21  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+92,dumm1,dumm2,dumm3)
      xs22   =  xkarl*averwt
      dxs22  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+80,dumm1,dumm2,dumm3)
      xs10   =  xkarl*averwt
      dxs10  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
      CALL gmonit(1,idyfs+81,dumm1,dumm2,dumm3)
      xs11   =  xkarl*averwt
      dxs11  =  xkarl*averwt*sqrt(erkarl**2+errela**2)
* ......................output winDOw b...............................
      WRITE(nout,bxope)
      WRITE(nout,bxtxt) '          model  output - winDOw b'
      WRITE(nout,bxtxt) '          final   state only      '
      WRITE(nout,bxtxt) 'analytical estimates of x-sections'
      WRITE(nout,bxtxt) 'units of born x-section           '
      prec = 1d-6
      ys00 =      brefkf( 300,prec)
      ys01 =      brefkf( 301,prec)
      ys02 =      brefkf( 302,prec)
      WRITE(nout,bxl1f) ys02,'x-section','o(alf2)','b1'
      WRITE(nout,bxl1f) ys01,'x-section','o(alf1)','b2'
      WRITE(nout,bxl1f) ys00,'x-section','o(alf0)','b3'
      WRITE(nout,bxl1f) ys02-ys01,'(o(alf2)-o(alf1))','/o(alf1)','b4'
      WRITE(nout,bxl1f) ys01-ys00,'(o(alf1)-o(alf0))','/o(alf0)','b5'
      ys20 = brefkf( 320,prec)
      ys21 = brefkf( 321,prec)
      ys22 = brefkf( 322,prec)
      ys10 = brefkf( 310,prec)
      ys11 = brefkf( 311,prec)
      WRITE(nout,bxl1f) ys20,'beta0          ','o(alf2)','b06'
      WRITE(nout,bxl1f) ys21,'     beta1     ','       ','b07'
      WRITE(nout,bxl1f) ys22,'          beta2','       ','b08'
      WRITE(nout,bxl1f) ys10,'beta0          ','o(alf1)','b09'
      WRITE(nout,bxl1f) ys11,'     beta1     ','       ','b10'
      WRITE(nout,bxclo)
* --------------------------------------------------------------------
* ......................output winDOw c...............................
      WRITE(nout,bxope)
      WRITE(nout,bxtxt) '             model  output - winDOw c'
      WRITE(nout,bxtxt) ' final   state  only'
      WRITE(nout,bxtxt) 'comparison of mc and analytical calc.'
      WRITE(nout,bxtxt) '(montecarlo  - analytical)/analytical'
      rs01 = xs01/ys00 -1
      rs02 = xs02/ys01 -1
      rs03 = xs03/ys02 -1
      WRITE(nout,bxl2f) rs03,dxs03/ys02,'x-section','o(alf2)','c1'
      WRITE(nout,bxl2f) rs02,dxs02/ys01,'x-section','o(alf1)','c2'
      WRITE(nout,bxl2f) rs01,dxs01/ys00,'x-section','o(alf0)','c3'
      WRITE(nout,bxtxt) 'beta contributions'
      rs20 = xs20/ys20 -1
      rs21 = xs21/ys21 -1
      rs22 = xs22/ys22 -1
      rs10 = xs10/ys10 -1
      rs11 = xs11/ys11 -1
      drs20 = dxs20/ys20
      drs21 = dxs21/ys21
      drs22 = dxs22/ys22
      drs10 = dxs10/ys10
      drs11 = dxs11/ys11
      WRITE(nout,bxl2f) rs20,drs20,'beta0          ','o(alf2)','c04'
      WRITE(nout,bxl2f) rs21,drs21,'     beta1     ','       ','c05'
      WRITE(nout,bxl2f) rs22,drs22,'          beta2','       ','c06'
      WRITE(nout,bxl2f) rs10,drs10,'beta0          ','o(alf1)','c07'
      WRITE(nout,bxl2f) rs11,drs11,'     beta1     ','       ','c08'
      WRITE(nout,bxclo)
* ...the END of the final-state report..............................
      END


*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------

      FUNCTION bremkf(key,erel)
*     *************************
* non-montecarlo integration of the v-distribution
* gauss method, change of variables with help of chbin1
* see vvdisb
* key= 1,2,3,...for various distributions
* key= 3 for mc generation, other for tests
* for keyfix=1, exeptionally, it provides integrand at vv=vvmax
* with born omitted
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / keyyfs / keyzet,keybrm,keyfix,keyred,keywgt
* COMMON keydst communicates only with vvdisb - integrand FUNCTION
      COMMON / keydst / keydis
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON / vvrec  / vvmin,vvmax,vv,beti
      SAVE
      EXTERNAL vvdisb
*
      keydis=key
      IF(keyfix .EQ. 0) THEN
         xborn  =borny(4d0*ene**2)
         prec=  xborn*erel
         xa= 0d0
         xb= 1d0
         CALL gausjd(vvdisb,xa,xb,prec,result)
         bremkf=result
      ELSE
         svar  = 4d0*ene**2
         bremkf= vvrho(keydis,svar,amel,vvmax,vvmin)
     $          /vvrho(     9,svar,amel,vvmax,vvmin)
      ENDIF
      END

      FUNCTION vvdisb(r)
*     ******************
* integrand for bremkf
* mapping xx => vv change  to improve on efficiency
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( fleps =1d-35)
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON / vvrec  / vvmin,vvmax,vv,beti
      COMMON / keydst / keydis
      SAVE
*
      keyd=keydis
      x = max(r,fleps**beti)
      alf=  beti
      bet=  1d0
* ...special cases
* ...monte carlo crude distr
      IF    (keyd .EQ. 1)  THEN
        bet=  -0.5d0
* ...yfs exponentiation beta0,1,2 contribs
      ELSEIF(keyd .EQ. 310)  THEN
        alf=  beti
      ELSEIF(keyd .EQ. 311)  THEN
        alf=  beti +1
      ELSEIF(keyd .EQ. 320)  THEN
        alf=  beti
      ELSEIF(keyd .EQ. 321)  THEN
        alf=  beti +1
      ELSEIF(keyd .EQ. 322)  THEN
        alf=  beti +2
* ...reference distr including dilatation factor damel
      ELSEIF(keyd .EQ. 12) THEN
        bet=  -0.5
      ENDIF
      CALL chbin1(x,alf,bet,vvmax,vv,rjac)
* born xsection
* note 1/(1-vv) factor because borny is in r-units
      svar   = 4d0*ene**2
      svar1  = svar*(1d0-vv)
      xborn  = borny(svar1)/(1d0-vv)
      vvdisb = vvrho(keyd,svar,amel,vv,vvmin) *rjac*xborn
      END

      FUNCTION vvrho(keydis,svar,amel,vv,vvmin)
*     *****************************************
*-------------------------------------------------------------
* convention for keydis
* pedagogical exercises
*     keydis   =  1      crude distribution for initial state mc
*     keydis   =  9      reference distr.  of yfs2 cpc paper
*     keydis   =  50-52  obsolete test distr. for yfs2 cpc paper
*     keydis   =  101    soft part yfs       first  order
*     keydis   =  102    soft part yfs       second order
*     keydis   =  105    hard non-exp.       first  order
*     keydis   =  106    hard non-exp.       second order
* total results
*     keydis   =  0 + r*100                  zero   order
*     keydis   =  1 + r*100                  first  order
*     keydis   =  2 + r*100                  second order
*     keydis   = 15     reference distr. of yfs paper
* beta contributions
*     keydis   = 10 + r*100      beta0       zero   order
*     keydis   = 11 + r*100      beta0       first  order
*     keydis   = 12 + r*100      beta1
*     keydis   = 20 + r*100      beta0       second order
*     keydis   = 21 + r*100      beta1
*     keydis   = 22 + r*100      beta2
*     r = 200 kuraev-fadin
*     r = 300 yfs (pragmatic)
*     r = 400 yfs single electron ll str. funct.
*-------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(pi= 3.1415926535897932d0, alfinv=137.03604d0)
      PARAMETER(alf1   = 1d0/pi/alfinv)
      PARAMETER(ceuler =0.57721566d0)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
*
      keyd = keydis
      bilg   = dlog(svar/amel**2)
      beti   = 2d0*alf1*(bilg-1d0)
*===================================================================
* ---------------------- keyd = 1 ----------------------------------
* ---- crude distribution in yfs2 initial state monte carlo --------
* ------------------------------------------------------------------
* dilat is related to dilatation jacobian in yfsgen
* damel is responsible for modIFication of photon ang. distribution
* see also weight wt=wt1 in   angbre
      IF(keyd .GE. 1 .AND. keyd .LT. 100) THEN
         dilat=1d0
         IF(vv .GT. vvmin) dilat=(1d0+1d0/sqrt(1d0-vv))/2d0
         beti2  = 2d0*alf1*bilg
         damel=1d0
         IF(vv .GT. vvmin) damel=beti2/beti*(vv/vvmin)**(beti2-beti)
*---------
         IF    (keyd .EQ. 1)  THEN
            distr= beti*vv**(beti-1d0)*dilat*damel
* ...reference distribution used in yfs2 paper --------------------
         ELSEIF(keyd .EQ.  9)  THEN
            distr= beti*vv**(beti-1d0)*(1+(1-vv)**2)/2
* basic reference distribution  xrefer=sigma-ref
         ELSEIF(keyd .EQ. 50) THEN
            distr= beti*vv**(beti-1d0)
* xrefer times damel
         ELSEIF(keyd .EQ. 51) THEN
            distr= beti*vv**(beti-1d0)*damel
* xrefer times dilatation factor dilat
         ELSEIF(keyd .EQ. 52) THEN
            distr= beti*vv**(beti-1d0)*dilat
         ENDIF
* ------------------------------------------------------------
* ------------- soft part only  -- yfs exponentiation
* ------------------------------------------------------------
      ELSEIF(keyd .GE.  101 .AND. keyd .LE. 102) THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac =exp(-ceuler*beti)/dpgamm(1d0+beti)
* ....first order
         IF(keyd   .EQ.  101) THEN
            dels = alf1*(bilg-1d0)
* ....second order
         ELSEIF(keyd   .EQ.  102) THEN
            dels = alf1*(bilg-1d0) +0.5d0*(alf1*bilg)**2
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*( 1d0+dels)
* ---------------------------------------------------------------
* ------------- hard part only - no exponentiation
* ---------------------------------------------------------------
      ELSEIF(keyd .GE.  105 .AND. keyd .LE. 106 ) THEN
         z=1d0-vv
* ....first  order
         IF(keyd   .EQ.  105) THEN
            distr =  beti*(1d0 +z**2)/2d0/vv
* ....second order
         ELSEIF(keyd   .EQ.  106) THEN
            dist1 =  beti*(1d0 +z**2)/2d0/vv
            delh2 =  (alf1*bilg)**2*(
     $             -(1d0+z*z)/vv*dlog(z)
     $              +(1d0+z)*(0.5d0*dlog(z)-2d0*dlog(vv))
     $              -2.5d0-0.5d0*z  )
            dist2 =  (alf1*bilg)**2*( 3 +4*dlog(vv))/vv +delh2
            distr =  dist1+dist2
         ENDIF
* ---------------------------------------------------------------
* ---------------------------------------------------------------
* ----- exponentiation following kuraiev fadin prescription
* ---------------------------------------------------------------
      ELSEIF(keyd .GE. 201 .AND. keyd .LE. 202) THEN
* first order  from gerrit burgers (polarization book)
         dz2 = pi**2/6d0
         z=1d0-vv
         IF(keyd   .EQ. 201) THEN
            delvs= alf1*(1.5d0*bilg +2d0*dz2-2d0)
            delh=  -alf1*(1d0+z)*(bilg-1d0)
         ELSEIF(keyd   .EQ. 202) THEN
            delvs= alf1*(1.5d0*bilg +2d0*dz2-2d0)
     $            +alf1**2*(9d0/8d0-2d0*dz2)*bilg**2
            delh=  -alf1*(1d0+z)*(bilg-1d0)
     $      +alf1**2*( -(1d0+z*z)/vv     *dlog(z)
     $              +(1d0+z)*(0.5d0*dlog(z)-2d0*dlog(vv))
     $              -2.5d0-0.5d0*z)*bilg**2
         ENDIF
         distr=  beti*vv**(beti-1d0)*( 1d0+delvs) +delh
* ---------------------------------------------------------------
* ------------- yfs ad-hoc exponentiation -----------------------
* ---------------------------------------------------------------
      ELSEIF(keyd .GE. 300 .AND. keyd .LE. 302) THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac =exp(-ceuler*beti)/dpgamm(1d0+beti)
* ....zero   order
         IF(keyd   .EQ. 300)  THEN
            dels = 0d0
            delh = -beti/4 *log(1-vv)
* ....first  order
         ELSEIF(keyd   .EQ. 301)  THEN
            dels = beti/2
            delh = vv*(-1 +vv/2)
     $      -beti/2*vv**2 - beti/4*(1-vv)**2*log(1-vv)
* ....second order
         ELSEIF(keyd   .EQ. 302)  THEN
            dels = alf1*(bilg-1d0) +0.5d0*(alf1*bilg)**2
            delh = vv*(-1d0+vv/2d0)
     $      +alf1*bilg*(-0.25d0*(4d0-6d0*vv+3d0*vv**2)*dlog(1d0-vv) -vv)
         ENDIF
         distr= gamfac*exp(delb)*beti*vv**(beti-1d0)*(1 +dels + delh )
*----------------------------------------------------------------
*-------------- yfs ad-hoc exponentiation -----------------------
*-------------contributions  from various beta's ----------------
*----------------------------------------------------------------
      ELSEIF(keyd .GE. 310 .AND. keyd .LE. 322)  THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/dpgamm(1d0+beti)
         soft   = 0d0
         delh   = 0d0
* ...beta0 first  order
         IF(keyd .EQ. 310) THEN
            soft = 1 + beti/2
            delh = -beti/4 *log(1-vv)
* ...beta1 first  order
         ELSEIF(keyd .EQ.  311)  THEN
            delh =
     $      vv*(-1d0+vv/2/(1+beti))*(1-0.5*beti*log(1-vv))
* ...beta0 second order
         ELSEIF(keyd .EQ. 320) THEN
            soft = 1 + beti/2  +beti**2/8
            delh = -beti/4 *log(1-vv)
* ...beta1 second order
         ELSEIF(keyd .EQ.  321)  THEN
            delh = vv*(-1+vv/2)
     $      -beti*vv/2 -beti*vv**2/4 +beti/8*(-2+6*vv-3*vv**2)*log(1-vv)
* ...beta2 second order
         ELSEIF(keyd .EQ.  322)  THEN
            delh =    beti*  vv**2/4d0
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
* -------------------------------------------------------------------
* -------------------yfs formula-------------------------------------
* -------------single fermion ll fragmentation ----------------------
* -------------------------------------------------------------------
      ELSEIF(keyd .GE. 400 .AND. keyd .LE. 402)  THEN
*&&&&&&  delb   = beti/4
         delb   = beti/4  +alf1*(-0.5d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/dpgamm(1d0+beti)
         soft   = 0d0
         delh   = 0d0
* ...zero   order
         IF(keyd .EQ. 400) THEN
            soft = 1
* ...first  order
         ELSEIF(keyd .EQ. 401) THEN
            soft = 1 + beti/2
            delh =
     $      vv*(-1d0+vv/2/(1+beti))
* ...second order
         ELSEIF(keyd .EQ. 402) THEN
            soft = 1 + beti/2  +beti**2/8
*           delh =
*    $      vv*(1+0.5*beti)*(-1d0+vv/2/(1+beti))
*    $      -0.50*beti*log(1-vv)*0.5*(1+(1-vv)**2)
*    $      + beti/4d0*vv*(vv +(1-vv/2)*dlog(1d0-vv))
            delh =
     $      vv*(-1d0+vv/2)
     $      + beti*(-vv/2 -(1+3*(1-vv)**2)/8*dlog(1d0-vv))
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
* -------------------------------------------------------------------
* -------------single fermion ll fragmentation ----------------------
* -------------contributions  from various beta's -------------------
* -------------------------------------------------------------------
      ELSEIF(keyd .GE. 400 .AND. keyd .LE. 422)  THEN
*&&&&&   delb   = beti/4
         delb   = beti/4  +alf1*(-0.5d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/dpgamm(1d0+beti)
         soft   = 0d0
         delh   = 0d0
* ...beta0 zero   order
         IF(keyd .EQ. 400) THEN
            soft = 1
* ...beta0 first  order
         ELSEIF(keyd .EQ. 410) THEN
            soft = 1 + beti/2
* ...beta1 first  order
         ELSEIF(keyd .EQ.  411)  THEN
            delh =
     $      vv*(-1d0+vv/2/(1+beti))
* ...beta0 second order
         ELSEIF(keyd .EQ. 420) THEN
            soft = 1 + beti/2  +beti**2/8
* ...beta1 second order
         ELSEIF(keyd .EQ.  421)  THEN
            delh =
     $      vv*(1+0.5*beti)*(-1d0+vv/2/(1+beti))
     $      -0.50*beti*log(1-vv)*0.5*(1+(1-vv)**2)
* ...beta2 second order
         ELSEIF(keyd .EQ.  422)  THEN
            delh = beti/4d0*vv*(vv +(1-vv/2)*dlog(1d0-vv))
         ELSE
            GOTO 900
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
* ----------------------------------------------------------------
* -------------single fermion ll fragmentation -------------------
* -------------i n f i n i t e   o r d e r -----------------------
* ----------------------------------------------------------------
      ELSEIF(keyd .EQ. 502)  THEN
         delb   = beti*0.75d0
         gamfac = exp(-ceuler*beti)/dpgamm(1d0+beti)
         soft   = 1
         delh =
     $   vv*(-1d0+vv/2)
     $   + beti*(-vv**2/4 -(1+3*(1-vv)**2)/8*dlog(1d0-vv))
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
      ELSE
         GOTO 900
      ENDIF
      vvrho = distr
      RETURN
 900  WRITE(6,*) ' ===--->  wrong keydis in vvrho',keyd
      STOP
      END

      FUNCTION brefkf(key,erel)
*     *************************
* non-montecarlo integration of the u-distribution final st.bremss.
* gauss method, change of variables with help of chbin1
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
* COMMON keydst communicates only uudisb - integrand FUNCTION
      COMMON / keydst / keydis
      SAVE   / keydst /
      EXTERNAL uudisb
*
      keydis=key
      prec= dabs( erel)
      xa= 0d0
      xb= 1d0
      CALL gausjd(uudisb,xa,xb,prec,result)
      brefkf=result
      END
      FUNCTION uudisb(r)
*     ******************
* integrand for brefkf
* mapping xx => uu change  to improve on efficiency
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( fleps =1d-35)
      PARAMETER( pi=3.1415926535897932d0,alfinv= 137.03604d0)
      PARAMETER( alf1=1d0/alfinv/pi)
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON / vvrec  / vvmin,vvmax,vv,beti
      COMMON / keydst / keydis
      SAVE / weking /,/ vvrec  /,/ keydst /
*
      keyd=keydis
      betf= 2/alfinv/pi*(dlog(4*ene**2/amfin**2)-1)
      x = max(r,fleps**betf)
      alf=  betf
      bet=  1d0
* ...yfs exponentiation beta0,1,2 contribs
      IF(keyd .EQ. 310) THEN
        alf=  betf
      ELSEIF(keyd .EQ. 311)  THEN
        alf=  betf +1
      ELSEIF(keyd .EQ. 320) THEN
        alf=  betf
      ELSEIF(keyd .EQ. 321)  THEN
        alf=  betf +1
      ELSEIF(keyd .EQ. 322)  THEN
        alf=  betf +2
      ENDIF
      CALL chbin1(x,alf,bet,vvmax,uu,rjac)
      svar   = 4d0*ene**2
      uudisb=uurho(keyd,svar,amfin,uu)*rjac
      END

      FUNCTION uurho(keydis,svar,amfin,uu)
*     ******************************************
*--------------------------------------------------------------
* the parametrization of the final state bremss. as in yfs3
* various types of the rho(u) distribution
*--------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(pi= 3.1415926535897932d0, alfinv=137.03604d0)
      PARAMETER(ceuler =0.57721566d0)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
*
      keyd   = keydis
      alf1   = 1d0/pi/alfinv
      sprim  = svar*(1-uu)
      bilg   = dlog(sprim/amfin**2)
      betf   = 2d0*alf1*(bilg-1d0)
* -------------yfs formula ------------------------------------------
      IF(keyd .GE. 300 .AND. keyd .LE. 302) THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac =exp(-ceuler*betf)/dpgamm(1d0+betf)
* ....zero   order
         IF(keyd   .EQ. 300)  THEN
            dels = 0d0
            delh = -betf/4 *log(1-uu)
* ....first  order
         ELSEIF(keyd   .EQ. 301)  THEN
            dels = betf/2
            delh = uu*(-1 +uu/2)
     $       +betf*(-uu**2/2 - 0.25d0*(1-uu)**2*log(1-uu) )
* final state specIFic contr.
     $       +betf*uu/2*(1-uu/2)*log(1-uu)
* ....second order
         ELSEIF(keyd   .EQ. 302)  THEN
            dels  = betf/2d0 +betf**2/8d0
            delh  = uu*(-1d0+uu/2d0)
     $        +betf*( -uu/2  +(2*uu-uu**2)/8*log(1-uu) )
         ENDIF
         distr= gamfac*exp(delb)*betf*uu**(betf-1d0)*(1 +dels +delh )
* -------------------------------------------------------------------
* -------------contributions  from various beta's -------------------
      ELSEIF(keyd .GE. 310 .AND. keyd .LE. 322)  THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac = exp(-ceuler*betf)/dpgamm(1d0+betf)
         soft   = 0d0
         delh   = 0d0
* ...beta0 first  order
         IF(keyd .EQ. 310) THEN
            soft = 1 + betf/2
            delh = -betf/4 *log(1-uu)
* ...beta1 first  order
         ELSEIF(keyd .EQ.  311)  THEN
            delh =
     $      uu*(-1d0+uu/2/(1+betf))*(1-0.5*betf*log(1-uu))
     $       -betf*uu/2*(-1+uu/2)*log(1-uu)
* ...beta0 second order
         ELSEIF(keyd .EQ. 320) THEN
            soft = 1 + betf/2  +betf**2/8
            delh = -betf/4 *log(1-uu)
* ...beta1 second order
         ELSEIF(keyd .EQ.  321)  THEN
            delh =
     $      uu*(1+0.5*betf)*(-1d0+uu/2/(1+betf))*(1-0.5*betf*log(1-uu))
     $       +0.25*betf*log(1-uu)*0.5*(1+(1-uu)**2)
* final state specIFic contrib.
     $       +betf*uu/2*(1-uu/2)*log(1-uu)
*           delh = uu*(-1d0+uu/2d0)
*    $       +betf*(-uu/2 -uu**2/4 +(2+2*uu-uu**2)/8d0*log(1-uu))
* ...beta2 second order
         ELSEIF(keyd .EQ.  322)  THEN
            delh =    betf*uu/2*( uu/2  -(1-uu/2)*log(1-uu) )
         ENDIF
         distr=  gamfac*exp(delb)*betf*uu**(betf-1d0)*(soft+delh)
      ELSE
         WRITE(nout,*) ' ===--->  wrong keydis in uurho'
         STOP
      ENDIF
      uurho = distr
      END
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------

