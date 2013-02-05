***********************************************************************
* library of utilities for yfs3 and bhlumi PROGRAMs
***********************************************************************

      SUBROUTINE chbin1(r,alf,bet,xmax,x,djac)
*     ****************************************
* This mapps variable r into x.
* to be employed in the integration (either ordinary or Monte Carlo)
* of distributions resambling
* the binomial distribution x**(alf-1)*(1-x)**(bet-1)
* with alf > 0 and  bet arbitrary.
* variable r is in (0,1) range and x is within (0,xmax) range.
* djac is jacobian factor d(x)/d(r).
* mapping is such that 1/djac is very CLOSE to
* binomial distribution x**(alf-1)*(1-x)**(bet-1).
* warning: mapping may fail very CLOSE to r=0. practiCALLy, one is
* recommended to obey: fleps**alf < r, where fleps = 1.d-30.
* problems may also arise for very small xmax ( below 1.d-12 ).
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      SAVE
*
      IF( alf .LE. 0d0 ) GOTO 900
      x0=(alf-1d0)/(alf+bet-2d0)
      IF(x0 .GT. xmax) x0=xmax
      x0= max(x0, 0d0)
      q1= 1d0/alf*x0**alf  *(1d0-x0)**(bet-1d0)
      q2= x0**(alf-1d0) /bet*((1d0-x0)**bet-(1d0-xmax)**bet)
      p1= q1/(q1+q2)
      IF( r .LE. p1 ) THEN
         x=x0*(r/p1)**(1d0/alf)
         dist= x**(alf-1d0)*(1d0-x0)**(bet-1d0)
      ELSE
         r1= (1d0-r)/(1d0-p1)
         x = (1d0-xmax)**bet + ((1d0-x0)**bet-(1d0-xmax)**bet)*r1
         x = 1d0 - x**(1d0/bet)
         dist= x0**(alf-1d0)*(1d0-x)**(bet-1d0)
      ENDIF
      djac=(q1+q2)/dist
      RETURN
  900 WRITE(*,*) ' ========= STOP in chbin1: wrong params'
      STOP
      END


      SUBROUTINE vesko1(mmode,funsko,xpar,ypar)
*     **********************************************
*======================================================================
*======================================================================
*===================== V E S K O 1 ====================================
*==================S. Jadach  September 1985===========================
*==================S. Jadach  November  1991===========================
*==================S. Jadach  May       1997===========================
*======================================================================
* One dimensional Monte Carlo  sampler. Version with weighted events!
*   Mode = -1,0,1 corresponds to initialization, production, finalization.
*   Funsko is the distribution to be generated.
*   Xpar is for input params
*   Ypar is for output params
*   Xpar and Ypar are passed down to funsko
***
* jlim1 is the number of entries in the equidistant latice which
* is formed in the first stage and jlim2 is the total maximum
* number of entries in the latice.
* For mild funsko jlim2=128 is enough.
*     **********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      SAVE
      DIMENSION xpar(*),ypar(*)
      PARAMETER(jlim1=64,jlim2=1024)
      DIMENSION xx(jlim2+1),yy(jlim2+1),zint(jlim2+1)
      EXTERNAL funsko
      DIMENSION drvec(10)
      DATA iwarm/0/
*
      mode=mmode
      IF(mode .EQ. -1) THEN
*     ===================
* initialisation part, see vinsko for more comments
         iniran=1
         iwarm=1
         wt=0.
         swt=0.
         sswt=0.
         nevs=0
* initialisation part, sampling distribution funsko
* and filling matrices xx,yy,zint etc.
         jmax=1
         xx(1)=0.
         xx(2)=1.
         yy(1)=funsko(xx(1),xpar,ypar)
         yy(2)=funsko(xx(2),xpar,ypar)
         IF(yy(1) .LT. 0.0 .OR. yy(2) .LT. 0.0) go to 999
         zint(1)=.5d0*(yy(2)+yy(1))*(xx(2)-xx(1))
*
         jdiv=1
         DO k=1,jlim2-1
            IF(jmax .LT. jlim1) THEN
*...    note that divsko1 increments jmax=jmax+1 in every CALL
               CALL divsko1(funsko,xpar,ypar,jdiv,jmax,xx,yy,zint)
               jdiv=jdiv+2
               IF(jdiv .GT. jmax) jdiv=1
            ELSE
               jdiv=1
               zmx=zint(1)
               DO j=1,jmax
                  IF(zmx .LT. zint(j)) THEN
                     zmx=zint(j)
                     jdiv=j
                  ENDIF
               ENDDO
               CALL divsko1(funsko,xpar,ypar,jdiv,jmax,xx,yy,zint)
            ENDIF
         ENDDO
*
*...  final administration, normalizing zint etc.
         zsum1=0.
         zsum =0.
         DO j=1,jmax
            zsum1=zsum1+zint(j)
            ymax= max( yy(j+1),yy(j))
            zint(j)=ymax*(xx(j+1)-xx(j))
            zsum=zsum+zint(j)
         ENDDO
         sum=0.
         DO j=1,jmax
            sum=sum+zint(j)
            zint(j)=sum/zsum
         ENDDO
         ypar(1)=  zsum
         ypar(2)=  zsum
         ypar(3)=  zsum

      ELSE IF(mode .EQ. 0) THEN
*     =======================
* generation part
         IF(iwarm .EQ. 0) GOTO 901
         CALL varran(drvec,1)
         rnumb = drvec(1)
         DO j=1,jmax
            jstop=j
            IF(zint(j) .GT. rnumb) GOTO 216
         ENDDO
 216     CONTINUE
         IF(jstop .EQ. 1) THEN
            d=rnumb/zint(1)
         ELSE
            d =(rnumb-zint(jstop-1))/(zint(jstop)-zint(jstop-1))
         ENDIF
         x=xx(jstop)*(1.d0 -d )+xx(jstop+1)*d
         fn=funsko(x,xpar,ypar)
         IF(fn .LT. 0.d0) GOTO 999
         yymax=max(yy(jstop+1),yy(jstop))
         wt=fn/yymax
         nevs=nevs+1
         swt=swt+wt
         sswt=sswt+wt*wt
         ypar(1)=  x
         ypar(2)=  fn
         ypar(3)=  wt
*
      ELSE IF(mode .EQ. 1) THEN
*     =======================
         cinteg=0d0
         errint=0d0
         IF(nevs .GT. 0) cinteg=zsum*swt/float(nevs)
         IF(nevs .GT. 0) errint=sqrt(sswt/swt**2-1.d0/float(nevs))
         ypar(1)=  cinteg
         ypar(2)=  errint
         ypar(3)=  zsum
*
      ELSE
*     ====
         GOTO  902
      ENDIF
*     =====
*
      RETURN
 901  WRITE(*,'(a)') ' **** STOP in vesko1, lack of initialisation'
      STOP
 902  WRITE(*,'(a)') ' **** STOP in vesko1, wrong mode '
      STOP
 999  WRITE(*,'(a)') ' **** STOP in vesk01, negative value of funsko '
      STOP
      END

      SUBROUTINE divsko1(funsko,xpar,ypar,jdiv,jmax,xx,yy,zint)
*     ********************************************************
* this routine belongs to vesko1 package
* it sudivides into two equal parts the interval
* (xx(jdiv),xx(jdiv+1))  in the 1-dim. latice
*     ************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      SAVE
      DIMENSION  xx(*),yy(*),zint(*)
      DIMENSION  xpar(*),ypar(*)
      EXTERNAL funsko
*
      xnew=.5d0*(xx(jdiv) +xx(jdiv+1))
      DO j=jmax,jdiv,-1
         xx(j+2)  =xx(j+1)
         yy(j+2)  =yy(j+1)
         zint(j+1)=zint(j)
      ENDDO
      xx(jdiv+1)= xnew
      yy(jdiv+1)= funsko(xnew,xpar,ypar)
      IF(yy(jdiv+1) .LT. 0.) GOTO 999
      zint(jdiv)  =.5d0*(yy(jdiv+1)+yy(jdiv)  )*(xx(jdiv+1)-xx(jdiv)  )
      zint(jdiv+1)=.5d0*(yy(jdiv+2)+yy(jdiv+1))*(xx(jdiv+2)-xx(jdiv+1))
      jmax=jmax+1
      RETURN
 999  CONTINUE
      WRITE(*,'(a)') ' **** STOP in divsko1, negative value of funsko '
      STOP
      END


      SUBROUTINE gausjad(fun,xpar,ypar,aa,bb,eeps,result)
*     ***************************************************
* Gauss-type integration by S. Jadach, Oct. 1990, June 1997
* this is non-adaptive (!!!!) unoptimized (!!!) integration subprogram.
*     *************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION wg(12),xx(12)
      EXTERNAL fun
      DIMENSION xpar(*),ypar(*)
      SAVE
      DATA wg
     $/0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0,
     $ 0.362683783378362d0, 0.027152459411754d0, 0.062253523938648d0,
     $ 0.095158511682493d0, 0.124628971255534d0, 0.149595988816577d0,
     $ 0.169156519395003d0, 0.182603415044924d0, 0.189450610455069d0/
      DATA xx
     $/0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0,
     $ 0.183434642495650d0, 0.989400934991650d0, 0.944575023073233d0,
     $ 0.865631202387832d0, 0.755404408355003d0, 0.617876244402644d0,
     $ 0.458016777657227d0, 0.281603550779259d0, 0.095012509837637d0/
      DATA itermx / 15/
*
      a  = aa
      b  = bb
      eps= eeps
      ndivi=1
* iteration over subdivisions terminated by precision requirement
      DO iter=1,itermx
         calk8  =0d0
         calk16 =0d0
*     sum over delta subintegrals
         DO k = 1,ndivi
            delta = (b-a)/ndivi
            x1    =  a + (k-1)*delta
            x2    =  x1+ delta
            xmidle= 0.5d0*(x2+x1)
            range = 0.5d0*(x2-x1)
            sum8 =0d0
            sum16=0d0
*     8- and 12-point   gauss integration over single delta subinterval
            DO i=1,12
               xplus= xmidle+range*xx(i)
               xminu= xmidle-range*xx(i)
               fplus=fun(xplus,xpar,ypar)
               fminu=fun(xminu,xpar,ypar)
               IF(i .LE. 4) THEN
                  sum8 =sum8  +(fplus+fminu)*wg(i)/2d0
               ELSE
                  sum16=sum16 +(fplus+fminu)*wg(i)/2d0
               ENDIF
            ENDDO
            calk8 = calk8 + sum8 *(x2-x1)
            calk16= calk16+ sum16*(x2-x1)
         ENDDO
         erabs = abs(calk16-calk8)
         erela = 0d0
         IF(calk16.ne.0d0) erela= erabs/abs(calk16)
*     WRITE(*,*) 'gausjad: calk8,calk16=',iter,calk8,calk16,erela
*     precision check to terminate integration
         IF(eeps .GT. 0d0) THEN
            IF(erabs .LT.  eps) GOTO 800
         ELSE
            IF(erela .LT.  eps) GOTO 800
         ENDIF
         ndivi=ndivi*2
      ENDDO
      WRITE(*,*) ' +++++ gausjad:  required precision to high!'
      WRITE(*,*) ' +++++ gausjad:  iter,erela=',iter,erela
  800 CONTINUE
      result = calk16
      END





      SUBROUTINE vesk2w(mode,funsko,x,y,wt)
*     *************************************
*=======================================================================
*=======================================================================
*=======================================================================
*===============TWO DIMENSIONAL SAMPLER VESK2W==========================
*=======================================================================
*=======================================================================
*=======================================================================
*                         vesk2w                                       c
*  general purpose routine to generate an arbitrary two DIMENSIONal    c
*  distribution supplied by user in a form of FUNCTION funsko(x,y)     c
*                 written november 1985                                c
*                    by S. Jadach                                      c
*                 last update:  07.nov.1990                            c
*                 version with weighted event....                      c
*======================================================================c
* vesko2 generates two DIMENSIONal distribution defined by arbitrary
* FUNCTION funsko(x,y) where x,y belong  to (0,1) range.
* the method consists in dividing unit plaquet into cells using
* sort of 'life-game' method in which the division of a cells is made
* (during initialisation) always for this cell which contains
* a maximum value of the integral over funsko in the cell.
* resulting cells contain (usually up to factor two) equal intergral
* value. the generation consists in choosing randomly  a cell
* according to its content and THEN in generating x,y within the cell.
* rejection method is applied at the END of the procedure in order to
* assure that x,y are distributed precisely according to funsko(x,y)
*                    PARAMETERs
* -/ mode = -1 initialisation, no (x,y) generated, CALL vesko2(-1,d1,d2)
*    has to be made prior  to generating first (x,y) pair
* -/ mode =  0 generation of (x,y) pair by CALL vesko2(0,x,y)
* -/ mode =  1 CALL vesko2(1,valint,errint) may be DOne after last
*    (x,y) was generated in order to obtain the value of the integral
*    valint and its error errint, integral is calculated using average
*    weights encoutered during generation phase
* -/ x,y  IF mode=-1 the they are dummy
*         IF mode= 0 the result of random generation according to
*                    FUNCTION funsko, x and y belong to (0,1)
*         IF mode= 1 x= value of integral and y=error (relative)
*                    wt = crude x-section
* ------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      SAVE
      PARAMETER( jlim1 = 64, jlim2 = 1000 , nout = 6 )
      COMMON / vesw2  / xx(jlim2,2),dx(jlim2,2),yy(jlim2,2,2)
     $  ,yymx(jlim2),zint(jlim2),zsum,lev(jlim2),jmax
      DOUBLE PRECISION drvec(100)
      EXTERNAL funsko
      DATA iwarm/77/

      IF(mode) 100,200,300
*...  initialisation part, see vinsko for more comments
  100 CALL vinskw(funsko)
      iwarm=0
      wt=0d0
      wtmax = 1d0
      wtmxx = wtmax
      nevov=0
      swt=0d0
      sswt=0d0
      nevs=0
*(((((((((((((
*     CALL hbook1(1, 16h wt-vesko2     $,75,0.0d0,1.5d0)
*     CALL hminim(1,0)
*     CALL hbook2(2,16h x-y vesko2    $, 64,0,1, 32,0,1,0)
*     CALL hscale(2)
*))))))))))))
      RETURN
*...
  200 CONTINUE
*...  generation part
      IF(iwarm .EQ. 77) go to 980
*c    IF(wt .GT. wtmax) THEN
*c      WRITE(6,*) ' vesko2: ev. overweighted, DOnt worry, wt=',wt
*c      wt=wt-wtmax
*c      nevov=nevov+1
*c    ELSE
        CALL varran(drvec,3)
        r = drvec(1)
        DO 215 j=1,jmax
        jstop=j
  215   IF(zint(j) .GT. r) GOTO 216
  216   CONTINUE
        xr=xx(jstop,1)+dx(jstop,1)*drvec(2)
        yr=xx(jstop,2)+dx(jstop,2)*drvec(3)
        fn=funsko(xr,yr)
        IF(fn .LT. 0.) GOTO 999
        yymax=yymx(jstop)
        wt=fn/yymax
        wtmxx = max(wtmxx,wt)
*c      IF(nevs .LE. (4*jlim2) .AND. wt .GT. wtmax) THEN
*c         wtmax=wt*1.1d0
*c         WRITE(6,*) ' vesko2: nevs, new wtmax= ',nevs,wtmax
*c      ENDIF
        nevs=nevs+1
        swt=swt+wt
        sswt=sswt+wt*wt
*((((((((((
*       CALL hfill(1,wt,0d0,1d0)
*))))))))))
*cc   ENDIF
*cc    CALL varran(drvec,1)
*cc    rn=drvec(1)
*cc   IF(wtmax*rn .GT. wt) GOTO 200
      x=xr
      y=yr
*((((((((((
*     CALL hfill(2,xr,yr)
*))))))))))
      RETURN
*...
  300 CONTINUE
* this is the value of the integral
      cinteg=zsum*swt/nevs
* and its error
      errint=sqrt(sswt/swt**2-1d0/nevs)
      x=cinteg
      y=errint
      wt=zsum
*((((((((((
*     CALL hprint(1)
*     CALL hdelet(1)
*     CALL hprint(2)
*     CALL hdelet(2)
      print 7000,nevs,nevov,wtmax,wtmxx
 7000 FORMAT(' vesk2w: nevs,nevov,wtmax,wtmxx= ',2i7,2f7.3)
*))))))))))
      RETURN
  980 WRITE(nout,9002)
 9002 FORMAT(' **** STOP in vesk2w, lack of initialisation   ')
      STOP
  999 WRITE(nout,9004)
 9004 FORMAT(' **** STOP in vesk2w, negative value of funsko ')
      STOP
      END

      SUBROUTINE vinskw(funsko)
*     *************************
* this routine belongs to vesko2 package
* jlim1 is the number of cells, division of the unit plaque into cells
* is made in the first stage.    jlim2 is the total maximum
* number of cells, note that DIMENSIONs of
* matrices in /veskoa/ should be at least jlim2
*     **********************************
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      SAVE
* ------------------------------------------------------------
      PARAMETER( jlim1 = 64, jlim2 = 1000 , nout = 6 )
      COMMON / vesw2  / xx(jlim2,2),dx(jlim2,2),yy(jlim2,2,2)
     $  ,yymx(jlim2),zint(jlim2),zsum,lev(jlim2),jmax
      EXTERNAL funsko

*...  initialisation part, sampling distribution funsko
*...  and filling matrices xx,yy,zint etc.
      jmax=1
      xx(1,1)=0d0
      xx(1,2)=0d0
      dx(1,1)=1d0
      dx(1,2)=1d0
      lev(1)=1
      sum=0d0
      DO 150 i=1,2
      DO 150 k=1,2
*... this is not elegant but simple
      yy(1,i,k)=funsko(xx(1,1)+(i-1.)*dx(1,1),xx(1,2)+(k-1.)*dx(1,2))
      IF(yy(1,i,k) .LT. 0.0) go to 999
  150 sum=sum+yy(1,i,k)
      zint(1)=sum*dx(1,1)*dx(1,2)/4d0

      jdiv=1
      DO 200 kk=1,jlim2-1
      IF(jmax .LT. jlim1) THEN
*...    note that divskw increments jmax=jmax+1 in every CALL
        CALL divskw(jdiv,funsko)
*(((((((((((
*      IF(jmax .EQ. jlim1) THEN
*      print 9900,jmax,(lev(i),i=1,jmax)
* 9900 FORMAT(///,' jmax...  lev lev lev lev lev',i10,/(24i5))
*      print 9901,((xx(jd,i),i=1,2),jd=1,jmax)
* 9901 FORMAT('  xx xx xx xx xx xx xx  ',/(10e12.5))
*      print 9902,((dx(jd,i),i=1,2),jd=1,jmax)
* 9902 FORMAT('  dx  dx dx dx dx dx ',/(10e12.5))
*      print 9903,(((yy(jd,i,k),i=1,2),k=1,2),jd=1,jmax)
* 9903 FORMAT('  yy  yy yy yy yy yy ',/(8e15.5))
*      print 9904,(zint(i),i=1,jmax)
* 9904 FORMAT('   zint zint zint zint ',/(10e12.5))
*      ENDIF
*))))))))))))
        jdiv=jdiv+2
        IF(jdiv .GT. jmax) jdiv=1
      ELSE
        jdiv=1
        zmx=zint(1)
        DO 180 j=1,jmax
        IF(zmx .LT. zint(j)) THEN
          zmx=zint(j)
          jdiv=j
        ENDIF
  180   CONTINUE
        CALL divskw(jdiv,funsko)
      ENDIF
  200 CONTINUE

*(((((((((((
*      jprn=64
*      print 9910,jmax,(lev(i),i=1,jmax)
* 9910 FORMAT(/,' jmax...  lev lev lev lev lev',i10,/(24i5))
*      IF(jmax .LE. jprn) print 9911,((xx(jd,i),i=1,2),jd=1,jmax)
* 9911 FORMAT('  xx xx xx xx xx xx xx  ',/(10e12.5))
*      IF(jmax .LE. jprn) print 9912,((dx(jd,i),i=1,2),jd=1,jmax)
* 9912 FORMAT('  dx  dx dx dx dx dx ',/(10e12.5))
*      IF(jmax .LE. jprn) print 9913,
*     $             (((yy(jd,i,k),i=1,2),k=1,2),jd=1,jmax)
* 9913 FORMAT('  yy  yy yy yy yy yy ',/(8e15.5))
*      IF(jmax .LE. jprn) print 9914,(zint(i),i=1,jmax)
* 9914 FORMAT('   zint zint zint zint ',/(10e12.5))
*     DO 902 j=1,jmax
*     z=1d0*j-.5d0
* 902 CALL hfill(202,z,zint(j))
*))))))))))))
*...  final administration, normalizing zint etc.
      zsum1=0d0
      zsum =0d0
      DO 260 j=1,jmax
      zsum1=zsum1+zint(j)
      ymax= 0d0
      DO 250 i=1,2
      DO 250 k=1,2
  250 ymax= max(ymax,yy(j,i,k))
      yymx(j)=ymax
      zint(j)=ymax*dx(j,1)*dx(j,2)
  260 zsum=zsum+zint(j)
*((((((((
      zr=zsum1/zsum
      print 7000,zr
 7000 FORMAT(' /////// zsum1/zsum= ',f20.8)
*)))))))))
      sum=0d0
      DO 240 j=1,jmax
      sum=sum+zint(j)
  240 zint(j)=sum/zsum
*(((((((((((
*     jprn=64
*     print 9932,jmax
*9932 FORMAT(/'=====jmax zint zint zint  ',i10)
*     IF(jmax .LE. jprn) print 9935,(zint(i),i=1,jmax)
*9935            FORMAT(10e12.5)
*     DO 901 j=2,jmax
* 901 CALL hfill(201,(zint(j)-zint(j-1))*jmax)
*     CALL hfill(201,zint(1)*jmax)
*))))))))))))
      RETURN
  999 WRITE(nout,9000)
 9000 FORMAT(' **** STOP in vinskw, negative value of funsko ')
      STOP
      END

      SUBROUTINE divskw(jd,funsko)
*     ****************************
* this routine belongs to vesko2 package
* it subdivides one cell (no. jd) into two equal size cells
*     **********************************
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      SAVE
* ------------------------------------------------------------
      PARAMETER( jlim1 = 64, jlim2 = 1000 , nout = 6 )
      COMMON / vesw2  / xx(jlim2,2),dx(jlim2,2),yy(jlim2,2,2)
     $  ,yymx(jlim2),zint(jlim2),zsum,lev(jlim2),jmax
      EXTERNAL funsko

*...  moove to make a hole for a new entry (one additional cell)
      DO 100 j=jmax,jd,-1
      zint(j+1)=zint(j)
      lev(j+1)=lev(j)
      DO 100 i=1,2
      xx(j+1,i)  =xx(j,i)
      dx(j+1,i)  =dx(j,i)
      DO 100 k=1,2
  100 yy(j+1,i,k)  =yy(j,i,k)
*...  create two new cells and store them
      ll= MOD(lev(jd),2)+1
      dx(jd,ll)=dx(jd,ll)/2d0
      dx(jd+1,ll)=dx(jd+1,ll)/2d0
      xx(jd+1,ll)=xx(jd,ll)+dx(jd,ll)
      IF(ll .EQ. 1) THEN
        DO 150 i=1,2
*... this is not elegant, probably could be DOne better
        yy(jd,2,i)=funsko(xx(jd,1)+dx(jd,1),xx(jd,2)+(i-1.)*dx(jd,2))
  150   yy(jd+1,1,i)=yy(jd,2,i)
      ELSE
        DO 152 i=1,2
        yy(jd,i,2)=funsko(xx(jd,1)+(i-1.)*dx(jd,1),xx(jd,2)+dx(jd,2))
  152   yy(jd+1,i,1)=yy(jd,i,2)
      ENDIF
*...  estimate the integrals over new cells resulting from division
      DO 220 jdv=jd,jd+1
      lev(jdv)=lev(jdv)+1
      sum=0d0
      DO 210 i=1,2
      DO 210 k=1,2
      IF(yy(jdv,i,k) .LT. 0.d0) go to 999
  210 sum=sum+yy(jdv,i,k)
  220 zint(jdv) =sum*dx(jdv,1)*dx(jdv,2)/4d0
      jmax=jmax+1
      RETURN
  999 WRITE(nout,9000)
 9000 FORMAT(' **** STOP in divskw, negative value of funsko ')
      STOP
      END


      SUBROUTINE gausjd(fun,aa,bb,eeps,result)
*     ****************************************
* gauss integration by s. jadach, oct. 90.
* this is non-adaptive (!!!!) unoptimized (!!!) integration subprogram.
*     *************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION wg(12),xx(12)
      COMMON / inout  / ninp,nout
      EXTERNAL fun
      SAVE /inout/,wg,xx,itermx
      DATA wg
     $/0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0,
     $ 0.362683783378362d0, 0.027152459411754d0, 0.062253523938648d0,
     $ 0.095158511682493d0, 0.124628971255534d0, 0.149595988816577d0,
     $ 0.169156519395003d0, 0.182603415044924d0, 0.189450610455069d0/
      DATA xx
     $/0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0,
     $ 0.183434642495650d0, 0.989400934991650d0, 0.944575023073233d0,
     $ 0.865631202387832d0, 0.755404408355003d0, 0.617876244402644d0,
     $ 0.458016777657227d0, 0.281603550779259d0, 0.095012509837637d0/
      DATA itermx / 15/
      eps=abs(eeps)
      a=aa
      b=bb
      ndivi=1
* iteration over subdivisions terminated by precision requirement
      DO 400 iter=1,itermx
      calk8  =0d0
      calk16 =0d0
* sum over delta subintegrals
      DO 200 k = 1,ndivi
      delta = (b-a)/ndivi
      x1    =  a + (k-1)*delta
      x2    =  x1+ delta
      xmidle= 0.5d0*(x2+x1)
      range = 0.5d0*(x2-x1)
      sum8 =0d0
      sum16=0d0
* 8- and 12-point   gauss integration over single delta subinterval
      DO 100 i=1,12
      xplus= xmidle+range*xx(i)
      xminu= xmidle-range*xx(i)
      fplus=fun(xplus)
      fminu=fun(xminu)
      IF(i .LE. 4) THEN
          sum8 =sum8  +(fplus+fminu)*wg(i)/2d0
      ELSE
          sum16=sum16 +(fplus+fminu)*wg(i)/2d0
      ENDIF
  100 CONTINUE
      calk8 = calk8 + sum8 *(x2-x1)
      calk16= calk16+ sum16*(x2-x1)
  200 CONTINUE
      erabs = abs(calk16-calk8)
      erela = 0d0
      IF(calk16.ne.0d0) erela= erabs/abs(calk16)
*     WRITE(6,*) 'gausjd: calk8,calk16=',iter,calk8,calk16,erela
* PRECISION check to terminate integration
      IF(eeps .GT. 0d0) THEN
        IF(erabs .LT.  eps) GOTO 800
      ELSE
        IF(erela .LT.  eps) GOTO 800
      ENDIF
  400 ndivi=ndivi*2
      WRITE(nout,*) ' +++++ gausjd:  required precision to high!'
      WRITE(nout,*) ' +++++ gausjd:  iter,erela=',iter,erela
  800 result= calk16
      END


      SUBROUTINE wmonit(mode,id,wt,wtmax,rn)
*     **************************************
* last correction 19 sept. 89
* utility PROGRAM for monitoring m.c. rejection weights.
* id is weight idendifier, maximum idmx (defined below).
* wt is weight, wtmax is maximum weight and rn is random number.
* IF(mode .EQ. -1) THEN
*          initalization IF entry id, other arguments are ignored
* ELSEIF(mode .EQ. 0) THEN
*          summing up weights etc. for a given event for entry id
*        - wt is current weight.
*        - wtmax is maximum weight used for couting overweighted
*          events with wt>wtmax.
*        - rn is random number used in rejection, it is used to
*          count no. of accepted (rn<wt/wtmax) and rejected
*          (wt>wt/wtmax) events,
*          IF ro rejection THEN put rn=0d0.
* ELSEIF(mode .EQ. 1) THEN
*          in this mode wmonit repports on accumulated statistics
*          and the inFORMATion is stored in COMMON /cmonit/
*        - averwt= average weight wt counting all event
*        - errela= relative error of averwt
*        - nevtot= total nimber of accounted events
*        - nevacc= no. of accepted events (rn<wt\wtmax)
*        - nevneg= no. of events with negative weight (wt<0)
*        - nevzer= no. of events with zero weight (wt .EQ. 0d0)
*        - nevove= no. of overweghted events (wt>wtmax)
*          and IF you DO not want to use cmonit THEN the value
*          the value of averwt is assigned to wt,
*          the value of errela is assigned to wtmax and
*          the value of wtmax  is assigned to rn in this mode.
* ELSEIF(modee .EQ. 2) THEN
*          all inFORMATion defined for entry id defined above
*          for mode=2 is just printed of unit nout
* ENDIF
* note that output repport (mode=1,2) is DOne dynamiCALLy just for a
* given entry id only and it may be repeated many times for one id and
* for various id's as well.
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      SAVE
      PARAMETER(idmx=100)
      COMMON / cmonit/ averwt,errela,nevtot,nevacc,nevneg,nevove,nevzer
      COMMON / inout  / ninp,nout
      INTEGER ntot(idmx),nacc(idmx),nneg(idmx),nove(idmx),nzer(idmx)
      DIMENSION swt(idmx),sswt(idmx),wwmx(idmx)
      DATA ntot /idmx* -1/  swt /idmx*   0d0/
      DATA sswt /idmx*0d0/ wwmx /idmx*-1d-20/
*
      IF(id .LE. 0 .OR. id .GT. idmx) THEN
           WRITE(nout,*) ' =====wmonit: wrong id',id
           STOP
      ENDIF
      IF(mode .EQ. -1) THEN
           ntot(id)=0
           nacc(id)=0
           nneg(id)=0
           nzer(id)=0
           nove(id)=0
           swt(id)   =0d0
           sswt(id)  =0d0
           wwmx(id)  = -1d-20
      ELSEIF(mode .EQ. 0) THEN
           IF(ntot(id) .LT. 0) THEN
              WRITE(nout,*) ' ==== warning from wmonit: '
              WRITE(nout,*) ' lack of initialization, id=',id
           ENDIF
           ntot(id)=ntot(id)+1
           swt(id)=swt(id)+wt
           sswt(id)=sswt(id)+wt**2
           wwmx(id)= max(wwmx(id),wt)
           IF(wt .EQ. 0d0)   nzer(id)=nzer(id)+1
           IF(wt .LT. 0d0)   nneg(id)=nneg(id)+1
           IF(wt .GT. wtmax)      nove(id)=nove(id)+1
           IF(rn*wtmax .LE. wt)   nacc(id)=nacc(id)+1
      ELSEIF(mode .EQ. 1) THEN
           IF(ntot(id) .LT. 0) THEN
              WRITE(nout,*) ' ==== warning from wmonit: '
              WRITE(nout,*) ' lack of initialization, id=',id
           ENDIF
           IF(ntot(id) .LE. 0 .OR. swt(id) .EQ. 0d0)  THEN
              averwt=0d0
              errela=0d0
           ELSE
              averwt=swt(id)/float(ntot(id))
              errela=sqrt(abs(sswt(id)/swt(id)**2-1d0/float(ntot(id))))
           ENDIF
           nevtot=ntot(id)
           nevacc=nacc(id)
           nevneg=nneg(id)
           nevzer=nzer(id)
           nevove=nove(id)
           wt=averwt
           wtmax=errela
           rn    =wwmx(id)
      ELSEIF(mode .EQ. 2) THEN
           IF(ntot(id) .LE. 0 .OR. swt(id) .EQ. 0d0)  THEN
              averwt=0d0
              errela=0d0
           ELSE
              averwt=swt(id)/float(ntot(id))
              errela=sqrt(abs(sswt(id)/swt(id)**2-1d0/float(ntot(id))))
              wwmax=wwmx(id)
           ENDIF
           WRITE(nout,1003) id, averwt, errela, wwmax
           WRITE(nout,1004) ntot(id),nacc(id),nneg(id),nove(id),nzer(id)
           wt=averwt
           wtmax=errela
           rn    =wwmx(id)
      ELSE
           WRITE(nout,*) ' =====wmonit: wrong mode',mode
           STOP
      ENDIF
 1003 FORMAT(
     $  ' =======================wmonit========================'
     $/,'   id           averwt         errela            wwmax'
     $/,    i5,           e17.7,         f15.9,           e17.7)
 1004 FORMAT(
     $  ' -----------------------------------------------------------'
     $/,'      nevtot      nevacc      nevneg      nevove      nevzer'
     $/,   5i12)
      END


      FUNCTION gaus(f,a,b,eeps)
*     *************************
* this is iterative integration procedure
* originates  probably from cern library
* it subdivides inegration range until required PRECISION is reached
* PRECISION is a difference from 8 and 16 point gauss itegr. result
* eeps positive treated as absolute PRECISION
* eeps negative treated as relative PRECISION
*     *************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION w(12),x(12)
      COMMON / inout  / ninp,nout
      EXTERNAL f
      DATA const /1.0d-19/
      SAVE     / inout/, const, w, x
      DATA w
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA x
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      eps=abs(eeps)
      delta=const*abs(a-b)
      gaus=0d0
      aa=a
    5 y=b-aa
      IF(abs(y)  .LE.  delta) RETURN
    2 bb=aa+y
      c1=0.5d0*(aa+bb)
      c2=c1-aa
      s8=0d0
      s16=0d0
      DO 1 i=1,4
      u=x(i)*c2
    1 s8=s8+w(i)*(f(c1+u)+f(c1-u))
      DO 3 i=5,12
      u=x(i)*c2
    3 s16=s16+w(i)*(f(c1+u)+f(c1-u))
      s8=s8*c2
      s16=s16*c2
      IF(eeps .LT. 0d0) THEN
        IF(abs(s16-s8)  .GT.  eps*abs(s16)) go to 4
      ELSE
        IF(abs(s16-s8)  .GT.  eps) go to 4
      ENDIF
      gaus=gaus+s16
      aa=bb
      go to 5
    4 y=0.5d0*y
      IF(abs(y)  .GT.  delta) GOTO 2
      WRITE(nout,7)
      gaus=0d0
      RETURN
    7 FORMAT(1x,36hgaus  ... too high accuracy required)
      END


      DOUBLE PRECISION FUNCTION dilogy(x)
*-------------------------------------------- remarks ---------------
* dilogarithm FUNCTION: dilog(x)=int( -ln(1-z)/z ) , 0 < z < x .
* this is the cernlib version.
*--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      z=-1.644934066848226d0
      IF(x  .LT. -1.d0) go to 1
      IF(x  .LE.  0.5d0) go to 2
      IF(x  .EQ.  1.d0) go to 3
      IF(x  .LE.  2.d0) go to 4
      z=3.289868133696453d0
    1 t=1.d0/x
      s=-0.5d0
      z=z-0.5d0*dlog(dabs(x))**2
      go to 5
    2 t=x
      s=0.5d0
      z=0.d0
      go to 5
    3 dilogy=1.644934066848226d0
      RETURN
    4 t=1.d0-x
      s=-0.5d0
      z=1.644934066848226d0-dlog(x)*dlog(dabs(t))
    5 y=2.666666666666667d0*t+0.666666666666667d0
      b=      0.000000000000001d0
      a=y*b  +0.000000000000004d0
      b=y*a-b+0.000000000000011d0
      a=y*b-a+0.000000000000037d0
      b=y*a-b+0.000000000000121d0
      a=y*b-a+0.000000000000398d0
      b=y*a-b+0.000000000001312d0
      a=y*b-a+0.000000000004342d0
      b=y*a-b+0.000000000014437d0
      a=y*b-a+0.000000000048274d0
      b=y*a-b+0.000000000162421d0
      a=y*b-a+0.000000000550291d0
      b=y*a-b+0.000000001879117d0
      a=y*b-a+0.000000006474338d0
      b=y*a-b+0.000000022536705d0
      a=y*b-a+0.000000079387055d0
      b=y*a-b+0.000000283575385d0
      a=y*b-a+0.000001029904264d0
      b=y*a-b+0.000003816329463d0
      a=y*b-a+0.000014496300557d0
      b=y*a-b+0.000056817822718d0
      a=y*b-a+0.000232002196094d0
      b=y*a-b+0.001001627496164d0
      a=y*b-a+0.004686361959447d0
      b=y*a-b+0.024879322924228d0
      a=y*b-a+0.166073032927855d0
      a=y*a-b+1.935064300869969d0
      dilogy=s*t*(a-b)+z
      END


      DOUBLE PRECISION FUNCTION dpgamm(z)
*     **********************************
* DOUBLE PRECISION gamma FUNCTION
      DOUBLE PRECISION z,z1,x,x1,x2,d1,d2,s1,s2,s3,pi,c(20),const
      SAVE c,pi,const
      DATA c( 1) / 8.3333333333333333333333333332d-02/
      DATA c( 2) /-2.7777777777777777777777777777d-03/
      DATA c( 3) / 7.9365079365079365079365079364d-04/
      DATA c( 4) /-5.9523809523809523809523809523d-04/
      DATA c( 5) / 8.4175084175084175084175084175d-04/
      DATA c( 6) /-1.9175269175269175269175269175d-03/
      DATA c( 7) / 6.4102564102564102564102564102d-03/
      DATA c( 8) /-2.9550653594771241830065359477d-02/
      DATA c( 9) / 1.7964437236883057316493849001d-01/
      DATA c(10) /-1.3924322169059011164274322169d+00/
      DATA c(11) / 1.3402864044168391994478951001d+01/
      DATA c(12) /-1.5684828462600201730636513245d+02/
      DATA c(13) / 2.1931033333333333333333333333d+03/
      DATA c(14) /-3.6108771253724989357173265219d+04/
      DATA c(15) / 6.9147226885131306710839525077d+05/
      DATA c(16) /-1.5238221539407416192283364959d+07/
      DATA c(17) / 3.8290075139141414141414141414d+08/
      DATA c(18) /-1.0882266035784391089015149165d+10/
      DATA c(19) / 3.4732028376500225225225225224d+11/
      DATA c(20) /-1.2369602142269274454251710349d+13/
      DATA pi    / 3.1415926535897932384626433832d+00/
      DATA const / 9.1893853320467274178032973641d-01/
      IF(z .GT. 5.75d 1)                                     GOTO  6666
      nn = z
      IF (z  -  dble(float(nn)))                 3,1,3
    1 IF (z     .LE.     0.d 0)                    GOTO 6667
      dpgamm = 1.d 0
      IF (z     .LE.     2.d 0)                    RETURN
      z1 = z
    2 z1 = z1  -  1.d 0
      dpgamm = dpgamm * z1
      IF (z1  -  2.d 0)                          61,61,2
    3 IF (dabs(z)     .LT.     1.d-29)             GOTO 60
      IF (z     .LT.     0.d 0)                    GOTO 4
      x  = z
      kk = 1
      GOTO 10
    4 x  = 1.d 0  -  z
      kk = 2
   10 x1 = x
      IF (x     .GT.     19.d 0)                   GOTO 13
      d1 = x
   11 x1 = x1  +  1.d 0
      IF (x1     .GE.     19.d 0)                  GOTO 12
      d1 = d1 * x1
      GOTO 11
   12 s3 = -dlog(d1)
      GOTO 14
   13 s3 = 0.d 0
   14 d1 = x1 * x1
      s1 = (x1  -  5.d-1) * dlog(x1)  -  x1  +  const
      DO 20                  k=1,20
      s2 = s1  +  c(k)/x1
      IF (dabs(s2  -  s1)     .LT.     1.d-28)     GOTO 21
      x1 = x1 * d1
   20 s1 = s2
   21 s3 = s3  +  s2
      GOTO (50,22),    kk
   22 d2 = dabs(z  -  nn)
      d1 = d2 * pi
      IF (d1     .LT.     1.d-15)                  GOTO 31
   30 x2 =  dlog(pi/dsin(d1))  -  s3
      GOTO 40
   31 x2 = -dlog(d2)
   40 mm = dabs(z)
      IF(x2       .GT.       1.74d2)                  GOTO 6666
      dpgamm = dexp(x2)
      IF (mm    .ne.    (mm/2) * 2)              RETURN
      dpgamm = -dpgamm
      RETURN
   50 IF(s3       .GT.       1.74d2)                  GOTO 6666
      dpgamm = dexp(s3)
      RETURN
 6666 print *, 2000
      RETURN
 6667 print *, 2001
      RETURN
   60 dpgamm = 0.d 0
      IF(dabs(z)    .LT.    1.d-77)   RETURN
      dpgamm = 1.d 0/z
   61 RETURN
 2000 FORMAT (/////, 2x, 32hdpgamm ..... argument too large., /////)
 2001 FORMAT (/////, 2x, 32hdpgamm ..... argument is a pole., /////)
      END




*=======================================================================
*=======================================================================
*=======================================================================
*==received: by dxmint.cern.ch (cernvax) (5.57/3.14)
*== id aa13405; wed, 23 jan 91 17:19:06 +0100
*==message-id: <9101231619.aa13405@dxmint.cern.ch>
*==received: by cernapo; wed, 23 jan 91 17:23:40 +0100
*==received: by apojames.cern.ch; wed, 23 jan 91 17:05:23 cet
*==date: wed, 23 jan 91 17:05:23 cet
*==from: james@cernapo.cern.ch (frederick james)
*==to: jadach@cernvm
*==subject: random generators
*==
*==      PROGRAM pseudoran
*==c  cpc # abtk                                           cpc # abtk
*==c         pseudorandom generator demonstration (test case)
*==      DIMENSION rvec(1000)
*==      DIMENSION veri(5), isd25(25)
*==c
*==c
*==c   ................................................
*==      WRITE(6,'(20x,a)') 'demonstration of pseudorandom generators'
*==      WRITE(6,'(20x,a)') 'machine/system: date:'
*==      WRITE(6,'(/20x,a/)') 'initialization and test of portability'
*==c   ................................................
*==c
*==c                   initialization and verification  ranmar
*==        DO 40 i9= 1, 20
*==   40   CALL ranmar(rvec,1000)
*==      CALL ranmar(rvec,5)
*==      DO 41 i= 1 ,5
*==   41 veri(i) = (4096.*rvec(i))*(4096.)
*==      WRITE(6,'(a,5f12.1/)') '  ranmar 20001  ',veri
*==c
*==c                   initialization and verification  ranecu
*==      CALL ranecu(rvec,1000)
*==      CALL ranecu(veri,5)
*==      DO 52 i= 1 ,5
*==   52 veri(i) = 4096.*(4096.*veri(i))
*==      WRITE(6,'(a,5f12.1/)') '  ranecu 1001   ',veri
*==c
*==c                   initialization and verification  rcarry
*==      CALL rcarry(rvec,1000)
*==      CALL rcarry(veri,5)
*==      DO 62 i= 1 ,5
*==   62 veri(i) = 4096.*(4096.*veri(i))
*==      WRITE(6,'(a,5f12.1/)') '  rcarry 1001   ',veri
*==c
*==      WRITE(6,'(//20x,a/)') 'test of repeatability'
*==c  .................................................
*==c                  verify restarting      ranmar
*==      WRITE(6,'(/a)') '   the next line should be repeated:'
*==      CALL rmarut(imar1,imar2,imar3)
*==      CALL ranmar(rvec,777)
*==      CALL ranmar(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       ranmar 1 ',veri
*==      CALL rmarin(imar1,imar2,imar3)
*==      CALL ranmar(rvec,777)
*==      CALL ranmar(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       ranmar 2 ',veri
*==c
*==c                  verify restarting      ranecu
*==      WRITE(6,'(/a)') '   the next line should be repeated:'
*==      CALL recuut(is1,is2)
*==      CALL ranecu(rvec,777)
*==      CALL ranecu(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       ranecu 1 ',veri
*==      CALL recuin(is1,is2)
*==      CALL ranecu(rvec,777)
*==      CALL ranecu(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       ranecu 2 ',veri
*==c
*==c                  verify restarting      rcarry
*==      WRITE(6,'(/a)') '   the next line should be repeated:'
*==      CALL rcarut(isd25)
*==      CALL rcarry(rvec,777)
*==      CALL rcarry(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       rcarry 1 ',veri
*==      CALL rcarin(isd25)
*==      CALL rcarry(rvec,777)
*==      CALL rcarry(veri,5)
*==      WRITE(6,'(a,5f12.9)') '       rcarry 2 ',veri
*==c
*==      STOP
*==      END
*=======================================================================
*=======================================================================
*=======================================================================
      SUBROUTINE marran(rvec,lenv)
* =======================s. jadach===================================
* == this commes from f. james, the name of ranmar is changed to   ==
* == marran in order to avoid interference with the version        ==
* == alREADy in use and the public library version (if present).   ==
* ==      this is the only modification !!!!                       ==
* ========================s. jadach==================================
* universal random number generator proposed by marsaglia and zaman
* in report fsu-scri-87-50
*        modified by f. james, 1988 and 1989, to generate a vector
*        of pseudorandom numbers rvec of length lenv, and to put in
*        the COMMON block everything needed to specify currrent state,
*        and to add input and output entry points marini, marout.
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*!!!  CALLing sequences for ranmar:                                  ++
*!!!      CALL ranmar (rvec, len)   RETURNs a vector rvec of len     ++
*!!!                   32-bit random floating point numbers between  ++
*!!!                   zero and one.                                 ++
*!!!      CALL marini(i1,n1,n2)   initializes the generator from one ++
*!!!                   32-bit INTEGER i1, and number counts n1,n2    ++
*!!!                  (for initializing, set n1=n2=0, but to restart ++
*!!!                    a previously generated sequence, use values  ++
*!!!                    output by marout)                            ++
*!!!      CALL marout(i1,n1,n2)   outputs the value of the original  ++
*!!!                  seed and the two number counts, to be used     ++
*!!!                  for restarting by initializing to i1 and       ++
*!!!                  skipping n2*100000000+n1 numbers.              ++
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION rvec(*)
      COMMON/raset1/u(97),c,i97,j97
      PARAMETER (modcns=1000000000)
      SAVE cd, cm, twom24, ntot, ntot2, ijkl
      DATA ntot,ntot2,ijkl/-1,0,0/
*
      IF (ntot  .GE.  0)  go to 50
*
*        default initialization. user has CALLed ranmar without marini.
      ijkl = 54217137
      ntot = 0
      ntot2 = 0
      kalled = 0
      go to 1
*
      entry      marini(ijklin, ntotin,ntot2n)
*         initializing routine for ranmar, may be CALLed before
*         generating pseudorandom numbers with ranmar. the input
*         values should be in the ranges:  0<=ijklin<=900 ooo ooo
*                                          0<=ntotin<=999 999 999
*                                          0<=ntot2n<<999 999 999!
* to get the standard values in marsaglia's paper, ijklin=54217137
*                                            ntotin,ntot2n=0
      ijkl = ijklin
      ntot = max(ntotin,0)
      ntot2= max(ntot2n,0)
      kalled = 1
*          always come here to initialize
    1 CONTINUE
      ij = ijkl/30082
      kl = ijkl - 30082*ij
      i = MOD(ij/177, 177) + 2
      j = MOD(ij, 177)     + 2
      k = MOD(kl/169, 178) + 1
      l = MOD(kl, 169)
      WRITE(6,'(a,5i10)')
     $'marran initialized: ij,kl,ijkl,ntot,ntot2=',ij,kl,ijkl,ntot,ntot2
      DO 2 ii= 1, 97
      s = 0.
      t = .5
      DO 3 jj= 1, 24
         m = MOD(MOD(i*j,179)*k, 179)
         i = j
         j = k
         k = m
         l = MOD(53*l+1, 169)
         IF (MOD(l*m,64)  .GE.  32)  s = s+t
    3    t = 0.5*t
    2 u(ii) = s
      twom24 = 1.0
      DO 4 i24= 1, 24
    4 twom24 = 0.5*twom24
      c  =   362436.*twom24
      cd =  7654321.*twom24
      cm = 16777213.*twom24
      i97 = 97
      j97 = 33
*       complete initialization by skipping
*            (ntot2*modcns + ntot) random numbers
      DO 45 loop2= 1, ntot2+1
      now = modcns
      IF (loop2  .EQ.  ntot2+1)  now=ntot
      IF (now  .GT.  0)  THEN
        WRITE(6,'(a,i15)') ' marini skipping over ',now
       DO 40 idum = 1, ntot
       uni = u(i97)-u(j97)
       IF (uni  .LT.  0.)  uni=uni+1.
       u(i97) = uni
       i97 = i97-1
       IF (i97  .EQ.  0)  i97=97
       j97 = j97-1
       IF (j97  .EQ.  0)  j97=97
       c = c - cd
       IF (c  .LT.  0.)  c=c+cm
   40  CONTINUE
      ENDIF
   45 CONTINUE
      IF (kalled  .EQ.  1)  RETURN
*
*          normal entry to generate lenv random numbers
   50 CONTINUE
      DO 100 ivec= 1, lenv
      uni = u(i97)-u(j97)
      IF (uni  .LT.  0.)  uni=uni+1.
      u(i97) = uni
      i97 = i97-1
      IF (i97  .EQ.  0)  i97=97
      j97 = j97-1
      IF (j97  .EQ.  0)  j97=97
      c = c - cd
      IF (c  .LT.  0.)  c=c+cm
      uni = uni-c
      IF (uni  .LT.  0.) uni=uni+1.
      rvec(ivec) = uni
*             replace exact zeros by uniform distr. *2**-24
         IF (uni  .EQ.  0.)  THEN
         zuni = twom24*u(2)
*             an exact zero here is very unlikely, but let's be safe.
         IF (zuni  .EQ.  0.) zuni= twom24*twom24
         rvec(ivec) = zuni
         ENDIF
  100 CONTINUE
      ntot = ntot + lenv
         IF (ntot  .GE.  modcns)  THEN
         ntot2 = ntot2 + 1
         ntot = ntot - modcns
         ENDIF
      RETURN
*           entry to output current status
      entry marout(ijklut,ntotut,ntot2t)
      ijklut = ijkl
      ntotut = ntot
      ntot2t = ntot2
      RETURN
      END
      SUBROUTINE carran(rvec,lenv)
*         add-and-carry random number generator proposed by
*         marsaglia and zaman in siam j. scientific and statistical
*             computing, to appear probably 1990.
*         modified with enhanced initialization by f. james, 1990
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*!!!  CALLing sequences for carran:                                  ++
*!!!      CALL carran (rvec, len)   RETURNs a vector rvec of len     ++
*!!!                   32-bit random floating point numbers between  ++
*!!!                   zero and one.                                 ++
*!!!      CALL carini(int)     initializes the generator from one    ++
*!!!                   32-bit INTEGER int                            ++
*!!!      CALL carres(ivec)    restarts the generator from vector    ++
*!!!                   ivec of 25 32-bit INTEGERs (see carout)       ++
*!!!      CALL carout(ivec)    outputs the current values of the 25  ++
*!!!                 32-bit INTEGER seeds, to be used for restarting ++
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION rvec(lenv)
      DIMENSION seeds(24), iseeds(24), isdext(25)
      PARAMETER (twop12=4096.)
      PARAMETER (itwo24=2**24, icons=2147483563)
      SAVE notyet, i24, j24, carry, seeds, twom24
      LOGICAL notyet
      DATA notyet/.true./
      DATA i24,j24,carry/24,10,0./
*
*              default initialization by multiplicative congruential
      IF (notyet) THEN
         notyet = .false.
         jseed = 314159265
         WRITE(6,'(a,i12)') ' carran default initialization: ',jseed
            twom24 = 1.
         DO 25 i= 1, 24
            twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         IF (jseed  .LT.  0)  jseed = jseed+icons
         iseeds(i) = MOD(jseed,itwo24)
   25    CONTINUE
         DO 50 i= 1,24
         seeds(i) = REAL(iseeds(i))*twom24
   50    CONTINUE
         i24 = 24
         j24 = 10
         carry = 0.
         IF (seeds(24)  .LT.  seeds(14)) carry = twom24
      ENDIF
*
*          the generator proper: "subtract-with-borrow",
*          as proposed by marsaglia and zaman,
*          florida state university, march, 1989
*
      DO 100 ivec= 1, lenv
      uni = seeds(i24) - seeds(j24) - carry
      IF (uni  .LT.  0.)  THEN
         uni = uni + 1.0
         carry = twom24
      ELSE
         carry = 0.
      ENDIF
      seeds(i24) = uni
      i24 = i24 - 1
      IF (i24  .EQ.  0)  i24 = 24
      j24 = j24 - 1
      IF (j24  .EQ.  0)  j24 = 24
      rvec(ivec) = uni
  100 CONTINUE
      RETURN
*           entry to input and float INTEGER seeds from previous run
      entry carres(isdext)
         twom24 = 1.
         DO 195 i= 1, 24
  195    twom24 = twom24 * 0.5
      WRITE(6,'(a)') ' full initialization of carran with 25 INTEGERs:'
      WRITE(6,'(5x,5i12)') isdext
      DO 200 i= 1, 24
      seeds(i) = REAL(isdext(i))*twom24
  200 CONTINUE
      carry = REAL(MOD(isdext(25),10))*twom24
      isd = isdext(25)/10
      i24 = MOD(isd,100)
      isd = isd/100
      j24 = isd
      RETURN
*                    entry to ouput seeds as INTEGERs
      entry carout(isdext)
      DO 300 i= 1, 24
         isdext(i) = int(seeds(i)*twop12*twop12)
  300 CONTINUE
      icarry = 0
      IF (carry  .GT.  0.)  icarry = 1
      isdext(25) = 1000*j24 + 10*i24 + icarry
      RETURN
*                    entry to initialize from one INTEGER
      entry carini(inseed)
      jseed = inseed
      WRITE(6,'(a,i12)') ' carran initialized from seed ',inseed
*      twom24 = 1.
         DO 325 i= 1, 24
           twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         IF (jseed  .LT.  0)  jseed = jseed+icons
         iseeds(i) = MOD(jseed,itwo24)
  325    CONTINUE
         DO 350 i= 1,24
         seeds(i) = REAL(iseeds(i))*twom24
  350    CONTINUE
         i24 = 24
         j24 = 10
         carry = 0.
         IF (seeds(24)  .LT.  seeds(14)) carry = twom24
      RETURN
      END

      SUBROUTINE ecuran(rvec,len)
*         random number generator given by l'ecuyer in
*            comm. acm vol 31, p.742, 1988
*            modified by f. james to RETURN a vector of numbers
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*!!!  CALLing sequences for ecuran:                                  ++
*!!!      CALL ecuran (rvec, len)   RETURNs a vector rvec of len     ++
*!!!                   32-bit random floating point numbers between  ++
*!!!                   zero and one.                                 ++
*!!!      CALL ecuini(i1,i2)    initializes the generator from two   ++
*!!!                   32-bit INTEGERs i1 and i2                     ++
*!!!      CALL ecuout(i1,i2)    outputs the current values of the    ++
*!!!                   two INTEGER seeds, to be used for restarting  ++
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION rvec(*)
      SAVE iseed1,iseed2
      DATA iseed1,iseed2 /12345,67890/
*
      DO 100 i= 1, len
      k = iseed1/53668
      iseed1 = 40014*(iseed1 - k*53668) - k*12211
      IF (iseed1  .LT.  0) iseed1=iseed1+2147483563
*
      k = iseed2/52774
      iseed2 = 40692*(iseed2 - k*52774) - k* 3791
      IF (iseed2  .LT.  0) iseed2=iseed2+2147483399
*
      iz = iseed1 - iseed2
      IF (iz  .LT.  1)  iz = iz + 2147483562
*
      rvec(i) = REAL(iz) * 4.656613e-10
  100 CONTINUE
      RETURN
*
      entry ecuini(is1,is2)
      iseed1 = is1
      iseed2 = is2
      RETURN
*
      entry ecuout(is1,is2)
      is1 = iseed1
      is2 = iseed2
      RETURN
      END

      SUBROUTINE varran(drvec,len)
*     ***************************
* switchable random number generator
* translation to DOUBLE PRECISION
*     ***************************
      COMMON / ranpar / keyrnd
      SAVE   / ranpar /
      DOUBLE PRECISION drvec(*)
      DIMENSION rvec(1000)
      IF(len .LT. 1 .OR. len .GT. 1000) GOTO 901
   10 CONTINUE
      IF(keyrnd .EQ. 1) THEN
         CALL marran(rvec,len)
      ELSEIF(keyrnd .EQ. 2) THEN
         CALL ecuran(rvec,len)
      ELSEIF(keyrnd .EQ. 3) THEN
         CALL carran(rvec,len)
      ELSE
         GOTO 902
      ENDIF
* random numbers 0 and 1 not accepted
      DO 30 i=1,len
      IF(rvec(i) .LE. 0e0 .OR. rvec(i) .GE. 1e0) THEN
        WRITE(6,*) ' +++++ varran: rvec=',rvec(i)
        GOTO 10
      ENDIF
      drvec(i)=rvec(i)
   30 CONTINUE
      RETURN
  901 WRITE(6,*) ' +++++ STOP in varran: len=',len
      STOP
  902 WRITE(6,*) ' +++++ STOP in varran: wrong keyrnd',keyrnd
      STOP
      END

      SUBROUTINE bostdq(mode,qq,pp,r)
*     *******************************
* boost along arbitrary axis (by ronald kleiss).
* p boosted into r  from actual frame to rest frame of q
* forth (mode = 1) or back (mode = -1).
* q must be a timelike, p may be arbitrary.
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER ( nout =6 )
      DIMENSION qq(*),pp(*),r(*)
      DIMENSION q(4),p(4)

      DO 10 k=1,4
      p(k)=pp(k)
   10 q(k)=qq(k)
      amq =dsqrt(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      IF    (mode .EQ. -1) THEN
         r(4) = (p(1)*q(1)+p(2)*q(2)+p(3)*q(3)+p(4)*q(4))/amq
         fac  = (r(4)+p(4))/(q(4)+amq)
      ELSEIF(mode .EQ.  1) THEN
         r(4) =(-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)+p(4)*q(4))/amq
         fac  =-(r(4)+p(4))/(q(4)+amq)
      ELSE
         WRITE(nout,*) ' ++++++++ wrong mode in boost3 '
         STOP
      ENDIF
      r(1)=p(1)+fac*q(1)
      r(2)=p(2)+fac*q(2)
      r(3)=p(3)+fac*q(3)
      END


* boost along x axis, exe=exp(eta), eta= hiperbolic velocity.
      SUBROUTINE bostd1(exe,pvec,qvec)
*     ********************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 pvec(4),qvec(4),rvec(4)
      DO 10 i=1,4
  10  rvec(i)=pvec(i)
      rpl=rvec(4)+rvec(1)
      rmi=rvec(4)-rvec(1)
      qpl=rpl*exe
      qmi=rmi/exe
      qvec(2)=rvec(2)
      qvec(3)=rvec(3)
      qvec(1)=(qpl-qmi)/2
      qvec(4)=(qpl+qmi)/2
      END

* boost along z axis, exe=exp(eta), eta= hiperbolic velocity.
      SUBROUTINE bostd3(exe,pvec,qvec)
*     ********************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 pvec(4),qvec(4),rvec(4)
      DO 10 i=1,4
  10  rvec(i)=pvec(i)
      rpl=rvec(4)+rvec(3)
      rmi=rvec(4)-rvec(3)
      qpl=rpl*exe
      qmi=rmi/exe
      qvec(1)=rvec(1)
      qvec(2)=rvec(2)
      qvec(3)=(qpl-qmi)/2
      qvec(4)=(qpl+qmi)/2
      END

      SUBROUTINE rotod1(ph1,pvec,qvec)
*     ********************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 pvec(4),qvec(4),rvec(4)
      phi=ph1
      cs=cos(phi)
      sn=sin(phi)
      DO 10 i=1,4
  10  rvec(i)=pvec(i)
      qvec(1)=rvec(1)
      qvec(2)= cs*rvec(2)-sn*rvec(3)
      qvec(3)= sn*rvec(2)+cs*rvec(3)
      qvec(4)=rvec(4)
      END

      SUBROUTINE rotod2(ph1,pvec,qvec)
*     ********************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 pvec(4),qvec(4),rvec(4)
      phi=ph1
      cs=cos(phi)
      sn=sin(phi)
      DO 10 i=1,4
  10  rvec(i)=pvec(i)
      qvec(1)= cs*rvec(1)+sn*rvec(3)
      qvec(2)=rvec(2)
      qvec(3)=-sn*rvec(1)+cs*rvec(3)
      qvec(4)=rvec(4)
      END

      SUBROUTINE rotod3(ph1,pvec,qvec)
*     ********************************
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 pvec(4),qvec(4),rvec(4)
      phi=ph1
      cs=cos(phi)
      sn=sin(phi)
      DO 10 i=1,4
  10  rvec(i)=pvec(i)
      qvec(1)= cs*rvec(1)-sn*rvec(2)
      qvec(2)= sn*rvec(1)+cs*rvec(2)
      qvec(3)=rvec(3)
      qvec(4)=rvec(4)
      END

      FUNCTION angfi(x,y)
*     *******************
* calculates angle in (0,2*pi) range out of x-y
*     ***********************
      IMPLICIT REAL*8(a-h,o-z)
      DATA pi /3.1415926535897932d0/

      IF(abs(y) .LT. abs(x)) THEN
        the=atan(abs(y/x))
        IF(x .LE. 0d0) the=pi-the
      ELSE
        the=acos(x/sqrt(x**2+y**2))
      ENDIF
      IF(y .LT. 0d0) the=2d0*pi-the
      angfi=the
      END

      SUBROUTINE dumpt(nunit,word,pp)
*     *******************************
      IMPLICIT REAL*8(a-h,o-z)
      CHARACTER*8 word
      REAL*8 pp(4)
      ams=pp(4)**2-pp(3)**2-pp(2)**2-pp(1)**2
      IF(ams .GT. 0.0) ams=sqrt(ams)
      WRITE(nunit,'(1x,a8,5(1x,f13.8))') word,(pp(i),i=1,4),ams
*======================================================================
*================end of yfslib=========================================
*======================================================================
      END
