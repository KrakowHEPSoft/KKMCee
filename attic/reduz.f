      SUBROUTINE reduz0(qq,p1,p2,pr1,pr2)
*     ***********************************
* reduction of momenta for beta0, second one
* i.e. we mapp:   p1,p2 ==> pr1,pr2
* such that  pr1+pr2 = qq
* resulting pri qri are in qq rest frame.
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( eps1 =1d-15)
      DIMENSION qq(4),p1(4),p2(4),pr1(4),pr2(4)
      DIMENSION pp(4),px1(4),px2(4),ppx(4)
*
      DO 20 k=1,4
 20   pp(k)=p1(k)+p2(k)
      IF((pp(1)**2+pp(2)**2+pp(3)**2)/pp(4)**2  .GT.  eps1) THEN
* transform all momenta to qq rest-frame
         CALL bostdq( 1,qq,p1 ,px1)
         CALL bostdq( 1,qq,p2 ,px2)
         CALL bostdq( 1,qq,pp ,ppx)
* transform all momenta to pp rest-frame
         CALL bostdq( 1,ppx,px1,px1)
         CALL bostdq( 1,ppx,px2,px2)
      ELSE
* DO nothing IF we are alREADy in pp rest-frame
         DO 23 k=1,4
            px1(k)=p1(k)
   23       px2(k)=p2(k)
      ENDIF
* construct reduced beam momenta pr1,pr2
* note: they are understood to be in qq rest-frame
      svar1 = qq(4)**2-qq(3)**2-qq(2)**2-qq(1)**2
      svar  = pp(4)**2-pp(3)**2-pp(2)**2-pp(1)**2
      vv    = 1d0 -svar1/svar
      IF(abs(vv) .GT.  eps1) THEN
         amel2=  p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2
         pr1(4)= sqrt(svar1)/2d0
         pr2(4)= pr1(4)
         pxmod = sqrt(px1(1)**2+px1(2)**2+px1(3)**2)
         prmod = sqrt(pr1(4)**2-amel2)
         DO 30 k=1,3
         pr1(k)= px1(k)/pxmod*prmod
 30      pr2(k)= px2(k)/pxmod*prmod
      ELSE
         DO 40 k=1,4
         pr1(k)= px1(k)
 40      pr2(k)= px2(k)
      ENDIF
      END

      SUBROUTINE reduz1(qq,p1,p2,ph,pr1,pr2,phr)
*     ******************************************
* reduction of 4-momenta for beta1
*           p1,p2,ph ==--> pr1,pr2,phr
* such that  pr1+pr2 = qq+phr
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( eps1 =1d-15)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
      DIMENSION qq(4), p1(4), p2(4), ph(4), pr1(4),pr2(4),phr(4)
      DIMENSION pp(4),qqk(4),ppx(4), ppk(4)
      DIMENSION px1(4),px2(4),phx(4)
*
      DO 20 k=1,4
      pp(k)   = p1(k)+p2(k)
      ppk(k)  = p1(k)+p2(k)-ph(k)
 20   qqk(k)  = qq(k)+ph(k)
      svar  =  pp(4)**2 -pp(3)**2 -pp(2)**2 -pp(1)**2
      svar1 =  qq(4)**2 -qq(3)**2 -qq(2)**2 -qq(1)**2
      ss1   = ppk(4)**2-ppk(3)**2-ppk(2)**2-ppk(1)**2
      ss2   = qqk(4)**2-qqk(3)**2-qqk(2)**2-qqk(1)**2
      IF((pp(1)**2+pp(2)**2+pp(3)**2)/pp(4)**2  .GT.  eps1) THEN
* transform all momenta to qq rest-frame
         CALL bostdq( 1,qq,p1 ,px1)
         CALL bostdq( 1,qq,p2 ,px2)
         CALL bostdq( 1,qq,ph ,phx)
         CALL bostdq( 1,qq,pp ,ppx)
* transform all momenta to pp rest-frame
         CALL bostdq( 1,ppx,px1,px1)
         CALL bostdq( 1,ppx,px2,px2)
         CALL bostdq( 1,ppx,phx,phx)
      ELSE
* DO nothing IF we are alREADy in pp rest-frame
         DO 23 k=1,4
            phx(k)=ph(k)
            px1(k)=p1(k)
   23       px2(k)=p2(k)
      ENDIF
* construct reduced beam momenta pr1,pr2
* note: they are understood to be in qq rest-frame
      vv2   = 1d0 - ss2/svar
      IF(abs(vv2) .GT.  eps1) THEN
         pk    =  (px1(4)+px2(4))*phx(4)
*cccc    xlam= sqrt(svar1/svar+(pk/svar)**2)+pk/svar
         xlam= sqrt(svar1/ss1)
         amel2=  p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2
         pxmod = sqrt(px1(1)**2+px1(2)**2+px1(3)**2)
         px1(4)= px1(4)*xlam
         px2(4)= px2(4)*xlam
*cc      prmod = sqrt(px1(4)**2-amel2)
         prmod =      px1(4)**2-amel2
         IF(prmod .LE. 0d0) WRITE(nout,*) ' reduz1: prmod=', prmod
         IF(prmod .LE. 0d0) WRITE(   6,*) ' reduz1: prmod=', prmod
         prmod = sqrt(abs(prmod))
         DO 30 k=1,3
         px1(k)= px1(k)/pxmod*prmod
 30      px2(k)= px2(k)/pxmod*prmod
         DO 31 k=1,4
 31      phx(k)= phx(k)*xlam
      ENDIF
* THEN, boost away the three-vector part of p1+p2-ph
* that is transform to qq rest frame
      DO 35 k=1,4
 35   pp(k)= px1(k)+px2(k)-phx(k)
      CALL bostdq( 1,pp,px1,pr1)
      CALL bostdq( 1,pp,px2,pr2)
      CALL bostdq( 1,pp,phx,phr)
      END

      SUBROUTINE reduz2(qq,p1,p2,ph1,ph2,pr1,pr2,ph1r,ph2r)
*     *****************************************************
* reduction for beta2
*           p1,p2,ph1,ph2 ==--> pr1,pr2,ph1r,ph2r
* such that  pr1+pr2 = ph1r+ph2r+qq
* input:  qq,p1,p2,ph1,ph2
* output: pr1,pr2,ph1r,ph2r
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( eps1 =1d-15)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
      DIMENSION qq(*), p1(*),  p2(*),  ph1(*),  ph2(*)
      DIMENSION        pr1(*), pr2(*), ph1r(*), ph2r(*)
      DIMENSION pp(4),qqk(4),ppk(4),ppx(4)
      DIMENSION px1(4),px2(4),ph1x(4),ph2x(4),sph(4)
*
      DO 20 k=1,4
      pp(k)   = p1(k)+p2(k)
      ppk(k)  = p1(k)+p2(k)-ph1(k)-ph2(k)
 20   qqk(k)  = qq(k)+ph1(k)+ph2(k)
      svar  =  pp(4)**2 -pp(3)**2 -pp(2)**2 -pp(1)**2
      svar1 =  qq(4)**2 -qq(3)**2 -qq(2)**2 -qq(1)**2
      ss1   = ppk(4)**2-ppk(3)**2-ppk(2)**2-ppk(1)**2
      ss2   = qqk(4)**2-qqk(3)**2-qqk(2)**2-qqk(1)**2
      IF((pp(1)**2+pp(2)**2+pp(3)**2)/pp(4)**2  .GT.  eps1) THEN
* transform all momenta to qq rest-frame
         CALL bostdq( 1,qq,p1 ,px1)
         CALL bostdq( 1,qq,p2 ,px2)
         CALL bostdq( 1,qq,ph1,ph1x)
         CALL bostdq( 1,qq,ph2,ph2x)
         CALL bostdq( 1,qq,pp ,ppx)
* transform all momenta to pp rest-frame
         CALL bostdq( 1,ppx,px1,px1)
         CALL bostdq( 1,ppx,px2,px2)
         CALL bostdq( 1,ppx,ph1x,ph1x)
         CALL bostdq( 1,ppx,ph2x,ph2x)
      ELSE
* DO nothing IF we are alREADy in pp rest-frame
         DO 23 k=1,4
            ph1x(k)=ph1(k)
            ph2x(k)=ph2(k)
            px1(k)=p1(k)
   23       px2(k)=p2(k)
      ENDIF
* construct reduced beam momenta pr1,pr2
* note: they are understood to be in qq rest-frame
      vv2   = 1d0 - ss2/svar
      IF(abs(vv2) .GT. 1d-6) THEN
* construct reduced beam momenta pr1,pr2
* start with dilatation of beams
         DO 24 k=1,4
         pp(k)  =  px1(k)+px2(k)
  24     sph(k) =  ph1(k)+ph2(k)
         pk     =  pp(4)*sph(4)
         sk2    =  sph(4)**2 -sph(3)**2 -sph(2)**2 -sph(1)**2
*ccc     xlam   =  sqrt((svar1-sk2)/svar+(pk/svar)**2)+pk/svar
         xlam   =  sqrt(svar1/ss1)
         amel2  =  p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2
         pxmod  =  sqrt(px1(1)**2+px1(2)**2+px1(3)**2)
         px1(4) =  px1(4)*xlam
         px2(4) =  px2(4)*xlam
*ccc     prmod  =  sqrt(px1(4)**2-amel2)
         prmod  =      px1(4)**2-amel2
         IF(prmod .LE. 0d0) WRITE(nout,*) ' reduz2: prmod=', prmod
         IF(prmod .LE. 0d0) WRITE(   6,*) ' reduz2: prmod=', prmod
         prmod  = sqrt(abs(prmod))
         DO 30 k=1,3
         px1(k) = px1(k)/pxmod*prmod
 30      px2(k) = px2(k)/pxmod*prmod
         DO 31 k=1,4
         ph1x(k)= ph1x(k)*xlam
 31      ph2x(k)= ph2x(k)*xlam
      ENDIF
* then, boost away the three-vector part of p1+p2-ph1-ph2
* that is transform to qq rest frame
      DO 35 k=1,4
 35   pp(k)= px1(k)+px2(k)-ph1x(k)-ph2x(k)
      CALL bostdq( 1,pp,px1,pr1)
      CALL bostdq( 1,pp,px2,pr2)
      CALL bostdq( 1,pp,ph1x,ph1r)
      CALL bostdq( 1,pp,ph2x,ph2r)
      END

      SUBROUTINE tralqq(mode,q,p,r)
*     *****************************
* boost along z axis to a frame where qq(3)=0
* and next along transverse direction of qq,
* forth (mode = 1) or back (mode = -1).
* q must be a timelike, p may be arbitrary.
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION q(*),p(*),r(*)
      DIMENSION ql(4),qt(4)
      ql(4)=q(4)
      ql(3)=q(3)
      ql(2)=0d0
      ql(1)=0d0
      CALL bostdq(1,ql,q,qt)
      IF(mode .EQ. 1) THEN
        CALL bostdq( 1,ql,p,r)
        CALL bostdq( 1,qt,r,r)
      ELSE
        CALL bostdq(-1,qt,p,r)
        CALL bostdq(-1,ql,r,r)
      ENDIF
      END
