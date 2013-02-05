


*=======================================================================
*=======================================================================
*    !!!!!! TRASHED !!!!!!
*=======================================================================
*=======================================================================



      SUBROUTINE dumpbt(nout)
*     ***********************
*     prints out information on beta's
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(npmx=100)
* spying on weights.........
      COMMON / betspy / beta00,beta01,beta02,beti01,beti02,
     $  beta10(npmx),beta11(npmx),beta20(npmx,npmx),
     $  sfacx(npmx),sfacy(npmx),
     $  betx10(npmx),betx11(npmx),bety10(npmx),bety11(npmx),
     $  betxx20(npmx,npmx),betxy20(npmx,npmx),betyy20(npmx,npmx),
     $  beti10(npmx),beti11(npmx),beti20(npmx,npmx),
     $  betf10(npmx),betf11(npmx),betf20(npmx,npmx)
      SAVE / betspy /
* ...........
      COMMON / momini / xf1(4),xf2(4),xphum(4),xphot(100,4),nphox
      COMMON / momfin / qf1(4),qf2(4),yphum(4),yphot(100,4),nphoy
      SAVE / momini /,/ momfin /

      WRITE(nout,*) '--------------------<dumbt>----------------------'
 2000 format(a30,3f15.8/(5f15.8))
* ...beta0
      WRITE(nout,2000) ' beta00, beta01/beta00, beta02/beta00'
     & ,beta00,beta01/beta00,beta02/beta00
* ...sfaci(i)
 2002 format(a16,4e16.7/(5e16.7))
      IF(nphox .GT. 0)WRITE(nout,2002) ' sfacx(i) ',(sfacx(i),i=1,nphox)
      IF(nphoy .GT. 0)WRITE(nout,2002) ' sfacy(i) ',(sfacy(i),i=1,nphoy)
* ...beta1   initial state ...............
      IF(nphox .GT. 0) WRITE(nout,2000) ' init. betx10(i)/beta00 ',
     $  ( betx10(i)/sfacx(i)/beta00  ,i=1,nphox)
      IF(nphox .GT. 0) WRITE(nout,2000) ' init. betx11(i)/beta00 ',
     $  ( betx11(i)/sfacx(i)/beta00  ,i=1,nphox)

* ...beta1   final state ...............
      IF(nphoy .GT. 0) WRITE(nout,2000) ' fin.  bety10(i)/beta00 ',
     $  ( bety10(i)/sfacy(i)/beta00  ,i=1,nphoy)
      IF(nphoy .GT. 0) WRITE(nout,2000) ' fin.  bety11(i)/beta00 ',
     $  ( bety11(i)/sfacy(i)/beta00  ,i=1,nphoy)

* ...beta2 pure initial
      IF(nphox .GE. 2) THEN
      WRITE(nout,*) ' init. state betxx20(i,j)'
      DO 100 j=2,nphox
      WRITE(nout,'(  8f12.8)')
     $ (betxx20(i,j)/sfacx(i)/sfacx(j)/beta00,i=1,j-1)
  100 CONTINUE
      ENDIF
* ...beta2 pure final
      IF(nphoy .GE. 2) THEN
      WRITE(nout,*) ' final state betyy20(i,j)'
      DO 110 j=2,nphoy
      WRITE(nout,'(  8f20.16)')
     $ (betyy20(i,j)/sfacy(i)/sfacy(j)/beta00,i=1,j-1)
  110 CONTINUE
      ENDIF
* ...beta2 initial/final
      IF(nphox .GE. 1  .AND.  nphoy .GE. 1) THEN
      WRITE(nout,*) ' in/fi state betxy20(i,j)'
      DO 120 j=1,nphoy
      WRITE(nout,'(  8f12.8)')
     $ (betxy20(i,j)/sfacx(i)/sfacy(j)/beta00,i=1,nphox)
  120 CONTINUE
      ENDIF
      END





      SUBROUTINE D_ini2(amel,y1,z1,y2,z2,pdp,gi1,gi2)
*     ***********************************************
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
* inline functions
      wm (a  )=     (1d0-a)**2
      wms(a,b)=     ((1d0-a)**2+(1d0-b)**2)
      wwm(am,a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
*
      ame2  = amel**2/(2d0*pdp)
      sfac1  =  2d0/(pdp*y1*z1)*wwm(ame2,y1,z1)
      sfac2  =  2d0/(pdp*y2*z2)*wwm(ame2,y2,z2)
      y1p= y1/(1d0-y2)
      z1p= z1/(1d0-z2)
      y2p= y2/(1d0-y1)
      z2p= z2/(1d0-z1)
      IF((y1+z1) .GT. (y2+z2)) THEN
        x1=wm (y1   )*wms(y2p,z2p) +wm (y1p    )*wms(y2,z2)
        x2=wm (   z1)*wms(y2p,z2p) +wm (    z1p)*wms(y2,z2)
      ELSE
        x1=wm (y2   )*wms(y1p,z1p) +wm (y2p    )*wms(y1,z1)
        x2=wm (   z2)*wms(y1p,z1p) +wm (    z2p)*wms(y1,z1)
      ENDIF
      gi1 = x1*sfac1*sfac2/8d0
      gi2 = x2*sfac1*sfac2/8d0
      END

      SUBROUTINE D_fin2(amfin,yy1,zz1,yy2,zz2,qdq,gf1,gf2)
*     ***********************************************
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
* inline functions
      wm (a  )=     (1d0-a)**2
      wms(a,b)=     ((1d0-a)**2+(1d0-b)**2)
      wwm(am,a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
*
      amf2 = amfin**2/(2d0*qdq)
      yy1p= yy1/(1d0+yy2)
      zz1p= zz1/(1d0+zz2)
      yy2p= yy2/(1d0+yy1)
      zz2p= zz2/(1d0+zz1)
      y1  = yy1/(1+yy1+zz1)
      y2  = yy2/(1+yy2+zz2)
      z1  = zz1/(1+yy1+zz1)
      z2  = zz2/(1+yy2+zz2)
      y1p = yy1p/(1+yy1p+zz1p)
      y2p = yy2p/(1+yy2p+zz2p)
      z1p = zz1p/(1+yy1p+zz1p)
      z2p = zz2p/(1+yy2p+zz2p)
      sfac1  =  2d0/(qdq*yy1*zz1)*wwm(amf2,y1,z1)
      sfac2  =  2d0/(qdq*yy2*zz2)*wwm(amf2,y2,z2)
      IF((y1+z1) .GT. (y2+z2)) THEN
* argument of wms should be yy2,zz2 etc.. ????
* not realy, yy,zz<1 violated!!!
        x1=wm (y1   )*wms(y2p,z2p) +wm (y1p    )*wms(y2,z2)
        x2=wm (   z1)*wms(y2p,z2p) +wm (    z1p)*wms(y2,z2)
      ELSE
        x1=wm (y2   )*wms(y1p,z1p) +wm (y2p    )*wms(y1,z1)
        x2=wm (   z2)*wms(y1p,z1p) +wm (    z2p)*wms(y1,z1)
      ENDIF
      gf1 = x1*sfac1*sfac2/8d0
      gf2 = x2*sfac1*sfac2/8d0
      END

      SUBROUTINE D_infi2(amel,amfin,y1,z1,yy2,zz2,
     $     pdp,qdq,gi1,gi2,gf1,gf2)
*     ***********************************************
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
* inline functions
      wwm(am,a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)

      sfac1  =  2d0/(pdp*y1*z1)
      ame2 = amel**2/(2d0*pdp)
      gi1   = 0.5d0*(1d0-y1)**2 *wwm(ame2,y1,z1) *sfac1
      gi2   = 0.5d0*(1d0-z1)**2 *wwm(ame2,y1,z1) *sfac1

      y2  = yy2/(1 +yy2+zz2)
      z2  = zz2/(1 +yy2+zz2)
      sfac2  =  2d0/(qdq*yy2*zz2)
      amf2 = amfin**2/(2d0*qdq)
      gf1   = 0.5d0*(1d0-y2)**2 *wwm(amf2,y2,z2) *sfac2
      gf2   = 0.5d0*(1d0-z2)**2 *wwm(amf2,y2,z2) *sfac2
      END


      SUBROUTINE D_bet0(xx,p1,p2,q1,q2,deli1,deli2,delf1,delf2)
*     ***************************************************************
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / keyyfs / keyzet,keybrm,keyfix,keyred,keywgt
      SAVE   / keyyfs /
      DIMENSION xx(*),p1(*),p2(*),q1(*),q2(*)
      CALL bvirt0(p1,p2,deli1,deli2)
      CALL bvirt0(q1,q2,delf1,delf2)
* ...initial/final state bremsstrahlung switches
      keybin  = MOD(keybrm,10)
      keybfi  = MOD(keybrm,100)/10
      deli1   = deli1*keybin
      deli2   = deli2*keybin
      delf1   = delf1*keybfi
      delf2   = delf2*keybfi
      END




      SUBROUTINE sfach0(p1,p2,ph,sfac0)
*     *********************************
* calculates soft factor for REAL soft photon.
*     *********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph(4),p1(4),p2(4)
      am2= p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      pk1= p1(4)*ph(4)-p1(1)*ph(1)-p1(2)*ph(2)-p1(3)*ph(3)
      pk2= p2(4)*ph(4)-p2(1)*ph(1)-p2(2)*ph(2)-p2(3)*ph(3)
      pp = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      sfac0  =  2d0*pp/pk1/pk2 -am2/pk1**2 -am2/pk2**2
      END


      SUBROUTINE sof1ini(p1,p2,ph,f1,f2)
*     *********************************
* calculates ingredients for REAL single photon diff. xsection
*     *****************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph(*),p1(*),p2(*)
      am2= p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      pk1= p1(4)*ph(4)-p1(1)*ph(1)-p1(2)*ph(2)-p1(3)*ph(3)
      pk2= p2(4)*ph(4)-p2(1)*ph(1)-p2(2)*ph(2)-p2(3)*ph(3)
      pp = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      b  =  pk2/pp
      a  =  pk1/pp
      sfac1  =  2d0*pp/pk1/pk2
      am = am2/(2d0*pp)
      wwm= 1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
      f1   = 0.5d0*(1d0-a)**2*wwm *sfac1
      f2   = 0.5d0*(1d0-b)**2*wwm *sfac1
      END

      SUBROUTINE sof1fin(p1,p2,ph,f1,f2)
*     *********************************
* final state now! but p <=> replacement kept
* calculates ingredients for REAL single photon dIFf. xsection
*     *****************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph(*),p1(*),p2(*)
      am2= p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      pk1= p1(4)*ph(4)-p1(1)*ph(1)-p1(2)*ph(2)-p1(3)*ph(3)
      pk2= p2(4)*ph(4)-p2(1)*ph(1)-p2(2)*ph(2)-p2(3)*ph(3)
      pp = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
*
      bb =  abs(pk2)/pp
      aa =  abs(pk1)/pp
      a  = aa/(1 +aa+bb)
      b  = bb/(1 +aa+bb)
      sfac1  =  2d0*pp/pk1/pk2
      am = am2/(2d0*pp)
      wwm= 1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
      f1   = 0.5d0*(1d0-a)**2*wwm *sfac1
      f2   = 0.5d0*(1d0-b)**2*wwm *sfac1
      END


      SUBROUTINE sof2ini(p1,p2,ph1,ph2,f1,f2)
*     **************************************
* new version by ela was
* calculates ingredients for REAL DOUBLE photon dIFf. xsection
*     *****************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph1(*),ph2(*),p1(*),p2(*)
*
      wm (a  )=     (1d0-a)**2
      wms(a,b)=     ((1d0-a)**2+(1d0-b)**2)
      wwm(a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
*
      pp = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      am2= p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      am = am2/(2d0*pp)
      b1 = (p2(4)*ph1(4)-p2(1)*ph1(1)-p2(2)*ph1(2)-p2(3)*ph1(3))/pp
      a1 = (p1(4)*ph1(4)-p1(1)*ph1(1)-p1(2)*ph1(2)-p1(3)*ph1(3))/pp
      b2 = (p2(4)*ph2(4)-p2(1)*ph2(1)-p2(2)*ph2(2)-p2(3)*ph2(3))/pp
      a2 = (p1(4)*ph2(4)-p1(1)*ph2(1)-p1(2)*ph2(2)-p1(3)*ph2(3))/pp
      sfac1  =  2d0/(pp*a1*b1)*wwm(a1,b1)
      sfac2  =  2d0/(pp*a2*b2)*wwm(a2,b2)
      a1p= a1/(1d0-a2)
      b1p= b1/(1d0-b2)
      a2p= a2/(1d0-a1)
      b2p= b2/(1d0-b1)
      IF((a1+b1) .GT. (a2+b2)) THEN
        x1=wm (a1   )*wms(a2p,b2p) +wm (a1p    )*wms(a2,b2)
        x2=wm (   b1)*wms(a2p,b2p) +wm (    b1p)*wms(a2,b2)
      ELSE
        x1=wm (a2   )*wms(a1p,b1p) +wm (a2p    )*wms(a1,b1)
        x2=wm (   b2)*wms(a1p,b1p) +wm (    b2p)*wms(a1,b1)
      ENDIF
      f1 = x1*sfac1*sfac2/8d0
      f2 = x2*sfac1*sfac2/8d0
*.. correction ela was november 1989................................
*.. this correction reconstructs properly DOUBLE collinear limit
*.. and affects below photon-fermion angle  <0.1 amel/ene
c      sfac1  =  2d0/(pp*a1*b1)
c      sfac2  =  2d0/(pp*a2*b2)
c      wwm1=1d0-wwm(a1,b1)
c      wwm2=1d0-wwm(a2,b2)
c      delt=(b2**2*a1**2+a2**2*b1**2)/(x1+x2)*2d0*
c     #  ( 1d0/(a1+a2)**2+1d0/(b1+b2)**2)
c      wminf=1d0-wwm1-wwm2+wwm1*wwm2*(1d0+delt)
c      f1 = x1*sfac1*sfac2/8d0*wminf
c      f2 = x2*sfac1*sfac2/8d0*wminf
*...END of correction............................................
      END

      SUBROUTINE sof2fin(q1,q2,ph1,ph2,f1,f2)
*     **************************************
* calculates ingredients for REAL DOUBLE photon dIFf. xsection
*     *****************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph1(*),ph2(*),q1(*),q2(*)
*
      wm (a  )=     (1d0-a)**2
      wms(a,b)=     ((1d0-a)**2+(1d0-b)**2)
      wwm(a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
*
      pp = q1(4)*q2(4)-q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)
      am2= q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      am = am2/(2d0*pp)
      bb1=abs(q2(4)*ph1(4)-q2(1)*ph1(1)-q2(2)*ph1(2)-q2(3)*ph1(3))/pp
      aa1=abs(q1(4)*ph1(4)-q1(1)*ph1(1)-q1(2)*ph1(2)-q1(3)*ph1(3))/pp
      bb2=abs(q2(4)*ph2(4)-q2(1)*ph2(1)-q2(2)*ph2(2)-q2(3)*ph2(3))/pp
      aa2=abs(q1(4)*ph2(4)-q1(1)*ph2(1)-q1(2)*ph2(2)-q1(3)*ph2(3))/pp
      aa1p= aa1/(1d0+aa2)
      bb1p= bb1/(1d0+bb2)
      aa2p= aa2/(1d0+aa1)
      bb2p= bb2/(1d0+bb1)
      a1  = aa1/(1+aa1+bb1)
      a2  = aa2/(1+aa2+bb2)
      b1  = bb1/(1+aa1+bb1)
      b2  = bb2/(1+aa2+bb2)
      a1p = aa1p/(1+aa1p+bb1p)
      a2p = aa2p/(1+aa2p+bb2p)
      b1p = bb1p/(1+aa1p+bb1p)
      b2p = bb2p/(1+aa2p+bb2p)
      sfac1  =  2d0/(pp*aa1*bb1)*wwm(a1,b1)
      sfac2  =  2d0/(pp*aa2*bb2)*wwm(a2,b2)
      IF((a1+b1) .GT. (a2+b2)) THEN
        x1=wm (a1   )*wms(a2p,b2p) +wm (a1p    )*wms(a2,b2)
        x2=wm (   b1)*wms(a2p,b2p) +wm (    b1p)*wms(a2,b2)
      ELSE
        x1=wm (a2   )*wms(a1p,b1p) +wm (a2p    )*wms(a1,b1)
        x2=wm (   b2)*wms(a1p,b1p) +wm (    b2p)*wms(a1,b1)
      ENDIF
      f1 = x1*sfac1*sfac2/8d0
      f2 = x2*sfac1*sfac2/8d0
      END


      SUBROUTINE sof2inix(p1,p2,ph1,ph2,f1,f2)
*     **************************************
* old version OBSOLETE
* calculates ingredients for REAL DOUBLE photon dIFf. xsection
*     *****************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION ph1(*),ph2(*),p1(*),p2(*)
*
      wm (a  )=     (1d0-a)**2
      wms(a,b)=     ((1d0-a)**2+(1d0-b)**2)
      wwm(a,b)=
     $   1d0-am*2d0*(1d0-a)*(1d0-b)/((1d0-a)**2+(1d0-b)**2)*(a/b+b/a)
*
      pp = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      am2= p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      am = am2/(2d0*pp)
      b1 = (p2(4)*ph1(4)-p2(1)*ph1(1)-p2(2)*ph1(2)-p2(3)*ph1(3))/pp
      a1 = (p1(4)*ph1(4)-p1(1)*ph1(1)-p1(2)*ph1(2)-p1(3)*ph1(3))/pp
      b2 = (p2(4)*ph2(4)-p2(1)*ph2(1)-p2(2)*ph2(2)-p2(3)*ph2(3))/pp
      a2 = (p1(4)*ph2(4)-p1(1)*ph2(1)-p1(2)*ph2(2)-p1(3)*ph2(3))/pp
      sfac1  =  2d0/(pp*a1*b1)*wwm(a1,b1)
      sfac2  =  2d0/(pp*a2*b2)*wwm(a2,b2)
      a1p= a1/(1d0-a2)
      b1p= b1/(1d0-b2)
      a2p= a2/(1d0-a1)
      b2p= b2/(1d0-b1)
      IF((a1+b1) .GT. (a2+b2)) THEN
        x1=wm (a1   )*wms(a2p,b2p) +wm (a1p    )*wms(a2,b2)
        x2=wm (   b1)*wms(a2p,b2p) +wm (    b1p)*wms(a2,b2)
      ELSE
        x1=wm (a2   )*wms(a1p,b1p) +wm (a2p    )*wms(a1,b1)
        x2=wm (   b2)*wms(a1p,b1p) +wm (    b2p)*wms(a1,b1)
      ENDIF
      f1 = x1*sfac1*sfac2/8d0
      f2 = x2*sfac1*sfac2/8d0
      END


      SUBROUTINE ntest0(qq,p1,p2,q1,q2)
*     *********************************
* ...testing redustion for beta0
*     *********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
      DIMENSION qq(*),p1(*),p2(*),q1(*),q2(*)
      DIMENSION pr1(4),pr2(4),qr1(4),qr2(4)
      DIMENSION ppr(4),qqr(4),qq1(4)
*
      WRITE(nout,*) '-----------ntest0---------==========='
      CALL reduz0(qq,q1,q2,qr1,qr2)
      CALL reduz0(qq,p1,p2,pr1,pr2)
      CALL bostdq(1,qq,qq,qq1)
      DO 10 k=1,4
      ppr(k)=pr1(k)+pr2(k)-qq1(k)
  10  qqr(k)=qr1(k)+qr2(k)-qq1(k)
      CALL dumpt(2,' ppr    ',ppr)
      CALL dumpt(2,' qqr    ',qqr)
      CALL dumpt(2,'  pr1   ',pr1)
      CALL dumpt(2,'  pr2   ',pr2)
      CALL dumpt(2,'  qr1   ',qr1)
      CALL dumpt(2,'  qr2   ',qr2)
      CALL gthet0(pr1,qr1,costh )
      WRITE(nout,'(a,3f20.12)') ' costh= ',costh
      END

      SUBROUTINE ntest1(qq,p1,p2,q1,q2,ph)
*     ************************************
* ...testing reduction for beta1
*     ************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
      DIMENSION qq(*),p1(*),p2(*),q1(*),q2(*),ph(*)
      DIMENSION pr1(4),pr2(4),qr1(4),qr2(4),phr(4)
      DIMENSION ppr(4),qqr(4),px1(4),px2(4),qq1(4),ppx(4)
*
      WRITE(nout,*) '-----------ntest1---------<<<<<<<<<<<'
      CALL reduz0(qq,q1,q2,qr1,qr2)
      CALL reduz1(qq,p1,p2,ph,pr1,pr2,phr)
      CALL reduz0(qq,p1,p2,px1,px2)
      CALL bostdq(1,qq,qq,qq1)
      DO 10 k=1,4
      ppx(k)=px1(k)+px2(k)        -qq1(k)
      ppr(k)=pr1(k)+pr2(k)-phr(k) -qq1(k)
  10  qqr(k)=qr1(k)+qr2(k)        -qq1(k)
      CALL dumpt(2,' ppx    ',ppx)
      CALL dumpt(2,' ppr    ',ppr)
      CALL dumpt(2,' qqr    ',qqr)
      CALL dumpt(2,'  pr1   ',pr1)
      CALL dumpt(2,'  pr2   ',pr2)
      CALL dumpt(2,'  px1   ',px1)
      CALL dumpt(2,'  px2   ',px2)
      CALL dumpt(2,'  qr1   ',qr1)
      CALL dumpt(2,'  qr2   ',qr2)
* single bremsstrahlung xsection
      CALL sof1ini(p1,p2,ph,gf1,gf2)
      CALL gthet1(pr1,pr2,qr1,costh1,costh2)
      CALL sfach0(p1,p2, ph ,sfacj)
      CALL gthet0(px1,qr1,costh)
      WRITE(nout,'(a,3f20.12)') ' ph(4)= ',ph(4)
      WRITE(nout,'(a,3f20.12)') ' costh= ',costh1-costh,costh2-costh
      WRITE(nout,'(a,3f20.12)') ' gf/sf= ',(gf1+gf2)/sfacj
      END

      SUBROUTINE ntest2(qq,p1,p2,q1,q2,ph1,ph2)
*     *****************************************
* ...testing reduction for beta2
*     ************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / inout  / ninp,nout
      SAVE   / inout  /
      DIMENSION qq(*),p1(*),p2(*),q1(*),q2(*),ph1(*),ph2(*)
      DIMENSION pr1(4),pr2(4),qr1(4),qr2(4),ph1r(4),ph2r(4)
      DIMENSION ppr(4),qqr(4),px1(4),px2(4),qq1(4)
*
      WRITE(nout,*) '-----------ntest2---------<<<<<<<<<<<'
      CALL reduz0(qq,q1,q2,qr1,qr2)
      CALL reduz2(qq,p1,p2,ph1,ph2,pr1,pr2,ph1r,ph2r)
      CALL bostdq(1,qq,qq,qq1)
      DO 10 k=1,4
      ppr(k)=pr1(k)+pr2(k)-ph1r(k)-ph2r(k) -qq1(k)
  10  qqr(k)=qr1(k)+qr2(k)                 -qq1(k)
      CALL dumpt(2,' ppr    ',ppr)
      CALL dumpt(2,' qqr    ',qqr)
      CALL dumpt(2,'  pr1   ',pr1)
      CALL dumpt(2,'  pr2   ',pr2)
      CALL dumpt(2,'  qr1   ',qr1)
      CALL dumpt(2,'  qr2   ',qr2)
* single bremsstrahlung xsection
      CALL sof2ini(p1,p2,ph1,ph2,gf1,gf2)
      CALL gthet1(pr1,pr2,qr1,costh1,costh2)
      CALL sfach0(p1,p2, ph1,sfac1)
      CALL sfach0(p1,p2, ph2,sfac2)
      CALL reduz0(qq,p1,p2,px1,px2)
      CALL gthet0(px1,qr1,costh)
      WRITE(nout,'(a,3f20.12)') 'phi(4)= ',ph1(4),ph2(4)
      WRITE(nout,'(a,3f20.12)') ' costh= ',costh1-costh,costh2-costh
      WRITE(nout,'(a,3f20.12)') ' gf/sf= ',(gf1+gf2)/sfac1/sfac2
      END


      SUBROUTINE gthet0(p1,q1,costh)
*     ******************************
* calculates costh between beam and final fermion
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION p1(*),q1(*)
      costh  = (p1(1)*q1(1) +p1(2)*q1(2) +p1(3)*q1(3))
     $            /sqrt((q1(1)**2 +q1(2)**2 +q1(3)**2)
     $                 *(p1(1)**2 +p1(2)**2 +p1(3)**2))
      END



      SUBROUTINE gthet1(p1,p2,q1,costh1,costh2)
*     *****************************************
* calculates costh1 and costh2 between beam amd final
* fermion momenta in final fermion rest frame q1(4)+q2(4)=0
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION p1(*),p2(*),q1(*)
      costh1 = (p1(1)*q1(1) +p1(2)*q1(2) +p1(3)*q1(3))
     $            /sqrt((q1(1)**2 +q1(2)**2 +q1(3)**2)
     $                 *(p1(1)**2 +p1(2)**2 +p1(3)**2))
      costh2 =-(p2(1)*q1(1) +p2(2)*q1(2) +p2(3)*q1(3))
     $            /sqrt((q1(1)**2 +q1(2)**2 +q1(3)**2)
     $                 *(p2(1)**2 +p2(2)**2 +p2(3)**2))
      END



      SUBROUTINE gthet3(p1,p2,q1,q2,cth11,cth12,cth21,cth22)
*     ***************************************************
* calculates costh1 and costh2 between beam amd final
* fermion momenta in Z resonance rest frame q1(4)+q2(4)=0
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION p1(*),p2(*),q1(*),q2(*)
      q1d=        sqrt(q1(1)**2 +q1(2)**2 +q1(3)**2)
      q2d=        sqrt(q2(1)**2 +q2(2)**2 +q2(3)**2)
      p1d=        sqrt(p1(1)**2 +p1(2)**2 +p1(3)**2)
      p2d=        sqrt(p2(1)**2 +p2(2)**2 +p2(3)**2)
      cth11 = (q1(1)*p1(1) +q1(2)*p1(2) +q1(3)*p1(3))/q1d/p1d
      cth12 =-(q1(1)*p2(1) +q1(2)*p2(2) +q1(3)*p2(3))/q1d/p2d
      cth21 =-(q2(1)*p1(1) +q2(2)*p1(2) +q2(3)*p1(3))/q2d/p1d
      cth22 = (q2(1)*p2(1) +q2(2)*p2(2) +q2(3)*p2(3))/q2d/p2d
      END
