      FUNCTION fpi2rnew(s)
c CMD-2 parametrization December 2001, gounaris sakurai type??? 
      implicit none
      complex*16 cI,cdels,cepdom,cbeta,cnorm,cfpis,
     &           cBWnu,cBWde,cBWGSr0,cBWGSr1,cBWom
      real *8 fpi2rnew,s,mp2,mp,mr2,mr,mo2,mo,e,xp,rats,beps,ppis,lgs,
     &        hs,pi,gr,mr4,ratm,bepm,ppim,lgm,hm,dhdsm,fs,d,grs,
     &        delta,beta,rBWnu,go,deltabs,deltarg
      real *8 mrho(2),grho(2)
      integer i
      data pi /3.141592653589793d0/
c Constants, Parameters
      mp     =139.56995D0
      mrho(1)=776.09d0          ! +/- 0.64 +/- 0.50    ***************
      grho(1)=144.46d0          ! +/- 1.33 +/- 0.80    ***************
      mrho(2)=1465.0d0          ! +/-  25.0
      grho(2)= 310.0d0          ! +/-  60.0            
      beta   =-.0695d0          ! +/-  .0053           ***************
      deltabs=1.57d-3           ! +/- 0.15d-3 +/- 0.05d-3     modulus of delta **********
      deltarg=12.6d0            ! +/- 3.7d0 +/- 0.2d0          phase angle in degrees *******
      mo     =782.71d0          ! +/- 0.08
      go     =  8.68d0          ! +/- 0.24
      cI   =DCMPLX(0.D0,1.D0)
      mp2  =mp*mp
      mr2  =mr*mr
      mo2  =mo*mo
c Auxiliary variables
      e    =dsqrt(s)
      xp   =4.d0*mp2
      rats =xp/s
      beps =dsqrt(1.d0-rats)
      ppis =0.5d0*e *beps
      lgs  =dlog((e +2.d0*ppis)/(2.d0*mp))
      hs   =2.d0/pi*ppis/e *lgs
      do i=1,2 
        mr   =mrho(i)
        gr   =grho(i)
        mr2  =mr *mr
        mr4  =mr2*mr2
        ratm =xp/mr2
        bepm =dsqrt(1.d0-ratm)
        ppim =0.5d0*mr*bepm
        lgm  =dlog((mr+2.d0*ppim)/(2.d0*mp))
        hm   =2.d0/pi*ppim/mr*lgm
        dhdsm=hm*(0.125d0/ppim**2-0.5d0/mr2)+0.5d0/pi/mr2
        fs   =gr*mr2/ppim**3*(ppis**2*(hs-hm)+(mr2-s)*ppim**2*dhdsm)
        d    =3.d0/pi*mp2/ppim**2*lgm+mr/2.d0/pi/ppim-mp2*mr/pi/ppim**3
        if (i.eq.1) then
          grs=gr*(ppis/ppim)**3*(mr/e)
        else 
          grs=gr
        endif
        rBWnu=mr2*(1.d0+d*gr/mr)
        cBWnu=DCMPLX(rBWnu,0.D0)
        cBWde=DCMPLX(mr2-s+fs,-mr*grs)
        if (i.eq.1) then
          cBWGSr0=cBWnu/cBWde
        else 
          cBWGSr1=cBWnu/cBWde
        endif
      enddo  
c delta argument in radians       
      deltarg=deltarg*2.d0*pi/360.d0
      cdels =DCMPLX(dcos(deltarg),dsin(deltarg))*
     &       DCMPLX(deltabs*s/mo2,0.d0)
      cBWom =DCMPLX(mo2,0.d0)/DCMPLX(mo2-s,-mo*go)
      cepdom=DCMPLX(1.d0,0.0d0)+cdels*cBWom
      cbeta =DCMPLX(beta,0.D0)
      cnorm =DCMPLX(1.d0+beta,0.D0)
      cfpis =(cBWGSr0*cepdom+cbeta*cBWGSr1)/cnorm
      fpi2rnew=cfpis*DCONJG(cfpis)
      fpi2rnew=abs(fpi2rnew)
      END

