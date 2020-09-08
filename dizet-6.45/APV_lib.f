      SUBROUTINE rokap0
     &    (tHmu2,rzm2,rtm2,rhm2,rwm2,DelRho,DelKap,DRh0Le,DRhoLe,DKapLe)
*-
      IMPLICIT NONE
*-      
      REAL*8 tHmu2,rzm2,rtm2,rhm2,rwm2
      REAL*8 ctw2,stw2,ctw4,rhw,rhz,rhw2,rhz2
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm
      REAL*8 DRfer0,DRfer,DRbos,DelRho,DelKap,DRh0Le,DRhoLe,DKapLe
      REAL*8 DREAL,DLOG
      COMPLEX*16 wm2,zm2,cz2,hm2
      COMPLEX*16 b0f,b0fwwz,b0fzww,b0fwwh,b0fzzh,b0fww0
      COMPLEX*16 SZZ_f,SWW_f,PZg_f
*- 
      wm2=DCMPLX(rwm2,-1d-20)
      zm2=DCMPLX(rzm2,-1d-20)
      hm2=DCMPLX(rhm2,-1d-20)
      cz2=DCMPLX(0d0 ,-1d-20)
*- 
      ctw2=rwm2/rzm2
      ctw4=ctw2**2
      stw2=1d0-ctw2
*-
      rhw = rhm2/rwm2
      rhz = rhm2/rzm2
      rhw2= rhw**2
      rhz2= rhz**2
*-
      Lnmuwm = DLOG(rwm2/tHmu2) 
      Lnmuzm = DLOG(rzm2/tHmu2)
      Lnmuhm = DLOG(rhm2/tHmu2)
*-
      b0fwwz = b0f(-rwm2,tHmu2,wm2,zm2)
      b0fzww = b0f(-rzm2,tHmu2,wm2,wm2)
      b0fwwh = b0f(-rwm2,tHmu2,wm2,hm2)
      b0fzzh = b0f(-rzm2,tHmu2,zm2,hm2)
      b0fww0 = b0f(-rwm2,tHmu2,wm2,cz2)
*-
      DRfer0=DREAL(SZZ_f(0d0)-SWW_f(0d0))/rwm2
      DelRho=3d0/4d0*(1d0/stw2*(Lnmuwm-Lnmuzm)
     &      +rhw/(1d0-rhw)*(Lnmuwm-Lnmuhm)
     &      -rhw/(1d0-rhz)*(Lnmuzm-Lnmuhm))-7d0/4d0+DRfer0
      DRh0Le=3d0/4d0*rtm2/rwm2
*-
      DRbos =DREAL((1d0/12d0/ctw4+4d0/3d0/ctw2-17d0/3d0-4d0*ctw2)
     &      *(b0fwwz-ctw2*b0fzww)
     &      +(1d0-1d0/3d0*rhw+1d0/12d0*rhw2)*b0fwwh
     &      -(1d0-1d0/3d0*rhz+1d0/12d0*rhz2)/ctw2*b0fzzh-4d0*stw2*b0fww0
     &      -1d0/12d0*(-(1d0/ctw4+6d0/ctw2-24d0+rhw)*Lnmuzm
     &      -rhw2*stw2*Lnmuhm
     &      +(14d0+1d0/ctw2+16d0*ctw2-48d0*ctw4+rhw)*Lnmuwm
     &      +1d0/ctw4+19d0/3d0/ctw2-22d0/3d0+stw2*rhw2))
*-
      DRfer =DREAL(SWW_f(rwm2)-SZZ_f(rzm2))/rwm2 
*-
      DelKap= ! -2d0/3d0*ve*(Lnmuzm-Lnmuem+1d0/6d0)
     &       +(1d0/6d0+7d0*ctw2)*Lnmuwm-8d0/9d0-2d0/3d0*ctw2
     &       -ctw2/stw2*(DRbos+DRfer)-DREAL(PZg_f(0d0))
      DRhoLe=DRbos+DRfer
      DKapLe=-ctw2/stw2*(DRbos+DRfer)
*-
      RETURN
      END

      SUBROUTINE INEETT(tHmu2A,tmA,wmA,zmA,hmA)
*-
      IMPLICIT NONE
*-
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      INTEGER*4 EWFFTR,GAMZTR,FERMTR,imf,TWOPOS
      REAL*8 PI,PI2,D2,D3,D5,CEILER
      REAL*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      REAL*8 tHmu2,tHmu2A
      REAL*8 AAv,vb,ab,qb,qe,AAe,ae,ae2,ve,ve2,anu,vnu,qe2,qt2,qeqt
      REAL*8 AAt,qt,qtm,at,at2,vt,vt2,vpae,vmae,vpau,vmau,de,dt,Nc
      REAL*8 rmf1,rmf2,rmf3,rmf,Qf,cf,I3f
      REAL*8 wm,wmA,zm,zmA,hm,hmA,tm,tmA,rwm2,rzm2,rhm2,rtm2
      REAL*8 ctw2,stw2,ctw4,ctw6
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,Rz,Rw  
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm
      REAL*8 DLOG
      COMPLEX*16 cz2,wm2,zm2,hm2,tm2,DCMPLX,Pgg_f
*-
      DIMENSION rmf1(12),rmf2(12),rmf3(12)
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-        nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
*-
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/TREATM/EWFFTR,GAMZTR,FERMTR
      COMMON/HISREM/TWOPOS
      COMMON/MATHCO/PI,PI2,D2,D3,D5,CEILER
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/tHScal/tHmu2
      COMMON/t_mass/tm,rtm2,tm2
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2
      COMMON/c_tw46/ctw4,ctw6
      COMMON/BOSrat/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,Rz,Rw
      COMMON/COUPLC/qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm
*-
      DATA rmf1/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &         .062d0,.083d0,1.50d0,.215d0,173.8d0,4.70d0/
      DATA rmf2/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &        .04851d0,.04851d0,1.50d0,.150d0,173.8d0,4.70d0/
*-             to reproduce ALQEDI=  128.88667 (Eidelman-Jegerlehner)
      DATA rmf3/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &        .075d0,.075d0,1.50d0,.250d0,173.8d0,4.70d0/
*-             to reproduce dispersive treatment of \Pi_{Z\gamma}(0)
      DATA Qf/0d0,-1d0,0d0,-1d0,0d0,-1d0,
     &        .666666666666667d0,-.333333333333333d0,
     &        .666666666666667d0,-.333333333333333d0,
     &        .666666666666667d0,-.333333333333333d0/
      DATA cf/6*1d0,6*3d0/
      DATA I3f/.5d0,-.5d0,.5d0,-.5d0,.5d0,-.5d0,.5d0,-.5d0,
     &         .5d0,-.5d0,.5d0,-.5d0/
      DATA alphaI/137.035989 5 D0/,CONHC/.389 379 66 D+9/
*-                                 conversion factor in pb's
*---------------------------------------- initialization --------------------
*-
      PI =4D0*DATAN(1D0)
      PI2=PI**2
      D2 =PI2/6D0
      D3 =1.2020569031596D0
      D5 =1.0369277551434D0
      CEILER=.577216D0      
*-
*-  Masses
*-
      rmf1(11)=tmA
      rmf2(11)=tmA
      rmf3(11)=tmA
*-
      if(FERMTR.eq.1) then
       do imf=1,12
        rmf(imf)=rmf1(imf)
       enddo
      elseif(FERMTR.eq.2) then
       do imf=1,12
        rmf(imf)=rmf2(imf)
       enddo
      elseif(FERMTR.eq.3) then
       do imf=1,12
        rmf(imf)=rmf3(imf)
       enddo
      endif
*-
      tm=tmA
      wm=wmA
      zm=zmA
      hm=hmA
*-  
      rwm2=wm**2
      rzm2=zm**2
      rtm2=tm**2
      rhm2=hm**2
*-
      cz2 =DCMPLX(0D0, -1D-20)
      wm2 =DCMPLX(rwm2,-1D-20)
      zm2 =DCMPLX(rzm2,-1D-20)
      tm2 =DCMPLX(rtm2,-1D-20)                 
      hm2 =DCMPLX(rhm2,-1D-20)                 
*-
*-  Mass ratios 
*-
      ctw2=rwm2/rzm2
      stw2=1d0-ctw2     
      ctw4=ctw2**2
      ctw6=ctw2**3
      rtz =rtm2/rzm2
      rtw =rtm2/rwm2
      rtw2=rtw**2
      rth =rtm2/rhm2
      rhz =rhm2/rzm2
      rhz2=rhz**2
      rhw =rhm2/rwm2
      rhw2=rhw**2
*-
      IF(TWOPOS.EQ.1) THEN
       GammaZ=2.499 776 D0 ! First  IPS
      ELSE
       GammaZ=2.499 538 D0 ! Second IPS
      ENDIF
*-
*-  For comparison with Hollik
*-
      if(GAMZTR.eq.0) GammaZ=0d0
*-
      if(IMOMS.eq.0) then 
       GFermi=pi/alphaI/sqrt(2d0)/stw2/ctw2/rzm2
      else
       GFermi=1.16637D-5
      endif      
*-
*-  Scales
*-
      tHmu2  = tHmu2A
      Lnmutm = DLOG(rtm2/tHmu2)
      Lnmuwm = DLOG(rwm2/tHmu2) 
      Lnmuzm = DLOG(rzm2/tHmu2)
      Lnmuhm = DLOG(rhm2/tHmu2)
*-
*-  Couplings
*-
      qt = 2d0/3
      qtm= qt
      qb =-1d0/3
      qe =-1d0
*-
      qe2=qe**2
      qt2=qt**2
      qeqt=qe*qt
*-
      anu=1d0/2
      vnu=1d0/2
*-
      ae =-1d0/2
      ve =-1d0/2-2d0*qe*stw2
      ae2=ae**2
      ve2=ve**2
      AAe=(3*ve2+ae2)*ae
*-
      at =1d0/2
      vt =1d0/2-2d0*qt*stw2
      at2=at**2
      vt2=vt**2
      AAt=(3*vt2+at2)*at
      AAv=(vt2+3*at2)*vt
*-
      ab =-1d0/2
      vb =-1d0/2-2d0*qb*stw2
*-
      vpae=ve+ae
      vmae=ve-ae
      vpau=vt+at
      vmau=vt-at
      de=-2d0*vmae
      dt=+2d0*vmau
      Nc=3
*-
      RETURN 
      END

      COMPLEX*16 FUNCTION SWW_f(s)
      IMPLICIT NONE
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,Qf,cf,I3f
      REAL*8 ctw2,stw2
      REAL*8 DLOG
      COMPLEX*16 bff,bfmdmu,b1mumd,b1mdmu,b0f,b0fsdu,b0p
      COMPLEX*16 mu2,md2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/c_stw2/ctw2,stw2
*-
      SWW_f=DCMPLX(0d0,0d0)
      do 1 i=1,6  
      mu2=DCMPLX((rmf(2*i-1))**2,-1d-20)
      md2=DCMPLX((rmf(2*i)  )**2,-1d-20)
      bfmdmu=bff(-s,tHmu2,md2,mu2)
      b0fsdu=b0f(-s,tHmu2,md2,mu2)
      if(s.ne.0d0) then
       b1mdmu=-1d0/2/s*                       
     &        (md2*(DLOG(DREAL(md2)/tHmu2)-1d0)                
     &        -mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)        
     &        +(md2-mu2+s)*b0fsdu)
       b1mumd=-1d0/2/s*
     &        (mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)                
     &        -md2*(DLOG(DREAL(md2)/tHmu2)-1d0)        
     &        +(mu2-md2+s)*b0fsdu) 
      else
       b1mdmu=-1d0/2*b0f(0d0,tHmu2,md2,mu2)
     &        +1d0/2*(md2-mu2)*b0p(0d0,md2,mu2)
       b1mumd=-1d0/2*b0f(0d0,tHmu2,mu2,md2)
     &        +1d0/2*(mu2-md2)*b0p(0d0,mu2,md2)
      endif        
      SWW_f =SWW_f - s*cf(2*i)*bfmdmu
     &       +cf(2*i-1)*mu2*b1mdmu+cf(2*i)*md2*b1mumd
 1    continue
      END

      COMPLEX*16 FUNCTION SZZ_f(s)
      IMPLICIT NONE
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,rmf2,Qf,cf,I3f
      REAL*8 af,af2,vf,vf2
      REAL*8 ctw2,stw2
      COMPLEX*16 bff,b0f,bffsff,b0fsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/c_stw2/ctw2,stw2
*-
      SZZ_f=DCMPLX(0d0,0d0)
      do 1 i=1,12
      rmf2=(rmf(i))**2
      mf2 =DCMPLX(rmf2,-1d-20)     
      bffsff = bff(-s,tHmu2,mf2,mf2)   
      b0fsff = b0f(-s,tHmu2,mf2,mf2)   
      af =I3f(i)                       
      vf =I3f(i)-2d0*Qf(i)*stw2   
      af2=af**2                        
      vf2=vf**2
      SZZ_f=SZZ_f+cf(i)*(-(vf2+af2)*s*bffsff-2d0*af2*rmf2*b0fsff)
 1    continue
      END

      COMPLEX*16 FUNCTION Pgg_f(s)
      IMPLICIT NONE
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,rmf2,Qf,cf,I3f
      COMPLEX*16 bff,bffsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-
      Pgg_f =DCMPLX(0d0,0d0)
      do 1 i=1,12  
      rmf2=(rmf(i))**2
      mf2=DCMPLX(rmf2,-1d-20)     
      bffsff=bff(-s,tHmu2,mf2,mf2)
      Pgg_f =Pgg_f+4d0*cf(i)*Qf(i)**2*bffsff
 1    continue
      END

      COMPLEX*16 FUNCTION PZg_f(s)
      IMPLICIT NONE
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,Qf,Qfm,cf,I3f
      REAL*8 ctw2,stw2 
      COMPLEX*16 bff,bffsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/c_stw2/ctw2,stw2
*-
      PZg_f =DCMPLX(0d0,0d0)
      do 1 i=1,12
      mf2=DCMPLX((rmf(i))**2,-1d-20)     
      bffsff=bff(-s,tHmu2,mf2,mf2)   
      Qfm   =abs(Qf(i))
      PZg_f =PZg_f+cf(i)*Qfm*(1-4*stw2*Qfm)*bffsff
 1    continue
      END

      COMPLEX*16 FUNCTION B0F(Q2,MU2,M12,M22)
*     ---------------------------------------
* general B0F; CalcPHEP group
*
      IMPLICIT NONE
      REAL*8 MU2,Q2,DREAL
      COMPLEX*16 M12,M22,SQR,SQS,ILQ,LOG,SQRT
*
      SQR=SQRT(Q2**2+2D0*Q2*(M12+M22)+(M12-M22)**2)
      SQS=Q2+M12+M22
*
      IF(ABS(Q2).NE.0D0) THEN
* B0F(Q2;MU2,M12,M22)
       IF(DREAL(M12).NE.0D0.AND.DREAL(M22).NE.0D0) THEN
        IF(DREAL((SQS+SQR)/SQR).LT.1D-7) THEN
          ILQ=-SQR/Q2*LOG((SQS-SQR)/(2D0*SQRT(M12*M22)))
        ELSE
          ILQ=+SQR/Q2*LOG((SQS+SQR)/(2D0*SQRT(M12*M22)))
        ENDIF
        B0F=2D0-LOG(SQRT(M12*M22)/MU2)+(M12-M22)/2D0/Q2*LOG(M12/M22)-ILQ
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        B0F=2D0-LOG(M22/MU2)-(1D0+M22/Q2)*LOG(1D0+Q2/M22)
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0F=2D0-LOG(M12/MU2)-(1D0+M12/Q2)*LOG(1D0+Q2/M12)
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0F=2D0-LOG(DCMPLX(Q2,-1D-20)/MU2)
       ELSE
        PRINT*,'NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP
       ENDIF
      ELSE
* B0F(0;MU2,M12,M22)
       IF((DREAL(M12).NE.DREAL(M22))
     &   .AND.DREAL(M22).NE.0D0.AND.DREAL(M12).NE.0D0) THEN
        B0F=1D0-LOG(M22/MU2)-M12/(M12-M22)*LOG(M12/M22)
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        B0F=1D0-LOG(M22/MU2)
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0F=1D0-LOG(M12/MU2)
       ELSEIF(DREAL(M12).EQ.DREAL(M22)) THEN
        B0F=-LOG(M12/MU2)
       ELSE
        PRINT*,'B0F(0;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP        
       ENDIF
      ENDIF
*
      RETURN    
      END

      COMPLEX*16 FUNCTION B0P(Q2,M12,M22)
*     -----------------------------------
* general B0P; CalcPHEP group
*
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Q2,DREAL
      COMPLEX*16 M12,M22,SD,B,XPL,XMI,FB0P,SUM,DCMPLX,LOG,SQRT
*
      IF(ABS(Q2).NE.0D0) THEN
* B0P(Q2;M12,M22)
       SD=SQRT(Q2**2+2D0*Q2*(M12+M22)+(M12-M22)**2)
       B=Q2-M22+M12
       XPL=(-B+SD)/(-2D0*Q2)
       XMI=(-B-SD)/(-2D0*Q2)
       IF(DREAL(M12).NE.0D0.AND.DREAL(M22).NE.0D0) THEN
        SUM=DCMPLX(0D0,0D0)
        DO N=1,2
          SUM=SUM+(-1)**N*(FB0P(N,XPL)-FB0P(N,XMI))
        ENDDO
        B0P=-1D0/SD*SUM
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        B0P=-1D0/Q2**2*(Q2-M22*LOG(1D0+Q2/M22))
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0P=-1D0/Q2**2*(Q2-M12*LOG(1D0+Q2/M12))
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0P=-1D0/Q2
       ELSE
        PRINT *,'NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP
       ENDIF        
      ELSE
* B0P(0;M12,M22)
       IF((DREAL(M12).NE.DREAL(M22))
     &   .AND.DREAL(M22).NE.0D0.AND.DREAL(M12).NE.0D0) THEN
        B0P=M12*M22/(M12-M22)**3*LOG(M12/M22)-(M12+M22)/2D0/(M12-M22)**2
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        B0P=-1D0/2D0/M22
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        B0P=-1D0/2D0/M12
       ELSEIF(DREAL(M12).EQ.DREAL(M22)) THEN
        B0P=-1D0/6D0/M12
       ELSE
        PRINT*,'B0P(0;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP        
       ENDIF
      ENDIF
*
      RETURN    
      END

      COMPLEX*16 FUNCTION FB0P(N,Y)
*
      IMPLICIT NONE
      INTEGER*4 N,L
      COMPLEX*16 Y,SUM,DCMPLX
*
      SUM=DCMPLX(0D0,0D0)
      DO L=1,N
        SUM=SUM+Y**(N-L)/L
      ENDDO
      FB0P=-Y**N*LOG(1D0-1D0/Y)-SUM
*
      RETURN
      END
      COMPLEX*16 FUNCTION BFF(Q2,MU2,M12,M22)
*     ---------------------------------------
* general B0F; CalcPHEP group
*
      IMPLICIT NONE
      REAL*8 MU2,Q2,DREAL
      COMPLEX*16 M12,M22,B,SQR,SQS,ILQ,LOG,SQRT
*
      SQR=SQRT(Q2**2+2D0*Q2*(M12+M22)+(M12-M22)**2)
      SQS=Q2+M12+M22
*
      IF(ABS(Q2).NE.0D0) THEN
* BFF(Q2;MU2,M12,M22)
       IF(DREAL(M12).NE.0D0.AND.DREAL(M22).NE.0D0) THEN
        IF(DREAL((SQS+SQR)/SQR).LT.1D-7) THEN
          ILQ=-SQR/Q2*LOG((SQS-SQR)/(2D0*SQRT(M12*M22)))
        ELSE
          ILQ=+SQR/Q2*LOG((SQS+SQR)/(2D0*SQRT(M12*M22)))
        ENDIF
        BFF=1D0/3D0*LOG(SQRT(M12*M22)/MU2)-5D0/9D0
     &     +2D0/3D0*((M12+M22)/Q2+(M12-M22)**2/Q2**2)
     &     +((M12**2-M22**2)/2D0/Q2**2
     &     +(M12-M22)**3/3D0/Q2**3)*LOG(M12/M22)
     &     +(1D0-(M12+M22)/Q2-2D0*(M12-M22)**2/Q2**2)/3D0*ILQ
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        BFF=1D0/3D0*LOG(M22/MU2)-5D0/9D0+2D0/3D0*(M22/Q2+M22**2/Q2**2)
     &   +(1D0-M22/Q2-2D0*M22**2/Q2**2)/3D0*(1D0+M22/Q2)*LOG(1D0+Q2/M22)
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        BFF=1D0/3D0*LOG(M12/MU2)-5D0/9D0+2D0/3D0*(M12/Q2+M12**2/Q2**2)
     &   +(1D0-M12/Q2-2D0*M12**2/Q2**2)/3D0*(1D0+M12/Q2)*LOG(1D0+Q2/M12)
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).EQ.0D0) THEN
        BFF=1D0/3D0*LOG(DCMPLX(Q2,-1D-20)/MU2)-5D0/9D0
       ELSE
        PRINT *,'NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP
       ENDIF
      ELSE
* BFF(0;MU2,M12,M22)
       IF((DREAL(M12).NE.DREAL(M22))
     &   .AND.DREAL(M22).NE.0D0.AND.DREAL(M12).NE.0D0) THEN
        B=M12/M22-1D0
        BFF=1D0/3D0*LOG(M22/MU2)
     &     +(1D0/3D0-1D0/B**2-2D0/3D0/B**3)*LOG(M12/M22)
     &     -5D0/18D0+2D0/3D0/B+2D0/3D0/B**2
       ELSEIF(DREAL(M12).EQ.0D0.AND.DREAL(M22).NE.0D0) THEN
        BFF=1D0/3D0*LOG(M22/MU2)-5D0/18D0
       ELSEIF(DREAL(M12).NE.0D0.AND.DREAL(M22).EQ.0D0) THEN
        BFF=1D0/3D0*LOG(M12/MU2)-5D0/18D0
       ELSEIF(DREAL(M12).EQ.DREAL(M22)) THEN
        BFF=1D0/3D0*LOG(M12/MU2)
       ELSE
        PRINT*,'BFF(0;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP        
       ENDIF
      ENDIF
*
      RETURN    
      END
