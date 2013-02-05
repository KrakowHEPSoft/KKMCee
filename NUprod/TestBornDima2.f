* Very simple program for e-e+ --> nu nubar  cross section
*     (g77 TestBornDima2.f; a.out)
      PROGRAM Main
*-------------------------------------------------------------
* Note that we define Born = 3/8 *dSigma/dCosTheta(CosTheta=0)
*-------------------------------------------------------------
* Output of this program is for nu_muon:
*  ========= test ============
*Test: CMSene =   91.1900  Born =     3.94526031817576
*Test: CMSene =  100.0000  Born =     0.08497608383105
*Test: CMSene =  140.0000  Born =     0.00377912901529
*Test: CMSene =  140.0000  Born =     0.00377912901529
*Test: CMSene =  189.0000  Born =     0.00116819286268
*Test: CMSene =  200.0000  Born =     0.00097868303435
*Test: CMSene =  206.0000  Born =     0.00089533050038
* Output for nu_electron:
* Comment out line BornT  = 0d0 !!!
*Test: CMSene =   91.1900  Born =     3.95460013827135
*Test: CMSene =  100.0000  Born =     0.13638590504141
*Test: CMSene =  140.0000  Born =     0.02138832243774
*Test: CMSene =  140.0000  Born =     0.02138832243774
*Test: CMSene =  189.0000  Born =     0.01272366028692
*Test: CMSene =  200.0000  Born =     0.01166075454654
*Test: CMSene =  206.0000  Born =     0.01114178091718
*-------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION   m_pi
      PARAMETER         (m_pi =3.1415926535897932d0)
      DOUBLE PRECISION Ve,Ae,Vf,Af,Qf,Qe,sinw2,costhe,CMSene,s,t
      DOUBLE PRECISION Gfermi,MZ,AlfInv,alfa,GZ,gnanob
      DOUBLE PRECISION Born,ff1,ff0,xe,xf,ye,yf
      DOUBLE PRECISION SqChiZ,ReChiZ,SqChiW,ReChiW,MW,BornT,BornS
      DOUBLE PRECISION RaZ
      INTEGER iEne
      WRITE(*,*) ' ========= test ============ '
*-------------------------------------------------------------
* input adjusted to Zfitter
*-------------------------------------------------------------
      Gfermi =   0.00001166370000d0 ! Last digit is 7 not 9!
      sinw2  =   0.22330592993199d0
      MZ     =  91.18700000000000d0
      GZ     =    2.4953427050359d0
      MW     =    80.363329699738d0
      AlfInv = 137.03598950000000d0 ! Thomson
***** AlfInv = 128.88524d0          ! at MZ
*-------------------------------------------------------------
      alfa   = 1d0/AlfInv
      gnanob =  389.37966d3
      Qe= -1d0                  !electron
      Qf=  0d0                  !neutrino
      Ve =  -1d0 -4*Qe*sinw2    !electron
      Ae =  -1d0                !electron
      Vf =  +1d0 -4*Qf*sinw2    !neutrino
      Af =  +1d0                !neutrino
      costhe = 0d0
      DO iEne =1,7
         IF(iEne.EQ.1)  CMSene = 91.19d0
         IF(iEne.EQ.2)  CMSene = 100d0
         IF(iEne.EQ.3)  CMSene = 140d0
         IF(iEne.EQ.5)  CMSene = 189d0
         IF(iEne.EQ.6)  CMSene = 200d0
         IF(iEne.EQ.7)  CMSene = 206d0
         s = CMSene**2
         xe= Ve**2 +Ae**2
         xf= Vf**2 +Af**2
         ye= 2*Ve*Ae
         yf= 2*Vf*Af
         RaZ  =(GFermi *s *AlfInv     )/( DSQRT(2d0) *8d0 *m_pi) !

         ReChiZ =(s-MZ**2)*MZ**2/((s-MZ**2)**2+(GZ*MZ)**2) *RaZ ! fixed width
         SqChiZ =          MZ**4/((s-MZ**2)**2+(GZ*MZ)**2) *RaZ**2 ! fixed width

* t-chanel
         t =-s*(1d0-costhe)/2d0
         ReChiW=      MW**2/(-t+MW**2)    *RaZ
         SqChiW=      MW**4/(-t+MW**2)**2 *RaZ**2
         BornT  = 16d0*SqChiW*(1d0+costhe)**2
     $            -8d0*ReChiZ*ReChiW*(Ae+Ve)*(1d0+costhe)**2
cc         BornT  = 0d0
* s-chanel
         ff0= Qe**2*Qf**2 +2*ReChiZ*Qe*Qf*Ve*Vf +SqChiZ*xe*xf
         ff1=             +2*ReChiZ*Qe*Qf*Ae*Af +SqChiZ*ye*yf
         BornS = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1 ! sigma=3/8*dsigma/dcostheta(0)

         Born = (BornS+BornT)*  4d0*m_pi*alfa**2/(3d0*s)  *gnanob
         WRITE(*,'(a,f10.4,a,f20.14)') 
     $        'Test: CMSene =',CMSene,'  Born = ', Born
      ENDDO
      END

