      PROGRAM Main
*--------------------------------------------
* Output of this program is:
* Test: CMSene =  189.  Born =   0.590062067
* Test: CMSene =  200.  Born =   0.494369669
* Test: CMSene =  206.  Born =   0.452277836
*---------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION   m_pi
      PARAMETER         (m_pi =3.1415926535897932d0)
      DOUBLE PRECISION Ve,Ae,Vf,Af,Qf,Qe,sinw2,costhe,CMSene,s
      DOUBLE PRECISION Gfermi,MZ,AlfInv,alfa,GZ,gnanob
      DOUBLE PRECISION Razy,chi2,rechi,Born,ff1,ff0,xe,xf,ye,yf
      INTEGER iEne
      WRITE(*,*) ' ========= test ============ '
      sinw2  = 0.22330592985051D0   ! more ditits added
      Gfermi = 1.16637d-5     ! With two-loop QED 
      MZ  = 91.187D0
      GZ  = 2.495343D0
      AlfInv = 137.0359895d0 ! Thomson
*      AlfInv = 128.88524d0   ! at MZ
      alfa   = 1d0/AlfInv
      gnanob =  389.37966d3
      Qe= -1d0
      Qf= -1d0
      Ve =  -1d0 -4*Qe*sinw2
      Ae =  -1d0
      Vf =  -1d0 -4*Qf*sinw2
      Af =  -1d0
      costhe = 0d0
      DO iEne =1,3
         IF(iEne.EQ.1)  CMSene = 189d0
         IF(iEne.EQ.2)  CMSene = 200d0
         IF(iEne.EQ.3)  CMSene = 206d0
         s = CMSene**2
*         Qe= 0d0                ! <-- this switches of Z-exchange
*         Qf= 0d0                ! <-- this switches of Z-exchange
***** rechi=(s-MZ**2)*s/((s-MZ**2)**2+(GZ*MZ)**2) *Razy    ! fixed width
***** chi2 =       s**2/((s-MZ**2)**2+(GZ*MZ)**2) *Razy**2 ! fixed width
         xe= Ve**2 +Ae**2
         xf= Vf**2 +Af**2
         ye= 2*Ve*Ae
         yf= 2*Vf*Af
         Razy =(GFermi *MZ**2 *AlfInv )/( DSQRT(2d0) *8d0 *m_pi)
         rechi=(s-MZ**2)*s/((s-MZ**2)**2+(GZ*s/MZ)**2)*Razy
         chi2 =       s**2/((s-MZ**2)**2+(GZ*s/MZ)**2)*Razy**2
         ff0= Qe**2*Qf**2 +2*rechi*Qe*Qf*Ve*Vf +chi2*xe*xf
         ff1=             +2*rechi*Qe*Qf*Ae*Af +chi2*ye*yf
         Born = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1 ! diff. x-sec. *3/8
         Born = Born*  4d0*m_pi*alfa**2/(3d0*s)  *gnanob*1d3
         WRITE(*,*) 'Test: CMSene =',CMSene,'  Born = ', Born
      ENDDO
      END


