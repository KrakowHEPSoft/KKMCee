*     Very simple program for e-e+ --> mu-mu+ cross section
*     (g77 TestBornDima.f; a.out)
      PROGRAM Main
*--------------------------------------------
* Output of this program is (for input from Dima):
*  ========= test ============
*Test: CMSene =   91.1900  Born =     2.00557303705243
*Test: CMSene =  100.0000  Born =     0.05208955209159
*Test: CMSene =  140.0000  Born =     0.00640805872136
*Test: CMSene =  140.0000  Born =     0.00640805872136
*Test: CMSene =  189.0000  Born =     0.00304924583604
*Test: CMSene =  200.0000  Born =     0.00268965550654
*Test: CMSene =  206.0000  Born =     0.00252119053337
*---------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION   m_pi
      PARAMETER         (m_pi =3.1415926535897932d0)
      DOUBLE PRECISION Ve,Ae,Vf,Af,Qf,Qe,sinw2,costhe,CMSene,s
      DOUBLE PRECISION Gfermi,MZ,AlfInv,alfa,GZ,gnanob,MW
      DOUBLE PRECISION Razy,chi2,rechi,Born,ff1,ff0,xe,xf,ye,yf
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
      Qe= -1d0
      Qf= -1d0
      Ve =  -1d0 -4*Qe*sinw2
      Ae =  -1d0
      Vf =  -1d0 -4*Qf*sinw2
      Af =  -1d0
      costhe = 0d0
      DO iEne =1,7
         IF(iEne.EQ.1)  CMSene = 91.19d0
         IF(iEne.EQ.2)  CMSene = 100d0
         IF(iEne.EQ.3)  CMSene = 140d0
         IF(iEne.EQ.5)  CMSene = 189d0
         IF(iEne.EQ.6)  CMSene = 200d0
         IF(iEne.EQ.7)  CMSene = 206d0
         s = CMSene**2
*         Qe= 0d0                ! <-- this switches of Z-exchange
*         Qf= 0d0                ! <-- this switches of Z-exchange
         xe= Ve**2 +Ae**2
         xf= Vf**2 +Af**2
         ye= 2*Ve*Ae
         yf= 2*Vf*Af
         Razy =(GFermi *MZ**2 *AlfInv )/( DSQRT(2d0) *8d0 *m_pi)

        rechi=(s-MZ**2)*s/((s-MZ**2)**2+(GZ*MZ)**2) *Razy    ! fixed width
        chi2 =       s**2/((s-MZ**2)**2+(GZ*MZ)**2) *Razy**2 ! fixed width
*         rechi=(s-MZ**2)*s/((s-MZ**2)**2+(GZ*s/MZ)**2)*Razy
*         chi2 =       s**2/((s-MZ**2)**2+(GZ*s/MZ)**2)*Razy**2

         ff0= Qe**2*Qf**2 +2*rechi*Qe*Qf*Ve*Vf +chi2*xe*xf
         ff1=             +2*rechi*Qe*Qf*Ae*Af +chi2*ye*yf
         Born = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1 ! diff. x-sec. *3/8
         Born = Born*  4d0*m_pi*alfa**2/(3d0*s)  *gnanob
         WRITE(*,'(a,f10.4,a,f20.14)') 
     $        'Test: CMSene =',CMSene,'  Born = ', Born
      ENDDO
      END


