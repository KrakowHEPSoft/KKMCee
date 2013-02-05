
c-------------- Couplings, masses and meson widths ------------------
      alpha =137.03599976d0     !   ---------- 1/alpha (QED)
      mpi = 0.13956995d0        ! GeV ---------- Charged pion mass
      al = 1.48d-3              ! -------------- Pion form factor parameters a and b: 
      be = -0.1473d0            ! -------------- Kuehn, Santamaria, ZPC48(1990)445
      mrho = 0.7726d0           ! GeV ---------- Rho mass
      mrhol = 1.46d0            ! GeV ---------- Rho' mass
      gammarho = 0.1437d0       ! GeV ---------- Total rho width 
      grhol = 0.31d0            ! GeV ---------- Rho' width
      momega = 0.78278d0        ! GeV ---------- Omega mass
      gomega = 8.68d-3          ! GeV ---------- Omega width
c --------------------------------------------------------------------
c     this is the main line for the call of the pion form factor
c     additional subroutines are:
c     - PionFormFactor2
c     - PionFormFactor
c     - BW
      pionFF = 4.d0*pi*alpha*PionFormFactor2(qq)
c --------------------------------------------------------------------

      double precision function PionFormFactor2(a)
      include 'phokhara_2.0.inc'       
      double precision a
      complex*16 PionFormFactor
       
      PionFormFactor2= PionFormFactor(a)*dconjg(PionFormFactor(a))
c
      return
      end

c --------------------------------------------------------------------

      complex*16 function PionFormFactor(a)
      include 'phokhara_2.0.inc'       
      double precision a
      complex*16 BW
       
      PionFormFactor = (BW(mrho,gammarho,a,1)
     &     *(1.D0+al*BW(momega,gomega,a,1))/(1.d0+al)
     &     +be*BW(mrhol,grhol,a,1))/(1.d0+be)
      return
      end

c --------------------------------------------------------------------

      complex*16 function BW(m,breite,x,k)
      include 'phokhara_2.0.inc'       
      integer k
      double precision m,breite,x,g
      complex *16 i

      if(breite.eq.gomega)then 
         g=breite
      else
         g=breite*m*m/x*(x-4.d0*mpi*mpi)**(1.5d0)/
     &     (m*m-4.d0*mpi*mpi)**(1.5d0)
      endif
      i=(0.d0,1.d0)
      if(k.eq.1)then
         BW=m*m/(m*m-x-i*sqrt(x)*g)
      else
         BW=m*m/(m*m-x+i*sqrt(x)*g)
      endif
      end

c --------------------------------------------------------------------


