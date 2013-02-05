*////////////////////////////////////////////////////////////////////////////
*//                                                                        //
*////////////////////////////////////////////////////////////////////////////
*
      INTEGER    imax
      PARAMETER (imax=1000)
      DOUBLE PRECISION   m_pi
      PARAMETER         (m_pi =3.1415926535897932d0)
*
      DOUBLE PRECISION    m_xpar,   m_CMSene,   m_amel,   m_amfin
      DOUBLE PRECISION    m_gnanob, m_Chini,  m_Chfin
      DOUBLE PRECISION    m_vvmin,  m_vvmax,  m_Xenph
      DOUBLE PRECISION    m_Zmass,  m_Zgamma, m_sinw2
      DOUBLE PRECISION    m_alfinv, m_alfinvMZ, m_GFermi
      DOUBLE PRECISION    m_svar1,  m_cMax,   m_cMin,     m_MfinMin
      DOUBLE PRECISION    m_eps1,   m_eps2,     m_ta1,    m_ta2
      INTEGER             m_KeyZet, m_KeyDis,   m_out
      INTEGER             m_KFini,  m_KFfin,    m_KeyFoB, m_KeyFSR
      INTEGER             m_Nchanel,m_KeyQCD
*
      COMMON  /c_Semalib/
     $    m_xpar(imax),               ! copy of input
     $    m_CMSene,                   ! CMS energy
     $    m_amel,                     ! beam  mass
     $    m_amfin,                    ! final mass
     $    m_Chini,                    ! beam  charge
     $    m_Chfin,                    ! final charge
     $    m_MfinMin,                  ! smallest final mass
     $    m_Zmass,                    ! Z mass
     $    m_Zgamma,                   ! Z width
     $    m_sinw2,                    ! EW mixing angle
     $    m_GFermi,                   ! Gmu Fermi
     $    m_alfinv,                   ! 1/alphaQED at q^2=0
     $    m_alfinvMZ,                 ! 1/alphaQED at MZ
     $    m_gnanob,                   ! GeV^2 -> nanobarns
     $    m_vvmin,                    ! minimum v
     $    m_vvmax,                    ! maximum v
     $    m_svar1,                    ! propagator sprim
     $    m_eps1,                     ! longit. polarization 1-st beam
     $    m_eps2,                     ! longit. polarization 2-nd beam
     $    m_ta1,                      ! longit. polarization 1-st final  ferm.
     $    m_ta2,                      ! longit. polarization 2-nd final  ferm.
     $    m_cMax,                     ! cos(theta_max)
     $    m_cMin,                     ! cos(theta_min)
     $    m_Xenph,                    ! Enhancement factor for Crude photon multiplicity
     $    m_KeyZet,                   ! type of Born
     $    m_KeyDis,                   ! type of bremsstrahlung
     $    m_KFini,                    ! KF flavour  code  1-st beam
     $    m_KFfin,                    ! KF flavour  code  1-st final fermion
     $    m_KeyFoB,                   ! togle swith for Forward or Backward
     $    m_KeyFSR,                   ! FSR on/off
     $    m_KeyQCD,                   ! QCD factor for quarks on/off
     $    m_Nchanel,                  ! No of open chanels
     $    m_out                       ! output file number
*
*////////////////////////////////////////////////////////////////////////////
*//                                                                        //
*////////////////////////////////////////////////////////////////////////////
 
