*////////////////////////////////////////////////////////////////////////////
*//                                                                        //
*////////////////////////////////////////////////////////////////////////////
*
      INTEGER    imax
      PARAMETER (imax=1000)
*
      DOUBLE PRECISION    m_xpar,   m_CMSene,   m_amel,   m_amfin
      DOUBLE PRECISION    m_Zmass,  m_vvmin,    m_vvmax,  m_Xenph
      INTEGER             m_KeyZet, m_KeyDis,   m_out
*
      COMMON  /c_Semalib/
     $    m_xpar(imax),               ! copy of input
     $    m_CMSene,                   ! CMS energy
     $    m_amel,                     ! electron mass
     $    m_amfin,                    ! final (muon) mass
     $    m_Zmass,                    ! Z mass
     $    m_vvmin,                    ! minimum v
     $    m_vvmax,                    ! maximum v
     $    m_Xenph,                    ! Enhancement factor for Crude photon multiplicity
     $    m_KeyZet,                   ! type of Born
     $    m_KeyDis,                   ! type of bremsstrahlung
     $    m_out                       ! output file number
*
*////////////////////////////////////////////////////////////////////////////
*//                                                                        //
*////////////////////////////////////////////////////////////////////////////
 
