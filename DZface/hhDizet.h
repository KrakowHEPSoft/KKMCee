* /////////////////////////////////////////////////////////////////////////////
* //                                                                         //
* //           Pseudo-Class hh_Dizet for Dizet interface                    //
* //                     to KKMC_hh. Allows multiple initial states.        //
* //                                                                         //
* /////////////////////////////////////////////////////////////////////////////
*//////////////////////////////////////////////////////////////////////////////
*//       Energy limits in the EW grid, w=sqrt(s) in GeV units.              //
*//////////////////////////////////////////////////////////////////////////////
*     340-point grid, only 80pt for NLC to be improved/tested in future
      INTEGER           m_poin1, m_poin2, m_poin3, m_poin4, m_poinG , m_poinQ !
      INTEGER           m_poTh1, m_poTh2, m_poTh3, m_poTh4
      PARAMETER(        m_poinG =   7 )                      ! No of EW formfactors
      PARAMETER(        m_poinQ =   4 )                      ! No of QCD corrections
*----------- Low energies and LEP1
      PARAMETER(        m_poin1 = 100 )                              ! LEP1 LOG(SQRT(s)) spacing
      PARAMETER(        m_poTh1 =   0 )                              ! Cost(heta) spacing
      DOUBLE PRECISION  m_WminLEP1,          m_WmaxLEP1              ! LEP1 basic range (m_WminLEP1,m_WmaxLEP1)
      PARAMETER(        m_WminLEP1=0.010d0,  m_WmaxLEP1= 95.000d0 )  ! LEP1 basic range (m_WminLEP1,m_WmaxLEP1)
*----------- Z resonance
      PARAMETER(        m_poin2 = 120 )                              ! Z range sqrt(s)    spacing
      PARAMETER(        m_poTh2 =  14 )                              ! =14 is overkill?
      DOUBLE PRECISION  m_WminZ, m_WmaxZ, m_WdelZ                    ! Z range (amz + m_WdelZ)
      PARAMETER(        m_WdelZ = 60.000d0)                          ! Z range (amz + m_WdelZ)
*----------- LEP2
      PARAMETER(        m_poTh3 =  30 )                              ! Overkill, bit lets kkep it
      PARAMETER(        m_poin3 = 145 )                              ! LEP2 interval sqrt(s)    spacing
      DOUBLE PRECISION  m_WmaxLEP2                                   ! LEP2 interval (m_WmaxLEP1,m_WmaxLEP2)
      PARAMETER(        m_WmaxLEP2  =240.001d0 )                     ! LEP2 interval (m_WmaxLEP1,m_WmaxLEP2)
*----------- Linear Colliders
      PARAMETER(        m_poin4 = 180 )                              ! NLC range sqrt(s)    spacing
      PARAMETER(        m_poTh4 =  14 )                              ! Cost(heta) spacing
      DOUBLE PRECISION  m_WmaxNLC                                    ! NLC range (m_WmaxLEP2,m_WmaxNLC)
      PARAMETER(        m_WmaxNLC  =8000.00d0 )                      ! NLC range (m_WmaxLEP2,m_WmaxNLC)
*//////////////////////////////////////////////////////////////////////////////
      DOUBLE COMPLEX     m_GSW
      DOUBLE PRECISION   m_QCDcor   ! obsolete
      DOUBLE PRECISION   m_QCDcorR
      COMMON /c_hhDizetG/
     $  m_GSW(    m_poinG),                    ! form-factors,   at the actual energy/angle
     $  m_QCDcorR(m_poinQ),                    ! QCD correction, at the actual energy/angle
     $  m_QCDcor                               ! obsolete!!!!
*//////////////////////////////////////////////////////////////////////////////
      DOUBLE COMPLEX    m_cyys, m_czzs, m_ctts, m_clcs
      DOUBLE PRECISION  m_syys, m_szzs, m_stts, m_slcs
      INTEGER m_ndisk
      DOUBLE PRECISION m_xpar
      COMMON /c_hhDizetF/
     &  m_cyys(m_poin1+1,          m_poinG,5,16), ! EW formfactor table
     &  m_czzs(m_poin2+1,m_poTh2+1,m_poinG,5,16), ! EW formfactor table
     &  m_ctts(m_poin3+1,m_poTh3+1,m_poinG,5,16), ! EW formfactor table
     &  m_clcs(m_poin4+1,m_poTh4+1,m_poinG,5,16), ! EW formfactor table
     $  m_syys(m_poin1+1,          m_poinQ,5,16), ! QCD correction,
     $  m_szzs(m_poin2+1,          m_poinQ,5,16), ! QCD correction,
     $  m_stts(m_poin3+1,          m_poinQ,5,16), ! QCD correction,
     $  m_slcs(m_poin3+1,          m_poinQ,5,16), ! QCD correction,
     &  m_xpar(10000),
     &  m_ndisk
      SAVE /c_hhDizetG/
      SAVE /c_hhDizetF/
