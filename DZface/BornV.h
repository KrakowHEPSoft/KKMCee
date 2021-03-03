*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  BornV                                  //
*//                                                                          //
*//  Purpose:                                                                //
*//  Provide Born angular distribution and integrated x-section              //
*//  as a function of s.                                                     //
*//                                                                          //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
*
*  Class members:
*
*//////////////////////////////////////////////////////////////////////////////
      DOUBLE PRECISION   m_pi
      PARAMETER         (m_pi =3.1415926535897932d0)
      DOUBLE PRECISION   m_fleps
****  PARAMETER (m_fleps = 1d-35)  ! original
****  PARAMETER (m_fleps = 1d-45)  ! enough???
      PARAMETER (m_fleps = 1d-100)  ! enough!!!
****      PARAMETER (m_fleps = 1d-200)  ! enough!!!
* EW parameters
      DOUBLE PRECISION   m_Gmu
      DOUBLE PRECISION   m_MZ,      m_amh,     m_amtop
      DOUBLE PRECISION   m_swsq,    m_gammz,   m_MW,    m_GammW
*
      DOUBLE PRECISION   m_CMSene,  m_XXXene,  m_HadMin, m_vvmin,  m_vvmax !
      DOUBLE PRECISION   m_AvMult,  m_YFSkon,  m_YFS_IR, m_alfinv, m_alfpi, m_Xenph !
      DOUBLE PRECISION   m_vv,      m_x1,      m_x2
      DOUBLE PRECISION   m_Qf,      m_T3f,     m_helic,  m_amferm, m_auxpar !
      DOUBLE PRECISION   m_gnanob,  m_xpar_input
      DOUBLE PRECISION   m_alfinvMZ,   m_alfQCDMZ
      INTEGER            m_IsGenerated, m_KFferm,  m_NCf
      INTEGER            m_KFini,       m_KFfin,   m_KeyINT,   m_KeyQCD,   m_KeyRes
      INTEGER            m_KeyElw,      m_KeyZet,  m_KeyWtm,   m_KeyFSR
      INTEGER            m_out,         m_ibox
*[[[ NEW!!!
      INTEGER            m_AntiQ, m_FoamMode

      COMMON /c_BornV/
     $  m_CMSene,                       ! Initial value of CMS energy
     $  m_XXXene,                       ! CMS energy after beamsstrahlung or beam spread
* -------------------- EVENT --------------------------
     $  m_x1,                           ! 1-z1 = x1 for first  beam(strahlung)
     $  m_x2,                           ! 1-z2 = x2 for second beam(strahlung)
     $  m_vv,                           ! v = 1-sprim/s
     $  m_AvMult,                       ! Average photon multiplicity CRude at given v
     $  m_YFSkon,                       ! YFS formfactor finite part
     $  m_YFS_IR,                       ! YFS formfactor IR part
* -----------------------------------------------------
     $  m_vvmin,                        ! minimum v, infrared cut
     $  m_vvmax,                        ! maximum v
     $  m_HadMin,                       ! minimum hadronization mass [GeV]
* Basic QED and QCD
     $  m_alfinv,                       ! 1/alphaQED, Thomson limit (Q^2=0)
     $  m_alfpi,                        ! alphaQED/pi
     $  m_Xenph,                        ! Enhancement factor for Crude photon multiplicity
* EW parameters
     $  m_MZ,                           ! Z mass
     $  m_amh,                          ! Higgs mass
     $  m_amtop,                        ! Top mass
     $  m_swsq,                         ! sin(thetaW)**2
     $  m_gammz,                        ! Z width
     $  m_MW,                           ! W mass
     $  m_GammW,                        ! W width
     $  m_Gmu,                          ! Fermi constant (from muon decay)
     $  m_alfinvMZ,                     ! 1/alphaQED at (Q^2=MZ^2)     DIZET
     $  m_alfQCDMZ,                     ! alphaQCD   at (Q^2=MZ^2)     DIZET
* Table of fermion paramerets, quarks (1->6) and leptons (11->16)
     $  m_KFferm(20),                   ! fermion KFcode (1->6) and (11->16)
     $  m_NCf(20),                      ! number of colours
     $  m_Qf(20),                       ! electric charge
     $  m_T3f(20),                      ! isospin, L-hand component
     $  m_helic(20),                    ! helicity or polarization
     $  m_amferm(20),                   ! fermion mass
     $  m_auxpar(20),                   ! auxiliary parameter
     $  m_IsGenerated(20),              ! Generation flag, only for SAN !!! 
* Normalization
     $  m_gnanob,                       ! GeV^(-2) to nanobarns
* Initial/final fermion types
     $  m_KFini,                        ! KF code of beam
     $  m_KFfin,                        ! KF code of final state
* Test switches
     $  m_KeyQCD,                       ! QCD FSR corr. switch
     $  m_KeyINT,                       ! ISR/FSR INTereference switch
     $  m_KeyFSR,                       ! ISR/FSR INTereference switch
     $  m_KeyElw,                       ! Type of Electrowak Library
     $  m_KeyZet,                       ! Z-boson on/off
     $  m_KeyWtm,                       ! Photon emission without mass terms
     $  m_KeyRes,                       ! experim. R for gamma* decays switch
     $  m_ibox,                         ! from DZface
     $  m_out,                          ! output unit for printouts
*[[[ NEW!!!
     $  m_FoamMode,                     ! Foam density control <0 for generation
     $  m_AntiQ,                        ! =1 if initial antiquark
     $  m_xpar_input(10000)
      SAVE /c_BornV/
*
* Formats for writing EW tables onto disk file.
      CHARACTER*80  m_fmt0, m_fmt1, m_fmt2
      PARAMETER (
     $  m_fmt0 ='(4g20.13)',                      ! Mz,Mt,Mh etc.
     $  m_fmt1 ='( a,  i4,  f11.5, i4,  f11.5 )', ! header
     $  m_fmt2 ='(6(g13.7,1x))'   )                    ! complex formfactors
*
*  Class procedures:
*
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  BornV                                 //
*//////////////////////////////////////////////////////////////////////////////
