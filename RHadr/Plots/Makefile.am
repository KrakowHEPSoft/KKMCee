AM_CPPFLAGS = -I$(top_srcdir)/RHadr $(all_includes) -I$(ROOTINCDIR)

AM_CXXFLAGS = -g -O2 -std=c++11 -lstdc++ 

METASOURCES = AUTO
bin_PROGRAMS = SkyLine Plot2b Plot2c Plot3a Plot3b Plot3m Plot3f

SkyLine_SOURCES = SkyLine.cxx
Plot2b_SOURCES  = Plot2b.cxx
Plot2c_SOURCES  = Plot2c.cxx
Plot3a_SOURCES  = Plot3a.cxx
Plot3b_SOURCES  = Plot3b.cxx
Plot3m_SOURCES  = Plot3m.cxx
Plot3f_SOURCES  = Plot3f.cxx

KK2f_Libs2= \
	$(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a  $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a  -lstdc++



SkyLine_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
SkyLine_LDFLAGS = -L$(ROOTLIBDIR)

Plot2b_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot2b_LDFLAGS  = -L$(ROOTLIBDIR)

Plot2c_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot2c_LDFLAGS  = -L$(ROOTLIBDIR)

Plot3b_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot3b_LDFLAGS  = -L$(ROOTLIBDIR)

Plot3a_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot3a_LDFLAGS  = -L$(ROOTLIBDIR)

Plot3m_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot3m_LDFLAGS  = -L$(ROOTLIBDIR)

Plot3f_LDADD =  $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot3f_LDFLAGS  = -L$(ROOTLIBDIR)


#=============================================================================
#                  User part
#-----------------------------------------------------------------------------
# The examples of plotings
#   (make SkyLine; ./Skyline )    this one will plot R(s)
#   (make Plot2b;  ./Plot2b )
#   (make Plot2c;  ./Plot2c )
#   etc...
#-----------------------------------------------------------------------------
# The other ploting ROOT scripts are activated as follows:
#   root Fig0.C      xcheck of normalization procedure
#   root Fig1.C      xsection(nphot), dsig/d(theta_gamma)
#   root Fig1b.C     Pions, pT and CoTheta 
#   root Fig2.C      Q2(pions) with various cuts
#   root Fig2b.C     Q2 distr. 2pi versus muons
#   root Fig2c.C     Q2 distr. NEW
#   root Fig3.C      Higer order effects		
#-----------------------------------------------------------------------------
# The input root files with histograms are in ../demoC/
# Note that variable RPATH for Root may have to be readjusted
#=============================================================================

