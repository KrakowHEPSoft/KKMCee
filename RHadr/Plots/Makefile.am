INCLUDES = -I$(top_srcdir)/RHadr $(all_includes) -I$(ROOTINCDIR)

METASOURCES = AUTO
bin_PROGRAMS = SkyLine Plot2b Plot2c Plot3a Plot3b Plot3m Plot3f

SkyLine_SOURCES = SkyLine.cxx
Plot2b_SOURCES  = Plot2b.cxx
Plot2c_SOURCES  = Plot2c.cxx
Plot3a_SOURCES  = Plot3a.cxx
Plot3b_SOURCES  = Plot3b.cxx
Plot3m_SOURCES  = Plot3m.cxx
Plot3f_SOURCES  = Plot3f.cxx

SkyLine_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(ROOTAUXLIBS) $(ROOTGLIBS)
SkyLine_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot2b_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot2b_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot2c_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot2c_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot3a_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot3a_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot3b_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot3b_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot3m_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot3m_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot3f_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot3f_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran


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

