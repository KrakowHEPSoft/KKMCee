INCLUDES = -I$(top_srcdir)/KKsem -I$(top_srcdir)/MaMar $(all_includes) \
	-I$(ROOTINCDIR)
METASOURCES = AUTO
bin_PROGRAMS = Plot1
Plot1_SOURCES = KKsem.cxx Plot1.cxx
Plot1_LDADD = $(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvTabs.a $(top_builddir)/tauola/libTauola.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a $(top_builddir)/KKsem/libKKsem.a \
	$(ROOTAUXLIBS) $(ROOTGLIBS)
Plot1_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran
