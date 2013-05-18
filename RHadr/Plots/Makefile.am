INCLUDES = -I$(top_srcdir)/RHadr $(all_includes) -I$(ROOTINCDIR)

METASOURCES = AUTO
bin_PROGRAMS = SkyLine Plot2b

SkyLine_SOURCES = SkyLine.cxx
Plot2b_SOURCES  = Plot2b.cxx

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

#----------------------------------------------------------------------
Plot2b-run: Plot2b
	./Plot2b
Plot2b-pdf: Plot2b-run
	(pdflatex Plot2b.tex; kpdf Plot2b.pdf&)
#----------------------------------------------------------------------

