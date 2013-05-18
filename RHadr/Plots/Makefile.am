INCLUDES = 
METASOURCES = AUTO
bin_PROGRAMS = SkyLine Plot2b

SkyLine_SOURCES = SkyLine.f
Plot2b_SOURCES  = Plot2b.f

SkyLine_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a

Plot2b_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a

#----------------------------------------------------------------------
Plot2b-run: Plot2b
	./Plot2b
Plot2b-pdf: Plot2b-run
	(pdflatex Plot2b.tex; kpdf Plot2b.pdf&)
#----------------------------------------------------------------------
