INCLUDES = 
METASOURCES = AUTO
bin_PROGRAMS = figinc figbeast

figinc_SOURCES = figinc.f readat.f
figbeast_SOURCES = figbeast.f readat.f

figbeast_LDADD = $(top_builddir)/KK2f/libKK2f.a\
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/glibk/libGlibk.a\
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a\
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a

figinc_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/dizet/libDizet.a

#----------------------------------------------------------------------
figinc-run: figinc
	./figinc
figinc-pdf: figinc-run
	(pdflatex figinc.tex; kpdf figinc.pdf&)
#----------------------------------------------------------------------
