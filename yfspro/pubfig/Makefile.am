INCLUDES = 
METASOURCES = AUTO

bin_PROGRAMS = chi_prepare chi_mcan afb_prepare afb_int delmu delmu_plot delhad bornex

chi_prepare_SOURCES = chi_prepare.f chi.h

chi_mcan_SOURCES = chi_mcan.f chi.h

afb_prepare_SOURCES = afb_prepare.f afb.h

afb_int_SOURCES = afb_int.f chi.h

delmu_SOURCES = delmu.f

delmu_plot_SOURCES = delmu_plot.f

delhad_SOURCES = delhad.f

bornex_SOURCES = bornex.f

bornex_LDADD = $(top_builddir)/KKsem/libKKsem.a $(top_builddir)/KK2f/libKK2f.a\
	$(top_builddir)/glibk/libGlibk.a $(top_builddir)/jetset/libJetset.a\
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a\
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a

chi_prepare_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a \
	$(top_builddir)/glibk/libGlibk.a

afb_prepare_LDADD = $(top_builddir)/KKsem/libKKsem.a \
	$(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a $(top_builddir)/jetset/libJetset.a \
	$(top_builddir)/bornv/libBornvTabs.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a

chi_mcan_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a\
	$(top_builddir)/glibk/libGlibk.a

afb_int_LDADD = $(top_builddir)/KKsem/libKKsem.a $(top_builddir)/KK2f/libKK2f.a\
	$(top_builddir)/glibk/libGlibk.a $(top_builddir)/jetset/libJetset.a\
	$(top_builddir)/bornv/libBornvTabs.a $(top_builddir)/tauola/libTauola.a\
	$(top_builddir)/photos/libPhotos.a

delhad_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a\
	$(top_builddir)/glibk/libGlibk.a

delmu_plot_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a\
	$(top_builddir)/glibk/libGlibk.a

delmu_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a\
	$(top_builddir)/glibk/libGlibk.a
