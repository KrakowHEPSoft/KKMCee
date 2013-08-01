
METASOURCES = AUTO

bin_PROGRAMS = chi_prepare chi_mcan afb_prepare afb_int delmu delmu_plot delhad bornex

chi_prepare_SOURCES = chi_prepare.f chi.h

chi_mcan_SOURCES    = chi_mcan.f chi.h

afb_prepare_SOURCES = afb_prepare.f afb.h

afb_int_SOURCES     = afb_int.f chi.h

delmu_SOURCES       = delmu.f

delmu_plot_SOURCES  = delmu_plot.f

delhad_SOURCES      = delhad.f

bornex_SOURCES      = bornex.f

bornex_LDADD = $(top_builddir)/KKsem/libKKsem.a $(top_builddir)/KK2f/libKK2f.a\
	$(top_builddir)/glibk/libGlibk.a $(top_builddir)/jetset/libJetset.a\
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/tauola/libTauola.a\
	$(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a

chi_prepare_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a \
	$(top_builddir)/glibk/libGlibk.a

afb_prepare_LDADD = $(top_builddir)/KKsem/libKKsem.a \
	$(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a $(top_builddir)/jetset/libJetset.a \
	$(top_builddir)/bornv/libBornvTabs.a $(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a

chi_mcan_LDADD = $(top_builddir)/yfspro/libfig/libSemaPlot.a \
	$(top_builddir)/glibk/libGlibk.a $(top_builddir)/KKsem/libKKsem.a

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

#######################################################################################################
###                    Private part of makefile
#######################################################################################################
#  Plots of Sigma and AFB
#     make afb_sig-ps  (sometimes twice!)
#     make afb_int-ps  (sometimes twice!)
#     make afb_sig-ps  afb_sig-ps with smaller table
#     make afb_ang-ps
#  other plots
#     make delta-ps
#  comparisons with semianalyticals
#     make chi_mcan-ps  (sometimes twice!)
############# postscript translation (one page) #################
XPS     = kghostview
SUFFIXES:  .eps .dvi .tex .txi .html .bib .blg .bbl .aux .txp
.txp.dvi:
	(latex $< )
#----------------------------------------------
.txp.eps:
	(latex $*.txp;    dvips -E  $*.dvi -o $*.eps)
#======================================================================================================
chi.hst:	chi_prepare
	(time ./chi_prepare )
BN0 = chi_mcan
CHI_EPS_LIST = $(BN0)-O2mca.eps  $(BN0)-O2dif.eps $(BN0)-O3dan.eps $(BN0)-O0dif.eps $(BN0)-O0tech.eps $(BN0)-O2mO1.eps $(BN0)-O3dO2.eps
chi_mcan-run: chi_mcan chi.hst
	(time ./chi_mcan; touch chi_mcan-run)
chi_mcan-eps:  $(CHI_EPS_LIST)
	touch chi_mcan-eps;
#--------------------------
$(BN0).dvi: chi.hst chi_mcan-run chi_mcan-eps	$(BN0).tex
	(latex   $(BN0).tex )
$(BN0)-dvi: 		$(BN0).dvi
	($(XDVI) $(BN0).dvi)
$(BN0).ps: 		$(BN0).dvi
	(dvips $(BN0).dvi -o )
$(BN0)-ps: 		$(BN0).ps
	($(XPS) $(BN0).ps)
#======================================================================================================
afb.hst:	afb_prepare
	(time ./afb_prepare )
#-----------------------------------------------------------------------------
BN4 = afb_int
AFB_EPS_LIST= $(BN4)-tab1.eps $(BN4)-tab2.eps $(BN4)-tab3.eps $(BN4)-G1.eps $(BN4)-G1x.eps $(BN4)-G1xxx.eps $(BN4)-afb1.eps  $(BN4)-afb2.eps $(BN4)-sig1.eps $(BN4)-Gsig.eps $(BN4)-Gsig.eps $(BN4)-GsigZF.eps $(BN4)-Gafb.eps $(BN4)-GafbZF.eps $(BN4)-afb1S.eps $(BN4)-sig1S.eps $(BN4)-afbHO.eps $(BN4)-sigHO.eps $(BN4)-com1.eps $(BN4)-com1x.eps $(BN4)-com1xxx.eps $(BN4)-AngMx.eps $(BN4)-comMx.eps $(BN4)-tabEWG1.eps $(BN4)-tabEWG2.eps
#-----------------------------------------------------------------------------
afb_int-run: afb_int afb.hst
	(time ./afb_int; touch afb_int-run)
afb_int-eps:  $(AFB_EPS_LIST)
	touch afb_int-eps;
#-----------------------------------------------------------------------------
$(BN4).dvi:	      afb.hst afb_int-run afb_int-eps	$(BN4).tex
	(latex   $(BN4).tex )
$(BN4).ps: 			$(BN4).dvi
	(dvips $(BN4).dvi -o )
$(BN4)-ps: 			$(BN4).ps
	($(XPS) $(BN4).ps&)
#-----------------------------------------------------------------------------
AFB_SIG = afb_sig
$(AFB_SIG).dvi:   afb.hst afb_int-run afb_int-eps $(AFB_SIG).tex
	(latex   $(AFB_SIG).tex )
$(AFB_SIG).ps: 		$(AFB_SIG).dvi
	(dvips $(AFB_SIG).dvi -o )
$(AFB_SIG)-ps: 		$(AFB_SIG).ps
	($(XPS) $(AFB_SIG).ps&)
#-----------------------------------------------------------------------------
AFB_SIG2 = afb_sig2
$(AFB_SIG2).dvi:  afb.hst afb_int-run afb_int-eps	 $(AFB_SIG2).tex
	(latex   $(AFB_SIG2).tex )
$(AFB_SIG2).ps: 		$(AFB_SIG2).dvi
	(dvips $(AFB_SIG2).dvi -o )
$(AFB_SIG2)-ps: 		$(AFB_SIG2).ps
	($(XPS) $(AFB_SIG2).ps&)
#-----------------------------------------------------------------------------
AFB_ANG = afb_ang
$(AFB_ANG).dvi:   afb.hst afb_int-run afb_int-eps	 $(AFB_ANG).tex
	(latex   $(AFB_ANG).tex )
$(AFB_ANG).ps: 		$(AFB_ANG).dvi
	(dvips $(AFB_ANG).dvi -o )
$(AFB_ANG)-ps: 		$(AFB_ANG).ps
	($(XPS) $(AFB_ANG).ps&)
#-----------------------------------------------------------------------------
#  make delta-ps
delmu.data:	delmu
	(time ./delmu )
	(more  delmu.data)
delhad.data:	delhad
	(time ./delhad)
	(more  delhad.data)
#--------
delmu-sigma1.txp delmu-sigma2.txp:	delmu.data delhad.data delmu_plot
	(time ./delmu_plot )
delmu-sigma1.eps:		delmu-sigma1.txp
delmu-sigma2.eps:		delmu-sigma2.txp
#-------
delta.dvi:		delta.tex delmu-sigma1.eps delmu-sigma2.eps
	(latex   delta.tex )
delta.ps: 		delta.dvi
	(dvips delta.dvi -o )
delta-ps: 		delta.ps
	($(XPS) delta.ps&)
#-----------------------------------------------------------------------------
myclean: clean
	(rm -f $(AFB_EPS_LIST) $(CHI_EPS_LIST) )
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap)
	(rm -f *.txp chi afb output-* *.hst *.data)
#-----------------------------------------------------------------------------
AM_FFLAGS = @MY_F77_FLAGS@
