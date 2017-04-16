INCLUDES = -I$(top_srcdir)/KKsem -I$(top_srcdir)/MaMar $(all_includes) \
	-I$(ROOTINCDIR)
METASOURCES = AUTO

bin_PROGRAMS = Plot1 Plot2 PlotFSR PlotAFB PlotAFB2 PlotIFI

########### KK2f MC libraries ###########
$(top_builddir)/KK2f/libKK2f.a:
	(cd $(top_builddir)/KK2f/; make)
$(top_builddir)/glibk/libGlibk.a:
	(cd $(top_builddir)/glibk/; make)
$(top_builddir)/jetset/libJetset.a:
	(cd $(top_builddir)/jetset/; make)
$(top_builddir)/bornv/libBornvDizet.a:
	(cd $(top_builddir)/bornv/; make)
$(top_builddir)/tauola/libTauola.a:
	(cd $(top_builddir)/tauola/; make)
$(top_builddir)/photos/libPhotos.a:
	(cd $(top_builddir)/photos/; make)
$(top_builddir)/dizet/libDizet.a:
	(cd $(top_builddir)/dizet/; make)
############### extras ##################
$(top_builddir)/KKsem/libKKsem.a:
	(cd $(top_builddir)/KKsem/; make)
$(top_builddir)/bornv/libBornvTabs.a:
	(cd $(top_builddir)/bornv/; make)
#########################################
	
KK2f_Libs2= \
	$(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvTabs.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a

Plot1_SOURCES = KKcol.cxx Plot1.cxx
Plot1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot1_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

PlotAFB_SOURCES = KKcol.cxx PlotAFB.cxx
PlotAFB_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
PlotAFB_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

PlotAFB2_SOURCES = KKcol.cxx PlotAFB2.cxx
PlotAFB2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
PlotAFB2_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot2_SOURCES = Plot2.cxx
Plot2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot2_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

PlotFSR_SOURCES = KKcol.cxx PlotFSR.cxx
PlotFSR_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam
PlotFSR_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

PlotIFI_SOURCES = KKcol.cxx PlotIFI.cxx
PlotIFI_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam
PlotIFI_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

######################### private #############################
SUFFIXES:  .ps .dvi .tex .txp .html .bib .blg .bbl .aux .fig .pdf .eps .jpg
#==============================================================
#===============================================
.txp.eps:
	(latex $*.txp; dvips -E $*.dvi -o $*.eps;)
.txp.pdf:
	( latex $*.txp; dvips -E $*.dvi -o $*.eps ; ps2pdf -dEPSCrop $*.eps ;)
#	($(FIG2PDF)  $*.fig $*.pdf )
#################
plot1: Plot1
	./Plot1
################# PRD88 Table I benchmark ##########
#		make BMarkPRD88-pdf
PlotIFI.phony: PlotIFI
	(./PlotIFI; touch PlotIFI.phony;)
#######
BM1 = BMarkPRD88
FigVprod.pdf:  PlotIFI.phony
FigVprod2.pdf: PlotIFI.phony
Table1.txp:    PlotIFI.phony
#######
$(BM1).pdf: $(BM1).tex Table1.pdf FigVprod.pdf
	pdflatex $(BM1).tex
$(BM1)-pdf: $(BM1).pdf
	evince  $(BM1).pdf
#######
klean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap *.phony)
	(rm -f *.txp chi afb output-* *.hst)
###############################################################
