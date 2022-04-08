#all_includes :=  -I$(top_srcdir)/SRCee $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS)    # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2 -lstdc++ # used while compiling and linking


METASOURCES = AUTO

bin_PROGRAMS = PlotBES Plot1 TabIIIPRD PlotNU PlotTau

############################################
Plot1_SOURCES = HisNorm.cxx Plot1.cxx
Plot1_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot1_LDFLAGS = -L$(ROOTLIBDIR)
############################################
PlotNU_SOURCES = HisNorm.cxx PlotNU.cxx
PlotNU_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotNU_LDFLAGS = -L$(ROOTLIBDIR)
############################################
PlotBES_SOURCES = HisNorm.cxx PlotBES.cxx
PlotBES_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotBES_LDFLAGS = -L$(ROOTLIBDIR)
############################################
TabIIIPRD_SOURCES = HisNorm.cxx TabIIIPRD.cxx
TabIIIPRD_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
TabIIIPRD_LDFLAGS = -L$(ROOTLIBDIR)
############################################
PlotTau_SOURCES = HisNorm.cxx PlotTau.cxx
PlotTau_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotTau_LDFLAGS = -L$(ROOTLIBDIR)

######################### private #############################
SUFFIXES:  .ps .dvi .tex .txp .html .bib .blg .bbl .aux .fig .pdf .eps .jpg
#==============================================================
#==============================================================
.txp.eps:
	(latex $*.txp; dvips -E $*.dvi -o $*.eps;)
.txp.pdf:
	( latex $*.txp; dvips -E $*.dvi -o $*.eps ; ps2pdf -dEPSCrop $*.eps ;)
#	($(FIG2PDF)  $*.fig $*.pdf )
#################
PlotBES-run: PlotBES
	./PlotBES
#################
PlotTau-run: PlotTau
	./PlotTau
#################
Plot1-run: Plot1
	./Plot1
#################
PlotNU-run: PlotNU
	./PlotNU
###################################################
#		make TabIIIPRD-pdf
TabIIIPRD.phony: TabIIIPRD
	(./TabIIIPRD; touch TabIIIPRD.phony;)
####
BM3 = TabIIIPRD
$(BM3).txp:   TabIIIPRD.phony
$(BM3).pdf: $(BM3).txp
#	pdflatex $(BM3).txp
$(BM3)-pdf: klean $(BM3).pdf
	evince  $(BM3).pdf
###################################################
###################################################
#      make Tab2-pdf
Tab2.phony: PlotNU
	(./PlotNU; touch PlotNU.phony;)
####
BM2 = Tab2
$(BM2).txp:   Tab2.phony
$(BM2).pdf: $(BM2).txp
#	pdflatex $(BM2).txp
$(BM2)-pdf: klean $(BM2).pdf
	evince  $(BM2).pdf
###################################################

###################################################
klean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap *.phony)
	(rm -f *.txp chi afb output-* *.hst)
###############################################################
