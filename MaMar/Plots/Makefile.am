INCLUDES = -I$(top_srcdir)/KKsem  \
    -I$(top_srcdir)/MaMar $(all_includes) \
	-I$(ROOTINCDIR)
METASOURCES = AUTO

bin_PROGRAMS = Plot1 Plot2 PlotZinv PlotFSR \
   PlotAFB PlotAFB2 PlotIFI PlotFoam0 PlotFoam1 PlotFoam2\
   PlotBorn PlotKKsem PlotPlay PlotALF

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
############################################
#  Tables (re)generated in flight.
#  To be xchecked if it works correctly!!!
############################################
KK2f_Libs2= \
	$(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a
############################################
#  From tables on the disk (original)
############################################
KK2f_Libs= \
	$(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvTabs.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a

Plot1_SOURCES = KKplot.cxx HisNorm.cxx Plot1.cxx
Plot1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot1_LDFLAGS = -L$(ROOTLIBDIR)

Plot2_SOURCES =  HisNorm.cxx Plot2.cxx
Plot2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot2_LDFLAGS = -L$(ROOTLIBDIR)

PlotZinv_SOURCES =  HisNorm.cxx PlotZinv.cxx
PlotZinv_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotZinv_LDFLAGS = -L$(ROOTLIBDIR)

PlotAFB_SOURCES = KKplot.cxx HisNorm.cxx PlotAFB.cxx
PlotAFB_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotAFB_LDFLAGS = -L$(ROOTLIBDIR)

PlotAFB2_SOURCES = KKplot.cxx HisNorm.cxx PlotAFB2.cxx
PlotAFB2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotAFB2_LDFLAGS = -L$(ROOTLIBDIR)

PlotFSR_SOURCES = KKplot.cxx HisNorm.cxx PlotFSR.cxx
PlotFSR_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam -lgfortran
PlotFSR_LDFLAGS = -L$(ROOTLIBDIR)

PlotIFI_SOURCES = KKplot.cxx HisNorm.cxx PlotIFI.cxx
PlotIFI_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotIFI_LDFLAGS = -L$(ROOTLIBDIR)

PlotFoam0_SOURCES = KKplot.cxx HisNorm.cxx HisReMake.cxx PlotFoam0.cxx
PlotFoam0_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam -lgfortran
PlotFoam0_LDFLAGS = -L$(ROOTLIBDIR)

PlotFoam1_SOURCES = KKplot.cxx HisNorm.cxx HisReMake.cxx PlotFoam1.cxx
PlotFoam1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam -lgfortran
PlotFoam1_LDFLAGS = -L$(ROOTLIBDIR)

PlotFoam2_SOURCES = KKplot.cxx HisNorm.cxx HisReMake.cxx PlotFoam2.cxx
PlotFoam2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam -lgfortran
PlotFoam2_LDFLAGS = -L$(ROOTLIBDIR)

PlotBorn_SOURCES = KKplot.cxx HisNorm.cxx PlotBorn.cxx
PlotBorn_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotBorn_LDFLAGS = -L$(ROOTLIBDIR) -no-install -g

PlotKKsem_SOURCES = KKplot.cxx HisNorm.cxx HisReMake.cxx PlotKKsem.cxx
PlotKKsem_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotKKsem_LDFLAGS = -L$(ROOTLIBDIR)

PlotPlay_SOURCES =  PlotPlay.cxx
PlotPlay_LDADD   =  $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotPlay_LDFLAGS = -L$(ROOTLIBDIR)

PlotALF_SOURCES = KKplot.cxx HisNorm.cxx PlotALF.cxx
PlotALF_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
PlotALF_LDFLAGS = -L$(ROOTLIBDIR)


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
plot1: Plot1
	./Plot1
#################
Plot1-run: Plot1
	./Plot1
PlotZinv-run: PlotZinv
	./PlotZinv
PlotIFI-run: PlotIFI
	./PlotIFI
PlotAFB-run: PlotAFB
	./PlotAFB
PlotAFB2-run: PlotAFB2
	./PlotAFB2
PlotFoam0-run: PlotFoam0
	./PlotFoam0
PlotFoam1-run: PlotFoam1
	./PlotFoam1
PlotFoam2-run: PlotFoam2
	./PlotFoam2
PlotBorn-run: PlotBorn
	./PlotBorn
PlotKKsem-run: PlotKKsem
	./PlotKKsem
PlotALF-run: PlotALF
	./PlotALF
################# PRD88 Table I benchmark ##########
#		make BMarkPRD88-pdf
PlotIFI.phony: PlotIFI
	(./PlotIFI; touch PlotIFI.phony;)
####### make BMarkPRD88-pdf
BM1 = BMarkPRD88
cFigOldBench.pdf:  PlotIFI.phony
cFigOldBench2.pdf: PlotIFI.phony
TabOldBench.txp:   PlotIFI.phony
#######
$(BM1).pdf: $(BM1).tex TabOldBench.pdf cFigOldBench.pdf
	pdflatex $(BM1).tex
$(BM1)-pdf: $(BM1).pdf
	evince  $(BM1).pdf
#######
klean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap *.phony)
	(rm -f *.txp chi afb output-* *.hst)
###############################################################
