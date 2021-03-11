all_includes :=  -I$(top_srcdir)/SRCee $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS)    # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2 -lstdc++ # used while compiling and linking


METASOURCES = AUTO

bin_PROGRAMS = Plot1 

############################################
Plot1_SOURCES = HisNorm.cxx Plot1.cxx
Plot1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
Plot1_LDFLAGS = -L$(ROOTLIBDIR)

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
Plot1-run: Plot1
	./Plot1

###################################################
klean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap *.phony)
	(rm -f *.txp chi afb output-* *.hst)
###############################################################
