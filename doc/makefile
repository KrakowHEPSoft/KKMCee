# make ICHEP2010.pdf
#==============================================================================
.SUFFIXES:  .ps .dvi .tex .txi .html .bib .blg .bbl .aux .fig .pdf .eps
#----------------------------------------------------------------------
FIG2PDF = fig2dev -L pdf -p 0
XPdf    = xpdf
#===============================================
.fig.pdf:
	($(FIG2PDF)  $*.fig $*.pdf )
#==============================================================================
all:	IFJtalk13.pdf
#==============================================================================
# List of pdf files exported from fig files
GRAFLIST = 

CLEAN_GRAFLIST:
	(rm -f $(GRAFLIST) )
#==============================================================================
fig_count:
	(cat $(TeX2).tex | grep includegraphics )
	(cat $(TeX2).tex | grep includegraphics | wc -l )
	(ls -altr | grep pdf | wc -l )
	(ls -altr | grep jpg | wc -l )
#==============================================================================
#==============================================================================
TeX2 =  IFJtalk13
$(TeX2).pdf: 	$(TeX2).tex $(GRAFLIST)
	(pdflatex	$(TeX2) )
$(TeX2)-pdf:	$(TeX2).pdf
	($(XPdf) $(TeX2).pdf& )
#=============================================================================
clean:
	(rm -f *.dvi *~ *.log  *.aux *_region_* core *.bbl *.blg  *.phony)
	(rm -f *.toc *.nav *.out *.snm *.bak)
#=============================================================================


