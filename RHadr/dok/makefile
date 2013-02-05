##############################################################################
# This is a version for presenting pdf file on the screen, use CTR-L
# You may need to adjust page size parameters in seminar.con
#	make LumLEP-pdf
##############################################################################
.SUFFIXES:  .eps .dvi .tex .txi .html .bib .blg .bbl .aux .txp
#---------------------------------------------------------------------------
XDVI    = xdvi      -bg white -fg black  -s 4 
XPSn    = gv -seascape -geometry +10+0 -magstep  0 -bg white -fg black
XPS     = $(XPSn)
XPdf    = acroread -geometry 800x675+10+0 
# print
XPRINT    = xprint -d 4-2-11
REDU      ='1:0@0.925(10mm,-1mm)'
#REDU      ='1:0@0.925(10mm,10mm)'
#---------------------------------------------------------------------------
#==================================================================================
#	make MBS-ps
TeX1 = MBS
##############  Creation of postscript, printing etc.
$(TeX1).dvi:	$(TeX1).tex seminar.con
	(latex  $(TeX1).tex)
$(TeX1)-dvi:	$(TeX1).dvi
	( $(XDVI) $(TeX1).dvi& )
$(TeX1).ps:	$(TeX1).dvi
	(dvips  $(TeX1).dvi -o )
$(TeX1)-ps:	$(TeX1).ps
	($(XPS) $(TeX1).ps& )
$(TeX1).pdf:	$(TeX1).ps
	(ps2pdf $(TeX1).ps $(TeX1).pdf  )
$(TeX1)-pdf:	$(TeX1).pdf
	($(XPdf) $(TeX1).pdf&)
#==================================================================================
#	make KKres-ps
TeX2 = KKres
##############  Creation of postscript, printing etc.
$(TeX2).dvi:	$(TeX2).tex seminar.con
	(latex  $(TeX2).tex)
$(TeX2)-dvi:	$(TeX2).dvi
	( $(XDVI) $(TeX2).dvi& )
$(TeX2).ps:	$(TeX2).dvi
	(dvips  $(TeX2).dvi -o )
$(TeX2)-ps:	$(TeX2).ps
	($(XPS) $(TeX2).ps& )
$(TeX2).pdf:	$(TeX2).ps
	(ps2pdf $(TeX2).ps $(TeX2).pdf  )
$(TeX2)-pdf:	$(TeX2).pdf
	($(XPdf) $(TeX2).pdf&)
#==================================================================================
clean:
	(rm -f *.o *.exe *~ *.dvi *.log *.blg *.bbl *.phony  *.aux dump.hst core lmap)
#==================================================================================
gtar:	clean
	(cd ..; gtar -cvzf LumLEP_2003feb27.tar.gz LumLEP)
##############################################################################

