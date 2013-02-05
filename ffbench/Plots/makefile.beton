##############################################################################
##############################################################################
#   Ploting various figures and tables
# 
#	make prod2-dvi		view results from prod2 run (BHLUM4)
#	make obis2-dvi		view results from obis2 run (OLDBIS)
#	make llog2-dvi		view results from llog2 run (LUMLOG)
#
#	make BeMaTabs-dvi	CPC benchmark
#	make BeMaFigs-dvi	CPC benchmark
#
#	make makfig4-dvi	th/96-156 plots, Missing O(alf3)LL
#	make makfig2-dvi	th/95-38 plots, Bhlumi-(Oldbis+Lumlog)
#	make Fig-MisO3-ps	again Missing O(alf3)LL
#
#	This should be invoked from ../ because ../prod3/bhl.hst may not exist
#       make makfig3-dvi	exercise on dummy parameters in MC
##############################################################################
##############################################################################
.SUFFIXES:  .eps .dvi .tex .txi .html .bib .blg .bbl .aux .txp
#----------------------------------------------------------------------
XDVI   = xdvi      -bg white -fg black  -s 4 
XPS    = ghostview -bg white -fg black  -magstep -1 
#----------------------------------------------------------------------
PTL1= ../../lib
EXTER_LIB1  = $(PTL1)/glib.a
#----------------------------------------------------------------------
$(EXTER_LIB1): $(PTL1)/*.f
	cd $(PTL1); make
#----------------------------------------------------------------------
########### collection of hp   flags
# FFLAGS  =   -K -g -C +T
# LDFLAGS =   -K -g -C +T
#--------------------------------------------------------
########### collection of ibm  flags 
# FFLAGS  =   -O   -C -qextchk  -bloadmap:lmap
# LDFLAGS =   -O   -C -qextchk  -bloadmap:lmap
# FFLAGS  =  -C -qextchk -qinitauto=FF -bloadmap:lmap -qflttrap=overflow:underflow:zerodivide:invalid:enable
# LDFLAGS =  -C -qextchk -qinitauto=FF -bloadmap:lmap -qflttrap=overflow:underflow:zerodivide:invalid:enable
########### actual flags
FFLAGS  = -O
LDFLAGS = -O
#----------------------------------------------------------------------
F77 = f77
LD  = $(F77)
#----------------------------------------------------------------------
##############################################################################
.f.o:
	$(F77) $(FFLAGS) -c $<
################## postscript translation (one page) ##################
.txp.dvi:
	(latex $< )
############## option E is for tight bounding on one ##################
.dvi.eps:
	(dvips -E $< -o $@;	$(XPS) $@ )
###	(dvips -E $< -o )
##############################################################################
#----------------------------------------------------------------------
TTL1 = BeMaTabs
objectTTL1 = $(TTL1).o plthis.o readat.o
$(TTL1).exe:			$(objectTTL1) $(EXTER_LIB1)
	$(LD) -o $@ $(LDFLAGS)	$(objectTTL1) $(EXTER_LIB1)
$(TTL1).tex: 	$(TTL1).exe ../prod2/bhl.hst ../obis2/bhl.hst ../llog2/bhl.hst
	(time   $(TTL1).exe )
$(TTL1).dvi:	$(TTL1).tex
	(time   $(TTL1).exe )
	(latex	$(TTL1).tex )
$(TTL1)-dvi:	$(TTL1).dvi
	($(XDVI) $(TTL1).dvi& )
$(TTL1).ps:	$(TTL1).dvi
	(dvips $(TTL1).dvi -o )
$(TTL1)-ps:	$(TTL1).ps
	($(XPS) $(TTL1).ps& )
#----------------------------------------------------------------------
TTL2 = BeMaFigs
objectTTL2 = $(TTL2).o plthis.o readat.o
$(TTL2).exe:			$(objectTTL2) $(EXTER_LIB1)
	$(LD) -o $@ $(LDFLAGS)	$(objectTTL2) $(EXTER_LIB1)
$(TTL2).tex: 	$(TTL2).exe ../prod2/bhl.hst ../obis2/bhl.hst ../llog2/bhl.hst
	(time   $(TTL2).exe )
$(TTL2).dvi:	$(TTL2).tex
	(time   $(TTL2).exe )
	(latex	$(TTL2).tex )
$(TTL2)-dvi:	$(TTL2).dvi
	($(XDVI) $(TTL2).dvi& )
$(TTL2).ps:	$(TTL2).dvi
	(dvips $(TTL2).dvi -o )
$(TTL2)-ps:	$(TTL2).ps
	($(XPS) $(TTL2).ps& )
#----------------------------------------------------------------------
#======================================================================
OBJECTS1 = makfig1.o plthis.o readat.o
makfig1.exe:                   $(OBJECTS1) $(EXTER_LIB1)  
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) $(EXTER_LIB1) 
#----------------------------------------------------------------------
makfig1: makfig1.exe
	(time ./makfig1.exe )
makfig1-dvi: makfig1
	(latex prod1.2M.Z-gamma.tex; $(XDVI)  prod1.2M.Z-gamma)
makfig1-ps: makfig1-dvi
	(dvips prod1.2M.Z-gamma.dvi -o;  $(XPS) prod1.2M.Z-gamma.ps& )
#======================================================================
OBJECTS2 = makfig2.o plthis.o readat.o
makfig2.exe:                   $(OBJECTS2) $(EXTER_LIB1)  
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXTER_LIB1) 
#----------------------------------------------------------------------
makfig2: makfig2.exe ../prod2/bhl.hst ../obis2/bhl.hst ../llog2/bhl.hst
	(time ./makfig2.exe )
makfig2-dvi: makfig2
	(latex fig-95-38; $(XDVI)  fig-95-38)
makfig2-ps: makfig2
	(dvips makfig2.dvi -o;  $(XPS) makfig2.ps& )
#======================================================================
OBJECTS3 = makfig3.o plthis.o readat.o
makfig3.exe:                   $(OBJECTS3) $(EXTER_LIB1)  
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS3) $(EXTER_LIB1) 
#----------------------------------------------------------------------
makfig3: makfig3.exe ../prod3/bhl.hst
	(time ./makfig3.exe )
makfig3-dvi: makfig3
	(latex fig-MC-params; $(XDVI)  fig-MC-params)
makfig3-ps: makfig3
	(dvips makfig3.dvi -o;  $(XPS) makfig3.ps& )
#======================================================================
OBJECTS4 = makfig4.o plthis.o readat.o
makfig4.exe:                   $(OBJECTS4) $(EXTER_LIB1)  
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS4) $(EXTER_LIB1) 
#----------------------------------------------------------------------
makfig4: makfig4.exe ../llog2/bhl.hst
	(time ./makfig4.exe )
makfig4-dvi: makfig4
	(latex fig-MisO3; $(XDVI)  fig-MisO3)
makfig4-ps: makfig4-dvi
	(dvips makfig4.dvi -o;  $(XPS) makfig4.ps& )
#----------------------------------------------------------------------
# Try 'make Fig-MisO3-ps'
#-----------------------
BN4 = MisOgam3
####
$(BN4)-txp-list =  $(BN4)-SICAL.txp $(BN4)-BARE1.txp $(BN4)-CALO2.txp \
	$(BN4)-SICAL2.txp  $(BN4)u-SICAL.txp $(BN4)u-BARE1.txp
$(BN4)-dvi-list = $(BN4)-SICAL.dvi $(BN4)-BARE1.dvi $(BN4)-CALO2.dvi \
	$(BN4)-SICAL2.dvi  $(BN4)u-SICAL.dvi $(BN4)u-BARE1.dvi
$(BN4)-eps-list = $(BN4)-SICAL.eps $(BN4)-BARE1.eps $(BN4)-CALO2.eps \
	$(BN4)-SICAL2.eps  $(BN4)u-SICAL.eps $(BN4)u-BARE1.eps
$($(BN4)-txp-list): makfig4.exe makfig4-dvi
$($(BN4)-dvi-list): 	$($(BN4)-txp-list)
$($(BN4)-eps-list): 	$($(BN4)-dvi-list)
makfig4-eps: $($(BN4)-eps-list)
Fig-MisO3.dvi: makfig4-eps Fig-MisO3.tex
	(latex Fig-MisO3 )
Fig-MisO3-dvi: Fig-MisO3.dvi
	($(XDVI) Fig-MisO3.dvi)
Fig-MisO3.ps: makfig4-eps Fig-MisO3.dvi  
	(dvips Fig-MisO3.dvi -o )
Fig-MisO3-ps: Fig-MisO3.ps
	($(XPS) Fig-MisO3.ps)
#======================================================================
#----------------------------------------------------------------------
Fp2 = prod2
OBJECFp2 = $(Fp2).o plthis.o readat.o
$(Fp2).exe:			$(OBJECFp2) $(EXTER_LIB1)
	$(LD) -o $@ $(LDFLAGS)	$(OBJECFp2) $(EXTER_LIB1)
$(Fp2).tex: 	$(Fp2).exe ../$(Fp2)/bhl.hst
	(time   $(Fp2).exe )
$(Fp2)-dvi: $(Fp2).tex
	(latex	$(Fp2).tex; 	$(XDVI) $(Fp2).dvi)
#----------------------------------------------------------------------
Fo2 = obis2
OBJECFo2 = $(Fo2).o plthis.o readat.o
$(Fo2).exe:			$(OBJECFo2) $(EXTER_LIB1)
	$(LD) -o $@ $(LDFLAGS)	$(OBJECFo2) $(EXTER_LIB1)
$(Fo2).tex: 	$(Fo2).exe ../$(Fo2)/bhl.hst
	(time   $(Fo2).exe )
$(Fo2)-dvi: $(Fo2).tex
	(latex	$(Fo2).tex; 	$(XDVI) $(Fo2).dvi)
#----------------------------------------------------------------------
Fl2 = llog2
OBJECFl2 = $(Fl2).o plthis.o readat.o
$(Fl2).exe:			$(OBJECFl2) $(EXTER_LIB1)
	$(LD) -o $@ $(LDFLAGS)	$(OBJECFl2) $(EXTER_LIB1)
$(Fl2).tex: 	$(Fl2).exe ../$(Fl2)/bhl.hst
	(time   $(Fl2).exe )
$(Fl2)-dvi: $(Fl2).tex
	(latex	$(Fl2).tex; 	$(XDVI) $(Fl2).dvi)
#----------------------------------------------------------------------
#======================================================================
# Miscelaneous, make ibm-to-hp
# Translates double slashes (IBM) to single slashes (HP) in fortran files
ibm-to-hp:
	(mv    ./makfig2.f ./makfig2.f.bak; \
	  	sed -e 's#\\\\#\\#g' ./makfig2.f.bak > ./makfig2.f)
	(mv    ./makfig4.f ./makfig4.f.bak; \
	  	sed -e 's#\\\\#\\#g' ./makfig4.f.bak > ./makfig4.f)
	(mv    ./BeMaTabs.f ./BeMaTabs.f.bak; \
	  	sed -e 's#\\\\#\\#g' ./BeMaTabs.f.bak > ./BeMaTabs.f)
	(mv    ./BeMaFigs.f ./BeMaFigs.f.bak; \
	  	sed -e 's#\\\\#\\#g' ./BeMaFigs.f.bak > ./BeMaFigs.f)
#======================================================================
clean:                                                     
	rm -f *.o *.a *~ *.exe dump.hst *.dvi *.log *.aux *.ps core
Clean:  clean                                   
	cd $(PTL1); make clean


