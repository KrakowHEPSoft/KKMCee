#======================================================================
#                    Makefile for everything
# Notes:
#	make makflag		<- upgrate of f77 compilation flags everywhere
#	source ./RootLogin	<- to enter Root enviroment
#	make rmain.exe		<- compile and link
#       make makprod		<- update production makefiles
#
#	Example of production cycle
#	make 172GeV.4J-start	<- production...
#	make 172GeV.4J-stop
#	make 172GeV.4J-cont
#	make 172GeV.4J-debug
#	(cd 172GeV.4J; make farm4-SP)		<- make subdirectory for batch
#       (cd 172GeV.4J/1/; llsubmit 172GeV.4J.cmd.1)	<- submit batch job
#
#	make fig                <-- view 8 figures from BE paper
#======================================================================
all: rmain.exe
#######################################################################
### version 0.8
### ROOTdir = $(ROOTDIR)
# version 0.9
ROOTdir = $(ROOTSYS)
#----------------------------------------------------------------------
.SUFFIXES:  .o .c .C .f .l  .s  .sh  .h  .a
#----------------------------------------------------------------------
#======================================================================
#        FORTRAN Compile Link Parameters
# fortran is not realy used in this directory but compilation
# flags are copied (with make makflag) to other directiories for debug
#----------------------------------------------------------------------
#--- collection of HP flags ----
# FFLAGS  =   -K -g -C +T
# LDFLAGS =   -K -g -C +T
#-------------------------------
F77  = xlf
FLD  = xlf
##FFLAGS   =  -g -qextname
##FLDFLAGS =  -g 
FFLAGS   =  -O -qextname
FLDFLAGS =  -O
#======================================================================
#  C++ Compile Link Parameters
#-----------------------------
#####  IBM AIX
CXX       = xlC
CXXFLAGS  = -g -qcheck=all -w -qnoro -qnoroconst -qmaxmem=-1 -I$(ROOTdir)/include
CLD       = xlC
CLDFLAGS  = -g -bloadmap:lmap
SLD       = /usr/lpp/xlC/bin/makeC++SharedLib -bloadmap:lmap
###CLDFLAGS  = -g -bloadmap:lmap -qcheck=all
#-----------------------------
#####  GNU c++ compiler
##CXX         =  g++
##CXXFLAGSS   =  -g
##CLD         =  g++
##CLDFLAGS    =  -g -Xlinker -bloadmap:lmap
##############
#======================================================================
#  Root 
ROOTLIBS      = -L$(ROOTdir)/lib -lRoot -lCint
LIBS          = $(ROOTLIBS) -lm -lxlf -lxlf90 
LIBSm         = $(ROOTLIBS) -lm
#======================================================================
# Compile all c++ files
#-----------------------------
.C.o:
	$(CXX) $(CXXFLAGS) -c $<
#-----------------------------
# Compile all fortran files
#-----------------------------
.f.o:
	$(F77) $(FFLAGS) -c $<
#======================================================================


#======================================================================
#       External Fortran Libraries
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#  Cern
CERNDIR       = /cern/pro
LIBCERN       = -L$(CERNDIR)/lib  -lpacklib -lpawlib -ljetset74
#======================================================================
#              Private libraries KORALW
#----------------------------------------------------------------------
EXTLIB0   = ../../korww/korww.a
$(EXTLIB0): ../../korww/*.f
	(cd ../../korww; make)
LIB0 = -L../../korww -lkorww
#---
EXTLIB1   = ../../glib/glib.a
$(EXTLIB1): ../../glib/*.f
	(cd ../../glib; make)
LIB1 = -L../../glib -lglib
#---
EXTLIB2   = ../../kwlund/kwlund.a
$(EXTLIB2): ../../kwlund/*.f
	(cd ../../kwlund; make)
LIB2 = -L../../kwlund -lkwlund
#---
EXTLIB3   = ../../model/model.a
$(EXTLIB3): ../../model/*.f
	(cd ../..//model; make)
LIB3 = -L../../model -lmodel
#---
EXTLIB4   = ../../tauola/lib/taulib.a
$(EXTLIB4): ../../tauola/lib/*.f
	(cd ../../tauola/lib; make)
LIB4 = -L../../tauola/lib -ltaulib
#---
EXTLIB5   = ../../interfaces/interface.a
$(EXTLIB5): ../../interfaces/*.f
	(cd ../../interfaces; make)
LIB5 = -L../../interfaces -linterface
#---
EXTLIB6   = ../../semian/wwanal.a
$(EXTLIB6): ../../semian/*.f
	(cd ../../semian; make)
LIB6 = -L../../semian -lwwanal
#---
EXTLIB7   = ../../ampli4f/ampli4f.a
$(EXTLIB7): ../../ampli4f/*.f
	(cd ../../ampli4f; make)
LIB7 = -L../../ampli4f -lampli4f
#---
EXTLIBS  = $(EXTLIB0) $(EXTLIB1) $(EXTLIB2) $(EXTLIB3) \
           $(EXTLIB4) $(EXTLIB5) $(EXTLIB6) $(EXTLIB7)
LIBSY = $(LIB1) $(LIB2) $(LIB3) $(LIB4) $(LIB5) $(LIB6) $(LIB7)
#======================================================================

ROBOL.o:	ROBOL.h
Semaph.o: 	Semaph.h
JetAnalyzer.o: 	JetAnalyzer.h
BEwtMaker.o: 	BEwtMaker.h
KoralwMaker.o:	KoralwMaker.h
KorEvent.o: 	KorEvent.h
PartLund.o: 	PartLund.h
VLorenz.o: 	VLorenz.h

###KORALWh = VLorenz.h PartLund.h KorEvent.h KoralwMaker.h
###KORALWo = VLorenz.o PartLund.o KorEvent.o KoralwMaker.o
KORALWh = VLorenz.h PartLund.h KorEvent.h KoralwMaker.h \
	BEwtMaker.h JetAnalyzer.h Semaph.h ROBOL.h
KORALWo = VLorenz.o PartLund.o KorEvent.o KoralwMaker.o \
	BEwtMaker.o JetAnalyzer.o Semaph.o ROBOL.o
#//////////////////////////////////////////////////////////////////////////////
#   creation of dictionary to be compiled and liked with the main program
#			make libKorw.so
#			make cleanso
KorwDict.C: $(KORALWh)
	 $(ROOTSYS)/bin/rootcint -f KorwDict.C -c $(KORALWh) LinkDef.h
#/////// here we also produce shared library for use under cint interpreter
libKorw.so: $(KORALWo) KorwDict.o
	$(SLD) -o libKorw.so $(ROOTLIBS) -lm -lxlf -lxlf90 -lcurses \
	-p 0 $(KORALWo) \
	partit1.o  ReaData.o  $(EXTLIBS) \
	KorwDict.o
	@echo " >>>>>>> libKorw.so done"
cleanso: clean
	rm -f libKorw.so KorwDict.*
#//////////////////////////////////////////////////////////////////////////////
#======================================================================
MAIN1 = mainVLorenz
####OBJ1  = $(MAIN1).o
####$(MAIN1).o: $(MAIN1).C
OBJ1  = $(MAIN1).o VLorenz.o
$(MAIN1).o: $(MAIN1).C VLorenz.h
$(MAIN1).exe:	$(OBJ1)  $(EXTLIBS)
		$(CLD) $(CLDFLAGS) $(OBJ1) $(LIBSm) -o $(MAIN1).exe
$(MAIN1): $(MAIN1).exe
	(time ./$(MAIN1).exe)
debug:	$(MAIN1).exe
	( xldb ./$(MAIN1).exe)
#======================================================================
MAIN3 = rmain
rmain.o: rmain.C $(KORALWh) partit1.h
OBJ3  = $(MAIN3).o partit1.o ReaData.o\
	$(KORALWo) KorwDict.o       #### WITH    dictionary
####	$(KORALWo)                  #### WITHOUT dictionary
$(MAIN3).exe:	$(OBJ3)  $(EXTLIBS)
		$(CLD) $(CLDFLAGS) $(OBJ3) $(LIBS) -o $(MAIN3).exe $(EXTLIBS)
#======================================================================
#-----------------------------------------------------------------------
#                   MC Production
#-----------------------------------------------------------------------
172GeV.4J.mini-start:	$(MAIN3).exe
	(cd ../run; make 172GeV.4J.mini-start)
172GeV.4J.mini-stop:
	(cd ../run; make 172GeV.4J.mini-stop)
172GeV.4J.mini-debug:	$(MAIN3).exe
	(cd ../run; make 172GeV.4J.mini-debug)
#-----------------------------------------------------------------------
172GeV.2J.mini-start:	$(MAIN3).exe
	(cd ../run; make 172GeV.2J.mini-start)
#-----------------------------------------------------------------------
172GeV.4J-start:	$(MAIN3).exe
	(cd ../run; make 172GeV.4J-start)
172GeV.4J-debug:	$(MAIN3).exe
	(cd ../run; make 172GeV.4J-debug)
172GeV.4J-stop:
	(cd ../run; make 172GeV.4J-stop)
172GeV.4J-cont:	$(MAIN3).exe
	(cd ../run; make 172GeV.4J-cont)
172GeV.4J-diff:
	(cd ../run; make 172GeV.4J;  diff koralw.output.100ev koralw.output)
#-----------------------------------------------------------------------
200GeV.4J-start:	$(MAIN3).exe
	(cd ../run; make 200GeV.4J-start)
200GeV.4J-stop:	$(MAIN3).exe
	(cd ../run; make 200GeV.4J-stop)
#-----------------------------------------------------------------------
200GeV.4J.mini-start:	$(MAIN3).exe
	(cd ../run; make 200GeV.4J.mini-start)
200GeV.4J.mini-stop:	$(MAIN3).exe
	(cd ../run; make 200GeV.4J.mini-stop)
#-----------------------------------------------------------------------
200GeV.2J.mini-start:	$(MAIN3).exe
	(cd ../run; make 200GeV.2J.mini-start)
#=============================================================================
# Update fortran compilation flags, adding lib prefix
makflag2:
	(cd ../../korww;   sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#korww.a#libkorww.a#'	makefile.templ > makefile )
	(cd ../../demo;    sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		makefile.templ > makefile )
	(cd ../../model;   sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#model.a#libmodel.a#'	makefile.templ > makefile )
	(cd ../../glib/;   sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#glib.a#libglib.a#'	makefile.templ > makefile )
	(cd ../../kwlund;  sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#kwlund.a#libkwlund.a#'	makefile.templ > makefile )
	(cd ../../tauola/lib; sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#taulib.a#libtaulib.a#'	makefile.templ > makefile )
	(cd ../../interfaces; sed -e 's#<FFLAGS>#$(FFLAGS)#' \
	    -e 's#interface.a#libinterface.a#'	makefile.templ > makefile )
	(cd ../../ampli4f; sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#ampli4f.a#libampli4f.a#'	makefile.templ > makefile )
	(cd ../../semian;  sed -e 's#<FFLAGS>#$(FFLAGS)#' \
		-e 's#wwanal.a#libwwanal.a#'	makefile.templ > makefile )
#=============================================================================
# Update fortran compilation flags, adding lib prefix
makflag:
	(cd ../../korww; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../demo;  sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../model; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../glib/; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../kwlund; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../tauola/lib; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../interfaces; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../ampli4f; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../../semian;  sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
#==============================================================================
makprod:
	(cd ../run; make makprod)
#==============================================================================
fig:
	(cd ../fig; make menu)
#==============================================================================
clean:
	(rm -f *.o *.a *~ *.exe core *.gz lmap *.out *.dvi *.aux *.log)
Clean:  clean
	(cd ../../;  make Clean)
	(rm -r ./dok)
CLEAN:  Clean
	(cd ../../;  make CLEAN)
#------------------------------------------------------------------------------
DIRKW = koralw-1.21-export
gtar:	clean
	(cd ../../../; \
	 gtar -cvzf $(DIRKW).BE.tar.gz $(DIRKW)/B.E.)
#==============================================================================
#                   End of makefile
#==============================================================================
