#=============================================================================
#	The example of ploting in batch program
#	this one will plot R(s), biut more complicated can be costructed
#		make SkyLine
#		make Plot2b		Many distribs...
#		make Plot2c		Many distribs...
#-----------------------------------------------------------------------------
#	The other ploting ROOT scripts are activated as follows:
#		root Fig0.C		xcheck of normalization procedure
#		root Fig1.C		xsection(nphot), dsig/d(theta_gamma)
#		root Fig1b.C		Pions, pT and CoTheta 
#		root Fig2.C		Q2(pions) with various cuts
#		root Fig2b.C		Q2 distr. 2pi versus muons
#		root Fig2c.C		Q2 distr. NEW
#		root Fig3.C		Higer order effects		
#-----------------------------------------------------------------------------
#	The input root files with histograms are in ../demoC/
#	Note that variable RPATH for Root has to be readjusted
#=============================================================================
.SUFFIXES:  .o .c .cxx .f .l  .s  .sh  .h  .a  .C
#-----------------------------------------------------------------------------
#  Root 
RPATH = /home/jadach/root
#  Root  libraries for linker
RLIB=   -L$(RPATH)/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
	-lPostscript -lMatrix -lPhysics -lm -ldl -lpthread -rdynamic
FLIB=   -lf2c
##======================================================================
#####  f77 compiler
# forbid underscore completely       -fno-second-underscore
# forbid double underscore only      -fno-underscoring
FFLAGS  =  -O -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash
F77 = g77
#  Fortran libraries for linker
#======================================================================
#####  c++ compiler
#CXXFLAGS    = -g -ansi -pedantic
#CXXFLAGS    = -O -Wall -fPIC -D_REENTRANT -I/home/jadach/root/include
CXXFLAGS    = -O  -ansi -Wall -fPIC -D_REENTRANT -I$(RPATH)/include
#CXXFLAGS    =  -g
CXX         =  g++
#CXX         =  insure
#CXX         =  codewizard
CLD         =  $(CXX)
#-----------------------------
# Compile all c++ files
#-----------------------------
.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<
.f.o:
	$(F77) $(FFLAGS) -c $<
.c.o:
	gcc -c $<
#================================================================================
SRCLIST= "../KK2f ../glibk ../jetset ../bornv "
#----------------------
EXT_LIB0  = ../../KKsem/KKsem.a
EXT_SRC0  = ../../KKsem/*.f  ../../KKsem/*.h
$(EXT_LIB0): $(EXT_SRC0)
	(cd ../../KKsem; make)
#############################################################
EXT_LIB1  = ../../KK2f/KK2f.a
EXT_SRC1  = ../../KK2f/*.f		../../KK2f/*.h
EXT_LIB2  = ../../glibk/glibk.a
EXT_SRC2  = ../../glibk/*.f	../../glibk/*.h
EXT_LIB3  = ../../jetset/jetset.a
EXT_SRC3  = ../../jetset/*.f
EXT_LIB4t  = ../../bornv/bornv-tabs.a
EXT_LIB4d  = ../../bornv/bornv-dizet.a
EXT_SRC4  = ../../bornv/*.f	../../bornv/*.h
EXT_LIB5  = ../../tauola/glib.a
EXT_SRC5  = ../../tauola/*.f
EXT_LIB6  = ../../photos/glib.a
EXT_SRC6  = ../../photos/photos.f	
EXT_LIB7  = ../../dizet/dizet.a
EXT_SRC7  = ../../dizet/*.f
###########
$(EXT_LIB1): $(EXT_SRC1)
	cd ../../KK2f; make
$(EXT_LIB2): $(EXT_SRC2)
	cd ../../glibk; make
$(EXT_LIB3): $(EXT_SRC3)
	cd ../../jetset; make
$(EXT_LIB4t): $(EXT_SRC4)
	cd ../../bornv; make
$(EXT_LIB4d): $(EXT_SRC4)
	cd ../../bornv; make
$(EXT_LIB5): $(EXT_SRC5)
	cd ../../tauola; make
$(EXT_LIB6): $(EXT_SRC6)
	cd ../../photos; make
$(EXT_LIB7): $(EXT_SRC7)
	cd ../../dizet; make
#===================================================================================================
#### Make your own choice of the two choices below, see comments on two EW modes above
EXT_LIBRARYt = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4t) $(EXT_LIB5) $(EXT_LIB6)
EXT_LIBRARYd = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4d) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7)
# For plots we actuaaly using these ones:
EXT_LIBRARY =  $(EXT_LIB4t) $(EXT_LIB1) $(EXT_LIB3)
EXT_LIBRARY3 = $(EXT_LIB0) $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4t) $(EXT_LIB5) $(EXT_LIB6)
#===================================================================================================
MAIN1 = SkyLine
OBJECTS1 = $(MAIN1).o
$(MAIN1).exe:                $(OBJECTS1) $(EXT_LIBRARY)
	g++ -o $(MAIN1).exe  $(OBJECTS1) $(EXT_LIBRARY)  $(RLIB) $(FLIB)
$(MAIN1):	$(MAIN1).exe
	./$(MAIN1).exe
#===================================================================================================
MAIN2 = Plot2b
OBJECTS2 = $(MAIN2).o
$(MAIN2).exe:                $(OBJECTS2) $(EXT_LIBRARY3)
	g++ -o $(MAIN2).exe  $(OBJECTS2) $(EXT_LIBRARY3)  $(RLIB) $(FLIB)
$(MAIN2):	$(MAIN2).exe
	./$(MAIN2).exe
#===================================================================================================
MAIN3 = Plot2c
OBJECTS3 = $(MAIN3).o
$(MAIN3).exe:                $(OBJECTS3) $(EXT_LIBRARY3)
	g++ -o $(MAIN3).exe  $(OBJECTS3) $(EXT_LIBRARY3)  $(RLIB) $(FLIB)
$(MAIN3):	$(MAIN3).exe
	./$(MAIN3).exe
#=====================================================================
clean:
	(rm -f *.o *.a *~ *.exe core *.gz lmap *.out *.dvi *.aux *.log *.so)
#==============================================================================
