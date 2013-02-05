#===============================================================================
#       make makflag		<==  in order to update f77 flags everywhere
#-------------------------------------------------------------------------------
#	make demoC-start	<==  basic demonstration program
#	make demoC-view		<==  preview histograms 
#===============================================================================
#  TWO WAYS OF ORGANIZING ELECTROWEAK LIBRARY INITIALIZATION.
#  If you want to READ EW formfactor from the disk (Dizet NOT linked) then
#  activate below   EXT_LIBRARY = EXT_LIBRARYa
#----------------------------------------------------------------------
line_count:
	(cat *.cxx *.h | wc -l )
#----------------------------------------------------------------------
.SUFFIXES:  .o .c .cxx .f .l  .s  .sh  .h  .a
#======================================================================
#RPATH = /afs/cern.ch/na49/library.4/ROOT/new/i386_redhat60/root
#RPATH = /afs/cern.ch/user/j/jadachsp/private/root
RPATH = /home/jadach/root
#  Root  libraries for linker
RLIB=   -L$(RPATH)/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
	-lPostscript -lMatrix -lPhysics -lm -ldl -lpthread -rdynamic
#  Fortran libraries for linker
FLIB=   -lg2c
#================================================================================
########### ALPHA  flags ###########
#FFLAGS  =  -O -w -c -static
#FFLAGS  =  -extend_source -qfixed=132 -qextname -O  -C -qextchk
########### collection of HP  flags ###########
### +es for 100 or more columns, +B for backslash in strings, +K  for global SAVE
#FFLAGS  =  -K +es +B -O
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
#FFLAGS  =  -qfixed=132 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=132 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable -bloadmap:lmap
#FFLAGS  =  -qfixed=132 -qextname  -C -qextchk -O
##======================================================================
########### linux flags ###########
#FFLAGS  =  -fno-automatic -ffixed-line-length-132 -fbackslash -C -O
#####  f77 compiler
# forbid underscore completely       -fno-second-underscore
# forbid double underscore only      -fno-underscoring
FFLAGS  =  -O -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash
F77 = g77
#======================================================================
#####  g++ compiler
# -shared -Wl,-soname
#CXXFLAGS    = -g -ansi -pedantic
#CXXFLAGS    = -O -Wall -fPIC -D_REENTRANT -I$(RPATH)/include
CXXFLAGS    = -O  -ansi -Wall -fPIC -D_REENTRANT -I$(RPATH)/include
CXX         =  g++
#CXX         =  insure
#CXX         =  codewizard
CLD         =  $(CXX)
#-------------------------------------------------------
# Compile all c++, f77 and c  files
.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<
.f.o:
	$(F77) $(FFLAGS) -c $<
.c.o:
	gcc -c $<
#================================================================================
SRCLIST= "../KK2f ../glibk ../jetset ../bornv "
#-----------------------------------------------------------------------
EXT_LIB1  = ../KK2f/KK2f.a
EXT_SRC1  = ../KK2f/*.f		../KK2f/*.h
EXT_LIB2  = ../glibk/glibk.a
EXT_SRC2  = ../glibk/*.f	../glibk/*.h
EXT_LIB3  = ../jetset/jetset.a
EXT_SRC3  = ../jetset/*.f
EXT_LIB4a = ../bornv/bornv-tabs.a
EXT_LIB4d = ../bornv/bornv-dizet.a
EXT_SRC4  = ../bornv/*.f	../bornv/*.h
EXT_LIB5  = ../tauola/glib.a
EXT_SRC5  = ../tauola/*.f
EXT_LIB6  = ../photos/glib.a
EXT_SRC6  = ../photos/photos.f	
EXT_LIB7  = ../dizet/dizet.a
EXT_SRC7  = ../dizet/*.f	../dizet/*.h	../bornv/*.h
###########
$(EXT_LIB1): $(EXT_SRC1)
	cd ../KK2f; make
$(EXT_LIB2): $(EXT_SRC2)
	cd ../glibk; make
$(EXT_LIB3): $(EXT_SRC3)
	cd ../jetset; make
$(EXT_LIB4a): $(EXT_SRC4)
	cd ../bornv; make
$(EXT_LIB4d): $(EXT_SRC4)
	cd ../bornv; make
$(EXT_LIB5): $(EXT_SRC5)
	cd ../tauola; make
$(EXT_LIB6): $(EXT_SRC6)
	cd ../photos; make
$(EXT_LIB7): $(EXT_SRC7)
	cd ../dizet; make
#===================================================================================================
MainP.o:  ROBOL.o ROBOL2.o KKMC.o 
KKMC.o:		KKMC.h PartLund.h
ROBOL.o:	ROBOL.h
ROBOL2.o:	ROBOL2.h
Semaph.o:	Semaph.h
PartLund.o: 	PartLund.h
#===================================================================================================
#### Make your own choice of the two choices below, see comments on two EW modes above
EXT_LIBRARYa = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4a) $(EXT_LIB5) $(EXT_LIB6)
EXT_LIBRARYd = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4d) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7)
# Note that benchmark ouputs are for EXT_LIBRARYd
EXT_LIBRARY = $(EXT_LIBRARYa)   # Dizet is NOT linked, look-up tables from disk are used
####EXT_LIBRARY = $(EXT_LIBRARYd)   # Dizet is linked, tables created in-flight during initialization
#-----------------------------------------------------------------------
# Link procedure
COMMAND1 = ProdMC.exe
OBJECTS1 = ProdMC.o
$(COMMAND1):                   $(OBJECTS1) $(EXT_LIBRARY)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) $(EXT_LIBRARY)
#-----------------------------------------------------------------------
COMMAND2 = MainP.exe
OBJECTS2 = MainP.o Semaph.o ROBOL.o ROBOL2.o KKMC.o PartLund.o
$(COMMAND2):                    $(OBJECTS2) $(EXT_LIBRARY)
	$(CLD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXT_LIBRARY) $(RLIB) $(FLIB)
#	$(CLD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXT_LIBRARY) $(RLIB) $(FLIB)
#-----------------------------------------------------------------------
demoC-start:	$(COMMAND2)
	(cd ./demoC;  make start )
demoC-cont:	$(COMMAND2)
	(cd ./demoC;  make cont )
demoC-stop:	$(COMMAND2)
	(cd ./demoC;  make stop )
demoC-debug:	$(COMMAND2)
	(cd ./demoC;  make debug )
demoC-view:
	(cd ./demoC;  root canvas.C)
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
	(sed -e 's#<main>#MainP#'   -e 's#<dset>#demoC#' -e 's#<class>#8nh#' ./makefile.production > ./demoC/makefile)
#-----------------------------------
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd ./Plots;	make clean)
	(cd ./dok;	make clean)
	(cd ./demoC;	make CLEAN)
#=============================================================================
# Update fortran compilation flags
makflag:
	(cd ../KK2f;   sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../jetset; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../tauola; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../photos; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../glibk;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../bornv;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../bornv/tabtest; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../dizet;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../KKsem;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#'  makefile.templ > makefile)
mbr:
	(cd ../bornv;  make mbr)
default:
	(cd ../dizet;  make default)
#-----------------------------------------------------------------------
Clean:	clean
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd ./Plots;	make clean)
	(cd ../glibk;	make clean)
	(cd ../KK2f;	make clean)
	(cd ../jetset;	make clean)
	(cd ../tauola;	make clean)
	(cd ../photos;	make clean)
	(cd ../bornv;	make clean)
	(cd ../dizet;	make clean)
	(cd ../bornv/tabtest; make clean)
#-----------------------------------------------------------------------
gtar: Clean
	(cd ../../; gtar -cvzf KK2f-all.tar.gz  KK2f-all )
	(cd ../../; chmod -w KK2f-all.tar.gz)
#-----------------------------------------------------------------------
