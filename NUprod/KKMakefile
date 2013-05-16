#===============================================================================
#       make makflag		<==  in order to update f77 flags everywhere
#-------------------------------------------------------------------------------
#	make NuAll-start
#	make ew-mu2
#	make ew-numu
#	make ew-nuel
#===============================================================================
#  TWO WAYS OF ORGANIZING ELECTROWEAK LIBRARY INITIALIZATION.
#  If you want to READ EW formfactor from the disk (Dizet NOT linked) then
#  activate below   EXT_LIBRARY = EXT_LIBRARYa
#-------------------------------------------------------------------------------
############
all:	demo-start
########### ALPHA  flags ###########
#FFLAGS  =  -O -w -c -static
#FFLAGS  =  -extend_source -qfixed=132 -qextname -O  -C -qextchk
########### linux flags ###########
FFLAGS  =  -fno-automatic -ffixed-line-length-132 -fbackslash -C -O
########### collection of HP  flags ###########
### +es for 100 or more columns, +B for backslash in strings, +K  for global SAVE
#FFLAGS  =  -K +es +B -O
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
#FFLAGS  =  -qfixed=132 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=132 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable -bloadmap:lmap
#FFLAGS  =  -qfixed=132 -qextname  -C -qextchk -O
###########
# under linux g77 is often used
F77 = g77
#F77 = f77 
LD  = $(F77)
LDFLAGS =
#-----------------------------------------------------------------------
SRCLIST= "../KK2f ../glibk ../jetset ../bornv "
#-----------------------------------------------------------------------
###########
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
EXT_SRC7  = ../dizet/*.f
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
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
#----------------------
EXT_LIB0  = ../KKsem/KKsem.a
EXT_SRC0  = ../KKsem/*.f  ../KKsem/*.h
$(EXT_LIB0): $(EXT_SRC0)
	(cd ../KKsem; make)
#===================================================================================================
#### Make your own choice of the two choices below, see comments on two EW modes above
EXT_LIBRARYa = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4a) $(EXT_LIB5) $(EXT_LIB6)
EXT_LIBRARYd = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4d) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7)
EXT_LIBRARY = $(EXT_LIBRARYa)   # Dizet is NOT linked, look-up tables from disk are used
###EXT_LIBRARY = $(EXT_LIBRARYd)   # Dizet is linked, tables created in-flight during initialization
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Link procedure
COMMAND2 = demo.exe
OBJECTS2 = demo.o
$(COMMAND2):                   $(OBJECTS2) $(EXT_LIBRARY)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXT_LIBRARY)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Link procedure
COMMAND1 = ProdNU.exe
OBJECTS1 = ProdNU.o
$(COMMAND1):                   $(OBJECTS1) $(EXT_LIBRARY)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) $(EXT_LIBRARY)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
demo-start:	$(COMMAND2)
	(cd ./demo;  make start )
#-----------------------------------------------------------------------
NuAll-start:	$(COMMAND1)
	(cd ./NuAll;  make start )
NuAll-stop:
	(cd ./NuAll;  make stop )
NuAll-start-debug:	$(COMMAND1)
	(cd ./NuAll;  make start-debug )
NuAll-cont-debug:	$(COMMAND1)
	(cp NuAll/pro.hst.110 NuAll/pro.hst )
	(cp NuAll/semaphore.110 NuAll/semaphore )
###	(cd ./NuAll;  make cont )
	(cd ./NuAll;  make cont-debug )
##############################################################################
# 	gmake bornx
COMMAND9 = bornx.exe
OBJECTS9 = bornx.o
$(COMMAND9): $(OBJECTS9)                   $(EXT_LIB0) $(EXT_LIBRARYa)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS9) $(EXT_LIB0) $(EXT_LIBRARYa)
#-----------------------------------------------------------------------
bornx:	$(COMMAND9)
	(cp bornx.input Semi/)
	(cd Semi; time ../$(COMMAND9) )
#-----------------------------------------------------------------------
ew-mu2:	$(COMMAND9)
	(cp input-ew-mu2 Semi/bornx.input)
	(cd Semi; time ../$(COMMAND9) )
#-----------------------------------------------------------------------
ew-numu:	$(COMMAND9)
	(cp input-ew-numu Semi/bornx.input)
	(cd Semi; time ../$(COMMAND9) )
#-----------------------------------------------------------------------
ew-nuel:	$(COMMAND9)
	(cp input-ew-nuel Semi/bornx.input)
	(cd Semi; time ../$(COMMAND9) )
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
	(sed -e 's#<main>#ProdNU#' -e 's#<dset>#NuAll#' -e 's#<class>#long#' ./makefile.production > ./NuAll/makefile)
	(sed -e 's#<main>#ProdNU#' -e 's#<dset>#NuEl#'  -e 's#<class>#long#' ./makefile.production > ./NuEl/makefile)
	(sed -e 's#<main>#ProdNU#' -e 's#<dset>#NuMu#'  -e 's#<class>#long#' ./makefile.production > ./NuMu/makefile)
#-----------------------------------
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
#=============================================================================
# Update fortran compilation flags
makflag:
	(cd ./Plots;   sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../KK2f;   sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../jetset; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../tauola; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../photos; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../glibk;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../bornv;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../bornv/tabtest; sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../dizet;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
	(cd ../KKsem;  sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#'  makefile.templ > makefile)
#-----------------------------------------------------------------------
Clean:
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
