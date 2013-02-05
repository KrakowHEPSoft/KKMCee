#===============================================================================
#       make makflag		<==  in order to update f77 flags everywhere
#-------------------------------------------------------------------------------
#	make 189GeV-start
#	make 200GeV-start
#	make 206GeV-start
#	make input-quark
#===============================================================================
#  TWO WAYS OF ORGANIZING ELECTROWEAK LIBRARY INITIALIZATION.
#  If you want to READ EW formfactor from the disk (Dizet NOT linked) then
#  activate below   EXT_LIBRARY = EXT_LIBRARYa
############
all:	demo-start
########### ALPHA  flags ###########
#FFLAGS  =  -O -w -c -static
#FFLAGS  =  -extend_source -qfixed=132 -qextname -O  -C -qextchk
########### linux flags ###########
##FFLAGS  =  -fno-automatic -ffixed-line-length-132 -fbackslash -C -O
FFLAGS  =  -O -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash
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
F77 = g77
LD  = $(F77)
LDFLAGS =
#-----------------------------------------------------------------------
SRCLIST= "../KK2f ../glibk ../jetset ../bornv "
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
ProdMain.o:	RobAll.h
#-----------------------------------------------------------------------
#################################################
.f.o:
	$(F77) $(FFLAGS) -c $<
#########################################
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
#######################################
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
#### Make your own choice of the two choices below, see comments on two EW modes above
EXT_LIBRARYa = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4a) $(EXT_LIB5) $(EXT_LIB6)
EXT_LIBRARYd = $(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4d) $(EXT_LIB5) $(EXT_LIB6) $(EXT_LIB7)
###EXT_LIBRARY = $(EXT_LIBRARYa)   # Dizet is NOT linked, look-up tables from disk are used
EXT_LIBRARY = $(EXT_LIBRARYd)   # Dizet is linked, tables created in-flight during initialization
#-----------------------------------------------------------------------
# Link procedure
COMMAND1 = ProdMain.exe
OBJECTS1 = ProdMain.o
$(COMMAND1):                   $(OBJECTS1) $(EXT_LIBRARY)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) $(EXT_LIBRARY)
#-----------------------------------------------------------------------
189GeV-start:	$(COMMAND1)
	(cd ./189GeV;  make start )
189GeV-batch:	$(COMMAND1)
	(cd ./189GeV;  make batch )
189GeV-stop:
	(cd ./189GeV;  make stop )
189GeV-cont:
	(cd ./189GeV;  make cont )
189GeV-plot:
	(cd ./Plots; make Plotter-ps)
189GeV-start-debug:	$(COMMAND1)
	(cd ./189GeV;  make start-debug )
189GeV-cont-debug:	$(COMMAND1)
	(cp 189GeV/pro.hst.382 189GeV/pro.hst )
	(cp 189GeV/semaphore.382 189GeV/semaphore )
	(cd ./189GeV;  make cont )
###	(cd ./189GeV;  make cont-debug )
#-----------------------------------------------------------------------
200GeV-start:	$(COMMAND1)
	(cd ./200GeV;  make start )
200GeV-batch:	$(COMMAND1)
	(cd ./200GeV;  make batch )
200GeV-stop:
	(cd ./200GeV;  make stop )
200GeV-cont:
	(cd ./200GeV;  make cont )
200GeV-plot:
	(cd ./Plots; make Plotter-ps)
#-----------------------------------------------------------------------
206GeV-start:	$(COMMAND1)
	(cd ./206GeV;  make start )
206GeV-batch:	$(COMMAND1)
	(cd ./206GeV;  make batch )
206GeV-stop:
	(cd ./206GeV;  make stop )
206GeV-start-debug:	$(COMMAND1)
	(cd ./206GeV;  make start-debug )
#-----------------------------------------------------------------------
# Link procedure
COMMAND2 = RealMain.exe
OBJECTS2 = RealMain.o
$(COMMAND2):                   $(OBJECTS2) $(EXT_LIBRARY)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXT_LIBRARY)
#-----------------------------------------------------------------------
189Rea-start:	$(COMMAND2)
	(cd ./189Rea;  make start )
189Rea-stop:
	(cd ./189Rea;  make stop )
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
	(sed -e 's#<main>#ProdMain#' -e 's#<dset>#189GeV#'   -e 's#<class>#8nh#' ./makefile.production > ./189GeV/makefile)
	(sed -e 's#<main>#ProdMain#' -e 's#<dset>#200GeV#'   -e 's#<class>#8nh#' ./makefile.production > ./200GeV/makefile)
	(sed -e 's#<main>#ProdMain#' -e 's#<dset>#206GeV#'   -e 's#<class>#8nh#' ./makefile.production > ./206GeV/makefile)
	(sed -e 's#<main>#RealMain#' -e 's#<dset>#189Rea#'   -e 's#<class>#8nh#' ./makefile.production > ./189Rea/makefile)
input-mu:
	(cp ./189GeV/Mu.input  ./189GeV/189GeV.input )
	(cp ./200GeV/Mu.input  ./200GeV/200GeV.input )
	(cp ./206GeV/Mu.input  ./206GeV/206GeV.input )
input-tau:
	(cp ./189GeV/Tau.input ./189GeV/189GeV.input )
	(cp ./200GeV/Tau.input ./200GeV/200GeV.input )
	(cp ./206GeV/Tau.input ./206GeV/206GeV.input )
input-quark:
	(cp ./189GeV/Quark.input ./189GeV/189GeV.input )
	(cp ./200GeV/Quark.input ./200GeV/200GeV.input )
	(cp ./206GeV/Quark.input ./206GeV/206GeV.input )
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
	(cd ../KKsem;   sed -e 's#<FFLAGS>#$(FFLAGS)#' -e 's#<F77>#$(F77)#' makefile.templ > makefile)
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
##############################################################################
default:
	(cd ../dizet/; make default)
YR:
	(cd ../dizet/; make YR)
#-----------------------------------------------------------------------
gtar: Clean
	(cd ../../; gtar -cvzf KK2f-all.tar.gz  KK2f-all )
	(cd ../../; chmod -w KK2f-all.tar.gz)
#-----------------------------------------------------------------------
