#-------------------------------------------------------------------------------
#	make Test_Rmat
#	make Test_SfacISR
#	make Test_SfacFSR
#	make Test_UV
#	make Test_BornXsec
#	make Test_IsrSingle
#	make Test_FsrSingle
#       make Test_NewSingle
#	make Test_MultiPhot
#	make Test_Virtual
#	make Test_Rules
#	make Test_DsigOverDtau
#-------------------------------------------------------------------------------
#	make Tau-start
#	make Tau-debug
#-------------------------------------------------------------------------------
########### collection of HP  flags ###########
### HP options: +se for 100 or more columns, +B for backslash in strings
###             +K  for global SAVE
#FFLAGS  =  -K +se +B
########### collection of IBM  flags ###########
# IBM options: fixed=100 for 100 or more columns,
#              qextname for C-style external names with underscore at the end
#FFLAGS  =  -qfixed=132 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=132 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable
###########
FFLAGS  =  -qfixed=132 -qextname  -C -qextchk -g
LDFLAGS =  -qextname -bloadmap:lmap
##
F77 = f77
LD  = $(F77)
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
EXT_LIB4  = ../bornv/bornv.a
EXT_SRC4  = ../bornv/*.f	../bornv/*.h
EXT_LIB5  = ../tauola/glib.a
EXT_SRC5  = ../tauola/*.f
EXT_LIB6  = ../photos/glib.a
EXT_SRC6  = ../photos/photos.f	
###########
.f.o:
	$(F77) $(FFLAGS) -c $<
$(EXT_LIB1): $(EXT_SRC1)
	cd ../KK2f; make
$(EXT_LIB2): $(EXT_SRC2)
	cd ../glibk; make
$(EXT_LIB3): $(EXT_SRC3)
	cd ../jetset; make
$(EXT_LIB4): $(EXT_SRC4)
	cd ../bornv; make
$(EXT_LIB5): $(EXT_SRC5)
	cd ../tauola; make
$(EXT_LIB6): $(EXT_SRC6)
	cd ../photos; make
###########################
ERW.f:		ERW.h
STest.f:	STest.h
###########################
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Link procedure
COMMAND1 = Main.exe
OBJECTS1 = Main.o ERW.o STest.o Weyl.o brem.o
$(COMMAND1): $(OBJECTS1)  \
		$(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS1) \
		$(EXT_LIB1) $(EXT_LIB2) $(EXT_LIB3) $(EXT_LIB4) $(EXT_LIB5) $(EXT_LIB6)
#-----------------------------------------------------------------------
Tau-start:	$(COMMAND1)
	(cd ./Tau;  make start )
Tau-stop:
	(cd ./Tau;  make stop )
Tau-cont:
	(cd ./Tau;  make cont )
Tau-plot:
	(cd ./profig; make figmix)
Tau-start-debug:	$(COMMAND1)
	(cd ./Tau;  make start-debug )
Tau-cont-debug:	$(COMMAND1)
	(cp Tau/pro.hst.382 Tau/pro.hst )
	(cp Tau/semaphore.382 Tau/semaphore )
###	(cd ./Tau;  make cont )
	(cd ./Tau;  make cont-debug )
#-----------------------------------------------------------------------
Test_SfacISR: 	$(COMMAND1)
	(cd ./Tau;  cp SfacISR.input Tau.input ; make test )
	(diff   Tau/SfacISR.output Tau/pro.output | more)
Test_SfacFSR: 	$(COMMAND1)
	(cd ./Tau;  cp SfacFSR.input Tau.input ; make test )
	(diff   Tau/SfacFSR.output Tau/pro.output | more)
#-----------------------------------------------------------------------
Test_Rmat: 	$(COMMAND1)
	(cd ./Tau;  cp Rmat.input Tau.input ; make test )
	(diff   Tau/Rmat.output Tau/pro.output | more )
#-----------------------------------------------------------------------
Test_UV: 	$(COMMAND1)
	(cd ./Tau;  cp UV.input Tau.input ; make test )
	(diff   Tau/UV.output Tau/pro.output | more)
#-----------------------------------------------------------------------
Test_BornXsec: 	$(COMMAND1)
	(cd ./Tau;  cp BornXsec.input Tau.input ; make test )
	(diff   Tau/BornXsec.output Tau/pro.output | more)
Test_BornXsec-cont-debug:	$(COMMAND1)
	(cd ./Tau;  cp BornXsec.input Tau.input )
	(cp Tau/pro.hst.1402   Tau/pro.hst )
	(cp Tau/semaphore.1402 Tau/semaphore )
	(cd ./Tau;  make cont )
###	(cd ./Tau;  make cont-debug )
#-----------------------------------------------------------------------
Test_IsrSingle: 	$(COMMAND1)
	(cd ./Tau;  cp IsrSingle.input Tau.input ; make test )
	(diff   Tau/IsrSingle.output Tau/pro.output | more)
Test_IsrSingle-debug: 	$(COMMAND1)
	(cd ./Tau;  cp IsrSingle.input Tau.input )
	(cd ./Tau;  make cont-debug )
Test_FsrSingle: 	$(COMMAND1)
	(cd ./Tau;  cp FsrSingle.input Tau.input ; make test )
	(diff   Tau/FsrSingle.output Tau/pro.output | more)
Test_NewSingle: 	$(COMMAND1)
	(cd ./Tau;  cp NewSingle.input Tau.input ; make test )
	(diff   Tau/NewSingle.output Tau/pro.output | more)
#-----------------------------------------------------------------------
Test_MultiPhot: 	$(COMMAND1)
	(cd ./Tau;  cp MultiPhot.input Tau.input ; make test )
	(diff   ./Tau/MultiPhot.output ./Tau/pro.output | more )
Test_MultiPhot-cont-debug:	$(COMMAND1)
	(cd ./Tau;  cp MultiPhot.input Tau.input )
	(cp Tau/pro.hst.1402 Tau/pro.hst )
	(cp Tau/semaphore.1402 Tau/semaphore )
###	(cd ./Tau;  make cont )
	(cd ./Tau;  make cont-debug )
#-----------------------------------------------------------------------
Test_Virtual: 	$(COMMAND1)
	(cd ./Tau;  cp Virtual.input Tau.input ; make test )
	(diff   Tau/Virtual.output Tau/pro.output | more)
#-----------------------------------------------------------------------
Test_Rules: 	$(COMMAND1)
	(cd ./Tau;  cp Rules.input Tau.input ; make test )
	(diff   Tau/Rules.output Tau/pro.output | more)
#-----------------------------------------------------------------------
Test_DsigOverDtau: 	$(COMMAND1)
	(cd ./Tau;  cp DsigOverDtau.input Tau.input ; make test )
	(diff   Tau/DsigOverDtau.output Tau/pro.output | more)
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
#-----------------------------------
	(sed -e 's#<main>#Main#' -e 's#<dset>#Tau#' -e 's#<class>#long#'\
		./makefile.production > ./Tau/makefile)
#-----------------------------------
clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
#=============================================================================
# Update fortran compilation flags
makflag:
	(cd ../KK2f;   sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../jetset; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../tauola; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../photos; sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../glibk;  sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../dizet;  sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../bornv;  sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
	(cd ../bornv/tabtest;  \
	    sed -e 's#<FFLAGS>#$(FFLAGS)#' makefile.templ > makefile)
#-----------------------------------------------------------------------
Clean:
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd ../glibk;	make clean)
	(cd ../KK2f;	make clean)
	(cd ../jetset;	make clean)
	(cd ../tauola;	make clean)
	(cd ../photos;	make clean)
	(cd ../bornv;	make clean)
	(cd ../dizet;	make clean)
	(cd ../bornv/tabtest; make clean)
#-----------------------------------------------------------------------

