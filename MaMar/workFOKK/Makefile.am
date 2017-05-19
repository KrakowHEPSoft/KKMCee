#=======================================================================
# Invokes all kind of scripts which setup (Start.C) and execute 
# program MAIN interactively and/or on the PC-farm (batch)
# Part edited by hand (Below part edited by Automake Manager od KDevelop)
#=======================================================================
MAIN=MainFOAM
DSET=workFOKK
## CERN  subCMC=bsub, CLASS= 8nm, 1nh, 1nd, 1nw
## Local qubCMC=qsub, CLASS= qunlimited
subCMD=qsub
CLASS=qunlimitted
##############
check_all:
	(cd ..; make )
start:  check_all
	rm -f ./pro.output ./histo.root ./mcgen.root
	$(ROOTEXEC)  -b -q -l ./Start.C
	time ../${MAIN}
stop:
	$(ROOTEXEC)  -b -q -l ../../MCdev/farming/Stop.C
cont:
	time ../${MAIN}
debug:  check_all
	ddd --directory="../" ../../bin/MainPr
#########################################################################
####      NEW UNIVERSAL VERSION tested on the local farm
#########################################################################
qfarm6:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start0.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qfarm24:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start0.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",24)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qfarm55:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start0.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",55)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qfarm40:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start0.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",40)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
submit6:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",6)' )
qsubmitall:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",200)' )
########################################################################
farm-gtar:
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
farm-stop:
	($(ROOTEXEC)  -b -q -l ../../MCdev/farming/StopFarm.C)
#-----------------------------------------------------------------------
run2:
	(cd ./farm/1/; csh ${DSET}.cmd.1)&
	(cd ./farm/2/; csh ${DSET}.cmd.2)&
#--------- cobmine (add) results from several subdirectories -----------
combine:
	($(ROOTSYS)/bin/hadd -f ./histo.root ./farm/*/histo.root)
	(cp -p ./farm/1/mcgen.root ./)
#-----------------------------------------------------------------------
#---- query
q-sem:
	(ls -altr ./farm/*/mcgen.root)
q-nev:
	($(ROOTEXEC)  -b -q -l ../../MCdev/farming/q-nev.C)
#-----------------------------------------------------------------------
farm-clean:
	(rm  -R -f farm *.bak )
#=======================================================================
#  Part edited by Automake Manager of KDvelelop
#=======================================================================
#INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/MCdev $(all_includes) \
#	-I$(ROOTINCDIR)
#METASOURCES = AUTO
#noinst_HEADERS = StartH.C
#========================================================================
#========================================================================
