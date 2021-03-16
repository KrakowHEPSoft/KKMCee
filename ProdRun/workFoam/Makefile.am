#=======================================================================
# Invokes all kind of scripts which setup (Start.C) and executes MAIN
#    make start; make stop
# Seting up input files/scripts and runing on the NZ42 PC-farm
# Farming under SLURM batch system on CC1
#    make SLfarm6; make SLsubmitall; ...
#=======================================================================
MAIN=MainFoam
DSET=workFoam
subCMD=qsub
CLASS=qunlimitted
##############
KKMCee_defaults:
	(ln -s ../../SRCee/KKMCee_defaults ./)

EWtables: KKMCee_defaults
	(cd $(top_builddir); make)
	cp ./$(DSET).input ./pro.input
	($(top_builddir)/DZface/TabMainC;)
##############
check_all: EWtables
	(cd $(top_srcdir); make)
start_n: check_all
	rm -f ./pro.output ./histo.root ./mcgen.root
	cp ./$(DSET).input ./pro.input
	$(ROOTEXEC)  -b -q -l ./Start.C
start:  start_n
	time ../${MAIN} &
stop:
	$(ROOTEXEC)  -b -q -l ../../MCdev/farming/Stop.C
cont:
	time ../${MAIN} &
debug:  check_all
	ddd --directory="../" ../../bin/MainPr
#########################################################################
####      NEW UNIVERSAL VERSION tested on the local farm
#########################################################################
####   version for CC1 cluster ####
SLfarm6:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm24:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",24)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm48:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",48)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm72:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",72)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm96:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",96)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLsubmitall:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmS.C("sbatch",200)' )
########################################################################
farm-gtar:
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
farm-stop:
	($(ROOTEXEC)  -b -q -l ../../MCdev/farming/StopFarm.C)
#--------- cobmine (add) results from several subdirectories -----------
combine:
#	($(ROOTSYS)/bin/hadd -f ./histo.root ./farm/*/histo.root)
	(hadd -f ./histo.root ./farm/*/histo.root)
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
#========================================================================
#========================================================================
