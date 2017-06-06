#=================================================================================
# Invokes all kind of scripts which setup (Start.C) and executes MAIN
#    make start; make stop
# Seting up input files/scripts and runing on the NZ42 PC-farm
#    make sfarm6; make ssubmit6; make q-sem; make q-nev
#    make farm-stop; make combine
# Farming under SLURM batch system on CC1
#    make SLfarm6; make SLsubmitall; ...
# https://www.rc.fas.harvard.edu/resources/documentation/convenient-slurm-commands/
#==================================================================================
MAIN=MainFOAM
DSET=workFOAM
## CERN  subCMC=bsub, CLASS= 8nm, 1nh, 1nd, 1nw
## Local qubCMC=qsub, CLASS= qunlimited
subCMD=qsub
CLASS=qunlimitted
##############
.KK2f_defaults:
	(ln -s ../../.KK2f_defaults ./)
check_all: .KK2f_defaults
	(cd ..; make )
start_n: check_all
	rm -f ./pro.output ./histo.root ./mcgen.root
	cp ./$(DSET).input ./pro.input
	$(ROOTEXEC)  -b -q -l ./Start.C
start:  start_n
	time ../$(MAIN)&
stop:
	$(ROOTEXEC)  -b -q -l ../../MCdev/farming/Stop.C
cont:
	time ../${MAIN} &
debug:  check_all
	ddd --directory="../" ../../bin/MainPr
#########################################################################
####      NEW UNIVERSAL VERSION tested on the local farm
#########################################################################
qfarm6:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qfarm40:	farm-clean check_all
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmQ.C("$(DSET)",40)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qsubmit6:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",6)' )
qqsubmitall:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",200)' )
#########################################################################
####   version for SLURM on NZ42 cluster ####
sfarm6:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmS.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
sfarm24:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmS.C("$(DSET)",24)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
ssubmit6:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmS.C("sbatch",6)' )
ssubmitall:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmS.C("sbatch",200)' )
#########################################################################
####   version for CC1 cluster ####
SLfarm6:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm84:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../../MCdev/farming/SetFarmSL.C("$(DSET)",86)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLsubmitall:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmS.C("sbatch",200)' )
#########################################################################
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
