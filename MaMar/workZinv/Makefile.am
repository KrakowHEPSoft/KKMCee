INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO
#-----------------------------------------------------------------------
#########  FARMING on NZ42 farm ###########
#  make prepare
#  (make CLEAN;   make farm24-KRK)
#  make SLsubmit
#  make q-semaphore   make combine
#-----------------------------------------------------------------------
MCMD=ProdZinv
DSET=workZinv
CLASS=qunlimitted
##-----------------------------------------------------------------------
SRCLIST= "../ ../../yfslib ../../bornv ../../glibk ../../KK2f"
#-----------------------------------------------------------------------
# Interactive:  run start stop and continue
#-----------------------------------------------------------------------
../${MCMD}:
	(cd ../; make)
#-----------------------------------------------------------------------
prepare:
	(ln -s ../../.KK2f_defaults ./.KK2f_defaults )
	(ln -s ../../dizet          ./dizet )
#-----------------------------------------------------------------------
start: ../${MCMD}
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 time ../${MCMD}			)&
#-----------------------------------------------------------------------
debug:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 ddd ../${MCMD}			)&
#-----------------------------------------------------------------------
backg:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 time ../${MCMD} > /dev/null	)&
#-----------------------------------------------------------------------
stop:
	(cp ../semaphore.stop ./semaphore )&
#-----------------------------------------------------------------------
cont:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         time ../${MCMD}			)&
#-----------------------------------------------------------------------
start-debug:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 xldb -I  ${SRCLIST} ../${MCMD}	)&
#-----------------------------------------------------------------------
cont-debug:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         xldb -I  ${SRCLIST} ../${MCMD}	)&
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#   Batch on farm:  initialize, submit and stop
#-----------------------------------------------------------------------
farm6-LSF:
	(csh ../farming/farm-init "1 2 3 4 5 6" \
	${MCMD} ${DSET} LSF ${CLASS} )
farm12-LSF:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MCMD} ${DSET} LSF ${CLASS} )
farm24-LSF:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	      13 14 15 16 17 18 19 20 21 22 23 24 " \
	${MCMD} ${DSET} LSF ${CLASS} )
farm12-SP:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MCMD} ${DSET} SP ${CLASS} )
farm24-SP:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	      13 14 15 16 17 18 19 20 21 22 23 24 " \
	${MCMD} ${DSET} SP ${CLASS} )
farm6-CSF:
	(csh ../farming/farm-init "1 2 3 4 5 6" \
	${MCMD} ${DSET} CSF ${CLASS} )
farm12-CSF:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MCMD} ${DSET} CSF ${CLASS} )
###-----------
farm6-KRK:
	(csh ../farming/farm-init "1 2 3 4 5 6" \
	${MCMD} ${DSET} KRK ${CLASS} )
farm12-KRK:
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MCMD} ${DSET} KRK ${CLASS} )
farm24-KRK:
	(csh ../farming/farm-init  "1 2 3 4 5 6 7 8 9 10 11 12 \
	      13 14 15 16 17 18 19 20 21 22 23 24 " \
	${MCMD} ${DSET} KRK ${CLASS} )
###-----------
farm6-CC1: CLEAN
	(csh ../farming/farm-init "1 2 3 4 5 6" \
	${MCMD} ${DSET} CC1 ${CLASS} )
farm12-CC1: CLEAN
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MCMD} ${DSET} CC1 ${CLASS} )
farm24-CC1: CLEAN
	(csh ../farming/farm-init "1 2 3 4 5 6 7 8 9 10 11 12  \
	                  13 14 15 16 17 18 19 20 21 22 23 24" \
	${MCMD} ${DSET} CC1 ${CLASS} )
#-----------------------------------------------------------------------
farm-gtar: prepare
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
bsubmit6:
	(csh ../farming/farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6")
bsubmit12:
	(csh ../farming/farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12")
bsubmit18:
	(csh ../farming/farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18")
bsubmit24:
	(csh ../farming/farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#######---------
qsubmit6:
	(csh ../farming/farm-qsubmit ${DSET} "1 2 3 4 5 6 ")
qsubmit12:
	(csh ../farming/farm-qsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
qsubmit18:
	(csh ../farming/farm-qsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 ")
qsubmit24:
	(csh ../farming/farm-qsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
llsubmit6:
	(csh ../farming/farm-llsubmit ${DSET} "1 2 3 4 5 6")
llsubmit12:
	(csh ../farming/farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
llsubmit18:
	(csh ../farming/farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18")
llsubmit24:
	(csh ../farming/farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#-----------------------------------------------------------------------
stop6:
	(csh ../farming/farm-stop   ${DSET} "1 2 3 4 5 6 ")
stop12:
	(csh ../farming/farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
stop18:
	(csh ../farming/farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 ")
stop24:
	(csh ../farming/farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#
cphst24:
	(csh ../farming/farm-cphst   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
########################################################################
SLsubmit:
	($(ROOTEXEC)  -b -q -l '../../MCdev/farming/SubmFarmS.C("sbatch",200)' )
combine:
	(hadd -f ./rmain.root ./farm/*/rmain.root)
########################################################################
#-----------------------------------------------------------------------
#--------- cobmine (add) results from several subdirectories -----------
Combine:
	(cd ../../glibk; make  -f KKMakefile combine.exe)
combine6:  combine
	(sed -e 's#?.hst#pro.hst#' -e 's#?#6#' ../combine.list \
	| ../../glibk/combine.exe; more pro.hst)
combine12:  combine
	(sed -e 's#?.hst#pro.hst#' -e 's#?#12#' ../combine.list \
	| ../../glibk/combine.exe; more pro.hst)
combine18:  combine
	(sed -e 's#?.hst#pro.hst#' -e 's#?#18#' ../combine.list \
	| ../../glibk/combine.exe; more pro.hst)
combine24:  combine
	(sed -e 's#?.hst#pro.hst#' -e 's#?#24#' ../combine.list \
	| ../../glibk/combine.exe; more pro.hst)
#-----------------------------------------------------------------------
#---- query
q-sem:
	(ls -altr farm/*/semaphore)
q-nev:
	(cat farm/*/semaphore | grep nev)
q-cont:
	(cat farm/*/semaphore | grep cont)
#-----------------------------------------------------------------------
clean:                                                     
	(rm -f *.o *.a *~ *.exe core *.err *.out)
CLEAN:  clean
	(rm  -R -f farm )
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Miscelaneous
#-----------------------------------------------------------------------
script-SP:
	(csh ../farming/farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MCMD} ${DSET} SP ${CLASS} )
script-SP850:
	(csh ../farming/farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MCMD} ${DSET} SP850 ${CLASS} )
#-----------------------------------------------------------------------
