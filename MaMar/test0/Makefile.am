AM_CPPFLAGS = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO

#-----------------------------------------------------------------------
#       Template for makefile in production directory demo
#-----------------------------------------------------------------------
MAIN=Prod0
DSET=test0
CLASS=long
#-----------------------------------------------------------------------
SRCLIST= "../ ../../yfslib ../../bornv ../../glibk ../../KK2f"
#-----------------------------------------------------------------------
# Interactive:  run start stop and continue
#-----------------------------------------------------------------------
../${MAIN}:
	(cd ..; make ../${MAIN})
prepare:
	(ln -s ../../.KK2f_defaults ./.KK2f_defaults )
	(ln -s ../../dizet          ./dizet )
#-----------------------------------------------------------------------
start: 
	(cd $(top_builddir); make -s)
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 time ../${MAIN}			)&
#-----------------------------------------------------------------------
debug:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 ddd ../${MAIN}			)&
#-----------------------------------------------------------------------
backg:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 time ../${MAIN} > /dev/null	)&
#-----------------------------------------------------------------------
stop:
	(cp ../semaphore.stop ./semaphore )&
#-----------------------------------------------------------------------
cont:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         time ../${MAIN}			)&
#-----------------------------------------------------------------------
start-debug:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 xldb -I  ${SRCLIST} ../${MAIN}	)&
#-----------------------------------------------------------------------
cont-debug:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         xldb -I  ${SRCLIST} ../${MAIN}	)&
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#   Batch on farm:  initialize, submit and stop
#-----------------------------------------------------------------------
farm6-LSF:
	(csh ../farm-init "1 2 3 4 5 6" \
	${MAIN} ${DSET} LSF ${CLASS} )
farm12-LSF:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MAIN} ${DSET} LSF ${CLASS} )
farm24-LSF:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	      13 14 15 16 17 18 19 20 21 22 23 24 " \
	${MAIN} ${DSET} LSF ${CLASS} )
#-----------------------------------------------------------------------
farm-gtar:
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
bsubmit6:
	(csh ../farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6")
bsubmit12:
	(csh ../farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12")
bsubmit18:
	(csh ../farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18")
bsubmit24:
	(csh ../farm-bsubmit ${DSET} ${CLASS} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#-----------------------------------------------------------------------
stop6:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 ")
stop12:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
stop18:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 ")
stop24:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#
cphst24:
	(csh ../farm-cphst   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#-----------------------------------------------------------------------
#--------- cobmine (add) results from several subdirectories -----------
combine:
	(cd ../../glibk; make combine.exe)
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
q-semaphore:
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
	(csh ../farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MAIN} ${DSET} SP ${CLASS} )
script-SP850:
	(csh ../farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MAIN} ${DSET} SP850 ${CLASS} )
#-----------------------------------------------------------------------
