#-----------------------------------------------------------------------
#       Template for makefile in production directory Tau
#-----------------------------------------------------------------------
MAIN=Main
DSET=Tau
CLASS=long
#-----------------------------------------------------------------------
SRCLIST= "../ ../../KK2f ../../bornv ../../jetset ../../glibk"
#-----------------------------------------------------------------------
# Interactive:  run start stop and continue
#-----------------------------------------------------------------------
test:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 time ../${MAIN}.exe			)
#-----------------------------------------------------------------------
start:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 time ../${MAIN}.exe			)&
#-----------------------------------------------------------------------
stop:
	(cp ../semaphore.stop ./semaphore )&
#-----------------------------------------------------------------------
cont:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         time ../${MAIN}.exe			)&
#-----------------------------------------------------------------------
start-debug:
	(cp ../iniseed/iniseed.0 ./iniseed	;\
	 cp ../semaphore.start   ./semaphore	;\
	 cp ./$(DSET).input ./pro.input		;\
	 xldb -I  ${SRCLIST} ../${MAIN}.exe	)&
#-----------------------------------------------------------------------
cont-debug:
	(cp semaphore  semaphore.bak		;\
	 cp pro.hst    pro.hst.bak		;\
	 cp pro.output pro.output.bak		;\
         xldb -I  ${SRCLIST} ../${MAIN}.exe	)&
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
clean:                                                     
	(rm -f *.o *.a *~ *.exe core *.err *.out)
CLEAN:  clean
	(rm  -R -f farm )
#-----------------------------------------------------------------------
