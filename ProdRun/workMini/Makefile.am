DSET = workMini

##############
KKMCee_defaults:
	(ln -s ../../SRCee/KKMCee_defaults ./)

EWtables: KKMCee_defaults
	(cd $(top_builddir); make)
	cp ./$(DSET).input ./pro.input
	($(top_builddir)/DZface/TabMainC;)

start: EWtables
	../MainMini
	
start2:
	(cd $(top_builddir); make)
	../MainPers
	