DSET = work1

##############
KKMCee_defaults:
	(ln -s ../../SRCee/KKMCee_defaults ./)

EWtables:
	(cd $(top_builddir); make)
	cp ./$(DSET).input ./pro.input
	($(top_builddir)/DZface/TabMainC;)
#	 cp -f DIZET-table1 DIZET-table1-T;)
	