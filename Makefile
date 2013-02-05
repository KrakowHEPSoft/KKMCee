#-----------------------------------------------------------------------
# for 	make gtar
DATE= 2003-06-02
# for	make export
VERSION=v.4.19.b
#-----------------------------------------------------------------------
Clean:
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe
	(cd ./glibk;	make clean)
	(cd ./KK2f;	make clean)
	(cd ./jetset;	make clean)
	(cd ./bornv;	make clean)
	(cd ./dizet;	make clean)
	(cd ./ffbench;	make clean)
	(cd ./photos;	make clean)
	(cd ./tauola;	make clean)
	(cd ./spinor;	make clean)
	(cd ./spinpro;	make Clean)
	(cd ./yfspro;	make Clean)
	(cd ./KKsem;	make Clean)
	(cd ./YRprod;	make Clean)
	(cd ./NUprod;	make Clean)
	(cd ./droot;	make clean)
	(cd ./MaMar;	make clean)
	(cd ./RHadr;	make clean)
#-----------------------------------------------------------------------
export: Clean
	(cd ../; \
	gtar -cvzf KKMC-$(VERSION)-export.tar.gz  \
		KK-all/Makefile  \
		KK-all/.KK2f_defaults \
		KK-all/TitlePage \
		KK-all/HighLights \
		KK-all/HowToStart \
		KK-all/HowToMore \
		KK-all/HowToCite \
		KK-all/HowToMore \
		KK-all/InstalationNotes\
		KK-all/ReleaseNotes\
		KK-all/Acknowledgements\
		KK-all/POOPrules\
		KK-all/KK2f\
		KK-all/glibk\
		KK-all/bornv\
		KK-all/dizet\
		KK-all/jetset\
		KK-all/photos\
		KK-all/tauola\
		KK-all/KKsem\
		KK-all/ffbench\
		KK-all/droot\
		)
	(chmod -w ../KKMC-$(VERSION)-export.tar.gz)
#-----------------------------------------------------------------------
gtar: Clean
	(cd ../; gtar -cvzf KK-all-$(DATE).tar.gz  KK-all )
	(chmod -w ../KK-all-$(DATE).tar.gz)
#-----------------------------------------------------------------------
