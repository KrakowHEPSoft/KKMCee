#!/bin/csh
# @ notify_user = Zbigniew.Was@cern.ch
# @ class = long
# @ preferences=(Feature=="lowio")
#######################################
unalias cp
unalias rm
unalias mv
#
#set dset=prod-kek1-161
#
# ======== WORK DIRECTORIES               
set pthfarm=/afs/cern.ch/user/w/wasm/private/amplitudy/trzynue
set pthlocal=$pthfarm
#
# ========= EXECUTION ========= 
#
cd $pthlocal
time trzynue_test.exe
#
# ========= COPYING   =========
#cp korw.output ${dset}.output
#cp korw.hst    histdump.hst
#
exit 0          
