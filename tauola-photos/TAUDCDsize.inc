      integer NMODE,NM1,NM2,NM3,NM4,NM5,NM6,NLT
      PARAMETER (NLT=2,NMODE=196,NM1=40,NM2=71,NM3=19,NM4=32,NM5=21,NM6=13)

C INFO FOR POSSIBLE CHANGES:

C it is semi automatic to increase  number of channels. 
C once one modifies/increase above declaration of   NM-i and NMODE 
C one has  to:

C ##A##  in file demo-xxx/taumain.f (or correspondintg interface file):

C 1.increase number of entries in:
C DATA KEYstrt-i, DATA NOPIK-i 

C 2. increase list of corresponding channels names, lines like
C IF(II.EQ.??) NAMES(I-NLT)='  TAU-  --> xxxxxxxixxxxxxxx   '  !  (sep 2004)
C have to be added

C ##B## in file tauola-c/TauolaStructs.h
C adjust accordingly declaration of matrix sizes used in C part of the code:  
C const int NLT=2, NMODE=196,NM1=40,NM2=71,NM3=19,NM4=32,NM5=21,NM6=13;

