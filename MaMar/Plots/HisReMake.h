#ifndef HisReMake_H
#define HisReMake_H
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TFile.h"
#include "TH2.h"

#include "HisReMake.h"
#include "HisNorm.h"
#include "KKplot.h"


void HisReMakeKKMC(  TFile *DiskFileA, int NbMax, int NbMax2);

void HisReMakeFoam35(TFile *DiskFileA, int NbMax, int NbMax2);

#endif
