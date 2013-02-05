/********************************************************************
* KorwDict.h
********************************************************************/
#ifdef __CINT__
#error KorwDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
extern "C" {
#define G__ANSIHEADER
#include "G__ci.h"
extern void G__cpp_setup_tagtableKorwDict();
extern void G__cpp_setup_inheritanceKorwDict();
extern void G__cpp_setup_typetableKorwDict();
extern void G__cpp_setup_memvarKorwDict();
extern void G__cpp_setup_globalKorwDict();
extern void G__cpp_setup_memfuncKorwDict();
extern void G__cpp_setup_funcKorwDict();
extern void G__set_cpp_environmentKorwDict();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "VLorenz.h"
#include "PartLund.h"
#include "KorEvent.h"
#include "KoralwMaker.h"
#include "BEwtMaker.h"
#include "JetAnalyzer.h"
#include "Semaph.h"
#include "ROBOL.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__KorwDictLN_TClass;
extern G__linked_taginfo G__KorwDictLN_ofstream;
extern G__linked_taginfo G__KorwDictLN_TObject;
extern G__linked_taginfo G__KorwDictLN_TNamed;
extern G__linked_taginfo G__KorwDictLN_VLorenz;
extern G__linked_taginfo G__KorwDictLN_PartLund;
extern G__linked_taginfo G__KorwDictLN_KorEvent;
extern G__linked_taginfo G__KorwDictLN_ifstream;
extern G__linked_taginfo G__KorwDictLN_KoralwMaker;
extern G__linked_taginfo G__KorwDictLN_TH1F;
extern G__linked_taginfo G__KorwDictLN_TNtuple;
extern G__linked_taginfo G__KorwDictLN_BEwtMaker;
extern G__linked_taginfo G__KorwDictLN_JetAnalyzer;
extern G__linked_taginfo G__KorwDictLN_Semaph;
extern G__linked_taginfo G__KorwDictLN_ROBOL;
