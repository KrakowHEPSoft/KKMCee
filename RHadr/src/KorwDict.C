/********************************************************
* KorwDict.C
********************************************************/
#include "KorwDict.h"

#ifdef G__MEMTEST
#undef malloc
#endif


extern "C" void G__set_cpp_environmentKorwDict() {
  G__add_compiledheader("TROOT.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("VLorenz.h");
  G__add_compiledheader("PartLund.h");
  G__add_compiledheader("KorEvent.h");
  G__add_compiledheader("KoralwMaker.h");
  G__add_compiledheader("BEwtMaker.h");
  G__add_compiledheader("JetAnalyzer.h");
  G__add_compiledheader("Semaph.h");
  G__add_compiledheader("ROBOL.h");
}
int G__cpp_dllrevKorwDict() { return(51111); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* VLorenz */
static int G__VLorenz_VLorenz_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   VLorenz *p=NULL;
   if(G__getaryconstruct()) p=new VLorenz[G__getaryconstruct()];
   else                    p=new VLorenz;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_VLorenz);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_VLorenz_1_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   VLorenz *p=NULL;
      p = new VLorenz(
(double)G__double(libp->para[0]),(double)G__double(libp->para[1])
,(double)G__double(libp->para[2]),(double)G__double(libp->para[3]));
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_VLorenz);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_print_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((VLorenz*)(G__getstructoffset()))->print();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_operatoroBcB_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      {
        double& obj=((VLorenz*)(G__getstructoffset()))->operator[]((int)G__int(libp->para[0]));
         result7->ref=(long)(&obj); result7->obj.d=(double)(obj);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_DeclFileName_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((VLorenz*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_DeclFileLine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((VLorenz*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_ImplFileName_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((VLorenz*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_ImplFileLine_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((VLorenz*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_Class_Version_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((VLorenz*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_Class_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((VLorenz*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_Dictionary_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((VLorenz*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_IsA_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((VLorenz*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_ShowMembers_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((VLorenz*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__VLorenz_Streamer_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((VLorenz*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__VLorenz_VLorenz_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash)
{
   VLorenz *p;
   if(1!=libp->paran) ;
   p=new VLorenz(*(VLorenz*)G__int(libp->para[0]));
   result7->obj.i = (long)p;
   result7->ref = (long)p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_VLorenz);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__VLorenz_wAVLorenz_6_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (VLorenz *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (VLorenz *)((G__getstructoffset())+sizeof(VLorenz)*i);
   else  delete (VLorenz *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* PartLund */
static int G__PartLund_PartLund_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   PartLund *p=NULL;
   if(G__getaryconstruct()) p=new PartLund[G__getaryconstruct()];
   else                    p=new PartLund;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_PartLund);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_PartLund_1_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   PartLund *p=NULL;
      p = new PartLund(
(int)G__int(libp->para[0]),(int)G__int(libp->para[1])
,(int)G__int(libp->para[2]),(int)G__int(libp->para[3])
,(int)G__int(libp->para[4]),(int)G__int(libp->para[5])
,(double)G__double(libp->para[6]),(double)G__double(libp->para[7])
,(double)G__double(libp->para[8]),(double)G__double(libp->para[9])
,(double)G__double(libp->para[10]),(double)G__double(libp->para[11])
,(double)G__double(libp->para[12]),(double)G__double(libp->para[13])
,(double)G__double(libp->para[14]),(double)G__double(libp->para[15]));
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_PartLund);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_print_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   switch(libp->paran) {
   case 1:
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->print((int)G__int(libp->para[0]));
      break;
   case 0:
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->print();
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_ListPrint_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->ListPrint();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_DeclFileName_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((PartLund*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_DeclFileLine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((PartLund*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_ImplFileName_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((PartLund*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_ImplFileLine_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((PartLund*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_Class_Version_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((PartLund*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_Class_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((PartLund*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_Dictionary_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_IsA_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((PartLund*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_ShowMembers_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__PartLund_Streamer_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((PartLund*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__PartLund_PartLund_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash)
{
   PartLund *p;
   if(1!=libp->paran) ;
   p=new PartLund(*(PartLund*)G__int(libp->para[0]));
   result7->obj.i = (long)p;
   result7->ref = (long)p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_PartLund);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__PartLund_wAPartLund_6_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (PartLund *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (PartLund *)((G__getstructoffset())+sizeof(PartLund)*i);
   else  delete (PartLund *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* KorEvent */
static int G__KorEvent_KorEvent_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   KorEvent *p=NULL;
   if(G__getaryconstruct()) p=new KorEvent[G__getaryconstruct()];
   else                    p=new KorEvent;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_KorEvent);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_Print_2_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->Print((long)G__int(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_GetPartonMass_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->GetPartonMass((double*)G__int(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_GetPartonAngle_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letdouble(result7,100,(double)((KorEvent*)(G__getstructoffset()))->GetPartonAngle((double*)G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_GetJetMass_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->GetJetMass((double*)G__int(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_GetJetAngle_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letdouble(result7,100,(double)((KorEvent*)(G__getstructoffset()))->GetJetAngle((double*)G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_DeclFileName_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((KorEvent*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_DeclFileLine_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((KorEvent*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_ImplFileName_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((KorEvent*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_ImplFileLine_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((KorEvent*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_Class_Version_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((KorEvent*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_Class_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((KorEvent*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_Dictionary_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_IsA_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((KorEvent*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_ShowMembers_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KorEvent_Streamer_6_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KorEvent*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__KorEvent_KorEvent_7_1(G__value *result7,char *funcname,struct G__param *libp,int hash)
{
   KorEvent *p;
   if(1!=libp->paran) ;
   p=new KorEvent(*(KorEvent*)G__int(libp->para[0]));
   result7->obj.i = (long)p;
   result7->ref = (long)p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_KorEvent);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__KorEvent_wAKorEvent_8_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (KorEvent *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (KorEvent *)((G__getstructoffset())+sizeof(KorEvent)*i);
   else  delete (KorEvent *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* KoralwMaker */
static int G__KoralwMaker_KoralwMaker_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   KoralwMaker *p=NULL;
   if(G__getaryconstruct()) p=new KoralwMaker[G__getaryconstruct()];
   else                    p=new KoralwMaker;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_ReadData_2_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->ReadData(libp->para[0].ref?*(long*)libp->para[0].ref:G__Mlong(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Initialize_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->Initialize();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Finalize_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->Finalize();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Generate_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->Generate(*(KorEvent*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_JetDefine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->JetDefine(*(KorEvent*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_LuGive_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->LuGive((char*)G__int(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_LuList_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->LuList(*(KorEvent*)libp->para[0].ref,(long)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_FortOpen_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->FortOpen(libp->para[0].ref?*(long*)libp->para[0].ref:G__Mlong(libp->para[0]),(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_FortClose_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->FortClose(libp->para[0].ref?*(long*)libp->para[0].ref:G__Mlong(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_VarRan_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->VarRan((double*)G__int(libp->para[0]),libp->para[1].ref?*(long*)libp->para[1].ref:G__Mlong(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_DeclFileName_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((KoralwMaker*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_DeclFileLine_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((KoralwMaker*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_ImplFileName_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((KoralwMaker*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_ImplFileLine_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((KoralwMaker*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Class_Version_6_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((KoralwMaker*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Class_7_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((KoralwMaker*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Dictionary_8_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_IsA_9_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((KoralwMaker*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_ShowMembers_0_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__KoralwMaker_Streamer_1_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((KoralwMaker*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__KoralwMaker_KoralwMaker_2_2(G__value *result7,char *funcname,struct G__param *libp,int hash)
{
   KoralwMaker *p;
   if(1!=libp->paran) ;
   p=new KoralwMaker(*(KoralwMaker*)G__int(libp->para[0]));
   result7->obj.i = (long)p;
   result7->ref = (long)p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__KoralwMaker_wAKoralwMaker_3_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (KoralwMaker *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (KoralwMaker *)((G__getstructoffset())+sizeof(KoralwMaker)*i);
   else  delete (KoralwMaker *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* BEwtMaker */
static int G__BEwtMaker_BEwtMaker_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   BEwtMaker *p=NULL;
   if(G__getaryconstruct()) p=new BEwtMaker[G__getaryconstruct()];
   else                    p=new BEwtMaker;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_SetModel_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((BEwtMaker*)(G__getstructoffset()))->SetModel((double)G__double(libp->para[0]),(long)G__int(libp->para[1])
,(double)G__double(libp->para[2]),(double)G__double(libp->para[3])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_SetRenorm_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((BEwtMaker*)(G__getstructoffset()))->SetRenorm((double)G__double(libp->para[0]),(double)G__double(libp->para[1])
,(double)G__double(libp->para[2]),(double)G__double(libp->para[3])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_MakeWeight_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((BEwtMaker*)(G__getstructoffset()))->MakeWeight(*(KorEvent*)libp->para[0].ref,libp->para[1].ref?*(long*)libp->para[1].ref:G__Mlong(libp->para[1])
,libp->para[2].ref?*(double*)libp->para[2].ref:G__Mdouble(libp->para[2]),libp->para[3].ref?*(double*)libp->para[3].ref:G__Mdouble(libp->para[3]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_BookLSP_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((BEwtMaker*)(G__getstructoffset()))->BookLSP(*(KorEvent*)libp->para[0].ref,(int)G__int(libp->para[1])
,(int)G__int(libp->para[2]),(double)G__double(libp->para[3])
,(double)G__double(libp->para[4]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_Q2pair_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letdouble(result7,100,(double)((BEwtMaker*)(G__getstructoffset()))->Q2pair((PartLund*)G__int(libp->para[0]),(PartLund*)G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_DeclFileName_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((BEwtMaker*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_DeclFileLine_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((BEwtMaker*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_ImplFileName_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((BEwtMaker*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_ImplFileLine_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((BEwtMaker*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_Class_Version_6_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((BEwtMaker*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_Class_7_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((BEwtMaker*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_Dictionary_8_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((BEwtMaker*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_IsA_9_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((BEwtMaker*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_ShowMembers_0_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((BEwtMaker*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BEwtMaker_Streamer_1_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((BEwtMaker*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__BEwtMaker_wABEwtMaker_2_2(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (BEwtMaker *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (BEwtMaker *)((G__getstructoffset())+sizeof(BEwtMaker)*i);
   else  delete (BEwtMaker *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* JetAnalyzer */
static int G__JetAnalyzer_JetAnalyzer_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   JetAnalyzer *p=NULL;
   if(G__getaryconstruct()) p=new JetAnalyzer[G__getaryconstruct()];
   else                    p=new JetAnalyzer;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_Book_2_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((JetAnalyzer*)(G__getstructoffset()))->Book(*(KorEvent*)libp->para[0].ref,(long)G__int(libp->para[1])
,(double)G__double(libp->para[2]),(double)G__double(libp->para[3]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_DeclFileName_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((JetAnalyzer*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_DeclFileLine_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((JetAnalyzer*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_ImplFileName_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((JetAnalyzer*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_ImplFileLine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((JetAnalyzer*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_Class_Version_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((JetAnalyzer*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_Class_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((JetAnalyzer*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_Dictionary_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((JetAnalyzer*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_IsA_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((JetAnalyzer*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_ShowMembers_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((JetAnalyzer*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JetAnalyzer_Streamer_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((JetAnalyzer*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__JetAnalyzer_JetAnalyzer_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash)
{
   JetAnalyzer *p;
   if(1!=libp->paran) ;
   p=new JetAnalyzer(*(JetAnalyzer*)G__int(libp->para[0]));
   result7->obj.i = (long)p;
   result7->ref = (long)p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__JetAnalyzer_wAJetAnalyzer_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (JetAnalyzer *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (JetAnalyzer *)((G__getstructoffset())+sizeof(JetAnalyzer)*i);
   else  delete (JetAnalyzer *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* Semaph */
static int G__Semaph_Semaph_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   Semaph *p=NULL;
   if(G__getaryconstruct()) p=new Semaph[G__getaryconstruct()];
   else                    p=new Semaph;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_Semaph);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_Initialize_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((Semaph*)(G__getstructoffset()))->Initialize(*(TString*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_ReadStatus_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((Semaph*)(G__getstructoffset()))->ReadStatus(*(TString*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_DeclFileName_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((Semaph*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_DeclFileLine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((Semaph*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_ImplFileName_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((Semaph*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_ImplFileLine_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((Semaph*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_Class_Version_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((Semaph*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_Class_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((Semaph*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_Dictionary_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((Semaph*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_IsA_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((Semaph*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_ShowMembers_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((Semaph*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Semaph_Streamer_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((Semaph*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__Semaph_wASemaph_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (Semaph *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (Semaph *)((G__getstructoffset())+sizeof(Semaph)*i);
   else  delete (Semaph *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* ROBOL */
static int G__ROBOL_ROBOL_0_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   ROBOL *p=NULL;
   if(G__getaryconstruct()) p=new ROBOL[G__getaryconstruct()];
   else                    p=new ROBOL;
      result7->obj.i = (long)p;
      result7->ref = (long)p;
      result7->type = 'u';
      result7->tagnum = G__get_linked_tagnum(&G__KorwDictLN_ROBOL);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Initialize_2_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->Initialize(libp->para[0].ref?*(long*)libp->para[0].ref:G__Mlong(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Production_3_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->Production(libp->para[0].ref?*(long*)libp->para[0].ref:G__Mlong(libp->para[0]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Finalize_4_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->Finalize();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_DeclFileName_5_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((ROBOL*)(G__getstructoffset()))->DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_DeclFileLine_6_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((ROBOL*)(G__getstructoffset()))->DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_ImplFileName_7_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,67,(long)((ROBOL*)(G__getstructoffset()))->ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_ImplFileLine_8_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,105,(long)((ROBOL*)(G__getstructoffset()))->ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Class_Version_9_0(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,115,(long)((ROBOL*)(G__getstructoffset()))->Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Class_0_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((ROBOL*)(G__getstructoffset()))->Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Dictionary_1_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->Dictionary();
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_IsA_2_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__letint(result7,85,(long)((ROBOL*)(G__getstructoffset()))->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_ShowMembers_3_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->ShowMembers(*(TMemberInspector*)libp->para[0].ref,(char*)G__int(libp->para[1]));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ROBOL_Streamer_4_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
      G__setnull(result7);
      ((ROBOL*)(G__getstructoffset()))->Streamer(*(TBuffer*)libp->para[0].ref);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
static int G__ROBOL_wAROBOL_5_1(G__value *result7,char *funcname,struct G__param *libp,int hash) {
   if(G__getaryconstruct())
     if(G__PVOID==G__getgvp())
       delete[] (ROBOL *)(G__getstructoffset());
     else
       for(int i=G__getaryconstruct()-1;i>=0;i--)
         delete (ROBOL *)((G__getstructoffset())+sizeof(ROBOL)*i);
   else  delete (ROBOL *)(G__getstructoffset());
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* VLorenz */

/* PartLund */

/* KorEvent */

/* KoralwMaker */

/* BEwtMaker */

/* JetAnalyzer */

/* Semaph */

/* ROBOL */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncKorwDict {
 public:
  G__Sizep2memfuncKorwDict() {p=&G__Sizep2memfuncKorwDict::sizep2memfunc;}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncKorwDict::*p)();
};

size_t G__get_sizep2memfuncKorwDict()
{
  G__Sizep2memfuncKorwDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceKorwDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_VLorenz))) {
     VLorenz *G__Lderived;
     G__Lderived=(VLorenz*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_VLorenz),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_PartLund))) {
     PartLund *G__Lderived;
     G__Lderived=(PartLund*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_PartLund),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_KorEvent))) {
     KorEvent *G__Lderived;
     G__Lderived=(KorEvent*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_KorEvent),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker))) {
     KoralwMaker *G__Lderived;
     G__Lderived=(KoralwMaker*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),G__get_linked_tagnum(&G__KorwDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker))) {
     BEwtMaker *G__Lderived;
     G__Lderived=(BEwtMaker*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),G__get_linked_tagnum(&G__KorwDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer))) {
     JetAnalyzer *G__Lderived;
     G__Lderived=(JetAnalyzer*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),G__get_linked_tagnum(&G__KorwDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_Semaph))) {
     Semaph *G__Lderived;
     G__Lderived=(Semaph*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_Semaph),G__get_linked_tagnum(&G__KorwDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_Semaph),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__KorwDictLN_ROBOL))) {
     ROBOL *G__Lderived;
     G__Lderived=(ROBOL*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_ROBOL),G__get_linked_tagnum(&G__KorwDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__KorwDictLN_ROBOL),G__get_linked_tagnum(&G__KorwDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableKorwDict() {

   /* Setting up typedef entry */
   G__search_typename2("Char_t",99,-1,0,
-1);
   G__setnewtype(-1,"Signed Character 1 byte",0);
   G__search_typename2("UChar_t",98,-1,0,
-1);
   G__setnewtype(-1,"Unsigned Character 1 byte",0);
   G__search_typename2("Short_t",115,-1,0,
-1);
   G__setnewtype(-1,"Signed Short integer 2 bytes",0);
   G__search_typename2("UShort_t",114,-1,0,
-1);
   G__setnewtype(-1,"Unsigned Short integer 2 bytes",0);
   G__search_typename2("Int_t",105,-1,0,
-1);
   G__setnewtype(-1,"Signed integer 4 bytes",0);
   G__search_typename2("UInt_t",104,-1,0,
-1);
   G__setnewtype(-1,"Unsigned integer 4 bytes",0);
   G__search_typename2("Seek_t",105,-1,0,
-1);
   G__setnewtype(-1,"File pointer",0);
   G__search_typename2("Long_t",108,-1,0,
-1);
   G__setnewtype(-1,"Signed long integer 8 bytes",0);
   G__search_typename2("ULong_t",107,-1,0,
-1);
   G__setnewtype(-1,"Unsigned long integer 8 bytes",0);
   G__search_typename2("Float_t",102,-1,0,
-1);
   G__setnewtype(-1,"Float 4 bytes",0);
   G__search_typename2("Double_t",100,-1,0,
-1);
   G__setnewtype(-1,"Float 8 bytes",0);
   G__search_typename2("Text_t",99,-1,0,
-1);
   G__setnewtype(-1,"General string",0);
   G__search_typename2("Bool_t",98,-1,0,
-1);
   G__setnewtype(-1,"Boolean (0=false, 1=true)",0);
   G__search_typename2("Byte_t",98,-1,0,
-1);
   G__setnewtype(-1,"Byte (8 bits)",0);
   G__search_typename2("Version_t",115,-1,0,
-1);
   G__setnewtype(-1,"Class version identifier",0);
   G__search_typename2("Option_t",99,-1,0,
-1);
   G__setnewtype(-1,"Option string",0);
   G__search_typename2("Ssiz_t",105,-1,0,
-1);
   G__setnewtype(-1,"String size",0);
   G__search_typename2("Real_t",102,-1,0,
-1);
   G__setnewtype(-1,"TVector and TMatrix element type",0);
   G__search_typename2("VoidFuncPtr_t",89,-1,0,
-1);
   G__setnewtype(-1,"pointer to void function",0);
   G__search_typename2("FreeHookFun_t",89,-1,0,
-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("ReAllocFun_t",81,-1,2,
-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("ReAllocCFun_t",81,-1,2,
-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("Axis_t",102,-1,0,
-1);
   G__setnewtype(-1,"Axis values type",0);
   G__search_typename2("Stat_t",100,-1,0,
-1);
   G__setnewtype(-1,"Statistics type",0);
   G__search_typename2("Font_t",115,-1,0,
-1);
   G__setnewtype(-1,"Font number",0);
   G__search_typename2("Style_t",115,-1,0,
-1);
   G__setnewtype(-1,"Style number",0);
   G__search_typename2("Marker_t",115,-1,0,
-1);
   G__setnewtype(-1,"Marker number",0);
   G__search_typename2("Width_t",115,-1,0,
-1);
   G__setnewtype(-1,"Line width",0);
   G__search_typename2("Color_t",115,-1,0,
-1);
   G__setnewtype(-1,"Color number",0);
   G__search_typename2("SCoord_t",115,-1,0,
-1);
   G__setnewtype(-1,"Screen coordinates",0);
   G__search_typename2("Coord_t",102,-1,0,
-1);
   G__setnewtype(-1,"Pad world coordinates",0);
   G__search_typename2("Angle_t",102,-1,0,
-1);
   G__setnewtype(-1,"Graphics angle",0);
   G__search_typename2("Size_t",102,-1,0,
-1);
   G__setnewtype(-1,"Attribute size",0);
   G__search_typename2("Double_t (*)(Double_t*, Double_t*)",81,-1,0,
-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* VLorenz */
static void G__setup_memvarVLorenz(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_VLorenz));
   { VLorenz *p; p=(VLorenz*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->m_comp)-(long)(p)),100,0,0,-1,-1,-1,1,"m_comp[4]=",0,"4-momentum, comp[0] is energy");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* PartLund */
static void G__setup_memvarPartLund(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_PartLund));
   { PartLund *p; p=(PartLund*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->previous)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,-1,1,"previous=",0,"pointer for constructing lists");
   G__memvar_setup((void*)((long)(&p->next)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,-1,1,"next=",0,"pointer for constructing lists");
   G__memvar_setup((void*)((long)(&p->m_lserial)-(long)(p)),105,0,0,-1,-1,-1,1,"m_lserial=",0,"Lund serial number");
   G__memvar_setup((void*)((long)(&p->m_kstatus)-(long)(p)),105,0,0,-1,-1,-1,1,"m_kstatus=",0,"status");
   G__memvar_setup((void*)((long)(&p->m_kflavor)-(long)(p)),105,0,0,-1,-1,-1,1,"m_kflavor=",0,"flavour");
   G__memvar_setup((void*)((long)(&p->m_kparent)-(long)(p)),105,0,0,-1,-1,-1,1,"m_kparent=",0,"parent");
   G__memvar_setup((void*)((long)(&p->m_kFirstChild)-(long)(p)),105,0,0,-1,-1,-1,1,"m_kFirstChild=",0,"first child");
   G__memvar_setup((void*)((long)(&p->m_kLastChild)-(long)(p)),105,0,0,-1,-1,-1,1,"m_kLastChild=",0,"last  child");
   G__memvar_setup((void*)((long)(&p->m_pmom)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_pmom=",0,"4-momentum, pmom[0] is energy [GeV]");
   G__memvar_setup((void*)((long)(&p->m_pmass)-(long)(p)),100,0,0,-1,-1,-1,1,"m_pmass=",0,"mass [GeV]");
   G__memvar_setup((void*)((long)(&p->m_Vertex)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_Vertex=",0,"vertex position [mm], V[0] is time in Lab");
   G__memvar_setup((void*)((long)(&p->m_LifeTime)-(long)(p)),100,0,0,-1,-1,-1,1,"m_LifeTime=",0,"lifetime [mm/sec]");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* KorEvent */
static void G__setup_memvarKorEvent(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_KorEvent));
   { KorEvent *p; p=(KorEvent*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->m_nphot)-(long)(p)),105,0,0,-1,-1,-1,1,"m_nphot=",0,"number of photons");
   G__memvar_setup((void*)((long)(&p->m_photmom)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_photmom[100]=",0,"list of photon momenta");
   G__memvar_setup((void*)((long)(&p->m_wminus)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_wminus=",0,"4-momentum of W-");
   G__memvar_setup((void*)((long)(&p->m_wplus)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_wplus=",0,"4-momentum of W+");
   G__memvar_setup((void*)((long)(&p->m_ferm1)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_ferm1=",0,"4-momentum of fermion 1");
   G__memvar_setup((void*)((long)(&p->m_ferm2)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_ferm2=",0,"4-momentum of fermion 2");
   G__memvar_setup((void*)((long)(&p->m_ferm3)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_ferm3=",0,"4-momentum of fermion 3");
   G__memvar_setup((void*)((long)(&p->m_ferm4)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,-1,1,"m_ferm4=",0,"4-momentum of fermion 4");
   G__memvar_setup((void*)((long)(&p->m_npart)-(long)(p)),105,0,0,-1,-1,-1,1,"m_npart=",0,"total number of particles in part[]");
   G__memvar_setup((void*)((long)(&p->m_part)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,-1,1,"m_part[20000]=",0,"all particles in Lund-like format");
   G__memvar_setup((void*)((long)(&p->m_njet)-(long)(p)),105,0,0,-1,-1,-1,1,"m_njet=",0,"total number of jets");
   G__memvar_setup((void*)((long)(&p->m_jet)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,-1,1,"m_jet[200]=",0,"jets in Lund-like format");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* KoralwMaker */
static void G__setup_memvarKoralwMaker(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker));
   { KoralwMaker *p; p=(KoralwMaker*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->NevTot)-(long)(p)),108,0,0,-1,-1,-1,1,"NevTot=",0,"total numer of events to be generated");
   G__memvar_setup((void*)((long)(&p->m_EvenCounter)-(long)(p)),108,0,0,-1,-1,-1,1,"m_EvenCounter=",0,"event serial counter");
   G__memvar_setup((void*)((long)(&p->m_PrintLimit)-(long)(p)),108,0,0,-1,-1,-1,1,"m_PrintLimit=",0,"print limit for debug");
   G__memvar_setup((void*)((long)(&p->xpar)-(long)(p)),100,0,0,-1,-1,-1,1,"xpar[10001]=",0,"xpar input/output aprameters of koralw");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* BEwtMaker */
static void G__setup_memvarBEwtMaker(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker));
   { BEwtMaker *p; p=(BEwtMaker*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->BEoutfile)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_ofstream),-1,-1,1,"BEoutfile=",0,"! output file for information printouts");
   G__memvar_setup((void*)((long)(&p->BEievent)-(long)(p)),108,0,0,-1,-1,-1,1,"BEievent=",0,"Serial number of event (for printouts)");
   G__memvar_setup((void*)((long)(&p->BEprintlevel)-(long)(p)),108,0,0,-1,-1,-1,1,"BEprintlevel=",0,"level of information printout");
   G__memvar_setup((void*)((long)(&p->BEprintlast)-(long)(p)),108,0,0,-1,-1,-1,1,"BEprintlast=",0,"number of events printed");
   G__memvar_setup((void*)((long)(&p->BErange)-(long)(p)),100,0,0,-1,-1,-1,1,"BErange=",0,"Range, maximum sqrt(Q2) range for a single pair [GeV]");
   G__memvar_setup((void*)((long)(&p->BEFuncType)-(long)(p)),108,0,0,-1,-1,-1,1,"BEFuncType=",0,"Function type =1,2 for Gausian,exponential");
   G__memvar_setup((void*)((long)(&p->BEpp)-(long)(p)),100,0,0,-1,-1,-1,1,"BEpp=",0,"p parameter");
   G__memvar_setup((void*)((long)(&p->BEradius)-(long)(p)),100,0,0,-1,-1,-1,1,"BEradius=",0,"R radius in fermi units");
   G__memvar_setup((void*)((long)(&p->BElambda)-(long)(p)),100,0,0,-1,-1,-1,1,"BElambda=",0,"lambda parameter for wt  (for wt renormalization)");
   G__memvar_setup((void*)((long)(&p->BEavewt)-(long)(p)),100,0,0,-1,-1,-1,1,"BEavewt=",0,"average weight   for wt  (for wt renormalization)");
   G__memvar_setup((void*)((long)(&p->BElambda2)-(long)(p)),100,0,0,-1,-1,-1,1,"BElambda2=",0,"lambda parameter for wt2 (for wt2 renormalization)");
   G__memvar_setup((void*)((long)(&p->BEavewt2)-(long)(p)),100,0,0,-1,-1,-1,1,"BEavewt2=",0,"average weight   for wt2 (for wt2 renormalization)");
   G__memvar_setup((void*)((long)(&p->hst_nclu)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_nclu=",0,"histo cluster multiplicity");
   G__memvar_setup((void*)((long)(&p->hst_unQ2)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_unQ2=",0,"histo Q2 of like-sign-pion-pair");
   G__memvar_setup((void*)((long)(&p->hst_wtQ2)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_wtQ2=",0,"histo Q2 of like-sign-pion-pair, wt");
   G__memvar_setup((void*)((long)(&p->hst_wt2Q2)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_wt2Q2=",0,"histo Q2 of like-sign-pion-pair, wt2");
   G__memvar_setup((void*)((long)(&p->hst_unQ3)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_unQ3=",0,"histo Q2 of like-sign-pion-triplet");
   G__memvar_setup((void*)((long)(&p->hst_wtQ3)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_wtQ3=",0,"histo Q2 of like-sign-pion-triplet, wt");
   G__memvar_setup((void*)((long)(&p->hst_wt2Q3)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_wt2Q3=",0,"histo Q2 of like-sign-pion-triplet, wt2");
   G__memvar_setup((void*)((long)(&p->ctuple2_counter)-(long)(p)),100,0,0,-1,-1,-1,1,"ctuple2_counter=",0,"counter of entries in ntuple");
   G__memvar_setup((void*)((long)(&p->ctuple2_max)-(long)(p)),100,0,0,-1,-1,-1,1,"ctuple2_max=",0,"maximum entries in ntuple");
   G__memvar_setup((void*)((long)(&p->ctuple2)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TNtuple),-1,-1,1,"ctuple2=",0,"ntuple Q2,wt,pn");
   G__memvar_setup((void*)((long)(&p->ctuple3_counter)-(long)(p)),100,0,0,-1,-1,-1,1,"ctuple3_counter=",0,"counter of entries in ntuple");
   G__memvar_setup((void*)((long)(&p->ctuple3_max)-(long)(p)),100,0,0,-1,-1,-1,1,"ctuple3_max=",0,"maximum entries in ntuple");
   G__memvar_setup((void*)((long)(&p->ctuple3)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TNtuple),-1,-1,1,"ctuple3=",0,"ntuple Q3,wt,pn");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* JetAnalyzer */
static void G__setup_memvarJetAnalyzer(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer));
   { JetAnalyzer *p; p=(JetAnalyzer*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->hst_pene)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_pene=",0,"parton minimum energy");
   G__memvar_setup((void*)((long)(&p->hst_jene)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_jene=",0,"jet minimum energy");
   G__memvar_setup((void*)((long)(&p->jtuple)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TNtuple),-1,-1,1,"jtuple=",0,"big ntuple with jet masses etc");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* Semaph */
static void G__setup_memvarSemaph(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_Semaph));
   { Semaph *p; p=(Semaph*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->InfileSema)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_ifstream),-1,-1,1,"InfileSema=",0,"! Input semaphore file");
   G__memvar_setup((void*)((long)(&p->InfileSeed)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_ifstream),-1,-1,1,"InfileSeed=",0,"! Input seed file");
   G__memvar_setup((void*)((long)(&p->OufileSema)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_ofstream),-1,-1,1,"OufileSema=",0,"! Output semaphore file");
   G__memvar_setup((void*)((long)(&p->ijklin)-(long)(p)),108,0,0,-1,-1,-1,1,"ijklin=",0,"for ranmar");
   G__memvar_setup((void*)((long)(&p->ntotin)-(long)(p)),108,0,0,-1,-1,-1,1,"ntotin=",0,"for ranmar");
   G__memvar_setup((void*)((long)(&p->ntot2n)-(long)(p)),108,0,0,-1,-1,-1,1,"ntot2n=",0,"for ranmar");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* ROBOL */
static void G__setup_memvarROBOL(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__KorwDictLN_ROBOL));
   { ROBOL *p; p=(ROBOL*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->KoralW)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),-1,-1,1,"KoralW=",0,"define MC generator maker (manager)");
   G__memvar_setup((void*)((long)(&p->BEfactory)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),-1,-1,1,"BEfactory=",0,"define BEwt maker (calculator)");
   G__memvar_setup((void*)((long)(&p->J4factory)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),-1,-1,1,"J4factory=",0,"define 4jet maker (analyzer)");
   G__memvar_setup((void*)((long)(&p->current_event)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__KorwDictLN_KorEvent),-1,-1,1,"current_event=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->hst_BEwt)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_BEwt=",0,"BE weight before rejection");
   G__memvar_setup((void*)((long)(&p->hst_BEwtAR)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TH1F),-1,-1,1,"hst_BEwtAR=",0,"BE weight After  Rejection");
   G__memvar_setup((void*)((long)(&p->KeyRej)-(long)(p)),108,0,0,-1,-1,-1,1,"KeyRej=",0,"rejection key (=0 for rejection ON)");
   G__memvar_setup((void*)((long)(&p->WtMax)-(long)(p)),100,0,0,-1,-1,-1,1,"WtMax=",0,"maximum rejection weight");
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarKorwDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncVLorenz(void) {
   /* VLorenz */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_VLorenz));
   G__memfunc_setup("VLorenz",720,G__VLorenz_VLorenz_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("VLorenz",720,G__VLorenz_VLorenz_1_0,105,G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,0,4,1,1,0,
"d - - 0 - en d - - 0 - px "
"d - - 0 - py d - - 0 - pz",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("print",557,G__VLorenz_print_3_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("operator[]",1060,G__VLorenz_operatoroBcB_4_0,100,-1,-1,1,1,1,1,0,"i - - 0 - index",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__VLorenz_DeclFileName_5_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__VLorenz_DeclFileLine_6_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__VLorenz_ImplFileName_7_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__VLorenz_ImplFileLine_8_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__VLorenz_Class_Version_9_0,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__VLorenz_Class_0_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__VLorenz_Dictionary_1_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__VLorenz_IsA_2_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__VLorenz_ShowMembers_3_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__VLorenz_Streamer_4_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic copy constructor
   G__memfunc_setup("VLorenz",720,G__VLorenz_VLorenz_5_1,(int)('i'),G__get_linked_tagnum(&G__KorwDictLN_VLorenz),-1,0,1,1,1,0,"u 'VLorenz' - 1 - -",(char*)NULL,(void*)NULL,0);
   // automatic destructor
   G__memfunc_setup("~VLorenz",846,G__VLorenz_wAVLorenz_6_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncPartLund(void) {
   /* PartLund */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_PartLund));
   G__memfunc_setup("PartLund",810,G__PartLund_PartLund_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("PartLund",810,G__PartLund_PartLund_1_0,105,G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,0,16,1,1,0,
"i - - 0 - lserial i - - 0 - kstatus "
"i - - 0 - kflavor i - - 0 - kparent "
"i - - 0 - kFirstChild i - - 0 - kLastChild "
"d - - 0 - px d - - 0 - py "
"d - - 0 - pz d - - 0 - en "
"d - - 0 - pmass d - - 0 - Vx "
"d - - 0 - Vy d - - 0 - Vz "
"d - - 0 - Vt d - - 0 - LifeTime",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("print",557,G__PartLund_print_3_0,121,-1,-1,0,1,1,1,0,"i - - 0 1 mode",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ListPrint",937,G__PartLund_ListPrint_4_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__PartLund_DeclFileName_5_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__PartLund_DeclFileLine_6_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__PartLund_ImplFileName_7_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__PartLund_ImplFileLine_8_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__PartLund_Class_Version_9_0,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__PartLund_Class_0_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__PartLund_Dictionary_1_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__PartLund_IsA_2_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__PartLund_ShowMembers_3_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__PartLund_Streamer_4_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic copy constructor
   G__memfunc_setup("PartLund",810,G__PartLund_PartLund_5_1,(int)('i'),G__get_linked_tagnum(&G__KorwDictLN_PartLund),-1,0,1,1,1,0,"u 'PartLund' - 1 - -",(char*)NULL,(void*)NULL,0);
   // automatic destructor
   G__memfunc_setup("~PartLund",936,G__PartLund_wAPartLund_6_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncKorEvent(void) {
   /* KorEvent */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_KorEvent));
   G__memfunc_setup("KorEvent",814,G__KorEvent_KorEvent_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_KorEvent),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Print",525,G__KorEvent_Print_2_0,121,-1,-1,0,1,1,1,0,"l - - 0 - level",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("GetPartonMass",1320,G__KorEvent_GetPartonMass_3_0,121,-1,-1,0,1,1,1,0,"D - - 0 - m",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("GetPartonAngle",1403,G__KorEvent_GetPartonAngle_4_0,100,-1,-1,0,1,1,1,0,"D - - 0 - angle",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("GetJetMass",983,G__KorEvent_GetJetMass_5_0,121,-1,-1,0,1,1,1,0,"D - - 0 - m",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("GetJetAngle",1066,G__KorEvent_GetJetAngle_6_0,100,-1,-1,0,1,1,1,0,"D - - 0 - angle",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__KorEvent_DeclFileName_7_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__KorEvent_DeclFileLine_8_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__KorEvent_ImplFileName_9_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__KorEvent_ImplFileLine_0_1,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__KorEvent_Class_Version_1_1,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__KorEvent_Class_2_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__KorEvent_Dictionary_3_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__KorEvent_IsA_4_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__KorEvent_ShowMembers_5_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__KorEvent_Streamer_6_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic copy constructor
   G__memfunc_setup("KorEvent",814,G__KorEvent_KorEvent_7_1,(int)('i'),G__get_linked_tagnum(&G__KorwDictLN_KorEvent),-1,0,1,1,1,0,"u 'KorEvent' - 1 - -",(char*)NULL,(void*)NULL,0);
   // automatic destructor
   G__memfunc_setup("~KorEvent",940,G__KorEvent_wAKorEvent_8_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncKoralwMaker(void) {
   /* KoralwMaker */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker));
   G__memfunc_setup("KoralwMaker",1120,G__KoralwMaker_KoralwMaker_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ReadData",758,G__KoralwMaker_ReadData_2_0,121,-1,-1,0,1,1,1,0,"l - - 1 - ntot",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Initialize",1042,G__KoralwMaker_Initialize_3_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Finalize",818,G__KoralwMaker_Finalize_4_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Generate",811,G__KoralwMaker_Generate_5_0,121,-1,-1,0,1,1,1,0,"u 'KorEvent' - 1 - E",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("JetDefine",878,G__KoralwMaker_JetDefine_6_0,121,-1,-1,0,1,1,1,0,"u 'KorEvent' - 1 - event",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("LuGive",588,G__KoralwMaker_LuGive_7_0,121,-1,-1,0,1,1,1,0,"C - - 0 - directive",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("LuList",605,G__KoralwMaker_LuList_8_0,121,-1,-1,0,2,1,1,0,
"u 'KorEvent' - 1 - event l - - 0 - Level",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("FortOpen",813,G__KoralwMaker_FortOpen_9_0,121,-1,-1,0,2,1,1,0,
"l - - 1 - lunit C - - 0 - filename",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("FortClose",913,G__KoralwMaker_FortClose_0_1,121,-1,-1,0,1,1,1,0,"l - - 1 - lunit",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("VarRan",586,G__KoralwMaker_VarRan_1_1,121,-1,-1,0,2,1,1,0,
"D - - 0 - rn l - - 1 - n",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__KoralwMaker_DeclFileName_2_1,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__KoralwMaker_DeclFileLine_3_1,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__KoralwMaker_ImplFileName_4_1,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__KoralwMaker_ImplFileLine_5_1,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__KoralwMaker_Class_Version_6_1,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__KoralwMaker_Class_7_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__KoralwMaker_Dictionary_8_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__KoralwMaker_IsA_9_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__KoralwMaker_ShowMembers_0_2,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__KoralwMaker_Streamer_1_2,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic copy constructor
   G__memfunc_setup("KoralwMaker",1120,G__KoralwMaker_KoralwMaker_2_2,(int)('i'),G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),-1,0,1,1,1,0,"u 'KoralwMaker' - 1 - -",(char*)NULL,(void*)NULL,0);
   // automatic destructor
   G__memfunc_setup("~KoralwMaker",1246,G__KoralwMaker_wAKoralwMaker_3_2,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncBEwtMaker(void) {
   /* BEwtMaker */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker));
   G__memfunc_setup("BEwtMaker",866,G__BEwtMaker_BEwtMaker_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("BEwtMaker",866,(G__InterfaceMethod)NULL,105,G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),-1,0,1,1,4,0,"u 'BEwtMaker' - 1 - org",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("SetModel",797,G__BEwtMaker_SetModel_3_0,105,G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),-1,0,4,1,1,0,
"d - - 0 - range l - - 0 - FuncType "
"d - - 0 - pp d - - 0 - radius",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("SetRenorm",927,G__BEwtMaker_SetRenorm_4_0,105,G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),-1,0,4,1,1,0,
"d - - 0 - lambda d - - 0 - avewt "
"d - - 0 - lambda2 d - - 0 - avewt2",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("CorFun1",638,(G__InterfaceMethod)NULL,100,-1,-1,0,4,1,4,0,
"d - - 0 - Q2 i - - 0 - mult "
"d - - 0 - Rf d - - 0 - pp",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("CorFun2",639,(G__InterfaceMethod)NULL,100,-1,-1,0,4,1,4,0,
"d - - 0 - Q2 i - - 0 - mult "
"d - - 0 - Rf d - - 0 - pp",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("WTcluster",941,(G__InterfaceMethod)NULL,121,-1,-1,0,3,1,4,0,
"i - - 0 - mult U 'PartLund' - 0 - first "
"d - - 1 - wt",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("BEchain2",700,(G__InterfaceMethod)NULL,121,-1,-1,0,4,1,4,0,
"u 'KorEvent' - 1 - event i - - 0 - kfP "
"i - - 1 - nP d - - 1 - wt",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("MakeWeight",998,G__BEwtMaker_MakeWeight_9_0,121,-1,-1,0,4,1,1,0,
"u 'KorEvent' - 1 - event l - - 1 - ntot "
"d - - 1 - wtot d - - 1 - wtot2",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("BookLSP",634,G__BEwtMaker_BookLSP_0_1,121,-1,-1,0,5,1,1,0,
"u 'KorEvent' - 1 - event i - - 0 - kfP "
"i - - 0 - ntot d - - 0 - wt "
"d - - 0 - wt2",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Q2pair",559,G__BEwtMaker_Q2pair_1_1,100,-1,-1,0,2,1,1,0,
"U 'PartLund' - 0 - first U 'PartLund' - 0 - second",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__BEwtMaker_DeclFileName_2_1,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__BEwtMaker_DeclFileLine_3_1,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__BEwtMaker_ImplFileName_4_1,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__BEwtMaker_ImplFileLine_5_1,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__BEwtMaker_Class_Version_6_1,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__BEwtMaker_Class_7_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__BEwtMaker_Dictionary_8_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__BEwtMaker_IsA_9_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__BEwtMaker_ShowMembers_0_2,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__BEwtMaker_Streamer_1_2,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic destructor
   G__memfunc_setup("~BEwtMaker",992,G__BEwtMaker_wABEwtMaker_2_2,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncJetAnalyzer(void) {
   /* JetAnalyzer */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer));
   G__memfunc_setup("JetAnalyzer",1129,G__JetAnalyzer_JetAnalyzer_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Book",395,G__JetAnalyzer_Book_2_0,121,-1,-1,0,4,1,1,0,
"u 'KorEvent' - 1 - event l - - 0 - ntot "
"d - - 0 - wt d - - 0 - wt2",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__JetAnalyzer_DeclFileName_3_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__JetAnalyzer_DeclFileLine_4_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__JetAnalyzer_ImplFileName_5_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__JetAnalyzer_ImplFileLine_6_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__JetAnalyzer_Class_Version_7_0,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__JetAnalyzer_Class_8_0,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__JetAnalyzer_Dictionary_9_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__JetAnalyzer_IsA_0_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__JetAnalyzer_ShowMembers_1_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__JetAnalyzer_Streamer_2_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic copy constructor
   G__memfunc_setup("JetAnalyzer",1129,G__JetAnalyzer_JetAnalyzer_3_1,(int)('i'),G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),-1,0,1,1,1,0,"u 'JetAnalyzer' - 1 - -",(char*)NULL,(void*)NULL,0);
   // automatic destructor
   G__memfunc_setup("~JetAnalyzer",1255,G__JetAnalyzer_wAJetAnalyzer_4_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncSemaph(void) {
   /* Semaph */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_Semaph));
   G__memfunc_setup("Semaph",606,G__Semaph_Semaph_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_Semaph),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Semaph",606,(G__InterfaceMethod)NULL,105,G__get_linked_tagnum(&G__KorwDictLN_Semaph),-1,0,1,1,4,0,"u 'Semaph' - 1 - org",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Initialize",1042,G__Semaph_Initialize_3_0,121,-1,-1,0,1,1,1,0,"u 'TString' - 1 - Semaphore",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ReadStatus",1024,G__Semaph_ReadStatus_4_0,121,-1,-1,0,1,1,1,0,"u 'TString' - 1 - Semaphore",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__Semaph_DeclFileName_5_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__Semaph_DeclFileLine_6_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__Semaph_ImplFileName_7_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__Semaph_ImplFileLine_8_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__Semaph_Class_Version_9_0,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__Semaph_Class_0_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__Semaph_Dictionary_1_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__Semaph_IsA_2_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__Semaph_ShowMembers_3_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__Semaph_Streamer_4_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic destructor
   G__memfunc_setup("~Semaph",732,G__Semaph_wASemaph_5_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncROBOL(void) {
   /* ROBOL */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__KorwDictLN_ROBOL));
   G__memfunc_setup("ROBOL",382,G__ROBOL_ROBOL_0_0,105,G__get_linked_tagnum(&G__KorwDictLN_ROBOL),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Initialize",1042,G__ROBOL_Initialize_2_0,121,-1,-1,0,1,1,1,0,"l - - 1 - NevTot",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Production",1063,G__ROBOL_Production_3_0,121,-1,-1,0,1,1,1,0,"l - - 1 - iEvent",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Finalize",818,G__ROBOL_Finalize_4_0,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileName",1145,G__ROBOL_DeclFileName_5_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("DeclFileLine",1152,G__ROBOL_DeclFileLine_6_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileName",1171,G__ROBOL_ImplFileName_7_0,67,-1,-1,0,0,1,1,1,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("ImplFileLine",1178,G__ROBOL_ImplFileLine_8_0,105,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class_Version",1339,G__ROBOL_Class_Version_9_0,115,-1,G__defined_typename("Version_t"),0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Class",502,G__ROBOL_Class_0_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("Dictionary",1046,G__ROBOL_Dictionary_1_1,121,-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__memfunc_setup("IsA",253,G__ROBOL_IsA_2_1,85,G__get_linked_tagnum(&G__KorwDictLN_TClass),-1,0,0,1,1,8,"",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("ShowMembers",1132,G__ROBOL_ShowMembers_3_1,121,-1,-1,0,2,1,1,0,
"u 'TMemberInspector' - 1 - insp C - - 0 - parent",(char*)NULL,(void*)NULL,1);
   G__memfunc_setup("Streamer",835,G__ROBOL_Streamer_4_1,121,-1,-1,0,1,1,1,0,"u 'TBuffer' - 1 - b",(char*)NULL,(void*)NULL,1);
   // automatic destructor
   G__memfunc_setup("~ROBOL",508,G__ROBOL_wAROBOL_5_1,(int)('y'),-1,-1,0,0,1,1,0,"",(char*)NULL,(void*)NULL,0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncKorwDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
extern "C" void G__cpp_setup_globalKorwDict() {

   /* Setting up global variables */
   G__resetplocal();


   G__resetglobalenv();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
extern "C" void G__cpp_setup_funcKorwDict() {
   G__lastifuncposition();


   G__resetifuncposition();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__KorwDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_ofstream = { "ofstream" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_VLorenz = { "VLorenz" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_PartLund = { "PartLund" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_KorEvent = { "KorEvent" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_ifstream = { "ifstream" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_KoralwMaker = { "KoralwMaker" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_TH1F = { "TH1F" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_TNtuple = { "TNtuple" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_BEwtMaker = { "BEwtMaker" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_JetAnalyzer = { "JetAnalyzer" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_Semaph = { "Semaph" , 99 , -1 };
G__linked_taginfo G__KorwDictLN_ROBOL = { "ROBOL" , 99 , -1 };

extern "C" void G__cpp_setup_tagtableKorwDict() {

   /* Setting up class,struct,union tag entry */
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_VLorenz),sizeof(VLorenz),-1,0,"VLorenz  class",G__setup_memvarVLorenz,G__setup_memfuncVLorenz);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_PartLund),sizeof(PartLund),-1,0,"PartLund  class, <--for dictionary",G__setup_memvarPartLund,G__setup_memfuncPartLund);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_KorEvent),sizeof(KorEvent),-1,0,"KorEvent  class, <--for dictionary",G__setup_memvarKorEvent,G__setup_memfuncKorEvent);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_KoralwMaker),sizeof(KoralwMaker),-1,0,"KorEvent  class, <--for dictionary",G__setup_memvarKoralwMaker,G__setup_memfuncKoralwMaker);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_BEwtMaker),sizeof(BEwtMaker),-1,0,(char*)NULL,G__setup_memvarBEwtMaker,G__setup_memfuncBEwtMaker);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_JetAnalyzer),sizeof(JetAnalyzer),-1,0,(char*)NULL,G__setup_memvarJetAnalyzer,G__setup_memfuncJetAnalyzer);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_Semaph),sizeof(Semaph),-1,0,(char*)NULL,G__setup_memvarSemaph,G__setup_memfuncSemaph);
   G__tagtable_setup(G__get_linked_tagnum(&G__KorwDictLN_ROBOL),sizeof(ROBOL),-1,0,(char*)NULL,G__setup_memvarROBOL,G__setup_memfuncROBOL);
}
extern "C" void G__cpp_setupKorwDict() {
  G__check_setup_version(51111,"G__cpp_setupKorwDict()");
  G__set_cpp_environmentKorwDict();
  G__cpp_setup_tagtableKorwDict();

  G__cpp_setup_inheritanceKorwDict();

  G__cpp_setup_typetableKorwDict();

  G__cpp_setup_memvarKorwDict();

  G__cpp_setup_memfuncKorwDict();
  G__cpp_setup_globalKorwDict();
  G__cpp_setup_funcKorwDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncKorwDict();
  return;
}
class G__cpp_setup_initKorwDict {
  public:
    G__cpp_setup_initKorwDict() { G__add_setup_func("KorwDict",&G__cpp_setupKorwDict); }
   ~G__cpp_setup_initKorwDict() { G__remove_setup_func("KorwDict"); }
};
G__cpp_setup_initKorwDict G__cpp_setup_initializerKorwDict;

//
// File generated by /home/jadach/lib/root/bin/rootcint at Mon Jul 27 12:28:59 1998.
// Do NOT change. Changes will be lost next time file is generated
//

#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, VLorenz *&obj)
{
   // Read a pointer to an object of class VLorenz.

   obj = (VLorenz *) buf.ReadObject(VLorenz::Class());
   return buf;
}

//______________________________________________________________________________
void VLorenz::Streamer(TBuffer &R__b)
{
   // Stream an object of class VLorenz.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.ReadStaticArray(m_comp);
   } else {
      R__b.WriteVersion(VLorenz::IsA());
      TObject::Streamer(R__b);
      R__b.WriteArray(m_comp, 4);
   }
}

//______________________________________________________________________________
void VLorenz::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class VLorenz.

   TClass *R__cl  = VLorenz::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "m_comp[4]", m_comp);
   TObject::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, PartLund *&obj)
{
   // Read a pointer to an object of class PartLund.

   obj = (PartLund *) buf.ReadObject(PartLund::Class());
   return buf;
}

//______________________________________________________________________________
void PartLund::Streamer(TBuffer &R__b)
{
   // Stream an object of class PartLund.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> previous;
      R__b >> next;
      R__b >> m_lserial;
      R__b >> m_kstatus;
      R__b >> m_kflavor;
      R__b >> m_kparent;
      R__b >> m_kFirstChild;
      R__b >> m_kLastChild;
      m_pmom.Streamer(R__b);
      R__b >> m_pmass;
      m_Vertex.Streamer(R__b);
      R__b >> m_LifeTime;
   } else {
      R__b.WriteVersion(PartLund::IsA());
      TObject::Streamer(R__b);
      R__b << previous;
      R__b << next;
      R__b << m_lserial;
      R__b << m_kstatus;
      R__b << m_kflavor;
      R__b << m_kparent;
      R__b << m_kFirstChild;
      R__b << m_kLastChild;
      m_pmom.Streamer(R__b);
      R__b << m_pmass;
      m_Vertex.Streamer(R__b);
      R__b << m_LifeTime;
   }
}

//______________________________________________________________________________
void PartLund::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class PartLund.

   TClass *R__cl  = PartLund::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "*previous", &previous);
   R__insp.Inspect(R__cl, R__parent, "*next", &next);
   R__insp.Inspect(R__cl, R__parent, "m_lserial", &m_lserial);
   R__insp.Inspect(R__cl, R__parent, "m_kstatus", &m_kstatus);
   R__insp.Inspect(R__cl, R__parent, "m_kflavor", &m_kflavor);
   R__insp.Inspect(R__cl, R__parent, "m_kparent", &m_kparent);
   R__insp.Inspect(R__cl, R__parent, "m_kFirstChild", &m_kFirstChild);
   R__insp.Inspect(R__cl, R__parent, "m_kLastChild", &m_kLastChild);
   m_pmom.ShowMembers(R__insp, strcat(R__parent,"m_pmom.")); R__parent[R__ncp] = 0;
   R__insp.Inspect(R__cl, R__parent, "m_pmass", &m_pmass);
   m_Vertex.ShowMembers(R__insp, strcat(R__parent,"m_Vertex.")); R__parent[R__ncp] = 0;
   R__insp.Inspect(R__cl, R__parent, "m_LifeTime", &m_LifeTime);
   TObject::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, KorEvent *&obj)
{
   // Read a pointer to an object of class KorEvent.

   obj = (KorEvent *) buf.ReadObject(KorEvent::Class());
   return buf;
}

//______________________________________________________________________________
void KorEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class KorEvent.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> m_nphot;
      int R__i;
      for (R__i = 0; R__i < 100; R__i++)
         m_photmom[R__i].Streamer(R__b);
      m_wminus.Streamer(R__b);
      m_wplus.Streamer(R__b);
      m_ferm1.Streamer(R__b);
      m_ferm2.Streamer(R__b);
      m_ferm3.Streamer(R__b);
      m_ferm4.Streamer(R__b);
      R__b >> m_npart;
      for (R__i = 0; R__i < 20000; R__i++)
         m_part[R__i].Streamer(R__b);
      R__b >> m_njet;
      for (R__i = 0; R__i < 200; R__i++)
         m_jet[R__i].Streamer(R__b);
   } else {
      R__b.WriteVersion(KorEvent::IsA());
      TObject::Streamer(R__b);
      R__b << m_nphot;
      int R__i;
      for (R__i = 0; R__i < 100; R__i++)
         m_photmom[R__i].Streamer(R__b);
      m_wminus.Streamer(R__b);
      m_wplus.Streamer(R__b);
      m_ferm1.Streamer(R__b);
      m_ferm2.Streamer(R__b);
      m_ferm3.Streamer(R__b);
      m_ferm4.Streamer(R__b);
      R__b << m_npart;
      for (R__i = 0; R__i < 20000; R__i++)
         m_part[R__i].Streamer(R__b);
      R__b << m_njet;
      for (R__i = 0; R__i < 200; R__i++)
         m_jet[R__i].Streamer(R__b);
   }
}

//______________________________________________________________________________
void KorEvent::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class KorEvent.

   TClass *R__cl  = KorEvent::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "m_nphot", &m_nphot);
   R__insp.Inspect(R__cl, R__parent, "m_photmom[100]", m_photmom);
   m_wminus.ShowMembers(R__insp, strcat(R__parent,"m_wminus.")); R__parent[R__ncp] = 0;
   m_wplus.ShowMembers(R__insp, strcat(R__parent,"m_wplus.")); R__parent[R__ncp] = 0;
   m_ferm1.ShowMembers(R__insp, strcat(R__parent,"m_ferm1.")); R__parent[R__ncp] = 0;
   m_ferm2.ShowMembers(R__insp, strcat(R__parent,"m_ferm2.")); R__parent[R__ncp] = 0;
   m_ferm3.ShowMembers(R__insp, strcat(R__parent,"m_ferm3.")); R__parent[R__ncp] = 0;
   m_ferm4.ShowMembers(R__insp, strcat(R__parent,"m_ferm4.")); R__parent[R__ncp] = 0;
   R__insp.Inspect(R__cl, R__parent, "m_npart", &m_npart);
   R__insp.Inspect(R__cl, R__parent, "m_part[20000]", m_part);
   R__insp.Inspect(R__cl, R__parent, "m_njet", &m_njet);
   R__insp.Inspect(R__cl, R__parent, "m_jet[200]", m_jet);
   TObject::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, KoralwMaker *&obj)
{
   // Read a pointer to an object of class KoralwMaker.

   obj = (KoralwMaker *) buf.ReadObject(KoralwMaker::Class());
   return buf;
}

//______________________________________________________________________________
void KoralwMaker::Streamer(TBuffer &R__b)
{
   // Stream an object of class KoralwMaker.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> NevTot;
      R__b >> m_EvenCounter;
      R__b >> m_PrintLimit;
      R__b.ReadStaticArray(xpar);
   } else {
      R__b.WriteVersion(KoralwMaker::IsA());
      TNamed::Streamer(R__b);
      R__b << NevTot;
      R__b << m_EvenCounter;
      R__b << m_PrintLimit;
      R__b.WriteArray(xpar, 10001);
   }
}

//______________________________________________________________________________
void KoralwMaker::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class KoralwMaker.

   TClass *R__cl  = KoralwMaker::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "NevTot", &NevTot);
   R__insp.Inspect(R__cl, R__parent, "m_EvenCounter", &m_EvenCounter);
   R__insp.Inspect(R__cl, R__parent, "m_PrintLimit", &m_PrintLimit);
   R__insp.Inspect(R__cl, R__parent, "xpar[10001]", xpar);
   TNamed::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, BEwtMaker *&obj)
{
   // Read a pointer to an object of class BEwtMaker.

   obj = (BEwtMaker *) buf.ReadObject(BEwtMaker::Class());
   return buf;
}

//______________________________________________________________________________
void BEwtMaker::Streamer(TBuffer &R__b)
{
   // Stream an object of class BEwtMaker.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> BEievent;
      R__b >> BEprintlevel;
      R__b >> BEprintlast;
      R__b >> BErange;
      R__b >> BEFuncType;
      R__b >> BEpp;
      R__b >> BEradius;
      R__b >> BElambda;
      R__b >> BEavewt;
      R__b >> BElambda2;
      R__b >> BEavewt2;
      R__b >> hst_nclu;
      R__b >> hst_unQ2;
      R__b >> hst_wtQ2;
      R__b >> hst_wt2Q2;
      R__b >> hst_unQ3;
      R__b >> hst_wtQ3;
      R__b >> hst_wt2Q3;
      R__b >> ctuple2_counter;
      R__b >> ctuple2_max;
      R__b >> ctuple2;
      R__b >> ctuple3_counter;
      R__b >> ctuple3_max;
      R__b >> ctuple3;
   } else {
      R__b.WriteVersion(BEwtMaker::IsA());
      TNamed::Streamer(R__b);
      R__b << BEievent;
      R__b << BEprintlevel;
      R__b << BEprintlast;
      R__b << BErange;
      R__b << BEFuncType;
      R__b << BEpp;
      R__b << BEradius;
      R__b << BElambda;
      R__b << BEavewt;
      R__b << BElambda2;
      R__b << BEavewt2;
      R__b << (TObject*)hst_nclu;
      R__b << (TObject*)hst_unQ2;
      R__b << (TObject*)hst_wtQ2;
      R__b << (TObject*)hst_wt2Q2;
      R__b << (TObject*)hst_unQ3;
      R__b << (TObject*)hst_wtQ3;
      R__b << (TObject*)hst_wt2Q3;
      R__b << ctuple2_counter;
      R__b << ctuple2_max;
      R__b << ctuple2;
      R__b << ctuple3_counter;
      R__b << ctuple3_max;
      R__b << ctuple3;
   }
}

//______________________________________________________________________________
void BEwtMaker::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class BEwtMaker.

   TClass *R__cl  = BEwtMaker::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "BEievent", &BEievent);
   R__insp.Inspect(R__cl, R__parent, "BEprintlevel", &BEprintlevel);
   R__insp.Inspect(R__cl, R__parent, "BEprintlast", &BEprintlast);
   R__insp.Inspect(R__cl, R__parent, "BErange", &BErange);
   R__insp.Inspect(R__cl, R__parent, "BEFuncType", &BEFuncType);
   R__insp.Inspect(R__cl, R__parent, "BEpp", &BEpp);
   R__insp.Inspect(R__cl, R__parent, "BEradius", &BEradius);
   R__insp.Inspect(R__cl, R__parent, "BElambda", &BElambda);
   R__insp.Inspect(R__cl, R__parent, "BEavewt", &BEavewt);
   R__insp.Inspect(R__cl, R__parent, "BElambda2", &BElambda2);
   R__insp.Inspect(R__cl, R__parent, "BEavewt2", &BEavewt2);
   R__insp.Inspect(R__cl, R__parent, "*hst_nclu", &hst_nclu);
   R__insp.Inspect(R__cl, R__parent, "*hst_unQ2", &hst_unQ2);
   R__insp.Inspect(R__cl, R__parent, "*hst_wtQ2", &hst_wtQ2);
   R__insp.Inspect(R__cl, R__parent, "*hst_wt2Q2", &hst_wt2Q2);
   R__insp.Inspect(R__cl, R__parent, "*hst_unQ3", &hst_unQ3);
   R__insp.Inspect(R__cl, R__parent, "*hst_wtQ3", &hst_wtQ3);
   R__insp.Inspect(R__cl, R__parent, "*hst_wt2Q3", &hst_wt2Q3);
   R__insp.Inspect(R__cl, R__parent, "ctuple2_counter", &ctuple2_counter);
   R__insp.Inspect(R__cl, R__parent, "ctuple2_max", &ctuple2_max);
   R__insp.Inspect(R__cl, R__parent, "*ctuple2", &ctuple2);
   R__insp.Inspect(R__cl, R__parent, "ctuple3_counter", &ctuple3_counter);
   R__insp.Inspect(R__cl, R__parent, "ctuple3_max", &ctuple3_max);
   R__insp.Inspect(R__cl, R__parent, "*ctuple3", &ctuple3);
   TNamed::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, JetAnalyzer *&obj)
{
   // Read a pointer to an object of class JetAnalyzer.

   obj = (JetAnalyzer *) buf.ReadObject(JetAnalyzer::Class());
   return buf;
}

//______________________________________________________________________________
void JetAnalyzer::Streamer(TBuffer &R__b)
{
   // Stream an object of class JetAnalyzer.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> hst_pene;
      R__b >> hst_jene;
      R__b >> jtuple;
   } else {
      R__b.WriteVersion(JetAnalyzer::IsA());
      TNamed::Streamer(R__b);
      R__b << (TObject*)hst_pene;
      R__b << (TObject*)hst_jene;
      R__b << jtuple;
   }
}

//______________________________________________________________________________
void JetAnalyzer::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class JetAnalyzer.

   TClass *R__cl  = JetAnalyzer::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "*hst_pene", &hst_pene);
   R__insp.Inspect(R__cl, R__parent, "*hst_jene", &hst_jene);
   R__insp.Inspect(R__cl, R__parent, "*jtuple", &jtuple);
   TNamed::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, Semaph *&obj)
{
   // Read a pointer to an object of class Semaph.

   obj = (Semaph *) buf.ReadObject(Semaph::Class());
   return buf;
}

//______________________________________________________________________________
void Semaph::Streamer(TBuffer &R__b)
{
   // Stream an object of class Semaph.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> ijklin;
      R__b >> ntotin;
      R__b >> ntot2n;
   } else {
      R__b.WriteVersion(Semaph::IsA());
      TNamed::Streamer(R__b);
      R__b << ijklin;
      R__b << ntotin;
      R__b << ntot2n;
   }
}

//______________________________________________________________________________
void Semaph::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class Semaph.

   TClass *R__cl  = Semaph::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   R__insp.Inspect(R__cl, R__parent, "ijklin", &ijklin);
   R__insp.Inspect(R__cl, R__parent, "ntotin", &ntotin);
   R__insp.Inspect(R__cl, R__parent, "ntot2n", &ntot2n);
   TNamed::ShowMembers(R__insp, R__parent);
}

//______________________________________________________________________________
TBuffer &operator>>(TBuffer &buf, ROBOL *&obj)
{
   // Read a pointer to an object of class ROBOL.

   obj = (ROBOL *) buf.ReadObject(ROBOL::Class());
   return buf;
}

//______________________________________________________________________________
void ROBOL::Streamer(TBuffer &R__b)
{
   // Stream an object of class ROBOL.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      KoralW.Streamer(R__b);
      BEfactory.Streamer(R__b);
      J4factory.Streamer(R__b);
      current_event.Streamer(R__b);
      R__b >> hst_BEwt;
      R__b >> hst_BEwtAR;
      R__b >> KeyRej;
      R__b >> WtMax;
   } else {
      R__b.WriteVersion(ROBOL::IsA());
      TNamed::Streamer(R__b);
      KoralW.Streamer(R__b);
      BEfactory.Streamer(R__b);
      J4factory.Streamer(R__b);
      current_event.Streamer(R__b);
      R__b << (TObject*)hst_BEwt;
      R__b << (TObject*)hst_BEwtAR;
      R__b << KeyRej;
      R__b << WtMax;
   }
}

//______________________________________________________________________________
void ROBOL::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
   // Inspect the data members of an object of class ROBOL.

   TClass *R__cl  = ROBOL::IsA();
   Int_t   R__ncp = strlen(R__parent);
   if (R__ncp || R__cl || R__insp.IsA()) { }
   KoralW.ShowMembers(R__insp, strcat(R__parent,"KoralW.")); R__parent[R__ncp] = 0;
   BEfactory.ShowMembers(R__insp, strcat(R__parent,"BEfactory.")); R__parent[R__ncp] = 0;
   J4factory.ShowMembers(R__insp, strcat(R__parent,"J4factory.")); R__parent[R__ncp] = 0;
   current_event.ShowMembers(R__insp, strcat(R__parent,"current_event.")); R__parent[R__ncp] = 0;
   R__insp.Inspect(R__cl, R__parent, "*hst_BEwt", &hst_BEwt);
   R__insp.Inspect(R__cl, R__parent, "*hst_BEwtAR", &hst_BEwtAR);
   R__insp.Inspect(R__cl, R__parent, "KeyRej", &KeyRej);
   R__insp.Inspect(R__cl, R__parent, "WtMax", &WtMax);
   TNamed::ShowMembers(R__insp, R__parent);
}

