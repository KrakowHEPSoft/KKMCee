//* BoX-FORMATs for nice and flexible outputs


#define BXOPE(file) file<<\
"********************************************************************************"<<endl<<\
"*                                                                              *"<<endl
#define BXTXT(file,text) file<<\
"*                   "<<setw(40)<<         text           <<"                   *"<<endl
#define BX1I(file,name,numb,text) file<<\
"* "<<setw(10)<<name<<" = "<<setw(10)<<numb<<" = "          <<setw(50)<<text<<" *"<<endl
#define BX1F(file,name,numb,text) file<<"* "<<setw(10)<<name<<\
" = "<<setw(15)<<setprecision(8)<<numb<<"   =    "<<setw(40)<<text<<" *"<<endl
#define BX2F(file,name,numb,err,text) file<<"* "<<setw(10)<<name<<\
" = "<<setw(15)<<setprecision(8)<<numb<<" +- "<<setw(15)<<setprecision(8)<<err<<\
                                                      "  = "<<setw(25)<<text<<" *"<<endl
#define BXCLO(file) file<<\
"*                                                                              *"<<endl<<\
"********************************************************************************"<<endl

