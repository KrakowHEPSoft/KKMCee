// 27.03.2014 JZ

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
 {
  fstream in;
  fstream out;

//commented lines alows, to make another file wih list og generated channel only
// if zero appears as number of efents in channel, that mean it's not initialized properly
/*  
  fstream chan;
  chan.open("channels.txt", ios::out);
 */
  in.open("../tauola.output", ios::in);
  out.open("out.txt",	 ios::out);
  if (in.good() == true) cout<<"plik otwarty"<<endl;
  else cout<<"pliku nie ma, lub niedostępny do odczytu"<<endl;
  if (out.good() == false) cout<<"nie mozna utworzyc pliku do zapisu"<<endl;
  string dane;
  int event_count=0;
  while (in.fail() == false)
  {
   getline(in, dane);
   if (dane == "                            Event listing (standard)")
      {
        for (int i=1; i<=8; i++) getline(in, dane);
//        for (int i=1; i<=18; i++) {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}
        while (dane != "") {out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);}
        out<<"\n";
        event_count++;
        
        for (int i=1; i<6; i++)
        {
         getline(in, dane);
         if (dane == " ***************************************************************************")
            {
             getline(in, dane);
             out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);
             out.write(& dane[0], dane.length()); out<<"\n"; getline(in, dane);
             out<<"\n";
            }
        }
      }
    if (dane == " *NCHAN    NOEVTS  PART.WIDTH     ERROR       ROUTINE    DECAY MODE        *")
     {
      for(int i=1; i<(event_count + 1) ;i++)
      {
       getline(in, dane);
        if(i==event_count)
          {
           out.write(& dane[0], dane.length());
//           chan.write(& dane[0], dane.length());
//           chan<<"\n";
           out<<"\n \n";
           out<<"________________________________________next channel____________________________________________________________";
	   out<<"\n";
	  }
      }
     }
  }
  in.close();
  out.close();
  return 0;
 }
