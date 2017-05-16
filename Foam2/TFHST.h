///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Primitive Class of 1-dimensional histograms                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef TFHST_1
#define TFHST_1
class TFHST{                            // histogram
  private:
    double m_Entries;                  // no of Entries
    int    m_Nb;                       // no of Bins
    double m_xmin;                     // (xmin,xmax) is histo range
    double m_xmax;                     // (xmin,xmax) is histo range
    double *m_Bin1;                    // Bin sum wt
    double *m_Bin2;                    // Bin sum wt*wt
    const static int m_chat =1;        // chatability (for debug)
  public:
    TFHST(void);                       // Constructor for empty histo
    TFHST(double, double, int );       // Constructor PRINCIPAL
    TFHST(TFHST &);                    // Copy constructor ***
    ~TFHST();                          // destructor
    //
    TFHST& operator =(TFHST&);         // Substitution operator 
    double& operator[](int);           // Returns bin content
    double& operator()(int);           // Returns bin error content
    //
    void Fill(double, double);         // Fill histogram
    double GetBinContent(int );        // Get content of bin nb
    double GetBinError(  int );        // Get error of bin nb
    double GetEntries(){ return m_Entries; }; // Get no. of entries
    int    GetNbin(){    return m_Nb; };      // Get no. of bin nb
    void Print(void);                  // Prints histogram
    void Print(char*);
    void Reset();                      // Reset bin content
};
#endif
