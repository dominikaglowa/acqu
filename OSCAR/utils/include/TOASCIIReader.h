// SVN Info: $Id: TOASCIIReader.h 1257 2012-07-26 15:33:13Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2007-2008
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOASCIIReader                                                        //
//                                                                      //
// This class provides the basic functionality for ASCII file reading.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef OSCAR_TOASCIIReader
#define OSCAR_TOASCIIReader

#include <fstream>
#include <cstring>

#include "TObject.h"
//#include "TError.h"


using namespace std;


class TOASCIIReader : public TObject
{

protected:
    ifstream* fFile;                                    // ASCII file to read          
    Char_t fFileName[256];                              // Name of the file 

    Int_t ReadLines(Bool_t inCallParser);
    virtual void ParseLine(const char* inLine) = 0;     // line parser method

public:
    TOASCIIReader() : fFile(0) { }
    TOASCIIReader(const char* inFileName);
    virtual ~TOASCIIReader();
    
    void ReadFile();
    virtual Int_t CountLines();                         // line counter method

    ClassDef(TOASCIIReader, 0)  // abstract ASCII file reader
};

#endif

