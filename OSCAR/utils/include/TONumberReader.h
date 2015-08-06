// SVN Info: $Id: TONumberReader.h 925 2011-05-29 03:25:14Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2007-2008
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TONumberReader                                                       //
//                                                                      //
// This class reads numbers arranged in multiple columns row by row     //
// and stores them in a list of arrays.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef OSCAR_TONumberReader
#define OSCAR_TONumberReader

#include <sstream>

#include "TMath.h"

#include "TOASCIIReader.h"
#include "TOArray2D.h"
#include "TOSUtils.h"


class TONumberReader : public TOASCIIReader
{

private:
    Int_t fNcol;                                   // Number of columns
    Int_t fNrow;                                   // Number of rows
    TOArray2D* fData;                              // 2-dim. array to store the numbers in

    virtual void ParseLine(const char* inLine);

public:
    TONumberReader() : TOASCIIReader(), fNcol(0), fNrow(0), fData(0) { }
    TONumberReader(const char* inFileName, Int_t inColumns); 
    virtual ~TONumberReader();
    
    inline TOArray2D* GetData() const { return fData; }
    Int_t GetNcol() const { return fNcol; }
    Int_t GetNrow() const { return fNrow; }
    Double_t* GetColumn(Int_t c) const { return fData->GetColumn(c); }
    Double_t* GetSum(Int_t c1, Int_t c2, Double_t w1 = 1.0, Double_t w2 = 1.0) const;
    Double_t* GetSum(Int_t c1, TONumberReader* r, Int_t c2, Double_t w1 = 1.0, Double_t w2 = 1.0) const;
    Double_t GetColumnAverage(Int_t c);
    Double_t GetColumnAverageWeighted(Int_t c, Double_t* outErr = 0);
    Double_t GetRowAverage(Int_t r);
    Double_t GetRowAverageWeighted(Int_t r, Double_t* outErr = 0);
    virtual Int_t CountLines();

    ClassDef(TONumberReader, 0)  // ASCII file number reader
};

#endif

