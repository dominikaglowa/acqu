// SVN Info: $Id: TCCalibType.cxx 913 2011-05-18 22:15:42Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibType                                                          //
//                                                                      //
// Calibration type class.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TCCalibType.h"

ClassImp(TCCalibType)


//______________________________________________________________________________
TCCalibType::TCCalibType(const Char_t* name, const Char_t* title)
    : TNamed(name, title)
{
    // Constructor.
    
    fData = new TList();
}

//______________________________________________________________________________
TCCalibType::~TCCalibType()
{
    // Destructor.

    if (fData) delete fData;
}

//______________________________________________________________________________
void TCCalibType::Print()
{
    // Print the content of this class.
    
    printf("CaLib Type Information\n");
    printf("Name           : %s\n", GetName());
    printf("Title          : %s\n", GetTitle());
    printf("CaLib data     : ");
    TIter next(fData);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        printf("%s,", d->GetName());
    }
    printf("\r\n");
}

