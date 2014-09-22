// SVN Info: $Id: TCCalibTargetPosition.h 768 2011-01-26 16:57:17Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibTargetPosition                                                //
//                                                                      //
// Calibration module for the target position.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCCALIBTARGETPOSITION_H
#define TCCALIBTARGETPOSITION_H

#include "TCanvas.h"
#include "TH2.h"
#include "TLine.h"

#include "TCCalib.h"
#include "TCFileManager.h"


class TCCalibTargetPosition : public TCCalib
{

private:
    Double_t fSigmaPrev;                // previous pi0 peak sigma
    TLine* fLine;                       // indicator line
    
    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

public:
    TCCalibTargetPosition();
    virtual ~TCCalibTargetPosition();

    virtual void Write();
    virtual void PrintValues();

    ClassDef(TCCalibTargetPosition, 0) // target position calibration
};

#endif

