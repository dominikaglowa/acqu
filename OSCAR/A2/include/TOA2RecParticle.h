// SVN Info: $Id: TOA2RecParticle.h 1257 2012-07-26 15:33:13Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2008-2010
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOA2RecParticle                                                      //
//                                                                      //
// Base class for reconstructed particles in the A2 setup.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef OSCAR_TOA2RecParticle
#define OSCAR_TOA2RecParticle

#include "TLorentzVector.h"
#include "TClass.h"
#include "TH1.h"

#include "TOGlobals.h"
#include "TOA2DetParticle.h"


class TOA2RecParticle : public TObject
{

protected:
    TLorentzVector* fP4;                            //-> 4-vector
    Int_t fNDetectedProducts;                       // number of detected decay products
    TOA2DetParticle** fDetectedProducts;            //[fNDetectedProducts] pointer to detected decay product list
    Double_t fChiSquare;                            // chi square value of reconstruction

protected:
    virtual void Calculate4Vector();

public:
    TOA2RecParticle(); 
    TOA2RecParticle(Int_t nParticle);
    TOA2RecParticle(const TOA2RecParticle& orig);
    virtual ~TOA2RecParticle();
 
    virtual Bool_t Reconstruct(Int_t nParticle, TOA2DetParticle** particleList);
    Double_t GetAverageTime(A2Detector_t detector) const;
    Int_t GetNDetectorHits(A2Detector_t detector) const;
    
    Bool_t CheckCoincidence(Double_t limitLowCB, Double_t limitHighCB, TH1* hCB,
                            Double_t limitLowTAPS, Double_t limitHighTAPS, TH1* hTAPS,
                            Double_t limitLowCBTAPS, Double_t limitHighCBTAPS, TH1* hCBTAPS);
    Bool_t CheckDetectorCoincidence(A2Detector_t detector, Double_t limitLow, Double_t limitHigh, TH1* h = 0);
    
    void CalculateConstrained4Vector(Double_t mass, TLorentzVector* outP4);

    Double_t GetTime() const;
    TLorentzVector* Get4Vector() const { return fP4; }
    Int_t GetNDetectedProducts() const { return fNDetectedProducts; }
    TOA2DetParticle** GetDetectedProducts() const { return fDetectedProducts; }
    TOA2DetParticle* GetDetectedProduct(Int_t index) const;
    Double_t GetChiSquare() const { return fChiSquare; }

    void Set4Vector(TLorentzVector* inP4) { if (fP4) *fP4 = *inP4; }
    void SetNDetectedProducts(Int_t inN) { fNDetectedProducts = inN; }
    void SetDetectedProducts(TOA2DetParticle** inProd) { fDetectedProducts = inProd; }

    void Clear(Option_t* option = "");
    virtual void Print(Option_t* option = "") const;
    
    ClassDef(TOA2RecParticle, 1)  // Reconstructed particles in the A2 setup
};


TOA2RecParticle operator+(const TOA2RecParticle& part1, const TOA2RecParticle& part2);

#endif

