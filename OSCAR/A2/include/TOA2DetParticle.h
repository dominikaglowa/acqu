// SVN Info: $Id: TOA2DetParticle.h 1345 2012-09-12 14:29:45Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2008
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOA2DetParticle                                                      //
//                                                                      //
// Class for reconstructed particles in the A2 setup.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef OSCAR_TOA2DetParticle
#define OSCAR_TOA2DetParticle

#include "TOA2BaseDetParticle.h"
#include "TOSUtils.h"


class TOA2DetParticle : public TOA2BaseDetParticle
{

protected:
    UInt_t* fClusterHits;                       //[fClusterSize] Cluster hit elements 
    Double_t* fClusterHitEnergies;              //[fClusterSize] Cluster hit energies
    Double_t* fClusterHitTimes;                 //[fClusterSize] Cluster hit times
    Double_t fErrorTheta;                       // error in theta angle
    Double_t fErrorPhi;                         // error in phi angle
    Double_t fErrorEnergy;                      // error in energy

public:
    TOA2DetParticle() : TOA2BaseDetParticle(),
                        fClusterHits(0), fClusterHitEnergies(0), fClusterHitTimes(0),
                        fErrorTheta(0), fErrorPhi(0), fErrorEnergy(0) { }

    TOA2DetParticle(const TOA2DetParticle& orig);
    virtual ~TOA2DetParticle();
 
    Double_t CalcChi2IMError(TOA2DetParticle* p);

    virtual void Print(Option_t* option = "") const;

    void SetClusterHits(Int_t size, UInt_t* hits);
    void SetClusterHitEnergies(Int_t size, Double_t* energies);
    void SetClusterHitTimes(Int_t size, Double_t* times);
    void SetErrorTheta(Double_t e) { fErrorTheta = e; }
    void SetErrorPhi(Double_t e) { fErrorPhi = e; }
    void SetErrorEnergy(Double_t e) { fErrorEnergy = e; }

    UInt_t* GetClusterHits() const { return fClusterHits; }
    Double_t* GetClusterHitEnergies() const { return fClusterHitEnergies; }
    Double_t* GetClusterHitTimes() const { return fClusterHitTimes; }
    Double_t GetErrorTheta() const { return fErrorTheta; }
    Double_t GetErrorPhi() const { return fErrorPhi; }
    Double_t GetErrorEnergy() const { return fErrorEnergy; }

    ClassDef(TOA2DetParticle, 1)  // Detected particles in the A2 setup
};

#endif

