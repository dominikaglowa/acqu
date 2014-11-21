//
// Offline analysis of gamma-deuterium interaction
//
// TA2GammaDeuterium
//
//

#ifndef __TA2GammaDeuterium_h__
#define __TA2GammaDeuterium_h__

#include "TAcquRoot.h"
#include "TAcquFile.h"
#include "TA2Physics.h"
#include "TA2Tagger.h"
#include "TA2Ladder.h"
#include "TA2Calorimeter.h"
#include "TA2CentralApparatus.h"
#include "TA2CalArray.h"
#include "TA2CylMwpc.h"
#include "TA2PlasticPID.h"
#include "TA2Particle.h"
#include "TA2Taps.h"

#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TChain.h"


// Default tagger offset

enum { ELadderOffset=0 };
static const Map_t kLadderKeys[]= {
  {"Offset:",        ELadderOffset},
  {NULL,             -1}
};


class TA2GammaDeuterium : public TA2Physics {
 protected:

   TA2Tagger* fTAGG;                    // Glasgow photon tagger
   TA2Ladder* fLADD;                    // Ladder
   TA2CentralApparatus* fCA;            // Central apparatus
   TA2CalArray* fCB;			// Crystal ball
   TA2CylMwpc* fMWPC;                   // Multi-wire chamber
   TA2PlasticPID* fPID;			// PID detector
   TA2Taps* fTAPS;               	// TAPS calorumeter

   TFile *fGammaDFile;                  // ROOT file to store data
   TTree *fGammaDRun;                   // ROOT tree to store run info
   TTree *fGammaDTree;                  // ROOT tree to store data

   UInt_t fRun;                         // run number

// Variables to store tagger info

   UInt_t fNTag;                        // number of tagger hits
   UInt_t fNMultTag;			// tagger hits multiplicity
   Int_t* fTagCh;                       // tagger hits
   Double_t* fTagTime;                  // tagger time
   Double_t* fTagEnergy;                // incident photon energy
   

// Variables to store counters info

   Double_t* fScalCurr;                 // Current scaler values
   Double_t* fScalAcc;                  // Accumulated scalers
   Double_t* fScalAccCorr;              // Corrected accumulated scalers
   UInt_t* fScalerIndex;                // Tagger indexes
   Int_t fScalerOffset;                 // Tagger indexes offset

// Varaibles to store CB / TAPS info

   UInt_t fMaxTAGG, fMaxCB, fMaxTAPS, fMaxPID; // max. size of arrays
   UInt_t fNCBcl, fNTAPScl;

   Double_t fBeamEnergy;		// beam energy
   Double_t* fCBEnergy;                 // photons enegry
   Double_t* fCBTheta;                  // theta angle / deg.
   Double_t* fCBPhi;                    // phi angle / deg.
   Double_t* fCBTime;                   // photons time
   Double_t* fCBRadius;			// cluster radius
   Double_t* fTAPSEnergy;               // protons energy
   Double_t* fTAPSTheta;                // theta angle / deg.
   Double_t* fTAPSPhi;                  // phi angle / deg.
   Double_t* fTAPSTime;                 // time
   Double_t* fTAPSRadius;		// cluster radius
   Int_t*    fTAPSVeto;			// Veto for charged particles

// Variables to store pid info

   UInt_t fNPid;                        // number of pid hits
   Int_t* fPidCh;                       // pid hits
   Double_t* fPidEnergy;                // pid energy
   Double_t* fPidPhi;                   // pid angle
   Double_t* fPidTime;                  // pid time

// Variables to store MWPC info

   UInt_t fMaxInter; 			// max. size of arrays

   UInt_t fNinter1;			// Number of hits for 1st chamber 
   Double_t* fInter1X;			// X position of hits ...
   Double_t* fInter1Y;			// Y position of hits ...
   Double_t* fInter1Z;   		// Z position of hits ...
   Double_t* fInter1Phi;		// Phi - angle of hits ...
   Int_t* fInter1Type;    		// hits type ...

   UInt_t fNinter2;			// Number of hits for 2st chamber 
   Double_t* fInter2X;			// X position of hits ...
   Double_t* fInter2Y;			// Y position of hits ...
   Double_t* fInter2Z;   		// Z position of hits ...
   Double_t* fInter2Phi;		// Phi - angle of hits ...
   Int_t* fInter2Type;    		// hits type ...

// MC variables 

   UInt_t fMaxPart;			// Max number of particles
   Int_t fNpart;			// Number of particles
   Double_t* fVertex;			// Initial vertex position
   Double_t* fEbeam;			// Px, Py, Pz of beam gamma 
   Int_t* fIdPart;			// Particle Geant4 code
   Double_t* fPartT;			// Tkin of particles
   Double_t* fPartX;			// Px of secondary partcile
   Double_t* fPartY;			// Px of secondary partcile
   Double_t* fPartZ;			// Px of secondary partcile

 public:

   TA2GammaDeuterium( const char*, TA2Analysis* );

   virtual ~TA2GammaDeuterium();
   virtual void LoadVariable();                   // load variables to plot/cut
   virtual void PostInit( );                      // detector initialization
   virtual void Reconstruct();                    // ask each event

   virtual void SetConfig(Char_t* , Int_t);       // set config
   virtual void CreateROOTFile();                 // create ROOT file
   virtual Char_t* CreateROOTName(const char*);   // create ROOT filename
   virtual void ProcessEnd();                     // finish part
   virtual TA2DataManager* CreateChild( const char*, Int_t ){ return NULL;}
   virtual void DecodeLadderINFO();               // get counters parameters
   virtual void ParseMisc(char* line);            // read original dat file
   virtual void GetMCINFO();			  // read initail MC INFO
   virtual void GetTaggerINFO();                  // get tagger parameters
   virtual void GetCrystalBallINFO();             // get CB parameters
   virtual void GetPidINFO();                     // get PID parameters
   virtual void GetMWPCINFO();                    // get MWPC parameters
   virtual void GetTAPSINFO();      		  // get TAPS parameters 

  ClassDef(TA2GammaDeuterium,1)
};

#endif

