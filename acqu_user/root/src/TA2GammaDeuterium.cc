//
// Offline analysis of gamma-deuterium interaction
//
// TA2GammaDeuterium
//
// Author : M. Martemyanov / April 2013
//

#include "TA2GammaDeuterium.h"

ClassImp(TA2GammaDeuterium)

//-----------------------------------------------------------------------------

TA2GammaDeuterium::TA2GammaDeuterium( const char* name, TA2Analysis* analysis )
  :TA2Physics( name, analysis ) {

  fTAGG = NULL; fLADD = NULL; fCA = NULL; fCB = NULL; fMWPC = NULL;

  fRun = 0; fNTag = 0; fNCBcl = 0; fNTAPScl = 0; fBeamEnergy = 0;

  fTagCh = NULL; fTagEnergy = NULL; fTagTime = NULL;

  fNpart = 0; fVertex = NULL; fEbeam = NULL; fIdPart = NULL; fPartT = NULL; 
  fPartX = NULL; fPartY = NULL; fPartZ = NULL;

  fCBEnergy = NULL; fCBTheta = NULL; fCBPhi = NULL; fCBTime = NULL; fCBRadius = NULL;

  fTAPSEnergy = NULL; fTAPSTheta = NULL; fTAPSPhi = NULL; fTAPSTime = NULL;
  fTAPSRadius = NULL; fTAPSVeto = NULL;

  fNPid = 0; fPidCh = NULL; fPidEnergy = NULL; fPidTime = NULL;
  fPidPhi = NULL;

  fScalCurr = NULL; fScalAcc = NULL; fScalAccCorr = NULL;
  fScalerIndex = NULL; fScalerOffset = 0;

  fNinter1 = 0; fInter1X = NULL; fInter1Y = NULL; fInter1Z = NULL; 
  fInter1Phi = NULL, fInter1Type = NULL;
  fNinter2 = 0; fInter2X = NULL; fInter2Y = NULL; fInter2Z = NULL; 
  fInter2Phi = NULL, fInter2Type = NULL;

}

//-----------------------------------------------------------------------------

TA2GammaDeuterium::~TA2GammaDeuterium() {

  delete fGammaDRun;
  delete fGammaDTree;
  delete fGammaDFile;

}

//---------------------------------------------------------------------------

void TA2GammaDeuterium::PostInit() {


//
// Create detector setup to be analysed
//

// Tagger part initialization

  fLADD = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild( "FPD");
  if ( !fLADD) PrintError( "", "<No Ladder class found>", EErrFatal);

  fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
  if (!fTAGG) PrintError( "", "<No Tagger class found>", EErrFatal);

  AddCmdList(kLadderKeys);
  fScalerIndex = fLADD->GetScalerIndex();

// Central appartus part (CB + MPWC + PID)

  fCA = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");
  if (!fCA) PrintError("","<No Central Apparatus class found>",EErrFatal);

  fCB = (TA2CalArray*)((TA2Analysis*)fParent)->GetGrandChild("NaI");  
  if (!fCB) PrintError("","<No CB class found>",EErrFatal); 

  fMWPC = (TA2CylMwpc*)((TA2Analysis*)fParent)->GetGrandChild("CylMWPC");
  if (!fMWPC) PrintError("","<No MWPC class found>",EErrFatal);

  fPID = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("PID");
  if (!fPID) PrintError("","<No PID class found>",EErrFatal);

// TAPS part initialization

  fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
  if (!fTAPS) PrintError( "", "<No TAPS class found>",EErrFatal);


//
// Open ROOT-file to store data
//


  CreateROOTFile();

  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------

void TA2GammaDeuterium::CreateROOTFile() {

//
// Create root file structure 
//

// Max. num. of particles in tagger


  if ( (fTAGG) && (fLADD) )  { 
    fNMultTag = fLADD->GetNMultihit();
    fMaxTAGG = fTAGG->GetMaxParticle();
    fBeamEnergy = fTAGG->GetBeamEnergy();
    if (fNMultTag > 0) fMaxTAGG = fMaxTAGG*fNMultTag;
  }

// For MC reject multi-tagger case (if accepted)

  if (! gAR->IsOnline() ) fNMultTag = -1;


// Max. num. of particles in CB

  if (fCB) fMaxCB = fCB->GetMaxCluster();
  else fMaxCB = 0;

// Max. num. of elements in PID
   if (fPID) fMaxPID = fPID->GetNelem();
   else fMaxPID = 0;

// Max. num. of intersections in MWPC
   if (fMWPC) fMaxInter = (UInt_t)fMWPC->GetMaxIntersect() + 1;
   else fMaxInter = 0; 

// Max. num. of particles in TAPS
   if (fTAPS) fMaxTAPS = fTAPS->GetMaxParticle();
   else fMaxTAPS = 0;


// MC parameters


  fMaxPart = 3;

  fVertex = new Double_t[3]; 
  fEbeam  = new Double_t[3]; 
  fIdPart = new Int_t[fMaxPart];
  fPartT  = new Double_t[fMaxPart];
  fPartX  = new Double_t[fMaxPart];
  fPartY  = new Double_t[fMaxPart];
  fPartZ  = new Double_t[fMaxPart];  
 
// Create arrays to be written

  fTagCh       = new Int_t[fMaxTAGG];

  fTagEnergy   = new Double_t[fMaxTAGG];
  fTagTime     = new Double_t[fMaxTAGG];

  fScalCurr    = new Double_t[fMaxTAGG];
  fScalAcc     = new Double_t[fMaxTAGG];
  fScalAccCorr = new Double_t[fMaxTAGG];

// CB cluster parameters

  fCBEnergy    = new Double_t[fMaxCB];
  fCBTheta     = new Double_t[fMaxCB];
  fCBPhi       = new Double_t[fMaxCB];
  fCBTime      = new Double_t[fMaxCB];
  fCBRadius    = new Double_t[fMaxCB];

// PID parameters

  fPidCh       = new Int_t[fMaxPID];
  fPidEnergy   = new Double_t[fMaxPID];
  fPidTime     = new Double_t[fMaxPID];
  fPidPhi      = new Double_t[fMaxPID];

// MWPC parameters

  fInter1X      = new Double_t[fMaxInter];
  fInter1Y      = new Double_t[fMaxInter];
  fInter1Z      = new Double_t[fMaxInter];
  fInter1Phi    = new Double_t[fMaxInter];
  fInter1Type   = new Int_t[fMaxInter];

  fInter2X      = new Double_t[fMaxInter];
  fInter2Y      = new Double_t[fMaxInter];
  fInter2Z      = new Double_t[fMaxInter];
  fInter2Phi    = new Double_t[fMaxInter];
  fInter2Type   = new Int_t[fMaxInter];

// TAPS cluster parameters

  fTAPSEnergy  = new Double_t[fMaxTAPS];
  fTAPSTheta   = new Double_t[fMaxTAPS];
  fTAPSPhi     = new Double_t[fMaxTAPS];
  fTAPSTime    = new Double_t[fMaxTAPS];
  fTAPSRadius  = new Double_t[fMaxTAPS];
  fTAPSVeto    = new Int_t[fMaxTAPS];

  if (gAR->IsOnline())
   fGammaDFile = new TFile(CreateROOTName("TA2GammaDeuterium"), "RECREATE", "", 3);
  else
   fGammaDFile = new TFile(CreateROOTName("TA2GammaDeuteriumMC"),"RECREATE", "", 3);


// Add additional histograms

// Add additional tree

   fGammaDRun = new TTree("TA2GammaDRun","");
   fGammaDRun->Branch("Run",            &fRun,          "Run/I");
   fGammaDRun->Branch("MaxTAGG",        &fMaxTAGG,      "MaxTAGG/I");
   fGammaDRun->Branch("MaxCB",          &fMaxCB,        "MaxCB/I");
   fGammaDRun->Branch("MaxInter",	&fMaxInter,	"MaxInter/I");
   fGammaDRun->Branch("MaxTAPS",        &fMaxTAPS,      "MaxTAPS/I");
   fGammaDRun->Branch("BeamEnergy",     &fBeamEnergy,   "BeamEnergy/D");

   fGammaDTree = new TTree("TA2GammaDTree","");

   if (! gAR->IsOnline()) {

     fGammaDTree->Branch("Npart",	&fNpart,	"Npart/I");
     fGammaDTree->Branch("Vertex",	fVertex,	"Vertex[3]/D");
     fGammaDTree->Branch("Ebeam",	fEbeam,		"Ebeam[3]/D");
     fGammaDTree->Branch("IdPart",	fIdPart,	"IdPart[Npart]/I");

     fGammaDTree->Branch("PartT",	fPartT,		"PartT[Npart]/D");	
     fGammaDTree->Branch("PartX",	fPartX,		"PartX[Npart]/D");	
     fGammaDTree->Branch("PartY",	fPartY,		"PartY[Npart]/D");
     fGammaDTree->Branch("PartZ",	fPartZ,		"PartZ[Npart]/D");  
 
   }

   fGammaDTree->Branch("NTag",          &fNTag,         "NTag/I");

   fGammaDTree->Branch("TagCh",         fTagCh,         "TagCh[NTag]/I");
   fGammaDTree->Branch("TagEnergy",     fTagEnergy,     "TagEnergy[NTag]/D");
   fGammaDTree->Branch("TagTime",       fTagTime,       "TagTime[NTag]/D");

   fGammaDTree->Branch("NCBcl",         &fNCBcl,        "NCBcl/I");
   fGammaDTree->Branch("CBEnergy",      fCBEnergy,      "CBEnergy[NCBcl]/D");
   fGammaDTree->Branch("CBTheta",       fCBTheta,       "CBTheta[NCBcl]/D");
   fGammaDTree->Branch("CBPhi",         fCBPhi,         "CBPhi[NCBcl]/D");
   fGammaDTree->Branch("CBTime",        fCBTime,        "CBTime[NCBcl]/D");
   fGammaDTree->Branch("CBRadius",      fCBRadius,      "CBRadius[NCBcl]/D");

   fGammaDTree->Branch("NPid",          &fNPid,         "NPid/I");
   fGammaDTree->Branch("PidCh",         fPidCh,         "PidCh[NPid]/I");
   fGammaDTree->Branch("PidEnergy",     fPidEnergy,     "PidEnergy[NPid]/D");
   fGammaDTree->Branch("PidPhi",        fPidPhi,        "PidPhi[NPid]/D");
   fGammaDTree->Branch("PidTime",       fPidTime,       "PidTime[NPid]/D");

   fGammaDTree->Branch("Ninter1",	&fNinter1,	"Ninter1/I");
   fGammaDTree->Branch("Inter1X",	fInter1X,	"Inter1X[Ninter1]/D");
   fGammaDTree->Branch("Inter1Y",	fInter1Y,	"Inter1Y[Ninter1]/D");
   fGammaDTree->Branch("Inter1Z",	fInter1Z,	"Inter1Z[Ninter1]/D");
   fGammaDTree->Branch("Inter1Phi",	fInter1Phi,	"Inter1Phi[Ninter1]/D");
   fGammaDTree->Branch("Inter1Type",	fInter1Type,	"Inter1Type[Ninter1]/I");

   fGammaDTree->Branch("Ninter2",	&fNinter2,	"Ninter2/I");
   fGammaDTree->Branch("Inter2X",	fInter2X,	"Inter2X[Ninter2]/D");
   fGammaDTree->Branch("Inter2Y",	fInter2Y,	"Inter2Y[Ninter2]/D");
   fGammaDTree->Branch("Inter2Z",	fInter2Z,	"Inter2Z[Ninter2]/D");
   fGammaDTree->Branch("Inter2Phi",	fInter2Phi,	"Inter2Phi[Ninter2]/D");
   fGammaDTree->Branch("Inter2Type",	fInter2Type,	"Inter2Type[Ninter2]/I");

   fGammaDTree->Branch("NTAPScl",       &fNTAPScl,      "NTAPScl/I");
   fGammaDTree->Branch("TAPSEnergy",    fTAPSEnergy,    "TAPSEnergy[NTAPScl]/D");
   fGammaDTree->Branch("TAPSTheta",     fTAPSTheta,     "TAPSTheta[NTAPScl]/D");
   fGammaDTree->Branch("TAPSPhi",       fTAPSPhi,       "TAPSPhi[NTAPScl]/D");
   fGammaDTree->Branch("TAPSTime",      fTAPSTime,      "TAPSTime[NTAPScl]/D");
   fGammaDTree->Branch("TAPSRadius",    fTAPSRadius,    "TAPSRadius[NTAPScl]/D");
   fGammaDTree->Branch("TAPSVeto",      fTAPSVeto,      "TAPSVeto[NTAPScl]/I");

  gROOT->cd();

}

//-----------------------------------------------------------------------------

Char_t* TA2GammaDeuterium::CreateROOTName(const char* spare) {

  enum {CSIZE = 200}; char buffer[CSIZE]; char strname[CSIZE];
  Char_t *CurrentFileName = new Char_t[CSIZE];

  char* split;

  if (gAR->IsOnline()) {
	split = strtok(strtok(gAR->GetFileName(),"."),"_");
   }
  else	{
        strcpy(strname, gAR->GetTreeFileList(0));
	split = strtok(strtok(strname,"."),"_");
  }

  while(split != NULL) {
   strcpy(buffer, split);
   split=strtok(NULL,"_");
  }

  fRun = atoi(buffer);
  sprintf(CurrentFileName,"%s_%04d.root",spare,fRun);
  
  return CurrentFileName;

}

//-----------------------------------------------------------------------------

void TA2GammaDeuterium::SetConfig(Char_t* line, Int_t key) {


 switch (key){
   case ELadderOffset:
    sscanf(line, "%d", &fScalerOffset);
   break;
   default:
    TA2Physics::SetConfig( line, key );
   break;
  }

}

// Parse line marked as Misc:

void TA2GammaDeuterium::ParseMisc(char* line) { }


//-----------------------------------------------------------------------------

inline void TA2GammaDeuterium::ProcessEnd() {

// Close ROOT file

// Add scaler standard information


   TList *listHist = (TList*)gROOT->GetList();

   fGammaDRun->Fill();
   fGammaDFile->cd();

   listHist->Write();
   listHist->Print();


   fGammaDRun->Print();
   fGammaDTree->Print();

//   fGammaDRun->Write();
//   fGammaDTree->Write();

    fGammaDFile->Write();

   char keyString[128];
   TIter next(fGammaDFile->GetListOfKeys());
   sprintf(keyString,"*;%d",((TKey*)next())->GetCycle()-1);
   cout << "TA2GammaDeuterium : keyString = " << keyString << endl;
   fGammaDFile->Delete(keyString);


   fGammaDFile->Close();


}

//_____________________________________________________________________________________


void TA2GammaDeuterium::LoadVariable(){

// Input name - variable pointer associations for any subsequent
// cut or histogram setup
// LoadVariable( "name", pointer-to-variable, type-spec );
// NB scaler variable pointers need the preceeding &
//    array variable pointers do not.
// type-spec ED prefix for a Double_t variable
//           EI prefix for an Int_t variable
// type-spec SingleX for a single-valued variable
//           MultiX  for a multi-valued variable

   TA2Physics::LoadVariable();

   TA2DataManager::LoadVariable("MaxTAGG",      &fMaxTAGG,      EISingleX);
   TA2DataManager::LoadVariable("MaxCB",        &fMaxCB,        EISingleX);
   TA2DataManager::LoadVariable("MaxTAPS",      &fMaxTAPS,      EISingleX);

   TA2DataManager::LoadVariable("MaxPID",       &fMaxPID,       EISingleX);
   TA2DataManager::LoadVariable("ScalCurr",     fScalCurr,      EDScalerX);
   TA2DataManager::LoadVariable("ScalAcc",      fScalAcc,       EDScalerX);
   TA2DataManager::LoadVariable("ScalAccCorr",  fScalAccCorr,   EDScalerX);

  return;

}

//
// Fucntion for each event
//

void TA2GammaDeuterium::Reconstruct() {
 

   DecodeLadderINFO();

   if ( (fTAGG) && ((Int_t)fTAGG->GetNparticle() > 0 ) ) GetTaggerINFO();    
   if ( (fCB)   && ((Int_t)fCB->GetNCluster() > 0) ) GetCrystalBallINFO();
   if ( (fPID)  && ((Int_t)fPID->GetNhits() > 0) )   GetPidINFO();
   if ( (fMWPC) && ((Int_t)fMWPC->GetNchamber() > 0) ) GetMWPCINFO();
   if ( (fTAPS) && ((Int_t)fTAPS->GetNparticle() > 0) ) GetTAPSINFO();

   if ( ! gAR->IsOnline() ) GetMCINFO();


// Store apparatus INFO

   fGammaDTree->Fill();

   fNTag = 0; fNCBcl = 0; fNTAPScl = 0; fNPid = 0; 
   fNinter1 = 0; fNinter2 = 0;

}

// Decode data coming from counters (seans March 2013)

void TA2GammaDeuterium::DecodeLadderINFO() {

   Double_t corrRatio = 1.;

   if ( (fLADD->IsScaler()) && (gAR->IsScalerRead()) &&
        (gAN->GetNEvent() > 1000) ) {

//     allTag = (Double_t)fScaler[191];
     Double_t corCB = (Double_t)fScaler[190];
     Double_t corTag = (Double_t)fScaler[186];
     if (corTag != 0) corrRatio = corCB/corTag;

     for (Int_t i = 0; i < (Int_t)fLADD->GetNelem(); i++) {
      fScalCurr[i] = (Double_t)fScaler[fScalerIndex[i] + fScalerOffset];
      fScalAcc[i] += fScalCurr[i];
      fScalAccCorr[i] += corrRatio*fScalCurr[i];

     }


   }

}

// Function to extract MC information 


void TA2GammaDeuterium::GetMCINFO() {

    fNpart = *(Int_t*)fEvent[EI_npart];

    Int_t* ipart = (Int_t*)fEvent[EI_idpart];  
    Float_t* vert = (Float_t*)fEvent[EI_vertex];
    Float_t* ptot = (Float_t*)fEvent[EI_plab];
    Float_t* etot = (Float_t*)fEvent[EI_elab];
    Float_t* beam = (Float_t*)fEvent[EI_beam];

    fVertex[0] = (Double_t)vert[0]; 
    fVertex[1] = (Double_t)vert[1]; 
    fVertex[2] = (Double_t)vert[2];

    fEbeam[0] = (Double_t)(beam[0]*beam[3]);
    fEbeam[1] = (Double_t)(beam[1]*beam[3]);
    fEbeam[2] = (Double_t)(beam[2]*beam[3]);

    TVector3 vDir;
    Float_t *dircos = (Float_t*)(fEvent[EI_dircos]);
    for (int i = 0; i < fNpart; i++) {
        fIdPart[i] = ipart[i];
        fPartT[i] = (Double_t)etot[i];
	vDir = TVector3(dircos + i*3);
        fPartX[i] = (Double_t)(ptot[i]*vDir[0]);
        fPartY[i] = (Double_t)(ptot[i]*vDir[1]);
        fPartZ[i] = (Double_t)(ptot[i]*vDir[2]);
    }

 }


// Function to store CB parameters

void TA2GammaDeuterium::GetTaggerINFO() {


  Int_t* FPChannel;  
  Double_t *FPTime;
  const Double_t* FPEnergy;  

  if (! fNMultTag) { 

   fNTag = fLADD->GetNhits();
   for (Int_t i = 0; i < (Int_t)fNTag; i++) {
    fTagCh[i] = fLADD->GetHits()[i];
    fTagEnergy[i] = fBeamEnergy - fLADD->GetEelecOR()[i];
    fTagTime[i] = fLADD->GetTimeOR()[i];  
 
   }
  
  }
 
  else {

   Int_t nStepSum = 0;
   for (Int_t i = 0; i < (Int_t)fNMultTag; i++) {

    FPChannel = fLADD->GetHitsM(i);
    FPEnergy = fLADD->GetECalibration();
    FPTime   = fLADD->GetTimeORM(i);  
 
    nStepSum = (Int_t)fLADD->GetNhitsM(i);
    for (Int_t j = 0; j < nStepSum; j++) {

     fTagCh[fNTag+j]     = FPChannel[j];
     fTagEnergy[fNTag+j] = fBeamEnergy - FPEnergy[FPChannel[j]];
     fTagTime[fNTag+j]   = FPTime[j];  
     if (i == 0) fTagTime[fNTag + j] = fLADD->GetTimeOR()[j];

    }

    fNTag += nStepSum;

   } 

  }
 

}

// Function to store selected CB parameters

void TA2GammaDeuterium::GetCrystalBallINFO() {

   Double_t* cTheta  = fCB->GetTheta();
   Double_t* cPhi    = fCB->GetPhi();
   Double_t* cEnergy = fCB->GetClEnergyOR();
   Double_t* cTime   = fCB->GetClTimeOR();
   Double_t* cRadius = fCB->GetClRadiusOR();

   fNCBcl = fCB->GetNCluster();
   for (Int_t i = 0; i < (Int_t)fNCBcl; i++) {
    fCBEnergy[i] = cEnergy[i];
    fCBTheta[i]  = cTheta[i];
    fCBPhi[i]    = cPhi[i];
    fCBTime[i]   = cTime[i];
    fCBRadius[i] = cRadius[i];
   }
  
}


// Function to store selected PID parameters

void TA2GammaDeuterium::GetPidINFO() {

   fNPid = fPID->GetNhits();
   Int_t* PIDHits = fPID->GetHits();

   for (Int_t i = 0; i < (Int_t)fNPid; i++) {

    fPidCh[i]     = PIDHits[i];
    fPidEnergy[i] = fPID->GetEnergy(PIDHits[i]);
    fPidTime[i]   = fPID->GetTime(PIDHits[i]);
    fPidPhi[i]    = fPID->GetPosition(PIDHits[i])->Z();

   }

}

// Function to store selected MWPC parameters

void TA2GammaDeuterium::GetMWPCINFO() {


  for (Int_t i = 0; i < (Int_t)fMWPC->GetNChamber(); i++) {

   if (i == 0) {
     fNinter1 = fMWPC->GetNinters(i);
     for (Int_t j = 0; j < (Int_t)fNinter1; j++) {
       const TA2MwpcIntersection *inter1 = fMWPC->GetInters(i,j);
       fInter1X[j] = inter1->GetPosition()->X();
       fInter1Y[j] = inter1->GetPosition()->Y();
       fInter1Z[j] = inter1->GetPosition()->Z();
       fInter1Phi[j] = inter1->GetPhi()*TMath::RadToDeg();
       fInter1Type[j] = (Int_t)inter1->GetType();
     }
    }

   if (i == 1) {
     fNinter2 = fMWPC->GetNinters(i);
     for (Int_t j = 0; j < (Int_t)fNinter2; j++) {
       const TA2MwpcIntersection *inter2 = fMWPC->GetInters(i,j);
       fInter2X[j] = inter2->GetPosition()->X();
       fInter2Y[j] = inter2->GetPosition()->Y();
       fInter2Z[j] = inter2->GetPosition()->Z();
       fInter2Phi[j] = inter2->GetPhi()* TMath::RadToDeg();
       fInter2Type[j] = (Int_t)inter2->GetType();
     }
    }

   }


}

// Function to store selected TAPS parameters

void TA2GammaDeuterium::GetTAPSINFO() {


// TAPS clusters from BaF2

   TA2ClusterDetector* fCl = fTAPS->GetCal();
   Bool_t* cVeto  = fTAPS->GetfIsVCharged();

   Double_t* cTheta  = fCl->GetTheta();
   Double_t* cPhi    = fCl->GetPhi();
   Double_t* cEnergy = fCl->GetClEnergyOR();
   Double_t* cTime   = fCl->GetClTimeOR();
   Double_t* cRadius = fCl->GetClRadiusOR(); 

   fNTAPScl = fCl->GetNCluster();

   for (Int_t i = 0; i < (Int_t)fNTAPScl; i++) {
    fTAPSEnergy[i] = cEnergy[i];
    fTAPSTheta[i]  = cTheta[i];
    fTAPSPhi[i]    = cPhi[i];
    fTAPSTime[i]   = cTime[i];
    fTAPSRadius[i] = cRadius[i];
    fTAPSVeto[i]   = (Int_t)cVeto[i];
   }

}
