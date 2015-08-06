// SVN Info: $Id: TOA2RecParticle.cxx 1257 2012-07-26 15:33:13Z werthm $

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


#include "TOA2RecParticle.h"

ClassImp(TOA2RecParticle)


//______________________________________________________________________________
TOA2RecParticle::TOA2RecParticle()
    : TObject()
{
    // Empty constructor.
    
    // do not save bit mask and unique ID
    Class()->IgnoreTObjectStreamer();
 
    fP4 = 0;

    fNDetectedProducts = 0;
    fDetectedProducts  = 0;
    fChiSquare         = 0.;
}

//______________________________________________________________________________
TOA2RecParticle::TOA2RecParticle(Int_t nParticle)
    : TObject()
{
    // Pseudo-non-empty constructor (nParticle is not really used here).
    
    // do not save bit mask and unique ID
    Class()->IgnoreTObjectStreamer();
 
    // create members
    fP4 = new TLorentzVector();

    fNDetectedProducts = 0;
    fDetectedProducts  = 0;
    fChiSquare         = 0.;
}

//______________________________________________________________________________
TOA2RecParticle::TOA2RecParticle(const TOA2RecParticle& orig)
    : TObject(orig)
{
    // Copy constructor.
    
    // do not save bit mask and unique ID
    Class()->IgnoreTObjectStreamer();
 
    // create members
    fP4 = new TLorentzVector(*orig.fP4);
    
    fNDetectedProducts = orig.fNDetectedProducts;
    fChiSquare         = orig.fChiSquare;
    
    fDetectedProducts = new TOA2DetParticle*[fNDetectedProducts];
    for (Int_t i = 0; i < fNDetectedProducts; i++) fDetectedProducts[i] = new TOA2DetParticle(*orig.GetDetectedProduct(i));
}

//______________________________________________________________________________
TOA2RecParticle::~TOA2RecParticle()
{
    // Destructor.
 
    Clear();
}

//______________________________________________________________________________
void TOA2RecParticle::Clear(Option_t* option)
{
    // Clear the memory.
    
    if (fP4) delete fP4;
    if (fDetectedProducts) 
    {
        for (Int_t i = 0; i < fNDetectedProducts; i++) delete fDetectedProducts[i];
        delete [] fDetectedProducts;
    }
}

//______________________________________________________________________________
Bool_t TOA2RecParticle::Reconstruct(Int_t nParticle, TOA2DetParticle** particleList)
{
    // Basic reconstruction of the particle from its detected decay products.
    // Get the mass of the detected decay products, calculate the 4-vectors
    // and combine them to this particle's 4-vector.
    //
    // Return kTRUE if the reconstruction was successful.
    //
    // Overwrite this method if you need a more sophisticated reconstruction.
    
    // check the number of particles
    if (!nParticle)
    {
        fNDetectedProducts = 0;
        fDetectedProducts  = 0;
        Error("Reconstruct", "Number of particles cannot be zero!");
        return kFALSE;
    }
    // copy the particles directly to the detected decay product list
    // and mark them as gammas
    else
    {
        fNDetectedProducts = nParticle;
        fDetectedProducts = new TOA2DetParticle*[fNDetectedProducts];
        for (Int_t i = 0; i < fNDetectedProducts; i++) 
        {
            fDetectedProducts[i] = new TOA2DetParticle(*particleList[i]);
            fDetectedProducts[i]->SetPDG_ID(kGamma);
        }

        // calculate the 4-vector
        Calculate4Vector();

        return kTRUE;
    }
}

//______________________________________________________________________________
void TOA2RecParticle::Calculate4Vector()
{
    // Basic method to calculate the 4-vector of the reconstructed particle
    // from the 4-vectors of the detected products.
    
    // clear particle 4-vector
    fP4->SetPxPyPzE(0., 0., 0., 0.);

    // loop over decay products
    for (Int_t i = 0; i < fNDetectedProducts; i++)
    {
        // determine mass of decay product
        Double_t mass;
        TParticlePDG* pdg_pa = TOGlobals::GetPDG_DB()->GetParticle(fDetectedProducts[i]->GetPDG_ID());
        
        // check if particle is in database
        if (pdg_pa) mass = pdg_pa->Mass();
        else
        {
            mass = 0.;
            Error("Calculate4Vector", "Particle with code %d was not found in the PDG database!", 
                  (Int_t)fDetectedProducts[i]->GetPDG_ID());
        }
        
        // build 4-vector of decay product
        TLorentzVector p4;
        fDetectedProducts[i]->Calculate4Vector(&p4, mass);

        // add 4-vector to total 4-vector
        *fP4 += p4;
    }
}

//______________________________________________________________________________
void TOA2RecParticle::CalculateConstrained4Vector(Double_t mass, TLorentzVector* outP4)
{
    // Calculate the 4-vector constrained to the mass 'mass' and save it in 'outP4.
    // If the particle was not reconstructed yet 0 is returned.
    
    if (!fP4->M()) 
    {
        Error("CalculateConstrained4Vector", "Cannot constrain 4-vector because particle mass is zero!");
        return;
    }

    /*
    //
    // gives the same result as the code below
    //
    Double_t ptot_old = (fP4->Vect()).Mag();
    Double_t energy = mass / fP4->M() * fP4->E();

    Double_t ptot_new = TMath::Sqrt(energy * energy - mass * mass);
    Double_t scale = ptot_new / ptot_old;

    fP4Const->SetPxPyPzE(fP4->Px() * scale, fP4->Py() * scale, fP4->Pz() * scale, energy);
    */
    
    Double_t scale = mass / fP4->M();
    outP4->SetPxPyPzE(fP4->Px()*scale, fP4->Py()*scale, fP4->Pz()*scale, fP4->E()*scale);
}

//______________________________________________________________________________
Double_t TOA2RecParticle::GetAverageTime(A2Detector_t detector) const
{
    // Return the average time of the decay products detected in the A2 detector
    // 'detector'. 
    //
    // If no decay products are set or no decay products are found int the 
    // specified detector 9999 is returned.

    if (fNDetectedProducts)
    {
        Double_t avr_time = 0.;
        Int_t n = 0;

        // loop over decay products
        for (Int_t i = 0; i < fNDetectedProducts; i++)
        {
            // check detector
            if (fDetectedProducts[i]->GetDetector() == detector) 
            {
                avr_time += fDetectedProducts[i]->GetTime();
                n++;
            }
        }

        if (!n) return 9999;
        else return avr_time / (Double_t) n;
    }
    else 
    {   
        Error("GetAverageTime", "No detected products defined!");
        return 9999;
    }
}

//______________________________________________________________________________
Int_t TOA2RecParticle::GetNDetectorHits(A2Detector_t detector) const
{
    // Return the number of detected decay products found in the detector 'detector'.

    Int_t n = 0;

    for (Int_t i = 0; i < fNDetectedProducts; i++)
    {
        if (fDetectedProducts[i]->GetDetector() == detector) n++;
    }
    
    return n;
}

//______________________________________________________________________________
Bool_t TOA2RecParticle::CheckCoincidence(Double_t limitLowCB, Double_t limitHighCB, TH1* hCB,
                                         Double_t limitLowTAPS, Double_t limitHighTAPS, TH1* hTAPS,
                                         Double_t limitLowCBTAPS, Double_t limitHighCBTAPS, TH1* hCBTAPS)
{
    // Check the coincidence of all detected particles. Use
    // - 'limitLowCB' / 'limitHighCB' as limits for coincident particles in CB
    // - 'limitLowTAPS' / 'limitHighTAPS' as limits for coincident particles in TAPS
    // - 'limitLowCBTAPS' / 'limitHighCBTAPS' as limits for the difference of the average 
    // times of particles in CB and TAPS, respectively.
    // If 'hCB', 'hTAPS' or 'hCBTAPS' are provided the time differences will be filled into these histograms.

    Bool_t coincOk = kTRUE;

    // check CB coincidence
    if (!CheckDetectorCoincidence(kCBDetector, limitLowCB, limitHighCB, hCB)) coincOk = kFALSE;

    // check TAPS coincidence
    if (!CheckDetectorCoincidence(kTAPSDetector, limitLowTAPS, limitHighTAPS, hTAPS)) coincOk = kFALSE;

    // check CB-TAPS coincidence if there is a hit in TAPS
    Int_t nTAPSHits = GetNDetectorHits(kTAPSDetector);
    if (nTAPSHits)
    {
        Double_t avrCBTime = GetAverageTime(kCBDetector);
        Double_t avrTAPSTime = GetAverageTime(kTAPSDetector);

        // the time values have to be added to remove the trigger time
        Double_t avrCBTAPSTime = avrCBTime + avrTAPSTime;
        hCBTAPS->Fill(avrCBTAPSTime);
        
        if (avrCBTAPSTime < limitLowCBTAPS || avrCBTAPSTime > limitHighCBTAPS) coincOk = kFALSE;
    }

    return coincOk;
}

//______________________________________________________________________________
Bool_t TOA2RecParticle::CheckDetectorCoincidence(A2Detector_t detector, Double_t limitLow, Double_t limitHigh, TH1* h)
{
    // Check if all combinations of particles detected in the detector 'detector' fulfill the coincidence
    // condition, i.e. their time difference lies between 'limitLow' and 'limitHigh'.
    // If 'h' is provided the time differences will be filled into this histogram.
    // 
    // Note: 
    // If no particles are found in the specified detector kTRUE is returned.

    Bool_t coincOk = kTRUE;

    // detected particle loop 1
    for (Int_t i = 0; i < fNDetectedProducts; i++)
    {
        // check the detector
        if (fDetectedProducts[i]->GetDetector() != detector) continue;

        // detected particle loop 2
        for (Int_t j = i+1; j < fNDetectedProducts; j++)
        {
            // check the detector
            if (fDetectedProducts[j]->GetDetector() != detector) continue;
    
            // calculate time difference
            Double_t timeDiff = fDetectedProducts[i]->GetTime() - fDetectedProducts[j]->GetTime();

            // fill histogram
            if (h) h->Fill(timeDiff);

            // check coincidence
            if (timeDiff < limitLow || timeDiff > limitHigh) coincOk = kFALSE;
        }
    }

    return coincOk;
}

//______________________________________________________________________________
Double_t TOA2RecParticle::GetTime() const
{
    // Return the time of the particle that is calculated by averaging the
    // times of the detected products. If there are products detected in TAPS
    // this time will be preferred because of its better resolution.
    
    // check if there are particles in TAPS
    if (GetNDetectorHits(kTAPSDetector)) return GetAverageTime(kTAPSDetector);
    else return GetAverageTime(kCBDetector);
}

//______________________________________________________________________________
TOA2DetParticle* TOA2RecParticle::GetDetectedProduct(Int_t index) const
{
    // Return the decay product with index 'index'. Return 0 if there is no
    // product with this index in the products list.

    // check index
    if (index < 0 || index >= fNDetectedProducts) 
    {
        Error("GetDetectedProduct", "Index %d out of bounds (size: %d)", index, fNDetectedProducts);
        return 0;
    }

    // return particle
    return fDetectedProducts[index];
}

//______________________________________________________________________________
void TOA2RecParticle::Print(Option_t* option) const
{
    // Print out the content of this class.

    printf("TOA2RecParticle content:\n");
    if (fP4) printf("4-vector               : %f %f %f %f %f\n", fP4->Px(), fP4->Py(),fP4->Pz(),fP4->E(), fP4->M());
    else printf("4-vector               : zero\n");
    printf("chi square (reconstr.) : %f\n", fChiSquare);

    printf("detected products      : %d\n", fNDetectedProducts);
    if (fDetectedProducts) 
    {   
        for (Int_t i = 0; i < fNDetectedProducts; i++) 
        {
            printf("\n");
            printf("product %d\n", i+1);
            fDetectedProducts[i]->Print();
        }
        printf("\n");
    }
}

TOA2RecParticle operator+(const TOA2RecParticle& part1, const TOA2RecParticle& part2)
{
    // Add two reconstructed particles.
    //
    // NOTE: Currently only the 4-vector and the list of detected particles are
    //       added.

    // create the object
    TOA2RecParticle out(0);

    // set the 4-vector
    TLorentzVector sum = *part1.Get4Vector() + *part2.Get4Vector();
    out.Set4Vector(&sum);
    
    //
    // set the detected products
    //

    // get the number of products, set the total number of products
    Int_t nPart1 = part1.GetNDetectedProducts();
    Int_t nPart2 = part2.GetNDetectedProducts();
    out.SetNDetectedProducts(nPart1+nPart2);
    
    // create the new product list and copy the products over there
    TOA2DetParticle** fDetectedProducts = new TOA2DetParticle*[nPart1+nPart2];
    for (Int_t i = 0; i < nPart1; i++) fDetectedProducts[i] = new TOA2DetParticle(*part1.GetDetectedProduct(i));
    for (Int_t i = 0; i < nPart2; i++) fDetectedProducts[nPart1+i] = new TOA2DetParticle(*part2.GetDetectedProduct(i));
    
    // set the list in the new object
    out.SetDetectedProducts(fDetectedProducts);

    return out;
}

