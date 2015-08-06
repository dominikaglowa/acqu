// SVN Info: $Id: MyHitClusterTAPS_t.cc 337 2008-11-26 12:28:11Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2008
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyHitClusterTAPS_t                                                   //
//                                                                      //
// My customized cluster reconstruction class for TAPS.                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "MyHitClusterTAPS_t.h"


// BaF2 crystal surface geometry
static const UInt_t TAPS_SurfPolyNPoints = 7;
static const Double_t TAPS_SurfOuterSideLength = 3.5;
static Double_t TAPS_SurfPolyPointsX[7] = { 3.4641014, 1.7320507, -1.7320507, -3.4641014,
                                           -1.7320507, 1.7320507,  3.4641014 };
static Double_t TAPS_SurfPolyPointsY[7] = { 0.0000000,  3.0000000, 3.0000000, 0.0000000,
                                           -3.0000000, -3.0000000, 0.0000000 };


//______________________________________________________________________________
MyHitClusterTAPS_t::MyHitClusterTAPS_t(Char_t* line, UInt_t index, Int_t sizefactor, TRandom* rg)
    : HitCluster_t(line, index, sizefactor)
{
    // Constructor.

    fRandGen = rg;
}

//______________________________________________________________________________
MyHitClusterTAPS_t::~MyHitClusterTAPS_t()
{
    // Destructor.
    
}

//______________________________________________________________________________
inline void MyHitClusterTAPS_t::SetRandomPosition(TVector3* in, TVector3* out)
{
    // Save a random hit position within the crystal surface of the
    // TAPS element with center 'in' to 'out'.
   
    Double_t xRand, yRand;

    // generate random point that fulfills polygon condition
    do
    {
        xRand = 2 * TAPS_SurfOuterSideLength * fRandGen->Rndm() - TAPS_SurfOuterSideLength;
        yRand = 2 * TAPS_SurfOuterSideLength * fRandGen->Rndm() - TAPS_SurfOuterSideLength;
    }
    while (!TMath::IsInside(xRand, yRand, TAPS_SurfPolyNPoints, TAPS_SurfPolyPointsX, TAPS_SurfPolyPointsY));

    out->SetX(in->X() + xRand);
    out->SetY(in->Y() + yRand);
    out->SetZ(in->Z());
}

//______________________________________________________________________________
void MyHitClusterTAPS_t::ClusterDetermine(TA2ClusterDetector* cl)
{
    // Overwrite default ClusterDetermine() for TAPS.
     
    // get detector pointers
    Double_t* energy   = cl->GetEnergy();
    Double_t* time     = cl->GetTime();
    UInt_t* hits       = cl->GetTempHits();
    TVector3** pos     = cl->GetPosition();
    UInt_t nhits       = cl->GetNhits();
    UInt_t max_cluster = cl->GetMaxCluster();

    // define local variables
    UInt_t i, j, k;
    Double_t w, w_tot;
    TVector3 det_pos[max_cluster];
    Double_t det_energy[max_cluster];

    // process central element
    if (fRandGen) SetRandomPosition(pos[fIndex], &det_pos[0]);
    else det_pos[0] = *(pos[fIndex]); 
    det_energy[0] = energy[fIndex];
    fEnergy = det_energy[0];
    if (cl->IsTime()) fTime = time[fIndex];       
    fSqrtEtot = 0.;                                 // not used
    fRadius = 0.;                                   // not used
    fHits[0] = fIndex;

    // Search neighbours
    for (i = 0, k = 1; i < nhits; i++)
    {
        if ((j = hits[i]) == ENullHit) continue;    // was previously counted
        if (j == fIndex)
        {
            hits[i] = ENullHit;
            continue;
        }                                           // already got center
        if (IsNeighbour(j))
        {                                           // a neighbour of the center?
            hits[i] = ENullHit;                     // so its not double counted
            det_pos[k] = *(pos[j]);
            det_energy[k] = energy[j];
            fEnergy += det_energy[k];
            fHits[k] = j;                           // add to cluster hits collection
            k++;
        }
    }
    fNhits = k;
    
    // If the iterate flag is TRUE check if any unused hits (ie NOT in the near
    // neighbour array) are connected to existing hits
    if (cl->IsIterate()) MoreNeighbours(cl, fNhits, det_pos, det_energy);

    fHits[fNhits] = EBufferEnd;                  // mark no more hits
    
    // Calculate weights and mean position
    w_tot = 0.;
    fMeanPosition->SetXYZ(0., 0., 0.);
    for (i = 0; i < fNhits; i++)
    {
        // calculate log. weight
        w = 5. + TMath::Log(det_energy[i] / fEnergy);
        if (w < 0)  w = 0.;
        w_tot += w;

        // update mean position
        *fMeanPosition += (w * det_pos[i]);
    }

    // Normalize mean position
    *fMeanPosition *= 1. / w_tot;

    // correct for the shower depth
    Double_t sh_dep = 2.05 * (TMath::Log(fEnergy / 12.7) + 1.2);    // shower depth
    if (sh_dep > 0)
    {
        Double_t x = fMeanPosition->X();                            // reconstructed x
        Double_t y = fMeanPosition->Y();                            // reconstructed y
        Double_t z = fMeanPosition->Z();                            // target-TAPS distance
        Double_t s = TMath::Sqrt(x*x + y*y + z*z);                  // flight path
        Double_t deltaX = x / (s / sh_dep + 1.);                    // x correction
        Double_t deltaY = y / (s / sh_dep + 1.);                    // y correction
        fMeanPosition->SetX(x - deltaX);                            // correct x-coord.
        fMeanPosition->SetY(y - deltaY);                            // correct y-coord.
    }

    // Calculate the rest of the cluster information
    fCentralFrac = energy[fIndex] / fEnergy;     
    fTheta = TMath::RadToDeg() * fMeanPosition->Theta();
    fPhi   = TMath::RadToDeg() * fMeanPosition->Phi();
}
   
//______________________________________________________________________________
void MyHitClusterTAPS_t::MoreNeighbours(TA2ClusterDetector* cl, UInt_t det_i, 
                                        TVector3* det_pos, Double_t* det_energy)
{
    // Overwrite default MoreNeighbours() for TAPS.

    UInt_t i, j, k, n;
    Double_t distApart;                          // distance between elements
    TVector3** pos = cl->GetPosition();          // element position array
    UInt_t* hits = cl->GetTempHits();            // hits array
    UInt_t nhits = cl->GetNhits();               // # hits
    Double_t* energy = cl->GetEnergy();          // array of energies
    UInt_t* nextIter;                            // -> relevent part fHits array
    UInt_t* hitsEnd = fHits + fMaxHits - 1;      // end-stop of fHits array
    UInt_t pNhits = 0;                           // # previously processed hits

    // Iterate round while new cluster members still found
    do
    {
        nextIter = fHits + fNhits;                          // -> start new cluster members
        for (i = pNhits, n = 0; i < fNhits; i++)            // loop existing clust hits
        {
            for (j = 0; j < nhits; j++)                     // loop total detector hits
            {
	        if ((k = hits[j]) == ENullHit) continue;    // hit already spoken for

	        // distance between cluster hit and potential new hit
	        // if its within limits add the new hit to the cluster
	        distApart = (*(pos[fHits[i]]) - *(pos[k])).Mag();
	        if ((distApart > cl->GetMinPosDiff()) && (distApart < cl->GetMaxPosDiff()) ) 
                {
	            hits[j] = ENullHit;                 // mark hit as spoken for
                    det_pos[det_i] = *(pos[k]);
                    det_energy[det_i] = energy[k];
                    fEnergy += det_energy[det_i];
	            *nextIter++ = k;                    // add index to cluster hits list
	            n++;                                // update # new clust members
	            det_i++;

                    // Check if space for futher cluster members, finish if not
	            if (nextIter >= hitsEnd)
                    {
	                fNhits += n;
	                return;
	            }
	        } 
            }
        }
    
        pNhits = fNhits;         // update previously processed hits
        fNhits += n;             // update total processed hits
    } while (n);                 // iterate while new cluster members found
}

