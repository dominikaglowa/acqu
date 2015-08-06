// SVN Info: $Id: MyHitClusterTAPS_t.h 337 2008-11-26 12:28:11Z werthm $

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


#ifndef __MyHitClusterTAPS_t_h__
#define __MyHitClusterTAPS_t_h__

#include "TA2ClusterDetector.h"
#include "HitCluster_t.h"


class MyHitClusterTAPS_t : public HitCluster_t
{

private:
    TRandom* fRandGen;                              // RG for random hit positions: activated if non-zero!

    void SetRandomPosition(TVector3* in, TVector3* out);

public:
    MyHitClusterTAPS_t(Char_t* line, UInt_t index, Int_t sizefactor = 1, TRandom* rg = 0);
    virtual ~MyHitClusterTAPS_t();
    
    void ClusterDetermine(TA2ClusterDetector* cl);
    void MoreNeighbours(TA2ClusterDetector* cl, UInt_t det_i, 
                                TVector3* det_pos, Double_t* det_energy);
};

#endif

