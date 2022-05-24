////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef CVN_TRAININGDATA_H
#define CVN_TRAININGDATA_H

#include <iostream>
#include <ostream>

#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/AssignLabels.h"

namespace cvn
{

  class TDNuInfo
  {
  public:
    TDNuInfo();
    
    void SetTruthInfo(float nuEnergy, float lepEnergy, float lepAngle, float weight);
    void SetRecoInfo(float nueEnergy, float numuEnergy, float nutauEnergy);

    // Set topology information separately to save having a large number of 
    // arguments in the constructor.
    void SetTopologyInformation(int pdg, int nproton, int npion,
                                int npizero, int nneutron, int toptype,
                                int toptypealt);

    friend std::ostream& operator<<(std::ostream& os, const TDNuInfo& td);

  private: 
    float    fNuEnergy;        ///< True energy of neutrino event
    float    fLepEnergy;       ///< True energy of outgoing lepton
    float    fLepAngle;       ///< True angle of outgoing lepton wrt neutrino
    float    fEventWeight;     ///< The event weight (norm * oscProb)
    
    float    fRecoNueEnergy;   ///< Reconstructed energy under nue hypothesis
    float    fRecoNumuEnergy;  ///< Reconstructed energy under numu hypothesis
    float    fRecoNutauEnergy; ///< Reconstructed energy under nutau hypothesis
 
    // If we are using topology information, store it here
    bool fUseTopology;
    int  fNuPDG;
    int  fNProton;
    int  fNPion;
    int  fNPizero;
    int  fNNeutron;
    int  fTopologyType;
    int  fTopologyTypeAlt;
 
  };
  
  /// \brief   The TrainingData objects contains a PixelMap and the
  ///          output class type, and any other bit that goes into the ANN

  template <class T> class TrainingData
  {

  public:
    TrainingData(){};
    TrainingData(const InteractionType& interaction,
                 const PixelMap& pMap,
                 const T info);

    unsigned int NOutput() const {return (unsigned int)kNIntType;};

    void FillOutputVector(float* output) const;

    InteractionType  fInt;     ///< Class of the event
    PixelMap fPMap;           ///< PixelMap for the event
    T fInfo;
  };

  typedef TrainingData<cvn::TDNuInfo> TrainingNuData;
} // end namespace

#endif // CVN_TRAININGDATA_H
//////////////////////////////////////////////////////////////////////////////
