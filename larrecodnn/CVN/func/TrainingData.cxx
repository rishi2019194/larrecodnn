////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com

//#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dunereco/CVN/func/TrainingData.h"

namespace cvn
{

  TDNuInfo::TDNuInfo():
    fNuEnergy(0.),
    fLepEnergy(0.),
    fLepAngle(0.),
    fRecoNueEnergy(0.),
    fRecoNumuEnergy(0.),
    fRecoNutauEnergy(0.),
    fEventWeight(1),
    fUseTopology(false),
    fNuPDG(0),
    fNProton(-1),
    fNPion(-1),
    fNPizero(-1),
    fNNeutron(-1),
    fTopologyType(-1),
    fTopologyTypeAlt(-1)
  { }

  void TDNuInfo::SetTruthInfo(float nuEnergy, float lepEnergy, float lepAngle, float weight){

    fNuEnergy = nuEnergy;
    fLepEnergy = lepEnergy;
    fLepAngle = lepAngle;
    fEventWeight = weight;
  }

  void TDNuInfo::SetRecoInfo(float nueEnergy, numuEnergy, nutauEnergy){

    fRecoNueEnergy = nueEnergy;
    fRecoNumuEnergy = numuEnergy;
    fRecoNutauEnergy = nutauEnergy;
  }

  void TDNuInfo::SetTopologyInformation(int pdg, int nproton,
    int npion, int npizero, int nneutron, int toptype,
    int toptypealt){

    fUseTopology = true;

    fNuPDG = pdg;
    fNProton = nproton;
    fNPion = npion;
    fNPizero = npizero;
    fNNeutron = nneutron;

    fTopologyType = toptype;
    fTopologyTypeAlt = toptypealt;

  }

  //----------------------------------------------------------------------

  TrainingData::TrainingData(const InteractionType& interaction,
                             const PixelMap& pMap,
                             const T info):
  fInt(interaction),
  fPMap(pMap),
  fInfo(info)
  {  }


  void TrainingData::FillOutputVector(float* output) const
  {
    for(unsigned int i = 0; i < kNIntType; ++i)
      output[i] = 0;

    output[fInt] = 1;
  }

} // end namespace cvn
////////////////////////////////////////////////////////////////////////
