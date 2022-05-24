////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com

//#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larrecodnn/CVN/func/TrainingData.h"

namespace cvn
{

  TDNuInfo::TDNuInfo():
    fNuEnergy(0.),
    fLepEnergy(0.),
    fLepAngle(0.),
    fEventWeight(1),
    fRecoNueEnergy(0.),
    fRecoNumuEnergy(0.),
    fRecoNutauEnergy(0.),
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

  void TDNuInfo::SetRecoInfo(float nueEnergy, float numuEnergy, float nutauEnergy){

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

  std::ostream& operator<<(std::ostream& os, const TDNuInfo& td)
  {
    // keep this order because the keras python scripts depend on this
    os << td.fNuEnergy << std::endl;
    os << td.fLepEnergy << std::endl;
    os << td.fRecoNueEnergy << std::endl;
    os << td.fRecoNumuEnergy << std::endl;
    os << td.fRecoNutauEnergy << std::endl;
    os << td.fEventWeight << std::endl;
    
    os << td.fNuPDG << std::endl;
    os << td.fNProton << std::endl;
    os << td.fNPion << std::endl;
    os << td.fNPizero << std::endl;
    os << td.fNNeutron << std::endl;

    os << td.fTopologyType << std::endl;
    os << td.fTopologyTypeAlt << std::endl;
    os << td.fLepAngle << std::endl;

    return os;
  }

  //----------------------------------------------------------------------

  template <class T> TrainingData<T>::TrainingData(const InteractionType& interaction,
                             const PixelMap& pMap,
                             const T info):
  fInt(interaction),
  fPMap(pMap),
  fInfo(info)
  {  }


  template <class T> void TrainingData<T>::FillOutputVector(float* output) const
  {
    for(unsigned int i = 0; i < kNIntType; ++i)
      output[i] = 0;

    output[fInt] = 1;
  }

  template class TrainingData<cvn::TDNuInfo>; 
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
