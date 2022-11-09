////////////////////////////////////////////////////////////////////////
// \file    Assignlabels.h
///\brief   Utility class for truth labels
///
// \author  Leigh Whitehead leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////
#ifndef LCVN_ASSIGNLABELS_H
#define LCVN_ASSIGNLABELS_H

#include "art/Framework/Principal/Handle.h"

#include "larrecodnn/CVN/func/InteractionType.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace lcvn {

  class AssignLabels {

  public:
    AssignLabels();

    InteractionType GetInteractionType(simb::MCNeutrino& truth) const;
    InteractionType GetInteractionTypeFromSlice(int nuPDG, bool nuCCNC, int nuMode) const;

    // Use the topology information
    void GetTopology(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits);
    void PrintTopology();
    unsigned short GetNProtons() const { return nProton; };
    unsigned short GetNPions() const { return nPion; };
    unsigned short GetNPizeros() const { return nPizero; };
    unsigned short GetNNeutrons() const { return nNeutron; };
    short GetPDG() const { return pdgCode; };
    unsigned short TauMode() const { return tauMode; };
    bool IsAntineutrino() const { return pdgCode < 0; };
    unsigned short GetTopologyType() const;
    unsigned short GetTopologyTypeAlt() const;

    // Get the pion interaction mode for ProtoDUNE specific code
    unsigned short GetProtoDUNEBeamInteractionType(const simb::MCParticle& particle) const;

  private:
    // Recursive function to get all hits from daughters of a neutral particle
    unsigned int GetNeutralDaughterHitsRecursive(const simb::MCParticle& particle) const;

    int GetProcessKey(std::string process) const;

    unsigned short nProton;
    unsigned short nPion;
    unsigned short nPizero;
    unsigned short nNeutron;
    short pdgCode;
    unsigned short tauMode;
  };
}

#endif // CVN_ASSIGNLABELS_H
