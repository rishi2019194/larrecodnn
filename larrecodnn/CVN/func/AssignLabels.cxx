#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/LArTrainingData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <iomanip>
#include <iostream>
#include <limits>

namespace lcvn {
  /// Default constructor
  AssignLabels::AssignLabels()
    : nProton(0), nPion(0), nPizero(0), nNeutron(0), pdgCode(0), tauMode(0)
  {}

  /// Get Interaction_t from pdg, mode and iscc.
  /// Setting pdg and mode to zero triggers cosmic ray
  InteractionType AssignLabels::GetInteractionType(simb::MCNeutrino& truth) const
  {

    int pdg = truth.Nu().PdgCode();
    bool iscc = truth.CCNC() == simb::kCC;
    int trueMode = truth.Mode();

    if (iscc) {
      if (abs(pdg) == 14) {
        switch (trueMode) {
        case simb::kQE: return kNumuQE; break;
        case simb::kRes: return kNumuRes; break;
        case simb::kDIS: return kNumuDIS; break;
        default: return kNumuOther;
        }
      }
      else if (abs(pdg) == 12) {
        switch (trueMode) {
        case simb::kQE: return kNueQE; break;
        case simb::kRes: return kNueRes; break;
        case simb::kDIS: return kNueDIS; break;
        default: return kNueOther;
        }
      }
      else if (abs(pdg) == 16) {
        switch (trueMode) {
        case simb::kQE: return kNutauQE; break;
        case simb::kRes: return kNutauRes; break;
        case simb::kDIS: return kNutauDIS; break;
        default: return kNutauOther;
        }
      }
    }
    else if (trueMode == simb::kNuElectronElastic) {
      return kNuElectronElastic;
    }

    return kNC;
  }

  InteractionType AssignLabels::GetInteractionTypeFromSlice(int pdg, bool iscc, int trueMode) const
  {

    if (iscc) {

      if (abs(pdg) == 14) {
        switch (trueMode) {
        case simb::kQE: return kNumuQE;
        case simb::kRes: return kNumuRes;
        case simb::kDIS: return kNumuDIS;
        default: return kNumuOther;
        }
      }
      else if (abs(pdg) == 12) {
        switch (trueMode) {
        case simb::kQE: return kNueQE;
        case simb::kRes: return kNueRes;
        case simb::kDIS: return kNueDIS;
        default: return kNueOther;
        }
      }
      else if (abs(pdg) == 16) {
        switch (trueMode) {
        case simb::kQE: return kNutauQE;
        case simb::kRes: return kNutauRes;
        case simb::kDIS: return kNutauDIS;
        default: return kNutauOther;
        }
      }
    }
    else if (trueMode == simb::kNuElectronElastic) {
      return kNuElectronElastic;
    }

    return kNC;
  }

  // This function uses purely the information from the neutrino generator to
  // find all of the final-state particles that contribute to the event.
  void AssignLabels::GetTopology(const art::Ptr<simb::MCTruth> truth,
                                 unsigned int nTopologyHits = 0)
  {

    const simb::MCNeutrino& nu = truth->GetNeutrino();

    // Get neutrino flavour
    if (nu.CCNC() == simb::kCC) { pdgCode = nu.Nu().PdgCode(); }
    else {
      pdgCode = 1;
    }

    // Get tau topology, if necessary
    tauMode = kNotNutau;
    if (abs(pdgCode) == 16) {
      tauMode = kNutauHad;
      for (int p = 0; p < truth->NParticles(); ++p) {
        if (truth->GetParticle(p).StatusCode() != 1) continue;
        int pdg = abs(truth->GetParticle(p).PdgCode());
        int parent = truth->GetParticle(p).Mother();
        while (parent > 0)
          parent = truth->GetParticle(parent).Mother();

        if (parent == 0) {
          if (pdg == 11) {
            tauMode = kNutauE;
            break;
          }
          else if (pdg == 13) {
            tauMode = kNutauMu;
            break;
          }
        }
      }
    }

    // Now we need to do some final state particle counting.
    //    unsigned int nParticle = truth.NParticles();

    nProton = 0;
    nPion = 0; // Charged pions, that is
    nPizero = 0;
    nNeutron = 0;

    // We need an instance of the backtracker to find the number of simulated hits for each track
    art::ServiceHandle<cheat::BackTrackerService> backTrack;
    art::ServiceHandle<cheat::ParticleInventoryService> partService;

    // Loop over all of the particles
    for (auto const thisPart : partService->MCTruthToParticles_Ps(truth)) {

      const simb::MCParticle& part = *thisPart;

      int pdg = part.PdgCode();

      // Make sure this is a final state particle
      if (part.StatusCode() != 1) { continue; }

      // Make sure this particle is a daughter of the neutrino
      if (part.Mother() != 0) { continue; }

      // GENIE has some fake particles for energy conservation - eg nuclear binding energy. Ignore these
      if (pdg > 2000000000) { continue; }

      // Also don't care about nuclear recoils
      if (pdg > 1000000) { continue; }

      // Find how many SimIDEs the track has
      unsigned int nSimIDE = backTrack->TrackIdToSimIDEs_Ps(part.TrackId()).size();

      // Check if we have more than 100 MeV of kinetic energy
      // float ke = part.E() - part.Mass();
      //    if( ke < 0.0){
      //      continue;
      //    }

      // Special case for pi-zeros since it is the decay photons and their pair produced electrons that deposit energy
      if (pdg == 111 || pdg == 2112) {
        // Decay photons
        for (int d = 0; d < part.NumberDaughters(); ++d) {
          nSimIDE += backTrack->TrackIdToSimIDEs_Ps(part.Daughter(d)).size();
        }
      }

      // Do we pass the number of hits cut?
      if (nSimIDE < nTopologyHits) { continue; }

      switch (abs(pdg)) {
      case 111: ++nPizero; break;
      case 211: ++nPion; break;
      case 2112: ++nNeutron; break;
      case 2212: ++nProton; break;
      default: break;
      }
    }

    std::cout << "Particle counts: " << nProton << ", " << nPion << ", " << nPizero << ", "
              << nNeutron << std::endl;
  }

  void AssignLabels::PrintTopology()
  {

    std::cout << "== Topology Information ==" << std::endl;

    std::cout << " - Neutrino PDG code = " << pdgCode << std::endl;

    std::cout << " - Number of protons (3 means >2) = " << nProton << std::endl;

    std::cout << " - Number of charged pions (3 means >2) = " << nPion << std::endl;

    std::cout << " - Number of pizeros (3 means >2) = " << nPizero << std::endl;

    std::cout << " - Number of neutrons (3 means >2) = " << nNeutron << std::endl;

    std::cout << " - Topology type is " << GetTopologyType() << std::endl;

    std::cout << " - Alternate topology type is " << GetTopologyTypeAlt() << std::endl;
  }

  unsigned short AssignLabels::GetTopologyType() const
  {

    if (abs(pdgCode) == 12) return kTopNue;
    if (abs(pdgCode) == 14) return kTopNumu;
    if (abs(pdgCode) == 16) {
      if (tauMode == kNutauE) return kTopNutauE;
      if (tauMode == kNutauMu) return kTopNutauMu;
      if (tauMode == kNutauHad) return kTopNutauHad;
    }
    if (pdgCode == 1) return kTopNC;
    throw std::runtime_error("Topology type not recognised!");
  }

  unsigned short AssignLabels::GetTopologyTypeAlt() const
  {

    if (abs(pdgCode) == 12) return kTopNueLike;
    if (abs(pdgCode) == 14) return kTopNumuLike;
    if (abs(pdgCode) == 16) {
      if (tauMode == kNutauE) return kTopNueLike;
      if (tauMode == kNutauMu) return kTopNumuLike;
      if (tauMode == kNutauHad) return kTopNutauLike;
    }
    if (pdgCode == 1) return kTopNCLike;
    throw std::runtime_error("Topology type not recognised!");
  }

  // Get the beam interaction mode for ProtoDUNE specific code
  unsigned short AssignLabels::GetProtoDUNEBeamInteractionType(
    const simb::MCParticle& particle) const
  {

    unsigned short baseProcess = std::numeric_limits<unsigned short>::max();

    // The first thing we can do is look at the process key
    std::string processName = particle.EndProcess();

    if (GetProcessKey(processName) > -1) {
      // Base process gives us a value from 0 to 44
      baseProcess = static_cast<unsigned int>(GetProcessKey(processName));
    }

    std::cout << "What interaction type, then? " << processName << std::endl;

    // In the case that we have an inelastic interaction, maybe we can do more.
    art::ServiceHandle<cheat::ParticleInventoryService> piService;

    unsigned int nPi0 = 0; // Pi-zeros
    unsigned int nPiM = 0; // Pi-minuses
    unsigned int nPiP = 0; // Pi-pluses
    unsigned int nNeu = 0; // Neutrons
    unsigned int nPro = 0; // Protons
    unsigned int nOth = 0; // Everything else

    for (int i = 0; i < particle.NumberDaughters(); ++i) {
      const simb::MCParticle* daughter = piService->TrackIdToParticle_P(particle.Daughter(i));
      switch (daughter->PdgCode()) {
      case 111: ++nPi0; break;
      case -211: ++nPiM; break;
      case 211: ++nPiP; break;
      case 2112: ++nNeu; break;
      case 2212: ++nPro; break;
      default: ++nOth; break;
      }
    }

    std::cout << "Base process = " << baseProcess << std::endl;
    std::cout << "Daughters = " << nPi0 << " pi0s, " << nPiM << " pi-s, " << nPiP << " pi+s, "
              << nNeu << " neutrons, " << nPro << " protons and " << nOth << " other particles."
              << std::endl;

    // If we have a pion with a pi0 in the final state we can flag it as charge exchange
    if (abs(particle.PdgCode()) == 211 && nPi0 == 1) {
      return 45; // First free value after those from the truth utility
    }
    else {
      return baseProcess;
    }
  }

  // Get process key.
  int lcvn::AssignLabels::GetProcessKey(std::string process) const
  {

    if (process.compare("primary") == 0) return 0;
    if (process.compare("hadElastic") == 0) return 1;
    if (process.compare("pi-Inelastic") == 0) return 2;
    if (process.compare("pi+Inelastic") == 0) return 3;
    if (process.compare("kaon-Inelastic") == 0) return 4;
    if (process.compare("kaon+Inelastic") == 0) return 5;
    if (process.compare("protonInelastic") == 0) return 6;
    if (process.compare("neutronInelastic") == 0) return 7;
    if (process.compare("kaon0SInelastic") == 0) return 8;
    if (process.compare("kaon0LInelastic") == 0) return 9;
    if (process.compare("lambdaInelastic") == 0) return 10;
    if (process.compare("omega-Inelastic") == 0) return 11;
    if (process.compare("sigma+Inelastic") == 0) return 12;
    if (process.compare("sigma-Inelastic") == 0) return 13;
    if (process.compare("sigma0Inelastic") == 0) return 14;
    if (process.compare("xi-Inelastic") == 0) return 15;
    if (process.compare("xi0Inelastic") == 0) return 16;
    if (process.compare("anti_protonInelastic") == 0) return 20;
    if (process.compare("anti_neutronInelastic") == 0) return 21;
    if (process.compare("anti_lambdaInelastic") == 0) return 22;
    if (process.compare("anti_omega-Inelastic") == 0) return 23;
    if (process.compare("anti_sigma+Inelastic") == 0) return 24;
    if (process.compare("anti_sigma-Inelastic") == 0) return 25;
    if (process.compare("anti_xi-Inelastic") == 0) return 26;
    if (process.compare("anti_xi0Inelastic") == 0) return 27;

    if (process.compare("Decay") == 0) return 30;
    if (process.compare("FastScintillation") == 0) return 31;
    if (process.compare("nKiller") == 0)
      return 32; // Remove unwanted neutrons: neutron kinetic energy threshold (default 0) or time limit for neutron track
    if (process.compare("nCapture") == 0) return 33; // Neutron capture

    if (process.compare("compt") == 0) return 40;                 // Compton Scattering
    if (process.compare("rayleigh") == 0) return 41;              // Rayleigh Scattering
    if (process.compare("phot") == 0) return 42;                  // Photoelectric Effect
    if (process.compare("conv") == 0) return 43;                  // Pair production
    if (process.compare("CoupledTransportation") == 0) return 44; //

    return -1;
  }

  // Recursive function to get all hits from daughters of a given neutral particle
  unsigned int AssignLabels::GetNeutralDaughterHitsRecursive(const simb::MCParticle& particle) const
  {

    unsigned int nSimIDEs = 0;

    // The backtrack and particle inventory service will be useful here
    art::ServiceHandle<cheat::BackTrackerService> backTrack;
    art::ServiceHandle<cheat::ParticleInventoryService> partService;

    for (int d = 0; d < particle.NumberDaughters(); ++d) {

      const simb::MCParticle* daughter = partService->TrackIdToParticle_P(particle.Daughter(d));
      unsigned int localSimIDEs = backTrack->TrackIdToSimIDEs_Ps(daughter->TrackId()).size();
      std::cout << "Got " << localSimIDEs << " hits from " << daughter->PdgCode() << std::endl;
      if (localSimIDEs == 0) localSimIDEs = GetNeutralDaughterHitsRecursive(*daughter);

      nSimIDEs += localSimIDEs;
    }

    return nSimIDEs;
  }

}
