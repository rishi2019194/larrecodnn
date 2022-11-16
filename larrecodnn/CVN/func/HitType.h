////////////////////////////////////////////////////////////////////////
// \file    HitType.h
///\brief   Defines an enumeration for cellhit classification
///
// \author psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LCVN_HITTYPE_H
#define LCVN_HITTYPE_H

namespace lcvn {

  enum HitType {
    kElectronHit,
    kMuonHit,
    kProtonHit,
    kNeutronHit,
    kPionHit,
    kPiZeroHit,
    kGammaHit,
    kOtherPDGhit,
    kUnknownHit,
    kEmptyHit
  };

}

#endif // CVN_HITTYPE_H
