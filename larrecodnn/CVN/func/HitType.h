////////////////////////////////////////////////////////////////////////
// \file    HitType.h
///\brief   Defines an enumeration for cellhit classification
///
// \author psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LCVN_HITTYPE_H
#define LCVN_HITTYPE_H

namespace lcvn {

  typedef enum HType {
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
  } HitType;

}

#endif // CVN_HITTYPE_H
