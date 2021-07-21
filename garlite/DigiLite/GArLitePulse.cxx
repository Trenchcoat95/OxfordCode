/** ****************************************************************************
* @file GArLitePulse.h
* @brief Definition of tracking plane pulses in GArLite
* @author federico.battisti@pmb.ox.ac.uk
* @see  GArLitePulse.cxx raw.h
* Created on the 20/27/21
*
* ****************************************************************************/

#include <limits>

#include "RawDataProducts/GArLitePulse.h"

// C/C++ standard libraries

namespace gar {
  namespace raw{

    //----------------------------------------------------------------------
    GArLitePulse::GArLitePulse()
    : fCharge(0.),
    fTime(0.),
    fX(0.),
    fY(0.),
    fZ(0.),
    fChannel(0)
    {}

    //----------------------------------------------------------------------
    GArLitePulse::GArLitePulse(float fCharge, double fTime, float x, float y, float z, fChannel_t fChannel)
    : fCharge(fCharge),
    fTime(fTime),
    fX(x),
    fY(y),
    fZ(z),
    fChannel(fChannel)
    {}

    } // namespace raw
  } // gar
    ////////////////////////////////////////////////////////////////////////
