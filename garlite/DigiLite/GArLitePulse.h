/** ****************************************************************************
 * @file GArLitePulse.h
 * @brief Definition of tracking plane pulses in GArLite
 * @author federico.battisti@pmb.ox.ac.uk
 * @see  GArLitePulse.cxx raw.h
 * Created on the 20/27/21
 *
 * ****************************************************************************/

#ifndef GAR_RAWDATA_GARLITEPULSE_H
#define GAR_RAWDATA_GARLITEPULSE_H

// C/C++ standard libraries
#include <stdint.h> // uint32_t
#include <cstdlib> // size_t
#include <vector>
#include "RawDataProducts/RawTypes.h"
#include "RtypesCore.h"  // ULong64_t

/// Raw data description and utilities
namespace gar {
  namespace raw {

    typedef long long int fChannel_t;

    class GArLitePulse {

    public:

        /// Default constructor: an empty raw digit with zeros put in for paraneters and an invalid channel
      GArLitePulse();

#ifndef __GCCXML__
    public:

      GArLitePulse(float fCharge, double fTime, float x, float y, float z, fChannel_t fChannel);

      //Copy constructor
      GArLitePulse(gar::raw::GArLitePulse const&) = default;

      /// Reference to the compressed ADC count vector
      float Charge() const;
      /// Timestmap
      double Time() const;
      /// X position
      float X() const;
      /// Y position
      float Y() const;
      /// Z position
      float Z() const;
      /// cellID
      fChannel_t Channel() const;


#endif // !__GCCXML__
    private:

      float fCharge; ///< Charge deposited on the Sipm
      double fTime; ///< time of the hit for the SiPm
      float fX; ///< x of the hit
      float fY; ///< y of the hit
      float fZ; ///< z of the hit
      fChannel_t fChannel;              ///< ID of the individual Sipm Channel

    }; // class GArLitePulse


  } // namespace raw
} // gar

//------------------------------------------------------------------------------
//--- inline implementation
//---
#ifndef __GCCXML__

inline float                                            gar::raw::GArLitePulse::Charge()        const {return fCharge;   }
inline double                                           gar::raw::GArLitePulse::Time()          const {return fTime;  }
inline float                                            gar::raw::GArLitePulse::X()             const {return fX;	 }
inline float                                            gar::raw::GArLitePulse::Y()             const {return fY;	 }
inline float                                            gar::raw::GArLitePulse::Z()             const {return fZ;	 }
inline gar::raw::fChannel_t                             gar::raw::GArLitePulse::Channel()       const {return fChannel;}

#endif // !__GCCXML__

#endif // GAR_RAWDATA_GARLITEPULSE_H

////////////////////////////////////////////////////////////////////////
