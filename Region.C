#include <istream>
#include <ostream>

#include "Region.H"

namespace cfdsim {

  bool operator==(const Region& lhs, const Region& rhs) {
    return lhs.interfaceName == rhs.interfaceName && \
    lhs.sNorm == rhs.sNorm;
  }
  
  bool operator <(const Region& r1, const Region& r2) {
    return r1.sNorm < r2.sNorm;
  }

  std::ostream& operator<<(std::ostream& os, const Region& r) {
    os << "{ " << r.interfaceName << " " << r.sNorm << " }";
    return os;
  }

  std::istream& operator>>(std::istream& is, Region& r) {
    std::string s;
    //double h;
    bool error = false;
    is >> s;
    if(s != "{") {
      error = true;
    }
    else {
      is >> s;
      r.interfaceName = s;
      is >> s;
      r.sNorm = stod(s);
      is >> s;
      if(s != "}") {
	error = true;
      }
    }
      
    if(error)
      is.setstate(std::ios::failbit);
    return is;
  }

  double Region::normalisedStreamwiseCoord() {
    return sNorm;
  }

}
