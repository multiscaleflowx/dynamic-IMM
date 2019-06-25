#ifndef REGION_H
#define REGION_H

#include <string>

namespace cfdsim {

  class Region {
  public:
    std::string interfaceName;
    double sNorm; // Normalised streamwise coordinate: 0 <= sNorm < 1.

    double normalisedStreamwiseCoord();
  };

  bool operator==(const Region&, const Region&);
  bool operator<(const Region&, const Region&);
  std::ostream& operator<<(std::ostream&, const Region&);
  std::istream& operator>>(std::istream&, Region&);

}

#endif
