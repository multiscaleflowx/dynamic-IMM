#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <set>          // std::set
#include <vector>       // std::vector

#include "CFDSim.H"

namespace cfdsim {

  int CFDSim::countOccurrencesOfSubstring(const std::string& str, const std::string& sub) {
    int occurrences = 0;
    std::string::size_type pos = 0;
    while ((pos = str.find(sub, pos )) != std::string::npos) {
      ++occurrences;
      pos += sub.length();
    }
    return occurrences;
  }

  void CFDSim::configureRegions() {
    std::cout << "Entering configureRegions" << std::endl;
    std::ifstream cmdFile("cmd");
    std::string cmd;
    getline(cmdFile, cmd);

    int count = countOccurrencesOfSubstring(cmd, std::string("mui-lmp"));
    std::cout << "count = " << count << std::endl;

    maxIndex = 0;
    if(count > 0) {
      int n = (count - 1) / 2;
      maxIndex = n * (n + 2);
    }
    std::cout << "maxIndex = " << maxIndex << std::endl;

    for(int j = 0; j < count; j++) {
      int index = (maxIndex - count) + j + 1;
      Region r;
      r.interfaceName = std::to_string(index); // Simulation index
      r.sNorm = double(j) / count;
      std::cout << "r = " << r << std::endl;
      regions.emplace_back(r);
      interfaceNames.emplace_back(r.interfaceName);
    }
    std::cout << "Leaving configureRegions" << std::endl;
  }

  /*
  std::vector<std::string> CFDSim::collectInterfaceIndices(const std::string& str) {
    std::string inMD("in.MD");
    std::vector<std::string> ifStrings;
    std::string::size_type pos = 0;
    while ((pos = str.find(inMD, pos )) != std::string::npos) {
      std::string ifStr = xxx;
      ifStrings.emplace_back(ifStr);
    }
    return ifStrings;
  }
  */

  /*
    The CFD class must implement a method to replace that given below.
    The implementation in this file is provided only to allow testing of the rest of the code.

    The function must return a bool indicating whether new simulations are needed or existing
    ones are no longer needed.

    This function populates two vectors: add and remove.
    When the function returns these 'add' must contain the indices of any new
    simulations/interfaces that are required and 'remove' must contain the indices of
    any simulations/interfaces that are no longer required.
  */
  bool CFDSim::changesRequired(int t, int iter) {
    std::cout << "Entering changesRequired" << std::endl;

    for(Region r : regions) {
      remove.emplace_back(r);
    }

    int count = interfaceNames.size();
    if(count == 0) count = 1;

    int numberOfRegionsToCreate = count + 2;
    for(int i = 0; i < numberOfRegionsToCreate; i++) {
      int index = maxIndex + 1 + i;
      Region r;
      r.interfaceName = std::to_string(index); // Simulation index
      r.sNorm = double(i) / numberOfRegionsToCreate;
      std::cout << "The next simulation will contain r = " << r << std::endl;
      add.emplace_back(r);
    }

    std::cout << "Leaving changesRequired" << std::endl;
    return true;
  }
}
