#ifndef CFDSIM_H
#define CFDSIM_H

#include <mui.h>
#include <vector>
#include <string>

namespace cfdsim {
  struct Region {
    std::string interfaceName;
    std::string ylo;
    std::string yhi;
  };

  class CFDSim {
    int iters;
    int startTime;
    std::set<std::string> inVars;
    std::set<std::string> outVars;
    std::set<Region> firstTimeRegions;
    std::vector<Region> regions;
    std::vector<Region> updatedRegions;
    std::vector<std::string> interfaceNames;
    //std::vector<std::string> updatedInterfaceNames;
    std::string cfdFileName;
    int maxIndex;
    std::vector<Region> add;
    std::vector<Region> remove;

  public:
    CFDSim(char*);

    void readConfigFile();

    void writeConfigFile(int);

    void writeCmdFile(std::vector<int>, std::string);

    void createInitialMDFile(Region, std::string);

    void createRestartedMDFile(Region, std::string);

    bool changeMDs(int t);

    int countNodes(MPI_Comm);

    bool calculateNodeDistribution(int, std::vector<int>&);

    void exchangeData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>&,
		      int);

    void run();

    /*
      This method must be reimplemented.
    */
    bool changesRequired();

    /*
      This method must be reimplemented.
    */
    void calculateOutputs(int,
			  std::map<std::string, std::map<std::string, double>>&);
  };
}

#endif
