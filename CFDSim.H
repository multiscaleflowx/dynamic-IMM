#ifndef CFDSIM_H
#define CFDSIM_H

#include "mui.h"
#include <vector>
#include <string>

#include "IOdictionary.H"
#include "polyMesh.H"
#include "Switch.H"
#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "SVD.H"
#include "SortableList.H"
#include "simpleMatrix.H"
#include "Time.H"

#include "Random.H"
#include "Pair.H"

#include "Region.H"
#include "Channel.H"

using namespace Foam;

namespace cfdsim {
  class CFDSim {
    double forceConversionFactor = 1.44e10;
    double massFlowRateConversionFactor = 1.661e-12;

    std::set<std::string> inVars;
    std::set<std::string> outVars;
    std::map<std::string, std::map<std::string, double>> outVarValues;
    std::map<std::string, std::map<std::string, std::vector<double>>> inVarValuesVec;
    std::map<std::string, std::map<std::string, std::vector<double>>> averages;
    std::vector<Region> regions;
    std::vector<Region> updatedRegions;
    std::vector<std::string> interfaceNames;
    std::vector<string> reorderedInterfaceNames;
    const std::string exeName;
    const std::string cfdFileName;
    const std::string caseName;
    int maxIndex;
    std::vector<Region> add;
    std::vector<Region> remove;
    Channel channel;
    std::ofstream ofs;

    label nIter_; // Number of iterations.
    scalar L_; // Length of macro channel.
    scalar F_; // Body force in Newtons

    double time_dt;
    double acceptableError;
    double tolerance;
    int64_t nStepsPerIter; // Total number of micro element steps per iteration.
    int32_t nEquilibration; // Number of equilibration steps.
    int32_t nStepsBetweenSamples;
    int32_t nSamples; // Number of samples to use when averaging.
    int64_t nMeasurement; // Number of measurement steps.
    double molMass;
    double lengthOfRegion;
    double widthOfRegion;

    scalarField s_; // collocation points, co-ordinates of micro elements
    scalarField mDot_; // mass flow rates mDot_i
    scalarField k_; // flow resistances k_i
    scalarField f_; // micro forces f_i
    scalarField f_old_; // previous micro forces f_old_i
    List< List< List<scalar> > > kIs_; // keep track of older flow resistance information 
    scalar rhoN_; // reference number density //*** check 
    scalarField phiCoeffs_; // Fourier phi coefficients
    scalarField phi_; // momentum flux corrector terms Phi_i

  public:
    CFDSim(const char*, const char*, const char*, const IOdictionary&);

    void initialise(int);

    int countOccurrencesOfSubstring(const std::string&, const std::string&);

    void configureRegions();

    //void readConfigFile();

    //void writeConfigFile();

    void writeCmdFile(std::vector<int>);

    void createMDForceFile(Region, std::string, std::string);
    
    void createInitialMDFile(Region,
			     std::string,
			     std::string,
			     std::string,
			     std::string,
			     std::string);

    void createRestartedMDFile(Region,
			       std::string,
			       std::string,
			       std::string,
			       std::string,
			       std::string);

    bool changeMDs(int, int);

    int countNodes(MPI_Comm);

    int countNodes();

    bool calculateNodeDistribution(int, std::vector<int>&);

    void solve(label, label);

    //void calculateBodyForces(label, label, std::map<std::string, std::map<std::string, double>>);

    void receiveData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>&,
		     int);

    void sendData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>&,
		  int);

    void calculateAverages();

    void clearAccumulators();

    void run();

    bool changesRequired(int, int);

    bool converged(int, double&);

    bool convergedOverall(bool, double, double);
  };
}

#endif
