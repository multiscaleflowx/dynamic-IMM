#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <set>          // std::set
#include <vector>       // std::vector

#include "CFDSim.H"
#include "common.H"

namespace cfdsim {

  std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
      throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
      result += buffer.data();
    }
    return result;
  }

  CFDSim::CFDSim(const char* exe,
		 const char* name,
		 const char* openfoamCase,
		 const IOdictionary& immDict): exeName(exe),cfdFileName(name), caseName(openfoamCase) {
    configureRegions();

    double length;
    if(!immDict.readIfPresent("length", length)) {
      haltMPMD("the keyword 'length' is not present in immDict.");
    }
    double hEnd;
    if(!immDict.readIfPresent("hEnd", hEnd)) {
      haltMPMD("the keyword 'hEnd' is not present in immDict.");
    }
    double hNeck;
    if(!immDict.readIfPresent("hNeck", hNeck)) {
      haltMPMD("the keyword 'hNeck' is not present in immDict.");
    }
    channel = Channel(length, hEnd, hNeck);

    if(!immDict.readIfPresent("rhoN", rhoN_)) {
      haltMPMD("the keyword 'rhoN' is not present in immDict.");
    }

    if(!immDict.readIfPresent("nSolverCalls", nIter_)) {
      haltMPMD("the keyword 'nSolverCalls' is not present in immDict.");
    }
    if(!immDict.readIfPresent("numberOfSamplesToAverageOver", nSamples)) {
      haltMPMD("the keyword 'numberOfSamplesToAverageOver' is not present in immDict.");
    }

    if(!immDict.readIfPresent("molMass", molMass)) {
      haltMPMD("the keyword 'molMass' is not present in immDict.");
    }

    if(!immDict.readIfPresent("lengthOfRegion", lengthOfRegion)) {
      haltMPMD("the keyword 'lengthOfRegion' is not present in immDict.");
    }
    if(!immDict.readIfPresent("widthOfRegion", widthOfRegion)) {
      haltMPMD("the keyword 'widthOfRegion' is not present in immDict.");
    }

    if(!immDict.readIfPresent("microTimeStep",  time_dt)) {
      haltMPMD("the keyword 'microTimeStep' is not present in immDict.");
    }
    if(!immDict.readIfPresent("numberOfEquilibrationSteps", nEquilibration)) {
      haltMPMD("the keyword 'numberOfEquilibrationSteps' is not present in immDict.");
    }
    if(!immDict.readIfPresent("numberOfStepsBetweenSamples", nStepsBetweenSamples)) {
      haltMPMD("the keyword 'numberOfStepsBetweenSamples' is not present in immDict.");
    }

    // The number of time steps used to compute the average values.
    nMeasurement = nSamples * nStepsBetweenSamples;

    wordList pushVars;
    if(!immDict.readIfPresent("push", pushVars)) {
      haltMPMD("the keyword 'push' is not present in immDict.");
    }
    forAll(pushVars, v) {
      word vn = pushVars[v];
      std::ostringstream out;
      out << vn;
      std::string outVar = out.str();
      outVars.insert(outVar);

      bool isSubDict = immDict.isDict(vn);
      if(!isSubDict) {
	std::string msg("the keyword '");
	msg.append(vn);
	msg.append("' does not denote a subdictionary of immDict");
	haltMPMD(msg.c_str());
      }
      const dictionary& subDict = immDict.subDict(vn);

      double conversionFactor;
      if(!subDict.readIfPresent("conversionFactor", conversionFactor)) {
	std::string msg("the keyword 'conversionFactor' is not present in the subdictionary ");
	msg.append(vn);
	msg.append(".");
	haltMPMD(msg.c_str());
      }
      outputConversionFactors[outVar] = conversionFactor;

      double initialValue;
      if(!subDict.readIfPresent("initialValue", initialValue)) {
	std::string msg("the keyword 'initialValue' is not present in the subdictionary ");
	msg.append(vn);
	msg.append(".");
	haltMPMD(msg.c_str());
      }
      initialValues[outVar] = initialValue;
      std::cout << "pushVar: " << outVar << ", conversionFactor = " << conversionFactor << ", initialValue = " << initialValue << std::endl;
    }

    wordList fetchVars(immDict.lookup("fetch"));
    forAll(fetchVars, v) {
      word vn = fetchVars[v];
      std::ostringstream in;
      in << vn;
      std::string inVar = in.str();
      inVars.insert(inVar);

      bool isSubDict = immDict.isDict(vn);
      if(!isSubDict) {
	std::string msg("the keyword '");
	msg.append(vn);
	msg.append("' does not denote a subdictionary of immDict");
	haltMPMD(msg.c_str());
      }
      const dictionary& subDict = immDict.subDict(vn);

      double conversionFactor;
      if(!subDict.readIfPresent("conversionFactor", conversionFactor)) {
	std::string msg("the keyword 'conversionFactor' is not present in the subdictionary ");
	msg.append(vn);
	msg.append(".");
	haltMPMD(msg.c_str());
      }
      inputConversionFactors[inVar] = conversionFactor;

      double acceptableError;
      if(!subDict.readIfPresent("acceptableError", acceptableError)) {
	std::string msg("the keyword 'acceptableError' is not present in the subdictionary ");
	msg.append(vn);
	msg.append(".");
	haltMPMD(msg.c_str());
      }
      acceptableErrors[inVar] = acceptableError;

      double tolerance;
      if(!subDict.readIfPresent("tolerance", tolerance)) {
	std::string msg("the keyword 'tolerance' is not present in the subdictionary ");
	msg.append(vn);
	msg.append(".");
	haltMPMD(msg.c_str());
      }
      tolerances[inVar] = tolerance;
      std::cout << "fetchVar: " << inVar << ", conversionFactor = " << conversionFactor << ", acceptableError = " << acceptableError << ", tolerance = " << tolerance << std::endl;
    }
  }

  void CFDSim::writeCmdFile(std::vector<int> nodeDistribution) {
    std::cout << "Entering writeCmdFile" << std::endl;
    std::ofstream outfile("cmd");
    outfile << "aprun -n 1 " << exeName << " -case " << caseName;
    int i = 0;
    for(auto r : updatedRegions) {
      int p = nodeDistribution[i] * 24;
      outfile << " :  -n " << p << " mui-lmp -in in.MD" << r.interfaceName;
      i++;
    }
    outfile << " > output" << nodeDistribution.size() << std::endl;
    outfile.close();
    std::cout << "Leaving writeCmdFile" << std::endl;
  }

  void CFDSim::createInitialMDFile(Region r,
				   std::set<std::string>& outVars,
				   std::string niterM1,
				   std::string nequib,
				   std::string nevery,
				   std::string nsteps) {
    std::cout << "Entering createInitialMDFile" << std::endl;
    std::ifstream infs;
    if(r.sNorm == 0.0) {
      infs = std::ifstream("templates/initial_end_MD_template");
      if(!infs.good()) {
	haltMPMD("the file 'templates/initial_end_MD_template' does not exist.");
      }
    }
    else {
      infs = std::ifstream("templates/initial_MD_template");
      if(!infs.good()) {
	haltMPMD("the file 'templates/initial_MD_template' does not exist.");
      }
    }
    std::ofstream outfs("in.MD" + r.interfaceName);
    std::string line;
    double h = channel.heightAt(r.sNorm);
    int32_t N = round(rhoN_*lengthOfRegion*widthOfRegion*h/molMass);
    h = h * 1.0e10; // Convert from metres to Angstroms.
    double length = lengthOfRegion * 1.0e10;
    double width = widthOfRegion * 1.0e10;
    std::cout << "Initial MD file: " << r << ", h = " << h << ", N = " << N << std::endl;
    std::ostringstream h_strs;
    h_strs << h;
    std::string h_str = h_strs.str();
    std::ostringstream N_strs;
    N_strs << N;
    std::string N_str = N_strs.str();
    std::ostringstream l_strs;
    l_strs <<length;
    std::string l_str = l_strs.str();
    std::ostringstream w_strs;
    w_strs << width;
    std::string w_str = w_strs.str();
    std::cout << "Initial MD file: " << r << ", h = " << h << ", N = " << N << std::endl;
    std::string atSymb("@");
    while(getline(infs, line)) {
      size_t index;
      for(std::string var : outVars) {
	double val =  outputConversionFactors[var] * initialValues[var];
	std::ostringstream val_strs;
	val_strs << val;
	std::string val_str = val_strs.str();
	std::string v = atSymb + var;
	std::cout << "Replacing " << v << " by " << val_str << std::endl;
	while((index = line.find(v)) != std::string::npos) {
	  line.replace(index, v.size(), val_str);
	}
      }
      while((index = line.find("@h")) != std::string::npos) {
	line.replace(index, 2, h_str);
      }
      while((index = line.find("@N")) != std::string::npos) {
	line.replace(index, 2, N_str);
      }
      while((index = line.find("@niterM1")) != std::string::npos) {
	line.replace(index, 8, niterM1);
      }
      while((index = line.find("@id")) != std::string::npos) {
	line.replace(index, 3, r.interfaceName);
      }
      while((index = line.find("@nequib")) != std::string::npos) {
	line.replace(index, 7, nequib);
      }
      while((index = line.find("@nevery")) != std::string::npos) {
	line.replace(index, 7, nevery);
      }
      while((index = line.find("@nsteps")) != std::string::npos) {
	line.replace(index, 7, nsteps);
      }
      while((index = line.find("@l")) != std::string::npos) {
        line.replace(index, 2, l_str);
      }
       while((index = line.find("@w")) != std::string::npos) {
        line.replace(index, 2, w_str);
      }
      outfs << line << std::endl;
    }
    std::cout << "Leaving createInitialMDFile" << std::endl;
  }

  void CFDSim::createRestartedMDFile(Region r,
				     std::set<std::string>& outVars,
				     std::string niterM1,
				     std::string nequib,
				     std::string nevery,
				     std::string nsteps) {
    std::cout << "Entering createRestartedMDFile" << std::endl;
    haltMPMD("this method has not been implemented.");
  }

  bool CFDSim::changeMDs(int t, int iter) {
    std::cout << "Entering changeMDs" << std::endl;
    bool change = false;

    updatedRegions.clear(); // This is necessary as changeMDs may be called more than once.
    // Make a copy of interfaceNames.
    for(auto r : regions) {
      updatedRegions.push_back(r);
    }

    if(changesRequired(t, iter)) {

      // Remove the interface names that are no longer needed.
      for(auto r : remove) {
	std::cout << "changeMDs: removing " << r << std::endl;
	auto it = std::find(updatedRegions.begin(), updatedRegions.end(), r);
	updatedRegions.erase(it);
      }

      // Those MDs that existed before this run and are still needed must have
      // restart files created for them.
      for(auto r : updatedRegions) {
	createRestartedMDFile(r,
			      outVars,
			      std::to_string(nIter_ - 1),
			      std::to_string(nEquilibration),
			      std::to_string(nStepsBetweenSamples),
			      std::to_string(nMeasurement));
      }

      // Add the new interface names and create initial restart files.
      for(auto r : add) {
	updatedRegions.emplace_back(r);
	createInitialMDFile(r,
			    outVars,
			    std::to_string(nIter_ - 1),
			    std::to_string(nEquilibration),
			    std::to_string(nStepsBetweenSamples),
			    std::to_string(nMeasurement));
      }

      change = true;
    }

    std::cout << "Leaving changeMDs: change = " << change << std::endl;
    return change;
  }

  int CFDSim::countNodes() {
    char* pVar;
    pVar = getenv ("PBS_JOBID");
    std::string command("qstat -a | grep ");
    command.append(pVar);
    command.append(" | awk '{print $6}'");
    const char* cmd = command.c_str();
    std::string result = exec(cmd);

    return std::stoi(result);
  }

  bool CFDSim::calculateNodeDistribution(int nNodes,
					 std::vector<int>& nodeDistribution) {
    std::cout << "Entering calculateNodeDistribution" << std::endl;
    int nMDs = nodeDistribution.size();
    std::cout << "calculateNodeDistribution: nMDs = " << nMDs << ", nMDNodes = " << nNodes << std::endl;
    if(nMDs > nNodes) { // Arbitrarily assign 1 node to each MD.
      for(auto& x : nodeDistribution) {
	x = 1;
      }
      std::cout << "Leaving calculateNodeDistribution" << std::endl;
      return false;
    }

    int baseNodes = nNodes / nMDs;
    int excess = nNodes - baseNodes*nMDs;
    std::cout << "calculateNodeDistribution: baseNodes = " << baseNodes << std::endl;
    std::cout << "calculateNodeDistribution: excess = " << excess << std::endl;
    for(auto& x : nodeDistribution) {
      if(excess) {
	x = baseNodes + 1;
	excess--;
      }
      else {
	x = baseNodes;
      }
      std::cout << "calculateNodeDistribution: x = " << x << std::endl;
    }

    std::cout << "Leaving calculateNodeDistribution" << std::endl;
    return true;
  }

  //TODO: order these correctly.
  void CFDSim::calculateAverages() {
    std::cout << "Entering calculateAverages" << std::endl;
    for(auto& ifn : interfaceNames) {
      for(std::string var : inVars) {
	double sum = 0.0;
	for(auto& val : inVarValuesVec[ifn][var]) {
	  sum += val;
	}
	std::cout << "ifn = " << ifn << ", var = " << var << ", average = " << (sum / nSamples) << std::endl;
	averages[ifn][var].push_back(sum / nSamples);
      }
    }
    std::cout << "Leaving calculateAverages" << std::endl;
  }

  void CFDSim::clearAccumulators() {
    std::cout << "Entering clearAccumulators" << std::endl;
    for(auto& ifn : interfaceNames) {
      for(std::string var : inVars) {
	inVarValuesVec[ifn][var].clear();
	std::cout << "cleared ifn = " << ifn << ", var = " << var << std::endl;
      }
    }
    std::cout << "Leaving clearAccumulators" << std::endl;
  }

  void CFDSim::receiveData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>& interfaces,
			   int step) {
    std::cout << "Entering receiveData" << std::endl;

    double t = step*time_dt;

    std::cout << "receiveData: step = " << step << ", t = " << t << std::endl;

    int i = 0;
    for(auto& interface : interfaces) {
      std::string ifn = reorderedInterfaceNames[i];

      interface->commit(t);  //TODO: understand this.

      interface->barrier(t);
      std::cout << "After barrier." << std::endl;
      for(std::string var : inVars) {
	double val = interface->fetch<double>(var);
	double inputValue = inputConversionFactors[var] * val;
	inVarValuesVec[ifn][var].push_back(inputValue);
	std::cout << "CFD in t = " << t << " " << ifn << ' ' << var << ' ' << inputValue << std::endl;
      }

      i++;
    }
    std::cout << "Leaving receiveData" << std::endl;
  }

  void CFDSim::sendData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>& interfaces,
			int step) {
    std::cout << "Entering sendData" << std::endl;

    double t = step*time_dt;

    std::cout << "sendData: step = " << step << ", t = " << t << std::endl;

    int i = 0;
    for(auto& interface : interfaces) {
      std::string ifn = reorderedInterfaceNames[i];

      for(std::string var : outVars) {
	double outputValue = outputConversionFactors[var] * outVarValues[ifn][var];
	std::cout << "outVarValues[" << ifn << "][" << var << "]" << outVarValues[ifn][var] << std::endl;
        std::cout << "CFD out t = " << t << " " << ifn << ' ' << var << ' ' << outputValue << std::endl;
        interface->push(var, outputValue);
      }
      interface->commit(t);

      i++;
    }
    std::cout << "Leaving sendData" << std::endl;
  }

  /*
    The following method is more complicated than it needs to be for the IMM application
    because it is a (not fully) simplified version of code that was developed with the
    understanding that it would interact with a dynamic CFD simulation, whereas it is
    now being used for a steady flow.
   */
  void CFDSim::run() {
    std::cout << "Entering run" << std::endl;

    std::cout << "Interface names:";
    for(std::string ifn : interfaceNames) {
      std::cout << " " << ifn;
    }
    std::cout << std::endl;

    std::cout << "CFD creating interfaces." << std::endl;
    auto interfaces = mui::create_uniface<mui::config_3d>(std::string("CFD"), interfaceNames);
    std::cout << "CFD interfaces created." << std::endl;

    std::hash<std::string> str_hash; // Hash function for strings.
    std::vector<int> hashes;
    std::set<int> hashSet;
    for(std::string ifn :interfaceNames) {
      int hash = (int)str_hash(ifn); // MUI has an implicit conversion from size_t to int.
      hashes.emplace_back(hash); // These are in the order in which the interface names are supplied.
      hashSet.insert(hash); // These are in lexical order of the hashes.
    }
    if(hashes.size() != hashSet.size()) {
      haltMPMD("not all of the hashed names are distinct.");
    }
    // It is easier to work with a sorted vector so create one from the sorted hash set.
    std::vector<int> sortedHashes(hashSet.begin(), hashSet.end());

    std::cout << "hashes:";
    for(int h : hashes) {
      printf( " %08X", h);
    }
    std::cout << std::endl;

    std::set<int>::iterator it = hashSet.begin();
    std::cout << "hashSet:";
    while (it != hashSet.end()) {
      printf( " %08X", (*it));
      it++;
    }
    std::cout << std::endl;

    std::cout << "sorted hashes:";
    for(int h : sortedHashes) {
      printf( " %08X", h);
    }
    std::cout << std::endl;

    // Sort the interface names to match the order in which the unifaces (interfaces) are created.
    for(int h : sortedHashes) {
      ptrdiff_t pos = std::find(hashes.begin(), hashes.end(), h) - hashes.begin();
      reorderedInterfaceNames.push_back(interfaceNames[pos]);
      std::cout << "pos = " << pos << ", ifn = " << interfaceNames[pos] << std::endl;
      printf( "hash = %08X\n", h);
    }

    std::cout << "sorted names:";
    for(std::string s : reorderedInterfaceNames) {
      std::cout << " " << s;
    }
    std::cout << std::endl;

    int nNodes;
    nNodes = countNodes();
    std::cout << "A: nNodes = " << nNodes << std::endl;
    int nMDNodes = nNodes-1;

    int NEVERY = nStepsBetweenSamples;

    bool halt = false;
    bool terminate = false;
    int64_t timestep = 1;
    int numberOfMDs = interfaceNames.size();
    std::map<std::string, double> previousEstimates;
    bool previousEstimatesExist = getPreviousEstimates(numberOfMDs, previousEstimates);

    std::string solver_output("solver_output");
    solver_output.append(std::to_string(numberOfMDs));
    ofs.open(solver_output); // Where to put a copy of the output from the solver.

    initialise(numberOfMDs);
    //TODO: replace this!
    //int i = 0;
    std::cout << "Set outputs before run." << std::endl;
    for(std::string ifn : interfaceNames) {
      for(std::string v : outVars) {
	outVarValues[ifn][v] = initialValues[v];
      }
    }
    std::cout << "Outputs set before run." << std::endl;

    std::cout << "run: nEquilibration = " << nEquilibration << std::endl;
    std::cout << "run: nMeasurement = " << nMeasurement << std::endl;

    bool hasConverged = (numberOfMDs == 0); // Do not want to test for convergence if the are no MDs.
    bool hasConvergedOverall = false;
    for(int iter_ = 0; iter_ < nIter_; iter_++) {
      int32_t sampleCount = 0;

      //sendData(interfaces, timestep);
      for(int i = 0; (i < nEquilibration) && (!terminate) && (!hasConverged) && (!hasConvergedOverall); i++) {
	//std::cout << "run:equib: timestep = " << timestep << std::endl;

	if((timestep % NEVERY) == 0) {
	  std::cout << "run:equib:NEVERY: timestep = " << timestep << std::endl;

	  receiveData(interfaces, timestep);

          sampleCount++;
          if(sampleCount == nSamples) {
	    calculateAverages();
	    setInputs();
	    clearAccumulators();
	    sampleCount = 0;
	  }

	  // The initial force is set in the LAMMPS script.
	  //sendData(interfaces, timestep+NEVERY);
	  sendData(interfaces, timestep);

          double val = double(false);
	  double max_val; // To keep MPI_Allreduce happy. It will be set to the value of val.
	  MPI_Allreduce(&val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  std::cout << "MPI_Allreduce: terminate at timestep " << timestep << " = " << max_val << ", halt = " << halt << std::endl;
	}
	timestep++;
      }

      std::cout << "run: may start measurement: iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
      std::map<std::string, double> means;
      for(int i = 0; (i < nMeasurement) && (!terminate); i++) {
	if(i == 0) {
	  std::cout << "run:start measurement: iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
	}
	if((timestep % NEVERY) == 0) {
	  std::cout << "run:measurement:NEVERY: timestep = " << timestep << ", hasConverged = " << hasConverged << ", iter_ = " << iter_ << std::endl;

	  receiveData(interfaces, timestep);

          sampleCount++;
	  //std::cout << "SAMPLE_COUNT = " << sampleCount << ", i = " << i << ", iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
          if(sampleCount == nSamples) {
	    //std::cout << "AVERAGE" << std::endl;
	    calculateAverages();
	    setInputs();
	    clearAccumulators();
	    sampleCount = 0;

	    std::cout << "SOLVE, iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
	    solve(iter_, numberOfMDs);
	    hasConverged = hasConverged || (converged(means) && (iter_ > 0));
	    hasConvergedOverall = (iter_ > 0) && convergedOverall(previousEstimatesExist, previousEstimates, means);

	    calculateOutputs(timestep);
          }

          sendData(interfaces, timestep);

	  if((i == (nMeasurement-1)) && (hasConverged || iter_ == (nIter_ - 1))) {
	    if(changeMDs(timestep, iter_)) {
	      int newNumberOfMDs = updatedRegions.size();
	      std::cout << "timestep = " << timestep << ": newNumberOfMDs = " << newNumberOfMDs << std::endl;
	      std::cout << "timestep = " << timestep << ": maxIndex = " << maxIndex << std::endl;
	      std::vector<int> nodeDistribution(newNumberOfMDs);
	      if((newNumberOfMDs == 0 ) && (numberOfMDs > 0)) {
		std::cout << "No MDs left, so must halt." << std::endl;
		halt = true;
	      }
	      else if(calculateNodeDistribution(nMDNodes, nodeDistribution)) {
		std::cout << "Change of configuration required for MDs, so must terminate." << std::endl;
	      }
	      else { // There is no possible node distribution.
		std::cout << "Too few nodes for new MDs, so must terminate." << std::endl;
		halt = true;
	      }

	      writeCmdFile(nodeDistribution);

	      if(halt) {
		std::system("mv cmd halt");
	      }
	      terminate = true;
	    }
	    if(hasConverged) {
	      terminate = true;
	    }
	    /*
	    if(numberOfMDs == 5) { // TODO: FOR TEST PUPOSES ONLY - REMOVE EVENTUALLY.
	      hasConvergedOverall = true;
	      ofs << "Asserting overall convergence for numberOfMDs = 5." << std::endl;
	    }
	    */
	    if(hasConvergedOverall) {
	      halt = true;
	      std::ifstream cmdFile("cmd");
	      if(cmdFile.good()) { // Call move only if the file exists.
		std::system("mv cmd halt");
	      }
	    }
	  }

	  double val = double(terminate);
	  double max_val; // To keep MPI_Allreduce happy. It will be set to the value of val.
	  MPI_Allreduce(&val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  //std::cout << "MPI_Allreduce: terminate at timestep " << timestep << " = " << max_val << ", halt = " << halt << std::endl;
	}
        timestep++;
	//if(terminate) {
	//  break;
	//}
      }
    }
    std::cout << "Leaving run" << std::endl;
  }

  bool CFDSim::converged(std::map<std::string, double>& means) {
    int nMicro_ = interfaceNames.size();

    if(nMicro_ == 0) {
      return true;
    }

    std::ofstream estimatesFileStream("estimates");
    bool hasConverged;
    for(std::string v : inVars) {
      double sum = 0.0;
      for(std::string ifn : interfaceNames) {
	sum += averages[ifn][std::string(v)].back();
      }
      double mean = sum / nMicro_;
      means[v] = mean;

      double variance = 0.0;
      for(std::string ifn : interfaceNames) {
	double valueForRegion = averages[ifn][std::string(v)].back();
	variance += pow(valueForRegion - mean, 2);
      }

      double standardError = std::sqrt(variance / nMicro_);

      double normalisedError = standardError / std::abs(mean);

      ofs << "For " << v << " the mean is " << mean << " and the normalised error is " << normalisedError << std::endl;

      estimatesFileStream << v << ' ' << mean << std::endl;

      hasConverged = hasConverged && normalisedError < acceptableErrors[v];
    }
    estimatesFileStream.close();

    ofs << "hasConverged = " << hasConverged << std::endl;
    return hasConverged;
  }

  bool CFDSim::convergedOverall(bool previousEstimatesExist,
				std::map<std::string, double>& previousEstimates,
				std::map<std::string, double>& means) {
    if(!previousEstimatesExist) {
      return false;
    }

    bool hasConvergedOverall;
    for(std::string v : inVars) {
      double previousEstimate = previousEstimates[v];
      double mean = means[v];
      std::cout << " For " << v << " previousEstimate = " << previousEstimate << ", mean = " << mean << std::endl;

      double normalisedChange = std::abs(previousEstimate - mean) / mean;

      ofs << "For " << v << " the overall normalised change is " << normalisedChange << std::endl;

      hasConvergedOverall = hasConvergedOverall && normalisedChange < tolerances[v];
    }

    ofs << "hasConvergedOverall = " << hasConvergedOverall << std::endl;
    return hasConvergedOverall;
  }

  bool CFDSim::getPreviousEstimates(int numberOfMDs,
				    std::map<std::string, double>& estimates) {
    std::cout << "Entering getPreviousEstimates." << std::endl;
    bool previousEstimatesExist = false;
    if(numberOfMDs != 0) {
      std::ifstream estimatesFileStream("estimates");
      previousEstimatesExist = estimatesFileStream.good();
      if(previousEstimatesExist) { // File exists.
	std::cout << "The file 'estimates' exists." << std::endl;
	std::string line;
	std::string name;
	double value;
	while(getline(estimatesFileStream, line)) {
	  std::istringstream estimatesStream(line);
	  estimatesStream >> name;
	  estimatesStream >> value;
	  std::cout << "name = " << name << ", value = " << value << std::endl;
	  estimates[name] = value;
	}
	estimatesFileStream.close();
      }
    }
    return previousEstimatesExist;
  }

}
