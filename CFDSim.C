#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <set>          // std::set
#include <vector>       // std::vector

#include "CFDSim.H"

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

namespace cfdsim {

  CFDSim::CFDSim(const char* exe,
		 const char* name,
		 const char* openfoamCase,
		 const IOdictionary& macroDict,
		 const IOdictionary& microDict): exeName(exe),cfdFileName(name), caseName(openfoamCase) {
    readConfigFile();

    const double length = readScalar(macroDict.lookup("length"));
    const double hEnd = readScalar(macroDict.lookup("hEnd"));
    const double hNeck = readScalar(macroDict.lookup("hNeck"));
    channel = Channel(length, hEnd, hNeck);

    std::string solver_output("solver_output");
    solver_output.append(std::to_string(startTime));
    ofs.open(solver_output); // Where to put a copy of the output from the solver.

    F_ = readScalar(macroDict.lookup("F"));
    std::cout << "Read in: F_ = " << F_ << std::endl;
    rhoN_ = readScalar(macroDict.lookup("rhoN"));

    nIter_= readLabel(macroDict.lookup("nSolverCalls"));
    nSamples = readLabel(macroDict.lookup("numberOfSamplesToAverageOver"));
    acceptableError = readScalar(macroDict.lookup("acceptableError"));
    tolerance = readScalar(macroDict.lookup("tolerance"));

    molMass = readScalar(microDict.lookup("molMass"));

    lengthOfRegion = readScalar(microDict.lookup("lengthOfRegion"));
    widthOfRegion = readScalar(microDict.lookup("widthOfRegion"));

    time_dt = readScalar(microDict.lookup("microTimeStep"));
    nEquilibration = readLabel(microDict.lookup("numberOfEquilibrationSteps"));
    nStepsBetweenSamples = readLabel(microDict.lookup("numberOfStepsBetweenSamples"));

    // The number of time steps used to compute an average flow rate.
    nMeasurement = nSamples * nStepsBetweenSamples;
  }

  void CFDSim::initialise(int nMicro_) {
    Info << nl << "Initialising macro solver" << nl << endl;

    s_.setSize(nMicro_, 0.0);
    mDot_.setSize(nMicro_, 0.0);
    k_.setSize(nMicro_, 1.0);
    f_.setSize(nMicro_, F_); // F_ is the correct initial value for f_.
    f_old_.setSize(nMicro_, 0.0);
    phi_.setSize(nMicro_, 0.0);
    phiCoeffs_.setSize(nMicro_, 0.0);
    kIs_.setSize(nMicro_);

    std::set<double> normalisedPositions;
    std::cout << "normalisedPositions:";
    for(Region r : regions) {
      normalisedPositions.insert(r.sNorm);
      std::cout << " " << r.sNorm;
    }
    std::cout << std::endl;

    // Compute the heights of the channel for the various regions.
    label i = 0;
    for(double s : normalisedPositions) {
      s_[i] = s*channel.length();
      std::cout << "s_[" << i << "] = " << s_[i] << std::endl;
      i++;
    }

    forAll(kIs_, i) {
      kIs_[i].setSize(nIter_);
      forAll(kIs_[i], j) {
	kIs_[i][j].setSize(3, 0.0);
      }
    }
    Info << nl << "Macro solver initialised" << nl << endl;
  }

  void CFDSim::readConfigFile() {
    std::cout << "Entering readConfigFile" << std::endl;
    std::ifstream infs(cfdFileName);
    std::string keyWord, line, s;

    infs >> keyWord;
    if(keyWord != "startAtIteration") {
      std::cout << "ERROR: '" << keyWord << "' found where 'startAtIteration' was expected." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    infs >> startTime;

    infs >> keyWord;
    if(keyWord != "push") {
      std::cout << "ERROR: '" << keyWord << "' found where 'push' was expected." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    getline(infs, line);
    std::istringstream str_stream2(line);
    while(str_stream2 >> s) {
      outVars.insert(s);
    }

    infs >> keyWord;
    if(keyWord != "fetch") {
      std::cout << "ERROR: '" << keyWord << "' found where 'fetch' was expected." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    getline(infs, line);
    std::istringstream str_stream3(line);
    while(str_stream3 >> s) {
      inVars.insert(s);
    }

    // Initially OpenFOAM will be run alone and no interfaces will be needed.
    // Once an MD simulation is required OpenFOAM will be restarted with interfaces
    // to the required MD simulations.
    infs >> keyWord;
    if(keyWord != "regions") {
      std::cout << "ERROR: '" << keyWord << "' found where 'regions' was expected." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    getline(infs, line);
    std::istringstream str_stream1(line);
    Region r;
    std::cout << "regions:";
    while(str_stream1 >> r) {
      regions.emplace_back(r);
      std::cout << " " << r;
      interfaceNames.emplace_back(r.interfaceName);
    }
    std::cout << std::endl;

    infs >> keyWord;
    if(keyWord != "maxIndex") {
      std::cout << "ERROR: '" << keyWord << "' found where 'maxIndex' was expected." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    infs >> maxIndex;
    if(interfaceNames.size() == 0) {
      if(maxIndex != 0) {
	std::cout << "ERROR: maxIndex '" << maxIndex << "' found where '0' was expected." << std::endl;
	// Make sure program terminates even in MPMD mode.
	MPI_Abort(MPI_COMM_WORLD, 999);
      }
    }

    infs.close();
    std::cout << "Leaving readConfigFile" << std::endl;
  }

  void CFDSim::writeConfigFile(int t) {
    std::cout << "Entering writeConfigFile" << std::endl;
    std::ifstream cfdInFS(cfdFileName);
    std::string newCfdFileName("new_"+cfdFileName);
    std::ofstream cfdOutFS(newCfdFileName);
    std::string cfdline;
    std::getline(cfdInFS, cfdline); // Discard and replace.
    cfdOutFS << "startAtIteration " << t+1 << std::endl;
    std::getline(cfdInFS, cfdline);
    cfdOutFS << cfdline <<std::endl; // Write push vars.
    std::getline(cfdInFS, cfdline);
    cfdOutFS << cfdline <<std::endl; // Write fetch vars.

    cfdOutFS << "regions";
    for(auto r : updatedRegions) {
      cfdOutFS << ' ' << r;
    }
    cfdOutFS << std::endl;

    cfdOutFS << "maxIndex " << maxIndex << std::endl;
    cfdOutFS.close();

    cfdInFS.close();

    std::system(("mv " + newCfdFileName + " " + cfdFileName).c_str());
    std::cout << "Leaving writeConfigFile" << std::endl;
  }

  void CFDSim::writeCmdFile(std::vector<int> nodeDistribution, std::string t) {
    std::cout << "Entering writeCmdFile" << std::endl;
    std::ofstream outfile("cmd");
    outfile << "aprun -n 1 " << exeName << " -case " << caseName;
    int i = 0;
    for(auto r : updatedRegions) {
      int p = nodeDistribution[i] * 24;
      outfile << " :  -n " << p << " mui-lmp -in in.MD" << r.interfaceName;
      i++;
    }
    outfile << " > output" << t << std::endl;
    outfile.close();
    std::cout << "Leaving writeCmdFile" << std::endl;
  }

  void CFDSim::createMDForceFile(Region r,
				 std::string f,
				 std::string a) {
    std::cout << "Entering createMDForceFile" << std::endl;
    std::ifstream infs("force_template");
    std::ofstream outfs("in.force" + r.interfaceName + "_" + a);
    std::string line;
    while(getline(infs, line)) {
      size_t index;
      while((index = line.find("@f")) != std::string::npos) {
	line.replace(index, 2, f);
      }
      while((index = line.find("@id")) != std::string::npos) {
	line.replace(index, 3, r.interfaceName);
      }
      outfs << line << std::endl;
    }
    outfs.close();
    std::cout << "Leaving createMDForceFile" << std::endl;
  }

  void CFDSim::createInitialMDFile(Region r,
				   std::string f,
				   std::string niterM1,
				   std::string nequib,
				   std::string nevery,
				   std::string nsteps,
				   std::string startTime) {
    std::cout << "Entering createInitialMDFile" << std::endl;
    std::ifstream infs("initial_MD_template");
    std::ofstream outfs("in.MD" + r.interfaceName);
    std::string line;
    double h = channel.heightAt(r.sNorm);
    int32_t N = round(rhoN_*lengthOfRegion*widthOfRegion*h/molMass);
    h = h * 1.0e10; // Convert from metres to Angstroms.
    double length = lengthOfRegion * 1.0e10;
    double width = widthOfRegion * 1.0e10;
    std::cout << "Initial MD file: " << r << ", h = " << h << ", N = " << N << ", f = " << f << std::endl;
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
    std::cout << "Initial MD file: " << r << ", h = " << h << ", N = " << N << ", f = " << f << std::endl;
    while(getline(infs, line)) {
      size_t index;
      while((index = line.find("@f")) != std::string::npos) {
	line.replace(index, 2, f);
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
      while((index = line.find("@t")) != std::string::npos) {
        line.replace(index, 2, startTime);
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
				     std::string f,
				     std::string niterM1,
				     std::string nequib,
				     std::string nevery,
				     std::string nsteps,
				     std::string startTime) {
    std::cout << "Entering createRestartedMDFile" << std::endl;
    std::ifstream infs("restarted_MD_template");
    std::ofstream outfs("in.MD" + r.interfaceName);
    std::string line;
    double h = channel.heightAt(r.sNorm);
    int32_t N = round(rhoN_*lengthOfRegion*widthOfRegion*h/molMass);
    h = h * 1.0e10; // Convert from metres to Angstroms.
    double length = lengthOfRegion * 1.0e10;
    double width = widthOfRegion * 1.0e10;
    std::cout << "Restarted MD file: " << r << ", h = " << h << ", N = " << N << ", f = " << f << std::endl;
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
    std::cout << "Restarted MD file: " << r << ", h = " << h << ", N = " << N << ", f = " << f << std::endl;
    while(getline(infs, line)) {
      size_t index;
      while((index = line.find("@f")) != std::string::npos) {
	line.replace(index, 2, f);
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
      while((index = line.find("@t")) != std::string::npos) {
        line.replace(index, 2, startTime);
      }
       while((index = line.find("@l")) != std::string::npos) {
        line.replace(index, 2, l_str);
      }
       while((index = line.find("@w")) != std::string::npos) {
        line.replace(index, 2, w_str);
      }
      outfs << line << std::endl;
    }
    std::cout << "Leaving createRestartedMDFile" << std::endl;
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
	double mdForce = forceConversionFactor * F_; // TODO: replace F_ by f_[i].
	std::ostringstream f_strs;
	f_strs << mdForce;
	std::string mdF_str = f_strs.str();

	createRestartedMDFile(r,
			      mdF_str,
			      std::to_string(nIter_ - 1),
			      std::to_string(nEquilibration),
			      std::to_string(nStepsBetweenSamples),
			      std::to_string(nMeasurement),
			      std::to_string(t));
      }

      // Add the new interface names and create initial restart files.
      for(auto r : add) {
	updatedRegions.emplace_back(r);

	double mdForce = forceConversionFactor * F_;
	std::cout << "Writing to MD file: F_ = " << F_ << std::endl;
	std::cout << "Writing to MD file: mdForce = " << mdForce << std::endl;
	std::ostringstream f_strs;
	f_strs << mdForce;
	std::string mdF_str = f_strs.str();

	createInitialMDFile(r,
			    mdF_str,
			    std::to_string(nIter_ - 1),
			    std::to_string(nEquilibration),
			    std::to_string(nStepsBetweenSamples),
			    std::to_string(nMeasurement),
			    std::to_string(t));

	//std::cout << "Creating force files for new simulations." << std::endl;
	//createMDForceFile(r, mdF_str, "0");	
      }

      writeConfigFile(t);

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

  int CFDSim::countNodes(MPI_Comm comm) {
    std::cout << "Entering countNodes" << std::endl;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int nameLength;
    MPI_Get_processor_name(processorName, &nameLength);

    int size;
    MPI_Comm_size(comm, &size);
    char processorNames[size*nameLength];
    int recvCounts[size];
    int displacements[size];
    for(int i = 0; i < size; i++) {
      recvCounts[i] = nameLength;
      displacements[i] = i*nameLength;
    }
    MPI_Allgatherv(processorName, nameLength, MPI_CHAR, processorNames,
		   recvCounts, displacements, MPI_CHAR, comm);

    std::set<std::string> names;
    for(int i = 0; i < size; i++) {
      char *name_p = processorNames + displacements[i];
      std::string name(name_p, nameLength);
      names.insert(name);
    }

    std::cout << "Leaving countNodes" << std::endl;
    return names.size();
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

  void CFDSim::solve(label iter_, label nMicro_) {
    scalar pi = constant::mathematical::pi;

    Info << nl << "iter:" << iter_
	 << " solving macro equation"
	 << endl;

    // 1. set k_i's from current MD simulation data

    if(iter_ == 0) {
      Info << nl << "applying initial flow resistance adjustments" << endl;

      for(label i = 0; i < nMicro_; i++) {
	ofs << "mDot_[" << i << "] = " << mDot_[i] << std::endl;
	if(mDot_[i] != 0.0) {
	  k_[i] = f_[i]/mDot_[i];

	  kIs_[i][iter_][0] = mDot_[i];
	  kIs_[i][iter_][1] = f_[i];
	  kIs_[i][iter_][2] = 0.1*mDot_[i]; //mDotError_[i] //*** need to pass flow rate error
	}
	ofs << "kIs_[" << i << "][" << iter_ << "][0] = " << kIs_[i][iter_][0] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][1] = " << kIs_[i][iter_][1] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][2] = " << kIs_[i][iter_][2] << std::endl;
	ofs << "k_[" << i << "] = " << k_[i] << std::endl;
      }
    }
    else {
      Info << nl << "Calculating k_i's" << nl << endl;
 
      for(label i = 0; i < nMicro_; i++) {
	kIs_[i][iter_][0] = mDot_[i];
	kIs_[i][iter_][1] = f_[i];
	kIs_[i][iter_][2] = 0.1*mDot_[i]; //mDotError_[i] //*** need to pass flow rate error

	// find k_i's
    
	scalar sumTop = 0.0;
	scalar sumBottom = 0.0;

	for (label j=0; j < nIter_; j++) {
	  if(j < iter_) {
	    scalar mDot = kIs_[i][j][0];
	    scalar force = kIs_[i][j][1];
	    scalar mDotError = kIs_[i][j][2];

	    sumTop += force*force/(mDotError*mDotError);
	    sumBottom += force*mDot/(mDotError*mDotError);
	  }
	}

	k_[i] = sumTop/sumBottom;

	ofs << "kIs_[" << i << "][" << iter_ << "][0] = " << kIs_[i][iter_][0] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][1] = " << kIs_[i][iter_][1] << std::endl;
	ofs << "kIs_[" << i << "][" << iter_ << "][2] = " << kIs_[i][iter_][2] << std::endl;
	ofs << "k_[" << i << "] = " << k_[i] << std::endl;
	//Calculate ki errors here //****
      }
    }

    // 2. set LU matrix (this is the hardest and trickiest part)

    std::cout << "Creating matrix" << std::endl;
    simpleMatrix<scalar> luMatrix(nMicro_+1, 0.0, 0.0);

    // column configuration

    // phiCoeff's... mbar

    label r = 0; // row counter

    label c = 0; // column pusher

    // overall pressure drop integral equation

    double L_ = channel.length();

    luMatrix[r][c+0] = L_;

    for(label j=1; j <= (nMicro_-1)/2; j++) {
      luMatrix[r][c+(2*j)-1] = (L_/(2*pi*j))*(1-Foam::cos(2*pi*j));
      luMatrix[r][c+(2*j)+1-1] = L_*(Foam::sin(2*pi*j))/(2*pi*j);
    }

    luMatrix.source()[r]=0.0;

    r++;

    // flow response equations

    for(label i = 0; i < nMicro_; i++) {
      c = nMicro_;

      luMatrix[r][c]=-rhoN_*k_[i];

      c = 0;

      luMatrix[r][c+0]=1;

      for(label j=1; j <= (nMicro_-1)/2; j++) {
	luMatrix[r][c+(2*j)-1] = Foam::sin(2*pi*s_[i]*j/L_);
	luMatrix[r][c+(2*j)+1-1] = Foam::cos(2*pi*s_[i]*j/L_);
      }

      luMatrix.source()[r]=-rhoN_*F_;

      r++;
    }

    //... done

    // 3. solve matrix
    scalarField M = luMatrix.LUsolve();

    // 4. set new phi's

    for(label i = 0; i < nMicro_; i++) {
      phiCoeffs_[i] = M[i];
    }

    for(label i = 0; i < nMicro_; i++) {
      phi_[i] = phiCoeffs_[0];

      for(label j=1; j <= (nMicro_-1)/2; j++) {
	phi_[i] += phiCoeffs_[(2*j)-1]*Foam::sin(2*pi*s_[i]*j/L_);
	phi_[i] += phiCoeffs_[(2*j)+1-1]*Foam::cos(2*pi*s_[i]*j/L_);
      }
    }

    //5. Set new forces

    List<scalar> mDotMacro_(nMicro_,0.0); //predicted Mass Flow Rates

    std::cout << "Setting new forces" << std::endl;
    for(label i = 0; i < nMicro_; i++) {
      f_old_[i] = f_[i];
      f_[i] = (F_ + phi_[i]/rhoN_);
      ofs << "f_[" << i << "] = " << f_[i] << " = " << F_ << " + " << phi_[i] << "/" << rhoN_ << std::endl;

      if(k_[i] != 0.0) {
	mDotMacro_[i] = f_[i]/k_[i];
      }
    }

    // next iteration

    //iter_++; // important

    Info << nl << "Statistics" << endl;
    Info << nl << "id \t s_i \t f_i \t mDotMicro_i \t mDotMacro_i \t k_i"
	 << "\t phi_i \t phiCoeffs_i \t"
	 << endl;

    for(label i = 0; i < nMicro_; i++) {
      Info << i << "\t" << s_[i] << "\t" << f_[i] << "\t" << mDot_[i]
	   << "\t" << mDotMacro_[i] << "\t" << k_[i] << "\t" << phi_[i]
	   << "\t" << phiCoeffs_[i]
	   << endl;
    }

    ofs << std::endl << "Statistics" << std::endl;
    ofs << std::endl << "id \t s_i \t f_i \t mDotMicro_i \t mDotMacro_i \t k_i"
        << "\t phi_i \t phiCoeffs_i \t"
        << std::endl;

    for(label i = 0; i < nMicro_; i++) {
      ofs << i << "\t" << s_[i] << "\t" << f_[i] << "\t" << mDot_[i]
	  << "\t" << mDotMacro_[i] << "\t" << k_[i] << "\t" << phi_[i]
	  << "\t" << phiCoeffs_[i]
	  << std::endl;
    }
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
	double cfdMassFlowRate = massFlowRateConversionFactor * val;
	inVarValuesVec[ifn][var].push_back(cfdMassFlowRate);
	std::cout << "CFD in t = " << t << " " << ifn << ' ' << var << ' ' << cfdMassFlowRate << std::endl;
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
	double mdForce = forceConversionFactor * outVarValues[ifn][var];
	std::cout << "outVarValues[" << ifn << "][" << var << "]" << outVarValues[ifn][var] << std::endl;
        std::cout << "CFD out t = " << t << " " << ifn << ' ' << var << ' ' << mdForce << std::endl;
        interface->push(var, mdForce);
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
      std::cout << "ERROR: not all of the hashed names are distinct." << std::endl;
      // Make sure program terminates even in MPMD mode.
      MPI_Abort(MPI_COMM_WORLD, 999);
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
    nNodes = countNodes(MPI_COMM_WORLD);
    std::cout << "B: nNodes = " << nNodes << std::endl;
    nNodes = countNodes();
    std::cout << "A: nNodes = " << nNodes << std::endl;
    int nMDNodes = nNodes-1;

    int NEVERY = nStepsBetweenSamples; //TODO: replace NEVERY?

    bool halt = false;
    bool terminate = false;
    int64_t timestep = startTime;
    int numberOfMDs = interfaceNames.size();
    bool flowRateHasBeenPreviouslyEstimated;
    double estimatedFlowRate;
    if(numberOfMDs != 0) {
      std::ifstream flowRateFileStream("estimatedFlowRate");
      flowRateHasBeenPreviouslyEstimated = flowRateFileStream.good();
      if(flowRateHasBeenPreviouslyEstimated) { // File exists.
	flowRateFileStream >> estimatedFlowRate;
	flowRateFileStream.close();
      }
    }
    else {
      flowRateHasBeenPreviouslyEstimated = false;
    }
    initialise(numberOfMDs);
    //TODO: replace this!
    //int i = 0;
    std::cout << "Set outputs before run." << std::endl;
    for(std::string ifn : interfaceNames) {
      for(std::string v : outVars) {
	outVarValues[ifn][v] = F_;
      }
      /*
      f_[i] = F_;
      f_old_[i] = F_;
      std::cout << "Before run f_[" << i << "] = " << f_[i] << std::endl;
      std::cout << "Before run f_old_[" << i << "] = " << f_old_[i] << std::endl;
      i++;
      */
    }
    std::cout << "Outputs set before run." << std::endl;

    std::cout << "run: nEquilibration = " << nEquilibration << std::endl;
    std::cout << "run: nMeasurement = " << nMeasurement << std::endl;

    bool hasConverged = (numberOfMDs == 0); // Do not want to test for convergence if the are no MDs.
    bool hasConvergedOverall = false;
    for(int iter_ = 0; iter_ < nIter_; iter_++) {
      int32_t sampleCount = 0;

      for(int i = 0; (i < nEquilibration) && (!terminate) && (!hasConverged) && (!hasConvergedOverall); i++) {
	//std::cout << "run:equib: timestep = " << timestep << std::endl;

	if((timestep % NEVERY) == 0) {
	  std::cout << "run:equib:NEVERY: timestep = " << timestep << std::endl;

	  receiveData(interfaces, timestep);

          sampleCount++;
          if(sampleCount == nSamples) {
	    calculateAverages();
	    // This relies on the regions being correctly ordered.
	    sort(regions.begin(), regions.end());
	    int32_t j = 0;
	    for(Region r : regions) {
	      std::cout << "Sorted: ifn = " << r.interfaceName << std::endl;
	      mDot_[j] = averages[r.interfaceName][std::string("mass_flow_x")].back();
	      std::cout << "equib: mDot_[" << j << "] = " << mDot_[j] << std::endl;
	      j++;
	    }
	    clearAccumulators();
	    sampleCount = 0;
	  }

	  // The initial force is set in the LAMMPS script.
	  sendData(interfaces, timestep);

          double val = double(false);
	  double max_val; // To keep MPI_Allreduce happy. It will be set to the value of val.
	  MPI_Allreduce(&val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  std::cout << "MPI_Allreduce: terminate at timestep " << timestep << " = " << max_val << ", halt = " << halt << std::endl;
	}
	timestep++;
      }

      double meanFlowRate;
      for(int i = 0; (i < nMeasurement) && (!terminate); i++) {
	//std::cout << "run:measurement: iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
	if((timestep % NEVERY) == 0) {
	  //std::cout << "run:measurement:NEVERY: timestep = " << timestep << ", hasConverged = " << hasConverged << ", iter_ = " << iter_ << std::endl;
	  if(hasConverged || iter_ == (nIter_ - 1)) {
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

	      writeCmdFile(nodeDistribution, std::to_string(timestep));

	      if(halt) {
		std::system("mv cmd halt");
	      }
	      terminate = true;
	    }
	    if(hasConverged) {
	      terminate = true;
	    }
	    hasConvergedOverall = convergedOverall(flowRateHasBeenPreviouslyEstimated, estimatedFlowRate, meanFlowRate);
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

	  receiveData(interfaces, timestep);

          sampleCount++;
	  //std::cout << "SAMPLE_COUNT = " << sampleCount << ", i = " << i << ", iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
          if(sampleCount == nSamples) {
	    //std::cout << "AVERAGE" << std::endl;
	    calculateAverages();
	    // This relies on the regions being correctly ordered.
	    sort(regions.begin(), regions.end());
	    int32_t j = 0;
	    for(Region r : regions) {
	      std::cout << "Sorted: ifn = " << r.interfaceName << std::endl;
	      mDot_[j] = averages[r.interfaceName][std::string("mass_flow_x")].back();
	      std::cout << "mDot_[" << j << "] = " << mDot_[j] << std::endl;
	      j++;
	    }
	    clearAccumulators();
	    sampleCount = 0;

	    std::cout << "SOLVE, iter_ = " << iter_ << ", timestep = " << timestep << std::endl;
	    solve(iter_, numberOfMDs);
	    hasConverged = converged(numberOfMDs, meanFlowRate);
	    /*
	    if(iter_ == 2) { // TODO: FOR TEST PUPOSES ONLY - REMOVE EVENTUALLY.
	      hasConverged = true;
	      ofs << "Asserting convergence at iter_ = 2." << std::endl;
	      std::cout << "Asserting convergence at iter_ = 2." << std::endl;
	    }
	    */

	    int32_t posn = 0;
	    for(std::string ifn : interfaceNames) {
	      for(std::string v : outVars) {
		if(v == std::string("force")) {
		  std::cout << "f_[" << posn << "] = " << f_[posn] << std::endl;
		  std::cout << "f_old_[" << posn << "] = " << f_old_[posn] << std::endl;
		  //double deltaF_ = f_[posn] - f_old_[posn];
		  //double deltaF_ = f_old_[posn] - f_[posn];
		  double deltaF_ = f_[posn];
		  outVarValues[ifn][v] = deltaF_;
		  std::cout << "Force delta set to " << deltaF_ << std::endl;
		}
		else {
		  std::cout << "Unexpected variable '" << v << "'." << std::endl;
		  MPI_Abort(MPI_COMM_WORLD, 999);
		}
	      }
	      posn++;
	    }
          }
	  //std::cout << "Outputs set" << std::endl;

          sendData(interfaces, timestep);

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

      /*
      // This relies on the regions being correctly ordered.
      std::cout << "Creating force files at end of iteration " << iter_ << std::endl;
      sort(regions.begin(), regions.end());
      int32_t j = 0;
      for(Region r : regions) {
	std::ostringstream strs;
	//strs << f_[j]; // TODO: reinstate this!
	strs << F_; // TODO: remove this!
	std::string f_str = strs.str();
	createMDForceFile(r, f_str, std::to_string(iter_ + 1));
	j++;
      }
      */
    }
    std::cout << "Leaving run" << std::endl;
  }

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
    std::ifstream agendaInFile("agenda.txt");
    std::ofstream agendaOutFile("new_agenda.txt");
    bool change = false;

    if(continue_counter <= 0) {
      // Get the first line of the file agenda.txt.
      std::string line;
      std::getline(agendaInFile, line);
      std::cout << "AGENDA: " << line << std::endl;
      std::istringstream iss(line);

      std::string token;
      iss >> token;
      if(token != std::string("continue")) {

	// Collect the descriptions of the new MD simulations to be created.
	if(token == std::string("add")) {
	  // Add the new regions.
	  while(iss >> token) {
	    if(token == std::string("subtract"))
	      break;
	    assert(token == "{");
	    iss >> token;
	    assert(std::stoi(token) > maxIndex);
	    Region r;
	    r.interfaceName = token; // Simulation index
	    iss >> token;
	    r.sNorm = std::stod(token);
	    iss >> token;
	    assert(token == "}");
	    add.emplace_back(r);
	    maxIndex++;
	    change = true;
	  }
	}

	while(iss >> token) {
	  // Collect the indices of the MD simulations that are no longer needed.
	  assert(token == "{");
	  iss >> token;
	  Region r;
	  r.interfaceName = token; // Simulation index
	  iss >> token;
	  r.sNorm = std::stod(token);
	  iss >> token;
	  assert(token == "}");
	  remove.emplace_back(r);
	  // Do not change max index as indices are never reused.
	  change = true;
	}

      }

      // Read how many times the program should continue with the current number of MD sims.
      iss >> continue_counter;
      std::cout << "continue_counter = " << continue_counter << std::endl;

      // Copy all but the first line of the file agenda.txt
      // to new_agenda.txt.
      while(std::getline(agendaInFile, line)) {
	agendaOutFile << line << std::endl;
      }
      std::cout << "Closing agenda.txt" << std::endl;
      agendaInFile.close();
      std::cout << "Closing new_agenda.txt" << std::endl;
      agendaOutFile.close();

      // Update agenda.txt.
      std::cout << "Moving new_agenda.txt" << std::endl;
      std::system("mv new_agenda.txt agenda.txt");
    }

    continue_counter--;
    //std::cout << "After decrement continue_counter = " << continue_counter << std::endl;

    std::cout << "Leaving changesRequired: change = " << change << std::endl;
    return change;
  }

  bool CFDSim::converged(int nMicro_, double& mean) {
    if(nMicro_ == 0) {
      return true;
    }

    double sum = 0.0;
    for(label i = 0; i < nMicro_; i++) {
      sum += mDot_[i];
    }

    mean = sum / nMicro_;

    double variance = 0.0;
    for(label i = 0; i < nMicro_; i++) {
      variance += pow(mDot_[i] - mean, 2);
    }

    double standardError = std::sqrt(variance / nMicro_);

    double normalisedError = standardError / std::abs(mean);

    ofs << "The mean is " << mean << " and the normalised error is " << normalisedError << std::endl;

    std::ofstream flowRateFileStream("estimatedFlowRate");
    flowRateFileStream << mean << std::endl;
    flowRateFileStream.close();

    bool hasConverged = normalisedError < acceptableError;
    ofs << "hasConverged = " << hasConverged << std::endl;
    return hasConverged;
  }

  bool CFDSim::convergedOverall(bool previousEstimateExists, double previousEstimate, double meanFlowRate) {
    if(!previousEstimateExists) {
      return false;
    }

    std::cout << "previousEstimate = " << previousEstimate << ", meanFlowRate = " << meanFlowRate << std::endl;

    double normalisedChange = std::abs(previousEstimate - meanFlowRate) / meanFlowRate;

    ofs << "The overall normalised change is " << normalisedChange << std::endl;

    bool hasConvergedOverall = normalisedChange < tolerance;
    ofs << "hasConvergedOverall = " << hasConvergedOverall << std::endl;
    return hasConvergedOverall;
  }
}
