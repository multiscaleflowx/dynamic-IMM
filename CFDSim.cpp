#include <map>
#include <string>
#include <istream>
#include <ostream>
#include "CFDSim.h"

namespace cfdsim {

  bool operator==(const Region& lhs, const Region& rhs) {
    return lhs.interfaceName == rhs.interfaceName && \
      lhs.ylo == rhs.ylo && \
      lhs.yhi == rhs.yhi;
  }
  
  bool operator <(const Region& r1, const Region& r2) {
    return r1.interfaceName < r2.interfaceName;
  }

  std::ostream& operator<<(std::ostream& os, const Region& r) {
    os << "{ " << r.interfaceName << " " << r.ylo << " " << r.yhi << " }";
    return os;
  }

  std::istream& operator>>(std::istream& is, Region& r) {
    std::string s;
    bool error = false;
    is >> s;
    if(s != "{") {
      error = true;
    }
    else {
      is >> s;
      r.interfaceName = s;
      is >> s;
      r.ylo = s;
      is >> s;
      r.yhi = s;
      is >> s;
      if(s != "}") {
	error = true;
      }
    }
      
    if(error)
      is.setstate(std::ios::failbit);
    return is;
  }

  CFDSim::CFDSim(char* name): cfdFileName(name) {
    readConfigFile();
  }

  void CFDSim::readConfigFile() {
    std::cout << "Entering readConfigFile" << std::endl;
    std::ifstream infs(cfdFileName);
    std::string keyWord, line, s;

    infs >> keyWord;
    assert(keyWord == "maxNumberOfIterations");
    infs >> iters;

    infs >> keyWord;
    assert(keyWord == "startAtIteration");
    infs >> startTime;

    infs >> keyWord;
    assert(keyWord == "push");
    getline(infs, line);
    std::istringstream str_stream2(line);
    while(str_stream2 >> s) {
      outVars.insert(s);
    }

    infs >> keyWord;
    assert(keyWord == "fetch");
    getline(infs, line);
    std::istringstream str_stream3(line);
    while(str_stream3 >> s) {
      inVars.insert(s);
    }

    // A firstTimeIndex corresponds to an MD simulation that is being run for the first time.
    // The simulation has already been started in parallel with OpenFOAM using an input file
    // (in.MD<firstTimeIndex>) which does not reference a restart file (as no restart file exists).
    // The next time the simulation is run a restart file will exist and the input file has to be
    // changed accordingly.
    infs >> keyWord;
    assert(keyWord == "firstTimeRegions");
    getline(infs, line);
    std::istringstream str_stream4(line);
    Region firstTimeRegion;
    while(str_stream4 >> firstTimeRegion) {
      firstTimeRegions.insert(firstTimeRegion);
    }

    // Initially OpenFOAM will be run alone and no interfaces will be needed.
    // Once an MD simulation is required OpenFOAM will be restarted with interfaces
    // to the required MD simulations.
    // The lower and upper bounds of the y direction for each
    // simulation region are specified.
    infs >> keyWord;
    assert(keyWord == "regions");
    getline(infs, line);
    std::istringstream str_stream1(line);
    Region r;
    while(str_stream1 >> r) {
      interfaceNames.emplace_back(r.interfaceName);
      regions.emplace_back(r);
    }

    infs >> keyWord;
    assert(keyWord == "maxIndex");
    infs >> maxIndex;
    // maxIndex should be set to 0 in the CFD config file if there are no MDs.
    if(interfaceNames.size() == 0)
      assert(maxIndex == 0);

    infs.close();
    std::cout << "Leaving readConfigFile" << std::endl;
  }

  void CFDSim::writeConfigFile(int t) {
    std::cout << "Entering writeConfigFile" << std::endl;
    std::ifstream cfdInFS(cfdFileName);
    std::string newCfdFileName("new_"+cfdFileName);
    std::ofstream cfdOutFS(newCfdFileName);
    std::string cfdline;
    std::getline(cfdInFS, cfdline);
    cfdOutFS << cfdline <<std::endl; // Write maximum number of iterations.
    std::getline(cfdInFS, cfdline); // Discard and replace.
    cfdOutFS << "startAtIteration " << t+1 << std::endl;
    std::getline(cfdInFS, cfdline);
    cfdOutFS << cfdline <<std::endl; // Write push vars.
    std::getline(cfdInFS, cfdline);
    cfdOutFS << cfdline <<std::endl; // Write fetch vars.
      
    cfdOutFS << "firstTimeRegions";
    for(auto r : add) {
      cfdOutFS << ' ' << r;
    }
    cfdOutFS << std::endl;

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
    outfile << "aprun -n 1 pseudoOpenfoam " << cfdFileName;
    int i = 0;
    for(auto r : updatedRegions) {
      int p = nodeDistribution[i] * 24;
      outfile << " :  -n " << p << " lmp_xc30 -in in.MD" << r.interfaceName;
      i++;
    }
    outfile << " > output" << t << std::endl;
    outfile.close();
    std::cout << "Leaving writeCmdFile" << std::endl;
  }

  void CFDSim::createInitialMDFile(Region r, std::string startTime) {
    std::cout << "Entering createInitialMDFile" << std::endl;
    std::ifstream infs("initial_MD_template");
    std::ofstream outfs("in.MD"+ r.interfaceName);
    std::string line;
    while(getline(infs, line)) {
      size_t index;
      while((index = line.find("$i")) != std::string::npos) {
	line.replace(index, 2, r.interfaceName);
      }
      while((index = line.find("$ylo")) != std::string::npos) {
	line.replace(index, 4, r.ylo);
      }
      while((index = line.find("$yhi")) != std::string::npos) {
	line.replace(index, 4, r.yhi);
      }
      while((index = line.find("$t")) != std::string::npos) {
	line.replace(index, 2, startTime);
      }
      outfs << line << std::endl;
    }
    std::cout << "Leaving createInitialMDFile" << std::endl;
  }

  void CFDSim::createRestartedMDFile(Region r, std::string startTime) {
    std::cout << "Entering createRestartedMDFile" << std::endl;
    std::ifstream infs("restarted_MD_template");
    std::ofstream outfs("in.MD"+ r.interfaceName);
    std::string line;
    while(getline(infs, line)) {
      size_t index;
      while((index = line.find("$i")) != std::string::npos) {
	line.replace(index, 2, r.interfaceName);
      }
      while((index = line.find("$ylo")) != std::string::npos) {
	line.replace(index, 4, r.ylo);
      }
      while((index = line.find("$yhi")) != std::string::npos) {
	line.replace(index, 4, r.yhi);
      }
      while((index = line.find("$t")) != std::string::npos) {
	line.replace(index, 2, startTime);
      }
      outfs << line << std::endl;
    }
    std::cout << "Leaving createRestartedMDFile" << std::endl;
  }

  bool CFDSim::changeMDs(int t) {
    std::cout << "Entering changeMDs" << std::endl;
    bool change = false;

    updatedRegions.clear(); // This is necessary as changeMDs may be called more than once.
    // Make a copy of interfaceNames.
    for(auto r : regions) {
      updatedRegions.push_back(r);
    }

    if(changesRequired()) {

      // Remove the interface names that are no longer needed.
      for(auto r : remove) {
	auto it = std::find(updatedRegions.begin(), updatedRegions.end(), r);
	updatedRegions.erase(it);
      }

      // Those MDs that existed before this run and are still needed must have
      // restart files created for them.
      for(auto r : updatedRegions) {
	createRestartedMDFile(r, std::to_string(t));
      }

      // Add the new interface names and create initial restart files.
      for(auto r : add) {
	updatedRegions.emplace_back(r);
	createInitialMDFile(r, std::to_string(t));
      }

      writeConfigFile(t);

      change = true;
    }

    std::cout << "Leaving changeMDs" << std::endl;
    return change;
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
    std::cout << "calculateNodeDistribution: nMDs = " << nMDs << ", nNodes = " << nNodes << std::endl;
    if(nMDs > nNodes) { // Arbitrarily assigne 1 node to each MD.
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

  void CFDSim::exchangeData(std::vector<std::unique_ptr<mui::uniface<mui::config_3d>>>& interfaces,
			    int step) {    
    std::cout << "Entering exchangeData" << std::endl;
    std::map<std::string, std::map<std::string, double>> inVarValues;
    std::map<std::string, std::map<std::string, double>> outVarValues;

    double time_dt = 0.003;
    double allowed_lag_time = 2*time_dt;
    double t = step*time_dt;

    calculateOutputs(t, outVarValues);

    int i = 0;
    for(auto& interface : interfaces) {
      std::string ifn = interfaceNames[i];

      // SEND.
      for(std::string var : outVars) {
	std::cout << "CFD out t = " << t << " " << ifn << ' ' << var << ' ' << outVarValues[ifn][var] << std::endl;
	interface->push(var, outVarValues[ifn][var]);
      }
      interface->commit(t);

      // RECEIVE.
      interface->barrier(t);
      for(std::string var : inVars) {
	inVarValues[ifn][var] = interface->fetch<double>(var);
	std::cout << "CFD in t = " << t << " " << ifn << ' ' << var << ' ' << inVarValues[ifn][var] << std::endl;
      }

      i++;
    }
    std::cout << "Leaving exchangeData" << std::endl;
  }

  void CFDSim::run() {
    std::cout << "Entering run" << std::endl;

    std::cout << "CFD creating interfaces." << std::endl;
    auto interfaces = mui::create_uniface<mui::config_3d>(std::string("CFD"), interfaceNames);
    std::cout << "CFD interfaces created." << std::endl;

    int nNodes = countNodes(MPI_COMM_WORLD);
    int nMDNodes = nNodes-1;

    int NEVERY = 10;

    bool halt = false;
    bool terminate = false;
    for(int t = startTime; t <= iters; t++) {
      if((t % NEVERY) == 0) {
	std::cout << "run: t = " << t << std::endl;
	int oldNumberOfMDs = interfaceNames.size();
	if(changeMDs(t)) {
	  int numberOfMDs = updatedRegions.size();
	  std::cout << "t = " << t << ": numberOfMDs = " << numberOfMDs << std::endl;
	  std::cout << "t = " << t << ": maxIndex = " << maxIndex << std::endl;
	  std::vector<int> nodeDistribution(numberOfMDs);
	  if((numberOfMDs == 0 ) && (oldNumberOfMDs > 0)) {
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

	  writeCmdFile(nodeDistribution, std::to_string(t));

	  if(halt) {
	    std::system("mv cmd halt");
	  }
	  terminate = true;
	}

	exchangeData(interfaces, t);
	double val = double(terminate);
	double max_val; // To keep MPI_Allreduce happy. It will be set to the value of val.
	MPI_Allreduce(&val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      }
      if(terminate) {
	break;
      }
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
  bool CFDSim::changesRequired() {
    std::cout << "Entering changesRequired" << std::endl;
    std::ifstream agendaInFile("agenda.txt");
    std::ofstream agendaOutFile("new_agenda.txt");
    bool change = false;

    // Get the first line of the file agenda.txt.
    std::string line;
    std::getline(agendaInFile, line);
    std::istringstream iss(line);

    std::string token;
    iss >> token;
    if(token != std::string("continue")) {

      // Collect the indices of the new MD simulations to be created.
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
	  r.ylo = token;
	  iss >> token;
	  r.yhi = token;
	  iss >> token;
	  assert(token == "}");
	  add.emplace_back(r);
	  maxIndex++;
	  change = true;
	}
      }

      while(iss >> token) {
	// Collect the indices of the MD simulations that are no longer needed.
	Region r;
	r.interfaceName = token; // Simulation index
	iss >> token;
	r.ylo = token;
	iss >> token;
	remove.emplace_back(r);
	// Do not change max index as indices are never reused.
	change = true;
      }

    }
  
    // Copy all but the first line of the file agenda.txt
    // to new_agenda.txt.
    while(std::getline(agendaInFile, line)) {
      agendaOutFile << line << std::endl;
    }
    agendaInFile.close();
    agendaOutFile.close();

    // Update agenda.txt.
    std::system("mv new_agenda.txt agenda.txt");

    std::cout << "Leaving changesRequired" << std::endl;
    return change;
  }

  /*
    The CFD class must implement a method to replace that given below.
    The implementation in this file is provided only to allow testing of the rest of the code.

    This function populates the map varValues. The keys of the map are the indices of the interfaces.
    The values are (variable name, value) pairs.
  */
  void CFDSim::calculateOutputs(int t,
				std::map<std::string, std::map<std::string, double>>& varValues) {
    std::cout << "Entering calculateOutputs" << std::endl;
    double val = 0.0;
    for(std::string ifn : interfaceNames) {
      for(std::string v : outVars) {
	varValues[ifn][v] = val;
	val = val + 1.0;
      }
    }
    std::cout << "Leaving calculateOutputs" << std::endl;
  }
}

int main(int argc, char ** argv) {
  MPI_Comm local_comm = mui::mpi_split_by_app();
  int size, rank, local_size, local_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(local_comm, &local_size);
  MPI_Comm_rank(local_comm, &local_rank);
  std::cout << "size, rank: " << size << " " << rank << std::endl;
  std::cout << "local size, local rank: " << local_size << " " << local_rank << std::endl;

  auto sim = cfdsim::CFDSim(argv[1]);
  sim.run();
  
  return 0;

}
