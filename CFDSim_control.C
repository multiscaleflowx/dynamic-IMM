#include <map>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <set>          // std::set
#include <vector>       // std::vector

#include "CFDSim.H"

namespace cfdsim {

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
      if((continue_counter != 0) && (continue_counter != (nSamples-1))){
	std::cout << "ERROR: (continue_counter = " << continue_counter <<") != ((nSamples-1) = " << (nSamples-1)  << ")" << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 999);
      }

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
}
