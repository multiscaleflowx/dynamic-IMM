#include <iostream>
#include <ostream>
#include <fstream>
#include <mpi.h>

#include "common.H"

namespace cfdsim {

  // To be used when there is an error.
  void haltMPMD(const char * msg) {
    std::string message("ERROR: ");
    message.append(msg);
    std::ofstream outfile("halt");
    outfile << message << std::endl;
    outfile.close();
    std::cout << message << std::endl << std::flush;
    std::cerr << message << std::endl << std::flush;
    MPI_Abort(MPI_COMM_WORLD, 999);
  }
}
