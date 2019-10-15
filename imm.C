#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "CFDSim.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char ** argv) {
  MPI_Comm local_comm = mui::mpi_split_by_app();
  int size, rank, local_size, local_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(local_comm, &local_size);
  MPI_Comm_rank(local_comm, &local_rank);
  std::cout << "size, rank: " << size << " " << rank << std::endl;
  std::cout << "local size, local rank: " << local_size << " " << local_rank << std::endl;
  std::cout << argv[0] << ": " << argv[1] << " " << argv[2] <<std::endl;

#include "setRootCase.H"
#include "createTime.H"

  IOdictionary immDict = (IOobject("immDict",
				   "",
				   runTime, 
				   IOobject::MUST_READ,
				   IOobject::NO_WRITE,
				   false
				   )
			  );

  auto sim = cfdsim::CFDSim(argv[0], "in.CFD", argv[2], immDict);
  sim.run();
  
  return 0;
}
