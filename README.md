AS yet, no real CFD program (such as OpenFOAM) is used. Instead there is a program which exchanges values with LAMMPS simulations and which decides when to change the number of simulations. These decisions are read from a file 'agenda'txt' at run time. This program is encoded in the files CFDSim.cpp and CFDSim.h.

To run this example (on ARCHER)
1) compile LAMMPS (from the mui-sequence branch) to produce an executable called mui-lmp.
2) Compile CFDSim.cpp to produce an executable called pseudoOpenfoam.
3) Create a directory in your work file system containing
    CP
    agenda.txt.orig
    cmd.orig
    in.CFD.orig
    initial_MD_template
    restarted_MD_template
    submit1.pbs
    submit2.pbs
    submit3.pbs
    pseudoOpenfoam
    mui-lmp
4) Run the commands
    source CP
    qsub -q short submit1.pbs
    mv halt cmd
    qsub -q short submit2.pbs
    mv halt cmd
    qsub -q short submit2.pbs
    mv halt cmd
    qsub -q short submit3.pbs
   noticing how files are created and overwritten.
   In particular notice how agenda.txt shrinks and the in.MD* files change.


CP copies agenda.txt.orig to agenda.txt, cmd.orig to cmd and in.CFD.orig to in.CFD. agenda.txt, cmd and in.CFD are updated by the sequence of command given in (4) above.

initial_MD_template and restarted_MD_template are templates for LAMMPS scripts used to either initiate a LAMMPS run or restart one. The names of the files created have the form in.MD<integer>.

Apart from the name of the job created, the submit*.pbs files differ only in the number of nodes that they request.


Compiling the Software on ARCHER
--------------------------------

git clone https://github.com/MxUI/MUI.git

Compiling LAMMPS
----------------

GIT_SSL_NO_VERIFY=true git clone https://git.ecdf.ed.ac.uk/multiscale/lammps.git

module swap PrgEnv-cray PrgEnv-gnu
module load fftw/3.3.4.5
module load cmake

cd lammps
git checkout mui-sequence
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=CC -DBUILD_MPI=ON -DBUILD_LIB=ON -DMUI_INCLUDE_DIR=../../test/MUI -DPKG_USER-MUI=ON ../cmake
make -j 8

Compiling CFDSim.cpp
--------------------

CC -g -O3 -std=c++11 -I<path_to_mui_headers> CFDSim.cpp -o pseudoOpenfoam
