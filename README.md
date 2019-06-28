In this software a steady state flow in a channel is modelled. The height of the channel (in the y direction) varies with the streamwise location (x direction). The channel is periodic in the z direction.

In brief, the software works as follows. There is a CFD component and two or more MD components that work together to generate a solution. The CFD component maintains an overview of the system which is a channel whilst each MD simulation models part of the channel.

In more detail, the CFD component of the simulator sends a body force to each of the MD simulations. (In the first instance the same body force is sent to each MD simulation.) After a number of time steps each MD simulation passes an instantaneous mass flow rate to the CFD component. For each of the selected locations in the channel a number of mass flow rates are accumulated and an average is computed. The average flow rates are then used to compute a new body force for each MD simulation (through the use of a linear solver) and process of computation starts again. This is carried out for a number of iterations that is specified at run time.

The CFD simulation is implemented using OpenFOAM whilst the MD simulations are implemented using LAMMPS.

In the job submission directory there must be two executables
1) the OpenFOAM executable cfdSim
2) the LAMMPS executable mui-lmp.
The OpenFOAM executable requires a directory called 'channel'  and a file called in.CFD to exist in the directory in which the simulation is run.

The simulation is controlled by a number of files as described below. The names of the files are given relative to the directory in which the simulation is run.

channel/macrodict sets the values of the constants
length
hEnd
hNeck
rho
F
nSolverCalls
numberOfSamplesToAverageOver

channel/microDict sets the values of the constants
molMass
lengthOfRegion
widthOfRegion
microTimeStep
numberOfEquilibrationSteps
numberOfStepsBetweenSamples

in.CFD set the values of
startAtIteration
push
fetch
regions
maxIndex
Here
1) startAtIteration specifies at what time step the MD simulations are to start,
2) push is followed by a list of names of variables the values of which are to be communicated to the MD simulations by CFD (in this case 'force'),
2) fetch is followed by a list of names of variables the values of which are to be communicated to CFD by the MD simulations (in this case 'mass_flow_x'),
3) regions is followed by a list of specifications of regions each of which to be simulated using MD, and
4) maxIndex specifies the highest index that has yet been used to label a MD simulation.
For example, in.CFD could have the form
startAtIteration 21
push force
fetch mass_flow_x
regions { 1 0 } { 2 0.5 }
maxIndex 2
This starts two MD simulations, labelled 1 and 2, of regions located at the streamwise positions 0 and 0.5. The locations are specified by normalised coordinates (i.e. coordinates in the range 0 to 1 with 0 representing the entrance to the channel and 1 the exit from it.

agenda.txt
As things currently stand the CFD program does not make run time decisions about when to start and stop MD simulations. Instead these decisions are encoded in this file.

At the end of a run the CFD simulation (if appropriate) writes out files which may be used to start a new run. The scripts required to run further MD simulations are constructed by editing template files whilst in.CFD is constructed directly. Also, a file is written which contains a suggested 'aprun' command. 

The shape of the channel can be changed by editing channel.C. The default shape is that of a channel that shrinks linearly in height until it reaches its minimum at the mid point and then expands to its original height. The precise dimensions are controlled by setting the values of 'length', 'hEnd' and 'hNeck'.

'rho' is the macroscopic density of the fluid whilst 'F' is an externally applied body force.

