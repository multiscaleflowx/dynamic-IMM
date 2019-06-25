In this software a steady state flow is modelled.

In brief, the software works as follows. There is a CFD component and two or more MD components that work together to generate a solution. The CFD component maintains an overview of the system which is a channel whilst each MD simulation models part of the channel.

In more detail, the CFD component of the simulator sends a body force to each of the MD simulations. (In the first instance the same body force is sent to each MD simulation.) After a number of time steps each MD simulation passes an instantaneous mass flow rate to the CFD component. For each of the selected locations in the channel a number of mass flow rates are accumulated and an average is computed. The average flow rates are then used to compute a new body force for each MD simulation (through the use of a linear solver) and process of computation starts again. This is carried out for a number of iterations that is specified at run time.

The OpenFOAM code requires a directory called 'channel' to exist in the directory in which the simulation is run. In this directory there must be two executables:
1) the OpenFOAM executable cfdSim
2) the LAMMPS executable mui-lmp.

The simulation is controlled by a number of files. The names of the files are given relative to the directory in which the simulation is run.

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

The code simulates a channel the height (in the y direction) of which varies with the streamwise location (x direction). The channel is periodic in the z direction. The shape of the channel can be changed by editing channel.C. The default shape is that of a channel that shrinks linearly in height until it reaches its minimum at the mid point and then expands to its original height. The precise dimensions are controlled by setting the values of 'length', 'hEnd' and 'hNeck'.

'rho' is the macroscopic density of the fluid whilst 'F' is an externally applied body force.

