# Airgun1D with shuttle integration
Code for source model used in *Influence of port opening dynamics on the acoustic signature of pneumatic marine seismic sources* (2023). This model invokes the quasi-one-dimensional Euler equations to model the gas in the pneumatic seismic source (airgun) coupled to a bubble and a shuttle constrained by the geometry of the airgun.

This repository extends the 1D airgun code of Leighton Watson et al. at [leighton-watson/AirGun1D](https://github.com/leighton-watson/AirGun1D) (more details within), which was used in [*What controls the initial peak of an air-gun source signature?*](https://library.seg.org/doi/full/10.1190/geo2018-0298.1) (2019). This version includes the effect of the shuttle on the gas discharge. The shuttle dynamics are controlled by the evolving gas pressures applied to either side of the shuttle. See the paper associated with this work for a schematic of the modeled source design, which is based on patent US9804280.

This repository uses sbplib (https://github.com/scicompuu/sbplib), a MATLAB library for summation-by-parts (SBP) finite-difference operators developed at Uppsala University. A copy of sbplib is included in this repository.

Field data used to generate a portion of the figures are not disclosed.

## Getting started

To run the model with default parameters, navigate to the root directory of this repository. The following lines of example code will run the code with default parameters, postprocess the solution into a human-readable struct, and then plot pressure for a specific time index.

```matlab
cd ./AirGun1D;           % Navigate to AirGun1D
nx = 40;                 % Number of grid points per meter of the 1D domain of the firing chamber
coupleToShuttle = true;  % Whether to control the port area by modeling the shuttle dynamics (true/false)
options = struct();      % Struct containing optional parameters; uses only default parameters when empty

% Run model and return a struct containing the solution and a struct containing the parameters used in the run
[solution, metadata] = airgunShuttleDeploy(nx, coupleToShuttle, options);

% (Postprocessing example) Postprocess solution at last time index i
i = size(solution.q, 1);
compute_all_variables = true; % Setting to compute all primitive variables in the 1D firing chamber

state = ...
    metadata.discretization.fullState( ...
    solution.q(:,i), ...       % Firing chamber 1D solution in outer product ordering, at time i
    solution.soln.x(i), ...    % Time vector for main ODE solve
    solution.bubble(:,i), ...  % Bubble state vector at time i
    solution.shuttle(:,i), ... % Shuttle state vector (set to dummy value [] if coupleToShuttle==false)
    ~coupleToShuttle, ...
    compute_all_variables);

% Access postprocessed data
eulerDomainStates = [state.eulerDomainStates];
bubbleStates =      [state.bubbleStates];
portStates =        [state.portStates];
shuttleStates =     [state.shuttleStates];

% Get 1D domain mesh in firing chamber
x = metadata.discretization.schm.u;

% Plot pressure
plot(x, [eulerDomainStates.p]);
```

### What's in my output?

* `solution` is a struct containing the numerical solution for all times, represented as matrices with rows corresponding to each degree of freedom in the model and columns corresponding to each time index.
* `metadata` contains data about the run, the parameters and initial conditions used, and the discretization object used to compute the solution.
* `state` is a struct of structs pertaining to each submodel (port, 1D Euler domain, bubble, shuttle); each field of `state` contains named variables corresponding to the variable in the paper (e.g., `p` for pressure, `M` for Mach number, etc.)

### How do I change model parameters?

* TODO

## Generating paper figures

The figures used in the associated publication are produced using code in [FigGen/](/AirGun1D/FigGen/), after running the model for the reference case using [runReferenceCase.m](/AirGun1D/Postprocess/runReferenceCase.m). Note that some figures have additional superimposed data, or require the field data to generate the plot.

## Simplified description of code structure

The wrapper [`airgunShuttleDeploy`](/AirGun1D/airgunShuttleDeploy.m
) calls [`runEulerCodeShuttleDual`](/AirGun1D/SBPSAT/runEulerCodeShuttleDual.m), which contains the code for the case where the port area is controlled by shuttle dynamics, and the case where the port area is set to the maximum port area at `t = 0`. The function `runEulerCodeShuttleDual` initializes the model and feeds the ODE right-hand-side (coming from the finite-difference discretization of the PDEs combined with the shuttle and bubble ODEs) to the ODE solver (`ode45` with maximum timestep given by the CFL condition). The ODE right-hand-side is defined in the class definition of [`DiscrAirgunShuttleMulti`](AirGun1D/SBPSAT/@DiscrAirgunShuttleMulti/DiscrAirgunShuttleMulti.m) . Useful functions of the state variables, such as the boundary flow state, are computed in [`fullState`](AirGun1D/SBPSAT/@DiscrAirgunShuttleMulti/fullState.m).
