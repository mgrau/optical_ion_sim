# optical_ion_sim
MATLAB class for simulating motion of ions trapped in optical potentials

This is a simple MATLAB class for simulating the motion of ions in optical potentials. Example usage might be

```
% initialize simulation with 50 ions
sim = ion_simulation(50)
% set initial positions of ions in a single 'pancake' of a 1D lattice
sim.init_pancake()
% run the Runge-Kutta simulation
sim.run()
% play a nice movie of the ion motion
sim.movie()
% plot the normal modes, assuming the final position is the equilibrium position
sim.normal_modes()
```
