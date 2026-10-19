# Molecular Dynamics Basics in MATLAB

Coursework for MNT 756: Molecular Modeling (Homework 1), North Dakota State University.

A small set of MATLAB simulations that implement molecular dynamics from scratch,
to explore how integration schemes and model parameters affect particle trajectories
and energy conservation.

## What it simulates

1. **Particle in a 2D circular box** – a single particle moving inside a circular boundary,
   kept in by a harmonic wall potential, U = ½·k_wall·(r − r_cut)², and integrated with a
   first-order (forward Euler) update. The trajectory is plotted as it runs.
2. **Particle in a 1D double-well potential** – integrated with the Verlet algorithm,
   x(t+dt) = 2x(t) − x(t−dt) + dt²·f(t)/m. Position, velocity, force and total energy
   are plotted to check energy conservation.
3. **Parameter study** – repeats the simulation with a heavier particle, a larger
   cutoff distance and a stiffer wall to see how each changes the dynamics.

## How to run

Open `twoD_single_particle_circle_MZ[final].m` in MATLAB and run it.
Initial conditions are random, so each run differs slightly.

## Files

- `twoD_single_particle_circle_MZ[final].m` – simulation script
- `w2.mlx` – MATLAB Live Script version
