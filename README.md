# N-body Simulation: Algorithm and Implementation Details

This document provides a comprehensive overview of the algorithmic implementation used in this N-body simulation codebase. The primary goal of the simulation is to evolve a system of particles under their mutual gravitational (or electrostatic) interaction, with optional support for grid-based field solvers.

The purpose of this code is to develop a mini app that enhances my understanding of the N-body problem and serves as a practical application for CINECA's Summer School 2026. My understanding has significantly improved after reading this article: https://arborjs.org/docs/barnes-hut.

At this point: the code is not fully optimized. 

---

## 1. Overview of the N-body Problem

The N-body problem concerns predicting the motion of N particles under mutual forces (typically gravity or Coulomb forces). The brute-force method computes all pairwise interactions, leading to O(N²) computational cost. More advanced algorithms, such as the Barnes-Hut tree code, can reduce this to O(N log N) by using hierarchical multipole expansions.

---

## 2. Code Structure and Initialization

### 2.1. Particle Initialization

Particles are initialized with positions and velocities, supporting different distributions:

- **Random Uniform Distribution:** Particles are placed randomly within a cubic box, with small random initial velocities.
- **Plummer Sphere:** Particles are distributed according to the Plummer model, a common approximation for spherical stellar systems. Positions are assigned using rejection sampling, and velocities are set to approximate virial equilibrium.

**Relevant code:**
```c
// Plummer sphere initialization (particle.c)
double radius = scale_radius / sqrt(pow(r1, -2.0/3.0) - 1.0);
// Random direction
double theta = acos(2.0 * r2 - 1.0);
double phi = 2.0 * M_PI * r3;
// Assign positions
sys->particles[i].x = radius * sin(theta) * cos(phi);
sys->particles[i].y = radius * sin(theta) * sin(phi);
sys->particles[i].z = radius * cos(theta);
// Virial velocities (simplified)
double v_escape = sqrt(2.0 * sys->n_particles / radius);
// ...
```

---

## 3. Main Simulation Loop

The simulation advances in discrete time steps. At each step, the following sequence is executed (see `main.c`):

1. **Reset Forces:** Set all particle forces to zero.
2. **Compute Gravitational Forces:**  
    - **Direct Method:** All pairwise forces are computed (`O(N²)`).
    - **Barnes-Hut Tree:** Space is recursively subdivided into octants, and distant groups of particles are approximated by their multipole expansion (`O(N log N)`).
3. **Grid-Based Forces (Optional):** If enabled, particle densities are mapped to a grid, and the Poisson equation is solved using FFT to obtain grid-based forces.
4. **Integrate Equations of Motion:** Particle positions and velocities are updated using either the leapfrog integrator (default) or an alternative, such as RK4.
5. **Output and Logging:** Particle states are periodically written to files for analysis and visualization.

**Relevant code:**
```c
for (step = 1; step <= n_steps; step++) {
    reset_forces(sys);
    if (params.use_tree) {
        build_tree(tree, sys);
        compute_gravity_tree(sys, tree, params.softening);
    } else {
        compute_gravity_direct(sys, params.softening);
    }
    if (params.use_grid) {
        particles_to_grid(sys, grid);
        solve_poisson_fft(grid);
        compute_electric_field(grid);
        grid_to_particles(sys, grid);
    }
    integrate_leapfrog(sys, params.dt);
    // Output, logging, etc.
}
```

---

## 4. Force Calculation

### 4.1. Direct Summation

For each particle, sum the force contribution from every other particle, applying softening to avoid singularities at short distance.

**Relevant code:**
```c
for (int i = 0; i < sys->n_particles; i++) {
    double fx = 0.0, fy = 0.0, fz = 0.0;
    for (int j = 0; j < sys->n_particles; j++) {
        if (i != j) {
            // Compute pairwise force with softening
        }
    }
    sys->particles[i].fx = fx;
    sys->particles[i].fy = fy;
    sys->particles[i].fz = fz;
}
```

### 4.2. Barnes-Hut Tree Algorithm

- **Tree Construction:** Particles are inserted into an octree, recursively subdividing space until each leaf contains a small number of particles.
- **Force Calculation:** For each particle, traverse the tree. For distant nodes, use the node's center of mass and total mass; for nearby nodes, recurse into children or compute directly.
- **Multipole Acceptance Criterion:** The ratio `s / r < θ` (where `s` is node size, `r` is distance, and θ is the opening angle) determines whether to approximate or recurse.

**Relevant code:**
```c
if (node_size / r < tree->theta) {
    // Use multipole approximation
    // Add force from this node's total mass at its center of mass
} else {
    // Recurse to children
}
```

---

## 5. Grid-Based Poisson Solver (Optional)

- **Density Assignment:** Particle masses are assigned to a 3D grid.
- **Poisson Equation:** The gravitational potential is computed via FFT.
- **Field Interpolation:** The grid-based force field is interpolated back to the particle positions.

---

## 6. Time Integration

### 6.1. Leapfrog Integrator (Default)

Symplectic and time-reversible, the leapfrog method alternates velocity and position updates:

1. Drift: `x_half = x + v * (dt/2)`
2. Kick: `v_new = v + a(x_half) * dt`
3. Drift: `x_new = x_half + v_new * (dt/2)`

### 6.2. Runge-Kutta 4 (Optional)

A higher-order, less commonly used method for N-body, provided for demonstration.

```c
void integrate_leapfrog(ParticleSystem *sys, double dt);
void integrate_rk4(ParticleSystem *sys, double dt);
```

---

## 7. Parallelization and GPU Support

- The code uses OpenACC pragmas for GPU acceleration in force computation and particle updates, where feasible.
- Tree-based force calculation is performed on CPU due to irregular memory access.

---

## 8. Output and Analysis

- Particle positions, velocities, and other properties are periodically written to text and binary files.
- Energy and other statistics are logged for post-analysis.
- A Python script (`nbodyplot.py`) is provided for 3D visualization of particle states.

---

## 9. References

- Barnes, J., & Hut, P. (1986). "A hierarchical O(N log N) force-calculation algorithm". Nature.
- Hockney, R. W., & Eastwood, J. W. (1988). "Computer Simulation Using Particles". CRC Press.

---

## 10. Summary

This implementation balances performance and accuracy by supporting both direct and tree-based force calculations, optional grid-based solvers, and efficient time integration.
---

**For further exploration, refer to the [source code](https://github.com/nzeal/nbody/tree/053ccd68796730cdec29429bf319bceefd868208/source) and [README](https://github.com/nzeal/nbody/blob/053ccd68796730cdec29429bf319bceefd868208/README.md).**
