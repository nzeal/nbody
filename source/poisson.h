// ===== poisson.h =====
#ifndef POISSON_H
#define POISSON_H

#include "particle.h"

typedef struct {
    int nx, ny, nz;
    double dx, dy, dz;
    double xmin, ymin, zmin;
    double *density;
    double *potential;
    double *field_x, *field_y, *field_z;
} Grid3D;

Grid3D* create_grid(int nx, int ny, int nz, double box_size);
void free_grid(Grid3D *grid);
void particles_to_grid(ParticleSystem *sys, Grid3D *grid);
void solve_poisson_fft(Grid3D *grid);
void compute_electric_field(Grid3D *grid);
void grid_to_particles(ParticleSystem *sys, Grid3D *grid);

#endif
