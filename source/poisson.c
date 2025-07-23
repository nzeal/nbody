// ===== poisson.c =====
#include "poisson.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

Grid3D* create_grid(int nx, int ny, int nz, double box_size) {
    Grid3D *grid = malloc(sizeof(Grid3D));
    
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    
    grid->dx = box_size / nx;
    grid->dy = box_size / ny;
    grid->dz = box_size / nz;
    
    grid->xmin = grid->ymin = grid->zmin = -box_size / 2.0;
    
    int total_cells = nx * ny * nz;
    
    grid->density = calloc(total_cells, sizeof(double));
    grid->potential = calloc(total_cells, sizeof(double));
    grid->field_x = calloc(total_cells, sizeof(double));
    grid->field_y = calloc(total_cells, sizeof(double));
    grid->field_z = calloc(total_cells, sizeof(double));
    
    return grid;
}

void free_grid(Grid3D *grid) {
    if (grid) {
        free(grid->density);
        free(grid->potential);
        free(grid->field_x);
        free(grid->field_y);
        free(grid->field_z);
        free(grid);
    }
}

void particles_to_grid(ParticleSystem *sys, Grid3D *grid) {
    // Clear density grid
    memset(grid->density, 0, grid->nx * grid->ny * grid->nz * sizeof(double));
    
    // Cloud-in-cell (CIC) assignment
    for (int p = 0; p < sys->n_particles; p++) {
        double x = sys->particles[p].x;
        double y = sys->particles[p].y;
        double z = sys->particles[p].z;
        double mass = sys->particles[p].mass;
        
        // Find grid position
        double fx = (x - grid->xmin) / grid->dx - 0.5;
        double fy = (y - grid->ymin) / grid->dy - 0.5;
        double fz = (z - grid->zmin) / grid->dz - 0.5;
        
        int ix = (int)floor(fx);
        int iy = (int)floor(fy);
        int iz = (int)floor(fz);
        
        double dx = fx - ix;
        double dy = fy - iy;
        double dz = fz - iz;
        
        // CIC weights
        double wx[2] = {1.0 - dx, dx};
        double wy[2] = {1.0 - dy, dy};
        double wz[2] = {1.0 - dz, dz};
        
        // Deposit mass to 8 neighboring cells
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    int gx = ix + i;
                    int gy = iy + j;
                    int gz = iz + k;
                    
                    // Periodic boundary conditions
                    gx = (gx + grid->nx) % grid->nx;
                    gy = (gy + grid->ny) % grid->ny;
                    gz = (gz + grid->nz) % grid->nz;
                    
                    int idx = gz * grid->nx * grid->ny + gy * grid->nx + gx;
                    grid->density[idx] += mass * wx[i] * wy[j] * wz[k];
                }
            }
        }
    }
}

void solve_poisson_fft(Grid3D *grid) {
    // Simplified Poisson solver using finite differences
    // In practice, you would use FFT for periodic boundaries
    
    double dx2 = grid->dx * grid->dx;
    double dy2 = grid->dy * grid->dy;
    double dz2 = grid->dz * grid->dz;
    
    // Simple Jacobi iteration (not optimal, but works)
    double *phi_new = malloc(grid->nx * grid->ny * grid->nz * sizeof(double));
    
    for (int iter = 0; iter < 100; iter++) {
        for (int i = 1; i < grid->nx - 1; i++) {
            for (int j = 1; j < grid->ny - 1; j++) {
                for (int k = 1; k < grid->nz - 1; k++) {
                    int idx = k * grid->nx * grid->ny + j * grid->nx + i;
                    
                    int idx_xp = k * grid->nx * grid->ny + j * grid->nx + (i+1);
                    int idx_xm = k * grid->nx * grid->ny + j * grid->nx + (i-1);
                    int idx_yp = k * grid->nx * grid->ny + (j+1) * grid->nx + i;
                    int idx_ym = k * grid->nx * grid->ny + (j-1) * grid->nx + i;
                    int idx_zp = (k+1) * grid->nx * grid->ny + j * grid->nx + i;
                    int idx_zm = (k-1) * grid->nx * grid->ny + j * grid->nx + i;
                    
                    phi_new[idx] = (grid->potential[idx_xp] + grid->potential[idx_xm]) / dx2 +
                                   (grid->potential[idx_yp] + grid->potential[idx_ym]) / dy2 +
                                   (grid->potential[idx_zp] + grid->potential[idx_zm]) / dz2 -
                                   4.0 * M_PI * grid->density[idx];
                    
                    phi_new[idx] /= (2.0/dx2 + 2.0/dy2 + 2.0/dz2);
                }
            }
        }
        
        // Copy back
        memcpy(grid->potential, phi_new, grid->nx * grid->ny * grid->nz * sizeof(double));
    }
    
    free(phi_new);
}

void compute_electric_field(Grid3D *grid) {
    // Compute E = -grad(phi) using finite differences
    for (int i = 1; i < grid->nx - 1; i++) {
        for (int j = 1; j < grid->ny - 1; j++) {
            for (int k = 1; k < grid->nz - 1; k++) {
                int idx = k * grid->nx * grid->ny + j * grid->nx + i;
                
                int idx_xp = k * grid->nx * grid->ny + j * grid->nx + (i+1);
                int idx_xm = k * grid->nx * grid->ny + j * grid->nx + (i-1);
                int idx_yp = k * grid->nx * grid->ny + (j+1) * grid->nx + i;
                int idx_ym = k * grid->nx * grid->ny + (j-1) * grid->nx + i;
                int idx_zp = (k+1) * grid->nx * grid->ny + j * grid->nx + i;
                int idx_zm = (k-1) * grid->nx * grid->ny + j * grid->nx + i;
                
                grid->field_x[idx] = -(grid->potential[idx_xp] - grid->potential[idx_xm]) / (2.0 * grid->dx);
                grid->field_y[idx] = -(grid->potential[idx_yp] - grid->potential[idx_ym]) / (2.0 * grid->dy);
                grid->field_z[idx] = -(grid->potential[idx_zp] - grid->potential[idx_zm]) / (2.0 * grid->dz);
            }
        }
    }
}

void grid_to_particles(ParticleSystem *sys, Grid3D *grid) {
    // Interpolate electric field from grid to particles
    for (int p = 0; p < sys->n_particles; p++) {
        double x = sys->particles[p].x;
        double y = sys->particles[p].y;
        double z = sys->particles[p].z;
        
        // Find grid position
        double fx = (x - grid->xmin) / grid->dx - 0.5;
        double fy = (y - grid->ymin) / grid->dy - 0.5;
        double fz = (z - grid->zmin) / grid->dz - 0.5;
        
        int ix = (int)floor(fx);
        int iy = (int)floor(fy);
        int iz = (int)floor(fz);
        
        double dx = fx - ix;
        double dy = fy - iy;
        double dz = fz - iz;
        
        // CIC weights
        double wx[2] = {1.0 - dx, dx};
        double wy[2] = {1.0 - dy, dy};
        double wz[2] = {1.0 - dz, dz};
        
        double ex = 0.0, ey = 0.0, ez = 0.0;
        
        // Interpolate from 8 neighboring cells
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    int gx = ix + i;
                    int gy = iy + j;
                    int gz = iz + k;
                    
                    // Periodic boundary conditions
                    gx = (gx + grid->nx) % grid->nx;
                    gy = (gy + grid->ny) % grid->ny;
                    gz = (gz + grid->nz) % grid->nz;
                    
                    int idx = gz * grid->nx * grid->ny + gy * grid->nx + gx;
                    double weight = wx[i] * wy[j] * wz[k];
                    
                    ex += weight * grid->field_x[idx];
                    ey += weight * grid->field_y[idx];
                    ez += weight * grid->field_z[idx];
                }
            }
        }
        
        // Add electric force (assuming unit charge)
        sys->particles[p].fx += ex * sys->particles[p].mass;
        sys->particles[p].fy += ey * sys->particles[p].mass;
        sys->particles[p].fz += ez * sys->particles[p].mass;
    }
}

