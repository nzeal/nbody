// ===== particle.h =====
#ifndef PARTICLE_H
#define PARTICLE_H

typedef struct {
    double x, y, z;     // Position
    double vx, vy, vz;  // Velocity
    double mass;        // Mass
    double fx, fy, fz;  // Force
} Particle;

typedef struct {
    int n_particles;
    Particle *particles;
    double box_size;
    double total_mass;
} ParticleSystem;

// Function declarations
ParticleSystem* init_particle_system(int n_particles, double box_size);
void free_particle_system(ParticleSystem *sys);
void init_random_particles(ParticleSystem *sys, unsigned int seed);
void init_plummer_sphere(ParticleSystem *sys, double scale_radius);
void reset_forces(ParticleSystem *sys);

#endif
