// ===== particle.c =====
#include "particle.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* 
 * ParticleSystem* init_particle_system(int n_particles, double box_size) {
    ParticleSystem *sys = malloc(sizeof(ParticleSystem));
    sys->n_particles = n_particles;
    sys->box_size = box_size;
    sys->total_mass = 0.0;
    
    // Allocate particles on both host and device
    sys->particles = malloc(n_particles * sizeof(Particle));
    
    #pragma acc enter data create(sys->particles[0:n_particles])
    
    return sys;
}

*/

ParticleSystem* init_particle_system(int n_particles, double box_size) {
    ParticleSystem *sys = malloc(sizeof(ParticleSystem));
    sys->n_particles = n_particles;
    sys->box_size = box_size;
    sys->total_mass = 0.0;

    // Allocate particle array
    sys->particles = malloc(n_particles * sizeof(Particle));

    // Make both the system and the particle array visible to the device
    #pragma acc enter data copyin(sys[0:1])
    #pragma acc enter data create(sys->particles[0:n_particles])

    return sys;
}


void free_particle_system(ParticleSystem *sys) {
    if (sys) {
        #pragma acc exit data delete(sys->particles[0:sys->n_particles])
        free(sys->particles);
        free(sys);
    }
}

void init_random_particles(ParticleSystem *sys, unsigned int seed) {
    srand(seed);
    double total_mass = 0.0;
    
    for (int i = 0; i < sys->n_particles; i++) {
        // Random positions in box
        sys->particles[i].x = ((double)rand() / RAND_MAX - 0.5) * sys->box_size;
        sys->particles[i].y = ((double)rand() / RAND_MAX - 0.5) * sys->box_size;
        sys->particles[i].z = ((double)rand() / RAND_MAX - 0.5) * sys->box_size;
        
        // Random velocities (small)
        sys->particles[i].vx = ((double)rand() / RAND_MAX - 0.5) * 0.1;
        sys->particles[i].vy = ((double)rand() / RAND_MAX - 0.5) * 0.1;
        sys->particles[i].vz = ((double)rand() / RAND_MAX - 0.5) * 0.1;
        
        // Equal masses
        sys->particles[i].mass = 1.0 / sys->n_particles;
        total_mass += sys->particles[i].mass;
        
        // Initialize forces to zero
        sys->particles[i].fx = sys->particles[i].fy = sys->particles[i].fz = 0.0;
    }
    
    sys->total_mass = total_mass;
    
    // Copy to device
    #pragma acc update device(sys->particles[0:sys->n_particles])
}

void init_plummer_sphere(ParticleSystem *sys, double scale_radius) {
    srand(42);
    double total_mass = 0.0;
    
    for (int i = 0; i < sys->n_particles; i++) {
        // Plummer sphere density profile
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        double r3 = (double)rand() / RAND_MAX;
        
        // Rejection sampling for Plummer sphere
        double radius = scale_radius / sqrt(pow(r1, -2.0/3.0) - 1.0);
        
        // Random direction
        double theta = acos(2.0 * r2 - 1.0);  // cos(theta) uniform in [-1,1]
        double phi = 2.0 * M_PI * r3;
        
        sys->particles[i].x = radius * sin(theta) * cos(phi);
        sys->particles[i].y = radius * sin(theta) * sin(phi);
        sys->particles[i].z = radius * cos(theta);
        
        // Virial equilibrium velocities (simplified)
        double v_escape = sqrt(2.0 * sys->n_particles / radius);
        double v_mag = v_escape * 0.5 * ((double)rand() / RAND_MAX);
        
        theta = acos(2.0 * ((double)rand() / RAND_MAX) - 1.0);
        phi = 2.0 * M_PI * ((double)rand() / RAND_MAX);
        
        sys->particles[i].vx = v_mag * sin(theta) * cos(phi);
        sys->particles[i].vy = v_mag * sin(theta) * sin(phi);
        sys->particles[i].vz = v_mag * cos(theta);
        
        sys->particles[i].mass = 1.0 / sys->n_particles;
        total_mass += sys->particles[i].mass;
        
        sys->particles[i].fx = sys->particles[i].fy = sys->particles[i].fz = 0.0;
    }
    
    sys->total_mass = total_mass;
    
    #pragma acc update device(sys->particles[0:sys->n_particles])
}

void reset_forces(ParticleSystem *sys) {
    #pragma acc parallel loop present(sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        sys->particles[i].fx = 0.0;
        sys->particles[i].fy = 0.0;
        sys->particles[i].fz = 0.0;
    }
}

