#ifndef IO_H
#define IO_H

#include "particle.h" // Include particle.h for Particle and ParticleSystem

// Function declarations
void write_particles_txt(ParticleSystem *sys, const char *outdir, const char *filename, int step);
void write_particles_binary(ParticleSystem *sys, const char *outdir, const char *filename, int step);
void write_energy_log(ParticleSystem *sys, const char *outdir, double time, double kinetic, double potential, const char *filename);
double compute_kinetic_energy(ParticleSystem *sys);
double compute_potential_energy(ParticleSystem *sys, double softening);
void print_system_stats(ParticleSystem *sys, int step, double time);

#endif
