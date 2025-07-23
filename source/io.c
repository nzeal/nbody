#include "io.h"
#include <stdio.h>
#include <math.h>
#include <sys/stat.h> // For mkdir

// Ensure the output directory exists
void ensure_output_directory(const char *outdir) {
    #ifdef _WIN32
        mkdir(outdir);
    #else
        mkdir(outdir, 0755); // Create directory with permissions in Unix-like systems
    #endif
}

void write_particles_txt(ParticleSystem *sys, const char *outdir, const char *filename, int step) {
    char full_filename[256];
    ensure_output_directory(outdir);
    snprintf(full_filename, sizeof(full_filename), "%s/%s_%04d.txt", outdir, filename, step);

    #pragma acc update host(sys->particles[0:sys->n_particles])

    FILE *file = fopen(full_filename, "w");
    if (!file) {
        printf("Error: Could not open file %s\n", full_filename);
        return;
    }

    fprintf(file, "x,y,z,vx,vy,vz,mass,fx,fy,fz\n");

    for (int i = 0; i < sys->n_particles; i++) {
        fprintf(file, "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                sys->particles[i].x, sys->particles[i].y, sys->particles[i].z,
                sys->particles[i].vx, sys->particles[i].vy, sys->particles[i].vz,
                sys->particles[i].mass,
                sys->particles[i].fx, sys->particles[i].fy, sys->particles[i].fz);
    }

    fclose(file);
    printf("Wrote particle data to %s\n", full_filename);
}

void write_particles_binary(ParticleSystem *sys, const char *outdir, const char *filename, int step) {
    char full_filename[256];
    ensure_output_directory(outdir);
    snprintf(full_filename, sizeof(full_filename), "%s/%s_%04d.bin", outdir, filename, step);

    #pragma acc update host(sys->particles[0:sys->n_particles])

    FILE *file = fopen(full_filename, "wb");
    if (!file) {
        printf("Error: Could not open file %s\n", full_filename);
        return;
    }

    fwrite(&sys->n_particles, sizeof(int), 1, file);
    fwrite(sys->particles, sizeof(Particle), sys->n_particles, file);

    fclose(file);
    printf("Wrote binary particle data to %s\n", full_filename);
}

void write_energy_log(ParticleSystem *sys, const char *outdir, double time, double kinetic,
                     double potential, const char *filename) {
    static int first_call = 1;
    char full_filename[256];
    ensure_output_directory(outdir);
    snprintf(full_filename, sizeof(full_filename), "%s/%s", outdir, filename);

    FILE *file = fopen(full_filename, first_call ? "w" : "a");
    if (!file) {
        printf("Error: Could not open energy log file %s\n", full_filename);
        return;
    }

    if (first_call) {
        fprintf(file, "time,kinetic_energy,potential_energy,total_energy\n");
        first_call = 0;
    }

    fprintf(file, "%.6e,%.6e,%.6e,%.6e\n", time, kinetic, potential, kinetic + potential);
    fclose(file);
    printf("Wrote energy log to %s\n", full_filename);
}

double compute_kinetic_energy(ParticleSystem *sys) {
    #pragma acc update host(sys->particles[0:sys->n_particles])

    double kinetic = 0.0;
    for (int i = 0; i < sys->n_particles; i++) {
        double v2 = sys->particles[i].vx * sys->particles[i].vx +
                   sys->particles[i].vy * sys->particles[i].vy +
                   sys->particles[i].vz * sys->particles[i].vz;
        kinetic += 0.5 * sys->particles[i].mass * v2;
    }
    return kinetic;
}

double compute_potential_energy(ParticleSystem *sys, double softening) {
    #pragma acc update host(sys->particles[0:sys->n_particles])

    double potential = 0.0;
    double eps2 = softening * softening;

    for (int i = 0; i < sys->n_particles; i++) {
        for (int j = i + 1; j < sys->n_particles; j++) {
            double dx = sys->particles[j].x - sys->particles[i].x;
            double dy = sys->particles[j].y - sys->particles[i].y;
            double dz = sys->particles[j].z - sys->particles[i].z;

            double r2 = dx*dx + dy*dy + dz*dz + eps2;
            double r = sqrt(r2);

            potential -= sys->particles[i].mass * sys->particles[j].mass / r;
        }
    }
    return potential;
}

void print_system_stats(ParticleSystem *sys, int step, double time) {
    double kinetic = compute_kinetic_energy(sys);
    double potential = compute_potential_energy(sys, 0.01);

    printf("Step %d, Time: %.3f, KE: %.6e, PE: %.6e, Total: %.6e\n",
           step, time, kinetic, potential, kinetic + potential);
}
