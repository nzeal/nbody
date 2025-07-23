#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strcpy and strcmp
#include <math.h>
#include <time.h>
#include "particle.h"
#include "tree.h"
#include "gravity.h"
#include "vlasov.h"
#include "poisson.h"
#include "io.h"

typedef struct {
    int n_particles;
    double box_size;
    double dt;
    double t_end;
    double softening;
    double theta;  // Barnes-Hut opening angle
    int output_freq;
    int use_tree;
    int use_grid;
    char output_prefix[256];
    char output_dir[256]; // Output directory
} SimulationParams;

void print_usage(const char *program_name) {
    printf("Usage: %s [options]\n", program_name);
    printf("Options:\n");
    printf("  -n <int>     Number of particles (default: 1024)\n");
    printf("  -L <double>  Box size (default: 10.0)\n");
    printf("  -dt <double> Time step (default: 0.01)\n");
    printf("  -T <double>  End time (default: 10.0)\n");
    printf("  -eps <double> Softening parameter (default: 0.01)\n");
    printf("  -theta <double> Barnes-Hut opening angle (default: 0.5)\n");
    printf("  -freq <int>  Output frequency (default: 10)\n");
    printf("  -tree        Use Barnes-Hut tree (default: direct)\n");
    printf("  -grid        Use grid-based Poisson solver\n");
    printf("  -o <string>  Output prefix (default: 'particles')\n");
    printf("  -outdir <string> Output directory (default: './data')\n");
    printf("  -h           Show this help\n");
}

void parse_arguments(int argc, char *argv[], SimulationParams *params) {
    // Set defaults
    params->n_particles = 10000;
    params->box_size = 100.0;
    params->dt = 0.01;
    params->t_end = 10.0;
    params->softening = 0.01;
    params->theta = 0.5;
    params->output_freq = 10;
    params->use_tree = 0;
    params->use_grid = 0;
    strcpy(params->output_prefix, "particles");
    strcpy(params->output_dir, "./data"); // Default output directory

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            params->n_particles = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-L") == 0 && i + 1 < argc) {
            params->box_size = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dt") == 0 && i + 1 < argc) {
            params->dt = atof(argv[++i]);
        } else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc) {
            params->t_end = atof(argv[++i]);
        } else if (strcmp(argv[i], "-eps") == 0 && i + 1 < argc) {
            params->softening = atof(argv[++i]);
        } else if (strcmp(argv[i], "-theta") == 0 && i + 1 < argc) {
            params->theta = atof(argv[++i]);
        } else if (strcmp(argv[i], "-freq") == 0 && i + 1 < argc) {
            params->output_freq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-tree") == 0) {
            params->use_tree = 1;
        } else if (strcmp(argv[i], "-grid") == 0) {
            params->use_grid = 1;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            strcpy(params->output_prefix, argv[++i]);
        } else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc) {
            strcpy(params->output_dir, argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            exit(0);
        }
    }
}

int main(int argc, char *argv[]) {
    SimulationParams params;
    parse_arguments(argc, argv, &params); // Fixed typo: ¶ms to params

    printf("=== Vlasov-Poisson N-body Simulation ===\n");
    printf("Particles: %d\n", params.n_particles);
    printf("Box size: %.2f\n", params.box_size);
    printf("Time step: %.4f\n", params.dt);
    printf("End time: %.2f\n", params.t_end);
    printf("Softening: %.4f\n", params.softening);
    printf("Method: %s\n", params.use_tree ? "Barnes-Hut Tree" : "Direct N-body");
    if (params.use_tree) {
        printf("Opening angle: %.2f\n", params.theta);
    }
    printf("Grid Poisson solver: %s\n", params.use_grid ? "Yes" : "No");
    printf("Output directory: %s\n", params.output_dir);
    printf("\n");

    // Initialize particle system
    ParticleSystem *sys = init_particle_system(params.n_particles, params.box_size);

    // Initialize with Plummer sphere
    init_plummer_sphere(sys, 1.0);

    // Initialize Barnes-Hut tree if needed
    BHTree *tree = NULL;
    if (params.use_tree) {
        tree = create_tree(params.theta);
    }

    // Initialize grid if needed
    Grid3D *grid = NULL;
    if (params.use_grid) {
        int grid_size = (int)cbrt(params.n_particles / 8);  // Rough heuristic
        grid_size = fmax(32, fmin(128, grid_size));
        grid = create_grid(grid_size, grid_size, grid_size, params.box_size);
        printf("Grid size: %d^3\n", grid_size);
    }

    // Main simulation loop
    double time = 0.0;
    int step = 0;
    int n_steps = (int)(params.t_end / params.dt);

    clock_t start_time = clock();

    // Initial output
    write_particles_txt(sys, params.output_dir, params.output_prefix, 0);
    print_system_stats(sys, 0, 0.0);

    printf("\nStarting simulation...\n");

    for (step = 1; step <= n_steps; step++) {
        // Reset forces
        reset_forces(sys);

        // Compute gravitational forces
        if (params.use_tree) {
            build_tree(tree, sys);
            compute_gravity_tree(sys, tree, params.softening);
        } else {
            compute_gravity_direct(sys, params.softening);
        }

        // Add grid-based forces if using Poisson solver
        if (params.use_grid) {
            particles_to_grid(sys, grid);
            solve_poisson_fft(grid);
            compute_electric_field(grid);
            grid_to_particles(sys, grid);
        }

        // Integrate equations of motion
        integrate_leapfrog(sys, params.dt);

        time += params.dt;

        // Output
        if (step % params.output_freq == 0) {
            write_particles_txt(sys, params.output_dir, params.output_prefix, step);
            print_system_stats(sys, step, time);

            // Write energy log (fixed: include time parameter)
            double ke = compute_kinetic_energy(sys);
            double pe = compute_potential_energy(sys, params.softening);
            write_energy_log(sys, params.output_dir, time, ke, pe, "energy_log.txt");
        }

        // Progress indicator
        if (step % (n_steps / 10) == 0) {
            printf("Progress: %d%% (Step %d/%d)\n",
                   (int)(100.0 * step / n_steps), step, n_steps);
        }
    }

    clock_t end_time = clock();
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nSimulation completed!\n");
    printf("Total time: %.2f seconds\n", cpu_time);
    printf("Time per step: %.4f seconds\n", cpu_time / n_steps);
    printf("Particle updates per second: %.0f\n",
           (double)params.n_particles * n_steps / cpu_time);

    // Final output
    write_particles_txt(sys, params.output_dir, params.output_prefix, step);
    write_particles_binary(sys, params.output_dir, params.output_prefix, step);

    // Cleanup
    free_particle_system(sys);
    if (tree) free_tree(tree);
    if (grid) free_grid(grid);

    return 0;
}
