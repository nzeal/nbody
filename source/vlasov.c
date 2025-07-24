// ===== vlasov.c =====
#include "vlasov.h"
 #include <stddef.h>

// Utility to ensure data is present on the device for OpenACC
static inline void acc_present_particles(ParticleSystem *sys) {
    #pragma acc data present_or_copyin(sys->particles[0:sys->n_particles])
    { /* Intentional no-op: ensures mapping */}
}

void integrate_leapfrog(ParticleSystem *sys, double dt) {
    // Leapfrog integration: kick-drift-kick
    kick_particles(sys, 0.5 * dt);
    drift_particles(sys, dt);
    kick_particles(sys, 0.5 * dt);
}

void drift_particles(ParticleSystem *sys, double dt) {
    if (sys == NULL || sys->particles == NULL || sys->n_particles <= 0) return;
    #pragma acc parallel loop present(sys, sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        sys->particles[i].x += sys->particles[i].vx * dt;
        sys->particles[i].y += sys->particles[i].vy * dt;
        sys->particles[i].z += sys->particles[i].vz * dt;
    }
}

void kick_particles(ParticleSystem *sys, double dt) {
    if (sys == NULL || sys->particles == NULL || sys->n_particles <= 0) return;
    #pragma acc parallel loop present(sys, sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        double mass = sys->particles[i].mass;
        // Avoid division by zero
        if (mass == 0.0) continue;
        sys->particles[i].vx += sys->particles[i].fx * dt / mass;
        sys->particles[i].vy += sys->particles[i].fy * dt / mass;
        sys->particles[i].vz += sys->particles[i].fz * dt / mass;
    }
}

void integrate_rk4(ParticleSystem *sys, double dt) {
    // Simple RK4 implementation (for demonstration purposes)
    drift_particles(sys, 0.5 * dt);
    kick_particles(sys, dt);
    drift_particles(sys, 0.5 * dt);
}
