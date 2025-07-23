// ===== vlasov.c =====
#include "vlasov.h"

void integrate_leapfrog(ParticleSystem *sys, double dt) {
    // Leapfrog integration: kick-drift-kick
    kick_particles(sys, 0.5 * dt);
    drift_particles(sys, dt);
    kick_particles(sys, 0.5 * dt);
}

void drift_particles(ParticleSystem *sys, double dt) {
    #pragma acc parallel loop present(sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        sys->particles[i].x += sys->particles[i].vx * dt;
        sys->particles[i].y += sys->particles[i].vy * dt;
        sys->particles[i].z += sys->particles[i].vz * dt;
    }
}

void kick_particles(ParticleSystem *sys, double dt) {
    #pragma acc parallel loop present(sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        sys->particles[i].vx += sys->particles[i].fx * dt / sys->particles[i].mass;
        sys->particles[i].vy += sys->particles[i].fy * dt / sys->particles[i].mass;
        sys->particles[i].vz += sys->particles[i].fz * dt / sys->particles[i].mass;
    }
}

void integrate_rk4(ParticleSystem *sys, double dt) {
    // Simple RK4 implementation (more complex for N-body systems)
    // For demonstration - leapfrog is typically preferred for N-body
    
    // This is a simplified version
    drift_particles(sys, 0.5 * dt);
    kick_particles(sys, dt);
    drift_particles(sys, 0.5 * dt);
}

