// ===== vlasov.h =====
#ifndef VLASOV_H
#define VLASOV_H

#include "particle.h"

void integrate_leapfrog(ParticleSystem *sys, double dt);
void integrate_rk4(ParticleSystem *sys, double dt);
void drift_particles(ParticleSystem *sys, double dt);
void kick_particles(ParticleSystem *sys, double dt);

#endif
