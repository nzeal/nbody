// ===== gravity.h =====
#ifndef GRAVITY_H
#define GRAVITY_H

#include "particle.h"
#include "tree.h"

void compute_gravity_direct(ParticleSystem *sys, double softening);
void compute_gravity_tree(ParticleSystem *sys, BHTree *tree, double softening);
void tree_force_calculation(TreeNode *node, int particle_idx, ParticleSystem *sys, 
                           BHTree *tree, double softening, double *fx, double *fy, double *fz);

#endif
