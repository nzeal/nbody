// ===== gravity.c =====
#include "gravity.h"
#include <math.h>

void compute_gravity_direct(ParticleSystem *sys, double softening) {
    double eps2 = softening * softening;
    
    #pragma acc parallel loop present(sys->particles[0:sys->n_particles])
    for (int i = 0; i < sys->n_particles; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        
        #pragma acc loop seq
        for (int j = 0; j < sys->n_particles; j++) {
            if (i != j) {
                double dx = sys->particles[j].x - sys->particles[i].x;
                double dy = sys->particles[j].y - sys->particles[i].y;
                double dz = sys->particles[j].z - sys->particles[i].z;
                
                double r2 = dx*dx + dy*dy + dz*dz + eps2;
                double r = sqrt(r2);
                double inv_r3 = 1.0 / (r * r2);
                
                double force = sys->particles[j].mass * inv_r3;
                
                fx += force * dx;
                fy += force * dy;
                fz += force * dz;
            }
        }
        
        sys->particles[i].fx = fx;
        sys->particles[i].fy = fy;
        sys->particles[i].fz = fz;
    }
}

void compute_gravity_tree(ParticleSystem *sys, BHTree *tree, double softening) {
    reset_forces(sys);
    
    // Note: Tree traversal is difficult to parallelize efficiently on GPU
    // For now, we'll do this on CPU and copy results to GPU
    #pragma acc update host(sys->particles[0:sys->n_particles])
    
    for (int i = 0; i < sys->n_particles; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        tree_force_calculation(tree->root, i, sys, tree, softening, &fx, &fy, &fz);
        sys->particles[i].fx = fx;
        sys->particles[i].fy = fy;
        sys->particles[i].fz = fz;
    }
    
    #pragma acc update device(sys->particles[0:sys->n_particles])
}

void tree_force_calculation(TreeNode *node, int particle_idx, ParticleSystem *sys,
                           BHTree *tree, double softening, double *fx, double *fy, double *fz) {
    if (!node || node->total_mass == 0.0) return;
    
    double dx = node->cm_x - sys->particles[particle_idx].x;
    double dy = node->cm_y - sys->particles[particle_idx].y;
    double dz = node->cm_z - sys->particles[particle_idx].z;
    
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);
    
    if (node->is_leaf) {
        // Direct interaction with all particles in leaf
        for (int i = 0; i < node->n_particles; i++) {
            int idx = node->particle_indices[i];
            if (idx != particle_idx) {
                dx = sys->particles[idx].x - sys->particles[particle_idx].x;
                dy = sys->particles[idx].y - sys->particles[particle_idx].y;
                dz = sys->particles[idx].z - sys->particles[particle_idx].z;
                
                r2 = dx*dx + dy*dy + dz*dz + softening*softening;
                r = sqrt(r2);
                double inv_r3 = 1.0 / (r * r2);
                
                double force = sys->particles[idx].mass * inv_r3;
                
                *fx += force * dx;
                *fy += force * dy;
                *fz += force * dz;
            }
        }
    } else {
        // Check Barnes-Hut criterion
        double node_size = fmax(fmax(node->xmax - node->xmin, 
                                    node->ymax - node->ymin),
                               node->zmax - node->zmin);
        
        if (node_size / r < tree->theta) {
            // Use multipole approximation
            r2 += softening * softening;
            r = sqrt(r2);
            double inv_r3 = 1.0 / (r * r2);
            
            double force = node->total_mass * inv_r3;
            
            *fx += force * dx;
            *fy += force * dy;
            *fz += force * dz;
        } else {
            // Recurse to children
            for (int i = 0; i < 8; i++) {
                if (node->children[i]) {
                    tree_force_calculation(node->children[i], particle_idx, sys, 
                                         tree, softening, fx, fy, fz);
                }
            }
        }
    }
}
