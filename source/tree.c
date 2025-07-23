// ===== tree.c =====
#include "tree.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

BHTree* create_tree(double theta) {
    BHTree *tree = malloc(sizeof(BHTree));
    tree->root = NULL;
    tree->max_particles = 0;
    tree->total_nodes = 0;
    tree->theta = theta;
    return tree;
}

void free_tree(BHTree *tree) {
    if (tree) {
        if (tree->root) {
            free_node(tree->root);
        }
        free(tree);
    }
}

TreeNode* create_node(double xmin, double xmax, double ymin, double ymax,
                     double zmin, double zmax, int level) {
    TreeNode *node = malloc(sizeof(TreeNode));
    
    node->xmin = xmin; node->xmax = xmax;
    node->ymin = ymin; node->ymax = ymax;
    node->zmin = zmin; node->zmax = zmax;
    
    node->cm_x = node->cm_y = node->cm_z = 0.0;
    node->total_mass = 0.0;
    
    for (int i = 0; i < 8; i++) {
        node->children[i] = NULL;
    }
    
    node->is_leaf = 1;
    node->n_particles = 0;
    node->particle_indices = NULL;
    node->level = level;
    
    return node;
}

void free_node(TreeNode *node) {
    if (node) {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                free_node(node->children[i]);
            }
        }
        if (node->particle_indices) {
            free(node->particle_indices);
        }
        free(node);
    }
}

void build_tree(BHTree *tree, ParticleSystem *sys) {
    // Free existing tree
    if (tree->root) {
        free_node(tree->root);
    }
    
    // Find bounding box
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = xmax = sys->particles[0].x;
    ymin = ymax = sys->particles[0].y;
    zmin = zmax = sys->particles[0].z;
    
    for (int i = 1; i < sys->n_particles; i++) {
        if (sys->particles[i].x < xmin) xmin = sys->particles[i].x;
        if (sys->particles[i].x > xmax) xmax = sys->particles[i].x;
        if (sys->particles[i].y < ymin) ymin = sys->particles[i].y;
        if (sys->particles[i].y > ymax) ymax = sys->particles[i].y;
        if (sys->particles[i].z < zmin) zmin = sys->particles[i].z;
        if (sys->particles[i].z > zmax) zmax = sys->particles[i].z;
    }
    
    // Expand box slightly to ensure all particles are inside
    double margin = 0.01 * fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin);
    xmin -= margin; xmax += margin;
    ymin -= margin; ymax += margin;
    zmin -= margin; zmax += margin;
    
    // Create root node
    tree->root = create_node(xmin, xmax, ymin, ymax, zmin, zmax, 0);
    tree->total_nodes = 1;
    
    // Insert all particles
    for (int i = 0; i < sys->n_particles; i++) {
        insert_particle(tree->root, i, sys);
    }
    
    // Compute center of mass for all nodes
    compute_center_of_mass(tree->root, sys);
}

void insert_particle(TreeNode *node, int particle_idx, ParticleSystem *sys) {
    if (node->is_leaf) {
        if (node->n_particles == 0) {
            // First particle in this node
            node->particle_indices = malloc(sizeof(int));
            node->particle_indices[0] = particle_idx;
            node->n_particles = 1;
        } else if (node->n_particles < MAX_PARTICLES_PER_NODE || 
                   node->level >= MAX_TREE_DEPTH) {
            // Add particle to existing leaf
            node->particle_indices = realloc(node->particle_indices, 
                                           (node->n_particles + 1) * sizeof(int));
            node->particle_indices[node->n_particles] = particle_idx;
            node->n_particles++;
        } else {
            // Need to subdivide
            subdivide_node(node);
            
            // Re-insert existing particles
            for (int i = 0; i < node->n_particles; i++) {
                int old_idx = node->particle_indices[i];
                
                // Determine which child this particle belongs to
                int child_idx = 0;
                double xmid = 0.5 * (node->xmin + node->xmax);
                double ymid = 0.5 * (node->ymin + node->ymax);
                double zmid = 0.5 * (node->zmin + node->zmax);
                
                if (sys->particles[old_idx].x > xmid) child_idx |= 1;
                if (sys->particles[old_idx].y > ymid) child_idx |= 2;
                if (sys->particles[old_idx].z > zmid) child_idx |= 4;
                
                insert_particle(node->children[child_idx], old_idx, sys);
            }
            
            // Insert new particle
            int child_idx = 0;
            double xmid = 0.5 * (node->xmin + node->xmax);
            double ymid = 0.5 * (node->ymin + node->ymax);
            double zmid = 0.5 * (node->zmin + node->zmax);
            
            if (sys->particles[particle_idx].x > xmid) child_idx |= 1;
            if (sys->particles[particle_idx].y > ymid) child_idx |= 2;
            if (sys->particles[particle_idx].z > zmid) child_idx |= 4;
            
            insert_particle(node->children[child_idx], particle_idx, sys);
            
            // Clean up leaf data
            free(node->particle_indices);
            node->particle_indices = NULL;
            node->n_particles = 0;
            node->is_leaf = 0;
        }
    } else {
        // Internal node - find correct child
        int child_idx = 0;
        double xmid = 0.5 * (node->xmin + node->xmax);
        double ymid = 0.5 * (node->ymin + node->ymax);
        double zmid = 0.5 * (node->zmin + node->zmax);
        
        if (sys->particles[particle_idx].x > xmid) child_idx |= 1;
        if (sys->particles[particle_idx].y > ymid) child_idx |= 2;
        if (sys->particles[particle_idx].z > zmid) child_idx |= 4;
        
        insert_particle(node->children[child_idx], particle_idx, sys);
    }
}

void subdivide_node(TreeNode *node) {
    double xmid = 0.5 * (node->xmin + node->xmax);
    double ymid = 0.5 * (node->ymin + node->ymax);
    double zmid = 0.5 * (node->zmin + node->zmax);
    
    // Create 8 children (octree)
    node->children[0] = create_node(node->xmin, xmid, node->ymin, ymid, 
                                   node->zmin, zmid, node->level + 1);
    node->children[1] = create_node(xmid, node->xmax, node->ymin, ymid, 
                                   node->zmin, zmid, node->level + 1);
    node->children[2] = create_node(node->xmin, xmid, ymid, node->ymax, 
                                   node->zmin, zmid, node->level + 1);
    node->children[3] = create_node(xmid, node->xmax, ymid, node->ymax, 
                                   node->zmin, zmid, node->level + 1);
    node->children[4] = create_node(node->xmin, xmid, node->ymin, ymid, 
                                   zmid, node->zmax, node->level + 1);
    node->children[5] = create_node(xmid, node->xmax, node->ymin, ymid, 
                                   zmid, node->zmax, node->level + 1);
    node->children[6] = create_node(node->xmin, xmid, ymid, node->ymax, 
                                   zmid, node->zmax, node->level + 1);
    node->children[7] = create_node(xmid, node->xmax, ymid, node->ymax, 
                                   zmid, node->zmax, node->level + 1);
}

void compute_center_of_mass(TreeNode *node, ParticleSystem *sys) {
    if (node->is_leaf) {
        // Compute center of mass for leaf node
        node->total_mass = 0.0;
        node->cm_x = node->cm_y = node->cm_z = 0.0;
        
        for (int i = 0; i < node->n_particles; i++) {
            int idx = node->particle_indices[i];
            double mass = sys->particles[idx].mass;
            
            node->total_mass += mass;
            node->cm_x += mass * sys->particles[idx].x;
            node->cm_y += mass * sys->particles[idx].y;
            node->cm_z += mass * sys->particles[idx].z;
        }
        
        if (node->total_mass > 0.0) {
            node->cm_x /= node->total_mass;
            node->cm_y /= node->total_mass;
            node->cm_z /= node->total_mass;
        }
    } else {
        // Compute center of mass for internal node
        node->total_mass = 0.0;
        node->cm_x = node->cm_y = node->cm_z = 0.0;
        
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                compute_center_of_mass(node->children[i], sys);
                
                double child_mass = node->children[i]->total_mass;
                if (child_mass > 0.0) {
                    node->total_mass += child_mass;
                    node->cm_x += child_mass * node->children[i]->cm_x;
                    node->cm_y += child_mass * node->children[i]->cm_y;
                    node->cm_z += child_mass * node->children[i]->cm_z;
                }
            }
        }
        
        if (node->total_mass > 0.0) {
            node->cm_x /= node->total_mass;
            node->cm_y /= node->total_mass;
            node->cm_z /= node->total_mass;
        }
    }
}

