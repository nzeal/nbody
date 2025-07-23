// ===== tree.h =====
#ifndef TREE_H
#define TREE_H

#include "particle.h"

#define MAX_PARTICLES_PER_NODE 1
#define MAX_TREE_DEPTH 20

typedef struct TreeNode {
    // Bounding box
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    // Center of mass and total mass
    double cm_x, cm_y, cm_z;
    double total_mass;
    
    // Tree structure
    struct TreeNode *children[8];  // Octree
    int is_leaf;
    int n_particles;
    int *particle_indices;
    
    // Tree level
    int level;
} TreeNode;

typedef struct {
    TreeNode *root;
    int max_particles;
    int total_nodes;
    double theta;  // Barnes-Hut opening angle criterion
} BHTree;

// Function declarations
BHTree* create_tree(double theta);
void free_tree(BHTree *tree);
TreeNode* create_node(double xmin, double xmax, double ymin, double ymax, 
                     double zmin, double zmax, int level);
void free_node(TreeNode *node);
void build_tree(BHTree *tree, ParticleSystem *sys);
void insert_particle(TreeNode *node, int particle_idx, ParticleSystem *sys);
void compute_center_of_mass(TreeNode *node, ParticleSystem *sys);
void subdivide_node(TreeNode *node);

#endif

