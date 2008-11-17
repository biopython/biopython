#include "Neighbor.h"

struct KDTree;

struct KDTree* KDTree_init(int dim, int bucket_size);
void KDTree_destroy(struct KDTree* tree);
int KDTree_set_data(struct KDTree* tree, float *coords, long int nr_points);
long int KDTree_get_count(struct KDTree* tree);
long int KDTree_neighbor_get_count(struct KDTree* tree);
int KDTree_search_center_radius(struct KDTree* tree, float *coord, float radius);
void KDTree_copy_indices(struct KDTree* tree, long *indices);
void KDTree_copy_radii(struct KDTree* tree, float *radii);
int KDTree_neighbor_search(struct KDTree* tree, float neighbor_radius, struct Neighbor** neighbors);
int KDTree_neighbor_simple_search(struct KDTree* tree, float radius, struct Neighbor** neighbors);
