#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
using namespace std;

#define INF 1000000
#undef NDEBUG


float KDTREE_dist(float *coord1, float *coord2, int dim);

class DataPoint
{
	private:
		long int _index;
		float *_coord;
	public:
		static int dim;
		// needed for vector & sort
		// T(); T(const T&); ~T(); T& operator=(const T&);
		friend int operator<(const DataPoint &self, const DataPoint &other);
		friend int operator==(const DataPoint &self, const DataPoint &other);
		// end needed for vector
		static int current_dim;
		void set_data(long int index, float *coord);
		float* get_coord(void);
		long int get_index(void);
};

class Node
{
	private:
		Node *_left;
		Node *_right;
		float _cut_value;
		int _cut_dim;
		long int _start, _end;
	public:
		Node(float cut_value, int cut_dim, long int start, long int end);
		~Node();
		void set_right_node(Node *node);
		void set_left_node(Node *node);
		Node *get_right_node(void);
		Node *get_left_node(void);
		long int get_index(void);
		float get_cut_value(void);
		int get_cut_dim(void);
		int is_leaf(void);
		int is_bucket(void);
		long int get_start(void);
		long int get_end(void);
};

class Region
{
	private:
		float *_left;
		float *_right;
	public:
		static int dim;
		Region(float *left=NULL, float *right=NULL);
		~Region();
		Region *intersect_right(float split_coord, int current_dim);
		Region *intersect_left(float split_coord, int current_dim);
		float *get_right(void);
		float *get_left(void);
		int encloses(float *coord);
		int test_intersection(Region *query_region, float radius=0);
		void print(void);
};

class KDTree
{
	private:
		vector<DataPoint> _data_point_list;
		vector<long int> _index_list;
		vector<float> _radius_list;
		vector<long int> _neighbor_index_list;
		vector<float> _neighbor_radius_list;
		Node *_root;
		Region *_query_region;
		long int _count;
		long int _neighbor_count;
		float _radius;
		float _radius_sq;
		float _neighbor_radius;
		float _neighbor_radius_sq;
		float *_center_coord;
		float *_coords;
		int _bucket_size;
		// Methods
		Node *_build_tree(long int offset_begin=0, long int offset_end=0, int depth=0);
		void _report_subtree(Node *node);
		void _report_point(long int index, float *coord);
		void _test_region(Node *node, Region *region, int depth); 
		void _set_query_region(float *left, float *right);
		void _add_point(long int index, float *coord);
		void _search(Region *region=NULL, Node *node=NULL, int depth=0);
		void _neighbor_search_pairs(Node *left, Region *left_region, Node *right, Region *right_region, int depth);
		void _neighbor_search(Node *root, Region *region, int depth);
		void _search_neighbors_between_buckets(Node *node1, Node *node2);
		void _search_neighbors_in_bucket(Node *node);
		void _test_neighbors(DataPoint &p1, DataPoint &p2);
	public:
		int dim;
		KDTree(int dim, int bucket_size);
		~KDTree();
		void set_data(float *coords, long int nr_points);
		// single neighbor search 
		void search_center_radius(float *coord, float radius);
		long int get_count(void);
		void copy_indices(long int *indices);
		void copy_radii(float *radii);
		// all neighbor search
		long int neighbor_get_count(void);
		void neighbor_search(float neighbor_radius);
		void neighbor_simple_search(float neighbor_radius);
		void neighbor_copy_indices(long int *indices);
		void neighbor_copy_radii(float *radii);
};




