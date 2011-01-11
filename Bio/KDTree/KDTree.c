#include <math.h>
#include <stdlib.h>
#include "KDTree.h"

#define INF 1000000

struct DataPoint;

struct Node;

struct Region;

struct Radius;

struct KDTree
{
    struct DataPoint* _data_point_list;
    int _data_point_list_size;
    struct Radius* _radius_list;
    struct Neighbor* _neighbor_list;
    struct Node *_root;
    struct Region *_query_region;
    long int _count;
    long int _neighbor_count;
    float _radius;
    float _radius_sq;
    float _neighbor_radius;
    float _neighbor_radius_sq;
    float *_center_coord;
    float *_coords;
    int _bucket_size;
    int dim;
};

/* DataPoint */

static int DataPoint_current_dim=0;

struct DataPoint
{
    long int _index;
    float *_coord;
};

static int compare(const void* self, const void* other)
{
    float a, b;
    const struct DataPoint* p = (const struct DataPoint*)self;
    const struct DataPoint* q = (const struct DataPoint*)other;

    a=p->_coord[DataPoint_current_dim];
    b=q->_coord[DataPoint_current_dim];

    if (a < b) return -1;
    if (a > b) return +1;
    return 0; 
}

static void DataPoint_sort(struct DataPoint* list, int n, int i)
{
    /* set sort dimension */
    DataPoint_current_dim=i;
    qsort(list, n, sizeof(struct DataPoint), compare);
}

/* Node */

struct Node
{
    struct Node *_left;
    struct Node *_right;
    float _cut_value;
    int _cut_dim;
    long int _start, _end;
};

static struct Node*
Node_create(float cut_value, int cut_dim, long int start, long int end)
{
    struct Node* node = malloc(sizeof(struct Node));
    if (node==NULL) return NULL;
    node->_left=NULL;
    node->_right=NULL;
    node->_cut_value=cut_value;
    node->_cut_dim=cut_dim;
    /* start and end index in _data_point_list */
    node->_start=start;
    node->_end=end;
    return node;
}

static void Node_destroy(struct Node* node)
{
    if(node==NULL) return;
    Node_destroy(node->_left);
    Node_destroy(node->_right);
    free(node);
}

static int Node_is_leaf(struct Node* node)
{
    if (node->_left==NULL && node->_right==NULL) return 1;
    else return 0;
}

/* Region */

static int Region_dim=3;

struct Region
{
    float *_left; 
    float *_right;
};

static struct Region* Region_create(const float *left, const float *right)
{
    struct Region* region = malloc(sizeof(struct Region));
    if(!region) return NULL;

    region->_left= malloc(Region_dim*sizeof(float));
    region->_right= malloc(Region_dim*sizeof(float));
    if (region->_left==NULL || region->_right==NULL)
    {
        if (region->_left) free(region->_left);
        if (region->_right) free(region->_right);
        free(region);
        return NULL;
    }

    if (left==NULL || right==NULL)
    {
        /* [-INF, INF] */
        int i;
        for (i=0; i<Region_dim; i++)
        {
            region->_left[i]=-INF;
            region->_right[i]=INF;
        }
    }
    else
    {
        int i;
        for (i=0; i<Region_dim; i++)
        {
            region->_left[i]=left[i];
            region->_right[i]=right[i];
        }
    }
    return region;
}

static void Region_destroy(struct Region* region)
{
    if(region==NULL) return;
    if(region->_left) free(region->_left);
    if(region->_right) free(region->_right);
    free(region);
}

static struct Region*
Region_create_intersect_left(struct Region* region, float split_coord, int current_dim)
{
    struct Region* p;
    const float value = region->_right[current_dim];
    region->_right[current_dim]=split_coord;
    p = Region_create(region->_left, region->_right);
    region->_right[current_dim] = value;
    return p;
}

static struct Region*
Region_create_intersect_right(struct Region* region, float split_coord, int current_dim)
{
    struct Region* p;
    const float value = region->_left[current_dim];
    region->_left[current_dim]=split_coord;
    p = Region_create(region->_left, region->_right);
    region->_left[current_dim]=value;
    return p;
}

static int
Region_test_intersect_left(struct Region* region, float split_coord, int current_dim)
{
    float l, r;

    r=region->_right[current_dim];
    l=region->_left[current_dim];

    if (split_coord<l) return -1;
    else if (split_coord<r) return 0; /* split point in interval */
    else return +1;
}

static int
Region_test_intersect_right(struct Region* region, float split_coord, int current_dim)
{
    float l, r;

    r=region->_right[current_dim];
    l=region->_left[current_dim];

    if (split_coord<=l) return -1;
    else if (split_coord<=r) return 0; /* split point in interval */
    else return +1;
}

static int Region_encloses(struct Region* region, float *coord)
{
    int i;
    for (i=0; i<Region_dim; i++)
    {
        if (!(coord[i]>=region->_left[i] && coord[i]<=region->_right[i]))
        {
            return 0;
        }
    }
    return 1;
}

static int
Region_test_intersection(struct Region* this_region, struct Region *query_region, float radius)
{
    int status=2;

    int i;
    for (i=0; i<Region_dim; i++)
    {
        float rs, rq, ls, lq;

        rs=this_region->_right[i];
        ls=this_region->_left[i];
        rq=query_region->_right[i];
        lq=query_region->_left[i];

        if (ls-rq>radius)
        {
            /* outside */
            return 0;
        }
        else if (lq-rs>radius)
        {
            /* outside */
            return 0;
        }
        else if (rs<=rq && ls>=lq)
        {
            /* inside (at least in dim i) */
            if (status > 2) status=2;
        }
        else
        {
            /* overlap (at least in dim i) */
            status=1;
        }
    }
    return status;
}

/* Radius */

struct Radius
{
    long int index;
    float value;
};

/* KDTree */

struct KDTree* KDTree_init(int dim, int bucket_size)
{
    struct KDTree* tree;

    tree = malloc(sizeof(struct KDTree));
    if (!tree) return NULL;

    tree->_center_coord= malloc(dim*sizeof(float));
    if (tree->_center_coord==NULL)
    {
        free(tree);
        return NULL;
    }

    tree->dim=dim;

    Region_dim=tree->dim;

    tree->_query_region=NULL;
    tree->_root=NULL;
    tree->_coords=NULL;
    tree->_radius_list = NULL;
    tree->_count=0;
    tree->_neighbor_count=0;
    tree->_neighbor_list = NULL;
    tree->_bucket_size=bucket_size;
    tree->_data_point_list = NULL;
    tree->_data_point_list_size = 0;

    return tree;
}

void KDTree_destroy(struct KDTree* tree)
{
    /* clean up KD tree */
    if (!tree) return;
    Node_destroy(tree->_root);
    Region_destroy(tree->_query_region);
    if (tree->_center_coord) free(tree->_center_coord);
    if (tree->_coords) free(tree->_coords);
    if (tree->_data_point_list) free(tree->_data_point_list);
    if (tree->_neighbor_list) free(tree->_neighbor_list);
    free(tree);
}

static struct Node *KDTree_build_tree(struct KDTree* tree, long int offset_begin, long int offset_end, int depth)
{
    int localdim;

    if (depth==0)
    {
        /* start with [begin, end+1] */
        offset_begin=0;
        offset_end=tree->_data_point_list_size;
        localdim=0;
    }
    else
    {
        localdim=depth%tree->dim;
    }

    if ((offset_end-offset_begin)<=tree->_bucket_size) 
    {
        /* leaf node */
        return Node_create(-1, localdim, offset_begin, offset_end);
    }
    else
    {
        long int offset_split;
        long int left_offset_begin, left_offset_end;
        long int right_offset_begin, right_offset_end;
        long int d;
        float cut_value;
        struct DataPoint data_point;
        struct Node *left_node, *right_node, *new_node;

        DataPoint_sort(tree->_data_point_list+offset_begin, offset_end-offset_begin, localdim);

        /* calculate index of split point */
        d=offset_end-offset_begin;
        offset_split=d/2+d%2;

        data_point=tree->_data_point_list[offset_begin+offset_split-1];
        cut_value = data_point._coord[localdim];

        /* create new node and bind to left & right nodes */
        new_node=Node_create(cut_value, localdim, offset_begin, offset_end);
        if (new_node==NULL) return NULL;

        /* left */
        left_offset_begin=offset_begin;
        left_offset_end=offset_begin+offset_split;
        left_node=KDTree_build_tree(tree, left_offset_begin, left_offset_end, depth+1);

        /* right */
        right_offset_begin=left_offset_end;
        right_offset_end=offset_end;
        right_node=KDTree_build_tree(tree, right_offset_begin, right_offset_end, depth+1);

        new_node->_left = left_node;
        new_node->_right = right_node;

        if (left_node==NULL || right_node==NULL)
        {
            Node_destroy(new_node);
            return NULL;
        }

        return new_node;
    }
}

static int KDTree_add_point(struct KDTree* tree, long int index, float *coord)
{
    int n;
    struct DataPoint* p;

    n = tree->_data_point_list_size;
    p = realloc(tree->_data_point_list, (n+1)*sizeof(struct DataPoint));
    if (p==NULL) return 0;

    p[n]._index = index;
    p[n]._coord = coord;

    tree->_data_point_list_size = n+1;
    tree->_data_point_list = p;

    return 1;
}

static float KDTree_dist(float *coord1, float *coord2, int dim)
{
    /* returns the SQUARE of the distance between two points */
    int i;
    float sum=0, dif=0;

    for (i=0; i<dim; i++)
    {
        dif=coord1[i]-coord2[i];
        sum+=dif*dif;
    }
    return sum;
}

static int KDTree_report_point(struct KDTree* tree, long int index, float *coord)
{
    float r;

    r=KDTree_dist(tree->_center_coord, coord, tree->dim);

    if (r<=tree->_radius_sq)
    {
        int n = tree->_count;
        struct Radius* p;

    p  = realloc(tree->_radius_list, (n+1)*sizeof(struct Radius));
        if (p==NULL)
        {
            return 0;
        }
        /* note use of sqrt - only calculated if necessary */
        p[n].index = index;
        p[n].value = sqrt(r);
        tree->_radius_list = p;
        tree->_count++;
    }
    return 1;
}

static int KDTree_report_subtree(struct KDTree* tree, struct Node *node)
{
    int ok;
    if (Node_is_leaf(node))
    {
        /* report point(s) */
        long int i;

        for (i=node->_start; i<node->_end; i++)
        {
            struct DataPoint data_point;
            data_point=tree->_data_point_list[i];
            ok = KDTree_report_point(tree, data_point._index, data_point._coord);
            if (!ok) return 0;
        }
    }
    else
    {
        /* find points in subtrees via recursion */
        ok = KDTree_report_subtree(tree, node->_left);
        if (!ok) return 0;
        ok = KDTree_report_subtree(tree, node->_right);
        if (!ok) return 0;
    }
    return 1;
}

long int KDTree_get_count(struct KDTree* tree)
{       
    return tree->_count;
}

long int KDTree_neighbor_get_count(struct KDTree* tree)
{       
    return tree->_neighbor_count;
}

static int KDTree_search(struct KDTree* tree, struct Region *region, struct Node *node, int depth);

static int KDTree_test_region(struct KDTree* tree, struct Node *node, struct Region *region, int depth)
{
    int ok;
    int intersect_flag;

    /* is node region inside, outside or overlapping 
     * with query region? */
    intersect_flag=Region_test_intersection(region, tree->_query_region, 0);

    if (intersect_flag==2)                
    {
        /* inside - extract points */
        ok = KDTree_report_subtree(tree, node);
        /* end of recursion -- get rid of region */
        Region_destroy(region);
        if (!ok) return 0;
    }
    else if (intersect_flag==1)
    {
        /* overlap - recursion */
        ok = KDTree_search(tree, region, node, depth+1);    
        /* search does cleanup of region */
        if (!ok) return 0;
    }
    else
    {
        /* outside - stop */

        /* end of recursion -- get rid of region */
        Region_destroy(region);    
    }
    return 1;
}

static int KDTree_search(struct KDTree* tree, struct Region *region, struct Node *node, int depth)
{
    int current_dim;
    int ok = 1;

    if(depth==0)
    {
        /* start with [-INF, INF] region */
        
        region = Region_create(NULL, NULL);
        if (region==NULL) return 0;

        /* start with root node */
        node=tree->_root;
    }

    current_dim=depth%tree->dim;

    if(Node_is_leaf(node))
    {
        long int i;

        for (i=node->_start; i<node->_end; i++)
        {
            struct DataPoint data_point;

            data_point=tree->_data_point_list[i];

            if (Region_encloses(tree->_query_region, data_point._coord))
            {
                /* point is enclosed in query region - report & stop */
                ok = KDTree_report_point(tree, data_point._index, data_point._coord);
            }
        }
    }
    else
    {
        struct Node *left_node, *right_node;
        struct Region *left_region, *right_region;
        int intersect_left, intersect_right;

        left_node=node->_left;

        /* LEFT HALF PLANE */

        /* new region */
        intersect_left=Region_test_intersect_left(region, node->_cut_value, current_dim);

        if(intersect_left==1)
        {
            left_region = Region_create(region->_left, region->_right);
            if (left_region)
                ok = KDTree_test_region(tree, left_node, left_region, depth);
            else
                ok = 0;
        }
        else if (intersect_left==0)
        {
            left_region = Region_create_intersect_left(region, node->_cut_value, current_dim);
            if (left_region)
                ok = KDTree_test_region(tree, left_node, left_region, depth);
            else
                ok = 0;
        }
        /* intersect_left is -1 if no overlap */

        /* RIGHT HALF PLANE */

        right_node=node->_right;

        /* new region */
        intersect_right=Region_test_intersect_right(region, node->_cut_value, current_dim);
        if (intersect_right==-1)
        {
            right_region = Region_create(region->_left, region->_right);
            /* test for overlap/inside/outside & do recursion/report/stop */
            if (right_region)
                ok = KDTree_test_region(tree, right_node, right_region, depth);
            else
                ok = 0;
        }
        else if (intersect_right==0)
        {
            right_region = Region_create_intersect_right(region, node->_cut_value, current_dim);
            /* test for overlap/inside/outside & do recursion/report/stop */
            if (right_region)
                ok = KDTree_test_region(tree, right_node, right_region, depth);
            else
                ok = 0;
        }
        /* intersect_right is +1 if no overlap */
    }

    Region_destroy(region);
    return ok;
}

int KDTree_set_data(struct KDTree* tree, float *coords, long int nr_points)
{
    long int i;
    int ok;

    Region_dim=tree->dim;

    /* clean up stuff from previous use */
    Node_destroy(tree->_root);
    if (tree->_coords) free(tree->_coords);
    if (tree->_radius_list)
    {
        free(tree->_radius_list);
        tree->_radius_list = NULL;
    }
    tree->_count=0;
    /* keep pointer to coords to delete it */
    tree->_coords=coords;
        
    for (i=0; i<nr_points; i++)
    {
        ok = KDTree_add_point(tree, i, coords+i*tree->dim);
        if (!ok) 
        {
            free(tree->_data_point_list);
            tree->_data_point_list = NULL;
            tree->_data_point_list_size = 0;
            return 0;
        }
    }

    /* build KD tree */
    tree->_root=KDTree_build_tree(tree, 0, 0, 0);
    if(!tree->_root) return 0;
    return 1;
}

int KDTree_search_center_radius(struct KDTree* tree, float *coord, float radius)
{
    int i;
    int dim = tree->dim;
    float* left = malloc(dim*sizeof(float));
    float* right = malloc(dim*sizeof(float));
    if (left==NULL || right==NULL)
    {
        if (left) free(left);
        if (right) free(right);
        return 0;
    }

    Region_dim=tree->dim;

    if (tree->_radius_list)
    {
        free(tree->_radius_list);
        tree->_radius_list = NULL;
    }
    tree->_count=0;

    tree->_radius=radius;
    /* use of r^2 to avoid sqrt use */
    tree->_radius_sq=radius*radius;

    for (i=0; i<tree->dim; i++)
    {
        left[i]=coord[i]-radius;
        right[i]=coord[i]+radius;
        /* set center of query */
        tree->_center_coord[i]=coord[i];
    }

    /* clean up! */
    if (coord) free(coord);

    Region_destroy(tree->_query_region);
    tree->_query_region= Region_create(left, right);

    free(left);
    free(right);

    if (!tree->_query_region) return 0;

    return KDTree_search(tree, NULL, NULL, 0);
}

void KDTree_copy_indices(struct KDTree* tree, long *indices)
{
    long int i;

    if (tree->_count==0) return;

    for(i=0; i<tree->_count; i++)
    {
        indices[i]=tree->_radius_list[i].index;
    }
}

void KDTree_copy_radii(struct KDTree* tree, float *radii)
{
    long int i;

    if (tree->_count==0) return;

    for(i=0; i<tree->_count; i++)
    {
        radii[i]=tree->_radius_list[i].value;
    }
}

static int KDTree_test_neighbors(struct KDTree* tree, struct DataPoint* p1, struct DataPoint* p2)
{
    float r;

    r=KDTree_dist(p1->_coord, p2->_coord, tree->dim);

    if(r<=tree->_neighbor_radius_sq)
    {
        /* we found a neighbor pair! */
        struct Neighbor* p;
        int n;
        n = tree->_neighbor_count;
        p = realloc(tree->_neighbor_list, (n+1)*sizeof(struct Neighbor));
        if (p==NULL) return 0;

        p[n].index1 = p1->_index;
        p[n].index2 = p2->_index;
        /* note sqrt */
        p[n].radius = sqrt(r);
        tree->_neighbor_list = p;
        tree->_neighbor_count++;
    }

    return 1;
}

static int KDTree_search_neighbors_between_buckets(struct KDTree* tree, struct Node *node1, struct Node *node2)
{
    long int i;
    int ok;

    for(i=node1->_start; i<node1->_end; i++)
    {
        struct DataPoint p1;
        long int j;

        p1=tree->_data_point_list[i];

        for (j=node2->_start; j<node2->_end; j++)
        {
            struct DataPoint p2;

            p2=tree->_data_point_list[j];

            ok = KDTree_test_neighbors(tree, &p1, &p2);
            if (!ok) return 0;
        }
    }
    return 1;
}

static int KDTree_neighbor_search_pairs(struct KDTree* tree, struct Node *down, struct Region *down_region, struct Node *up, struct Region *up_region, int depth)
{
    int down_is_leaf, up_is_leaf;
    int localdim;
    int ok = 1;

    /* if regions do not overlap - STOP */
    if (!down || !up || !down_region || !up_region) 
    {
        /* STOP */
        return ok;
    }
    
    if (Region_test_intersection(down_region, up_region, tree->_neighbor_radius)==0)
    {
        /* regions cannot contain neighbors */
        return ok;
    }

    /* dim */
    localdim=depth%tree->dim;

    /* are they leaves? */
    up_is_leaf=Node_is_leaf(up);
    down_is_leaf=Node_is_leaf(down);

    if (up_is_leaf && down_is_leaf)
    {
        /* two leaf nodes */
        ok = KDTree_search_neighbors_between_buckets(tree, down, up);
    }
    else
    {
        /* one or no leaf nodes */

        struct Node *up_right, *up_left, *down_left, *down_right;
        struct Region *up_left_region = NULL;
        struct Region *up_right_region = NULL;
        struct Region *down_left_region = NULL;
        struct Region *down_right_region = NULL;  
        
        if (down_is_leaf)
        {
            down_left=down;
            /* make a copy of down_region */
            down_left_region= Region_create(down_region->_left, down_region->_right);
            if (down_left_region==NULL) ok = 0;
            down_right=NULL;
            down_right_region=NULL;
        }
        else
        {
            float cut_value;
            int intersect;

            cut_value=down->_cut_value;

            down_left=down->_left;
            down_right=down->_right;
            intersect=Region_test_intersect_left(down_region, cut_value, localdim);

            if(intersect==1)
            {
                down_left_region = Region_create(down_region->_left, down_region->_right);
                if (down_left_region==NULL) ok = 0;
            }
            else if(intersect==0)
            {
                down_left_region = Region_create_intersect_left(down_region, cut_value, localdim);
                if (down_left_region==NULL) ok = 0;
            }
            else if(intersect==-1)
            /* intersect is -1 if no overlap */
            {
                down_left_region = NULL;
            }

            intersect=Region_test_intersect_right(down_region, cut_value, localdim); 
            if(intersect==-1)
            {
                down_right_region = Region_create(down_region->_left, down_region->_right);
                if (down_right_region==NULL) ok = 0;
            }
            else if(intersect==0)
            {
                down_right_region = Region_create_intersect_right(down_region, cut_value, localdim);
                if (down_right_region==NULL) ok = 0;
            }
            else if(intersect==+1)
            {
                down_right_region = NULL;
            }
        }

        if (up_is_leaf)
        {
            up_left=up;
            /* make a copy of up_region */
            up_left_region= Region_create(up_region->_left, up_region->_right);
            if (up_left_region==NULL) ok = 0;
            up_right=NULL;
            up_right_region=NULL;
        }
        else
        {
            float cut_value;
            int intersect;

            cut_value=up->_cut_value;

            up_left=up->_left;
            up_right=up->_right;
            intersect=Region_test_intersect_left(up_region, cut_value, localdim);

            if(intersect==1)
            {
                up_left_region = Region_create(up_region->_left, up_region->_right);
                if (up_left_region==NULL) ok = 0;
            }
            else if(intersect==0)
            {
                up_left_region = Region_create_intersect_left(up_region, cut_value, localdim);
                if (up_left_region==NULL) ok = 0;
            }
            else if(intersect==-1)
            /* intersect is -1 if no overlap */
            {
                up_left_region = NULL;
            }

            intersect=Region_test_intersect_right(up_region, cut_value, localdim);
            if(intersect==-1)
            {
                up_right_region = Region_create(up_region->_left, up_region->_right);
                if (up_right_region==NULL) ok = 0;
            }
            else if(intersect==0)
            {
                up_right_region = Region_create_intersect_right(up_region, cut_value, localdim);
                if (up_right_region==NULL) ok = 0;
            }
            else if(intersect==+1)
            /* intersect is +1 if no overlap */
            {
                up_right_region = NULL;
            }
        }

        if (ok)
            ok = KDTree_neighbor_search_pairs(tree, up_left, up_left_region, down_left, down_left_region, depth+1);
        if (ok)
            ok = KDTree_neighbor_search_pairs(tree, up_left, up_left_region, down_right, down_right_region, depth+1);
        if (ok)
            ok = KDTree_neighbor_search_pairs(tree, up_right, up_right_region, down_left, down_left_region, depth+1);
        if (ok)
            ok = KDTree_neighbor_search_pairs(tree, up_right, up_right_region, down_right, down_right_region, depth+1);

        Region_destroy(down_left_region);
        Region_destroy(down_right_region);
        Region_destroy(up_left_region);
        Region_destroy(up_right_region);
    }
    return ok;
}

static int KDTree_search_neighbors_in_bucket(struct KDTree* tree, struct Node *node)
{
    long int i;
    int ok;

    for(i=node->_start; i<node->_end; i++)
    {
        struct DataPoint p1;
        long int j;

        p1=tree->_data_point_list[i];

        for (j=i+1; j<node->_end; j++)
        {
            struct DataPoint p2;

            p2=tree->_data_point_list[j];

            ok = KDTree_test_neighbors(tree, &p1, &p2);
            if (!ok) return 0;
        }
    }
    return 1;
}

static int KDTree__neighbor_search(struct KDTree* tree, struct Node *node, struct Region *region, int depth)
{
    struct Node *left, *right;
    struct Region *left_region = NULL;
    struct Region *right_region = NULL;
    int localdim;
    int intersect;
    float cut_value;
    int ok = 1;

    localdim=depth%tree->dim;

    left=node->_left;
    right=node->_right;

    cut_value = node->_cut_value;

    /* planes of left and right nodes */
    intersect=Region_test_intersect_left(region, cut_value, localdim);
    if(intersect==1)
    {
        left_region = Region_create(region->_left, region->_right);
        if (!left_region) ok = 0;
    }
    else if(intersect==0)
    {
        left_region = Region_create_intersect_left(region, cut_value, localdim);
        if (!left_region) ok = 0;
    }
    else if(intersect==-1)
    /* intersect is -1 if no overlap */
    {
        left_region = NULL;
    }

    intersect=Region_test_intersect_right(region, cut_value, localdim);
    if(intersect==-1)
    {
        right_region = Region_create(region->_left, region->_right);
        if (!right_region) ok = 0;
    }
    else if(intersect==0)
    {
        right_region = Region_create_intersect_right(region, cut_value, localdim);
        if (!right_region) ok = 0;
    }
    else if(intersect==+1)
    /* intersect is +1 if no overlap */
    {
        right_region = NULL;
    }

    if (ok)
    {
        if (!Node_is_leaf(left))
        {
            /* search for pairs in this half plane */
            ok = KDTree__neighbor_search(tree, left, left_region, depth+1);
        }
        else
        {
            ok = KDTree_search_neighbors_in_bucket(tree, left);
        }
    }
    
    if (ok)
    {
        if (!Node_is_leaf(right))
        {
            /* search for pairs in this half plane */
            ok = KDTree__neighbor_search(tree, right, right_region, depth+1);
        }
        else
        {
            ok = KDTree_search_neighbors_in_bucket(tree, right);
        }
    }

    /* search for pairs between the half planes */
    if (ok)
    {
        ok = KDTree_neighbor_search_pairs(tree, left, left_region, right, right_region, depth+1);
    }

    /* cleanup */
    Region_destroy(left_region);
    Region_destroy(right_region);

    return ok;
}

int
KDTree_neighbor_search(struct KDTree* tree, float neighbor_radius,
                       struct Neighbor** neighbors)
{
    long int i;
    int ok;
    Region_dim=tree->dim;

    if(tree->_neighbor_list)
    {
        free(tree->_neighbor_list);
        tree->_neighbor_list = NULL;
    }
    tree->_neighbor_count=0;
    /* note the use of r^2 to avoid use of sqrt */
    tree->_neighbor_radius=neighbor_radius;
    tree->_neighbor_radius_sq=neighbor_radius*neighbor_radius;

    if (Node_is_leaf(tree->_root))
    {
        /* this is a boundary condition */
        /* bucket_size>nr of points */
        ok = KDTree_search_neighbors_in_bucket(tree, tree->_root);
    }
    else
    {
        /* "normal" situation */
        struct Region *region;
        /* start with [-INF, INF] */
        region= Region_create(NULL, NULL);
        if (!region) return 0;
        ok = KDTree__neighbor_search(tree, tree->_root, region, 0);
        Region_destroy(region);
    }
    if (!ok) return 0;

    *neighbors = NULL;
    for (i = 0; i < tree->_neighbor_count; i++)
    {
        struct Neighbor* neighbor = malloc(sizeof(struct Neighbor));
        if (!neighbor)
        {
            while(1)
            {
                neighbor = *neighbors;
                if (!neighbor) return 0;
                *neighbors = neighbor->next;
                free(neighbor);
            }
        }
        *neighbor = tree->_neighbor_list[i];
        neighbor->next = *neighbors;
        *neighbors = neighbor;
    }

    return 1;
}

int
KDTree_neighbor_simple_search(struct KDTree* tree, float radius,
                              struct Neighbor** neighbors)
{
    long int i;
    int ok = 1;

    Region_dim=tree->dim;

    tree->_neighbor_radius=radius;
    tree->_neighbor_radius_sq=radius*radius;

    tree->_neighbor_count=0;
    if (tree->_neighbor_list)
    {
        free(tree->_neighbor_list);
        tree->_neighbor_list = NULL;
    }

    DataPoint_sort(tree->_data_point_list, tree->_data_point_list_size, 0);

    for (i=0; i<tree->_data_point_list_size; i++)
    {
        float x1;
        long int j;
        struct DataPoint p1;

        p1=tree->_data_point_list[i];
        x1=p1._coord[0];

        for (j=i+1; j<tree->_data_point_list_size; j++)
        {
            struct DataPoint p2;
            float x2;

            p2=tree->_data_point_list[j];
            x2=p2._coord[0];

            if (fabs(x2-x1)<=radius)
            {
                ok = KDTree_test_neighbors(tree, &p1, &p2);
                if (!ok) break;
            }
            else
            {
                break;
            }
        }
    }

    if (!ok) return 0;

    *neighbors = NULL;

    for (i = 0; i < tree->_neighbor_count; i++)
    {
        struct Neighbor* neighbor = malloc(sizeof(struct Neighbor));
        if (!neighbor)
        {
            while(1)
            {
                neighbor = *neighbors;
                if (!neighbor) return 0;
                *neighbors = neighbor->next;
                free(neighbor);
            }
        }
        *neighbor = tree->_neighbor_list[i];
        neighbor->next = *neighbors;
        *neighbors = neighbor;
    }

    return 1;
}
