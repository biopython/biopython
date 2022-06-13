/* Copyright 2018-2020 by Michiel de Hoon.  All rights reserved.
 *
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */

#include <math.h>
#include <stdlib.h>
#include <Python.h>

#define INF 1000000

#define DIM 3 /* three spatial dimensions */

/* DataPoint */

static int DataPoint_current_dim = 0;

typedef struct
{
    long int _index;
    double _coord[DIM];
} DataPoint;

static int compare(const void* self, const void* other)
{
    const DataPoint* p = self;
    const DataPoint* q = other;
    const double a = p->_coord[DataPoint_current_dim];
    const double b = q->_coord[DataPoint_current_dim];
    if (a < b) return -1;
    if (a > b) return +1;
    return 0;
}

static void DataPoint_sort(DataPoint* list, int n, int i)
{
    /* set sort dimension */
    DataPoint_current_dim = i;
    qsort(list, n, sizeof(DataPoint), compare);
}

/* Point */

typedef struct {
    PyObject_HEAD
    long int index;
    double radius;
} Point;

static int
Point_init(Point *self, PyObject *args, PyObject *kwds)
{
    int index;
    double radius = 0.0;
    static char *kwlist[] = {"index", "radius", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|d", kwlist,
                                     &index, &radius))
        return -1;

    self->index = index;
    self->radius = radius;
    return 0;
}

static PyObject*
Point_repr(Point* self)
{
    return PyUnicode_FromFormat("%ld: %g", self->index, self->radius);
}

static char Point_index__doc__[] =
"index";

static PyObject*
Point_getindex(Point* self, void* closure)
{
    return PyLong_FromLong(self->index);
}

static char Point_radius__doc__[] = "the radius";

static PyObject*
Point_getradius(Point* self, void* closure)
{
    const double value = self->radius;
    return PyFloat_FromDouble(value);
}

static PyGetSetDef Point_getset[] = {
    {"index", (getter)Point_getindex, NULL, Point_index__doc__, NULL},
    {"radius", (getter)Point_getradius, NULL, Point_radius__doc__, NULL},
    {NULL}  /* Sentinel */
};

static char Point_doc[] =
"A single point; attributes are index and radius.\n";

static PyTypeObject PointType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Point",                   /* tp_name*/
    sizeof(Point),             /* tp_basicsize*/
    0,                         /* tp_itemsize*/
    0,                         /* tp_dealloc*/
    0,                         /* tp_print*/
    0,                         /* tp_getattr*/
    0,                         /* tp_setattr*/
    0,                         /* tp_compare*/
    (reprfunc)Point_repr,      /* tp_repr*/
    0,                         /* tp_as_number*/
    0,                         /* tp_as_sequence*/
    0,                         /* tp_as_mapping*/
    0,                         /* tp_hash */
    0,                         /* tp_call*/
    0,                         /* tp_str*/
    0,                         /* tp_getattro*/
    0,                         /* tp_setattro*/
    0,                         /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /* tp_flags*/
    Point_doc,                 /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    Point_getset,              /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Point_init,      /* tp_init */
};

/* Neighbor */

typedef struct {
    PyObject_HEAD
    long int index1;
    long int index2;
    double radius;
} Neighbor;

static int
Neighbor_init(Neighbor *self, PyObject *args, PyObject *kwds)
{
    int index1, index2;
    double radius = 0.0;
    static char *kwlist[] = {"index1", "index2", "radius", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ii|d", kwlist,
                                     &index1, &index2, &radius))
        return -1;
    self->index1 = index1;
    self->index2 = index2;
    self->radius = radius;

    return 0;
}

static PyObject*
Neighbor_repr(Neighbor* self)
{
    return PyUnicode_FromFormat("(%ld, %ld): %g",
                                self->index1, self->index2, self->radius);
}

static char Neighbor_index1__doc__[] =
"index of the first neighbor";

static PyObject*
Neighbor_getindex1(Neighbor* self, void* closure)
{
    return PyLong_FromLong(self->index1);
}

static char Neighbor_index2__doc__[] =
"index of the second neighbor";

static PyObject*
Neighbor_getindex2(Neighbor* self, void* closure)
{
    return PyLong_FromLong(self->index2);
}

static char Neighbor_radius__doc__[] = "the radius";

static PyObject*
Neighbor_getradius(Neighbor* self, void* closure)
{
    const double value = self->radius;
    return PyFloat_FromDouble(value);
}

static PyGetSetDef Neighbor_getset[] = {
    {"index1", (getter)Neighbor_getindex1, NULL, Neighbor_index1__doc__, NULL},
    {"index2", (getter)Neighbor_getindex2, NULL, Neighbor_index2__doc__, NULL},
    {"radius", (getter)Neighbor_getradius, NULL, Neighbor_radius__doc__, NULL},
    {NULL}  /* Sentinel */
};

static char Neighbor_doc[] =
"A neighbor pair; attributes are index1, index2, and radius.\n";

static PyTypeObject NeighborType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Neighbor",                /* tp_name*/
    sizeof(Neighbor),          /* tp_basicsize*/
    0,                         /* tp_itemsize*/
    0,                         /* tp_dealloc*/
    0,                         /* tp_print*/
    0,                         /* tp_getattr*/
    0,                         /* tp_setattr*/
    0,                         /* tp_compare*/
    (reprfunc)Neighbor_repr,   /* tp_repr*/
    0,                         /* tp_as_number*/
    0,                         /* tp_as_sequence*/
    0,                         /* tp_as_mapping*/
    0,                         /* tp_hash */
    0,                         /* tp_call*/
    0,                         /* tp_str*/
    0,                         /* tp_getattro*/
    0,                         /* tp_setattro*/
    0,                         /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /* tp_flags*/
    Neighbor_doc,              /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    Neighbor_getset,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Neighbor_init,   /* tp_init */
};

/* Node */

typedef struct Node
{
    struct Node *_left;
    struct Node *_right;
    double _cut_value;
    int _cut_dim;
    long int _start, _end;
} Node;

static Node*
Node_create(double cut_value, int cut_dim, long int start, long int end)
{
    Node* node = PyMem_Malloc(sizeof(Node));
    if (node == NULL) return NULL;
    node->_left = NULL;
    node->_right = NULL;
    node->_cut_value = cut_value;
    node->_cut_dim = cut_dim;
    /* start and end index in _data_point_list */
    node->_start = start;
    node->_end = end;
    return node;
}

static void Node_destroy(Node* node)
{
    if (node == NULL) return;
    Node_destroy(node->_left);
    Node_destroy(node->_right);
    PyMem_Free(node);
}

static int Node_is_leaf(Node* node)
{
    if (node->_left == NULL && node->_right == NULL) return 1;
    else return 0;
}

/* Region */

typedef struct
{
    double _left[DIM];
    double _right[DIM];
} Region;

static Region* Region_create(const double *left, const double *right)
{
    int i;
    Region* region = PyMem_Malloc(sizeof(Region));
    if (region == NULL) return NULL;

    if (left == NULL || right == NULL)
    {
        /* [-INF, INF] */
        for (i = 0; i < DIM; i++) {
            region->_left[i] = -INF;
            region->_right[i] = INF;
        }
    }
    else
    {
        for (i = 0; i < DIM; i++) {
            region->_left[i] = left[i];
            region->_right[i] = right[i];
        }
    }
    return region;
}

static void Region_destroy(Region* region)
{
    if (region) PyMem_Free(region);
}

static int Region_encloses(Region* region, double *coord)
{
    int i;
    for (i = 0; i < DIM; i++)
    {
        if (!(coord[i] >= region->_left[i] && coord[i] <= region->_right[i]))
        {
            return 0;
        }
    }
    return 1;
}

static int
Region_test_intersect_left(Region* region, double split_coord, int current_dim)
{
    const double r = region->_right[current_dim];
    const double l = region->_left[current_dim];
    if (split_coord < l) return -1;
    else if (split_coord < r) return 0; /* split point in interval */
    else return +1;
}

static int
Region_test_intersect_right(Region* region, double split_coord, int current_dim)
{
    const double r = region->_right[current_dim];
    const double l = region->_left[current_dim];
    if (split_coord <= l) return -1;
    else if (split_coord <= r) return 0; /* split point in interval */
    else return +1;
}

static int
Region_test_intersection(Region* this_region, Region *query_region, double radius)
{
    int status = 2;

    int i;
    for (i = 0; i < DIM; i++)
    {
        double rs = this_region->_right[i];
        double ls = this_region->_left[i];
        double rq = query_region->_right[i];
        double lq = query_region->_left[i];

        if (ls-rq > radius)
        {
            /* outside */
            return 0;
        }
        else if (lq-rs > radius)
        {
            /* outside */
            return 0;
        }
        else if (rs <= rq && ls>=lq)
        {
            /* inside (at least in dim i) */
            if (status > 2) status = 2;
        }
        else
        {
            /* overlap (at least in dim i) */
            status = 1;
        }
    }
    return status;
}

static Region*
Region_create_intersect_left(Region* region, double split_coord, int current_dim)
{
    Region* p;
    const double value = region->_right[current_dim];
    region->_right[current_dim] = split_coord;
    p = Region_create(region->_left, region->_right);
    region->_right[current_dim] = value;
    return p;
}

static Region*
Region_create_intersect_right(Region* region, double split_coord, int current_dim)
{
    Region* p;
    const double value = region->_left[current_dim];
    region->_left[current_dim] = split_coord;
    p = Region_create(region->_left, region->_right);
    region->_left[current_dim] = value;
    return p;
}

/* Radius */

typedef struct
{
    long int index;
    double value;
} Radius;

/* KDTree */

typedef struct {
    PyObject_HEAD
    DataPoint* _data_point_list;
    int _data_point_list_size;
    Node *_root;
    int _bucket_size;
    /* The following are temporary variables used during a search only. */
    double _radius;
    double _radius_sq;
    double _neighbor_radius;
    double _neighbor_radius_sq;
    double _center_coord[DIM];
} KDTree;

static double KDTree_dist(double *coord1, double *coord2)
{
    /* returns the SQUARE of the distance between two points */
    int i;
    double sum = 0, dif = 0;

    for (i = 0; i < DIM; i++) {
        dif = coord1[i]-coord2[i];
        sum += dif*dif;
    }
    return sum;
}

static int
KDTree_report_point(KDTree* self, DataPoint* data_point, PyObject* points)
{
    int ok;
    long int index = data_point->_index;
    double *coord = data_point->_coord;
    const double r = KDTree_dist(self->_center_coord, coord);
    if (r <= self->_radius_sq)
    {
        Point* point;
        point = (Point*) PointType.tp_alloc(&PointType, 0);
        if (!point) return 0;
        point->index = index;
        point->radius = sqrt(r); /* note sqrt */
        ok = PyList_Append(points, (PyObject*)point);
        Py_DECREF(point);
        if (ok == -1) return 0;
    }
    return 1;
}

static int
KDTree_test_neighbors(KDTree* self, DataPoint* p1, DataPoint* p2, PyObject* neighbors)
{
    int ok;
    const double r = KDTree_dist(p1->_coord, p2->_coord);
    if (r <= self->_neighbor_radius_sq)
    {
        /* we found a neighbor pair! */
        Neighbor* neighbor;
        long int index1, index2;
        neighbor = (Neighbor*) NeighborType.tp_alloc(&NeighborType, 0);
        if (!neighbor) return 0;
        index1 = p1->_index;
        index2 = p2->_index;
        if (index1 < index2) {
            neighbor->index1 = index1;
            neighbor->index2 = index2;
        }
        else {
            neighbor->index1 = index2;
            neighbor->index2 = index1;
        }
        neighbor->radius = sqrt(r); /* note sqrt */
        ok = PyList_Append(neighbors, (PyObject*)neighbor);
        Py_DECREF(neighbor);
        if (ok == -1) return 0;
    }

    return 1;
}

static int
KDTree_search_neighbors_in_bucket(KDTree* self, Node *node, PyObject* neighbors)
{
    long int i;
    int ok;

    for (i = node->_start; i < node->_end; i++)
    {
        DataPoint p1;
        long int j;

        p1 = self->_data_point_list[i];

        for (j = i+1; j < node->_end; j++) {
            DataPoint p2 = self->_data_point_list[j];
            ok = KDTree_test_neighbors(self, &p1, &p2, neighbors);
            if (!ok) return 0;
        }
    }
    return 1;
}

static int KDTree_search_neighbors_between_buckets(KDTree* self, Node *node1, Node *node2, PyObject* neighbors)
{
    long int i;
    int ok;

    for (i = node1->_start; i < node1->_end; i++)
    {
        DataPoint p1;
        long int j;

        p1 = self->_data_point_list[i];

        for (j = node2->_start; j < node2->_end; j++)
        {
            DataPoint p2 = self->_data_point_list[j];
            ok = KDTree_test_neighbors(self, &p1, &p2, neighbors);
            if (!ok) return 0;
        }
    }
    return 1;
}

static int KDTree_neighbor_search_pairs(KDTree* self, Node *down, Region *down_region, Node *up, Region *up_region, int depth, PyObject* neighbors)
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

    if (Region_test_intersection(down_region, up_region, self->_neighbor_radius)== 0)
    {
        /* regions cannot contain neighbors */
        return ok;
    }

    /* dim */
    localdim = depth % DIM;

    /* are they leaves? */
    up_is_leaf = Node_is_leaf(up);
    down_is_leaf = Node_is_leaf(down);

    if (up_is_leaf && down_is_leaf)
    {
        /* two leaf nodes */
        ok = KDTree_search_neighbors_between_buckets(self, down, up, neighbors);
    }
    else
    {
        /* one or no leaf nodes */

        Node *up_right, *up_left, *down_left, *down_right;
        Region *up_left_region = NULL;
        Region *up_right_region = NULL;
        Region *down_left_region = NULL;
        Region *down_right_region = NULL;

        if (down_is_leaf)
        {
            down_left = down;
            /* make a copy of down_region */
            down_left_region = Region_create(down_region->_left, down_region->_right);
            if (down_left_region == NULL) ok = 0;
            down_right = NULL;
            down_right_region = NULL;
        }
        else
        {
            double cut_value;
            int intersect;

            cut_value = down->_cut_value;

            down_left = down->_left;
            down_right = down->_right;
            intersect = Region_test_intersect_left(down_region, cut_value, localdim);
            switch (intersect) {
                case 1:
                    down_left_region = Region_create(down_region->_left, down_region->_right);
                    if (down_left_region == NULL) ok = 0;
                    break;
                case 0:
                    down_left_region = Region_create_intersect_left(down_region, cut_value, localdim);
                    if (down_left_region == NULL) ok = 0;
                    break;
                case -1: /* intersect is -1 if no overlap */
                    down_left_region = NULL;
                    break;
            }

            intersect = Region_test_intersect_right(down_region, cut_value, localdim);
            switch (intersect) {
                case -1:
                    down_right_region = Region_create(down_region->_left, down_region->_right);
                    if (down_right_region == NULL) ok = 0;
                    break;
                case 0:
                    down_right_region = Region_create_intersect_right(down_region, cut_value, localdim);
                    if (down_right_region == NULL) ok = 0;
                    break;
                case +1:
                    down_right_region = NULL;
                    break;
            }
        }

        if (up_is_leaf)
        {
            up_left = up;
            /* make a copy of up_region */
            up_left_region = Region_create(up_region->_left, up_region->_right);
            if (up_left_region == NULL) ok = 0;
            up_right = NULL;
            up_right_region = NULL;
        }
        else
        {
            double cut_value;
            int intersect;

            cut_value = up->_cut_value;

            up_left = up->_left;
            up_right = up->_right;
            intersect = Region_test_intersect_left(up_region, cut_value, localdim);

            switch (intersect) {
                case 1:
                    up_left_region = Region_create(up_region->_left, up_region->_right);
                    if (up_left_region == NULL) ok = 0;
                    break;
                case 0:
                    up_left_region = Region_create_intersect_left(up_region, cut_value, localdim);
                    if (up_left_region == NULL) ok = 0;
                    break;
                case -1: /* intersect is -1 if no overlap */
                    up_left_region = NULL;
                    break;
            }

            intersect = Region_test_intersect_right(up_region, cut_value, localdim);
            switch (intersect) {
                case -1:
                    up_right_region = Region_create(up_region->_left, up_region->_right);
                    if (up_right_region == NULL) ok = 0;
                    break;
                case 0:
                    up_right_region = Region_create_intersect_right(up_region, cut_value, localdim);
                    if (up_right_region == NULL) ok = 0;
                    break;
                case +1: /* intersect is +1 if no overlap */
                    up_right_region = NULL;
                    break;
            }
        }

        if (ok)
            ok = KDTree_neighbor_search_pairs(self, up_left, up_left_region, down_left, down_left_region, depth+1, neighbors);
        if (ok)
            ok = KDTree_neighbor_search_pairs(self, up_left, up_left_region, down_right, down_right_region, depth+1, neighbors);
        if (ok)
            ok = KDTree_neighbor_search_pairs(self, up_right, up_right_region, down_left, down_left_region, depth+1, neighbors);
        if (ok)
            ok = KDTree_neighbor_search_pairs(self, up_right, up_right_region, down_right, down_right_region, depth+1, neighbors);

        Region_destroy(down_left_region);
        Region_destroy(down_right_region);
        Region_destroy(up_left_region);
        Region_destroy(up_right_region);
    }
    return ok;
}

static int KDTree_neighbor_search(KDTree* self, Node *node, Region *region, int depth, PyObject* neighbors)
{
    Node *left, *right;
    Region *left_region = NULL;
    Region *right_region = NULL;
    int localdim;
    int intersect;
    double cut_value;
    int ok = 1;

    localdim = depth % DIM;

    left = node->_left;
    right = node->_right;

    cut_value = node->_cut_value;

    /* planes of left and right nodes */
    intersect = Region_test_intersect_left(region, cut_value, localdim);
    switch (intersect) {
        case 1:
            left_region = Region_create(region->_left, region->_right);
            if (!left_region) ok = 0;
            break;
        case 0:
            left_region = Region_create_intersect_left(region, cut_value, localdim);
            if (!left_region) ok = 0;
            break;
        case -1: /* intersect is -1 if no overlap */
            left_region = NULL;
            break;
    }

    intersect = Region_test_intersect_right(region, cut_value, localdim);
    switch (intersect) {
        case -1:
            right_region = Region_create(region->_left, region->_right);
            if (!right_region) ok = 0;
            break;
        case 0:
            right_region = Region_create_intersect_right(region, cut_value, localdim);
            if (!right_region) ok = 0;
            break;
        case +1: /* intersect is +1 if no overlap */
            right_region = NULL;
            break;
    }

    if (ok)
    {
        if (!Node_is_leaf(left))
        {
            /* search for pairs in this half plane */
            ok = KDTree_neighbor_search(self, left, left_region, depth+1, neighbors);
        }
        else
        {
            ok = KDTree_search_neighbors_in_bucket(self, left, neighbors);
        }
    }

    if (ok)
    {
        if (!Node_is_leaf(right))
        {
            /* search for pairs in this half plane */
            ok = KDTree_neighbor_search(self, right, right_region, depth+1, neighbors);
        }
        else
        {
            ok = KDTree_search_neighbors_in_bucket(self, right, neighbors);
        }
    }

    /* search for pairs between the half planes */
    if (ok)
    {
        ok = KDTree_neighbor_search_pairs(self, left, left_region, right, right_region, depth+1, neighbors);
    }

    /* cleanup */
    Region_destroy(left_region);
    Region_destroy(right_region);

    return ok;
}

static Node *
KDTree_build_tree(KDTree* self, long int offset_begin, long int offset_end, int depth)
{
    int localdim;

    if (depth == 0)
    {
        /* start with [begin, end+1] */
        offset_begin = 0;
        offset_end = self->_data_point_list_size;
        localdim = 0;
    }
    else
    {
        localdim = depth % DIM;
    }

    if ((offset_end-offset_begin) <= self->_bucket_size)
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
        double cut_value;
        DataPoint data_point;
        Node *left_node, *right_node, *new_node;

        DataPoint_sort(self->_data_point_list+offset_begin, offset_end-offset_begin, localdim);

        /* calculate index of split point */
        d = offset_end-offset_begin;
        offset_split = d/2+d%2;

        data_point = self->_data_point_list[offset_begin+offset_split-1];
        cut_value = data_point._coord[localdim];

        /* create new node and bind to left & right nodes */
        new_node = Node_create(cut_value, localdim, offset_begin, offset_end);
        if (new_node == NULL) return NULL;

        /* left */
        left_offset_begin = offset_begin;
        left_offset_end = offset_begin+offset_split;
        left_node = KDTree_build_tree(self, left_offset_begin, left_offset_end, depth+1);

        /* right */
        right_offset_begin = left_offset_end;
        right_offset_end = offset_end;
        right_node = KDTree_build_tree(self, right_offset_begin, right_offset_end, depth+1);

        new_node->_left = left_node;
        new_node->_right = right_node;

        if (left_node == NULL || right_node == NULL)
        {
            Node_destroy(new_node);
            return NULL;
        }

        return new_node;
    }
}

static int KDTree_report_subtree(KDTree* self, Node *node, PyObject* points)
{
    int ok;
    if (Node_is_leaf(node)) {
        /* report point(s) */
        long int i;
        for (i = node->_start; i < node->_end; i++) {
            ok = KDTree_report_point(self, &self->_data_point_list[i], points);
            if (!ok) return 0;
        }
    }
    else {
        /* find points in subtrees via recursion */
        ok = KDTree_report_subtree(self, node->_left, points);
        if (!ok) return 0;
        ok = KDTree_report_subtree(self, node->_right, points);
        if (!ok) return 0;
    }
    return 1;
}

static int
KDTree_search(KDTree* self, Region *region, Node *node, int depth, Region* query_region, PyObject* points);

static int KDTree_test_region(KDTree* self, Node *node, Region *region, int depth, Region* query_region, PyObject* points)
{
    int ok;
    int intersect_flag;

    /* is node region inside, outside or overlapping
     * with query region? */
    intersect_flag = Region_test_intersection(region, query_region, 0);

    switch (intersect_flag) {
        case 2:
            /* inside - extract points */
            ok = KDTree_report_subtree(self, node, points);
            /* end of recursion -- get rid of region */
            Region_destroy(region);
            break;
        case 1:
            /* overlap - recursion */
            ok = KDTree_search(self, region, node, depth+1, query_region, points);
            /* search does cleanup of region */
            break;
        default:
            /* outside - stop */
            ok = 1;
            /* end of recursion -- get rid of region */
            Region_destroy(region);
            break;
    }
    return ok;
}

static int
KDTree_search(KDTree* self, Region *region, Node *node, int depth, Region* query_region, PyObject* points)
{
    int current_dim;
    int ok = 1;

    if (depth == 0)
    {
        /* start with [-INF, INF] region */

        region = Region_create(NULL, NULL);
        if (region == NULL) return 0;

        /* start with root node */
        node = self->_root;
    }

    current_dim = depth % DIM;

    if (Node_is_leaf(node)) {
        long int i;
        DataPoint* data_point;
        for (i = node->_start; i < node->_end; i++) {
            data_point = &self->_data_point_list[i];
            if (Region_encloses(query_region, data_point->_coord)) {
                /* point is enclosed in query region - report & stop */
                ok = KDTree_report_point(self, data_point, points);
            }
        }
    }
    else {
        Node *left_node, *right_node;
        Region *left_region, *right_region;
        int intersect_left, intersect_right;

        left_node = node->_left;

        /* LEFT HALF PLANE */

        /* new region */
        intersect_left = Region_test_intersect_left(region, node->_cut_value, current_dim);

        switch (intersect_left) {
            case 1:
                left_region = Region_create(region->_left, region->_right);
                if (left_region)
                    ok = KDTree_test_region(self, left_node, left_region, depth, query_region, points);
                else
                    ok = 0;
                break;
            case 0:
                left_region = Region_create_intersect_left(region, node->_cut_value, current_dim);
                if (left_region)
                    ok = KDTree_test_region(self, left_node, left_region, depth, query_region, points);
                else
                    ok = 0;
                break;
            case -1:
                /* intersect_left is -1 if no overlap */
                break;
        }

        /* RIGHT HALF PLANE */

        right_node = node->_right;

        /* new region */
        intersect_right = Region_test_intersect_right(region, node->_cut_value, current_dim);
        switch (intersect_right) {
            case -1:
                right_region = Region_create(region->_left, region->_right);
                /* test for overlap/inside/outside & do recursion/report/stop */
                if (right_region)
                    ok = KDTree_test_region(self, right_node, right_region, depth, query_region, points);
                else
                    ok = 0;
                break;
            case 0:
                right_region = Region_create_intersect_right(region, node->_cut_value, current_dim);
                /* test for overlap/inside/outside & do recursion/report/stop */
                if (right_region)
                    ok = KDTree_test_region(self, right_node, right_region, depth, query_region, points);
                else
                    ok = 0;
                break;
            case +1:
                /* intersect_right is +1 if no overlap */
                break;
        }
    }

    Region_destroy(region);
    return ok;
}

/* Python interface */

static void
KDTree_dealloc(KDTree* self)
{
    Node_destroy(self->_root);
    if (self->_data_point_list) PyMem_Free(self->_data_point_list);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
KDTree_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
    int bucket_size = 1;
    double* coords;
    Py_ssize_t n, i, j;
    PyObject *obj;
    const int flags = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    Py_buffer view;
    KDTree* self;
    DataPoint* data_point_list;
    double value;

    if (!PyArg_ParseTuple(args, "O|i:KDTree_new" , &obj, &bucket_size))
        return NULL;

    if (bucket_size <= 0) {
        PyErr_SetString(PyExc_ValueError, "bucket size should be positive");
        return NULL;
    }

    if (PyObject_GetBuffer(obj, &view, flags) == -1) return NULL;
    if (view.itemsize != sizeof(double)) {
        PyBuffer_Release(&view);
        PyErr_SetString(PyExc_RuntimeError,
                        "coords array has incorrect data type");
        return NULL;
    }
    if (view.ndim != 2 || view.shape[1] != 3) {
        PyBuffer_Release(&view);
        PyErr_SetString(PyExc_ValueError, "expected a Nx3 numpy array");
        return NULL;
    }
    n = view.shape[0];

    data_point_list = PyMem_Malloc(n*sizeof(DataPoint));
    if (data_point_list == NULL) {
        /* KDTree_dealloc will deallocate data already stored in KDTree */
        PyBuffer_Release(&view);
        return PyErr_NoMemory();
    }

    coords = view.buf;
    for (i = 0; i < n; i++) {
        data_point_list[i]._index = i;
        for (j = 0; j < DIM; j++, coords++) {
            value = *coords;
            if (value <= -1e6 || value >= 1e6) {
                PyMem_Free(data_point_list);
                PyBuffer_Release(&view);
                PyErr_SetString(PyExc_ValueError,
                    "coordinate values should lie between -1e6 and 1e6");
                return NULL;
            }
            data_point_list[i]._coord[j] = value;
        }
    }
    PyBuffer_Release(&view);

    /* build KD tree */
    self = (KDTree*)type->tp_alloc(type, 0);
    if (!self) {
        PyMem_Free(data_point_list);
        return NULL;
    }
    self->_bucket_size = bucket_size;
    self->_data_point_list = data_point_list;
    self->_data_point_list_size = n;

    self->_root = KDTree_build_tree(self, 0, 0, 0);
    if (!self->_root) {
        Py_DECREF(self);
        return PyErr_NoMemory();
    }
    return (PyObject*)self;
}

PyDoc_STRVAR(PyKDTree_search__doc__,
"Search all points within the given radius of center.\n\
\n\
Arguments:\n\
 - center: NumPy array of size 3.\n\
 - radius: float>0\n\
\n\
Returns a list of Point objects; each neighbor has an attribute\n\
index corresponding to the index of the point, and an attribute\n\
radius with the radius between them.");


static PyObject*
PyKDTree_search(KDTree* self, PyObject* args)
{
    PyObject *obj;
    double radius;
    long int i;
    double *coords;
    const int flags = PyBUF_ND | PyBUF_C_CONTIGUOUS;
    Py_buffer view;
    double left[DIM];
    double right[DIM];
    Region* query_region = NULL;
    PyObject* points = NULL;

    if (!PyArg_ParseTuple(args, "Od:search", &obj, &radius))
        return NULL;

    if (radius <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    if (PyObject_GetBuffer(obj, &view, flags) == -1) return NULL;
    if (view.itemsize != sizeof(double)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "coords array has incorrect data type");
        goto exit;
    }
    if (view.ndim != 1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "coords array must be one-dimensional");
        goto exit;
    }
    if (view.shape[0] != DIM) {
        PyErr_SetString(PyExc_RuntimeError,
                        "coords array dimension must be 3");
        goto exit;
    }
    coords = view.buf;

    self->_radius = radius;
    /* use of r^2 to avoid sqrt use */
    self->_radius_sq = radius*radius;

    for (i = 0; i < DIM; i++)
    {
        left[i] = coords[i] - radius;
        right[i] = coords[i] + radius;
        /* set center of query */
        self->_center_coord[i] = coords[i];
    }

    query_region = Region_create(left, right);

    if (!query_region) {
        PyErr_NoMemory();
        goto exit;
    }

    points = PyList_New(0);
    if (!points) goto exit;

    if (!KDTree_search(self, NULL, NULL, 0, query_region, points)) {
        PyErr_NoMemory();
        Py_DECREF(points);
        points = NULL;
        goto exit;
    }

exit:
    if (query_region) Region_destroy(query_region);
    PyBuffer_Release(&view);
    return points;
}

PyDoc_STRVAR(PyKDTree_neighbor_search__doc__,
"All fixed neighbor search.\n\
\n\
Find all point pairs that are within radius of each other.\n\
\n\
Arguments:\n\
 - radius: float (>0)\n\
\n\
Returns a list of Neighbor objects; each neighbor has attributes\n\
index1, index2 corresponding to the indices of the point pair,\n\
and an attribute radius with the radius between them.");


static PyObject*
PyKDTree_neighbor_search(KDTree* self, PyObject* args)
{
    int ok = 0;
    double radius;
    PyObject* neighbors;

    if (!PyArg_ParseTuple(args, "d:neighbor_search", &radius))
        return NULL;

    if (radius <= 0) {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    neighbors = PyList_New(0);

    /* note the use of r^2 to avoid use of sqrt */
    self->_neighbor_radius = radius;
    self->_neighbor_radius_sq = radius*radius;

    if (Node_is_leaf(self->_root)) {
        /* this is a boundary condition */
        /* bucket_size > nr of points */
        ok = KDTree_search_neighbors_in_bucket(self, self->_root, neighbors);
    }
    else {
        /* "normal" situation */
        /* start with [-INF, INF] */
        Region *region = Region_create(NULL, NULL);
        if (region) {
            ok = KDTree_neighbor_search(self, self->_root, region, 0, neighbors);
            Region_destroy(region);
        }
    }
    if (!ok) {
        Py_DECREF(neighbors);
        return PyErr_NoMemory();
    }
    return neighbors;
}

PyDoc_STRVAR(PyKDTree_neighbor_simple_search__doc__,
"All fixed neighbor search (for testing purposes only).\n\
\n\
Find all point pairs that are within radius of each other, using a simple\n\
but slow algorithm. This function is provided to be able to verify the\n\
correctness of fast algorithm using the KD Tree for testing purposes.\n\
\n\
Arguments:\n\
 - radius: float (>0)\n\
\n\
Returns a list of Neighbor objects; each neighbor has attributes\n\
index1, index2 corresponding to the indices of the point pair,\n\
and an attribute radius with the radius between them.");

static PyObject*
PyKDTree_neighbor_simple_search(KDTree* self, PyObject* args)
{
    int ok;
    double radius;
    PyObject* neighbors;
    Py_ssize_t i;

    if (!PyArg_ParseTuple(args, "d:neighbor_simple_search", &radius))
        return NULL;

    if (radius <= 0) {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    neighbors = PyList_New(0);
    if (!neighbors) return NULL;

    self->_neighbor_radius = radius;
    self->_neighbor_radius_sq = radius*radius;

    DataPoint_sort(self->_data_point_list, self->_data_point_list_size, 0);

    for (i = 0; i < self->_data_point_list_size; i++) {
        double x1;
        long int j;
        DataPoint p1;

        p1 = self->_data_point_list[i];
        x1 = p1._coord[0];

        for (j = i+1; j < self->_data_point_list_size; j++) {
            DataPoint p2 = self->_data_point_list[j];
            double x2 = p2._coord[0];
            if (fabs(x2-x1) <= radius)
            {
                ok = KDTree_test_neighbors(self, &p1, &p2, neighbors);
                if (!ok) return PyErr_NoMemory();
            }
            else
            {
                break;
            }
        }
    }
    return neighbors;
}

static PyMethodDef KDTree_methods[] = {
    {"search",
     (PyCFunction)PyKDTree_search,
      METH_VARARGS,
      PyKDTree_search__doc__},
    {"neighbor_search",
     (PyCFunction)PyKDTree_neighbor_search,
      METH_VARARGS,
      PyKDTree_neighbor_search__doc__},
    {"neighbor_simple_search",
     (PyCFunction)PyKDTree_neighbor_simple_search,
      METH_VARARGS,
      PyKDTree_neighbor_simple_search__doc__},
    {NULL}  /* Sentinel */
};

PyDoc_STRVAR(KDTree_doc,
"KDTree(coordinates, bucket_size=1) -> new KDTree\n\
\n\
Create a new KDTree object for the given coordinates and bucket size,\n\
where coordinates is an Nx3 NumPy array (N being the number of points).\n\
\n\
The KDTree data structure can be used for neighbor searches (find all\n\
points within a radius of a given point) and for finding all point pairs\n\
within a certain radius of each other.\n\
\n\
Reference:\n\
\n\
Computational Geometry: Algorithms and Applications\n\
Second Edition\n\
Mark de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf\n\
published by Springer-Verlag\n\
2nd rev. ed. 2000.\n\
ISBN: 3-540-65620-0\n\
\n\
The KD tree data structure is described in chapter 5, pg. 99.\n\
\n\
The following article made clear to me that the nodes should\n\
contain more than one point (this leads to dramatic speed\n\
improvements for the \"all fixed radius neighbor search\", see\n\
below):\n\
\n\
JL Bentley, \"K-d trees for semidynamic point sets,\" in Sixth Annual\n\
ACM Symposium on Computational Geometry, vol. 91. San Francisco, 1990\n\
\n\
This KD implementation also performs an \"all fixed radius neighbor search\",\n\
i.e. it can find all point pairs in a set that are within a certain radius\n\
of each other. As far as I know the algorithm has not been published.");


static PyTypeObject KDTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "C KDTree",                  /*tp_name*/
    sizeof(KDTree),              /*tp_basicsize*/
    0,                           /*tp_itemsize*/
    (destructor)KDTree_dealloc,  /*tp_dealloc*/
    0,                           /*tp_print*/
    0,                           /*tp_getattr*/
    0,                           /*tp_setattr*/
    0,                           /*tp_compare*/
    0,                           /*tp_repr*/
    0,                           /*tp_as_number*/
    0,                           /*tp_as_sequence*/
    0,                           /*tp_as_mapping*/
    0,                           /*tp_hash */
    0,                           /*tp_call*/
    0,                           /*tp_str*/
    0,                           /*tp_getattro*/
    0,                           /*tp_setattro*/
    0,                           /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,          /*tp_flags*/
    KDTree_doc,                  /* tp_doc */
    0,                           /* tp_traverse */
    0,                           /* tp_clear */
    0,                           /* tp_richcompare */
    0,                           /* tp_weaklistoffset */
    0,                           /* tp_iter */
    0,                           /* tp_iternext */
    KDTree_methods,              /* tp_methods */
    NULL,                        /* tp_members */
    0,                           /* tp_getset */
    0,                           /* tp_base */
    0,                           /* tp_dict */
    0,                           /* tp_descr_get */
    0,                           /* tp_descr_set */
    0,                           /* tp_dictoffset */
    0,                           /* tp_init */
    0,                           /* tp_alloc */
    (newfunc)KDTree_new,         /* tp_new */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

PyDoc_STRVAR(module_doc,
"KDTree implementation for fast neighbor searches in 3D structures.\n\n\
This module implements three objects: KDTree, Point, and Neighbor. Refer \
to their docstrings for more documentation on usage and implementation."
);

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "kdtrees",
        module_doc,
        -1,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit_kdtrees(void)
{
  PyObject *module;

  PointType.tp_new = PyType_GenericNew;
  NeighborType.tp_new = PyType_GenericNew;

  if (PyType_Ready(&KDTreeType) < 0)
      return NULL;
  if (PyType_Ready(&PointType) < 0)
      return NULL;
  if (PyType_Ready(&NeighborType) < 0)
      return NULL;

  module = PyModule_Create(&moduledef);
  if (module == NULL) return NULL;

  Py_INCREF(&KDTreeType);
  if (PyModule_AddObject(module, "KDTree", (PyObject*) &KDTreeType) < 0) {
      Py_DECREF(module);
      Py_DECREF(&KDTreeType);
      return NULL;
  }

  Py_INCREF(&PointType);
  if (PyModule_AddObject(module, "Point", (PyObject*) &PointType) < 0) {
      Py_DECREF(module);
      Py_DECREF(&PointType);
      return NULL;
  }

  Py_INCREF(&NeighborType);
  if (PyModule_AddObject(module, "Neighbor", (PyObject*) &NeighborType) < 0) {
      Py_DECREF(module);
      Py_DECREF(&NeighborType);
      return NULL;
  }

  return module;
}
