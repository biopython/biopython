#include "_KDTree.h"

float KDTREE_dist(float *coord1, float *coord2, int dim)
{
	int i;
	float sum=0, dif=0;

	for (i=0; i<dim; i++)
	{
		dif=coord1[i]-coord2[i];
		sum+=dif*dif;
	}
	return sqrt(sum);
}


// DataPoint

int DataPoint::current_dim=0;
int DataPoint::dim=3;

int operator<(const DataPoint &self, const DataPoint &other)
{
	float a, b;

	a=self._coord[DataPoint::current_dim];
	b=other._coord[DataPoint::current_dim];

	return a<b; 
}

int operator==(const DataPoint &self, const DataPoint &other)
{
	float a, b;

	a=self._coord[DataPoint::current_dim];
	b=other._coord[DataPoint::current_dim];
	
	return a==b;
}

void DataPoint::set_data(long int index, float *coord)
{
	_index=index;
	_coord=coord;
}

float *DataPoint::get_coord(void)
{
	return _coord;
}

long int DataPoint::get_index(void)
{
	return _index;
}

// Node

Node::Node(DataPoint &data_point)
{
	_left=NULL;
	_right=NULL;
	_coord=data_point.get_coord();
	_index=data_point.get_index();
}

Node::~Node()
{
	delete _left;
	delete _right;
}

void Node::set_left_node(Node *node)
{
	_left=node;
}

void Node::set_right_node(Node *node)
{
	_right=node;
}

Node *Node::get_left_node(void)
{
	return _left;
}

Node *Node::get_right_node(void)
{
	return _right;
}

long int Node::get_index(void)
{
	return _index;
}

float *Node::get_coord(void)
{
	return _coord;
}

int Node::is_leaf(void)
{
	if (_left==NULL && _right==NULL)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Region

int Region::dim=3;

Region::Region(float *left, float *right)
{
	_left=new float[Region::dim];
	_right=new float[Region::dim];

	int i;
	for (i=0; i<Region::dim; i++)
	{
		_left[i]=left[i];
		_right[i]=right[i];
	}
}

Region::~Region()
{
	delete [] _left;
	delete [] _right;
}

Region *Region::intersect_right(float *split_coord, int current_dim)
{
	float new_left[Region::dim];
	int i;

	for (i=0; i<Region::dim; i++)
	{
		new_left[i]=_left[i];
	}

	new_left[current_dim]=split_coord[current_dim];

	return new Region(new_left, _right);
}

Region *Region::intersect_left(float *split_coord, int current_dim)
{
	float new_right[Region::dim];
	int i;

	for (i=0; i<Region::dim; i++)
	{
		new_right[i]=_right[i];
	}

	new_right[current_dim]=split_coord[current_dim];

	return new Region(_left, new_right);
}

int Region::encloses(float *coord)
{
	int i;
	for (i=0; i<Region::dim; i++)
	{
		if (!(coord[i]>=_left[i] && coord[i]<=_right[i]))
		{
			return 0;
		}
	}
	return 1;
}

float *Region::get_left(void)
{
	return _left;
}

float *Region::get_right(void)
{
	return _right;
}

int Region::test_intersection(Region *query_region)
{
	int status=2;

	int i;
	for (i=0; i<Region::dim; i++)
	{
		float rs, rq, ls, lq;

		rs=_right[i];
		ls=_left[i];
		rq=query_region->get_right()[i];
		lq=query_region->get_left()[i];

		if (ls>rq)
		{
			// outside
			return 0;
		}
		else if (rs<lq)
		{
			// outside
			return 0;
		}
		else if (rs<=rq && ls>=lq)
		{
			// inside (at least in dim i)
			status=min(status, 2);
		}
		else
		{
			// overlap (at least in dim i)    
			status=1;
		}
	}
	return status;
}

// KDTree

int KDTree::dim=3;

KDTree::KDTree(int dim)
{
	// set dimension
	KDTree::dim=dim;
	DataPoint::dim=dim;
	Region::dim=dim;

	_center_coord=new float[KDTree::dim];
	_query_region=NULL;
	_root=NULL;
	_coords=NULL;
	_count=0;
	_neighbor_count=0;
}

KDTree::~KDTree()
{
	// clean up KD tree
	delete _root;
	delete _query_region;
	delete [] _center_coord;
	delete [] _coords;
}

Node *KDTree::_build_tree(long int offset_begin, long int offset_end, int depth)
{

	if (depth==0)
	{
		// start with [begin, end+1[
		offset_begin=0;
		offset_end=_data_point_list.size();
	}

	if ((offset_end-offset_begin)==1) 
	{
		// leaf node
		return new Node(_data_point_list[offset_begin]);
	}
	else
	{
		long int offset_split;
		long int left_offset_begin, left_offset_end;
		long int right_offset_begin, right_offset_end;
		long int d;

		Node *left_node, *right_node, *new_node;

		// set sort dimension
		DataPoint::current_dim=depth%KDTree::dim;

		sort(_data_point_list.begin()+offset_begin, _data_point_list.begin()+offset_end);

		// calculate index of split point
		d=offset_end-offset_begin;
		offset_split=d/2+d%2;

		
		// create new node and bind to left & right nodes
		new_node=new Node(_data_point_list[offset_begin+offset_split-1]);

		// left
		left_offset_begin=offset_begin;
		left_offset_end=offset_begin+offset_split;
		left_node=_build_tree(left_offset_begin, left_offset_end, depth+1);

		// right
		right_offset_begin=offset_begin+offset_split;
		right_offset_end=offset_end;
		right_node=_build_tree(right_offset_begin, right_offset_end, depth+1);

		new_node->set_left_node(left_node);
		new_node->set_right_node(right_node);

		return new_node;
	}
}

void KDTree::_add_point(long int index, float *coord)
{
	DataPoint data_point;

	data_point.set_data(index, coord);
	// add to list of points
	_data_point_list.push_back(data_point);
}

void KDTree::_set_query_region(float *left, float *right)
{
	delete _query_region;
	_query_region=new Region(left, right);
}

void KDTree::_search(Region *region, Node *node, int depth)
{
	int current_dim;

	if(depth==0)
	{
		float left[KDTree::dim];
		float right[KDTree::dim];

		// start with [-INF, INF] region
		int i;
		for (i=0; i<KDTree::dim; i++)
		{
			left[i]=-INF;
			right[i]=INF;
		}
		
		region=new Region(left, right);

		// start with root node
		node=_root;
	}

	current_dim=depth%KDTree::dim;

	if(node->is_leaf())
	{
		if (_query_region->encloses(node->get_coord()))
		{
			// point is enclosed in query region - report & stop
			_report_point(node->get_index(), node->get_coord());
		}
	}
	else
	{
		Node *left_node, *right_node;
		Region *left_region, *right_region;

		left_node=node->get_left_node();

		// LEFT HALF PLANE

		// new region
		left_region=region->intersect_left(node->get_coord(), current_dim);

		_test_region(left_node, left_region, depth);

		// RIGHT HALF PLANE

		right_node=node->get_right_node();

		// new region
		right_region=region->intersect_right(node->get_coord(), current_dim);

		// test for overlap/inside/outside & do recursion/report/stop
		_test_region(right_node, right_region, depth);
	}

	delete region;
}

void KDTree::_test_region(Node *node, Region *region, int depth)
{
	int intersect_flag;

	// is node region inside, outside or overlapping 
	// with query region?
	intersect_flag=region->test_intersection(_query_region);

	if (intersect_flag==2)				
	{
		// inside - extract points
		
		_report_subtree(node);
		// end of recursion
		// get rid of region
		delete region;
	}
	else if (intersect_flag==1)			
	{
		// overlap - recursion
		
		_search(region, node, depth+1);	
		// search does cleanup of region
	}
	else
	{
		// outside - stop
		
		// end of recursion
		// get rid of region
		delete region;	
	}
}

void KDTree::_report_subtree(Node *node)
{
	if (node->is_leaf())
	{
		// report point
		_report_point(node->get_index(), node->get_coord());
	}
	else
	{
		// find points in subtrees via recursion 
		_report_subtree(node->get_left_node());
		_report_subtree(node->get_right_node());
	}
}

void KDTree::_report_point(long int index, float *coord)
{
	float r;

	r=KDTREE_dist(_center_coord, coord, KDTree::dim);

	if (r<=_radius)
	{
		_index_list.push_back(index);
		_radius_list.push_back(r);
		_count++;
	}
}
		
void KDTree::set_data(float *coords, long int nr_points)
{
	long int i;

	// clean up stuff from previous use
	delete _root;
	delete [] _coords;
	_index_list.clear();
	_radius_list.clear();
	_count=0;
	// keep pointer to coords to delete it
	_coords=coords;
		
	for (i=0; i<nr_points; i++)
	{
		_add_point(i, coords+i*KDTree::dim);
	}

	// build KD tree
	_root=_build_tree();
	// we can get rid of the vector now
	_data_point_list.clear();
}
	
void KDTree::search_center_radius(float *coord, float radius)
{
	int i;
	float left[KDTree::dim], right[KDTree::dim];

	_index_list.clear();
	_radius_list.clear();
	_count=0;

	_radius=radius;

	for (i=0; i<KDTree::dim; i++)
	{
		left[i]=coord[i]-radius;
		right[i]=coord[i]+radius;
		// set center of query
		_center_coord[i]=coord[i];
	}

	// clean up! 
	delete [] coord;

	_set_query_region(left, right);

	_search();
}

long int KDTree::get_count(void)
{
	return _count;
}

void KDTree::copy_indices(long *indices)
{
	long int i;

	if (_count==0)
	{
		return;
	}

	for(i=0; i<_count; i++)
	{
		indices[i]=_index_list[i];
	}
}

void KDTree::copy_radii(float *radii)
{
	long int i;

	if (_count==0)
	{
		return;
	}

	for(i=0; i<_count; i++)
	{
		radii[i]=_radius_list[i];
	}
}

void KDTree::copy_neighbors(long int *indices)
{
	long int i;

	if (_neighbor_count==0)
	{
		return;
	}

	for(i=0; i<_neighbor_count*2; i++)
	{
		indices[i]=_neighbor_index_list[i];
	}
}

long int KDTree::get_neighbor_count(void)
{
	return _neighbor_count;
}

void KDTree::neighbor_search(float neighbor_radius)
{
	_neighbor_index_list.clear();
	_neighbor_radius_list.clear();
	_neighbor_radius=neighbor_radius;
	_neighbor_count=0;

	_neighbor_search(_root, 0);
}

void KDTree::_neighbor_search(Node *root, int depth)
{
	Node *left, *right;

	left=root->get_left_node();
	right=root->get_right_node();

	if (!left->is_leaf())
	{
		// search for pairs in this half plane
		_neighbor_search(left, depth+1);
	}
	
	if (!right->is_leaf())
	{
		// search for pairs in this half plane
		_neighbor_search(right, depth+1);
	}

	// search for pairs between the half planes
	_neighbor_search_pairs(left, right, depth+1);
}

void KDTree::_neighbor_search_pairs(Node *left, Node *right, int depth)
{
	int left_is_leaf, right_is_leaf;
	int dim, near=0;
	float xl, xr, d;

	// dim
	dim=depth%KDTree::dim;

	// are the nodes within radius along dim?
	xl=(left->get_coord())[dim];
	xr=(right->get_coord())[dim];

	// difference along dim
	d=xr-xl;

	if ((fabsf(d)<=_neighbor_radius))
	{
		near=1;
	}

	// order nodes: left should be <= right in this dim
	if (d<0)
	{
		// right<left in this dimension so switch nodes
		Node *n;
		n=right;
		right=left;
		left=n;
	}

	// are they leaves?
	left_is_leaf=left->is_leaf();
	right_is_leaf=right->is_leaf();

	if (left_is_leaf && right_is_leaf)
	{
		// two leaf nodes
		
		// test if two points are within radius 
		if (near==1)
		{
			// the two leaves MIGHT be within _neighbor_radius
			
			float r;
			long int li, ri;

			li=left->get_index();
			ri=right->get_index();

			r=KDTREE_dist(left->get_coord(), right->get_coord(), KDTree::dim);

			if(r<=_neighbor_radius)
			{
				// we found a neighbor pair!
				_neighbor_index_list.push_back(li);
				_neighbor_index_list.push_back(ri);
				_neighbor_count++;
			}
		}
		else
		{
			// the two leaves are not within _neighbor_radius
			return;
		}
	}
	else
	{
		// one or no leaf nodes
		
		if (left_is_leaf)
		{
			if (near)
			{
				_neighbor_search_pairs(left, right->get_right_node(), depth+1);
				_neighbor_search_pairs(left, right->get_left_node(), depth+1);
			}
			else
			{
				_neighbor_search_pairs(left, right->get_left_node(), depth+1);
			}
		}
		else if (right_is_leaf)
		{
			if (near)
			{
				_neighbor_search_pairs(left->get_right_node(), right, depth+1);
				_neighbor_search_pairs(left->get_left_node(), right, depth+1);
			}
			else
			{
				_neighbor_search_pairs(left->get_right_node(), right, depth+1);
			}
		}
		else
		{
			if (near)
			{
				// test all 4 combinations
				_neighbor_search_pairs(left->get_left_node(), right->get_left_node(), depth+1);
				_neighbor_search_pairs(left->get_left_node(), right->get_right_node(), depth+1);
				_neighbor_search_pairs(left->get_right_node(), right->get_left_node(), depth+1);
				_neighbor_search_pairs(left->get_right_node(), right->get_right_node(), depth+1);
			}
			else
			{
				// only test 3 "overlapping" nodes
				_neighbor_search_pairs(left->get_left_node(), right->get_left_node(), depth+1);
				// NOT this one since they are >R apart, so this is were we save time!
				//_neighbor_search_pairs(left->get_left_node(), right->get_right_node(), depth+1);
				_neighbor_search_pairs(left->get_right_node(), right->get_left_node(), depth+1);
				_neighbor_search_pairs(left->get_right_node(), right->get_right_node(), depth+1);
			}
		}
	}
	
}


