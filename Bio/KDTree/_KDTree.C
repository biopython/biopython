#include "_KDTree.h"

float KDTREE_dist(float *coord1, float *coord2, int dim)
{
	// returns the SQUARE of the distance between two points
	int i;
	float sum=0, dif=0;

	for (i=0; i<dim; i++)
	{
		dif=coord1[i]-coord2[i];
		sum+=dif*dif;
	}
	return sum;
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

Node::Node(float cut_value, int cut_dim, long int start, long int end)
{
	_left=NULL;
	_right=NULL;
	_cut_value=cut_value;
	_cut_dim=cut_dim;
	// start and end index in _data_point_list
	_start=start;
	_end=end;
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

long int Node::get_start(void)
{
	return _start;
}

long int Node::get_end(void)
{
	return _end;
}

float Node::get_cut_value(void)
{
	return _cut_value;
}

int Node::get_cut_dim(void)
{
	return _cut_dim;
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

int Node::is_bucket(void)
{
	if (_start==_end+1)
	{
		// Node contains a single point
		return 0;
	}
	else
	{
		// Node contains several points
		return 1;
	}
}
// Region

int Region::dim=3;

Region::Region(float *left, float *right)
{
	_left=new float[Region::dim];
	_right=new float[Region::dim];

	if (left==NULL || right==NULL)
	{
		// [-INF, INF]
		
		int i;

		for (i=0; i<Region::dim; i++)
		{
			_left[i]=-INF;
			_right[i]=INF;
		}
	}
	else
	{
		int i;
		for (i=0; i<Region::dim; i++)
		{
			_left[i]=left[i];
			_right[i]=right[i];
		}
	}
}

Region::~Region()
{
	delete [] _left;
	delete [] _right;
}

Region *Region::intersect_right(float split_coord, int current_dim)
{
	float l, r;

	r=_right[current_dim];
	l=_left[current_dim];

	if (split_coord<=l)
	{
		// split point lies to the left
		return new Region(_left, _right);
	}
	else
	{
		if (split_coord<=r)
		{
			// split point in interval
			// adjust left
			int i;
			float new_left[Region::dim];

			for (i=0; i<Region::dim; i++)
			{
				new_left[i]=_left[i];
			}
			new_left[current_dim]=split_coord;
			return new Region(new_left, _right);
		}
		else
		{
			// interval lies to the left of split point

			return NULL;
		}
	}
}

Region *Region::intersect_left(float split_coord, int current_dim)
{
	float l, r;

	r=_right[current_dim];
	l=_left[current_dim];

	if (split_coord<l)
	{
		// nothing to the left
		return NULL;
	}
	else
	{
		if (split_coord<r)
		{
			// split point in interval
			// adjust right
			int i;
			float new_right[Region::dim];

			for (i=0; i<Region::dim; i++)
			{
				new_right[i]=_right[i];
			}
			new_right[current_dim]=split_coord;
			return new Region(_left, new_right);
		}
		else
		{
			return new Region(_left, _right);
		}
	}
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

int Region::test_intersection(Region *query_region, float radius)
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

		if (ls-rq>radius)
		{
			// outside
			return 0;
		}
		else if (lq-rs>radius)
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

KDTree::KDTree(int dim, int bucket_size)
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
	_bucket_size=bucket_size;
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
	int dim;

	if (depth==0)
	{
		// start with [begin, end+1[
		offset_begin=0;
		offset_end=_data_point_list.size();
		dim=0;
	}
	else
	{
		dim=depth%KDTree::dim;
	}

	if ((offset_end-offset_begin)<=_bucket_size) 
	{
		// leaf node

		return new Node(-1, dim, offset_begin, offset_end);
	}
	else
	{
		long int offset_split;
		long int left_offset_begin, left_offset_end;
		long int right_offset_begin, right_offset_end;
		long int d;
		float cut_value;
		DataPoint data_point;
		Node *left_node, *right_node, *new_node;

		// set sort dimension
		DataPoint::current_dim=dim;

		// sort method sorts [first, last[
		sort(_data_point_list.begin()+offset_begin, _data_point_list.begin()+offset_end);

		// calculate index of split point
		d=offset_end-offset_begin;
		offset_split=d/2+d%2;

		data_point=_data_point_list[offset_begin+offset_split-1];
		cut_value=(data_point.get_coord())[dim];
		
		// create new node and bind to left & right nodes
		new_node=new Node(cut_value, dim, offset_begin, offset_end);

		// left
		left_offset_begin=offset_begin;
		left_offset_end=offset_begin+offset_split;
		left_node=_build_tree(left_offset_begin, left_offset_end, depth+1);

		// right
		right_offset_begin=left_offset_end;
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
		// start with [-INF, INF] region
		
		region=new Region();

		// start with root node
		node=_root;
	}

	current_dim=depth%KDTree::dim;

	if(node->is_leaf())
	{
		long int i;

		for (i=node->get_start(); i<node->get_end(); i++)
		{
			DataPoint data_point;

			data_point=_data_point_list[i];

			if (_query_region->encloses(data_point.get_coord()))
			{
				// point is enclosed in query region - report & stop
				_report_point(data_point.get_index(), data_point.get_coord());
			}
		}
	}
	else
	{
		Node *left_node, *right_node;
		Region *left_region, *right_region;

		left_node=node->get_left_node();

		// LEFT HALF PLANE

		// new region
		left_region=region->intersect_left(node->get_cut_value(), current_dim);

		_test_region(left_node, left_region, depth);

		// RIGHT HALF PLANE

		right_node=node->get_right_node();

		// new region
		right_region=region->intersect_right(node->get_cut_value(), current_dim);

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
		// report point(s)
		long int i;

		for (i=node->get_start(); i<node->get_end(); i++)
		{
			DataPoint data_point;
			data_point=_data_point_list[i];
			_report_point(data_point.get_index(), data_point.get_coord());
		}
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

	if (r<=_radius_sq)
	{
		_index_list.push_back(index);
		// note use of sqrt - only calculated if necessary
		_radius_list.push_back(sqrt(r));
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
}

	
void KDTree::search_center_radius(float *coord, float radius)
{
	int i;
	float left[KDTree::dim], right[KDTree::dim];

	_index_list.clear();
	_radius_list.clear();
	_count=0;

	_radius=radius;
	// use of r^2 to avoid sqrt use
	_radius_sq=radius*radius;

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

void KDTree::neighbor_copy_indices(long int *indices)
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

void KDTree::neighbor_copy_radii(float *radii)
{
	long int i;

	if (_neighbor_count==0)
	{
		return;
	}

	for(i=0; i<_neighbor_count; i++)
	{
		radii[i]=_neighbor_radius_list[i];
	}
}

long int KDTree::neighbor_get_count(void)
{
	return _neighbor_count;
}

void KDTree::neighbor_search(float neighbor_radius)
{
	Region *region;

	_neighbor_index_list.clear();
	_neighbor_radius_list.clear();
	// note the use of r^2 to avoid use of sqrt
	_neighbor_radius=neighbor_radius;
	_neighbor_radius_sq=neighbor_radius*neighbor_radius;
	_neighbor_count=0;

	// start with [-INF, INF]
	region=new Region();

	if (_root->is_leaf())
	{
		// this is a boundary condition
		// bucket_size>nr of points
		_search_neighbors_in_bucket(_root);
	}
	else
	{
		// "normal" situation
		_neighbor_search(_root, region, 0);
	}
	delete region;
}

void KDTree::_neighbor_search(Node *node, Region *region, int depth)
{
	Node *left, *right;
	Region *left_region, *right_region; 
	int dim;

	dim=depth%KDTree::dim;

	left=node->get_left_node();
	right=node->get_right_node();

	// planes of left and right nodes
	left_region=region->intersect_left(node->get_cut_value(), dim);
	right_region=region->intersect_right(node->get_cut_value(), dim);

	if (!left->is_leaf())
	{
		// search for pairs in this half plane
		_neighbor_search(left, left_region, depth+1);
	}
	else
	{
		_search_neighbors_in_bucket(left);
	}
	
	if (!right->is_leaf())
	{
		// search for pairs in this half plane
		_neighbor_search(right, right_region, depth+1);
	}
	else
	{
		_search_neighbors_in_bucket(right);
	}

	// search for pairs between the half planes
	_neighbor_search_pairs(left, left_region, right, right_region, depth+1);

	// cleanup
	delete left_region;
	delete right_region;
}

void KDTree::_test_neighbors(DataPoint &p1, DataPoint &p2)
{
	float r;

	r=KDTREE_dist(p1.get_coord(), p2.get_coord(), KDTree::dim);

	if(r<=_neighbor_radius_sq)
	{
		// we found a neighbor pair!
		_neighbor_index_list.push_back(p1.get_index());
		_neighbor_index_list.push_back(p2.get_index());
		// note sqrt
		_neighbor_radius_list.push_back(sqrt(r));
		_neighbor_count++;
	}
}

void KDTree::_search_neighbors_in_bucket(Node *node)
{
	long int i;

	for(i=node->get_start(); i<node->get_end(); i++)
	{
		DataPoint p1;
		long int j;

		p1=_data_point_list[i];

		for (j=i+1; j<node->get_end(); j++)
		{
			DataPoint p2;

			p2=_data_point_list[j];

			_test_neighbors(p1, p2);
		}
	}
}

void KDTree::_search_neighbors_between_buckets(Node *node1, Node *node2)
{
	long int i;

	for(i=node1->get_start(); i<node1->get_end(); i++)
	{
		DataPoint p1;
		long int j;

		p1=_data_point_list[i];

		for (j=node2->get_start(); j<node2->get_end(); j++)
		{
			DataPoint p2;

			p2=_data_point_list[j];

			_test_neighbors(p1, p2);
		}
	}
}

void KDTree::_neighbor_search_pairs(Node *down, Region *down_region, 
		Node *up, Region *up_region, int depth)
{
	int down_is_leaf, up_is_leaf;
	int dim;

	// if regions do not overlap - STOP
	if (!down || !up || !down_region || !up_region) 
	{
		// STOP
		return;
	}
	
	if (down_region->test_intersection(up_region, _neighbor_radius)==0)
	{
		// regions cannot contain neighbors
		return;
	}

	// dim
	dim=depth%KDTree::dim;

	// are they leaves?
	up_is_leaf=up->is_leaf();
	down_is_leaf=down->is_leaf();

	if (up_is_leaf && down_is_leaf)
	{
		// two leaf nodes
		_search_neighbors_between_buckets(down, up);
	}
	else
	{
		// one or no leaf nodes

		Node *up_right, *up_left, *down_left, *down_right;
		Region *up_left_region, *up_right_region, 
			*down_left_region, *down_right_region;  
		
		if (down_is_leaf)
		{
			down_left=down;
			// make a copy of down_region
			down_left_region=new Region(down_region->get_left(), down_region->get_right());
			down_right=NULL;
			down_right_region=NULL;
		}
		else
		{
			float cut_value;

			cut_value=down->get_cut_value();

			down_left=down->get_left_node();
			down_right=down->get_right_node();
			down_left_region=down_region->intersect_left(cut_value, dim); 
			down_right_region=down_region->intersect_right(cut_value, dim); 
		}

		if (up_is_leaf)
		{
			up_left=up;
			// make a copy of up_region
			up_left_region=new Region(up_region->get_left(), up_region->get_right());
			up_right=NULL;
			up_right_region=NULL;
		}
		else
		{
			float cut_value;

			cut_value=up->get_cut_value();

			up_left=up->get_left_node();
			up_right=up->get_right_node();
			up_left_region=up_region->intersect_left(cut_value, dim); 
			up_right_region=up_region->intersect_right(cut_value, dim); 
		}

		_neighbor_search_pairs(up_left, up_left_region, down_left, down_left_region, depth+1);
		_neighbor_search_pairs(up_left, up_left_region, down_right, down_right_region, depth+1);
		_neighbor_search_pairs(up_right, up_right_region, down_left, down_left_region, depth+1);
		_neighbor_search_pairs(up_right, up_right_region, down_right, down_right_region, depth+1);

		delete down_left_region;
		delete down_right_region;
		delete up_left_region;
		delete up_right_region;
	}
}

void KDTree::neighbor_simple_search(float radius)
{
	long int i;

	_neighbor_radius=radius;
	_neighbor_radius_sq=radius*radius;

	_neighbor_count=0;
	_neighbor_index_list.clear();
	_neighbor_radius_list.clear();

	DataPoint::current_dim=0;
	sort(_data_point_list.begin(), _data_point_list.end());

	for (i=0; i<_data_point_list.size(); i++)
	{
		float x1;
		long int j;
		DataPoint p1;

		p1=_data_point_list[i];
		x1=p1.get_coord()[0];

		for (j=i+1; j<_data_point_list.size(); j++)
		{
			DataPoint p2;
			float x2;

			p2=_data_point_list[j];
			x2=p2.get_coord()[0];

			if (fabs(x2-x1)<=radius)
			{
				_test_neighbors(p1, p2);
			}
			else
			{
				break;
			}
		}
	}
}



