#include "kdtree.hh"
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

KDTree::KDTree(Grid g_):
g(g_)
{
	double min[2];
	double max[2];
	min[0]=min[1]=(-1)*g.rad;
	max[0]=max[1]=g.rad;
	nodes.push_back(Node());
	grow_tree(0,min, max);
	search_table=make_table();
}

unsigned int KDTree::make_table(){
	unsigned int a =0;
	a = (a<<2)+1;
	a = (a<<2)+3;
	a = (a<<2)+2;

	a = (a<<2)+1;
	a = (a<<2);
	a = (a<<2)+2;
	a = (a<<2)+3;

	a = (a<<2)+2;
	a = (a<<2)+1;
	a = (a<<2)+3;
	a = (a<<2);


	a = (a<<2)+3;
	a = (a<<2);
	a = (a<<2)+2;
	a = (a<<2)+1;

	return a;
}

/* Utility to calcuate squared distance*/
double KDTree::dsqr(double ax, double ay, double bx, double by){
	return (bx-ax)*(bx-ax)+(by-ay)*(by-ay);
}

/* Return minimal distance from a point to a node, return zero if point is in the node*/
double KDTree::pt_to_node(double x, double y, Node &target){
	double dist =0;
	double minx = target.min[0];
	double miny = target.min[1];
	double maxx = target.max[0];
	double maxy = target.max[1];
	if(y>maxy){
		if (x>maxx) dist = dsqr(x,y,maxx, maxy);
		else if( x<=maxx && x>=minx) dist = (y-maxy)*(y-maxy);
		else dist = dsqr(x,y,minx, maxy);
	}
	else if (y<miny){
		if(x>maxx) dist = dsqr(x,y,maxx, miny);
		else if( x<=maxx && x>=minx) dist = (miny-y)*(miny-y);
		else dist = dsqr(x,y,minx, miny);
	}
	else {
		if(x<=maxx && x>=minx) dist=0;
		else if (x>maxx) dist= (x-maxx)*(x-maxx);
		else dist = (minx-x)*(minx-x);
	}
	return dist;
}

/* Create a 1D list of particles from member **ba in a dld simulation*/
void KDTree::flatten(){
	int tot_blc = g.nb*g.nb;
	//iterate over all particles, and put their addresses in the list
	for(int i=0; i<tot_blc;i++)
		for(int j=0; j<g.co[i]; j++)
			particles.push_back(&(g.ba[i][j]));
}

/* A sort on this list of particles along a given axis to find a median particle, report median location use nth_element in algorithm package*/
void KDTree::custom_sort(int a, int b, int &mid, int axis){
	//insert code
}


/* Tally the current section of the vector particles, if it gets to more than 3, break and set hm=true, else do nothing. */
bool KDTree::has_minimum(int n, int m){
	bool hm =false;
	//insert code
	return hm;
}

/* Check whether a node in quadtree has minmal dimension*/
bool KDTree::has_minimum(double *min_, double *max_){
	bool hm = false;
	double diff_x = max_[0] - min_[0];
	double diff_y = max_[1] - min_[1];
	if(diff_x > g.dx && diff_y > g.dx) hm = true;
	return hm;
}

/* Build a quadtree, taking advantage of dld simulation square Grid structure.
 dld_sim makes sure that m, n are even number */
void KDTree::grow_tree(int node_id, double *min, double *max){
	//printf("node id=%d\n", node_id);
	//printf("size of vector<Node> nodes %d\n", int (nodes.size()));
	Node &current = nodes[node_id];
	current.min[0] =min[0];
	current.min[1] =min[1];
	current.max[0] =max[0];
	current.max[1] =max[1];
	int nx = int ((min[0]+g.rad)/g.dx);if(nx<0) nx=0;if(nx>=g.nb) nx=g.nb-1;
	int ny = int ((min[1]+g.rad)/g.dx);if(ny<0) ny=0;if(ny>=g.nb) ny=g.nb-1;
	current.blc = nx+g.nb*ny;
	if (has_minimum(min, max)) {
		current.type = 1;
		//generate the bounding boxes

		double midx = 0.5*(max[0]+min[0]);
		double midy = 0.5*(max[1]+min[1]);
		current.mid[0] = midx;
		current.mid[1] = midy;
		//upper left
		double c0min[2];
		c0min[0]=min[0];
		c0min[1]=midy;
		double c0max[2];
		c0max[0] =midx;
		c0max[1] =max[1];
		//upper right
		double c1min[2];
		c1min[0] = midx; c1min[1] = midy;
		double c1max[2];
		c1max[0] =max[0]; c1max[1]=max[1];
		//lower right
		double c2min[2];
		c2min[0] = midx;
		c2min[1] = min[1];
		double c2max[2];
		c2max[0] =max[0];
		c2max[1] =midy;
		//lower left
		double c3min[2];
		c3min[0] =min[0]; c3min[1]=min[1];
		double c3max[2];
		c3max[0]=midx;
		c3max[1]=midy;

		int ul = current.cl[0]=int (nodes.size());
		int ur = current.cl[1]=int (nodes.size()) +1;
		int lr = current.cl[2]=int (nodes.size()) +2;
		int ll = current.cl[3]=int (nodes.size()) +3;
		//printf("ul, ur ll, lr = %d %d %d %d \n",ul,ur,ll,lr);
		nodes.push_back(Node());
		nodes.push_back(Node());
		nodes.push_back(Node());
		nodes.push_back(Node());
		grow_tree(ul, c0min, c0max);
		grow_tree(ur, c1min, c1max);
		grow_tree(lr, c2min, c2max);
		grow_tree(ll, c3min, c3max);
		//delete c0min; delete c0max; delete c1min; delete c1max;
		//delete c2min; delete c2max; delete c3min; delete c3max;

	}
	else {
		current.type = 0;
		return;
	}
}

 /* Build a more generic (k=2)dtree, using a list of particles. Later can change that to a template.
void KDTree::grow_tree(int a, int b, int depth){
	//Make sure the axis to split on alternates with each split
	int axis = depth % k;
	int m=0;
	// Sort and find the median of ba_flat, along a given axis
	custom_sort(a,b,m, axis);
	// insert code to determine the min[2] max[2] of the current node
	// basically finding the bouding box of the current set of points
	// set the children node index
	// then pass the portion of the particle array to them for them to split
	//nodes.push_back(current);
	if(has_minimum(a, m)) grow_tree(a, m, depth+1);
	if(has_minimum(m, b)) grow_tree(m, b, depth+1);

}*/

int * KDTree::sort_kids(Nextkid *kids){
	int *kids_in_order = new int[4];
	for(int i=1;i<4;i++){
		int insert=i;
		for(int j=i-1; j>=0; j--){
			if(kids[i].dsq<kids[j].dsq) insert = j;
		}
		Nextkid tmp = kids[i];
		for(int k=i; k>insert; k--) {
			kids[k] = kids[k-1];
		}
		kids[insert]=tmp;
	}
	for(int i=0;i<4;i++) kids_in_order[i] = kids[i].kid_id;
	return kids_in_order;
}

void KDTree::search_next(double x, double y, Node &root, int &fst, int &sec, int &thd, int &forth){
	double midx = 0.5*(root.min[0]+root.max[0]);
	double midy = 0.5*(root.min[1]+root.max[1]);
	if(x >= midx and y >= midy) {
		fst=1; forth=3;
		if(pt_to_node(x, y, nodes[root.cl[0]]) < pt_to_node(x, y, nodes[root.cl[2]])) {sec=0; thd=2;}
		else {sec=2;thd=0;}
	}
	else if (x >=midx and y<midy){
		fst=2; forth=0;
		if(pt_to_node(x,y,nodes[root.cl[1]]) < pt_to_node(x,y,nodes[root.cl[3]])) {sec=1; thd=3;}
		else {sec=3; thd=1;}
	}
	else if (x <midx and y>= midy){
		fst=0; forth=2;
		if(pt_to_node(x,y,nodes[root.cl[1]]) < pt_to_node(x,y,nodes[root.cl[3]])) {sec=1; thd=3;}
		else {sec=3; thd=1;}
	}
	else {
		fst=3; forth=1;
		if(pt_to_node(x, y, nodes[root.cl[0]]) < pt_to_node(x, y, nodes[root.cl[2]])) {sec=0; thd=2;}
					else {sec=2;thd=0;}
	}
}
/* perform a search to the next nearest node */
void KDTree::search(double x, double y, int node_id, double &min_d, double &sec_d){
	//distance to node
	Node& root = nodes[node_id];
	//printf("searching node id (%g %g) (%g %g), node type =%d\n", root.min[0], root.min[1], root.max[0], root.max[1], root.type);
	double d = pt_to_node(x,y,root);
	if(sec_d <d) return;
	else {
		// check if it's a leaf node (type = 1 means it's not a leaf node)
		if(root.type){
			//search 4 children in ascending order in distance
			//using sort with custom sort object
			/*
			std::vector<double> dist;
			std::vector<int> index;
			for(int i=0; i<4; i++){
				index.push_back(i);
				dist.push_back(pt_to_node(x, y, nodes[root.cl[i]]));
			}
			std::sort(index.begin(), index.end(), sort_index(dist));
			for(int i =0; i<4; i++){
				search(x, y, root.cl[index[i]], min_d, sec_d);
			}
			*/

			/*
			//Using my custom sorter
			//Ran out of memory for cluster radius 1400
			Nextkid *kids = new Nextkid[4];
			for (int i=0; i<4; i++) kids[i] = Nextkid(root.cl[i],pt_to_node(x, y, nodes[root.cl[i]]));
			int *kio = sort_kids(kids); //note that kio here stores the actual indeces of children node, not the rank of the children
			delete [] kids;
			for (int i=0; i<4; i++) search(x, y, kio[i],min_d, sec_d);
			*/

			/*
			//Using search_next, hybrid between a table and 2 pt_to_node calculations per node
			int kio []={0,0,0,0};//kio = kids in order, rank of the children
			search_next(x,y,root,kio[0], kio[1],kio[2],kio[3]);
			for (int i=0; i<4; i++) search(x, y, root.cl[kio[i]],min_d, sec_d);
			*/
			int order;
			double a =root.mid[0]; double b=root.mid[1];
			if(x<a and y>=b) order =255&(search_table>>24);
			else if(x>=a and y>=b) order=255&(search_table>>16);
			else if(x>=a and y<b) order=255&(search_table>>8);
			else order=255&search_table;
			for(int i=0; i<4; i++){
				int id=3&(order>>(6-2*i));
				search(x,y, root.cl[id], min_d, sec_d);
			}
		}
		else{
			int s=root.blc;
			for(int i=0; i<g.co[s]; i++){
				particle *bb=&g.ba[s][i];
				double rsq=dsqr(x,y,bb->x, bb->y);
				if(rsq<min_d){
					sec_d=min_d;
					min_d=rsq;
					nearest = bb;
				}
				else {
					if(rsq<sec_d) sec_d=rsq;
				}
			}
			//printf("%g, %g\n", min_d, sec_d);
		}
	}
}
particle* KDTree::closest(double x, double y, double &min_d, double &sec_d){
	if (nodes.empty()) exit(1);
	search(x, y, 0, min_d, sec_d);
	return nearest;
}
