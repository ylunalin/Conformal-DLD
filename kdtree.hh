#ifndef KDTREE_HH
#define KDTREE_HH
#include "particle.hh"
#include <vector>
#include <cmath>
#include <algorithm>

/** Define grid structure */
struct Grid {
	/**total number of blocks along one axis*/
	int nb;
	/** Lower and Upper bound of the simulation grid */
	double rad;
	/** grid unit size */
	double dx;
	/** counter array in a dld simulation*/
	int *co;
	/** particle storage in  a dld simulation*/
	particle ** ba;
	Grid(){}
	Grid(double rad_, double dx_):
	nb(1), rad(rad_), dx(dx_) {}

	Grid(int nb_, int *co_, double rad_, double dx_, particle ** ba_) :
	nb(nb_), rad(rad_), dx(dx_), co(co_), ba(ba_) {}
};

/** Define a node by its coordinate, as well as pointers to its daughter nodes.*/
struct Node {
	bool type;//1 = non-leaf node, 0 = leaf-node
	//block number on the grid
	int blc;
	//children nodes, upper left, upper right, lower right, lower left
	int cl[4];
	double min[2];
	double max[2];
	double mid[2];
	Node():
	type(0), blc(0)
	{
		for(int i=0; i<4;i++) cl[i]=0;
		for(int i=0; i<2;i++) {min[i]=0; max[i]=0; mid[i]=0;}
	}
};

class KDTree{
	public:
	KDTree(){}
	KDTree(Grid g);
	particle* closest(double x, double y, double &min_d, double &sec_d);

	//private:
	/** Dimension */
	static const int k=2;
	/** Minimal number of particles in each leaf*/
	static const int min_np = 10;
	struct Nextkid{
		int kid_id;
		double dsq;
		Nextkid(){}
		Nextkid(int id, double d):
		kid_id(id), dsq(d){}
	};
	// compare two distances
	struct sort_index{
		std::vector<double> dist;
		bool operator()(int i, int j) {
			return (dist[i] < dist[j]);
		}
		sort_index(std::vector<double> dist_):
			dist(dist_) {}
	};

	unsigned int search_table;
	particle * nearest;
	std::vector<particle *> particles;
	std::vector<Node> nodes;
	Grid g;
	double dsqr(double ax, double ay, double bx, double by);
	double pt_to_node(double x, double y, Node &target);
	unsigned int make_table();
	void flatten();
	void custom_sort(int a, int b, int &min, int axis);
	bool has_minimum(int n, int m);
	bool has_minimum(double *min, double *max);
	void grow_tree(int node_id, double *min, double *max);
	//void grow_tree(int a, int b, int depth);
	int * sort_kids(Nextkid *kids);
	void search_next(double x, double y, Node &n, int &i, int &j, int &k, int &l);
	void search(double x, double y, int node_id, double &min_d, double &sec_d);
};
#endif
