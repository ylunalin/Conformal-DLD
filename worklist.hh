#ifndef DLD_SIM_WORKLIST_HH
#define DLD_SIM_WORKLIST_HH

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

//add template of the min, max, min_d type
template <typename T>
struct Box{
	//Type 0 box hasn't been marked or searched
	//Type 1, 2, 3... box has been marked to be searched number correspond to the index of thresholds, starting from 1
	int type;
	T min[2];
	T max[2];
	T min_d;

	Box (T *bmin, T *bmax):
	type(0){
		min[0]=bmin[0];
		min[1]=bmin[1];
		max[0]=bmax[0];
		max[1]=bmax[1];
	}
	Box():
	type(0){
		min[0]=0;
		min[1]=0;
		max[0]=0;
		max[1]=0;
	}
};

//Since kdtree and dld_sim both have grid struct, use SGrid to distinguish to denote kdtree's
//Grid implement Box template struct with integer
struct SGrid{
	int *min;
	int *max;
	std::vector< Box<int> > boxes;
	SGrid (int *gmin, int *gmax):
		min(gmin), max(gmax)
	{
		int bmin[2], bmax[2];
		for(int j=min[1]; j<max[1]; j++){
			bmin[1] = j; bmax[1] = j+1;
			for(int i=min[0]; i<max[0]; i++){
				bmin[0]=i; bmax[0]=i+1;
				Box<int> dummy(bmin,bmax);
				boxes.push_back(dummy);
			}
		}
	}
};

class wl_gen{
	public:
	std::vector<unsigned long int> worklist;
	wl_gen(int *min, int *max, int* grid_center);
	wl_gen(int *min, int *max, std::vector<int> thr);
	wl_gen(int *min, int *max, int max_thr, int ncell);

	private:
	/*member data*/
	SGrid g;
	int ncell;
	//after user provide a double value of maximum search threshold
	//it's converted to an integer, representing the threshold distance squared
	//in unit of the cell length squared,
	int max_thr;
	Box<int> walker;
	std::vector< Box<double> > walker_cells;
	std::vector<int> threshold;
	std::vector<int> user_def_thr;

	/*functions*/
	void reset();

	template <typename T>
	T dsq(T x1, T y1, T x2, T y2);

	template <typename T>
	void box_dist(Box<T> b, bool subgrid);

	void est_threshold();
	void est_threshold(int *min, int *max, int index);
	void encoder(bool usr_def_thr, bool subgrid, int nwl);

	//A custom sort object to sort the boxes first by type,
	//within the same type, by distance to the grid center
	struct weighted_sort{
		int grid_center[2];
		SGrid *grid;
		Box<int> *walker;
		bool operator()(Box<int> b1, Box<int> b2){
			//Calculate the dist from a box to walker, weighted by the distance from walker to grid center
			double weight = 1.0 / ((grid->min[0]-grid->max[0])*(grid->min[0]-grid->max[0])+(grid->min[1]-grid->max[1])*(grid->min[1]-grid->max[1]));
			double pt1[2] ={0.5*(double(b1.min[0]) +double(b1.max[0])),0.5*(double(b1.min[1]) +double(b1.max[1])) };
			double pt2[2] = {0.5*(double(b2.min[0]) +double(b2.max[0])),0.5*(double(b2.min[1]) +double(b2.max[1])) };
			double pt_walker[2] = {0.5*(double(walker->min[0]) +double(walker->max[0])),0.5*(double(walker->min[1]) +double(walker->max[1])) };
			double d1 = (pt1[0]-grid_center[0])*(pt1[0]-grid_center[0])+ (pt1[1]- grid_center[1])*(pt1[1]- grid_center[1]);
			double d2 = (pt2[0]-grid_center[0])*(pt2[0]-grid_center[0])+ (pt2[1]- grid_center[1])*(pt2[1]- grid_center[1]);
			d1 *=weight; d2 *=weight;
			d1 +=(pt1[0]-pt_walker[0])*(pt1[0]-pt_walker[0])+(pt1[1]-pt_walker[1])*(pt1[1]-pt_walker[1]);
			d2 +=(pt2[0]-pt_walker[0])*(pt2[0]-pt_walker[0])+(pt2[1]-pt_walker[1])*(pt2[1]-pt_walker[1]);
			bool comp=false;
			if(b1.type<b2.type) comp=true;
			else if  (b1.type==b2.type) {
				if(d1<d2) comp=true;
			}
			return comp;
		}
		weighted_sort(int *gc, SGrid *g, Box<int> *w):
		grid(g), walker(w)
		{
			grid_center[0]=gc[0];
			grid_center[1]=gc[1];
		}
	};

	//A custom sort object to sort the boxes first by type
	///within the same type, by box index sequentially
	struct seq_sort{
		SGrid *grid;
		bool operator()(Box<int> b1, Box<int> b2){
			bool comp=false;
							if(b1.type<b2.type) comp=true;
			else if (b1.type == b2.type) {
				int s1 = b1.min[0]+grid->max[0]*b1.min[1];
				int s2 = b2.min[0]+grid->max[0]*b2.min[1];
				if (s1<s2) comp=true;
			}
			return comp;

		}
		seq_sort(SGrid *g):grid(g)
		{}
	};

	//A custom sort by shortest distance from box to walker cell
	struct dist_sort{
		SGrid *grid;
		bool operator()(Box<int> b1, Box<int> b2){
			bool comp=false;
			if (b1.min_d<b2.min_d) comp = true;
			else if(b1.min_d==b2.min_d) {
									int s1 = b1.min[0]+grid->max[0]*b1.min[1];
									int s2 = b2.min[0]+grid->max[0]*b2.min[1];
				if (s1<s2) comp=true;
			}
			return comp;
		}
		dist_sort(SGrid *g):grid(g){}
	};
};//wl_gen

#endif
