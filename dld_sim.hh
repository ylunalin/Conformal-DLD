#ifndef DLA_SIM_HH
#define DLA_SIM_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <vector>
#include <limits>
#include <fstream>
#include "gsl/gsl_rng.h"
#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif

#include "common.hh"
#include "particle.hh"
#include "kdtree.hh"
#include "worklist.hh"
#include "sim_params.hh"

//typedef complex<double> Comp;
/** The width of particle */
static const double sfac=0.50;
/** The initial memory per each simulation block. */
const int init_region_memory=32;
/** The maximum memory that can be allocated to a particular simulation. */
const int max_region_memory=1048576;
/** The maximum number of blocs in one dimension.*/
const int min_blc_size=17;
class dld_sim {
	public:
		//Constructor and Destructor
		dld_sim(sim_params *spars);
		dld_sim();
		~dld_sim();

		//Dissolution and Aggregation
		/* Perform one single dissolution simulation. */
		void dissolve();
		/** A hybrid model of dld and dla. */
		int dl2a();
		/** Repeatedly dissolve, clean, reseed clusters. For radial geometry */
		void rapid_dissolve(double &x, double &y);

		// Helper function
		/* Use kdtree structure (TODO find out why it's so memory demanding.) to search for nearest neighbor.*/
		void test_search(bool swp);
		/*  Set the random number generator with a seed. */
		inline void set_rnd(int s){ gsl_rng_set(rando, s); }
		/* For multi realization simulations. */
		void clean(bool dla=false);
		/* Calculates the number of particles in the cluster were it be created.*/
		void calculate_cluster_size();

	//private:
		/** Flat geometry flag */
		bool is_flat;
		/** Total number of realizations */
		int tot_nr;
		/** The number of blocks in the x direction. */
		int nh;
		/** The number of blocks in the y direction. */
		int nv;
		/** The total number of blocks. */
		int mn;
		/** Random number generator seed */
		int seed;
		/** A counter for the total number of particles inserted */
		int n_part;
		/** A counter for total number of particles dissolved. */
		int n_dissd;
		/** Number of worklists */
		int nocell;
		/** Maximum sweeping distance in grid unit length */
		int max_swp_thr;
		/** the row of the lowest particle on flat interface*/
		int blclow;
		/** the length of the adld boundary shape array*/
		int adld_len;
		/** the tangent(tilt) calculation numerator */
		int par_n;
		/** the tangent(tilt) calculation denominator */
		int par_d;
		/** radius of the simulation cluster*/
		double rad;
		/** the coordinate of the righ edge, in a flat interface with periodic B.C.*/
		double redge;
		/** height of the simulation canvas, if not a square, could be change by flat()*/
		double height;
		/** total width */
		double twidth;
		/** True tilting angle after running continued fraction*/
		double true_tilt;
		/** length of a unit block of particles*/
		double ul;
		/** Tolerance of error for continued fraction angle estimate*/
		double cftol;
		/** Inverse of total number of realizations. */
		double tot_nrsp;
		/* 1 / (nocell*nocell), to normalize worklist threshold */
		double cell_normalizer;
		/* 1 / cell_length = nocell /dx */
		double dspcell;
		/** The length of a block in the x direction. */
		double dx;
		/** A small number to avoid numerical over/under flow*/
		double epsilon;
		/** The inverse length of a block in the x direction. */
		double xsp;
		/** Current radius of the cluster, sqaured*/
		double mrsq;
		/** Current radius of the cluster*/
		double mr;
		/** The distance above mr to generate the walkers */
		double abmr;
		/** Y coordinate of the lowest particle on flat interface*/
		double ylow;
		/** Y coordinate of the ceiling, will increase as blocks shift up */
		double yref;
		/** X coordinate of the lowest particle on flat interface*/
		double xlow;
		/** An array containing the number of particles per block. */
		int *co;
		/** An array containing the memory allocation per block. */
		int *mem;
		/** Holds the first position of each worklist in the precomputed worklist */
		int *wlpos;
		/** The name of the directory to store the output. */
		char *filename;
		/** Random number generator from gsl.*/
		gsl_rng *rando;
		/** A temporary character buffer used to create output
		 * filenames. */
		char *buf;
		/** Holds complex numbers the traces an adld shape*/
		Comp * shape;
		/** An array of particle information per block. */
		particle **ba;
		/** User defined threshold values */
		std::vector<int> thr;
		/** A quad- or kd- tree */
		KDTree mytree;
		/** precomputed worklist */
		std::vector<unsigned long int> wl;
		/** A simulation parameter object */
		sim_params *spars_p;
		/* A whole bunch of error msgs.*/
		static const char * LOC_NAN;
		static const char * TREE_NO_PBC;
		static const char * MINDIST_NO_PBC;
		static const char * NO_PBC;
		static const char * INSIDE_ERR;
		static const char * EQUI_DIST;
		static const char * FAR;
		static const char * TOO_FAR;
		static const char * NNB_SEARCH_FAIL;
		static const char * DIV_LINE;
		static const char * XSTEP_OUT_OF_BOUND;
		static const char * SPECIAL_SEARCH;
		static const char * HEX_ERR;
		static const char * SURVEY_ERR;
		static const char * TOO_FEW;
		static const char * BIN_MAT_ERR;
		static const char * INIT_FILE_ERR;
		static const char * PARAM_NOT_SET;
		// Constructor functions
		/** Used by the constructor to set block edge and size given radius*/
		void get_sim_params();
		void set_blc_size();
		void set_blc_size(double ** params, int* terms, int length, double ur, double dtheta);
		void set_blc_size(double tilt);
		void set_blc_size(int n);

		// House keeping

		// Helper funcs for all geometry
		void write(const char * k);
		void write(int k);
        void write_ADLD_shape(){
            sprintf(buf,"%s/shape.dat", spars_p->sim_dir);
            FILE *sh = safe_fopen(buf, "w");
            for (int j =0; j<adld_len; j++){
                fprintf(sh,"%.6f %.6f\n", real(shape[j]), imag(shape[j]));
            }
            fclose(sh);
        }
		int total_particles();
		void check_dup();
		particle * put(double x,double y);
		bool init_from_file(char* filename);
		void take(particle *pp);
		void add_region_memory(int s);
		void halve_region_memory(int s);
		void shrink_r();
		particle last();
		double bias_angle(double r);
		// Is this only good for non-periodic boundary condition?
		inline void int_box(double x,double y,double r,int &li,int &ui,int &lj,int &uj) {
			int ax = (-1)*rad;
			int ay = ax;
			li=int((x-ax-r)*xsp);if(li<0) li=0;if(li>=nh) li=nh-1;
			ui=int((x-ax+r)*xsp);if(ui<0) ui=0;if(ui>=nh) ui=nh-1;
			lj=int((y-ay-r)*xsp);if(lj<0) lj=0;if(lj>=nv) lj=nv-1;
			uj=int((y-ay+r)*xsp);if(uj<0) uj=0;if(uj>=nv) uj=nv-1;
		}
		inline double dis_sq(double gx,double gy,double qx,double qy) {
			double ex=gx-qx,ey=gy-qy;
			return ex*ex+ey*ey;
		}
		inline double rnd() {
			return gsl_rng_uniform(rando);
		}
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return double(clock())/CLOCKS_PER_SEC;}
#endif
		inline void msg(const char* func_name, const char* message){
			printf("%s: %s\n", func_name, message);
		}
		inline void fatal_error(const char* func_name, const char* message){
			printf("%s: Fatal error causes an exit. %s\n", func_name, message);
			exit(1);
		}
		inline void pos_of_rw(double dsq, bool &inside, bool &on_edge, bool &outside, bool &far_away){
			inside = on_edge= outside = far_away = false;

			inside = (dsq<1.0 - 4.0*epsilon);
			if(inside) return;

			on_edge = (dsq <=(1.+ 2*epsilon)*(1.+ 2*epsilon));
			if(on_edge) return;

			outside = (dsq<=9.0);
			if(outside) return;

			far_away = true;
		}

		// flat geometry
		bool ex_neighbor(double x, double y, short dir);
		void set_lowest(particle *pp);
		bool diss_enough();
		void pbcremap( double &x, double &y);
		// Put down particles in a cluster - periodic boundary condition, flat interface: here
		bool cont_frac(double deg, int& n, int& d, int par_lim);
		void shiftup(bool dupcheck,bool verbose);
		void flat(bool verbose);

		// ADLD shape
		struct vertex_pair{
			double v1[2];
			double v2[2];
			bool single;
			bool graze;
			vertex_pair(double x1, double y1, double x2, double y2):
			single(0), graze(0)
			{ v1[0] = x1; v1[1] =y1; v2[0]=x2; v2[1] = y2;};

			vertex_pair(double x1, double y1):
			single(1), graze(0)
			{ v1[0] = x1; v1[1] =y1; v2[0] =0; v2[1]= 0;};
			vertex_pair(){};

		};
		void adld_outline(Comp *&result, double ** params, int* terms, int length, double ur, double dtheta);
		void adld_cluster();

		// Circular cluster grid creation
		void square(double r);
		void hex(double r);


		// Dissolution and Aggregation
		// Common tools
		/** Repeatedly dissolve, clean, reseed clusters. For flat geometry. */
		particle * add_particle();
		particle * minus_particle();
		particle * min_distance(double x, double y, double &dsq, double &sec_d, double &rot_angle, bool verbose);
		particle * min_distance(double x, double y, int s, int stepx, int stepy, double &dsq, double &sec_d, bool verbose);

		//--Pacman
		bool test_isolation(particle *target);
		double pacman_tilt(double x, double y, particle *target);
		Comp pacman(double dsq, double sec_d, double rot_angle, particle * target, bool &eaten);

		//--Sweep
		bool exists(int s, int x, int y);
		particle * sweep(double &x, double y, double &dsq, double &sec_d, double &rot_angle,bool verbose);
		void wl_pos();
		void decoder(unsigned long int encoded, int &x_, int &y_, double &thr_);

		//--QuadTree
		Grid extract_grid();
		particle * treesearch(double x, double y, KDTree qt, double &dsq, double &sec_d, double &rot_angle);

		//--Stats
		void set_pdf_field(double **f, int k);
		double rgyration();
		void bump(double *cdf, double y, double h, double dh, int avail_n);
		double * survey(double dh, int &n, bool init = false);
		void active_zone_dist(int&nN, particle *p);
};
#endif
