/** A header file for the sim_params class, which contains
 * simulation parameters for the bac_sim class. Part of it
 * is adapted from fileinfo class in Chris's code.
 *
 * Author    : Y Luna Lin
 * Email     : y.lin2783@gmail.com
 * Date      : June, 2017
 */
#ifndef DLD_SIM_PARAMS
#define DLD_SIM_PARAMS

#include "omp.h"
#include "common.hh"
#include <cstring>
#include <vector>

// Parameters that can be set by dld_sim class are marked with          = D (dld)
// Parameters that need to be set by users, and have no default values  = U (user via config)
// Parameters that can be overwritten by users							= O (over via config)
// Parameters that are hard coded constants								= H (hardcoded)
// Other parameters														= M (misc)
class sim_params{
	public:
	// O booleans
	bool aggre;
	bool write;
	bool sweep;
	bool tree;
	bool verbose;
	bool local_verbose;
	bool hexa;
	bool dupcheck;

	// M booleans
	bool common_set;
	bool adld_set;
	bool aflat_set;
	bool nflat_set;
	bool output_set;
	bool dl2a_set;

	// H common between all simulations
	static const int safe_limit = 8192;
	static const int CIRCLE = 0;
	static const int ADLD = 1;
	static const int FLAT1 = 2;
	static const int FLAT2 = 3;

	// ===== common
	// write-out framerate
	int f_rate;//DO
	// used in all flat sims, and dla, dl2a sims
	int target;//DO
	int omp_thr;//O
	// number of cells to divide the center cell, for worklist
	int nocell;//O
	// number of realizations to do
	int nr;//U
	// random number generator seed
	int randseed; //O
	// simulation type: CIRCLE, ADLD, FLAT1, FLAT2
	int simtype; //U
	// length of N vector
	int N_len; //U

	// ===== radial
	// for radial geometry the total number of particles at the beginning
	int tb_static; //D


	// ===== flat
	// number of particle to discard at the beginning
	int burnin; //DO
	// length of discretized height
	int f_len; //D
	// number of frames to output
	int num_frame; //D
	// number of particles on the side
	int siden; //U

	// ===== dla
	int start_fr; //U
	int end_fr;//U

	// ===== adld
	int length;//U

	// ===== common
	double rad;//U

	// ===== flat
	// the true width of the simulation after cluster initiation
	double true_width;//D
	// height interval for CDF
	double dh;//U
	// the user input height of simulation block
	double height;//U
	// angle we input
	double tilt;//U
	// angle computed by dld
	double angle;//U
	// continued fraction calculation tolerance
	double cftol;//U

	// ===== dla
	// the probability of aggregation
	double p_aggre;//U

	// ===== adld
	double dtheta;//U

	// ===== common
	// arrays needed for active zone calculation
	int * N;///U through manager
	// if set in circular sims, output intermediate results
	// if set in flat sims, output interface CDF at every f_rate
	double **f;//U throguh manager
	// a pointer to store actiuve zone data from a particular realization
	double *r_nr;//U through manager
	double *r0_nr;//U through manager
	double *x_nr;//U through manager
	double *y_nr;//U through manager
	// simulation directory
	char * sim_dir;//U through manager

	// ===== dla
	// the data array to write radius of gyration stats to
	double *rg_data;//U through manager
	// the filename root of files to initializae dla clusters
	char * init_filename;//U

	// ===== adld
	int * terms;//U
	double ** params;//U
    double aspect_ratio;

	// Constructor
	sim_params():
	omp_thr(1), nocell(1), nr(-1), randseed(142857), tb_static(-1)
	{
		false_all();
		f = NULL;
		sim_dir = NULL;
		N = NULL;
		r_nr = NULL;
		rg_data = NULL;
		init_filename = NULL;
		terms = NULL;
		params = NULL;
        aspect_ratio = 1.0;
	}

	// Functions to set parameters
	void die(const char *msg);
	void false_all();
	void set_simdir(char *sd);
	void set_simtype(int t);
	void allow_growth();
	void set_rad(double r);
	void set_omp_thr(int t);
	void set_nr(int n);
	void reduce_thr_set();
	void set_common(double rad_, int omp_thr_, int nr_, int nocell_, int N_len_, int *N_);
	void set_common(double rad_, int omp_thr_, int nr_, int nocell_);
	inline void set_N(int N_len_, int *N_){
		N_len = N_len_;
		N = N_;
	}
	void set_tilt(double tilt_);
	void set_cftol(double cftol_);
	void set_siden(int n);
	void set_aflat(double tilt_, double cftol_, double dh_);
	void set_nflat(int n, double dh_);
	void set_flat(int n, double cftol_, double dh_);
	void set_flat_array(double **f_);
	void set_circle_array(double **f_);
	void set_output_params(int tb);
	void set_output_params(double twidth, int nv, double dx, double angle_);
	void set_dl2a(int sf, int ed, int target_, double agg_prob, double *f_, char* init_fn);
	void set_adld(double **params_, int *terms_, int l, double t, double dtheta);
	void set_prn(int *N_, int N_len_);
	void set_prn_r(double ** r_, double **r0_, double **x_, double **y_, int n);
	// Functions to parse config files
	// Copied from Chris's fileinfo.hh/.cc in ADLD file directory

	private:
		/** Tests to see if two strings are equal.
		 * \param[in] p1 a pointer to the first string.
		 * \param[in] p2 a pointer to the second string.
		 * \return True if they are equal, false otherwise. */
		inline bool se(const char *p1,const char *p2) {
			return strcmp(p1,p2)==0;
		}
		/** Finds the next token in a string and interprets it as a
		 * double precision floating point number. If none is availble,
		 * it gives an error message.
		 * \param[in] ln the current line number. */
		inline double next_double(int ln) {
			return atof(next_token(ln));
		}
		/** Finds the next token in a string, interprets it as a double
		 * precision floating point number, and checks that there are
		 * no subsequent values.
		 * \param[in] ln the current line number. */
		inline double final_double(int ln) {
			double temp=next_double(ln);
			check_no_more(ln);
			return temp;
		}
		/** Finds the next token in a string and interprets it as an
		 * integer. If none is availble, it gives an error message.
		 * \param[in] ln the current line number. */
		inline int next_int(int ln) {
			return atoi(next_token(ln));
		}
		/** Finds the next token in a string, interprets it as an
		 * integer, and checks that there are no subsequent values.
		 * \param[in] ln the current line number. */
		inline int final_int(int ln) {
			double temp=next_int(ln);
			check_no_more(ln);
			return temp;
		}
		char* next_token(int ln);
		void invalid_error(const char *cs,int ln);
		void coeff_error(const char* cs,int ln);
		void check_allocated(int ln);
		void check_no_more(int ln);
};
#endif
