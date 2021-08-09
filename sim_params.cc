#include "sim_params.hh"

// Functions of sim_params class

void sim_params::die(const char * msg){
		printf("sim_params:: %s\n", msg);
		exit(1);
	}

// functions that are being used by other codes
void sim_params::false_all(){
	write = sweep = tree = verbose = local_verbose = hexa = dupcheck = false;
	aggre = common_set = adld_set = aflat_set = nflat_set = output_set = dl2a_set = false;
}
void sim_params::set_simdir(char * sd){
	sim_dir = sd;
}
void sim_params::set_simtype(int t){
	// 0 circle
	// 1 adld shapes
	// 2 flat geometry, use angle to initialize substrate
	// 3 flat geometry, use number of particle to initialize
	if(t<4 && t>=0) simtype = t;
	else die("sim type code ouf of range. 0-circle, 1-adld shapes, 2-aflat, 3-nflat.");
}
void sim_params::allow_growth(){
	aggre = true;
}
void sim_params::set_rad(double r){
	rad=r;
	randseed=(int) rad;
}
void sim_params::set_omp_thr(int t){
	omp_thr = t;
	reduce_thr_set();
}
void sim_params::set_nr(int n){
	nr=n;
	reduce_thr_set();
}
void sim_params::reduce_thr_set(){
	if(nr!=-1 && (simtype == CIRCLE || simtype == ADLD) && omp_thr>nr) omp_thr = nr;
	if(omp_thr<=0) omp_thr =1;
#ifdef _OPENMP
	omp_set_num_threads(omp_thr);
#endif
}
void sim_params::set_common(double rad_, int omp_thr_, int nr_, int nocell_, int N_len_, int *N_){
	set_rad(rad_);
	nr = nr_;
	omp_thr = omp_thr_;
	reduce_thr_set();
	nocell = nocell_;
	set_prn(N_, N_len_);
	common_set = true;
}
void sim_params::set_common(double rad_, int omp_thr_, int nr_, int nocell_){
	set_rad(rad_);
	nr = nr_;
	omp_thr = omp_thr_;
	reduce_thr_set();
	nocell = nocell_;
	N = NULL;
	common_set = true;
}
void sim_params::set_tilt(double tilt_){
	tilt = tilt_;
}
void sim_params::set_cftol(double cftol_){
	cftol = cftol_;
}
void sim_params::set_siden(int n){
	siden = n;
}
void sim_params::set_aflat(double tilt_, double cftol_, double dh_){
	if(simtype!=FLAT1) die("set_aflat() simtype mismatch.");
	set_tilt(tilt_);
	set_cftol(cftol_);
	dh = dh_;
	aflat_set = true;
}
void sim_params::set_nflat(int n, double dh_){
	if(simtype!=FLAT2) die("set_nflat() simtype mismatch.");
	set_siden(n);
	dh  = dh_;
	nflat_set = true;
}
void sim_params::set_flat(int n, double cftol_, double dh_){
	if(simtype==FLAT1) set_aflat( (double) n, cftol_, dh_);
	else if(simtype==FLAT2) set_nflat( n, dh_);
	else die("set_flat() simtype mismatch.");
}
// The following two functions used by managers most likely
void sim_params::set_flat_array(double **f_){
	if(simtype!=FLAT1 && simtype!=FLAT2) die("set_flat_array() simtype mismatch.");
	f = f_;
	output_set = true;
}
void sim_params::set_circle_array(double **f_){
	if(simtype!=CIRCLE) die("set_circle_array() simtype mismatch.");
	f = f_;
	output_set = true;
}

void sim_params::set_output_params(int tb){
	tb_static = tb;
}
void sim_params::set_output_params(double twidth, int nv, double dx, double angle_){
	true_width = twidth;
	f_rate = int(twidth);
	// compute burnin and target
	burnin = f_rate*f_rate;
	target = burnin*2;
	num_frame = f_rate + 1;
	// compute f_len (nelem)
	f_len = int (std::ceil( (nv*dx)/dh ));
	angle = angle_/pi * 180.;
}
void sim_params::set_dl2a(int sf, int ed, int target_, double agg_prob, double *f_, char* init_fn){
	if(!aggre) die("setting params for dl2a, but aggregation isn't turned on.");
	start_fr = sf;
	end_fr = ed;
	p_aggre = agg_prob;
	init_filename = init_fn;
	rg_data = f_;
	target = target_;
	dl2a_set = true;
}
void sim_params::set_adld(double **params_, int *terms_, int l, double t, double asp){
	if(simtype!=ADLD) die("set_adld() simtype mismatch.");
	params = params_;
	terms = terms_;
	length = l;
	dtheta = t;
    aspect_ratio = asp;
	adld_set = true;
}
// set parameter for active zone data collection
// number of N values, and the N values
// TODO could be passed in the same way as params files for ADLD shapes
void sim_params::set_prn(int *N_, int N_len_){
	N_len = N_len_;
	N = N_;
}
// set the pointer to the r data array, N_len points for any given realization
void sim_params::set_prn_r(double ** r_, double **r0_, double **x_, double **y_, int n){
	r_nr = r_[n];
	r0_nr = r0_[n];
	x_nr = x_[n];
	y_nr = y_[n];
}
