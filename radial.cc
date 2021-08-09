/**
################################################################
Author: Y. Luna Lin
Email: y_lin@g.harvard.edu

Description: A program that creates circular cluster, or a shape
			 that is defined by the configuration file (see ADLD
			 paper on using Laurent series to define shapes).

			 The program then runs dissolution until only one
			 particle is left.

			 At the end, the averaged collapse coordinates are
			 recorded in "diss_pts" file; active_zone statistics
			 are recorded in "azone_N%d" file, identified by
			 value of N; and a configuration file for povray is
			 generated.
Parallel: Each thread owns its dld_sim object. And when the
		  threads divide up the realizations, it cleans its
		  simulation object for the next run.
 ###############################################################
*/
#include "dld_sim.hh"
#include "sim_manager.hh"

double ellipse_bn (const double aspect) {
    return (aspect-1) / (aspect+1);
}

int main(int argc, char *argv[]){
    double aspect = 1;
	int def_growth = false;
	int def_N_len = 30;
	int def_nr = 100;
	int def_st = 0;
#ifdef _OPENMP
	int def_thr = 10;
#else
	int def_thr = 1;
#endif
	double def_rad = 100;

    double ** params = new double* [2];
    params[0] = new double[2];
    params[1] = new double[2];
    params[0][0] = 1.0;
    params[0][1] = 0;
    params[1][1] = 0.;
    int length = 2;
    int * terms = new int[length];
    terms[0] = 1;
    terms[1] = -1;
    double dtheta = 0.01;


	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-h")==0){
			puts("Usage: ./radial <switches>");
            puts(" -a      \taspect ratio");
			puts(" -g      \tallows aggregation. Default false.");
			puts(" -h      \tprint this message.");
			puts(" -n <nr> \tset number of runs. Default 100.");
			puts(" -p <st> \tset simulation type. Default 0 (CIRCLE), 1 = ellipse.");
			puts(" -r <rad>\tset width of simulation. Default 100.");
			puts(" -t <thr>\tset number of OpenMP threads. Default 10.");
			puts(" -z <cfg>\tgive name of the active zone config file. TO BE IMPLEMENTED.");
			exit(0);
		}
		else if(strcmp(argv[i], "-n")==0) def_nr = (atof(argv[i+1]));
		else if(strcmp(argv[i], "-p")==0) def_st = (atof(argv[i+1]));
		else if(strcmp(argv[i], "-t") ==0) def_thr = (atof(argv[i+1]));
		else if(strcmp(argv[i], "-g") ==0) def_growth = true;
		else if(strcmp(argv[i], "-r") == 0) def_rad = (atof(argv[i+1]));
		else if(strcmp(argv[i], "-z")==0){
			// open the file for a list of N's, at which point to record prob(r,N)
			// set_prn(len, N);
		}
		else if(strcmp(argv[i], "-a") == 0) aspect = (atof(argv[i+1]));
	}
    if(def_st == sim_params::ADLD) {
        params[1][0] = ellipse_bn(aspect);
    }

	// creating a team of sim_params
	sim_params * spars= new sim_params[def_thr];
	for(int i=0;i<def_thr;i++){
		spars[i] = sim_params();
		spars[i].write=true;
		spars[i].sweep=true;
		// Set default values for common simulation parameters
		// argument list: radius, number of threads, number of realizations, number of workcells
		spars[i].set_common(def_rad, def_thr, def_nr, 1);
		spars[i].set_simtype(def_st);
		if(def_growth) spars[i].allow_growth();
        if(def_st == sim_params::ADLD) {
            spars[i].set_adld(params, terms, length, dtheta, aspect);
        }
	}

	// use the first spars pointer to get test cluter size
	dld_sim ** dlds = new dld_sim*[def_thr];
	for(int i=0;i<def_thr;i++){
		dlds[i] = new dld_sim(spars+i);
		dlds[i]->calculate_cluster_size();

	}

    double ntot_par = spars[0].tb_static;
	// create evenly spaced active zone sample points
	std::vector<int> ns;
	double start = log(2);
	double end = log(ntot_par);
	int numpts =5* (int) round(end+0.5);
	double incr = (end-start) / numpts;

    ns.push_back(ntot_par-1);
	int tn_prev =0, tn=2, cn=1;
	while(tn<ntot_par){
		if(tn!=tn_prev) ns.push_back(ntot_par - tn);
		tn_prev = tn;
		tn = int( exp(start + cn*incr));
		cn++;
	}

	def_N_len = ns.size();
	int * def_N = new int [def_N_len];
	for(int i=0;i<def_N_len;i++){
		def_N[i] = ns[def_N_len-1-i];
	}

	for(int i=0;i<def_thr;i++) {
		spars[i].set_N(def_N_len, def_N);
	}
	//=============================================================================================
	// Finish configuring sim_params //

	sim_manager manager(spars);
	//spars.set_circle(manager.interim);

	// Calculate where the left edge of the padded domain is
	manager.tic();
#pragma omp parallel for
	for (int k=0; k<spars->nr; k++){
#ifdef _OPENMP
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif
		double x=0, y=0;
		manager.checkin(k);
		manager.set_prn_r(k,tid);
		dlds[tid]->set_rnd(spars->randseed+k);
		dlds[tid]->rapid_dissolve(x, y);
		manager.record_clps(x,y);
		dlds[tid]->clean();
	}
	manager.toc();
	manager.end_report();
	manager.write_clps_dist();
	manager.write_povconfig();
    dlds[0]->write_ADLD_shape();
	manager.record_active_zone();

	// clean up
	for(int i=0;i<def_thr;i++) delete dlds[i];
	delete [] dlds;
	delete [] spars;
	delete [] def_N;
	delete [] params[0];
	delete [] params[1];
	delete [] params;
	delete [] terms;
}
