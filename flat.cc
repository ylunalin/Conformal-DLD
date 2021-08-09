/**
################################################################
Author: Y. Luna Lin
Email: y_lin@g.harvard.edu

Description: A program that creates flat periodic clusters and
			 runs dissolution until only a given number of
			 particles have been dissolved.
			 This process is repeated for nr times,
			 and statistics avergaed.
Parallel: Simulations for each angle has its own dld object,
		  which is cleaned at the end of a run for the next run.
 ###############################################################
*/
#include <sys/types.h>
#include <sys/stat.h>
#include "dld_sim.hh"
#include "sim_manager.hh"

int main(int argc, char *argv[]){
	int def_N_len = 1;
	int def_nr = 100;
	int max_n = 72;
	int incr_n = 3;
	int def_st = 3; // FLAT2
	int * def_N = new int [def_N_len];
	def_N[0] = 1000;

#ifdef _OPENMP
	int def_thr = 10;
#else
	int def_thr = 1;
#endif
	double def_rad = 100.0, def_dh = 0.01, def_cftol = 0.05;

	// Get the information to set up the parameter objects and manager objects
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-h")==0){
			puts("Usage: ./nflat <switches>");
			puts(" -c <inc>\tset the increment interval for tilt. Default 5.");
			puts(" -d <dh> \tset height interval to survey CDF. Default 0.01.");
			puts(" -f <tol>\tset continued fraction tolerance. Default 0.05.");
			puts(" -g      \tallows aggregation. Default false.");
			puts(" -h      \tprint this message.");
			puts(" -m <max>\tset maximum values for tilt (in degrees or number of particles). Default 50.");
			puts(" -n <nr> \tset number of runs. Default 100.");
			puts(" -p <st> \tset simulation type. Default 3 (FLAT2).");
			puts(" -r <rad>\tset width of simulation. Default 100.");
			puts(" -t <thr>\tset number of OpenMP threads. Default 10.");
			puts(" -z <cfg>\tgive name of the active zone config file. TO BE IMPLEMENTED.");
			exit(0);
		}
		else if(strcmp(argv[i], "-c")==0) incr_n = atoi(argv[i+1]);
		else if(strcmp(argv[i], "-d")==0) def_dh = atof(argv[i+1]);
		else if(strcmp(argv[i], "-f")==0) def_cftol = atof(argv[i+1]);
		else if(strcmp(argv[i], "-m")==0) max_n = atoi(argv[i+1]);
		else if(strcmp(argv[i], "-p")==0) def_st = (atoi(argv[i+1]));
		else if(strcmp(argv[i],"-t")==0) def_thr = atoi(argv[i+1]);
		else if(strcmp(argv[i], "-r") == 0) def_rad = (atof(argv[i+1]));
		else if(strcmp(argv[i], "-n")==0) def_nr = (atoi(argv[i+1]));
	}

	sim_params * spars = new sim_params[def_thr];
	sim_manager ** managers = new sim_manager*[def_thr];

	for(int i=0;i<def_thr;i++){
		managers[i] = NULL;
		spars[i] = sim_params();
		spars[i].write=true;
		spars[i].sweep=true;
		// Set default values for common simulation parameters
		// argument list: radius, number of threads, number of realizations, number of workcells
		spars[i].set_common(def_rad, def_thr, def_nr, 1);
		spars[i].set_simtype(def_st);
	}

	// Get the other command line customizations
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-g")==0){
			for(int j=0;j<def_thr;j++) spars[j].allow_growth();
		}
		else if(strcmp(argv[i], "-z")==0){
			// open the file for a list of N's, at which point to record prob(r,N)
			// set_prn(len, N);
		}
	}

	//=============================================================================================
	// Finish configuring sim_params and sim_manager

	int dn = int(max_n/incr_n);
	double twtime = 0;
#ifdef _OPENMP
		twtime = omp_get_wtime();
#else
		twtime = double (clock()) /CLOCKS_PER_SEC;
#endif
#pragma omp parallel for
	for(int i = 0; i<=dn; i++){
		int cn = max_n+i*incr_n;
		int tid = omp_get_thread_num();
		sim_params *sparp = spars+tid;
		sim_params &spar = *sparp;
		// important: set the numerator particle in parameters, flipp nflat_set to true
		spar.set_flat(cn, def_cftol, def_dh);

		// create the simulation, so that some system specific parameters are calculated and fed to sparp
		dld_sim dld(sparp);

		// create random sample points for active zone calc
		int N_len = (spar.target - spar.burnin) / spar.true_width;
		int * N = new int[N_len];
		//printf("burnin %d\n", spar.burnin);
		for(int n=0;n<N_len;n++){
			N[n] = spar.burnin + n*spar.true_width + (int) ( dld.rnd()*spar.true_width);
			//printf("%d N %d\n", n, N[n]);
		}
		spar.set_N(N_len, N);

		// create or reset the manager to have the right data array, output subdirectory, and update params
		if( managers[tid] == NULL) {
			managers[tid] = new sim_manager(spars+tid);
		}
		else managers[tid]->reset();
		sim_manager &manager = *managers[tid];

		manager.tic();
		for (int k=0; k<manager.nr; k++){
			double x=0, y=0;
			manager.checkin(k);
			manager.set_prn_r(k);
			dld.set_rnd((spar.randseed) + manager.nr*i + k);
			dld.rapid_dissolve(x, y);
			dld.clean();
		}
		manager.toc();

		manager.end_report();
		manager.record_active_zone();
		manager.record_cdf();
	}
#ifdef _OPENMP
		twtime = omp_get_wtime() - twtime;
#else
		twtime = double (clock()) /CLOCKS_PER_SEC -twtime;
#endif
	// clean up
	delete [] managers;
	delete [] spars;
	delete [] def_N;
	int hrs = int(twtime/3600);
	int mins = int((twtime - hrs*3600)/60);
	double secs = twtime - hrs*3600 - mins*60;
	printf("Program terminated normally. Took total %d hours %d mins %3.2f seconds\n", hrs, mins, secs);
}
