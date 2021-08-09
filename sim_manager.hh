#ifndef SIM_MANAGER_HH
#define SIM_MANAGER_HH

#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif
#include "sim_params.hh"

enum unit {
	none = '\0',
	kilo = 'k',
	mil = 'm'
};

class sim_manager{
	public:
	int nr;
	int grid_dim;
	int master_dim1;
	int master_dim2;
	int output_width;
	int ledge;
	int checkpt;
	int simtype;
	int tb_static;

	double padd;
	double nrsp;
	double padded_r;
	double time;
	double ** interim;
	/** Collapse point for cluster; CDF for flat */
	double ** master;
	/** Active zone distance, wrt center of mass */
	double ** r;
	/** Active zone distance, wrt origin */
	double ** r0;
	/** Active zone particle position, x and y */
	double **x;
	double **y;

	sim_params *spars;
	char *sim_dir;
	char *sub_sim_dir;
	char * buffer;

	sim_manager():
	interim(NULL), master(NULL), r(NULL), r0(NULL), x(NULL), y(NULL),
	spars(NULL), sim_dir(NULL), sub_sim_dir(NULL), buffer(NULL){}

	sim_manager(sim_params *spars_):
	spars(spars_), buffer(new char[256])
	{
		setup_manager();
	}

	~sim_manager(){
		if(x!=NULL){
			for (int i=0;i<nr;i++) delete [] x[i];
			delete [] x;
		}
		if(y!=NULL){
			for (int i=0;i<nr;i++) delete [] y[i];
			delete [] y;
		}
		if(r0!=NULL){
			for (int i=0;i<nr;i++) delete [] r0[i];
			delete [] r0;
		}
		if(r!=NULL){
			for (int i=0;i<nr;i++) delete [] r[i];
			delete [] r;
		}

		if(master!=NULL){
			for (int i=0;i<master_dim1;i++) delete [] master[i];
			delete [] master;
		}

		if(interim!=NULL){
			for (int i=0;i<4;i++) delete [] interim[i];
			delete [] interim;
		}

		if(sub_sim_dir!=NULL) delete [] sub_sim_dir;
		if(sim_dir!=NULL) delete [] sim_dir;
		if(buffer!=NULL) delete [] buffer;
	}

	inline void setup_manager(){
		padd = 4; // somewhat arbitrary
		double rad = spars->rad;
		padded_r = rad+padd;

		simtype = spars->simtype;

		int N_len = spars->N_len;
		nr = spars->nr;
		nrsp = 1.0 / (double) nr;
		if(nr==0) nrsp = 0;

		// set dimension of the master array
		// for radial dissolutions, it's the collapse point stats
		// for flat dissoltuions, it's the CDF stats
		grid_dim = ((int) rad)*2+1;
		int grid_dim_sq = grid_dim*grid_dim;
		master_dim1 = master_dim2 = int ((rad+3*padd))*2+1;
		if(simtype == spars->FLAT1 || simtype == spars->FLAT2) {
			master_dim2 = spars->f_len;
			master_dim1 = spars->num_frame;
		}
		else if(simtype == spars->CIRCLE || simtype == spars->ADLD){
			spars->num_frame = 4;
			spars->f_len = grid_dim_sq;
			// get staic cluster size from DLD
			tb_static = spars->tb_static;
		}

		// values needed to update the stat array and checkin
		ledge = -1*((int) (rad + 3*padd));
		if(nr<=100) checkpt = nr;
		else if(nr<=10000) checkpt = int(nr/100);
		else checkpt = int(nr/1000);
		if (checkpt<1) checkpt=1;



		sim_dir = new char[256];
		sub_sim_dir = new char[256];
		int nr_short = nr;
		unit u = none;
		if(nr>=1000 && nr<1000000) {
			nr_short = int(nr/1000);
			u = kilo;
		}
		else if(nr>=1000000){
			nr_short = int(nr/1000000);
			u = mil;
		}

		if(simtype == spars->CIRCLE) sprintf(sim_dir, "circle_r%.0f_nr%d%c",rad,nr_short,u);
		else if (simtype == spars->ADLD) sprintf(sim_dir, "adld_asp%.1f_r%.0f_nr%d%c",spars->aspect_ratio,rad,nr_short,u);
		else if (simtype == spars->FLAT1) sprintf(sim_dir, "aflat_w%.0f_nr%d%c",rad,nr_short,u);
		else if (simtype == spars->FLAT2) sprintf(sim_dir, "nflat_w%.0f_nr%d%c",rad,nr_short,u);
		else sprintf(sim_dir, "unknown_shape%.0f_nr%d",rad,nr);
		mkdir(sim_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

		// if necessary make subdirectory
		if(simtype == spars->FLAT1) {
			sprintf(sub_sim_dir, "%s/t%d", sim_dir, int(spars->tilt));
			mkdir(sub_sim_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
			spars->set_simdir(sub_sim_dir);
		}
		else if(simtype == spars->FLAT2) {
			sprintf(sub_sim_dir, "%s/n%d", sim_dir, spars->siden);
			mkdir(sub_sim_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
			spars->set_simdir(sub_sim_dir);
		}
		else {
			sprintf(sub_sim_dir, "%s", sim_dir);
			spars->set_simdir(sim_dir);
		}

		interim = new double *[4];
		for(int i=0;i<4;i++){
			interim[i] = new double [grid_dim_sq];
			for (int j=0; j<grid_dim_sq; j++) interim[i][j]=0;
		}

		master = new double *[master_dim1];
		for(int i=0;i<master_dim1;i++){
			master[i] = new double [master_dim2];
			for(int j=0;j<master_dim2;j++) master[i][j]=0;
		}

		// set CDF array for flat simulations
		if(simtype == spars->FLAT1 || simtype == spars->FLAT2){
			spars->set_flat_array(master);
			sprintf(buffer, "%s/size_angle", sim_dir);
			FILE *sa = safe_fopen(buffer, "w");
			fprintf(sa, "#input width\tinput n or tilt\ttrue width\ttrue tilt\n");
			fclose(sa);
			record_size_angle();
		}
		// not setting for radial simulations, left that to user

		r = new double *[nr];
		r0 = new double *[nr];
		x = new double *[nr];
		y = new double *[nr];
		for(int i=0;i<nr;i++){
			r[i] = new double[N_len];
			r0[i] = new double[N_len];
			x[i] = new double[N_len];
			y[i] = new double[N_len];
			for(int j=0;j<N_len;j++) {
				r[i][j] = 0;
				r0[i][j] = 0;
				x[i][j] = 0;
				y[i][j] = 0;
			}
		}

		// set active zone array
		set_prn_r(0);

		start_update();
	}
	void reset(){
		// TODO here r, r0, x, y must change as well since N would change.
		// as the angle changes, we reset the data array, output directory, and give an update
		if(simtype!=spars->FLAT1 &&simtype!=spars->FLAT2){
			puts("sim_manager:: Only flat simulations need to reset.");
			exit(1);
		}
		// must use old parameters to delete
		for(int i=0;i<master_dim1;i++) delete [] master[i];
		delete [] master;
		// get new parameters
		master_dim1 = spars->num_frame;
		master_dim2 = spars->f_len;

		// make new output subdirectory
		if(simtype == spars->FLAT1) {
			sprintf(sub_sim_dir, "%s/t%d", sim_dir, int(spars->tilt));
			mkdir(sub_sim_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
		}
		else if(simtype == spars->FLAT2) {
			sprintf(sub_sim_dir, "%s/n%d", sim_dir, spars->siden);
			mkdir(sub_sim_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
		}
		// reset data array
		master = new double *[master_dim1];
		for(int i=0;i<master_dim1;i++){
			master[i] = new double [master_dim2];
			for(int j=0;j<master_dim2;j++) master[i][j]=0;
		}
		spars->set_flat_array(master);

		int N_len = spars->N_len;
		for(int i=0;i<nr;i++) {
			double *temp_r = new double [N_len];
			delete [] r[i];
			r[i] = temp_r;

			double *temp_r0 = new double [N_len];
			delete [] r0[i];
			r0[i] = temp_r0;

			double *temp_x = new double [N_len];
			delete [] x[i];
			x[i] = temp_x;

			double *temp_y = new double [N_len];
			delete [] y[i];
			y[i] = temp_y;

			for(int j=0; j<N_len; j++){
				r[i][j] = 0;
				r0[i][j] = 0;
				x[i][j] = 0;
				y[i][j] = 0;
			}
		}

		set_prn_r(0);

		record_size_angle();
		start_update();
	}

#ifdef _OPENMP
	inline void tic() {
		time = omp_get_wtime();
	}
#else
	inline void tic(){
		time = double (clock()) /CLOCKS_PER_SEC;
	}
#endif

#ifdef _OPENMP
	inline void toc() {
		time = omp_get_wtime()-time;
	}
#else
	inline void toc(){
		time = double (clock()) /CLOCKS_PER_SEC - time;
	}
#endif

#ifdef _OPENMP
	inline void checkin(int k){
		if( (k+1)%checkpt ==0) printf("sim_manager:: thread #%d, realization %d.\n", omp_get_thread_num(), k+1);
	}
#else
	inline void checkin(int k){
		if( (k+1)%checkpt ==0) printf("sim_manager:: realization %d.\n", k+1);
	}
#endif

	void record_clps(double x, double y){
		int xpos = int(x) - ledge;
        int ypos = int(y) - ledge;
#pragma omp atomic
		master[ypos][xpos] += nrsp;
	}
	void write_clps_dist(){
		// Write the collaspe point distribution to file
		int lo=-1, hi=-1;
		int jmax=-1, imax = -1;
		for (int i=0;i<master_dim1;i++){
			lo = hi = -1;
			for(int j=0;j<master_dim2;j++){
				if( master[i][j] !=0 && lo == -1) lo = j;
				else if( master[i][j] ==0 && lo!=-1 && hi==-1) hi = j;
			}
			int curmax = hi-lo;
			if(curmax > jmax) jmax = curmax;
		}
		for (int j=0;j<master_dim2;j++){
			lo = hi = -1;
			for(int i=0;i<master_dim1;i++){
				if( master[i][j] !=0 && lo == -1) lo = i;
				else if( master[i][j] ==0 && lo!=-1 && hi==-1) hi = i;
			}
			int curmax = hi-lo;
			if(curmax > imax) imax = curmax;
		}
		if(imax<jmax) imax = jmax;
		output_width = int(imax/2) +2; // add 1 to account for rounding, add another 1 for safety
		sprintf(buffer,"%s/diss_pts", sim_dir);
		FILE * fg = safe_fopen(buffer, "w");
		// output centered at (0,0)
		for (int i =0; i<master_dim1; i++){
			if(i>-1*ledge-output_width && i< -1*ledge+output_width){
				for (int j=0; j<master_dim2; j++){
					if (j >-1*ledge-output_width && j <-1*ledge+output_width)
					fprintf(fg, "%.6f ", master[i][j]);
				}
				fprintf(fg, "\n");
			}
		}
		fclose(fg);
	}

	void write_povconfig(){
		sprintf(buffer,"%s/pov_config",sim_dir);
		FILE * fh = safe_fopen(buffer, "w");
		// XXX how does povray see xy dimension?
		fprintf(fh, "%d\n%d\n%d\n%d\n", master_dim1-1, master_dim2-1, 0, 0);
		fclose(fh);
	}

	void report_time(){
#ifdef _OPENMP
		int tid = omp_get_thread_num();
		if (simtype==spars->FLAT1 || simtype == spars->FLAT2) sprintf(buffer, "sim_manager:: manager %d\n", tid);
		else sprintf(buffer, "sim_manager:: ");
#else
		sprintf(buffer, "sim_manager:: ");
#endif
		if(time<=60.0) printf("%sDone in %3.1f seconds.\n", buffer, time);
		else if(time>60.0 && time<=3600.0){
			int mins = int(time/60.0);
			double seconds = time - mins*60.0;
			printf("%sDone in %d mins %3.1f seconds.\n", buffer, mins, seconds);
		}
		else if(time>3600){
			int hours = int(time/3600.0);
			int mins = int( (time - hours*3600.0) /60.0);
			double seconds = time - mins*60.0 - hours*3600.0;
			printf("%sDone in %d hours %d mins %3.1f seconds.\n", buffer, hours, mins, seconds);
		}
	}

	// set the pointers to use to record p(r,N) at a particular realization
	void set_prn_r(int n){
		spars->set_prn_r(r,r0,x,y,n);
	}

	void set_prn_r(int n, int ind){
		spars[ind].set_prn_r(r,r0,x,y,n);
	}

	void record_active_zone(){
		int N_len = spars-> N_len;
		double ** writeout = new double *[N_len];
		for(int i = 0;i<N_len;i++){
			writeout[i] = new double [nr];
			for(int j=0;j<nr;j++) {
				writeout[i][j] = r[j][i];
			}
		}
		// write out active zone distance wrt center of mass
		for(int i=0;i<N_len;i++){
			if(simtype==spars->CIRCLE || simtype == spars->ADLD) sprintf(buffer, "%s/azone_N%d", sub_sim_dir, tb_static-spars->N[i]);
			else sprintf(buffer, "%s/azone_N%d", sub_sim_dir, spars->N[i]);
			//printf("Filename %s\n", buffer);
			FILE * fh = safe_fopen(buffer, "wb");
			fwrite(writeout[i], sizeof(double), nr, fh);
			fclose(fh);
		}

		// write out active zone distance wrt origin
		if(simtype==spars->CIRCLE || simtype == spars->ADLD){
			for(int i = 0;i<N_len;i++){
				for(int j=0;j<nr;j++) {
					writeout[i][j] = r0[j][i];
				}
			}
			for(int i=0;i<N_len;i++){
				sprintf(buffer, "%s/azone_O%d", sub_sim_dir, tb_static-spars->N[i]);
				//printf("Filename %s\n", buffer);
				FILE * fh = safe_fopen(buffer, "wb");
				fwrite(writeout[i], sizeof(double), nr, fh);
				fclose(fh);
			}
		}

		// write out active zone position x and y
		for(int i = 0;i<N_len;i++){
			for(int j=0;j<nr;j++) {
				writeout[i][j] = x[j][i];
			}
		}
		for(int i=0;i<N_len;i++){
			if(simtype==spars->CIRCLE || simtype == spars->ADLD) sprintf(buffer, "%s/azone_X%d", sub_sim_dir, tb_static-spars->N[i]);
			else sprintf(buffer, "%s/azone_X%d", sub_sim_dir, spars->N[i]);
			//printf("Filename %s\n", buffer);
			FILE * fh = safe_fopen(buffer, "wb");
			fwrite(writeout[i], sizeof(double), nr, fh);
			fclose(fh);
		}

		// write out active zone position x and y
		for(int i = 0;i<N_len;i++){
			for(int j=0;j<nr;j++) {
				writeout[i][j] = y[j][i];
			}
		}
		for(int i=0;i<N_len;i++){
			if(simtype==spars->CIRCLE || simtype == spars->ADLD) sprintf(buffer, "%s/azone_Y%d", sub_sim_dir, tb_static-spars->N[i]);
			else sprintf(buffer, "%s/azone_Y%d", sub_sim_dir, spars->N[i]);
			//printf("Filename %s\n", buffer);
			FILE * fh = safe_fopen(buffer, "wb");
			fwrite(writeout[i], sizeof(double), nr, fh);
			fclose(fh);
		}

		for(int i=0; i<N_len;i++) delete [] writeout[i];
		delete [] writeout;
	}
	void record_size_angle(){
		// record in sub_sim_dir
#pragma omp critical
{
		sprintf(buffer, "%s/size_angle", sim_dir);
		FILE *sa = safe_fopen(buffer, "a");
		printf("true_angle %f\n", spars->angle);
		if(simtype == spars->FLAT2) fprintf(sa, "%.0f %d %.6f %.6f\n", spars->rad, spars->siden, spars->true_width, spars->angle);
		else if(simtype == spars->FLAT1) fprintf(sa, "%.0f %.0f %.6f %.6f\n", spars->rad, spars->tilt, spars->true_width, spars->angle);
		fclose(sa);
}
	}
	void record_cdf(){
		// as the angle changes, we reset the data array, output directory, and give an update
		if(simtype!=spars->FLAT1 &&simtype!=spars->FLAT2){
			puts("sim_manager:: Only flat simulations need to record cdf.");
			exit(1);
		}
		double dh = spars->dh;
		for(int k=0; k<master_dim1; k++){
			if(simtype == spars->FLAT1) sprintf(buffer, "%s/w.%d", sub_sim_dir,int(spars->tilt));
			else if(simtype == spars->FLAT2) sprintf(buffer, "%s/w.%d", sub_sim_dir, spars->siden);

			FILE *cdf = safe_fopen(buffer, "w");
			// TODO record nonzero might be faulty
			for (int l=0;l<master_dim2; l++){
				if(master[k][l]!=0.0) fprintf(cdf, "%g %1.10g\n", l*dh, master[k][l]*nrsp);
			}
			fclose(cdf);
		}
	}
	void start_update(){
		if(nr == 0 ) { puts("sim_manager:: nr set to zero, nothing will be done."); return; }

        printf("Output directory: %s\n", spars->sim_dir);
		if(simtype == spars->CIRCLE)
            printf("sim_manager:: circular simulation x%d runs;\n\tradius %.0f, padded radius = %.0f, padded total width = %.d.", nr, spars->rad, padded_r, master_dim1);
        else if(simtype == spars->ADLD) {
            printf("sim_manager:: ADLD cluster (ellipse) simulation x%d runs;\n\twidth%.0f, aspect_ratio %.1f\n", nr, spars->rad, (spars->params[1][0] + 1)/(1 - spars->params[1][0]) );
        }
		else if(simtype == spars->FLAT2)
            printf("sim_manager:: nflat simulation x%d runs;\n\tinput width %.0f true width %.6f; input numerator particle %d, true angle %.6f.\n", nr, spars->rad, spars->true_width, spars->siden, spars->angle);
		else puts("sim_manager:: update has not been implemented for this sim type.");
	}
	void end_report(){
		if(nr==0) return;
		report_time();
	}
};//end of class sim_manager

#endif
